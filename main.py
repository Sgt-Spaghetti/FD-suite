'''
  |----------------------------------|
  |  Written by Leonardo Cherin, UCL |
  |    Distributed as free software  | 
  |       with the GPL3 license      |
  |----------------------------------|
'''

import tkinter as tk
from tkinter import *
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
import h5py
import os

# A few useful global variables for the lifetime of the program
class GLOBALVARS():
	def __init__(self) -> None:
		self.all_files: list = []
		self.files_awaiting_selection: list = []
		self.selected_files: list = []
		self.groups = []
		self.active_file = None
		self.output_directory: str = ""
		self.graph_image = None
		self.extension_speed_um_s = 0.5
		self.frame_rate = 80
GLOBALVARS = GLOBALVARS()

# Create the "FD" class, to store all opened FD curves
class FD():
	def __init__(self, filepath: str) -> None:
		self.filepath: str = filepath
		self.name: str = os.path.splitext(os.path.basename(filepath))[0]
		self.group = None
		self.force_units: str = "pN"
		self.extension_units: str = "um"
		self.time_units: str = "sec"
		self.plot_time = False
		self.xmin: float = 0
		self.xmax: float = 0
		self.ymin: float = 0
		self.ymax: float = 0
		self.trimmed = False
		self.has_fit: bool = False
		self.baseline = False
		self.reference = False
		self.is_baseline_subtracted = False
		self.is_force_scaled = False
		self.current_plotted_trimmed: str = None # Holds "extension" or "retraction" keywords

		self.dataframe_extension = pd.DataFrame({"Force_Extension": [], "Distance_Extension": [], "Time_Extension": []})
		self.dataframe_retraction = pd.DataFrame({"Force_Retraction": [], "Distance_Retraction": [], "Time_Retraction": []})
		self.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": [], "Fit_Distance_Extension": []})
		self.fit_dataframe_retraction = pd.DataFrame({"Fit_Force_Retraction": [], "Fit_Distance_Retraction": []})

		self.fit_parameters = pd.DataFrame({"Lp_ext": [], "Lc_ext": [], "S_ext": [], "F0_ext": [],"Lp_ret": [], "Lc_ret": [], "S_ret": [], "F0_ret": []})

		self.initialise_attributes()

	def initialise_attributes(self) -> None:
		# NOTE: files are kept in RAM. If opening thousands this
		# might cause an issue, but it is highly unlikely
		raw_data = h5py.File(self.filepath, 'r')
		force_data: list[float] = raw_data["Force LF"]["Force 2x"]["Value"]
		distance_data: list[float] = raw_data["Distance"]["Distance 1"]["Value"]
		distance_to_time_conversion: list[float] = [(i/GLOBALVARS.frame_rate) for i in range(len(distance_data))]
	
		# Compute derivatives of the raw data, useful for data trimming
		derivatives: list = self.differentiate_savgol(distance_to_time_conversion, force_data, 2*GLOBALVARS.frame_rate, 2)
		first_derivative = [derivatives[0], derivatives[1]]
		second_derivative = [derivatives[0], derivatives[2]]

		# Initialise the starting core dataset, and a dummy "processed" dataset onto which any corrections (eg baseline subtraction) will be applied
		self.dataframe = pd.DataFrame({"Force": force_data, "Distance": distance_data, "Time": distance_to_time_conversion})
		self.processed_dataframe = pd.DataFrame({"Processed_Force": force_data, "Processed_Distance": distance_data, "Processed_Time":  distance_to_time_conversion})
		self.first_derivative_dataframe= pd.DataFrame({"First_Derivative": first_derivative[1], "Time": first_derivative[0], "Distance": distance_data})
		self.second_derivative_dataframe= pd.DataFrame({"Second_Derivative": second_derivative[1], "Time": second_derivative[0], "Distance": distance_data})

	# The central function responsible for plotting the curves, depending on what "state" the curve is in
	# For example, if it has been trimmed into extension / retraction curves, or fit with the eOdjik model
	def plot(self, expand_graph = False, extension_or_retraction = None) -> None:
		# If the framerate has been changed since last function call, update the global variable and adjust axis.
		if float(entry_frame_rate.get()) != GLOBALVARS.frame_rate:
			GLOBALVARS.frame_rate = float(entry_frame_rate.get())
			distance_to_time_conversion: list[float] = [(i/GLOBALVARS.frame_rate) for i in range(len(self.dataframe["Distance"]))]
			self.dataframe["Time"] = distance_to_time_conversion
		
		selected = False	
		for f in GLOBALVARS.selected_files:
			if self.name == f.name:
				selected = True
		if self.trimmed == False:
			if selected == True:
				self.plot_scatter_processed_data(self.xmin, self.xmax, self.ymin, self.ymax)
			else:
				self.plot_scatter_raw_data()
		else:
			if variable_checkbutton_view.get() == True:
				self.plot_scatter_processed_data(self.xmin, self.xmax, self.ymax)
			else:
				if self.has_fit == True:
					self.plot_fit_data(self.current_plotted_trimmed)
				else:
					self.plot_trimmed_data(self.current_plotted_trimmed)

		if expand_graph == True:
			plt.show()
			plt.close()
		else:
			plt.close()
			

	def plot_scatter_raw_data(self) -> None:
		window.update()
		window.update_idletasks()
		width: int = canvas_graph_display.winfo_width()
		height: int = canvas_graph_display.winfo_height()
		plt.figure(figsize=(width/100, height/100))
		plt.title(self.name)
		if self.plot_time == True:
			plt.scatter(self.dataframe["Time"], self.dataframe["Force"], s=0.1)
			plt.xlabel("Time (s)")
		else:
			plt.scatter(self.dataframe["Distance"], self.dataframe["Force"], s=0.1)
			plt.xlabel("Distance (um)")
		plt.ylabel("Force (pN)")
		plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_RAW_SCATTER.png"))
		plt.savefig("TEMP_PLOT.png")


	# Plot the data which might have been baseline subtracted or pixel corrected
	# Therefore, use the "processed_dataframe" dataset.
	def plot_scatter_processed_data(self, xmin=0, xmax=0, ymin=0, ymax=0) -> None:
		window.update()
		window.update_idletasks()
		width: int = canvas_graph_display.winfo_width()
		height: int = canvas_graph_display.winfo_height()
		
		fig, ax = plt.subplots(3,1,figsize=(width/100, height/100))
		fig.suptitle(self.name)
		if self.plot_time == True:
			ax[0].scatter(self.processed_dataframe["Processed_Time"], self.processed_dataframe["Processed_Force"], s=0.1)
			ax[1].plot(self.first_derivative_dataframe["Time"], self.first_derivative_dataframe["First_Derivative"])
			ax[2].plot(self.second_derivative_dataframe["Time"], self.second_derivative_dataframe["Second_Derivative"])
			ax[2].set_xlabel("Time (s)")
			if xmin > int((min(self.dataframe["Time"])*100)+0.5)/100:
				ax[0].axvline(xmin, color="r")
				ax[1].axvline(xmin, color="r")
				ax[2].axvline(xmin, color="r")
			if xmax > int((min(self.dataframe["Time"])*100)+0.5)/100:
				ax[0].axvline(xmax, color="g")
				ax[1].axvline(xmax, color="g")
				ax[2].axvline(xmax, color="g")
		else:
			ax[0].scatter(self.processed_dataframe["Processed_Distance"], self.processed_dataframe["Processed_Force"], s=0.1)
			ax[1].plot(self.first_derivative_dataframe["Distance"], self.first_derivative_dataframe["First_Derivative"])
			ax[2].plot(self.second_derivative_dataframe["Distance"], self.second_derivative_dataframe["Second_Derivative"])
			ax[2].set_xlabel("Distance(um)")
			if xmin > int((min(self.processed_dataframe["Processed_Distance"])*100)+0.5)/100:
				ax[0].axvline(xmin, color="r")
				ax[1].axvline(xmin, color="r")
				ax[2].axvline(xmin, color="r")
			if xmax > int((min(self.processed_dataframe["Processed_Distance"])*100)+0.5)/100:
				ax[0].axvline(xmax, color="g")
				ax[1].axvline(xmax, color="g")
				ax[2].axvline(xmax, color="g")
		if ymax > 0:
			ax[0].axhline(ymax, c="black")
		if ymin > 0:
			ax[0].axhline(ymin, c="blue")

		ax[1].axhline(0, c="black")
		ax[2].axhline(0, c="black")
		ax[0].set_ylabel("Force (pN)")
		ax[1].set_ylabel("dy/dx")
		ax[2].set_ylabel("ddy/dx")
		plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_PROCESSED_SCATTER.png"))
		plt.savefig("TEMP_PLOT.png")


	# Plot data which has been trimmed into extension or retraction curves
	def plot_trimmed_data(self, extension_or_retraction=None) -> None:
		if self.trimmed == True:
			extension_or_retraction = self.current_plotted_trimmed
			window.update()
			window.update_idletasks()
			width: int = canvas_graph_display.winfo_width()
			height: int = canvas_graph_display.winfo_height()
			plt.figure(figsize=(width/100, height/100))
			plt.title(self.name)
			if extension_or_retraction == "extension":
				if self.plot_time == True:
					plt.scatter(self.dataframe_extension["Time_Extension"], self.dataframe_extension["Force_Extension"], s=0.1)
					plt.xlabel("Time (s)")
				else:
					plt.scatter(self.dataframe_extension["Distance_Extension"], self.dataframe_extension["Force_Extension"], s=0.1)
					plt.xlabel("Distance (um)")
			else:
				if self.plot_time == True:
					plt.scatter(self.dataframe_retraction["Time_Retraction"], self.dataframe_retraction["Force_Retraction"], s=0.1)
					plt.xlabel("Time (s)")
				else:
					plt.scatter(self.dataframe_retraction["Distance_Retraction"], self.dataframe_retraction["Force_Retraction"], s=0.1)
					plt.xlabel("Distance (um)")

			plt.ylabel("Force (pN)")
			plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_TRIMMED_SCATTER.png"))
			plt.savefig("TEMP_PLOT.png")
			
	def plot_fit_data(self, extension_or_retraction=None) -> None:
		if self.trimmed == True:
			extension_or_retraction = self.current_plotted_trimmed
			window.update()
			window.update_idletasks()
			width: int = canvas_graph_display.winfo_width()
			height: int = canvas_graph_display.winfo_height()
			plt.figure(figsize=(width/100, height/100))
			plt.title(self.name)
			if extension_or_retraction == "extension":
				plt.scatter(self.dataframe_extension["Distance_Extension"], self.dataframe_extension["Force_Extension"], s=0.1)
				plt.plot(self.fit_dataframe_extension["Fit_Distance_Extension"], self.fit_dataframe_extension["Fit_Force_Extension"],c="r")
				plt.xlabel("Distance (um)")
			else:
				plt.scatter(self.dataframe_retraction["Distance_Retraction"], self.dataframe_retraction["Force_Retraction"], s=0.1)
				plt.plot(self.fit_dataframe_retraction["Fit_Distance_Retraction"], self.fit_dataframe_retraction["Fit_Force_Retraction"],c="r")
				plt.xlabel("Distance (um)")

			plt.ylabel("Force (pN)")
			plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_FIT.png"))
			plt.savefig("TEMP_PLOT.png")

	# Data is differentiated by applying a Savitzky-Golay filter to the raw noisy data,
	# And extracting the first and second derivatives directly from the differentiation of
	# the polynomial coefficients returned buy the Savitzky-Golay fit.
	# This is performed in the time domain, which guarantees evenly spaced data by
	# 1/framerate for the optical trap's camera system.
	# We will use a window size of 2*framerate to create a pseudo-2Hz lowpass filter, and we
	# will fit the window to a second degree polynomial.

	def differentiate_savgol(self, xdata, ydata, windowsize=160, degree=2):
		first_derivative = scipy.signal.savgol_filter(ydata, window_length=windowsize, polyorder=degree, mode="nearest", deriv=1)
		second_derivative = scipy.signal.savgol_filter(ydata, window_length=windowsize, polyorder=degree, mode="nearest", deriv=2)
		return [xdata, first_derivative, second_derivative]

	def subtract_baseline(self, baseline_x, baseline_y):
		if self.baseline == False:

			# Linear interpolation from baseline to datapoint xpos, fiding ypos to then subtract
			def apply_baseline(xdata, ydata):
				corrected_force = []
				for index, value in enumerate(xdata):
					corrected = False
					for i in range(len(baseline_x)-1):
						if baseline_x[i] <= value and baseline_x[i+1] > value:
							distance_before = baseline_x[i]
							distance_after = baseline_x[i+1]
							force_before = baseline_y[i]
							force_after = baseline_y[i+1]
							gradient = (force_after - force_before) / (distance_after - distance_before)
							corrected = True
					
					if corrected:
						corrected_force.append(ydata[index] - (force_before + (gradient *(value - distance_before))))
					else:
						corrected_force.append(ydata[index])
				return corrected_force

			self.is_baseline_subtracted = True

			# Baseline subtract every data point, from the full "processed_data" curve.
			corrected_force = apply_baseline(self.processed_dataframe["Processed_Distance"], self.processed_dataframe["Processed_Force"])
			self.processed_dataframe = pd.DataFrame({"Processed_Force": corrected_force, "Processed_Distance": self.processed_dataframe["Processed_Distance"], "Processed_Time": self.processed_dataframe["Processed_Time"]})
			if self.trimmed == True:
				# Subtract from the extension curve if applicable
				corrected_force = apply_baseline(self.dataframe_extension["Distance_Extension"], self.dataframe_extension["Force_Extension"])
				self.dataframe_extension = pd.DataFrame({"Force_Extension": corrected_force, "Distance_Extension": self.dataframe_extension["Distance_Extension"], "Time_Extension": self.dataframe_extension["Time_Extension"]})
				# Subtract from the retraction curve if applicable
				corrected_force = apply_baseline(self.dataframe_retraction["Distance_Retraction"], self.dataframe_retraction["Force_Retraction"])
				self.dataframe_extension = pd.DataFrame({"Force_Retraction": corrected_force, "Distance_Retraction": self.dataframe_retraction["Distance_Retraction"], "Time_Retraction": self.dataframe_retraction["Time_Retraction"]})

			
			# Compute derivatives of the raw data, useful for data trimming
			derivatives: list = self.differentiate_savgol(self.processed_dataframe["Processed_Time"], self.processed_dataframe["Processed_Force"], 2*GLOBALVARS.frame_rate, 2)
			first_derivative = [derivatives[0], derivatives[1]]
			second_derivative = [derivatives[0], derivatives[2]]

			self.first_derivative_dataframe= pd.DataFrame({"First_Derivative": first_derivative[1], "Time": first_derivative[0], "Distance": self.processed_dataframe["Processed_Distance"]})
			self.second_derivative_dataframe= pd.DataFrame({"Second_Derivative": second_derivative[1], "Time": second_derivative[0], "Distance": self.processed_dataframe["Processed_Distance"]})


		replot_canvas()


# Handle loading a folder of h5 files into the program
# NOTE: files are kept in RAM. If opening thousands this
# might cause an issue, but it is highly unlikely
def open_folder() -> list[str]:
	folder_path: str = filedialog.askdirectory()
	h5_files = []
	if folder_path != "":
		files: list[str] = os.listdir(folder_path)
		try:
			GLOBALVARS.output_directory = os.path.join(folder_path,"FDPLOT_OUTPUT")
			os.mkdir(GLOBALVARS.output_directory)
		except OSError as e:
			for f in os.listdir(GLOBALVARS.output_directory):
				os.remove(os.path.join(GLOBALVARS.output_directory,f))
			os.rmdir(GLOBALVARS.output_directory)
			os.mkdir(GLOBALVARS.output_directory)
		listbox_all_h5_files.delete(0, tk.END)
		for file in files:
			if os.path.splitext(file)[1] == ".h5":
				h5_files.append(file)

		sorted_names = sorted(h5_files)
		for file in sorted_names:
			FD_curve = FD(os.path.join(folder_path,file))
			GLOBALVARS.all_files.append(FD_curve)
			listbox_all_h5_files.insert(tk.END,FD_curve.name)	
		
		
# Handle selections in the first listbox
def all_h5_listbox_select(event) -> None:
	
	GLOBALVARS.files_awaiting_selection = []
	if len(listbox_all_h5_files.curselection()) > 0:
		# get the indecies of the selected files in the listbox
		indecies_highlighted_files_names: list[int] = listbox_all_h5_files.curselection()
		# if there is only one selection, plot it
		if len(indecies_highlighted_files_names) == 1:
			highlighted_file_name: str = listbox_all_h5_files.get(indecies_highlighted_files_names)
			# loop through all the curve objects, compare names to find the
			# matching one. Tell it to plot itself!
			for curve in GLOBALVARS.all_files:
				if highlighted_file_name == curve.name:
					GLOBALVARS.active_file = curve

					if curve.plot_time == False:
						scale_select_max_time.configure(from_ = min(curve.dataframe["Distance"]), to = max(curve.dataframe["Distance"]))
						scale_select_min_time.configure(from_ = min(curve.dataframe["Distance"]), to = max(curve.dataframe["Distance"]))
					else:
						scale_select_max_time.configure(from_ = min(curve.dataframe["Time"]), to = max(curve.dataframe["Time"]))
						scale_select_min_time.configure(from_ = min(curve.dataframe["Time"]), to = max(curve.dataframe["Time"]))
					
					GLOBALVARS.files_awaiting_selection.append(curve)
					scale_select_max_time.set(curve.xmax)
					scale_select_min_time.set(curve.xmin)
					curve.plot()
					update_canvas()

		else: # We have a multifile selection
			highlighted_files: list[str] = [listbox_all_h5_files.get(i) for i in indecies_highlighted_files_names]
			for curve in GLOBALVARS.all_files:
				for highlighted_file in highlighted_files:
					if highlighted_file == curve.name:
						GLOBALVARS.files_awaiting_selection.append(curve)
						if curve.name == highlighted_files[-1]:
							GLOBALVARS.active_file = curve
							if curve.plot_time == False:
								scale_select_max_time.configure(from_ = min(curve.dataframe["Distance"]), to = max(curve.dataframe["Distance"]))
								scale_select_min_time.configure(from_ = min(curve.dataframe["Distance"]), to = max(curve.dataframe["Distance"]))
							else:
								scale_select_max_time.configure(from_ = min(curve.dataframe["Time"]), to = max(curve.dataframe["Time"]))
								scale_select_min_time.configure(from_ = min(curve.dataframe["Time"]), to = max(curve.dataframe["Time"]))
							scale_select_max_time.set(0)
							scale_select_min_time.set(0)
							curve.plot()
							update_canvas()
	else:
		GLOBALVARS.files_awaiting_selection = []
						
	
# Handle selections in the second listbox
def all_selected_listbox_select(event) -> None:
	if len(listbox_all_selected_files.curselection()) > 0:
		index_highlighted_file_name: list[int] = listbox_all_selected_files.curselection()
		highlighted_file_name = listbox_all_selected_files.get(index_highlighted_file_name)
		# There is only one selection, plot it
		# loop through all selected curve objects, compare names to find the
		# matching one. Tell it to plot itself!
		for curve in GLOBALVARS.selected_files:
			if highlighted_file_name == curve.name:
				GLOBALVARS.active_file = curve
				if curve.plot_time == False:
					scale_select_max_time.configure(from_ = min(curve.dataframe["Distance"]), to = max(curve.dataframe["Distance"]))
					scale_select_min_time.configure(from_ = min(curve.dataframe["Distance"]), to = max(curve.dataframe["Distance"]))
				else:
					scale_select_max_time.configure(from_ = min(curve.dataframe["Time"]), to = max(curve.dataframe["Time"]))
					scale_select_min_time.configure(from_ = min(curve.dataframe["Time"]), to = max(curve.dataframe["Time"]))
				scale_select_max_time.set(curve.xmax)
				scale_select_min_time.set(curve.xmin)
				curve.plot()
				update_canvas()

def toggle_time() -> None:
	if GLOBALVARS.active_file != None:
		max_time: float = max(GLOBALVARS.active_file.dataframe["Time"])
		min_time: float = min(GLOBALVARS.active_file.dataframe["Time"])
		max_d: float = max(GLOBALVARS.active_file.dataframe["Distance"])
		min_d: float = min(GLOBALVARS.active_file.dataframe["Distance"])
		if GLOBALVARS.active_file.plot_time == False:
			GLOBALVARS.active_file.plot_time = True
			scale_select_max_time.configure(from_ = min_time, to = max_time)
			scale_select_min_time.configure(from_ = min_time, to = max_time)
			GLOBALVARS.active_file.xmin = (((GLOBALVARS.active_file.xmin-min_d)/(max_d-min_d))*(max_time-min_time)) + min_time
			GLOBALVARS.active_file.xmax = (((GLOBALVARS.active_file.xmax-min_d)/(max_d-min_d))*(max_time-min_time)) + min_time

		else:
			GLOBALVARS.active_file.plot_time = False
			scale_select_max_time.configure(from_ = min_d, to = max_d)
			scale_select_min_time.configure(from_ = min_d, to = max_d)
			GLOBALVARS.active_file.xmin = ((GLOBALVARS.active_file.xmin-min_time)/(max_time-min_time))*(max_d-min_d) + min_d
			GLOBALVARS.active_file.xmax = ((GLOBALVARS.active_file.xmax-min_time)/(max_time-min_time))*(max_d-min_d) + min_d

		GLOBALVARS.active_file.ymax = float(entry_ymax.get())
		scale_select_max_time.set(GLOBALVARS.active_file.xmax)
		scale_select_min_time.set(GLOBALVARS.active_file.xmin)
		GLOBALVARS.active_file.plot()
		replot_canvas()

def slider_max_release(event) -> None:
	if GLOBALVARS.active_file != None:
		GLOBALVARS.active_file.xmax = scale_select_max_time.get()
		GLOBALVARS.active_file.ymin = float(entry_ymin.get())
		GLOBALVARS.active_file.ymax = float(entry_ymax.get())
		entry_xmax.delete(0, tk.END)
		entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax))
		GLOBALVARS.active_file.plot()
		update_canvas()

def slider_min_release(event) -> None:
	if GLOBALVARS.active_file != None:
		GLOBALVARS.active_file.xmin = scale_select_min_time.get()
		GLOBALVARS.active_file.ymin = float(entry_ymin.get())
		GLOBALVARS.active_file.ymax = float(entry_ymax.get())
		entry_xmin.delete(0, tk.END)
		entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin))
		replot_canvas()

def update_trim_settings() -> None:
	if GLOBALVARS.active_file != None:
		GLOBALVARS.active_file.xmin = float(entry_xmin.get())
		GLOBALVARS.active_file.xmax = float(entry_xmax.get())
		GLOBALVARS.active_file.ymin = float(entry_ymin.get())
		GLOBALVARS.active_file.ymax = float(entry_ymax.get())
		replot_canvas()

def update_canvas() -> None:
	GLOBALVARS.graph_image = tk.PhotoImage(file="TEMP_PLOT.png")
	canvas_graph_display.create_image(0,0,image=GLOBALVARS.graph_image, anchor="nw")

def replot_canvas(expanded_graph = False) -> None:
	if GLOBALVARS.active_file != None:
		GLOBALVARS.active_file.plot(expanded_graph)
		GLOBALVARS.graph_image = tk.PhotoImage(file="TEMP_PLOT.png")
		canvas_graph_display.create_image(0,0,image=GLOBALVARS.graph_image, anchor="nw")

		if GLOBALVARS.active_file.has_fit == True:
			if GLOBALVARS.active_file.current_plotted_trimmed == "extension":
				Lp_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lp_ext"][0])))
				Lc_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lc_ext"][0])))
				S_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["S_ext"][0])))
				F0_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["F0_ext"][0])))
			else:
				Lp_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lp_ret"][0])))
				Lc_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lc_ret"][0])))
				S_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["S_ret"][0])))
				F0_display.config(text=str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["F0_ret"][0])))

def window_resize(event) -> None:
	replot_canvas()

def update_optic_settings() -> None:
	replot_canvas()
	if GLOBALVARS.active_file != None:
		GLOBALVARS.active_file.ymax = float(entry_ymax.get())
		if GLOBALVARS.active_file.plot_time == True:
			scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Time"]), to = max(GLOBALVARS.active_file.dataframe["Time"]))
			scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Time"]), to = max(GLOBALVARS.active_file.dataframe["Time"]))
		else:
			scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Distance"]), to = max(GLOBALVARS.active_file.dataframe["Distance"]))
			scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Distance"]), to = max(GLOBALVARS.active_file.dataframe["Distance"]))

def add_selected_curves() -> None:
	for curve in GLOBALVARS.files_awaiting_selection:
		GLOBALVARS.selected_files.append(curve)
		listbox_all_selected_files.insert(tk.END, curve.name)
	GLOBALVARS.files_awaiting_selection = []

def deselect_curves() -> None:
	# The listbox only has curve names, not the actual object
	current_selection = listbox_all_selected_files.get(listbox_all_selected_files.curselection())
	old_selection = GLOBALVARS.selected_files
	updated_selection = []
	for curve in old_selection:
		if curve.name != current_selection:
			updated_selection.append(curve)
	GLOBALVARS.selected_files = updated_selection
	listbox_all_selected_files.delete(0,tk.END)
	for curve in updated_selection:
		listbox_all_selected_files.insert(tk.END, curve.name)

def group_curves() -> None:
	pass

def auto_trim_data() -> None:
	if GLOBALVARS.active_file in GLOBALVARS.selected_files:
		
		inflection_point = 0
		max_force = 0
		force_cap = float(entry_ymax.get())
		for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])-1):
			if GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] > max_force and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] > GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i+1]:	
				max_force = GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i]
				inflection_point = i
		force_ext = []
		dist_ext= []
		time_ext = []
		force_ret = []
		dist_ret = []
		time_ret = []

		for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][0:inflection_point])):
			if GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] < force_cap:
				force_ext.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
				time_ext.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
				dist_ext.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
		for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][inflection_point+1:])):
			if GLOBALVARS.active_file.processed_dataframe["Processed_Force"][inflection_point+i+1] < force_cap:
				force_ret.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][inflection_point+i+1])
				time_ret.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][inflection_point+i+1])
				dist_ret.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][inflection_point+i+1])

		GLOBALVARS.active_file.dataframe_extension = pd.DataFrame({"Force_Extension": force_ext, "Distance_Extension": dist_ext, "Time_Extension": time_ext})
		GLOBALVARS.active_file.dataframe_retraction = pd.DataFrame({"Force_Retraction": force_ret, "Distance_Retraction": dist_ret, "Time_Retraction": time_ret})

		GLOBALVARS.active_file.trimmed=True
		GLOBALVARS.active_file.current_plotted_trimmed="extension"
		GLOBALVARS.active_file.plot("extension")
		replot_canvas()

def manual_trim_data() -> None:
	if GLOBALVARS.active_file in GLOBALVARS.selected_files:
		trimmed_force = []
		trimmed_time = []
		trimmed_dist= []
		if str(variable_radio_buttons.get()) == "extension":
			if GLOBALVARS.active_file.plot_time == True:
				for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])):
					if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] >= GLOBALVARS.active_file.xmin and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmax and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax:
						trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
						trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
						trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
				GLOBALVARS.active_file.dataframe_extension = pd.DataFrame({"Force_Extension": [], "Distance_Extension": [], "Time_Extension": []})
				GLOBALVARS.active_file.dataframe_extension["Force_Extension"] = trimmed_force
				GLOBALVARS.active_file.dataframe_extension["Time_Extension"] = trimmed_time
				GLOBALVARS.active_file.dataframe_extension["Distance_Extension"] = trimmed_dist

			else:
				for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])):
					if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] >= GLOBALVARS.active_file.xmin and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] <= GLOBALVARS.active_file.xmax and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax:
						trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
						trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
						trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
				GLOBALVARS.active_file.dataframe_extension = pd.DataFrame({"Force_Extension": [], "Distance_Extension": [], "Time_Extension": []})
				GLOBALVARS.active_file.dataframe_extension["Force_Extension"] = trimmed_force
				GLOBALVARS.active_file.dataframe_extension["Time_Extension"] = trimmed_time
				GLOBALVARS.active_file.dataframe_extension["Distance_Extension"] = trimmed_dist

			GLOBALVARS.active_file.trimmed=True
			GLOBALVARS.active_file.current_plotted_trimmed="extension"
			GLOBALVARS.active_file.plot("extension")
			replot_canvas()

		elif str(variable_radio_buttons.get()) == "retraction":
			if GLOBALVARS.active_file.plot_time == True:
				for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])):
					if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] >= GLOBALVARS.active_file.xmin and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmax and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax:
						trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
						trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
						trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
				GLOBALVARS.active_file.dataframe_retraction = pd.DataFrame({"Force_Retraction": [], "Distance_Retraction": [], "Time_Retraction": []})
				GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"] = trimmed_force
				GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"] = trimmed_time
				GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"] = trimmed_dist
			else:
				for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])):
					if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] >= GLOBALVARS.active_file.xmin and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] <= GLOBALVARS.active_file.xmax and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax:
						trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
						trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
						trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
				GLOBALVARS.active_file.dataframe_retraction = pd.DataFrame({"Force_Retraction": [], "Distance_Retraction": [], "Time_Retraction": []})
				GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"] = trimmed_force
				GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"] = trimmed_time
				GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"] = trimmed_dist

			GLOBALVARS.active_file.trimmed=True
			GLOBALVARS.active_file.current_plotted_trimmed="retraction"
			GLOBALVARS.active_file.plot("retraction")
			replot_canvas()

def radio_button_select() -> None:
	if variable_checkbutton_view.get() != True:
		if GLOBALVARS.active_file.trimmed == True:
			GLOBALVARS.active_file.current_plotted_trimmed = variable_radio_buttons.get()
			replot_canvas()
	else:
		replot_canvas()


def fit_eOdijk_F0(data_x, data_y) -> list:
	def eOdjik_force_offset(d, Lp=50, Lc=16.5, S=1500, F0=0, kT=4.11):
		output = (2*(Lp*Lc*S*d - Lp*S*(Lc**2)))/(3*Lp*(Lc**2)) - \
		   (-16*(Lp**2)*(S**2)*(d**2)*(Lc**2) + \
		    32*(Lp**2)*(S**2)*d*(Lc**3) - 16*(Lp**2)*(S**2)*(Lc**4))/ \
		     (24*Lp*(Lc**2)*(-8*(Lp**3)*(S**3)*(d**3)*(Lc**3) + \
		       24*(Lp**3)*(S**3)*(d**2)*(Lc**4) - \
			     24*(Lp**3)*(S**3)*d*(Lc**5) + \
		       27*kT*(Lp**2)*(S**2)*(Lc**6) + 8*(Lp**3)*(S**3)*(Lc**6) + \
		       3*np.sqrt(3)* \
			np.sqrt(-16*kT*(Lp**5)*(S**5)*(d**3)*(Lc**9) + \
			  48*kT*(Lp**5)*(S**5)*(d**2)*(Lc**10) - \
				   48*kT*(Lp**5)*(S**5)*d*(Lc**11) + \
			  27*(kT**2)*(Lp**4)*(S**4)*(Lc**12) + \
				   16*kT*(Lp**5)*(S**5)*(Lc**12)))**(1/3)) + \
		      (1/(6*Lp*(Lc**2)))* \
		     (-8*(Lp**3)*(S**3)*(d**3)*(Lc**3) + \
		     24*(Lp**3)*(S**3)*(d**2)*(Lc**4) - \
		     24*(Lp**3)*(S**3)*d*(Lc**5) + \
			  27*kT*(Lp**2)*(S**2)*(Lc**6) + \
		     8*(Lp**3)*(S**3)*(Lc**6) + \
		     3*np.sqrt(3)* \
		      np.sqrt(-16*kT*(Lp**5)*(S**5)*(d**3)*(Lc**9) + \
			48*kT*(Lp**5)*(S**5)*(d**2)*(Lc**10) - \
				48*kT*(Lp**5)*(S**5)*d*(Lc**11) + \
			27*(kT**2)*(Lp**4)*(S**4)*(Lc**12) + \
				16*kT*(Lp**5)*(S**5)*(Lc**12)))**(1/3)
		return output + F0

	parameters, covariance = scipy.optimize.curve_fit(eOdjik_force_offset,np.array(data_x),np.array(data_y), (float(Lp_entry.get()), float(Lc_entry.get()), float(S_entry.get()), float(F0_entry.get())))
	errors = np.sqrt(np.diag(covariance))
	predicted_force = eOdjik_force_offset(np.array(data_x), *parameters)
	residuals = np.array(data_y) - predicted_force
	return [parameters, errors, predicted_force, residuals]


def fit() -> None:
	if GLOBALVARS.active_file in GLOBALVARS.selected_files:
		if GLOBALVARS.active_file.trimmed == True:
			if GLOBALVARS.active_file.current_plotted_trimmed == "extension":
				fit_result = fit_eOdijk_F0(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"], GLOBALVARS.active_file.dataframe_extension["Force_Extension"])
				GLOBALVARS.active_file.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": [], "Fit_Distance_Extension": []})
				GLOBALVARS.active_file.fit_dataframe_extension["Fit_Force_Extension"] = fit_result[2]
				GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"] = GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]
				GLOBALVARS.active_file.fit_parameters["Lp_ext"] = [fit_result[0][0], fit_result[1][0]]
				GLOBALVARS.active_file.fit_parameters["Lc_ext"] = [fit_result[0][1], fit_result[1][1]]
				GLOBALVARS.active_file.fit_parameters["S_ext"] = [fit_result[0][2], fit_result[1][2]]
				GLOBALVARS.active_file.fit_parameters["F0_ext"] = [fit_result[0][3], fit_result[1][3]]

			else:
				fit_result = fit_eOdijk_F0(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"], GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"])
				GLOBALVARS.active_file.fit_dataframe_retraction = pd.DataFrame({"Fit_Force_Retraction": [], "Fit_Distance_Retraction": []})
				GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Force_Retraction"] = fit_result[2]
				GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"] = GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]
				GLOBALVARS.active_file.fit_parameters["Lp_ret"] = [fit_result[0][0], fit_result[1][0]]
				GLOBALVARS.active_file.fit_parameters["Lc_ret"] = [fit_result[0][1], fit_result[1][1]]
				GLOBALVARS.active_file.fit_parameters["S_ret"] = [fit_result[0][2], fit_result[1][2]]
				GLOBALVARS.active_file.fit_parameters["F0_ret"] = [fit_result[0][3], fit_result[1][3]]

			GLOBALVARS.active_file.has_fit = True
			GLOBALVARS.active_file.plot()
			replot_canvas()

def mark_baseline() -> None:
	GLOBALVARS.active_file.baseline = True

def mark_reference() -> None:
	GLOBALVARS.active_file.reference = True

def baseline_correction():
	minimum = 100
	maximum = 0
	densist_curve = 0
	for curve in GLOBALVARS.selected_files:
		if curve.baseline == True:
			data = curve.dataframe["Distance"]
			if len(data) > densist_curve:
				densist_curve = len(data)
			if min(data) < minimum:
				minimum = min(data)
			if max(data) > maximum:
				maximum = max(data)
	baseline_x = np.linspace(minimum, maximum, densist_curve)
	baseline_y = []
	last_added = 0
	for i in range(len(baseline_x) - 1):
		points = []
		for curve in GLOBALVARS.selected_files:
			if curve.baseline == True:
				for index, value in enumerate(curve.dataframe["Distance"]):
					if value >= baseline_x[i] and value < baseline_x[i+1]:
						points.append(curve.dataframe["Force"].iloc[index])
		if points != []:
			mean = np.mean(points)
			baseline_y.append(mean)
			last_added = mean
		else:
			baseline_y.append(last_added)
	
	for curve in GLOBALVARS.selected_files:
		curve.subtract_baseline(baseline_x[:-1], baseline_y)

	replot_canvas()

# For every reference curve, assume the curve is trimmed nicely to have one extension ready for fitting
# fit every reference curve, and average each Lc. From the average Lc, convert each curve to Lc space.
# Find the force at critical Lc from the full curve (equivalent to 22.15um in lambda), the force at this
# critical Lc should be 110pN. Get a correcton factor for each curve. Average the correction factors for
# each curve, then apply the averaged correction factors to every selected curve (via processed_data).
def auto_force_scale() -> None:

	reference_curves = []
	contour_lengths = []
	for curve in GLOBALVARS.selected_files:
		if curve.reference == True:
			reference_curves.append(curve)

	for curve in reference_curves:
		if curve.trimmed == False: # if the curve is not nicely trimmed already, attempt to do it automatically. Copy of the function before.
			if curve in GLOBALVARS.selected_files:	
				inflection_point = 0
				max_force = 0
				force_cap = float(entry_ymax.get())
				for i in range(len(curve.processed_dataframe["Processed_Time"])-1):
					if curve.processed_dataframe["Processed_Force"][i] > max_force and curve.processed_dataframe["Processed_Force"][i] > curve.processed_dataframe["Processed_Force"][i+1]:	
						max_force = curve.processed_dataframe["Processed_Force"][i]
						inflection_point = i
				force_ext = []
				dist_ext= []
				time_ext = []
				force_ret = []
				dist_ret = []
				time_ret = []

				for i in range(len(curve.processed_dataframe["Processed_Force"][0:inflection_point])):
					if curve.processed_dataframe["Processed_Force"][i] < force_cap:
						force_ext.append(curve.processed_dataframe["Processed_Force"][i])
						time_ext.append(curve.processed_dataframe["Processed_Time"][i])
						dist_ext.append(curve.processed_dataframe["Processed_Distance"][i])
				for i in range(len(curve.processed_dataframe["Processed_Force"][inflection_point+1:])):
					if curve.processed_dataframe["Processed_Force"][inflection_point+i+1] < force_cap:
						force_ret.append(curve.processed_dataframe["Processed_Force"][inflection_point+i+1])
						time_ret.append(curve.processed_dataframe["Processed_Time"][inflection_point+i+1])
						dist_ret.append(curve.processed_dataframe["Processed_Distance"][inflection_point+i+1])

				curve.dataframe_extension = pd.DataFrame({"Force_Extension": force_ext, "Distance_Extension": dist_ext, "Time_Extension": time_ext})
				curve.dataframe_retraction = pd.DataFrame({"Force_Retraction": force_ret, "Distance_Retraction": dist_ret, "Time_Retraction": time_ret})

				curve.trimmed=True
				curve.current_plotted_trimmed="extension"

		if curve.has_fit == False: # if the curve is not already fit, then fit the trimmed data.
			fit_result = fit_eOdijk_F0(curve.dataframe_extension["Distance_Extension"], curve.dataframe_extension["Force_Extension"])
			curve.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": [], "Fit_Distance_Extension": []})
			curve.fit_dataframe_extension["Fit_Force_Extension"] = fit_result[2]
			curve.fit_dataframe_extension["Fit_Distance_Extension"] = curve.dataframe_extension["Distance_Extension"]
			curve.fit_parameters["Lp_ext"] = [fit_result[0][0], fit_result[1][0]]
			curve.fit_parameters["Lc_ext"] = [fit_result[0][1], fit_result[1][1]]
			curve.fit_parameters["S_ext"] = [fit_result[0][2], fit_result[1][2]]
			curve.fit_parameters["F0_ext"] = [fit_result[0][3], fit_result[1][3]]

			curve.has_fit = True
		contour_lengths.append(curve.fit_parameters["Lc_ext"][0])

	mean_Lc = np.mean(contour_lengths)
	force_corrections = []
	
	for curve in reference_curves:
		Lc_space = curve.processed_dataframe["Processed_Distance"] / mean_Lc # Normalise to Lc space
		critical_Lc = 1.342424 # lambda Lc=16.5um. 22.15um = 1.342424x Lc. 1.342424 LC = 110 pN.
		force_at_critical_lc = 1
		for i in range(len(Lc_space)-1):
			if Lc_space[i] <= critical_Lc and Lc_space[i+1] > critical_Lc:
				force_at_critical_lc = (curve.processed_dataframe["Processed_Force"][i] + curve.processed_dataframe["Processed_Force"][i+1])/2
		force_correction_factor = 110/force_at_critical_lc
		force_corrections.append(force_correction_factor)

	mean_force_correction = np.mean(force_corrections)

	# apply the force correction to every selected curve.
	for curve in GLOBALVARS.selected_files:
		curve.processed_dataframe["Processed_Force"] = curve.processed_dataframe["Processed_Force"] * mean_force_correction
		curve.is_force_scaled = True

	replot_canvas()
			

def expand_graph() -> None:
	replot_canvas(True)

def export_data() -> None:
	fit_params: dict = {"File": [], "Lp-e": [], "Lc-e": [], "S-e": [], "F0-e": [], "Lp-r": [], "Lc-r": [], "S-r": [], "F0-r": []}
	for file in GLOBALVARS.selected_files:
		output_data = pd.concat([file.dataframe, file.processed_dataframe, file.dataframe_extension, file.dataframe_retraction, file.first_derivative_dataframe, file.second_derivative_dataframe, file.fit_dataframe_extension, file.fit_dataframe_retraction, file.fit_parameters, pd.DataFrame({"is_force_scaled": [file.is_force_scaled]})], axis = 1)
		output_data.to_csv(os.path.join(GLOBALVARS.output_directory, file.name+".csv"), index=False)

		if file.has_fit == True:

			fit_params["File"].append(file.name)
			fit_params["Lp-e"].append(file.fit_parameters["Lp_ext"][0])
			fit_params["Lc-e"].append(file.fit_parameters["Lc_ext"][0])
			fit_params["S-e"].append(file.fit_parameters["S_ext"][0])
			fit_params["F0-e"].append(file.fit_parameters["F0_ext"][0])
			fit_params["Lp-r"].append(file.fit_parameters["Lp_ret"][0])
			fit_params["Lc-r"].append(file.fit_parameters["Lc_ret"][0])
			fit_params["S-r"].append(file.fit_parameters["S_ret"][0])
			fit_params["F0-r"].append(file.fit_parameters["F0_ret"][0])
	out = pd.DataFrame(fit_params)
	out.to_csv(os.path.join(GLOBALVARS.output_directory,"FIT_PARAMETERS.csv"), index=False)

'''
 |------------------|
 |  GUI management  |
 |------------------|
'''



window = tk.Tk()
window.title("FD-curve suite")
window.columnconfigure(0, weight=1)
window.columnconfigure(1, weight=3)
window.rowconfigure(0, weight=0)
window.rowconfigure(1, weight=3)
window.rowconfigure(2, weight=1)


variable_radio_buttons = tk.StringVar()
variable_checkbutton_view = tk.BooleanVar()

# Create and set all top menubar options
menubar = Menu(window)
window.config(menu=menubar)

# Create all the menubar options
file_menu = Menu(menubar)
file_menu.add_command(label='Open Folder',command=open_folder)
file_menu.add_command(label='Export to CSV',command=export_data)
file_menu.add_command(label='Exit',command=window.destroy)


# Create the selection options
selection_menu = Menu(menubar)
selection_menu.add_command(label='Select Highlighted Curves',command=add_selected_curves)
selection_menu.add_command(label='Deselect Highlighted Curves',command=deselect_curves)

# Create the selection options
calibration_menu = Menu(menubar)
selection_menu.add_command(label='Mark Curve as Baseline',command=mark_baseline)
selection_menu.add_command(label='Mark Curve as Reference',command=mark_reference)
calibration_menu.add_command(label='Subtract Baseline',command=baseline_correction)
calibration_menu.add_command(label='Auto Force Scale',command=auto_force_scale)
calibration_menu.add_command(label='Calculate Supercoiling Density',command=baseline_correction)

# Add the dropdowns to the menubar
menubar.add_cascade(label="File",menu=file_menu)
menubar.add_cascade(label="Selection",menu=selection_menu)
menubar.add_cascade(label="Calibration",menu=calibration_menu)

frame_title_manager = tk.Frame(master=window)
frame_title_manager.grid(row=0, column=0, columnspan=2)
frame_file_managers = tk.Frame(master=window)
frame_file_managers.grid(row=1, rowspan=2, column=0, sticky=[tk.N,tk.E, tk.S, tk.W])
frame_file_managers.columnconfigure(0, weight=1)
frame_file_managers.columnconfigure(1, weight=0)
frame_file_managers.columnconfigure(2, weight=1)
frame_file_managers.columnconfigure(3, weight=0)
frame_file_managers.rowconfigure(0, weight=0)
frame_file_managers.rowconfigure(1, weight=1)
frame_file_managers.rowconfigure(2, weight=0)
frame_file_managers.rowconfigure(3, weight=0)
frame_file_managers.rowconfigure(4, weight=0)

frame_graphing_windows = tk.Frame(master=window)
frame_graphing_windows.grid(row=1, column=1, sticky=[tk.N, tk.E, tk.S, tk.W])
frame_input_buttons = tk.Frame(master=window)
frame_input_buttons.grid(row=2, column=1)
frame_graphing_windows.columnconfigure(0, weight=0)
frame_graphing_windows.columnconfigure(1, weight=1)
frame_graphing_windows.columnconfigure(2, weight=1)
frame_graphing_windows.columnconfigure(3, weight=1)
frame_graphing_windows.columnconfigure(4, weight=1)
frame_graphing_windows.columnconfigure(5, weight=1)
frame_graphing_windows.columnconfigure(6, weight=1)
frame_graphing_windows.columnconfigure(7, weight=1)
frame_graphing_windows.columnconfigure(8, weight=0)
frame_graphing_windows.rowconfigure(0, weight=1)
frame_graphing_windows.rowconfigure(1, weight=0)
frame_graphing_windows.rowconfigure(2, weight=0)
frame_graphing_windows.rowconfigure(3, weight=0)

title_label = tk.Label(master=frame_title_manager, text="DEMO")
title_label.pack()

frame_optical_settings = tk.Frame(master=frame_file_managers)
frame_optical_settings.grid(row=4, column=0, columnspan=3, sticky=[tk.E, tk.W])
frame_optical_settings.rowconfigure(0, weight=0)
frame_optical_settings.columnconfigure(0, weight=0)
frame_optical_settings.columnconfigure(1, weight=0)
frame_optical_settings.columnconfigure(2, weight=0)
frame_optical_settings.columnconfigure(3, weight=1)
frame_optical_settings.columnconfigure(4, weight=0)
frame_optical_settings.columnconfigure(5, weight=0)
frame_optical_settings.columnconfigure(6, weight=0)
frame_optical_settings.columnconfigure(7, weight=1)
frame_optical_settings.columnconfigure(8, weight=0)

# File selection listboxes

## Create listboxes
tk.Label(master=frame_file_managers, text="All H5 Files:").grid(row=0, column=0)
listbox_all_h5_files = tk.Listbox(master=frame_file_managers,selectmode=tk.MULTIPLE)
listbox_all_h5_files.bind('<<ListboxSelect>>', all_h5_listbox_select)
listbox_all_h5_files.grid(row=1, column=0, sticky=[tk.N,tk.E, tk.S, tk.W])
## associates scrollbars
scrollbar_v_all_h5 = Scrollbar(master=frame_file_managers, orient=VERTICAL, command=listbox_all_h5_files.yview)
scrollbar_v_all_h5.grid(row=1, column=1, sticky=[tk.N, tk.S])
scrollbar_h_all_h5 = Scrollbar(master=frame_file_managers, orient=HORIZONTAL, command=listbox_all_h5_files.xview)
scrollbar_h_all_h5.grid(row=2, column=0, sticky=[tk.E, tk.W])
## Second listbox
tk.Label(master=frame_file_managers, text="Selected Files:").grid(row=0, column=2)
listbox_all_selected_files = tk.Listbox(master=frame_file_managers)
listbox_all_selected_files.bind('<<ListboxSelect>>', all_selected_listbox_select)
listbox_all_selected_files.grid(row=1, column=2, sticky=[tk.N,tk.E, tk.S, tk.W])
## associated scrollbars
scrollbar_v_selected_h5 = Scrollbar(master=frame_file_managers, orient=VERTICAL, command=listbox_all_selected_files.yview)
scrollbar_v_selected_h5.grid(row=1, column=3, sticky=[tk.N, tk.S])
scrollbar_h_selected_h5 = Scrollbar(master=frame_file_managers, orient=HORIZONTAL, command=listbox_all_selected_files.xview)
scrollbar_h_selected_h5.grid(row=2, column=2, sticky=[tk.E, tk.W])
## Attach scrollbars to listboxes
listbox_all_h5_files.config(yscrollcommand=scrollbar_v_all_h5.set,xscrollcommand=scrollbar_h_all_h5.set)
listbox_all_selected_files.config(yscrollcommand=scrollbar_v_selected_h5.set,xscrollcommand=scrollbar_h_selected_h5.set)

# Set key optical settings used in the session
tk.Label(master=frame_optical_settings, text="FPS: ").grid(row=0, column=0, sticky=[tk.E])
entry_frame_rate = tk.Entry(master=frame_optical_settings, width=3)
entry_frame_rate.insert(0, str(GLOBALVARS.frame_rate))
entry_frame_rate.grid(row=0, column=1, sticky=[tk.W])
tk.Label(master=frame_optical_settings, text="Hz").grid(row=0, column=2, sticky=tk.W)
tk.Label(master=frame_optical_settings, text="Extension Speed: ").grid(row=0, column=4, sticky=[tk.E])
entry_extension_speed = tk.Entry(master=frame_optical_settings, width=4)
entry_extension_speed.insert(0, str(GLOBALVARS.extension_speed_um_s))
entry_extension_speed.grid(row=0, column=5, sticky=tk.W)
tk.Label(master=frame_optical_settings, text="um/s").grid(row=0, column=6, sticky=tk.W)
tk.Button(master=frame_optical_settings, text="Update", command=update_optic_settings).grid(row=0, column=8)


# Canvas to display graphs
canvas_graph_display = tk.Canvas(master=frame_graphing_windows, bg="#856ff8")
canvas_graph_display.grid(row=0, column=0, columnspan=9, sticky=[tk.N, tk.E, tk.S, tk.W])

# FD-time cutoff scrollers
tk.Label(master=frame_graphing_windows, text="Xmin: ").grid(row=1, column=0, sticky=tk.S)
scale_select_min_time = Scale(master=frame_graphing_windows, orient=HORIZONTAL, resolution=0.01)
scale_select_min_time.grid(row=1, column=1, columnspan=8, sticky=[tk.E, tk.W])
tk.Label(master=frame_graphing_windows, text="Xmax: ").grid(row=2, column=0, sticky=tk.S)
scale_select_min_time.bind("<ButtonRelease-1>", slider_min_release)
scale_select_max_time = Scale(master=frame_graphing_windows, orient=HORIZONTAL, resolution=0.01)
scale_select_max_time.grid(row=2, column=1, columnspan=8,sticky=[tk.E, tk.W])
scale_select_max_time.bind("<ButtonRelease-1>", slider_max_release)

tk.Label(master=frame_graphing_windows, text="Xmin:").grid(row=3, column=0, sticky=tk.E)
entry_xmin = tk.Entry(master=frame_graphing_windows,width=6)
entry_xmin.insert(0,"0")
entry_xmin.grid(row=3, column=1, sticky=tk.W)
tk.Label(master=frame_graphing_windows, text="Xmax:").grid(row=3, column=2, sticky=tk.E)
entry_xmax = tk.Entry(master=frame_graphing_windows,width=6)
entry_xmax.insert(0,"0")
entry_xmax.grid(row=3,column=3, sticky=tk.W)


tk.Label(master=frame_graphing_windows, text="Ymin:").grid(row=3, column=4, sticky=tk.E)
entry_ymin = tk.Entry(master=frame_graphing_windows,width=6)
entry_ymin.insert(0,"-5")
entry_ymin.grid(row=3, column=5, sticky=tk.W)
tk.Label(master=frame_graphing_windows, text="Ymax:").grid(row=3, column=6, sticky=tk.E)
entry_ymax = tk.Entry(master=frame_graphing_windows,width=6)
entry_ymax.insert(0,"30")
entry_ymax.grid(row=3, column=7, sticky=tk.W)

tk.Button(master=frame_graphing_windows, text="Update", command=update_trim_settings).grid(row=3, column=8)

# Buttons to set correction factors
button_1 = Button(master=frame_input_buttons, text="Toggle Time", command=toggle_time)
button_1.grid(row=0,column=0)
button_2 = Button(master=frame_input_buttons, text="Auto Trim", command=auto_trim_data)
button_2.grid(row=0,column=1)
button_3 = Button(master=frame_input_buttons, text="Manual Trim", command=manual_trim_data)
button_3.grid(row=0,column=2)
tk.Checkbutton(master=frame_input_buttons, text="Display Full", variable=variable_checkbutton_view, command=radio_button_select).grid(row=1, column=0)
tk.Radiobutton(master=frame_input_buttons, text="Set Extension", value = "extension", variable=variable_radio_buttons, command=radio_button_select).grid(row=1, column=1)
tk.Radiobutton(master=frame_input_buttons, text="Set Retraction", value="retraction", variable=variable_radio_buttons, command=radio_button_select).grid(row=1, column = 2)

Enlarge_Button = Button(master=frame_input_buttons, text="Expand Graph", command=expand_graph)
Enlarge_Button.grid(row=0, column=3)


tk.Label(master=frame_input_buttons, text="Lp:").grid(row=2, column=0, sticky=tk.E)
Lp_entry = tk.Entry(master=frame_input_buttons,width=6)
Lp_entry.insert(0,"50")
Lp_entry.grid(row=2, column=1)
Lp_display = tk.Label(master=frame_input_buttons, text="")
Lp_display.grid(row=2, column=2, sticky=tk.EW)
tk.Label(master=frame_input_buttons, text="Lc:").grid(row=3, column=0, sticky=tk.E)
Lc_entry = tk.Entry(master=frame_input_buttons,width=6)
Lc_entry.insert(0,"16.5")
Lc_entry.grid(row=3, column=1)
Lc_display = tk.Label(master=frame_input_buttons, text="")
Lc_display.grid(row=3, column=2, sticky=tk.EW)
tk.Label(master=frame_input_buttons, text="S:").grid(row=4, column=0, sticky=tk.E)
S_entry = tk.Entry(master=frame_input_buttons,width=6)
S_entry.insert(0,"1500")
S_entry.grid(row=4, column=1)
S_display = tk.Label(master=frame_input_buttons, text="")
S_display.grid(row=4, column=2, sticky=tk.EW)
tk.Label(master=frame_input_buttons, text="F0:").grid(row=5, column=0, sticky=tk.E)
F0_entry = tk.Entry(master=frame_input_buttons,width=6)
F0_entry.insert(0,"0")
F0_entry.grid(row=5, column=1)
F0_display = tk.Label(master=frame_input_buttons, text="")
F0_display.grid(row=5, column=2, sticky=tk.EW)

Fit_Button = Button(master=frame_input_buttons, text="Fit", command=fit)
Fit_Button.grid(row=8, column=1)



canvas_graph_display.bind("<Configure>", window_resize)
window.mainloop()
