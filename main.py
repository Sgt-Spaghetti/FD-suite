'''
  |----------------------------------|
  |  Written by Leonardo Cherin, UCL |
  |    Distributed as free software  | 
  |       with the GPL3 license      |
  |----------------------------------|

No AI was used at any point in development
'''

import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
import h5py
import os

# A few useful global variables for the lifetime of the program
class GLOBALVARS():
	def __init__(self) -> None:
		self.all_files: list = [] # Every h5 in the folder
		self.files_awaiting_selection: list = [] # files multi-selected
		self.selected_files: list = [] # files confirmed to process
		self.active_file = None # The currently focused file
		self.output_directory: str = ""
		self.graph_image = None
		self.extension_speed_um_s: float = 0.15
		self.frame_rate: int = 100
		self.baseline_curve = pd.DataFrame({"Force": [], "Distance": []})
		self.show_first_deriv = True
		self.show_second_deriv = False

GLOBALVARS = GLOBALVARS()

# Create the "FD" class, to store all opened FD curves
class FD():
	def __init__(self, filepath: str) -> None:
		self.filepath: str = filepath
		self.name: str = os.path.splitext(os.path.basename(filepath))[0]
		self.plot_time: bool = False

		# Bounding boxes for discrete extension curve
		#                                   D,T,OD,OT
		self.xmin_e: list[float] = [0,0, 0, 0]
		self.xmax_e: list[float] = [0,0, 0, 0]
		self.ymin_e: float = -5
		self.ymax_e: float = 125
		# Bounding boxes for discrete retraction curve
		self.xmin_r: list[float] = [0,0,0,0]
		self.xmax_r: list[float] = [0,0,0,0]
		self.ymin_r: float = -5
		self.ymax_r: float = 125
		# Bounding boxes for sub-area of extension curve to fit to
		self.xmin_f_e: list[float, float] = [0,0,0,0]
		self.xmax_f_e: list[float, float] = [0,0,0,0]
		self.ymin_f_e: float = -5
		self.ymax_f_e: float = 30
		# Bounding boxes for sub-area of retraction curve to fit to
		self.xmin_f_r: list[float, float] = [0,0,0,0]
		self.xmax_f_r: list[float, float] = [0,0,0,0]
		self.ymin_f_r: float = -5
		self.ymax_f_r: float = 30

		self.sigma_e: float = -100
		self.sigma_r: float = 0
		self.trimmed_e = False
		self.trimmed_r = False
		self.trimmed_f_e = False
		self.trimmed_f_r = False
		self.has_fit_e: bool = False
		self.has_fit_r: bool = False
		self.baseline: bool = False
		self.reference: bool = False
		self.is_baseline_subtracted: bool = False
		self.is_force_scaled: bool = False

		# Trimmed area of the curve that is one discrete extension curve
		self.dataframe_extension = pd.DataFrame({"Force_Extension": [], "Distance_Extension": [], "Time_Extension": []})
		# Trimmed area of the curve that is one discrete retraction curve
		self.dataframe_retraction = pd.DataFrame({"Force_Retraction": [], "Distance_Retraction": [], "Time_Retraction": []})
		# sub-area of extension curve used for fitting
		self.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": [], "Fit_Distance_Extension": [], "Fit_Time_Extension": []})
		# sub-area of retraction curve used for fitting
		self.fit_dataframe_retraction = pd.DataFrame({"Fit_Force_Retraction": [], "Fit_Distance_Retraction": [], "Fit_Time_Retraction": []})
	
		# Once fit, evaluate the WLC model and "precalculate" a linegraph for the fit, convienient for plotting with external programs
		self.precalculated_fit_extension = pd.DataFrame({"Precalculated_Fit_Force_Extension": [], "Precalculated_Fit_Distance_Extension": [], "Precalculated_Fit_Time_Extension": []})
		self.precalculated_fit_retraction = pd.DataFrame({"Precalculated_Fit_Force_Retraction": [], "Precalculated_Fit_Distance_Retraction": [], "Precalculated_Fit_Time_Retraction": []})

		self.fit_parameters = pd.DataFrame({"Lp_ext": [], "Lc_ext": [], "S_ext": [], "F0_ext": [],"Lp_ret": [], "Lc_ret": [], "S_ret": [], "F0_ret": []})

		# Critical force
		self.fc_e: float = 0
		self.fc_r: float = 0

		self.initialise_attributes()

	def initialise_attributes(self) -> None:
		# NOTE: files are kept in RAM. If opening thousands this
		# might cause an issue, but it is highly unlikely
		raw_data = h5py.File(self.filepath, 'r')
		force_data: list[float] = raw_data["Force LF"][y_variable_combo.get()]["Value"]
		distance_data: list[float] = raw_data["Distance"][x_variable_combo.get()]["Value"]
		#distance_to_time_conversion: list[float] = [(i/GLOBALVARS.frame_rate) for i in range(len(distance_data))]
		time_data: list[float] = raw_data["Force LF"][y_variable_combo.get()]["Timestamp"]
		zero_time = time_data[0]
		distance_to_time_conversion: list[float] = (time_data - zero_time)/1000000000

		# Compute derivatives of the raw data, useful for data trimming
		derivatives: list = self.differentiate_savgol(distance_to_time_conversion, force_data, 0.75*GLOBALVARS.frame_rate/GLOBALVARS.extension_speed_um_s, 2)
		first_derivative = [derivatives[0], derivatives[1]]
		second_derivative = [derivatives[0], derivatives[2]]

		# Initialise the starting core dataset, and a dummy "processed" dataset onto which any corrections (eg baseline subtraction) will be applied
		self.dataframe = pd.DataFrame({"Force": force_data, "Distance": distance_data, "Time": distance_to_time_conversion})
		self.processed_dataframe = pd.DataFrame({"Processed_Force": force_data, "Processed_Distance": distance_data, "Processed_Time":  distance_to_time_conversion})
		self.first_derivative_dataframe= pd.DataFrame({"First_Derivative": first_derivative[1], "Time": first_derivative[0], "Distance": distance_data})
		self.second_derivative_dataframe= pd.DataFrame({"Second_Derivative": second_derivative[1], "Time": second_derivative[0], "Distance": distance_data})

		self.trimmed_first_derivative_dataframe_extension = pd.DataFrame({"Trimmed_First_Derivative": [], "Trimmed_Time": [], "Trimmed_Distance": []})
		self.trimmed_first_derivative_dataframe_retraction = pd.DataFrame({"Trimmed_First_Derivative": [], "Trimmed_Time": [], "Trimmed_Distance": []})

	# The central function responsible for plotting the curves, depending on what "state" the curve is in
	# For example, if it has been trimmed into extension / retraction curves, or fit with the eOdjik model
	def plot(self, expand_graph = False) -> None:

		selected = False	
		for f in GLOBALVARS.selected_files:
			if self.name == f.name:
				selected = True

		if self.trimmed_e == False and self.trimmed_r == False: # There is no trimmed data to plot
			if selected == True: # If the curve is selected
				self.plot_scatter_processed_data()
			else: # It is just a "all files" preview
				self.plot_scatter_raw_data()
		else: # we have a trim to plot
			if variable_radio_buttons_view.get() == "full": # We want to see the full curve regardless
				self.plot_scatter_processed_data()
			else: # We are looking for a particular trim!
				if variable_checkbutton_view_fit.get() == True: # We actually want to see the fit
					self.plot_fit_data()
				else:
					self.plot_trimmed_data()

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
			plt.xlabel("Distance (\u03bcm)")
		plt.ylabel("Force (pN)")
		plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_RAW_SCATTER.png"))
		plt.savefig("TEMP_PLOT.png")


	# Plot the data which might have been baseline subtracted or pixel corrected
	# Therefore, use the "processed_dataframe" dataset.
	def plot_scatter_processed_data(self) -> None:
		window.update()
		window.update_idletasks()
		width: int = canvas_graph_display.winfo_width()
		height: int = canvas_graph_display.winfo_height()
	
		if self.plot_time == True:
			if GLOBALVARS.show_first_deriv == True and GLOBALVARS.show_second_deriv == True:
				fig, ax = plt.subplots(3,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax[0].scatter(self.processed_dataframe["Processed_Time"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax[1].plot(self.first_derivative_dataframe["Time"], self.first_derivative_dataframe["First_Derivative"])
				ax[2].plot(self.second_derivative_dataframe["Time"], self.second_derivative_dataframe["Second_Derivative"])
				ax[2].set_xlabel("Time (s)")
				if self.xmin_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_e[1], color="r")
					ax[1].axvline(self.xmin_e[1], color="r")
					ax[2].axvline(self.xmin_e[1], color="r")
				if self.xmax_e[1]> int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_e[1], color="g")
					ax[1].axvline(self.xmax_e[1], color="g")
					ax[2].axvline(self.xmax_e[1], color="g")
				if self.ymax_e > 0:
					ax[0].axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax[0].axhline(self.ymin_e, c="blue")
				if self.xmin_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_r[1], color="r", linestyle="--")
					ax[1].axvline(self.xmin_r[1], color="r", linestyle="--")
					ax[2].axvline(self.xmin_r[1], color="r", linestyle="--")
				if self.xmax_r[1]> int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_r[1], color="g", linestyle="--")
					ax[1].axvline(self.xmax_r[1], color="g", linestyle="--")
					ax[2].axvline(self.xmax_r[1], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax[0].axhline(self.ymax_r, c="black", linestyle="--")
				if self.ymin_r > 0:
					ax[0].axhline(self.ymin_r, c="blue", linestyle="--")

				ax[0].set_ylabel("Force (pN)")
				ax[1].axhline(0, c="black")
				ax[1].set_ylabel("dy/dx (pN/sm)")
				ax[2].axhline(0, c="black")
				ax[2].set_ylabel("ddy/dx (pN/s$^2$)")
			elif GLOBALVARS.show_first_deriv == True and GLOBALVARS.show_second_deriv == False: 
				fig, ax = plt.subplots(2,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax[0].scatter(self.processed_dataframe["Processed_Time"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax[1].plot(self.first_derivative_dataframe["Time"], self.first_derivative_dataframe["First_Derivative"])
				ax[1].set_xlabel("Time (s)")
				if self.xmin_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_e[1], color="r")
					ax[1].axvline(self.xmin_e[1], color="r")
				if self.xmax_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_e[1], color="g")
					ax[1].axvline(self.xmax_e[1], color="g")
				if self.ymax_e > 0:
					ax[0].axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax[0].axhline(self.ymin_e_e, c="blue")
				if self.xmin_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_r[1], color="r", linestyle="--")
					ax[1].axvline(self.xmin_r[1], color="r", linestyle="--")
				if self.xmax_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_r[1], color="g", linestyle="--")
					ax[1].axvline(self.xmax_r[1], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax[0].axhline(self.ymax_r, c="black", linestyle="--")
				if self.ymin_r > 0:
					ax[0].axhline(self.ymin_r, c="blue", linestyle="--")

				ax[0].set_ylabel("Force (pN)")
				ax[1].axhline(0, c="black")
				ax[1].set_ylabel("dy/dx (pN/s)")
			elif GLOBALVARS.show_first_deriv == False and GLOBALVARS.show_second_deriv == True: 
				fig, ax = plt.subplots(2,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax[0].scatter(self.processed_dataframe["Processed_Time"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax[1].plot(self.second_derivative_dataframe["Time"], self.second_derivative_dataframe["Second_Derivative"])
				ax[1].set_xlabel("Time (s)")
				if self.xmin_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_e[1], color="r")
					ax[1].axvline(self.xmin_e[1], color="r")
				if self.xmax_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_e[1], color="g")
					ax[1].axvline(self.xmax_e[1], color="g")
				if self.ymax_e > 0:
					ax[0].axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax[0].axhline(self.ymin_e, c="blue")
				if self.xmin_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_r[1], color="r", linestyle="--")
					ax[1].axvline(self.xmin_r[1], color="r", linestyle="--")
				if self.xmax_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_r[1], color="g", linestyle="--")
					ax[1].axvline(self.xmax_r[1], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax[0].axhline(self.ymax_r[1], c="black", linestyle="--")
				if self.ymin_r > 0:
					ax[0].axhline(self.ymin_r, c="blue", linestyle="--")

				ax[0].set_ylabel("Force (pN)")
				ax[1].axhline(0, c="black")
				ax[1].set_ylabel("ddy/ddx (pN/s$^2$)")
			else:
				fig, ax = plt.subplots(1,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax.scatter(self.processed_dataframe["Processed_Time"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax.set_xlabel("Time (s)")
				if self.xmin_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax.axvline(self.xmin_e[1], color="r")
				if self.xmax_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax.axvline(self.xmax_e[1], color="g")
				if self.ymax_e > 0:
					ax.axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax.axhline(self.ymin_e, c="blue")
				if self.xmin_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax.axvline(self.xmin_r[1], color="r", linestyle="--")
				if self.xmax_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
					ax.axvline(self.xmax_r[1], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax.axhline(self.ymax_r, c="black", linestyle="--")
				if self.ymin_r > 0:
					ax.axhline(self.ymin_r, c="blue", linestyle="--")

				ax.set_ylabel("Force (pN)")
		else:
			if GLOBALVARS.show_first_deriv == True and GLOBALVARS.show_second_deriv == True:
				fig, ax = plt.subplots(3,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax[0].scatter(self.processed_dataframe["Processed_Distance"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax[1].plot(self.first_derivative_dataframe["Distance"], self.first_derivative_dataframe["First_Derivative"])
				ax[2].plot(self.second_derivative_dataframe["Distance"], self.second_derivative_dataframe["Second_Derivative"])
				ax[2].set_xlabel("Distance (\u03bcm)")
				if self.xmin_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_e[0], color="r")
					ax[1].axvline(self.xmin_e[0], color="r")
					ax[2].axvline(self.xmin_e[0], color="r")
				if self.xmax_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_e[0], color="g")
					ax[1].axvline(self.xmax_e[0], color="g")
					ax[2].axvline(self.xmax_e[0], color="g")
				if self.ymax_e > 0:
					ax[0].axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax[0].axhline(self.ymin_e, c="blue")
				if self.xmin_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_r[0], color="r", linestyle="--")
					ax[1].axvline(self.xmin_r[0], color="r", linestyle="--")
					ax[2].axvline(self.xmin_r[0], color="r", linestyle="--")
				if self.xmax_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_r[0], color="g", linestyle="--")
					ax[1].axvline(self.xmax_r[0], color="g", linestyle="--")
					ax[2].axvline(self.xmax_r[0], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax[0].axhline(self.ymax_r, c="black", linestyle="--")
				if self.ymin_r > 0:
					ax[0].axhline(self.ymin_r, c="blue", linestyle="--")

				ax[0].set_ylabel("Force (pN)")
				ax[1].axhline(0, c="black")
				ax[1].set_ylabel("dy/dx (pN/\u03bcm)")
				ax[2].axhline(0, c="black")
				ax[2].set_ylabel("ddy/dx (pN/\u03bcm$^2$)")
			elif GLOBALVARS.show_first_deriv == True and GLOBALVARS.show_second_deriv == False: 
				fig, ax = plt.subplots(2,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax[0].scatter(self.processed_dataframe["Processed_Distance"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax[1].plot(self.first_derivative_dataframe["Distance"], self.first_derivative_dataframe["First_Derivative"])
				ax[1].set_xlabel("Distance (\u03bcm)")
				if self.xmin_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_e[0], color="r")
					ax[1].axvline(self.xmin_e[0], color="r")
				if self.xmax_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_e[0], color="g")
					ax[1].axvline(self.xmax_e[0], color="g")
				if self.ymax_e > 0:
					ax[0].axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax[0].axhline(self.ymin_e, c="blue")
				if self.xmin_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_r[0], color="r", linestyle="--")
					ax[1].axvline(self.xmin_r[0], color="r", linestyle="--")
				if self.xmax_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_r[0], color="g", linestyle="--")
					ax[1].axvline(self.xmax_r[0], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax[0].axhline(self.ymax_r, c="black", linestyle="--")
				if self.ymin_r > 0:
					ax[0].axhline(self.ymin_r, c="blue", linestyle="--")

				ax[0].set_ylabel("Force (pN)")
				ax[1].axhline(0, c="black")
				ax[1].set_ylabel("dy/dx (pN/\u03bcm)")
			elif GLOBALVARS.show_first_deriv == False and GLOBALVARS.show_second_deriv == True: 
				fig, ax = plt.subplots(2,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax[0].scatter(self.processed_dataframe["Processed_Distance"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax[1].plot(self.second_derivative_dataframe["Distance"], self.second_derivative_dataframe["Second_Derivative"])
				ax[1].set_xlabel("Distance (\u03bcm)")
				if self.xmin_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_e[0], color="r")
					ax[1].axvline(self.xmin_e[0], color="r")
				if self.xmax_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_e[0], color="g")
					ax[1].axvline(self.xmax_e[0], color="g")
				if self.ymax_e > 0:
					ax[0].axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax[0].axhline(self.ymin_e, c="blue")
				if self.xmin_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmin_r[0], color="r", linestyle="--")
					ax[1].axvline(self.xmin_r[0], color="r", linestyle="--")
				if self.xmax_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax[0].axvline(self.xmax_r[0], color="g", linestyle="--")
					ax[1].axvline(self.xmax_r[0], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax[0].axhline(self.ymax_r, c="black", linestyle="--")
				if self.ymin_r > 0:
					ax[0].axhline(self.ymin_r, c="blue", linestyle="--")

				ax[0].set_ylabel("Force (pN)")
				ax[1].axhline(0, c="black")
				ax[1].set_ylabel("ddy/ddx (pN/\u03bcm$^2$)")
			else:
				fig, ax = plt.subplots(1,1,figsize=(width/100, height/100))
				fig.suptitle(self.name)
				ax.scatter(self.processed_dataframe["Processed_Distance"], self.processed_dataframe["Processed_Force"], s=0.1)
				ax.set_xlabel("Distance (\u03bcm)")
				if self.xmin_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax.axvline(self.xmin_e[0], color="r")
				if self.xmax_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax.axvline(self.xmax_e[0], color="g")
				if self.ymax_e > 0:
					ax.axhline(self.ymax_e, c="black")
				if self.ymin_e > 0:
					ax.axhline(self.ymin_e, c="blue")
				if self.xmin_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax.axvline(self.xmin_r[0], color="r", linestyle="--")
				if self.xmax_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
					ax.axvline(self.xmax_r[0], color="g", linestyle="--")
				if self.ymax_r > 0:
					ax.axhline(self.ymax_r, c="black", linestyle="--")
				if self.ymin_r > 0:
					ax.axhline(self.ymin_r, c="blue", linestyle="--")

				ax.set_ylabel("Force (pN)")

		plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_PROCESSED_SCATTER.png"))
		plt.savefig("TEMP_PLOT.png")


	# Plot data which has been trimmed into extension or retraction curves
	def plot_trimmed_data(self) -> None:
		window.update()
		window.update_idletasks()
		width: int = canvas_graph_display.winfo_width()
		height: int = canvas_graph_display.winfo_height()
		if GLOBALVARS.show_first_deriv == False:
			fig, ax = plt.subplots(2,1,figsize=(width/100, height/100))
			fig.suptitle(self.name)
			if variable_radio_buttons_view.get() == "extension":
				if self.plot_time == True:
					ax.scatter(self.dataframe_extension["Time_Extension"], self.dataframe_extension["Force_Extension"], s=0.1)
					if self.xmin_f_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax.axvline(self.xmin_f_e[1], color="r")
					if self.xmax_f_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax.axvline(self.xmax_f_e[1], color="g")
					if self.ymax_f_e > 0:
						ax.axhline(self.ymax_f_e, c="black")
					if self.ymin_f_e > 0:
						ax.axhline(self.ymin_f_e, c="blue")
					ax.set_xlabel("Time (s)")
				else:
					plt.scatter(self.dataframe_extension["Distance_Extension"], self.dataframe_extension["Force_Extension"], s=0.1)
					if self.xmin_f_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax.axvline(self.xmin_f_e[0], color="r")
					if self.xmax_f_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax.axvline(self.xmax_f_e[0], color="g")
					if self.ymax_f_e > 0:
						ax.axhline(self.ymax_f_e, c="black")
					if self.ymin_f_e > 0:
						ax.axhline(self.ymin_f_e, c="blue")
					ax.set_xlabel("Distance (\u03bcm)")

			elif variable_radio_buttons_view.get() == "retraction":
				if self.plot_time == True:
					ax.scatter(self.dataframe_retraction["Time_Retraction"], self.dataframe_retraction["Force_Retraction"], s=0.1)
					if self.xmin_f_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax.axvline(self.xmin_f_r[1], color="r")
					if self.xmax_f_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax.axvline(self.xmax_f_r[1], color="g")
					if self.ymax_f_r > 0:
						ax.axhline(self.ymax_f_r, c="black")
					if self.ymin_f_r > 0:
						ax.axhline(self.ymin_f_r, c="blue")
					ax.set_xlabel("Time (s)")
				else:
					ax.scatter(self.dataframe_retraction["Distance_Retraction"], self.dataframe_retraction["Force_Retraction"], s=0.1)
					if self.xmin_f_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax.axvline(self.xmin_f_r[0], color="r")
					if self.xmax_f_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax.axvline(self.xmax_f_r[0], color="g")
					if self.ymax_f_r > 0:
						ax.axhline(self.ymax_f_r, c="black")
					if self.ymin_f_r > 0:
						ax.axhline(self.ymin_f_r, c="blue")
					ax.set_xlabel("Distance (\u03bcm)")
			ax.set_ylabel("Force (pN)")

		else:
			fig, ax = plt.subplots(2,1,figsize=(width/100, height/100))
			fig.suptitle(self.name)
			if variable_radio_buttons_view.get() == "extension":
				if self.plot_time == True:
					ax[0].scatter(self.dataframe_extension["Time_Extension"], self.dataframe_extension["Force_Extension"], s=0.1)
					ax[1].plot(self.trimmed_first_derivative_dataframe_extension["Trimmed_Time"], self.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"])
					if self.xmin_f_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax[0].axvline(self.xmin_f_e[1], color="r")
						ax[1].axvline(self.xmin_f_e[1], color="r")
					if self.xmax_f_e[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax[0].axvline(self.xmax_f_e[1], color="g")
						ax[1].axvline(self.xmax_f_e[1], color="g")
					if self.ymax_f_e > 0:
						ax[0].axhline(self.ymax_f_e, c="black")
					if self.ymin_f_e > 0:
						ax[0].axhline(self.ymin_f_e, c="blue")
					ax[1].set_xlabel("Time (s)")
					ax[1].set_ylabel("dy/dx (pN/s)")
					ax[0].set_ylabel("Force (pN)")
				else:
					ax[0].scatter(self.dataframe_extension["Distance_Extension"], self.dataframe_extension["Force_Extension"], s=0.1)
					ax[1].plot(self.trimmed_first_derivative_dataframe_extension["Trimmed_Distance"], self.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"])
					if self.xmin_f_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax[0].axvline(self.xmin_f_e[0], color="r")
						ax[1].axvline(self.xmin_f_e[0], color="r")
					if self.xmax_f_e[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax[0].axvline(self.xmax_f_e[0], color="g")
						ax[1].axvline(self.xmax_f_e[0], color="g")
					if self.ymax_f_e > 0:
						ax[0].axhline(self.ymax_f_e, c="black")
					if self.ymin_f_e > 0:
						ax[0].axhline(self.ymin_f_e, c="blue")
					ax[1].set_xlabel("Distance (\u03bcm)")
					ax[1].set_ylabel("dy/dx (pN/\u03bcm)")
					ax[0].set_ylabel("Force (pN)")

			elif variable_radio_buttons_view.get() == "retraction":
				if self.plot_time == True:
					ax[0].scatter(self.dataframe_retraction["Time_Retraction"], self.dataframe_retraction["Force_Retraction"], s=0.1)
					ax[1].plot(self.trimmed_first_derivative_dataframe_retraction["Trimmed_Time"], self.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"])
					if self.xmin_f_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax[0].axvline(self.xmin_f_r[1], color="r")
						ax[1].axvline(self.xmin_f_r[1], color="r")
					if self.xmax_f_r[1] > int((min(self.dataframe["Time"])*100)+0.5)/100:
						ax[0].axvline(self.xmax_f_r[1], color="g")
						ax[1].axvline(self.xmax_f_r[1], color="g")
					if self.ymax_f_r > 0:
						ax[0].axhline(self.ymax_f_r, c="black")
					if self.ymin_f_r > 0:
						ax[0].axhline(self.ymin_f_r, c="blue")
					ax[1].set_xlabel("Time (s)")
					ax[1].set_ylabel("dy/dx (pN/s)")
					ax[0].set_ylabel("Force (pN)")
				else:
					ax[0].scatter(self.dataframe_retraction["Distance_Retraction"], self.dataframe_retraction["Force_Retraction"], s=0.1)
					ax[1].plot(self.trimmed_first_derivative_dataframe_retraction["Trimmed_Distance"], self.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"])
					if self.xmin_f_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax[0].axvline(self.xmin_f_r[0], color="r")
						ax[1].axvline(self.xmin_f_r[0], color="r")
					if self.xmax_f_r[0] > int((min(self.dataframe["Distance"])*100)+0.5)/100:
						ax[0].axvline(self.xmax_f_r[0], color="g")
						ax[1].axvline(self.xmax_f_r[0], color="g")
					if self.ymax_f_r > 0:
						ax[0].axhline(self.ymax_f_r, c="black")
					if self.ymin_f_r > 0:
						ax[0].axhline(self.ymin_f_r, c="blue")
					ax[1].set_xlabel("Distance (\u03bcm)")
					ax[1].set_ylabel("dy/dx (pN/\u03bcm)")
					ax[0].set_ylabel("Force (pN)")

		plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_TRIMMED_SCATTER.png"))
		plt.savefig("TEMP_PLOT.png")
			
	def plot_fit_data(self) -> None:
		window.update()
		window.update_idletasks()
		width: int = canvas_graph_display.winfo_width()
		height: int = canvas_graph_display.winfo_height()
		plt.figure(figsize=(width/100, height/100))
		plt.title(self.name)
		if variable_checkbutton_view_fit.get() == True and variable_radio_buttons_view.get() == "extension":
			if self.plot_time == False:
				plt.scatter(self.fit_dataframe_extension["Fit_Distance_Extension"], self.fit_dataframe_extension["Fit_Force_Extension"],s=0.1)
				if self.has_fit_e == True:
					plt.plot(self.precalculated_fit_extension["Precalculated_Fit_Distance_Extension"], self.precalculated_fit_extension["Precalculated_Fit_Force_Extension"],c="r")
				plt.xlabel("Distance (\u03bcm)")
			else:
				plt.scatter(self.fit_dataframe_extension["Fit_Time_Extension"], self.fit_dataframe_extension["Fit_Force_Extension"],s=0.1)
				if self.has_fit_e == True:
					plt.plot(self.precalculated_fit_extension["Precalculated_Fit_Time_Extension"], self.precalculated_fit_extension["Precalculated_Fit_Force_Extension"],c="r")
				plt.xlabel("Time(s)")
		elif variable_checkbutton_view_fit.get() == True and variable_radio_buttons_view.get() == "retraction":
			if self.plot_time == False:
				plt.scatter(self.fit_dataframe_retraction["Fit_Distance_Retraction"], self.fit_dataframe_retraction["Fit_Force_Retraction"],s=0.1)
				if self.has_fit_r == True:
					plt.plot(self.precalculated_fit_retraction["Precalculated_Fit_Distance_Retraction"], self.precalculated_fit_retraction["Precalculated_Fit_Force_Retraction"],c="r")
				plt.xlabel("Distance (\u03bcm)")
			else:
				plt.scatter(self.fit_dataframe_retraction["Fit_Time_Retraction"], self.fit_dataframe_retraction["Fit_Force_Retraction"],s=0.1)
				if self.has_fit_r == True:
					plt.plot(self.precalculated_fit_retraction["Precalculated_Fit_Time_Retraction"], self.precalculated_fit_retraction["Precalculated_Fit_Force_Retraction"],c="r")
				plt.xlabel("Time(s)")

		plt.ylabel("Force (pN)")
		plt.savefig(os.path.join(GLOBALVARS.output_directory, self.name + "_FIT.png"))
		plt.savefig("TEMP_PLOT.png")

	# Data is differentiated by applying a Savitzky-Golay filter to the raw noisy data,
	# And extracting the first and second derivatives directly from the differentiation of
	# the polynomial coefficients returned buy the Savitzky-Golay fit.
	# This is performed in the time domain, which guarantees evenly spaced data by
	# 1/framerate for the optical trap's camera system.
	# We will use a window size of 0.5*framerate to create a pseudo-0.5Hz lowpass filter, and we
	# will fit the window to a second degree polynomial.

	def differentiate_savgol(self, xdata, ydata, windowsize=20, degree=2):
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
			if self.trimmed_e == True:
				# Subtract from the extension curve if applicable
				corrected_force = apply_baseline(self.dataframe_extension["Distance_Extension"], self.dataframe_extension["Force_Extension"])
				self.dataframe_extension = pd.DataFrame({"Force_Extension": corrected_force, "Distance_Extension": self.dataframe_extension["Distance_Extension"], "Time_Extension": self.dataframe_extension["Time_Extension"]})
				# Subtract from the retraction curve if applicable
				corrected_force = apply_baseline(self.dataframe_retraction["Distance_Retraction"], self.dataframe_retraction["Force_Retraction"])
				self.dataframe_extension = pd.DataFrame({"Force_Retraction": corrected_force, "Distance_Retraction": self.dataframe_retraction["Distance_Retraction"], "Time_Retraction": self.dataframe_retraction["Time_Retraction"]})

			
			# Compute derivatives of the raw data, useful for data trimming
			derivatives: list = self.differentiate_savgol(self.processed_dataframe["Processed_Time"], self.processed_dataframe["Processed_Force"], 0.75*GLOBALVARS.frame_rate/GLOBALVARS.extension_speed_um_s, 2)
			first_derivative = [derivatives[0], derivatives[1]]
			second_derivative = [derivatives[0], derivatives[2]]

			self.first_derivative_dataframe= pd.DataFrame({"First_Derivative": first_derivative[1], "Time": first_derivative[0], "Distance": self.processed_dataframe["Processed_Distance"]})
			self.second_derivative_dataframe= pd.DataFrame({"Second_Derivative": second_derivative[1], "Time": second_derivative[0], "Distance": self.processed_dataframe["Processed_Distance"]})

		replot_canvas()

	def reset(self) -> None:
		self.xmin: float = 0
		self.xmax: float = 0
		self.ymin: float = 0
		self.ymax: float = 0
		self.trimmed: bool = False
		self.has_fit_e: bool = False
		self.baseline: bool = False
		self.reference: bool = False
		self.is_baseline_subtracted: bool = False
		self.is_force_scaled: bool = False
		self.current_plotted_trimmed: str = None # Holds "extension" or "retraction" keywords

		self.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": [], "Fit_Distance_Extension": []})
		self.fit_dataframe_retraction = pd.DataFrame({"Fit_Force_Retraction": [], "Fit_Distance_Retraction": []})

		self.fit_parameters = pd.DataFrame({"Lp_ext": [], "Lc_ext": [], "S_ext": [], "F0_ext": [],"Lp_ret": [], "Lc_ret": [], "S_ret": [], "F0_ret": []})
		self.fc_e: float = 0
		self.fc_r: float = 0

def on_startup() -> None:
	frame_graphing_windows.grid_remove()
	frame_input_buttons.grid_remove()

	frame_session_initialisation.grid(row=1, column=1, padx=2, ipadx=2)

	label_initialisation_blurb.grid(row=0, column=0,columnspan = 3, sticky=tk.NSEW, pady=10)
	label_initialisation_framerate.grid(row=1, column=0, sticky=tk.E)
	label_initialisation_extension_speed.grid(row=2, column=0, sticky=tk.E)
	label_initialisation_extension_speed_units.grid(row=2, column=2, sticky=tk.W)
	label_initialisation_framerate_units.grid(row=2, column=1, sticky=tk.W)

	# Set key optical settings used in the session for savgol filter smoothing
	entry_frame_rate.insert(0, str(GLOBALVARS.frame_rate))
	entry_frame_rate.grid(row=1, column=1, sticky=[tk.W])
	entry_extension_speed.insert(0, str(GLOBALVARS.extension_speed_um_s))
	entry_extension_speed.grid(row=2, column=1, sticky=tk.W)

	# Combobox for selecting whether to use Force 2x / Trap 2 etc
	x_variable_combo["values"] = ["Distance 1", "Distance 2"]
	x_variable_combo.current(0)
	x_variable_combo.grid(row=3, column=0, pady=10)

	y_variable_combo["values"] = ["Force 2x", "Force 2y", "Trap 2"]
	y_variable_combo.current(2)
	y_variable_combo.grid(row=3, column=2, pady=10)

# Handle loading a folder of h5 files into the program
# NOTE: files are kept in RAM. If opening thousands this
# might cause an issue, but it is highly unlikely
def open_folder() -> list[str]:
	folder_path: str = filedialog.askdirectory()
	h5_files: list = []
	if folder_path != "":
		initiate_session()
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
	
def initiate_session() -> None:	
	frame_session_initialisation.grid_remove()

	frame_graphing_windows.grid()
	frame_input_buttons.grid()
	
	GLOBALVARS.frame_rate = float(entry_frame_rate.get())
	GLOBALVARS.extension_speed_um_s = float(entry_extension_speed.get())

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
					scale_select_max_time.set(0)
					scale_select_min_time.set(0)
					#curve.plot()
					replot_canvas()

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
							#curve.plot()
							replot_canvas()
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
				scale_select_max_time.set(0)
				scale_select_min_time.set(0)
			
				if curve.trimmed_e == False and curve.trimmed_r == False:
					variable_radio_buttons_view.set("full")
					variable_checkbutton_set_fit.set(False)
					variable_checkbutton_view_fit.set(False)

				#curve.plot()
				update_trim_entries_ui()
				replot_canvas()

def toggle_time() -> None:
	if GLOBALVARS.active_file != None:
		if (GLOBALVARS.active_file.trimmed_e == False and GLOBALVARS.active_file.trimmed_r == False) or variable_radio_buttons_view.get() == "full":
			max_time: float = max(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])
			min_time: float = min(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])
			max_d: float = max(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])
			min_d: float = min(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])

			if GLOBALVARS.active_file.plot_time == False:
				GLOBALVARS.active_file.plot_time = True
				scale_select_max_time.configure(from_ = min_time, to = max_time)
				scale_select_min_time.configure(from_ = min_time, to = max_time)

				active_index_xmin_e = True
				active_index_xmax_e = True
				active_index_xmin_r = True
				active_index_xmax_r = True
				index_xmin_e = 0
				index_xmax_e = 0
				index_xmin_r = 0
				index_xmax_r = 0
				for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])-1):
					if active_index_xmin_e == True:
						if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] <= GLOBALVARS.active_file.xmin_e[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i+1] > GLOBALVARS.active_file.xmin_e[0]:
							index_xmin_e = i
							active_index_xmin_e = False
					if active_index_xmax_e == True:
						if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] <= GLOBALVARS.active_file.xmax_e[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i+1] > GLOBALVARS.active_file.xmax_e[0]:
							index_xmax_e = i
							active_index_xmax_e = False
					if active_index_xmin_r == True:
						if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] >= GLOBALVARS.active_file.xmin_r[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i+1] < GLOBALVARS.active_file.xmin_r[0]:
							index_xmin_r = i
							#active_index_xmin_r = False
					if active_index_xmax_r == True:
						if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] >= GLOBALVARS.active_file.xmax_r[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i+1] < GLOBALVARS.active_file.xmax_r[0]:
							index_xmax_r = i
							#active_index_xmax_r = False

			
				if GLOBALVARS.active_file.xmin_e[0] != GLOBALVARS.active_file.xmin_e[2]:
					GLOBALVARS.active_file.xmin_e[1] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmin_e]
					GLOBALVARS.active_file.xmin_e[3] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmin_e]
				if GLOBALVARS.active_file.xmax_e[0] != GLOBALVARS.active_file.xmax_e[2]:
					GLOBALVARS.active_file.xmax_e[1] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmax_e]
					GLOBALVARS.active_file.xmax_e[3] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmax_e]
				# Exchange xmin and xmax to retain continuity in the time dimension when going from distance to time
				if GLOBALVARS.active_file.xmin_r[0] != GLOBALVARS.active_file.xmin_r[2]:
					GLOBALVARS.active_file.xmin_r[1] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmax_r]
					GLOBALVARS.active_file.xmin_r[3] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmax_r]
				if GLOBALVARS.active_file.xmax_r[0] != GLOBALVARS.active_file.xmax_r[2]:
					GLOBALVARS.active_file.xmax_r[1] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmin_r]
					GLOBALVARS.active_file.xmax_r[3] = GLOBALVARS.active_file.processed_dataframe["Processed_Time"][index_xmin_r]

				if variable_radio_buttons.get == "extension":
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_e[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_e[1])
				elif variable_radio_buttons.get == "retraction":
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_r[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_r[1])

				update_trim_entries_ui()

			else:
				GLOBALVARS.active_file.plot_time = False
				scale_select_max_time.configure(from_ = min_d, to = max_d)
				scale_select_min_time.configure(from_ = min_d, to = max_d)

				active_index_xmin_e = True
				active_index_xmax_e = True
				active_index_xmin_r = True
				active_index_xmax_r = True
				index_xmin_e = 0
				index_xmax_e = 0
				index_xmin_r = 0
				index_xmax_r = 0
				for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])-1):
					if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmin_e[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i+1] > GLOBALVARS.active_file.xmin_e[1]:
						if active_index_xmin_e == True:
							index_xmin_e = i
							active_index_xmin_e = False
					if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmax_e[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i+1] > GLOBALVARS.active_file.xmax_e[1]:
						if active_index_xmax_e == True:
							index_xmax_e = i
							active_index_xmax_e = False
					if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmin_r[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i+1] > GLOBALVARS.active_file.xmin_r[1]:
						if active_index_xmin_r == True:
							index_xmin_r = i
							#active_index_xmin_r = False
					if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmax_r[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i+1] > GLOBALVARS.active_file.xmax_r[1]:
						if active_index_xmax_r == True:
							index_xmax_r = i
							#active_index_xmax_r = False

				if GLOBALVARS.active_file.xmin_e[1] != GLOBALVARS.active_file.xmin_e[3]:
					GLOBALVARS.active_file.xmin_e[0] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmin_e]
					GLOBALVARS.active_file.xmin_e[2] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmin_e]
				if GLOBALVARS.active_file.xmax_e[1] != GLOBALVARS.active_file.xmax_e[3]:
					GLOBALVARS.active_file.xmax_e[0] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmax_e]
					GLOBALVARS.active_file.xmax_e[2] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmax_e]
				if GLOBALVARS.active_file.xmin_r[1] != GLOBALVARS.active_file.xmin_r[3]:
					GLOBALVARS.active_file.xmin_r[0] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmax_r]
					GLOBALVARS.active_file.xmin_r[2] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmax_r]
				if GLOBALVARS.active_file.xmax_r[1] != GLOBALVARS.active_file.xmax_r[3]:
					GLOBALVARS.active_file.xmax_r[0] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmin_r]
					GLOBALVARS.active_file.xmax_r[2] = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][index_xmin_r]

				if variable_radio_buttons.get == "extension":
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_e[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_e[0])
				elif variable_radio_buttons.get == "retraction":
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_r[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_r[0])

				update_trim_entries_ui()

		else:
			if variable_radio_buttons_view.get() == "extension":
				max_time: float = max(GLOBALVARS.active_file.dataframe_extension["Time_Extension"])
				min_time: float = min(GLOBALVARS.active_file.dataframe_extension["Time_Extension"])
				max_d: float = max(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"])
				min_d: float = min(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"])

				if GLOBALVARS.active_file.plot_time == False:
					GLOBALVARS.active_file.plot_time = True
					scale_select_max_time.configure(from_ = min_time, to = max_time)
					scale_select_min_time.configure(from_ = min_time, to = max_time)
					active_index_xmin_f_e = True
					active_index_xmax_f_e = True
					index_xmin_f_e = 0
					index_xmax_f_e = 0
					for i in range(len(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"])-1):
						if GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i] <= GLOBALVARS.active_file.xmin_f_e[0] and GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i+1] > GLOBALVARS.active_file.xmin_f_e[0]:
							if active_index_xmin_f_e == True:
								index_xmin_f_e = i
								active_index_xmin_f_e = False
						if GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i] <= GLOBALVARS.active_file.xmax_f_e[0] and GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i+1] > GLOBALVARS.active_file.xmax_f_e[0]:
							if active_index_xmax_f_e == True:
								index_xmax_f_e = i
								active_index_xmax_f_e = False
					GLOBALVARS.active_file.xmin_f_e[1] = GLOBALVARS.active_file.dataframe_extension["Time_Extension"][index_xmin_f_e]
					GLOBALVARS.active_file.xmax_f_e[1] = GLOBALVARS.active_file.dataframe_extension["Time_Extension"][index_xmax_f_e]
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_e[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_e[1])

					update_trim_entries_ui()

				else:
					GLOBALVARS.active_file.plot_time = False
					scale_select_max_time.configure(from_ = min_d, to = max_d)
					scale_select_min_time.configure(from_ = min_d, to = max_d)
					active_index_xmin_f_e = True
					active_index_xmax_f_e = True
					index_xmin_f_e = 0
					index_xmax_f_e = 0
					for i in range(len(GLOBALVARS.active_file.dataframe_extension["Time_Extension"])-1):
						if GLOBALVARS.active_file.dataframe_extension["Time_Extension"].iloc[i] <= GLOBALVARS.active_file.xmin_f_e[1] and GLOBALVARS.active_file.dataframe_extension["Time_Extension"].iloc[i+1] > GLOBALVARS.active_file.xmin_f_e[1]:
							if active_index_xmin_f_e == True:
								index_xmin_f_e = i
								active_index_xmin_f_e = False
						if GLOBALVARS.active_file.dataframe_extension["Time_Extension"].iloc[i] <= GLOBALVARS.active_file.xmax_f_e[1] and GLOBALVARS.active_file.dataframe_extension["Time_Extension"].iloc[i+1] > GLOBALVARS.active_file.xmax_f_e[1]:
							if active_index_xmax_f_e == True:
								index_xmax_f_e = i
								active_index_xmax_f_e = False
					GLOBALVARS.active_file.xmin_f_e[0] = GLOBALVARS.active_file.dataframe_extension["Distance_Extension"].iloc[index_xmin_f_e]
					GLOBALVARS.active_file.xmax_f_e[0] = GLOBALVARS.active_file.dataframe_extension["Distance_Extension"].iloc[index_xmax_f_e]
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_e[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_e[0])

					update_trim_entries_ui()


			elif variable_radio_buttons_view.get() == "retraction":
				max_time: float = max(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"])
				min_time: float = min(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"])
				max_d: float = max(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"])
				min_d: float = min(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"])

				if GLOBALVARS.active_file.plot_time == False:
					GLOBALVARS.active_file.plot_time = True
					scale_select_max_time.configure(from_ = min_time, to = max_time)
					scale_select_min_time.configure(from_ = min_time, to = max_time)
					active_index_xmin_f_r = True
					active_index_xmax_f_r = True
					index_xmin_f_r = 0
					index_xmax_f_r = 0
					for i in range(len(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"])-1):
						if GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i] >= GLOBALVARS.active_file.xmin_f_r[0] and GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i+1] < GLOBALVARS.active_file.xmin_f_r[0]:
							if active_index_xmin_f_r == True:
								index_xmin_f_r = i
								#active_index_xmin_f_r = False
						if GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i] >= GLOBALVARS.active_file.xmax_f_r[0] and GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i+1] < GLOBALVARS.active_file.xmax_f_r[0]:
							if active_index_xmax_f_r == True:
								index_xmax_f_r = i
								#active_index_xmax_f_r = False
					GLOBALVARS.active_file.xmin_f_r[1] = GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][index_xmax_f_r]
					GLOBALVARS.active_file.xmax_f_r[1] = GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][index_xmin_f_r]
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_r[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_r[1])

					update_trim_entries_ui()

				else:
					GLOBALVARS.active_file.plot_time = False
					scale_select_max_time.configure(from_ = min_d, to = max_d)
					scale_select_min_time.configure(from_ = min_d, to = max_d)
					active_index_xmin_f_r = True
					active_index_xmax_f_r = True
					index_xmin_f_r = 0
					index_xmax_f_r = 0
					for i in range(len(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"])-1):
						if GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i] <= GLOBALVARS.active_file.xmin_f_r[1] and GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i+1] > GLOBALVARS.active_file.xmin_f_r[1]:
							if active_index_xmin_f_r == True:
								index_xmin_f_r = i
								#active_index_xmin_f_r = False
						if GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i] <= GLOBALVARS.active_file.xmax_f_r[1] and GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i+1] > GLOBALVARS.active_file.xmax_f_r[1]:
							if active_index_xmax_f_r == True:
								index_xmax_f_r = i
								#active_index_xmax_f_r = False
					GLOBALVARS.active_file.xmin_f_r[0] = GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][index_xmax_f_r]
					GLOBALVARS.active_file.xmax_f_r[0] = GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][index_xmin_f_r]
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_r[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_r[0])

					update_trim_entries_ui()

		replot_canvas()


def slider_max_release(event) -> None:
	if GLOBALVARS.active_file != None:
		if GLOBALVARS.active_file.plot_time == False:
			index = 0
		else:
			index = 1
		if variable_radio_buttons.get() == "extension":
			if variable_checkbutton_set_fit.get() == False:
				GLOBALVARS.active_file.xmax_e[index] = scale_select_max_time.get()
				GLOBALVARS.active_file.ymin_e = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_e = float(entry_ymax.get())
				entry_xmax.delete(0, tk.END)
				entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_e[index]))
			else:
				GLOBALVARS.active_file.xmax_f_e[index] = scale_select_max_time.get()
				GLOBALVARS.active_file.ymin_f_e = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_f_e = float(entry_ymax.get())
				entry_xmax.delete(0, tk.END)
				entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_f_e[index]))
				
		elif variable_radio_buttons.get() == "retraction":
			if variable_checkbutton_set_fit.get() == False:
				GLOBALVARS.active_file.xmax_r[index] = scale_select_max_time.get()
				GLOBALVARS.active_file.ymin_r = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_r = float(entry_ymax.get())
				entry_xmax.delete(0, tk.END)
				entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_r[index]))
			else:
				GLOBALVARS.active_file.xmax_f_r[index] = scale_select_max_time.get()
				GLOBALVARS.active_file.ymin_f_r = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_f_r = float(entry_ymax.get())
				entry_xmax.delete(0, tk.END)
				entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_f_r[index]))

		replot_canvas()

def slider_min_release(event) -> None:
	if GLOBALVARS.active_file != None:
		if GLOBALVARS.active_file.plot_time == False:
			index = 0
		else:
			index = 1

		if variable_radio_buttons.get() == "extension":
			if variable_checkbutton_set_fit.get() == False:
				GLOBALVARS.active_file.xmin_e[index] = scale_select_min_time.get()
				GLOBALVARS.active_file.ymin_e = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_e = float(entry_ymax.get())
				entry_xmin.delete(0, tk.END)
				entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_e[index]))
			else:
				GLOBALVARS.active_file.xmin_f_e[index] = scale_select_min_time.get()
				GLOBALVARS.active_file.ymin_f_e = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_f_e = float(entry_ymax.get())
				entry_xmin.delete(0, tk.END)
				entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_f_e[index]))

		elif variable_radio_buttons.get() == "retraction":
			if variable_checkbutton_set_fit.get() == False:
				GLOBALVARS.active_file.xmin_r[index] = scale_select_min_time.get()
				GLOBALVARS.active_file.ymin_r = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_r = float(entry_ymax.get())
				entry_xmin.delete(0, tk.END)
				entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_r[index]))
			else:
				GLOBALVARS.active_file.xmin_f_r[index] = scale_select_min_time.get()
				GLOBALVARS.active_file.ymin_f_r = float(entry_ymin.get())
				GLOBALVARS.active_file.ymax_f_r = float(entry_ymax.get())
				entry_xmin.delete(0, tk.END)
				entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_f_r[index]))
				
		replot_canvas()

def update_trim_settings() -> None:
	if GLOBALVARS.active_file != None:
		if variable_radio_buttons.get() == "extension":
			if variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == False:
					GLOBALVARS.active_file.xmin_e[0] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_e[0] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_e = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_e = float(entry_ymax.get())
				else:
					GLOBALVARS.active_file.xmin_e[1] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_e[1] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_e = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_e = float(entry_ymax.get())
			else:
				if GLOBALVARS.active_file.plot_time == False:
					GLOBALVARS.active_file.xmin_f_e[0] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_f_e[0] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_f_e = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_f_e = float(entry_ymax.get())
				else:
					GLOBALVARS.active_file.xmin_f_e[1] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_f_e[1] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_f_e = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_f_e = float(entry_ymax.get())
		elif variable_radio_buttons.get() == "retraction":
			if variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == False:
					GLOBALVARS.active_file.xmin_r[0] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_r[0] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_r = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_r = float(entry_ymax.get())
				else:
					GLOBALVARS.active_file.xmin_r[1] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_r[1] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_r = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_r = float(entry_ymax.get())
			else:
				if GLOBALVARS.active_file.plot_time == False:
					GLOBALVARS.active_file.xmin_f_r[0] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_f_r[0] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_f_r = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_f_r = float(entry_ymax.get())
				else:
					GLOBALVARS.active_file.xmin_f_r[1] = float(entry_xmin.get())
					GLOBALVARS.active_file.xmax_f_r[1] = float(entry_xmax.get())
					GLOBALVARS.active_file.ymin_f_r = float(entry_ymin.get())
					GLOBALVARS.active_file.ymax_f_r = float(entry_ymax.get())
		
		scale_select_max_time.set(float(entry_xmax.get()))
		scale_select_min_time.set(float(entry_xmin.get()))

		replot_canvas()

def update_trim_entries_ui() -> None:
	if GLOBALVARS.active_file != None:
		if variable_radio_buttons.get() == "extension":
			if variable_checkbutton_set_fit.get() == False:
				entry_xmin.delete(0, tk.END)
				entry_xmax.delete(0, tk.END)
				entry_ymin.delete(0, tk.END)
				entry_ymax.delete(0, tk.END)
				if GLOBALVARS.active_file.plot_time == False:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_e[0]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_e[0]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_e))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_e))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_e[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_e[0])
				else:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_e[1]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_e[1]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_e))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_e))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_e[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_e[1])
			else:
				entry_xmin.delete(0, tk.END)
				entry_xmax.delete(0, tk.END)
				entry_ymin.delete(0, tk.END)
				entry_ymax.delete(0, tk.END)
				if GLOBALVARS.active_file.plot_time == False:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_f_e[0]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_f_e[0]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_f_e))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_f_e))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_e[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_e[0])
				else:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_f_e[1]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_f_e[1]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_f_e))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_f_e))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_e[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_e[1])
		else:
			if variable_checkbutton_set_fit.get() == False:
				entry_xmin.delete(0, tk.END)
				entry_xmax.delete(0, tk.END)
				entry_ymin.delete(0, tk.END)
				entry_ymax.delete(0, tk.END)
				if GLOBALVARS.active_file.plot_time == False:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_r[0]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_r[0]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_r))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_r))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_r[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_r[0])
				else:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_r[1]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_r[1]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_r))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_r))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_r[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_r[1])
			else:
				entry_xmin.delete(0, tk.END)
				entry_xmax.delete(0, tk.END)
				entry_ymin.delete(0, tk.END)
				entry_ymax.delete(0, tk.END)
				if GLOBALVARS.active_file.plot_time == False:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_f_r[0]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_f_r[0]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_f_r))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_f_r))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_r[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_r[0])
				else:
					entry_xmin.insert(0, str(GLOBALVARS.active_file.xmin_f_r[1]))
					entry_xmax.insert(0, str(GLOBALVARS.active_file.xmax_f_r[1]))
					entry_ymin.insert(0, str(GLOBALVARS.active_file.ymin_f_r))
					entry_ymax.insert(0, str(GLOBALVARS.active_file.ymax_f_r))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_r[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_r[1])
				
				

def update_canvas() -> None:
	GLOBALVARS.graph_image = tk.PhotoImage(file="TEMP_PLOT.png")
	canvas_graph_display.create_image(0,0,image=GLOBALVARS.graph_image, anchor="nw")

def replot_canvas(expanded_graph = False) -> None:
	if GLOBALVARS.active_file != None:
		GLOBALVARS.active_file.plot(expanded_graph)
		GLOBALVARS.graph_image = tk.PhotoImage(file="TEMP_PLOT.png")
		canvas_graph_display.create_image(0,0,image=GLOBALVARS.graph_image, anchor="nw")

		if GLOBALVARS.active_file.sigma_e != -100 and variable_radio_buttons_view.get() == "extension":
			sigma_display.config(state="normal")
			sigma_display.delete(0,tk.END)
			sigma_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.sigma_e)))
			sigma_display.config(state="readonly")
		elif GLOBALVARS.active_file.sigma_r != -100 and variable_radio_buttons_view.get() == "retraction":
			sigma_display.config(state="normal")
			sigma_display.delete(0,tk.END)
			sigma_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.sigma_r)))
			sigma_display.config(state="readonly")

		Lp_display.config(state="normal")
		Lc_display.config(state="normal")
		S_display.config(state="normal")
		F0_display.config(state="normal")

		Lp_display.delete(0,tk.END)
		Lc_display.delete(0,tk.END)
		S_display.delete(0,tk.END)
		F0_display.delete(0,tk.END)

		if GLOBALVARS.active_file.has_fit_e == True and variable_radio_buttons_view.get() == "extension":
				Lp_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lp_ext"][0])))
				Lc_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lc_ext"][0])))
				S_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["S_ext"][0])))
				F0_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["F0_ext"][0])))
		if GLOBALVARS.active_file.has_fit_r == True and variable_radio_buttons_view.get() == "retraction":
				Lp_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lp_ret"][0])))
				Lc_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["Lc_ret"][0])))
				S_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["S_ret"][0])))
				F0_display.insert(0,str(float("%.4g" % GLOBALVARS.active_file.fit_parameters["F0_ret"][0])))

		Lp_display.config(state="readonly")
		Lc_display.config(state="readonly")
		S_display.config(state="readonly")
		F0_display.config(state="readonly")


def window_resize(event) -> None:
	replot_canvas()

def update_optic_settings() -> None:

	scale_factor = float(entry_frame_rate.get())/GLOBALVARS.frame_rate
	GLOBALVARS.frame_rate = float(entry_frame_rate.get())
	for curve in GLOBALVARS.all_files:
		curve.dataframe = pd.DataFrame({"Force": curve.dataframe["Force"], "Distance": curve.dataframe["Distance"], "Time": [i/scale_factor for i in curve.dataframe["Time"]]})
		curve.first_derivative_dataframe = pd.DataFrame({"First_Derivative": curve.first_derivative_dataframe["First_Derivative"], "Distance": curve.first_derivative_dataframe["Distance"], "Time": [i/scale_factor for i in curve.first_derivative_dataframe["Time"]]})
		curve.second_derivative_dataframe = pd.DataFrame({"Second_Derivative": curve.second_derivative_dataframe["Second_Derivative"], "Distance": curve.second_derivative_dataframe["Distance"], "Time": [i/scale_factor for i in curve.second_derivative_dataframe["Time"]]})
		curve.processed_dataframe = pd.DataFrame({"Processed_Force": curve.processed_dataframe["Processed_Force"], "Processed_Distance": curve.processed_dataframe["Processed_Distance"], "Processed_Time": [i/scale_factor for i in curve.processed_dataframe["Processed_Time"]]})
		curve.dataframe_extension = pd.DataFrame({"Force_Extension": curve.dataframe_extension["Force_Extension"], "Distance_Extension": curve.dataframe_extension["Distance_Extension"], "Time_Extension": [i/scale_factor for i in curve.dataframe_extension["Time_Extension"]]})
		curve.dataframe_retraction = pd.DataFrame({"Force_Retraction": curve.dataframe_retraction["Force_Retraction"], "Distance_Retraction": curve.dataframe_retraction["Distance_Retraction"], "Time_Retraction": [i/scale_factor for i in curve.dataframe_retraction["Time_Retraction"]]})

	if GLOBALVARS.active_file != None:
		GLOBALVARS.active_file.ymax = float(entry_ymax.get())
		if GLOBALVARS.active_file.plot_time == True:
			scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Time"]), to = max(GLOBALVARS.active_file.dataframe["Time"]))
			scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Time"]), to = max(GLOBALVARS.active_file.dataframe["Time"]))
		else:
			scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Distance"]), to = max(GLOBALVARS.active_file.dataframe["Distance"]))
			scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe["Distance"]), to = max(GLOBALVARS.active_file.dataframe["Distance"]))
	replot_canvas()

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
		else:
			curve.reset()
	GLOBALVARS.selected_files = updated_selection
	listbox_all_selected_files.delete(0,tk.END)
	for curve in updated_selection:
		listbox_all_selected_files.insert(tk.END, curve.name)
		if curve.baseline == True:
			listbox_all_selected_files.itemconfig(tk.END, fg = "red")	
		elif curve.reference == True:
			listbox_all_selected_files.itemconfig(tk.END, fg = "green")
			

def toggle_first_derivative() -> None:
	if GLOBALVARS.show_first_deriv == True:
		GLOBALVARS.show_first_deriv = False
	else:
		GLOBALVARS.show_first_deriv = True
	replot_canvas()
def toggle_second_derivative() -> None:
	if GLOBALVARS.show_second_deriv == True:
		GLOBALVARS.show_second_deriv = False
	else:
		GLOBALVARS.show_second_deriv = True
	replot_canvas()


def manual_trim_data() -> None:
	if GLOBALVARS.active_file in GLOBALVARS.selected_files:

		# Delete any fits associated with extension or retraction curves
		if variable_radio_buttons.get() == "extension":
			if GLOBALVARS.active_file.has_fit_e == True:
				GLOBALVARS.active_file.has_fit_e = False
				GLOBALVARS.active_file.trimmed_f_e = False
				GLOBALVARS.active_file.precalculated_fit_extension= pd.DataFrame({"Precalculated_Fit_Force_Extension": [], "Precalculated_Fit_Distance_Extension": [], "Precalculated_Fit_Time_Extension": []})
				GLOBALVARS.active_file.fit_parameters = pd.DataFrame({"Lp_ext": [0,0], "Lc_ext": [0,0], "S_ext": [0,0],"F0_ext": [0,0],
				"Lp_ret": GLOBALVARS.active_file.fit_parameters["Lp_ret"],
				"Lc_ret": GLOBALVARS.active_file.fit_parameters["Lc_ret"],
				"S_ret": GLOBALVARS.active_file.fit_parameters["S_ret"],
				"F0_ret": GLOBALVARS.active_file.fit_parameters["F0_ret"]})
		elif variable_radio_buttons.get() == "retraction":
			if GLOBALVARS.active_file.has_fit_r == True:
				GLOBALVARS.active_file.has_fit_r = False
				GLOBALVARS.active_file.trimmed_f_r = False
				GLOBALVARS.active_file.precalculated_fit_retraction= pd.DataFrame({"Precalculated_Fit_Force_Retraction": [], "Precalculated_Fit_Distance_Retraction": [], "Precalculated_Fit_Time_Retraction": []})
				GLOBALVARS.active_file.fit_parameters = pd.DataFrame({"Lp_ret": [0,0], "Lc_ret": [0,0], "S_ret": [0,0],"F0_ret": [0,0],
				"Lp_ext": GLOBALVARS.active_file.fit_parameters["Lp_ext"],
				"Lc_ext": GLOBALVARS.active_file.fit_parameters["Lc_ext"],
				"S_ext": GLOBALVARS.active_file.fit_parameters["S_ext"],
				"F0_ext": GLOBALVARS.active_file.fit_parameters["F0_ext"]})
		
		trimmed_force = []
		trimmed_time = []
		trimmed_dist = []
		trimmed_first_deriv_force = []
		trimmed_first_deriv_distance = []
		trimmed_first_deriv_time = []

		# If we are trimming an extension curve...
		if variable_radio_buttons.get() == "extension":
			variable_radio_buttons_view.set("extension")
			# but we are not setting a fit range,
			if variable_checkbutton_set_fit.get() == False:
				# if we have plotted time, use the processed time as the x axis range
				if GLOBALVARS.active_file.plot_time == True:
					for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])):
						if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] >= GLOBALVARS.active_file.xmin_e[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmax_e[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax_e and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] >= GLOBALVARS.active_file.ymin_e:
							trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
							trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
							trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
							trimmed_first_deriv_force.append(GLOBALVARS.active_file.first_derivative_dataframe["First_Derivative"][i])
							trimmed_first_deriv_distance.append(GLOBALVARS.active_file.first_derivative_dataframe["Distance"][i])
							trimmed_first_deriv_time.append(GLOBALVARS.active_file.first_derivative_dataframe["Time"][i])
					GLOBALVARS.active_file.dataframe_extension = pd.DataFrame({"Force_Extension": [], "Distance_Extension": [], "Time_Extension": []})
					GLOBALVARS.active_file.dataframe_extension["Force_Extension"] = trimmed_force
					GLOBALVARS.active_file.dataframe_extension["Time_Extension"] = trimmed_time
					GLOBALVARS.active_file.dataframe_extension["Distance_Extension"] = trimmed_dist
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension = pd.DataFrame({"Trimmed_First_Derivative": [], "Trimmed_Distance": [], "Trimmed_Time": []})
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"] = trimmed_first_deriv_force
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension["Trimmed_Distance"] = trimmed_first_deriv_distance
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension["Trimmed_Time"] = trimmed_first_deriv_time

					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]))

				else: # we are trimming in distance, use processed distance
					for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])):
						if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] >= GLOBALVARS.active_file.xmin_e[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] <= GLOBALVARS.active_file.xmax_e[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax_e and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] >= GLOBALVARS.active_file.ymin_e:
							trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
							trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
							trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
							trimmed_first_deriv_force.append(GLOBALVARS.active_file.first_derivative_dataframe["First_Derivative"][i])
							trimmed_first_deriv_distance.append(GLOBALVARS.active_file.first_derivative_dataframe["Distance"][i])
							trimmed_first_deriv_time.append(GLOBALVARS.active_file.first_derivative_dataframe["Time"][i])
					GLOBALVARS.active_file.dataframe_extension = pd.DataFrame({"Force_Extension": [], "Distance_Extension": [], "Time_Extension": []})
					GLOBALVARS.active_file.dataframe_extension["Force_Extension"] = trimmed_force
					GLOBALVARS.active_file.dataframe_extension["Time_Extension"] = trimmed_time
					GLOBALVARS.active_file.dataframe_extension["Distance_Extension"] = trimmed_dist
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension= pd.DataFrame({"Trimmed_First_Derivative": [], "Trimmed_Distance": [], "Trimmed_Time": []})
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"] = trimmed_first_deriv_force
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension["Trimmed_Distance"] = trimmed_first_deriv_distance
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension["Trimmed_Time"] = trimmed_first_deriv_time
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]))

				GLOBALVARS.active_file.trimmed_e = True
				variable_checkbutton_set_fit.set(True)
				update_trim_entries_ui()
				replot_canvas()
			else: # we are trimming the extension curve to designate a propper fit range, therefore use dataframe_extension not processed_data
				if GLOBALVARS.active_file.plot_time == True:
					for i in range(len(GLOBALVARS.active_file.dataframe_extension["Time_Extension"])):
						if GLOBALVARS.active_file.dataframe_extension["Time_Extension"][i] >= GLOBALVARS.active_file.xmin_f_e[1] and GLOBALVARS.active_file.dataframe_extension["Time_Extension"][i] <= GLOBALVARS.active_file.xmax_f_e[1] and GLOBALVARS.active_file.dataframe_extension["Force_Extension"][i] <= GLOBALVARS.active_file.ymax_f_e and GLOBALVARS.active_file.dataframe_extension["Force_Extension"][i] >= GLOBALVARS.active_file.ymin_f_e:
							trimmed_force.append(GLOBALVARS.active_file.dataframe_extension["Force_Extension"][i])
							trimmed_dist.append(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i])
							trimmed_time.append(GLOBALVARS.active_file.dataframe_extension["Time_Extension"][i])
					GLOBALVARS.active_file.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": [], "Fit_Distance_Extension": []})
					GLOBALVARS.active_file.fit_dataframe_extension["Fit_Force_Extension"] = trimmed_force
					GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"] = trimmed_dist
					GLOBALVARS.active_file.fit_dataframe_extension["Fit_Time_Extension"] = trimmed_time
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Time_Extension"]), to = max(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Time_Extension"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Time_Extension"]), to = max(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Time_Extension"]))

				else:
					for i in range(len(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"])):
						if GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i] >= GLOBALVARS.active_file.xmin_f_e[0] and GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i] <= GLOBALVARS.active_file.xmax_f_e[0] and GLOBALVARS.active_file.dataframe_extension["Force_Extension"][i] <= GLOBALVARS.active_file.ymax_f_e and GLOBALVARS.active_file.dataframe_extension["Force_Extension"][i] >= GLOBALVARS.active_file.ymin_f_e:
							trimmed_force.append(GLOBALVARS.active_file.dataframe_extension["Force_Extension"][i])
							trimmed_dist.append(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"][i])
							trimmed_time.append(GLOBALVARS.active_file.dataframe_extension["Time_Extension"][i])
					GLOBALVARS.active_file.fit_dataframe_extension= pd.DataFrame({"Fit_Force_Extension": [], "Fit_Distance_Extension": []})
					GLOBALVARS.active_file.fit_dataframe_extension["Fit_Force_Extension"] = trimmed_force
					GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"] = trimmed_dist
					GLOBALVARS.active_file.fit_dataframe_extension["Fit_Time_Extension"] = trimmed_time
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"]), to = max(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"]), to = max(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"]))

				GLOBALVARS.active_file.trimmed_f_e = True
				variable_checkbutton_view_fit.set(True)
				update_trim_entries_ui()
				replot_canvas()
				

		elif variable_radio_buttons.get() == "retraction":
			variable_radio_buttons_view.set("retraction")
			if variable_checkbutton_set_fit.get() == False: # We are trimming to create a retraction curve, but not specifying the fit range.
				if GLOBALVARS.active_file.plot_time == True:
					for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Time"])):
						if GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] >= GLOBALVARS.active_file.xmin_r[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i] <= GLOBALVARS.active_file.xmax_r[1] and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax_r:
							trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
							trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
							trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
							trimmed_first_deriv_force.append(GLOBALVARS.active_file.first_derivative_dataframe["First_Derivative"][i])
							trimmed_first_deriv_distance.append(GLOBALVARS.active_file.first_derivative_dataframe["Distance"][i])
							trimmed_first_deriv_time.append(GLOBALVARS.active_file.first_derivative_dataframe["Time"][i])
					GLOBALVARS.active_file.dataframe_retraction = pd.DataFrame({"Force_Retraction": [], "Distance_Retraction": [], "Time_Retraction": []})
					GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"] = trimmed_force
					GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"] = trimmed_time
					GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"] = trimmed_dist
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction= pd.DataFrame({"Trimmed_First_Derivative": [], "Trimmed_Distance": [], "Trimmed_Time": []})
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"] = trimmed_first_deriv_force
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction["Trimmed_Distance"] = trimmed_first_deriv_distance
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction["Trimmed_Time"] = trimmed_first_deriv_time

					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]))
				else:
					for i in range(len(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])):
						if GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] >= GLOBALVARS.active_file.xmin_r[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i] <= GLOBALVARS.active_file.xmax_r[0] and GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i] <= GLOBALVARS.active_file.ymax_r:
							trimmed_force.append(GLOBALVARS.active_file.processed_dataframe["Processed_Force"][i])
							trimmed_time.append(GLOBALVARS.active_file.processed_dataframe["Processed_Time"][i])
							trimmed_dist.append(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"][i])
							trimmed_first_deriv_force.append(GLOBALVARS.active_file.first_derivative_dataframe["First_Derivative"][i])
							trimmed_first_deriv_distance.append(GLOBALVARS.active_file.first_derivative_dataframe["Distance"][i])
							trimmed_first_deriv_time.append(GLOBALVARS.active_file.first_derivative_dataframe["Time"][i])
					GLOBALVARS.active_file.dataframe_retraction = pd.DataFrame({"Force_Retraction": [], "Distance_Retraction": [], "Time_Retraction": []})
					GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"] = trimmed_force
					GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"] = trimmed_time
					GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"] = trimmed_dist
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction= pd.DataFrame({"Trimmed_First_Derivative": [], "Trimmed_Distance": [], "Trimmed_Time": []})
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"] = trimmed_first_deriv_force
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction["Trimmed_Distance"] = trimmed_first_deriv_distance
					GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction["Trimmed_Time"] = trimmed_first_deriv_time

					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]))

				GLOBALVARS.active_file.trimmed_r = True
				variable_checkbutton_set_fit.set(True)
				update_trim_entries_ui()
				replot_canvas()
			else: # we are trimming the retraction curve to create the fit range!
				if GLOBALVARS.active_file.plot_time == True:
					for i in range(len(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"])):
						if GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i] >= GLOBALVARS.active_file.xmin_f_r[1] and GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i] <= GLOBALVARS.active_file.xmax_f_r[1] and GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"][i] <= GLOBALVARS.active_file.ymax_f_r and GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"][i] >= GLOBALVARS.active_file.ymin_f_r:
							trimmed_force.append(GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"][i])
							trimmed_dist.append(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i])
							trimmed_time.append(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i])
					GLOBALVARS.active_file.fit_dataframe_retraction = pd.DataFrame({"Fit_Force_Retraction": [], "Fit_Distance_Retraction": [],"Fit_Time_Retraction": [] })
					GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Force_Retraction"] = trimmed_force
					GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"] = trimmed_dist
					GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Time_Retraction"] = trimmed_time
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Time_Retraction"]), to = max(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Time_Retraction"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Time_Retraction"]), to = max(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Time_Retraction"]))

				else:
					for i in range(len(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"])):
						if GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i] >= GLOBALVARS.active_file.xmin_f_r[0] and GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i] <= GLOBALVARS.active_file.xmax_f_r[0] and GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"][i] <= GLOBALVARS.active_file.ymax_f_r and GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"][i] >= GLOBALVARS.active_file.ymin_f_r:
							trimmed_force.append(GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"][i])
							trimmed_dist.append(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"][i])
							trimmed_time.append(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"][i])
					GLOBALVARS.active_file.fit_dataframe_retraction= pd.DataFrame({"Fit_Force_Retraction": [], "Fit_Distance_Retraction": [],"Fit_Time_Retraction": []})
					GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Force_Retraction"] = trimmed_force
					GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"] = trimmed_dist
					GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Time_Retraction"] = trimmed_time
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"]), to = max(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"]), to = max(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"]))

				GLOBALVARS.active_file.trimmed_f_r = True
				variable_checkbutton_view_fit.set(True)
				update_trim_entries_ui()
				replot_canvas()



def radio_button_select() -> None:
	if GLOBALVARS.active_file != None:
		if variable_radio_buttons_view.get() != "full":
			if GLOBALVARS.active_file.trimmed_e == True and variable_radio_buttons_view.get() == "extension" and variable_checkbutton_set_fit.get() == False:
				variable_checkbutton_set_fit.set(True)
				variable_checkbutton_view_fit.set(True)

			elif GLOBALVARS.active_file.trimmed_r == True and variable_radio_buttons_view.get() == "retraction" and variable_checkbutton_set_fit.get() == False:
				variable_checkbutton_set_fit.set(True)
				variable_checkbutton_view_fit.set(True)

			elif GLOBALVARS.active_file.trimmed_e == False and variable_radio_buttons_view.get() == "extension" and variable_checkbutton_set_fit.get() == True:
				variable_checkbutton_set_fit.set(False)
				variable_checkbutton_view_fit.set(False)
			elif GLOBALVARS.active_file.trimmed_r == False and variable_radio_buttons_view.get() == "retraction" and variable_checkbutton_set_fit.get() == True:
				variable_checkbutton_set_fit.set(False)
				variable_checkbutton_view_fit.set(False)

			if GLOBALVARS.active_file.trimmed_e == True and variable_radio_buttons_view.get() == "extension" and variable_checkbutton_set_fit.get() == True:
				if GLOBALVARS.active_file.plot_time == False:
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_e[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_e[0])
				else:
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]), to = max(GLOBALVARS.active_file.dataframe_extension["Time_Extension"]))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_e[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_e[1])


				update_trim_entries_ui()
				replot_canvas()

			elif GLOBALVARS.active_file.trimmed_r == True and variable_radio_buttons_view.get() == "retraction" and variable_checkbutton_set_fit.get() == True:
				if GLOBALVARS.active_file.plot_time == False:
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_r[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_r[0])
				else:
					scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]))
					scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]), to = max(GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]))
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_f_r[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_f_r[1])


				update_trim_entries_ui()
				replot_canvas()
		else:
			variable_checkbutton_set_fit.set(False)
			variable_checkbutton_view_fit.set(False)
			if GLOBALVARS.active_file.plot_time == False:
				scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]), to = max(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]))
				scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]), to = max(GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]))
			else:
				scale_select_max_time.configure(from_ = min(GLOBALVARS.active_file.processed_dataframe["Processed_Time"]), to = max(GLOBALVARS.active_file.processed_dataframe["Processed_Time"]))
				scale_select_min_time.configure(from_ = min(GLOBALVARS.active_file.processed_dataframe["Processed_Time"]), to = max(GLOBALVARS.active_file.processed_dataframe["Processed_Time"]))

			if variable_radio_buttons.get() == "extension":
				if GLOBALVARS.active_file.plot_time == False:
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_e[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_e[0])
				else:
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_e[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_e[1])
			else:
				if GLOBALVARS.active_file.plot_time == False:
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_r[0])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_r[0])
				else:
					scale_select_max_time.set(GLOBALVARS.active_file.xmax_r[1])
					scale_select_min_time.set(GLOBALVARS.active_file.xmin_r[1])

			update_trim_entries_ui()
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


	# Allow parameter "fixing" by definding a new lambda function with that as a constant.
	parameters_to_fix: list[bool] = [variable_lp_fix.get(),variable_lc_fix.get(),variable_s_fix.get(),variable_f0_fix.get()]

	if parameters_to_fix[0] == True:
		fixed_lp = lambda d, Lp, Lc, S, F0: eOdjik_force_offset(d, float(Lp_entry.get()), Lc, S, F0)
	else:
		fixed_lp = lambda d, Lp, Lc, S, F0: eOdjik_force_offset(d, Lp, Lc, S, F0)
	if parameters_to_fix[1] == True:
		fixed_lc = lambda d, Lp, Lc, S, F0: fixed_lp(d, Lp, float(Lc_entry.get()), S, F0)
	else:
		fixed_lc = lambda d, Lp, Lc, S, F0: fixed_lp(d, Lp, Lc, S, F0)
	if parameters_to_fix[2] == True:
		fixed_s = lambda d, Lp, Lc, S, F0: fixed_lc(d, Lp, Lc, float(S_entry.get()), F0)
	else:
		fixed_s = lambda d, Lp, Lc, S, F0: fixed_lc(d, Lp, Lc, S, F0)
	if parameters_to_fix[3] == True:
		fixed_f0 = lambda d, Lp, Lc, S, F0: fixed_s(d, Lp, Lc, S, float(F0_entry.get()))
	else:
		fixed_f0 = lambda d, Lp, Lc, S, F0: fixed_s(d, Lp, Lc, S, F0)

	parameters, covariance = scipy.optimize.curve_fit(fixed_f0,np.array(data_x),np.array(data_y), (float(Lp_entry.get()), float(Lc_entry.get()), float(S_entry.get()), float(F0_entry.get())))
	errors = np.sqrt(np.diag(covariance))
	predicted_force = eOdjik_force_offset(np.array(data_x), *parameters)
	residuals = np.array(data_y) - predicted_force
	return [parameters, errors, predicted_force, residuals]


def fit() -> None:
	if GLOBALVARS.active_file in GLOBALVARS.selected_files:
		if GLOBALVARS.active_file.trimmed_f_e == True and variable_radio_buttons_view.get() == "extension":
				fit_result = fit_eOdijk_F0(GLOBALVARS.active_file.fit_dataframe_extension["Fit_Distance_Extension"], GLOBALVARS.active_file.fit_dataframe_extension["Fit_Force_Extension"])
				GLOBALVARS.active_file.precalculated_fit_extension= pd.DataFrame({"Precalculated_Fit_Force_Extension": [], "Precalculated_Fit_Distance_Extension": [], "Precalculated_Fit_Time_Extension": []})

				GLOBALVARS.active_file.precalculated_fit_extension["Precalculated_Fit_Force_Extension"] = fit_result[2]
				GLOBALVARS.active_file.precalculated_fit_extension["Precalculated_Fit_Distance_Extension"] = GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]
				GLOBALVARS.active_file.precalculated_fit_extension["Precalculated_Fit_Time_Extension"] = GLOBALVARS.active_file.fit_dataframe_extension["Fit_Time_Extension"]
				GLOBALVARS.active_file.fit_parameters["Lp_ext"] = [fit_result[0][0], fit_result[1][0]]
				GLOBALVARS.active_file.fit_parameters["Lc_ext"] = [fit_result[0][1], fit_result[1][1]]
				GLOBALVARS.active_file.fit_parameters["S_ext"] = [fit_result[0][2], fit_result[1][2]]
				GLOBALVARS.active_file.fit_parameters["F0_ext"] = [fit_result[0][3], fit_result[1][3]]

				GLOBALVARS.active_file.has_fit_e = True

		elif GLOBALVARS.active_file.trimmed_f_r == True and variable_radio_buttons_view.get() == "retraction":
				fit_result = fit_eOdijk_F0(GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"], GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Force_Retraction"])
				GLOBALVARS.active_file.precalculated_fit_retraction = pd.DataFrame({"Precalculated_Fit_Force_Retraction": [], "Precalculated_Fit_Distance_Retraction": [], "Precalculated_Fit_Time_Retraction": []})

				GLOBALVARS.active_file.precalculated_fit_retraction["Precalculated_Fit_Force_Retraction"] = fit_result[2]
				GLOBALVARS.active_file.precalculated_fit_retraction["Precalculated_Fit_Distance_Retraction"] = GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Distance_Retraction"]
				GLOBALVARS.active_file.precalculated_fit_retraction["Precalculated_Fit_Time_Retraction"] = GLOBALVARS.active_file.fit_dataframe_retraction["Fit_Time_Retraction"]
				GLOBALVARS.active_file.fit_parameters["Lp_ret"] = [fit_result[0][0], fit_result[1][0]]
				GLOBALVARS.active_file.fit_parameters["Lc_ret"] = [fit_result[0][1], fit_result[1][1]]
				GLOBALVARS.active_file.fit_parameters["S_ret"] = [fit_result[0][2], fit_result[1][2]]
				GLOBALVARS.active_file.fit_parameters["F0_ret"] = [fit_result[0][3], fit_result[1][3]]

				GLOBALVARS.active_file.has_fit_r = True

		replot_canvas()

def mark_baseline() -> None:
	if GLOBALVARS.active_file.baseline == False:
		GLOBALVARS.active_file.baseline = True
		if len(listbox_all_selected_files.curselection()) > 0:
			index_highlighted_file_name: list[int] = listbox_all_selected_files.curselection()
			listbox_all_selected_files.itemconfig(index_highlighted_file_name, fg = "red")	
	else:
		GLOBALVARS.active_file.baseline = False
		if len(listbox_all_selected_files.curselection()) > 0:
			index_highlighted_file_name: list[int] = listbox_all_selected_files.curselection()
			listbox_all_selected_files.itemconfig(index_highlighted_file_name, fg = "black")	

def mark_reference() -> None:
	if GLOBALVARS.active_file.reference ==  False:
		GLOBALVARS.active_file.reference = True
		if len(listbox_all_selected_files.curselection()) > 0:
			index_highlighted_file_name: list[int] = listbox_all_selected_files.curselection()
			listbox_all_selected_files.itemconfig(index_highlighted_file_name, fg = "green")	
	else:
		GLOBALVARS.active_file.reference = False
		if len(listbox_all_selected_files.curselection()) > 0:
			index_highlighted_file_name: list[int] = listbox_all_selected_files.curselection()
			listbox_all_selected_files.itemconfig(index_highlighted_file_name, fg = "black")	

def baseline_correction():
	if len(GLOBALVARS.baseline_curve["Distance"]) > 1:
		for curve in GLOBALVARS.selected_files:
			curve.subtract_baseline(GLOBALVARS.baseline_curve["Distance"], GLOBALVARS.baseline_curve["Force"])
			recalculate_derivatives(curve)

	replot_canvas()

def load_baseline() -> None:
	file_path = filedialog.askopenfilename(title="Select Baseline CSV", filetypes=[("CSV", ('*.csv')), ("All files", "*.*")])	
	GLOBALVARS.baseline_curve = pd.read_csv(file_path)

def save_baseline() -> None:
	if len(GLOBALVARS.baseline_curve["Distance"]) > 1:
		GLOBALVARS.baseline_curve.to_csv(os.path.join(GLOBALVARS.output_directory, "BASELINE.csv"), index=False)

def calculate_baseline() -> None:
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
	GLOBALVARS.baseline_curve = pd.DataFrame({"Force": baseline_y, "Distance": baseline_x[:-1]})	

def view_baseline() -> None:
	if len(GLOBALVARS.baseline_curve["Distance"]) > 1:
		plt.plot(GLOBALVARS.baseline_curve["Distance"], GLOBALVARS.baseline_curve["Force"], )
		plt.xlabel("Distance (\u03bcm)")
		plt.ylabel("Force (pN)")
		plt.show()
		plt.close()

# if the curve is not nicely trimmed already, attempt to do it automatically!
# NOTE: may not be correct for unusual FD curves. It is always best to trim manually
def auto_split_and_trim(curve) -> None:
	inflection_point = 0
	max_force = 0
	for i in range(len(curve.processed_dataframe["Processed_Time"])-1):
		if curve.processed_dataframe["Processed_Force"][i] > max_force and curve.processed_dataframe["Processed_Force"][i] > curve.processed_dataframe["Processed_Force"][i+1]:	
			max_force = curve.processed_dataframe["Processed_Force"][i]
			inflection_point = i

	force_ext = []
	dist_ext= []
	time_ext = []
	first_deriv_ext = []
	force_ret = []
	dist_ret = []
	time_ret = []
	first_deriv_ret = []

	for i in range(len(curve.processed_dataframe["Processed_Force"][0:inflection_point])):
		force_ext.append(curve.processed_dataframe["Processed_Force"][i])
		time_ext.append(curve.processed_dataframe["Processed_Time"][i])
		dist_ext.append(curve.processed_dataframe["Processed_Distance"][i])
		first_deriv_ext.append(curve.first_derivative_dataframe["First_Derivative"][i])

	for i in range(len(curve.processed_dataframe["Processed_Force"][inflection_point+1:])):
		force_ret.append(curve.processed_dataframe["Processed_Force"][inflection_point+i+1])
		time_ret.append(curve.processed_dataframe["Processed_Time"][inflection_point+i+1])
		dist_ret.append(curve.processed_dataframe["Processed_Distance"][inflection_point+i+1])
		first_deriv_ret.append(curve.first_derivative_dataframe["First_Derivative"][inflection_point+i+1])

	curve.dataframe_extension = pd.DataFrame({"Force_Extension": force_ext, "Distance_Extension": dist_ext, "Time_Extension": time_ext})
	curve.dataframe_retraction = pd.DataFrame({"Force_Retraction": force_ret, "Distance_Retraction": dist_ret, "Time_Retraction": time_ret})
	curve.trimmed_first_derivative_dataframe_extension = pd.DataFrame({"Trimmed_First_Derivative": first_deriv_ext, "Trimmed_Distance": dist_ext, "Trimmed_Time": time_ext})
	curve.trimmed_first_derivative_dataframe_retraction = pd.DataFrame({"Trimmed_First_Derivative": first_deriv_ret, "Trimmed_Distance": dist_ret, "Trimmed_Time": time_ret})

	curve.trimmed_e = True
	curve.trimmed_r = True

	curve.xmin_e[0] = curve.dataframe_extension["Distance_Extension"].iloc[0]
	curve.xmin_e[1] = curve.dataframe_extension["Time_Extension"].iloc[0]
	curve.xmin_r[0] = curve.dataframe_retraction["Distance_Retraction"].iloc[-1]
	curve.xmin_r[1] = curve.dataframe_retraction["Time_Retraction"].iloc[-1]

	curve.xmax_e[0] = curve.dataframe_extension["Distance_Extension"].iloc[-1]
	curve.xmax_e[1] = curve.dataframe_extension["Time_Extension"].iloc[-1]
	curve.xmax_r[0] = curve.dataframe_retraction["Distance_Retraction"].iloc[0]
	curve.xmax_r[1] = curve.dataframe_retraction["Time_Retraction"].iloc[0]

	curve.ymax_e = max_force
	curve.ymax_r = max_force

	# Now, from the extension curves, trim them to Fc for future fitting!
	# Use the peak of the first derivative
	fc_inflection_point_ext = 0
	fc_ext = 0
	max_first_deriv_ext = 0
	fc_inflection_point_ret = 0
	fc_ret = 0
	max_first_deriv_ret = 0

	for i in range(len(curve.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"])):
		if curve.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"][i] > max_first_deriv_ext and curve.dataframe_extension["Force_Extension"][i] < 80 and curve.dataframe_extension["Force_Extension"][i] > 4: # shouldn't be above 80pN or below 4pN anyway
			fc_inflection_point_ext = i
			max_first_deriv_ext = curve.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"][i]
			fc_ext = curve.dataframe_extension["Force_Extension"][i]
	for i in range(len(curve.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"])):
		if abs(curve.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"][i]) > max_first_deriv_ret and curve.dataframe_retraction["Force_Retraction"][i] < 80 and curve.dataframe_retraction["Force_Retraction"][i] > 4:
			fc_inflection_point_ret = i
			max_first_deriv_ret = abs(curve.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"][i])
			fc_ret = curve.dataframe_retraction["Force_Retraction"][i]
	
	curve.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": curve.dataframe_extension["Force_Extension"][0:fc_inflection_point_ext],
						      "Fit_Distance_Extension": curve.dataframe_extension["Distance_Extension"][0:fc_inflection_point_ext],
						      "Fit_Time_Extension": curve.dataframe_extension["Time_Extension"][0:fc_inflection_point_ext]})

	curve.fit_dataframe_retraction = pd.DataFrame({"Fit_Force_Retraction": curve.dataframe_retraction["Force_Retraction"][fc_inflection_point_ret : -1],
						      "Fit_Distance_Retraction": curve.dataframe_retraction["Distance_Retraction"][fc_inflection_point_ret : -1],
						      "Fit_Time_Retraction": curve.dataframe_retraction["Time_Retraction"][fc_inflection_point_ret : -1]})

	curve.trimmed_f_e = True
	curve.trimmed_f_r = True

	curve.xmin_f_e[0] = curve.fit_dataframe_extension["Fit_Distance_Extension"].iloc[0]
	curve.xmin_f_e[1] = curve.fit_dataframe_extension["Fit_Time_Extension"].iloc[0]
	curve.xmin_f_r[0] = curve.fit_dataframe_retraction["Fit_Distance_Retraction"].iloc[-1]
	curve.xmin_f_r[1] = curve.fit_dataframe_retraction["Fit_Time_Retraction"].iloc[-1]

	curve.xmax_f_e[0] = curve.fit_dataframe_extension["Fit_Distance_Extension"].iloc[-1]
	curve.xmax_f_e[1] = curve.fit_dataframe_extension["Fit_Time_Extension"].iloc[-1]
	curve.xmax_f_r[0] = curve.fit_dataframe_retraction["Fit_Distance_Retraction"].iloc[0]
	curve.xmax_f_r[1] = curve.fit_dataframe_retraction["Fit_Time_Retraction"].iloc[0]

	curve.ymax_f_e = fc_ext
	curve.ymax_f_r = fc_ret

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
		if curve.trimmed_e == False and curve.trimmed_r == False:
			if curve in GLOBALVARS.selected_files:	
				auto_split_and_trim(curve)

		if curve.has_fit_e == False and curve.trimmed_e == True and len(curve.dataframe_extension["Distance_Extension"]) > 30: # if the curve is not already fit, then fit the trimmed data.
			fc_inflection_point_ext = 0
			fc_ext = 0
			max_first_deriv_ext = 0

			for i in range(len(curve.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"])):
				if curve.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"][i] > max_first_deriv_ext and curve.dataframe_extension["Force_Extension"][i] < 80 and curve.dataframe_extension["Force_Extension"][i] > 4: # shouldn't be above 80pN or below 4pN anyway
					fc_inflection_point_ext = i
					max_first_deriv_ext = curve.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"][i]
					fc_ext = curve.dataframe_extension["Force_Extension"][i]

			curve.fit_dataframe_extension = pd.DataFrame({"Fit_Force_Extension": curve.dataframe_extension["Force_Extension"][0:fc_inflection_point_ext],
								      "Fit_Distance_Extension": curve.dataframe_extension["Distance_Extension"][0:fc_inflection_point_ext],
								      "Fit_Time_Extension": curve.dataframe_extension["Time_Extension"][0:fc_inflection_point_ext]})

			fit_result = fit_eOdijk_F0(curve.fit_dataframe_extension["Fit_Distance_Extension"], curve.fit_dataframe_extension["Fit_Force_Extension"])
			curve.precalculated_fit_extension = pd.DataFrame({"Precalculated_Fit_Force_Extension": fit_result[2], "Precalculated_Fit_Distance_Extension": curve.fit_dataframe_extension["Fit_Distance_Extension"], "Precalculated_Fit_Time_Extension": curve.fit_dataframe_extension["Fit_Time_Extension"]})

			curve.fit_parameters["Lp_ext"] = [fit_result[0][0], fit_result[1][0]]
			curve.fit_parameters["Lc_ext"] = [fit_result[0][1], fit_result[1][1]]
			curve.fit_parameters["S_ext"] = [fit_result[0][2], fit_result[1][2]]
			curve.fit_parameters["F0_ext"] = [fit_result[0][3], fit_result[1][3]]

			curve.xmin_f_e[0] = curve.fit_dataframe_extension["Fit_Distance_Extension"].iloc[0]
			curve.xmin_f_e[1] = curve.fit_dataframe_extension["Fit_Time_Extension"].iloc[0]

			curve.xmax_f_e[0] = curve.fit_dataframe_extension["Fit_Distance_Extension"].iloc[-1]
			curve.xmax_f_e[1] = curve.fit_dataframe_extension["Fit_Time_Extension"].iloc[-1]

			curve.ymax_f_e = fc_ext

			curve.has_fit_e = True
			curve.trimmed_f_e = True

		if curve.has_fit_r == False and curve.trimmed_r == True and len(curve.dataframe_retraction["Distance_Retraction"]) > 30: # if the curve is not already fit, then fit the trimmed data.

			fc_inflection_point_ret = 0
			fc_ret = 0
			max_first_deriv_ret = 0

			for i in range(len(curve.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"])):
				if abs(curve.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"][i]) > max_first_deriv_ret and curve.dataframe_retraction["Force_Retraction"][i] < 80 and curve.dataframe_retraction["Force_Retraction"][i] > 4:
					fc_inflection_point_ret = i
					max_first_deriv_ret = abs(curve.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"][i])
					fc_ret = curve.dataframe_retraction["Force_Retraction"][i]

			curve.fit_dataframe_retraction = pd.DataFrame({"Fit_Force_Retraction": curve.dataframe_retraction["Force_Retraction"][fc_inflection_point_ret : -1],
								      "Fit_Distance_Retraction": curve.dataframe_retraction["Distance_Retraction"][fc_inflection_point_ret : -1],
								      "Fit_Time_Retraction": curve.dataframe_retraction["Time_Retraction"][fc_inflection_point_ret : -1]})

			fit_result = fit_eOdijk_F0(curve.fit_dataframe_retraction["Fit_Distance_Retraction"], curve.fit_dataframe_retraction["Fit_Force_Retraction"])
			curve.precalculated_fit_retraction = pd.DataFrame({"Precalculated_Fit_Force_Retraction": fit_result[2], "Precalculated_Fit_Distance_Retraction": curve.fit_dataframe_retraction["Fit_Distance_Retraction"], "Precalculated_Fit_Time_Retraction": curve.fit_dataframe_retraction["Fit_Time_Retraction"]})

			curve.fit_parameters["Lp_ret"] = [fit_result[0][0], fit_result[1][0]]
			curve.fit_parameters["Lc_ret"] = [fit_result[0][1], fit_result[1][1]]
			curve.fit_parameters["S_ret"] = [fit_result[0][2], fit_result[1][2]]
			curve.fit_parameters["F0_ret"] = [fit_result[0][3], fit_result[1][3]]

			curve.xmin_f_r[0] = curve.fit_dataframe_retraction["Fit_Distance_Retraction"].iloc[-1]
			curve.xmin_f_r[1] = curve.fit_dataframe_retraction["Fit_Time_Retraction"].iloc[-1]
			curve.xmax_f_r[0] = curve.fit_dataframe_retraction["Fit_Distance_Retraction"].iloc[0]
			curve.xmax_f_r[1] = curve.fit_dataframe_retraction["Fit_Time_Retraction"].iloc[0]
			curve.ymax_f_r = fc_ret

			curve.has_fit_r = True
			curve.trimmed_f_r = True

		if curve.has_fit_e == True:
			contour_lengths.append(curve.fit_parameters["Lc_ext"][0])
		if curve.has_fit_r == True:
			contour_lengths.append(curve.fit_parameters["Lc_ret"][0])

	mean_Lc = np.mean(contour_lengths)
	force_corrections = []
	critical_Lc = 1.342424 # lambda Lc=16.5um. 22.15um = 1.342424x Lc. 1.342424 LC = 110 pN.
	
	for curve in reference_curves:

		if curve.trimmed_e == True and len(curve.dataframe_extension["Distance_Extension"]) > 30:
			Lc_space_ext = curve.dataframe_extension["Distance_Extension"] / mean_Lc # Normalise to Lc space
			for i in range(len(Lc_space_ext)-1):
				if Lc_space_ext.iloc[i] <= critical_Lc and Lc_space_ext.iloc[i+1] > critical_Lc:
					# Take an average of 5 surrounding data points
					force_at_critical_lc_ext = (curve.dataframe_extension["Force_Extension"].iloc[i] + curve.dataframe_extension["Force_Extension"].iloc[i+1]+curve.dataframe_extension["Force_Extension"].iloc[i+2]+curve.dataframe_extension["Force_Extension"].iloc[i-1]+curve.dataframe_extension["Force_Extension"].iloc[i-2] )/5

					force_correction_factor_ext = 110/force_at_critical_lc_ext
					force_corrections.append(force_correction_factor_ext)

		if curve.trimmed_r == True and len(curve.dataframe_retraction["Distance_Retraction"]) > 30:
			Lc_space_ret = curve.dataframe_retraction["Distance_Retraction"] / mean_Lc # Normalise to Lc space
			for i in range(len(Lc_space_ret)-1):
				if Lc_space_ret.iloc[i] <= critical_Lc and Lc_space_ret.iloc[i+1] > critical_Lc:
					# Take an average of 5 surrounding data points
					force_at_critical_lc_ret = (curve.dataframe_retraction["Force_Retraction"].iloc[i] + curve.dataframe_retraction["Force_Retraction"].iloc[i+1]+curve.dataframe_retraction["Force_Retraction"].iloc[i+2]+curve.dataframe_retraction["Force_Retraction"].iloc[i-1]+curve.dataframe_retraction["Force_Retraction"].iloc[i-2] )/5

					force_correction_factor_ret = 110/force_at_critical_lc_ret
					force_corrections.append(force_correction_factor_ret)
	if len(force_corrections) > 0:
		mean_force_correction = np.mean(force_corrections)
	else:
		mean_force_correction = 1

	# apply the force correction to every selected curve.
	for curve in GLOBALVARS.selected_files:
		curve.processed_dataframe["Processed_Force"] = curve.processed_dataframe["Processed_Force"] * mean_force_correction
		if curve.trimmed_e == True and len(curve.dataframe_extension["Distance_Extension"]) > 30:
			curve.dataframe_extension["Force_Extension"] = curve.dataframe_extension["Force_Extension"] * mean_force_correction
		if curve.trimmed_r == True and len(curve.dataframe_retraction["Distance_Retraction"]) > 30:
			curve.dataframe_retraction["Force_Retraction"] = curve.dataframe_retraction["Force_Retraction"] * mean_force_correction
		if curve.trimmed_e == True and curve.trimmed_f_e == True and len(curve.fit_dataframe_extension["Fit_Distance_Extension"]) > 30:
			curve.fit_dataframe_extension["Fit_Force_Extension"] = curve.fit_dataframe_extension["Fit_Force_Extension"] * mean_force_correction
		if curve.trimmed_r == True and curve.trimmed_f_r == True and len(curve.fit_dataframe_retraction["Fit_Distance_Retraction"]) > 30:
			curve.fit_dataframe_retraction["Fit_Force_Retraction"] = curve.fit_dataframe_retraction["Fit_Force_Retraction"] * mean_force_correction
		if curve.has_fit_e == True:
			curve.precalculated_fit_extension["Precalculated_Fit_Force_Extension"] = curve.precalculated_fit_extension["Precalculated_Fit_Force_Extension"] * mean_force_correction
		if curve.has_fit_r == True:
			curve.precalculated_fit_retraction["Precalculated_Fit_Force_Retraction"] = curve.precalculated_fit_retraction["Precalculated_Fit_Force_Retraction"] * mean_force_correction
		curve.first_derivative_dataframe["First_Derivative"] = curve.first_derivative_dataframe["First_Derivative"] * mean_force_correction
		curve.second_derivative_dataframe["Second_Derivative"] = curve.second_derivative_dataframe["Second_Derivative"] * mean_force_correction
		curve.is_force_scaled = True
	
	update_trim_entries_ui()
	replot_canvas()
			

def expand_graph() -> None:
	replot_canvas(True)

def export_data() -> None:
	fit_params: dict = {"File": [], "Lp-e": [], "Lc-e": [], "S-e": [], "F0-e": [], "Lp-r": [], "Lc-r": [], "S-r": [], "F0-r": [], "Fc-e": [], "Fc-r": [], "sigma_e" : []}
	for file in GLOBALVARS.selected_files:
		output_data = pd.concat([file.dataframe, file.processed_dataframe, file.dataframe_extension, file.dataframe_retraction, file.first_derivative_dataframe, file.second_derivative_dataframe, file.fit_dataframe_extension, file.fit_dataframe_retraction, file.precalculated_fit_extension, file.precalculated_fit_retraction, file.fit_parameters, pd.DataFrame({"Fc_e": [file.fc_e]}), pd.DataFrame({"Fc_r": [file.fc_r]}), pd.DataFrame({"sigma_e": [file.sigma_e]}), pd.DataFrame({"is_baseline": [file.baseline]}), pd.DataFrame({"is_reference": [file.reference]}), pd.DataFrame({"is_baseline_subtracted": [file.is_baseline_subtracted]}), pd.DataFrame({"is_force_scaled": [file.is_force_scaled]}),pd.DataFrame({"x_variable": [x_variable_combo.get()]}) ,pd.DataFrame({"y_variable": [y_variable_combo.get()]})], axis = 1)
		output_data.to_csv(os.path.join(GLOBALVARS.output_directory, file.name+".csv"), index=False)

		if file.has_fit_e == True:

			fit_params["File"].append(file.name)
			fit_params["Lp-e"].append(file.fit_parameters["Lp_ext"][0])
			fit_params["Lc-e"].append(file.fit_parameters["Lc_ext"][0])
			fit_params["S-e"].append(file.fit_parameters["S_ext"][0])
			fit_params["F0-e"].append(file.fit_parameters["F0_ext"][0])
			fit_params["Lp-r"].append(file.fit_parameters["Lp_ret"][0])
			fit_params["Lc-r"].append(file.fit_parameters["Lc_ret"][0])
			fit_params["S-r"].append(file.fit_parameters["S_ret"][0])
			fit_params["F0-r"].append(file.fit_parameters["F0_ret"][0])
			fit_params["Fc-e"].append(file.fc_e)
			fit_params["Fc-r"].append(file.fc_r)
			fit_params["sigma_e"].append(file.sigma_e)

	out = pd.DataFrame(fit_params)
	out.to_csv(os.path.join(GLOBALVARS.output_directory,"FIT_PARAMETERS.csv"), index=False)

def auto_fc() -> None:
	if GLOBALVARS.active_file != None:
		if  variable_radio_buttons_view.get() == "full":
			if variable_radio_buttons.get() == "extension" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					xmin = GLOBALVARS.active_file.xmin_e[1]
					xmax = GLOBALVARS.active_file.xmax_e[1]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Time"]
				else:
					xmin = GLOBALVARS.active_file.xmin_e[0]
					xmax = GLOBALVARS.active_file.xmax_e[0]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]
			elif variable_radio_buttons.get() == "retraction" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					xmin = GLOBALVARS.active_file.xmin_r[1]
					xmax = GLOBALVARS.active_file.xmax_r[1]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Time"]
				else:
					xmin = GLOBALVARS.active_file.xmin_r[0]
					xmax = GLOBALVARS.active_file.xmax_r[0]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]
		
			force_data = GLOBALVARS.active_file.processed_dataframe["Processed_Force"]
			first_deriv_force_data = GLOBALVARS.active_file.first_derivative_dataframe["First_Derivative"]

		elif variable_radio_buttons_view.get() == "extension" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_e == True:
			if GLOBALVARS.active_file.plot_time == True:
				xmin = GLOBALVARS.active_file.xmin_f_e[1]
				xmax = GLOBALVARS.active_file.xmax_f_e[1]
				x_data = GLOBALVARS.active_file.dataframe_extension["Time_Extension"]
			else:
				xmin = GLOBALVARS.active_file.xmin_f_e[0]
				xmax = GLOBALVARS.active_file.xmax_f_e[0]
				x_data = GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]
			force_data = GLOBALVARS.active_file.dataframe_extension["Force_Extension"]
			first_deriv_force_data = GLOBALVARS.active_file.trimmed_first_derivative_dataframe_extension["Trimmed_First_Derivative"]
		elif variable_radio_buttons_view.get() == "retraction" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_r == True:
			if GLOBALVARS.active_file.plot_time == True:
				xmin = GLOBALVARS.active_file.xmin_f_r[1]
				xmax = GLOBALVARS.active_file.xmax_f_r[1]
				x_data = GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]
			else:
				xmin = GLOBALVARS.active_file.xmin_f_r[0]
				xmax = GLOBALVARS.active_file.xmax_f_r[0]
				x_data = GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]
			force_data = GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"]
			first_deriv_force_data = GLOBALVARS.active_file.trimmed_first_derivative_dataframe_retraction["Trimmed_First_Derivative"]

		ixmin = 0
		ixmax = len(x_data)
		for i in range(len(x_data)-1):
			if x_data[i] < xmin and x_data[i+1] >= xmin:
				ixmin = i
			if x_data[i] < xmax and x_data[i+1] >= xmax :
				ixmax = i
		fc = 0
		ifc = 0
		for i in range(ixmin, ixmax):
			if abs(first_deriv_force_data[i]) >= fc:
				fc = abs(first_deriv_force_data[i])
				ifc = i

		ymax = force_data[ifc]

		if first_deriv_force_data[ifc] > 0:
			xmax = x_data[ifc]
			GLOBALVARS.active_file.fc_e = force_data[ifc]
		else:
			xmin = x_data[ifc]
			GLOBALVARS.active_file.fc_r = force_data[ifc]

		
		if  variable_radio_buttons_view.get() == "full":
			if variable_radio_buttons.get() == "extension" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					GLOBALVARS.active_file.xmin_e[1] = xmin
					GLOBALVARS.active_file.xmax_e[1] = xmax
				else:
					GLOBALVARS.active_file.xmin_e[0] = xmin
					GLOBALVARS.active_file.xmax_e[0] = xmax
				GLOBALVARS.active_file.ymax_e = ymax
			elif variable_radio_buttons.get() == "retraction" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					GLOBALVARS.active_file.xmin_r[1] = xmin
					GLOBALVARS.active_file.xmax_r[1] = xmax
				else:
					GLOBALVARS.active_file.xmin_r[0] = xmin
					GLOBALVARS.active_file.xmax_r[0] = xmax
				GLOBALVARS.active_file.ymax_r = ymax
		
		elif variable_radio_buttons_view.get() == "extension" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_e == True:
			if GLOBALVARS.active_file.plot_time == True:
				GLOBALVARS.active_file.xmin_f_e[1] = xmin
				GLOBALVARS.active_file.xmax_f_e[1] = xmax
			else:
				GLOBALVARS.active_file.xmin_f_e[0] = xmin
				GLOBALVARS.active_file.xmax_f_e[0] = xmax
			GLOBALVARS.active_file.ymax_f_e = ymax

		elif variable_radio_buttons_view.get() == "retraction" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_r == True:
			if GLOBALVARS.active_file.plot_time == True:
				GLOBALVARS.active_file.xmin_f_r[1] = xmin
				GLOBALVARS.active_file.xmax_f_r[1] = xmax
			else:
				GLOBALVARS.active_file.xmin_f_r[0] = xmin
				GLOBALVARS.active_file.xmax_f_r[0] = xmax
			GLOBALVARS.active_file.ymax_f_r = ymax
		update_trim_entries_ui()
		replot_canvas()


def options_popup() -> None:

	def update_options_settings() -> None:
		GLOBALVARS.frame_rate = float(entry_popup_frame_rate.get())
		GLOBALVARS.extension_speed_um_s = float(entry_popup_extension_speed.get())

	popup_window = tk.Toplevel()
	popup_window.rowconfigure(0,weight=1)
	popup_window.columnconfigure(0,weight=1)

	popup_frame = tk.Frame(master=popup_window)
	popup_frame.grid(row=0, column=0)

	label_popup_framerate = tk.Label(master = popup_frame, text="Frame Rate:")
	label_popup_extension_speed = tk.Label(master = popup_frame, text="Extension Speed:")
	label_popup_extension_speed_units = tk.Label(master= popup_frame, text="\u03bcm/s")
	label_popup_framerate_units = tk.Label(master= popup_frame, text="Hz")
	entry_popup_frame_rate = tk.Entry(master= popup_frame, width=3)
	entry_popup_extension_speed = tk.Entry(master= popup_frame, width=4)
	popup_x_variable_combo = ttk.Combobox(master= popup_frame, state="readonly")
	popup_y_variable_combo = ttk.Combobox(master= popup_frame, state="readonly")


	label_popup_framerate.grid(row=0, column=0, sticky=tk.E)
	label_popup_extension_speed.grid(row=1, column=0, sticky=tk.E)
	label_popup_extension_speed_units.grid(row=1, column=2, sticky=tk.W)
	label_popup_framerate_units.grid(row=0, column=2, sticky=tk.W)

	# Set key optical settings used in the session for savgol filter smoothing
	entry_popup_frame_rate.insert(0, str(GLOBALVARS.frame_rate))
	entry_popup_frame_rate.grid(row=0, column=1)
	entry_popup_extension_speed.insert(0, str(GLOBALVARS.extension_speed_um_s))
	entry_popup_extension_speed.grid(row=1, column=1)

	# Combobox for selecting whether to use Force 2x / Trap 2 etc
	popup_x_variable_combo["values"] = ["Distance 1", "Distance 2"]
	popup_x_variable_combo.current(0)
	popup_x_variable_combo.grid(row=2, column=0, pady=10)

	popup_y_variable_combo["values"] = ["Force 2x", "Force 2y", "Trap 2"]
	popup_y_variable_combo.current(2)
	popup_y_variable_combo.grid(row=2, column=2, pady=10)

	popup_update_button = tk.Button(master=popup_frame,text="Update", command=update_options_settings)
	popup_update_button.grid(row=3, column=1, pady=10)


def recalculate_derivatives(curve=None) -> None:
	if GLOBALVARS.active_file != None and curve==None:
		# Compute derivatives of the processed data, useful for data trimming
		derivatives: list = GLOBALVARS.active_file.differentiate_savgol(GLOBALVARS.active_file.processed_dataframe["Processed_Time"], GLOBALVARS.active_file.processed_dataframe["Processed_Force"], 0.75*GLOBALVARS.frame_rate/GLOBALVARS.extension_speed_um_s, 2)
		first_derivative = [derivatives[0], derivatives[1]]
		second_derivative = [derivatives[0], derivatives[2]]

		GLOBALVARS.active_file.first_derivative_dataframe = pd.DataFrame({"First_Derivative": first_derivative[1], "Time": first_derivative[0], "Distance": GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]})
		GLOBALVARS.active_file.second_derivative_dataframe = pd.DataFrame({"Second_Derivative": second_derivative[1], "Time": second_derivative[0], "Distance": GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]})

		replot_canvas()

	elif curve != None:
		derivatives: list = curve.differentiate_savgol(curve.processed_dataframe["Processed_Time"], curve.processed_dataframe["Processed_Force"], 0.75*GLOBALVARS.frame_rate/GLOBALVARS.extension_speed_um_s, 2)
		first_derivative = [derivatives[0], derivatives[1]]
		second_derivative = [derivatives[0], derivatives[2]]

		curve.first_derivative_dataframe = pd.DataFrame({"First_Derivative": first_derivative[1], "Time": first_derivative[0], "Distance": curve.processed_dataframe["Processed_Distance"]})
		curve.second_derivative_dataframe = pd.DataFrame({"Second_Derivative": second_derivative[1], "Time": second_derivative[0], "Distance": curve.processed_dataframe["Processed_Distance"]})

def supercoiling_density_estimation() -> None:

	# Length at 70pN:
	# f(x) = -1.0414x + 1.1
	# R2 = 0.9975
	# X shifted to start at 0,0, forced intercept at 0,0, equation shifted back
	def formula(relative_length_via_lc) -> float:
		return -1.04141820566195*relative_length_via_lc+1.10074780083851

	reference_curves = []
	contour_lengths = []	
	
	for curve in GLOBALVARS.selected_files:
		if curve.reference == True:
			reference_curves.append(curve)

	for curve in reference_curves:
		if curve.trimmed_e == False and curve.trimmed_r == False:
			if curve in GLOBALVARS.selected_files:	
				auto_split_and_trim(curve)

		if curve.trimmed_e == True and curve.has_fit_e == False and len(curve.dataframe_extension["Distance_Extension"]) > 30: # if the curve is not already fit, then fit the trimmed data.
			fit_result = fit_eOdijk_F0(curve.fit_dataframe_extension["Fit_Distance_Extension"], curve.fit_dataframe_extension["Fit_Force_Extension"])
			curve.precalculated_fit_extension = pd.DataFrame({"Precalculated_Fit_Force_Extension": fit_result[2], "Precalculated_Fit_Distance_Extension": curve.fit_dataframe_extension["Fit_Distance_Extension"], "Precalculated_Fit_Time_Extension": curve.fit_dataframe_extension["Fit_Time_Extension"]})

			curve.fit_parameters["Lp_ext"] = [fit_result[0][0], fit_result[1][0]]
			curve.fit_parameters["Lc_ext"] = [fit_result[0][1], fit_result[1][1]]
			curve.fit_parameters["S_ext"] = [fit_result[0][2], fit_result[1][2]]
			curve.fit_parameters["F0_ext"] = [fit_result[0][3], fit_result[1][3]]

			curve.has_fit_e = True

		if curve.trimmed_r == True and curve.has_fit_r == False and len(curve.dataframe_retraction["Distance_Retraction"]) > 30: # if the curve is not already fit, then fit the trimmed data.
			fit_result = fit_eOdijk_F0(curve.fit_dataframe_retraction["Fit_Distance_Retraction"], curve.fit_dataframe_retraction["Fit_Force_Retraction"])
			curve.precalculated_fit_retraction = pd.DataFrame({"Precalculated_Fit_Force_Retraction": fit_result[2], "Precalculated_Fit_Distance_Retraction": curve.fit_dataframe_retraction["Fit_Distance_Retraction"], "Precalculated_Fit_Time_Retraction": curve.fit_dataframe_retraction["Fit_Time_Retraction"]})

			curve.fit_parameters["Lp_ret"] = [fit_result[0][0], fit_result[1][0]]
			curve.fit_parameters["Lc_ret"] = [fit_result[0][1], fit_result[1][1]]
			curve.fit_parameters["S_ret"] = [fit_result[0][2], fit_result[1][2]]
			curve.fit_parameters["F0_ret"] = [fit_result[0][3], fit_result[1][3]]

			curve.has_fit_r = True

		if curve.has_fit_e == True:
			contour_lengths.append(curve.fit_parameters["Lc_ext"][0])
		if curve.has_fit_r == True:
			contour_lengths.append(curve.fit_parameters["Lc_ret"][0])

	mean_Lc = np.mean(contour_lengths)

	ref_sigma_values = []

	for curve in reference_curves:

		if curve.trimmed_e == True and curve.dataframe_extension["Force_Extension"].iloc[-1] >= 70:
			force_index = 0
			for i in range(len(curve.dataframe_extension["Force_Extension"])):
				if curve.dataframe_extension["Force_Extension"].iloc[i] >= 70 and curve.dataframe_extension["Force_Extension"].iloc[i-1] < 70:
					force_index = i
			ref_sigma_values.append(formula((curve.dataframe_extension["Distance_Extension"][force_index] / mean_Lc)))

		if curve.trimmed_r == True and curve.dataframe_retraction["Force_Retraction"].iloc[0] >= 70:
			force_index = 0
			for i in range(len(curve.dataframe_retraction["Force_Retraction"])):
				if curve.dataframe_retraction["Force_Retraction"].iloc[i] <= 70 and curve.dataframe_retraction["Force_Retraction"].iloc[i-1] > 70:
					force_index = i
			ref_sigma_values.append(formula((curve.dataframe_retraction["Distance_Retraction"].iloc[force_index] / mean_Lc)))

	offset_factor = np.mean(ref_sigma_values) # Offset, sigma 0 should give 0

	for curve in GLOBALVARS.selected_files:
		if curve.trimmed_e == True and curve.dataframe_extension["Force_Extension"].iloc[-1] >= 70:
			force_index = 0
			for i in range(len(curve.dataframe_extension["Force_Extension"])):
				if curve.dataframe_extension["Force_Extension"].iloc[i] >= 70 and curve.dataframe_extension["Force_Extension"].iloc[i-1] < 70:
					force_index = i
			distance_at_70 = curve.dataframe_extension["Distance_Extension"].iloc[force_index]
			curve.sigma_e = formula((distance_at_70 / mean_Lc))-offset_factor

		if curve.trimmed_r == True and curve.dataframe_retraction["Force_Retraction"].iloc[0] >= 70:
			force_index = 0
			for i in range(len(curve.dataframe_retraction["Force_Retraction"])):
				if curve.dataframe_retraction["Force_Retraction"].iloc[i] <= 70 and curve.dataframe_retraction["Force_Retraction"].iloc[i-1] > 70:
					force_index = i
			distance_at_70 = curve.dataframe_retraction["Distance_Retraction"].iloc[force_index]
			curve.sigma_r = formula((distance_at_70 / mean_Lc))-offset_factor

		# Otherwise, try to crudly and temporarily split
		if curve.trimmed_e == False and curve.trimmed_r == False: # Do not save this temporary trim
			inflection_point = 0
			max_force = 0
			for i in range(len(curve.processed_dataframe["Processed_Time"])-1):
				if curve.processed_dataframe["Processed_Force"].iloc[i] > max_force and curve.processed_dataframe["Processed_Force"].iloc[i] > curve.processed_dataframe["Processed_Force"].iloc[i+1]:	
					max_force = curve.processed_dataframe["Processed_Force"].iloc[i]
					inflection_point = i
			
			force_ext = curve.processed_dataframe["Processed_Force"][0:inflection_point]
			time_ext = curve.processed_dataframe["Processed_Time"][0:inflection_point]
			dist_ext = curve.processed_dataframe["Processed_Distance"][0:inflection_point]
			force_ret = curve.processed_dataframe["Processed_Force"][inflection_point:-1]
			time_ret = curve.processed_dataframe["Processed_Time"][inflection_point:-1]
			dist_ret = curve.processed_dataframe["Processed_Distance"][inflection_point:-1]

			if force_ext.iloc[-1] >= 70:
				for i in range(len(force_ext)):
					if force_ext.iloc[i] >= 70 and force_ext.iloc[i-1] < 70:
						force_index = i
				distance_at_70 = dist_ext.iloc[force_index]
				curve.sigma_e = formula((distance_at_70 / mean_Lc))-offset_factor

			if force_ret.iloc[0] >= 70:
				for i in range(len(force_ret)):
					if force_ret.iloc[i] <= 70 and force_ret.iloc[i-1] > 70:
						force_index = i
				distance_at_70 = dist_ret.iloc[force_index]
				curve.sigma_r = formula((distance_at_70 / mean_Lc))-offset_factor
	
	replot_canvas()	
		

def plot_distance_time() -> None:
	if GLOBALVARS.active_file != None:
		plt.xlabel("Time (s)")
		plt.ylabel("Distance (\u03bcm)")
		plt.plot(GLOBALVARS.active_file.processed_dataframe["Processed_Time"], GLOBALVARS.active_file.processed_dataframe["Processed_Distance"])
		plt.show()
		plt.close()

def find_peak() -> None:
	if GLOBALVARS.active_file != None:
		if  variable_radio_buttons_view.get() == "full":
			if variable_radio_buttons.get() == "extension" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					xmin = GLOBALVARS.active_file.xmin_e[1]
					xmax = GLOBALVARS.active_file.xmax_e[1]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Time"]
				else:
					xmin = GLOBALVARS.active_file.xmin_e[0]
					xmax = GLOBALVARS.active_file.xmax_e[0]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]
			elif variable_radio_buttons.get() == "retraction" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					xmin = GLOBALVARS.active_file.xmin_r[1]
					xmax = GLOBALVARS.active_file.xmax_r[1]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Time"]
				else:
					xmin = GLOBALVARS.active_file.xmin_r[0]
					xmax = GLOBALVARS.active_file.xmax_r[0]
					x_data = GLOBALVARS.active_file.processed_dataframe["Processed_Distance"]
		
			force_data = GLOBALVARS.active_file.processed_dataframe["Processed_Force"]

		elif variable_radio_buttons_view.get() == "extension" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_e == True:
			if GLOBALVARS.active_file.plot_time == True:
				xmin = GLOBALVARS.active_file.xmin_f_e[1]
				xmax = GLOBALVARS.active_file.xmax_f_e[1]
				x_data = GLOBALVARS.active_file.dataframe_extension["Time_Extension"]
			else:
				xmin = GLOBALVARS.active_file.xmin_f_e[0]
				xmax = GLOBALVARS.active_file.xmax_f_e[0]
				x_data = GLOBALVARS.active_file.dataframe_extension["Distance_Extension"]
			force_data = GLOBALVARS.active_file.dataframe_extension["Force_Extension"]
		elif variable_radio_buttons_view.get() == "retraction" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_r == True:
			if GLOBALVARS.active_file.plot_time == True:
				xmin = GLOBALVARS.active_file.xmin_f_r[1]
				xmax = GLOBALVARS.active_file.xmax_f_r[1]
				x_data = GLOBALVARS.active_file.dataframe_retraction["Time_Retraction"]
			else:
				xmin = GLOBALVARS.active_file.xmin_f_r[0]
				xmax = GLOBALVARS.active_file.xmax_f_r[0]
				x_data = GLOBALVARS.active_file.dataframe_retraction["Distance_Retraction"]
			force_data = GLOBALVARS.active_file.dataframe_retraction["Force_Retraction"]

		ixmin = 0
		ixmax = len(x_data)
		for i in range(len(x_data)-1):
			if x_data[i] < xmin and x_data[i+1] >= xmin:
				ixmin = i
			if x_data[i] < xmax and x_data[i+1] >= xmax:
				ixmax = i
		peak = 0
		ipeak = 0
		for i in range(ixmin, ixmax):
			if abs(force_data[i]) >= peak:
				peak = abs(force_data[i])
				ipeak = i

		ymax = force_data[ipeak]

		if (ixmax-ipeak) < (ixmax-ixmin)/2:
			xmax = x_data[ipeak]
		else:
			xmin = x_data[ipeak]

		
		if  variable_radio_buttons_view.get() == "full":
			if variable_radio_buttons.get() == "extension" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					GLOBALVARS.active_file.xmin_e[1] = xmin
					GLOBALVARS.active_file.xmax_e[1] = xmax
				else:
					GLOBALVARS.active_file.xmin_e[0] = xmin
					GLOBALVARS.active_file.xmax_e[0] = xmax
				GLOBALVARS.active_file.ymax_e = ymax
			elif variable_radio_buttons.get() == "retraction" and variable_checkbutton_set_fit.get() == False:
				if GLOBALVARS.active_file.plot_time == True:
					GLOBALVARS.active_file.xmin_r[1] = xmin
					GLOBALVARS.active_file.xmax_r[1] = xmax
				else:
					GLOBALVARS.active_file.xmin_r[0] = xmin
					GLOBALVARS.active_file.xmax_r[0] = xmax
				GLOBALVARS.active_file.ymax_r = ymax
		
		elif variable_radio_buttons_view.get() == "extension" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_e == True:
			if GLOBALVARS.active_file.plot_time == True:
				GLOBALVARS.active_file.xmin_f_e[1] = xmin
				GLOBALVARS.active_file.xmax_f_e[1] = xmax
			else:
				GLOBALVARS.active_file.xmin_f_e[0] = xmin
				GLOBALVARS.active_file.xmax_f_e[0] = xmax
			GLOBALVARS.active_file.ymax_f_e = ymax

		elif variable_radio_buttons_view.get() == "retraction" and variable_checkbutton_view_fit.get() == False and GLOBALVARS.active_file.trimmed_r == True:
			if GLOBALVARS.active_file.plot_time == True:
				GLOBALVARS.active_file.xmin_f_r[1] = xmin
				GLOBALVARS.active_file.xmax_f_r[1] = xmax
			else:
				GLOBALVARS.active_file.xmin_f_r[0] = xmin
				GLOBALVARS.active_file.xmax_f_r[0] = xmax
			GLOBALVARS.active_file.ymax_f_r = ymax
		update_trim_entries_ui()
		replot_canvas()

def plot_reference_curves_overlaid() -> None:
	for curve in GLOBALVARS.selected_files:
		if curve.reference == True:
			plt.plot(curve.processed_dataframe["Processed_Distance"], curve.processed_dataframe["Processed_Force"])
	plt.xlabel("Distance (\u03bcm)")
	plt.ylabel("Force (pN)")
	plt.show()
	plt.close()

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
variable_radio_buttons_view = tk.StringVar()
variable_checkbutton_view_fit = tk.BooleanVar()
variable_checkbutton_set_fit = tk.BooleanVar()

variable_lp_fix = tk.BooleanVar()
variable_lc_fix = tk.BooleanVar()
variable_s_fix = tk.BooleanVar()
variable_f0_fix = tk.BooleanVar()

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
selection_menu.add_command(label='Mark Curve as Baseline',command=mark_baseline)
selection_menu.add_command(label='Mark Curve as Reference',command=mark_reference)
selection_menu.add_command(label='Recalculate Derivatives',command=recalculate_derivatives)

# Create the selection options
calibration_menu = Menu(menubar)
calibration_menu.add_command(label='View Baseline',command=view_baseline)
calibration_menu.add_command(label='Load Baseline',command=load_baseline)
calibration_menu.add_command(label='Save Baseline',command=save_baseline)
calibration_menu.add_command(label='Calculate Baseline',command=calculate_baseline)
calibration_menu.add_command(label='Subtract Baseline',command=baseline_correction)
calibration_menu.add_command(label='Auto Force Scale',command=auto_force_scale)
calibration_menu.add_command(label='Calculate Supercoiling Density',command=supercoiling_density_estimation)

view_menu = Menu(menubar)
view_menu.add_command(label='Distance Vs Time',command=plot_distance_time)
view_menu.add_command(label='Overlay References',command=plot_reference_curves_overlaid)
view_menu.add_command(label='Toggle First Derivative',command=toggle_first_derivative)
view_menu.add_command(label='Toggle Second Derivative',command=toggle_second_derivative)

# Add the dropdowns to the menubar
menubar.add_cascade(label="File",menu=file_menu)
menubar.add_cascade(label="Selection",menu=selection_menu)
menubar.add_cascade(label="Calibration",menu=calibration_menu)
menubar.add_cascade(label="View",menu=view_menu)
menubar.add_command(label="Options",command=options_popup)

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

frame_graphing_windows = tk.Frame(master=window, borderwidth=1, relief="solid")
frame_graphing_windows.grid(row=1, column=1, sticky=[tk.N, tk.S, tk.E, tk.W], padx=2, ipadx=2)
frame_graphing_windows.columnconfigure(0, weight=0)
frame_graphing_windows.columnconfigure(1, weight=1)
frame_graphing_windows.rowconfigure(0, weight=1)
frame_graphing_windows.rowconfigure(1, weight=0)
frame_graphing_windows.rowconfigure(2, weight=0)
frame_graphing_windows.rowconfigure(3, weight=0)

frame_input_buttons = tk.Frame(master=window)
frame_input_buttons.grid(row=2, column=1)

frame_graph_settings = tk.Frame(master=frame_graphing_windows, borderwidth=2, relief="groove")
frame_graph_settings.grid(row=3, column=0, columnspan=2, sticky=tk.S, pady=4, ipady=4)
frame_graph_settings.columnconfigure(0, weight=0)
frame_graph_settings.columnconfigure(1, weight=1)
frame_graph_settings.columnconfigure(2, weight=1)
frame_graph_settings.columnconfigure(3, weight=1)
frame_graph_settings.columnconfigure(4, weight=1)
frame_graph_settings.columnconfigure(5, weight=1)
frame_graph_settings.columnconfigure(6, weight=1)
frame_graph_settings.columnconfigure(7, weight=1)
frame_graph_settings.columnconfigure(8, weight=1)
frame_graph_settings.columnconfigure(8, weight=0)
frame_graph_settings.rowconfigure(0, weight=1)


frame_session_initialisation = tk.Frame(master=window, borderwidth=1, relief = "flat")
label_initialisation_blurb = tk.Label(master = frame_session_initialisation, text="Please initialise key session variables:")
label_initialisation_framerate = tk.Label(master = frame_session_initialisation, text="Frame Rate:")
label_initialisation_extension_speed = tk.Label(master = frame_session_initialisation, text="Extension Speed:")
label_initialisation_extension_speed_units = tk.Label(master=frame_session_initialisation, text="\u03bcm/s")
label_initialisation_framerate_units = tk.Label(master=frame_session_initialisation, text="Hz")
entry_frame_rate = tk.Entry(master=frame_session_initialisation, width=3)
entry_extension_speed = tk.Entry(master=frame_session_initialisation, width=4)
x_variable_combo = ttk.Combobox(master=frame_session_initialisation, state="readonly")
y_variable_combo = ttk.Combobox(master=frame_session_initialisation, state="readonly")



title_label = tk.Label(master=frame_title_manager, text="DEMO")
title_label.pack()

# File selection listboxes

## Create listboxes
tk.Label(master=frame_file_managers, text="All H5 Files:").grid(row=0, column=0)
listbox_all_h5_files = tk.Listbox(master=frame_file_managers,selectmode=tk.EXTENDED)
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

# Canvas to display graphs
canvas_graph_display = tk.Canvas(master=frame_graphing_windows, bg="#856ff8")
canvas_graph_display.grid(row=0, column=0, columnspan=2, sticky=[tk.N, tk.E, tk.S, tk.W])

# FD-time cutoff scrollers
tk.Label(master=frame_graphing_windows, text="Xmin: ").grid(row=1, column=0, sticky=tk.SW, padx=4)
scale_select_min_time = Scale(master=frame_graphing_windows, orient=HORIZONTAL, resolution=0.01)
scale_select_min_time.grid(row=1, column=1, sticky=[tk.E, tk.W, tk.S])
tk.Label(master=frame_graphing_windows, text="Xmax: ").grid(row=2, column=0, sticky=tk.SW, padx=4)
scale_select_min_time.bind("<ButtonRelease-1>", slider_min_release)
scale_select_max_time = Scale(master=frame_graphing_windows, orient=HORIZONTAL, resolution=0.01)
scale_select_max_time.grid(row=2, column=1,sticky=[tk.E, tk.W, tk.S])
scale_select_max_time.bind("<ButtonRelease-1>", slider_max_release)

tk.Label(master=frame_graph_settings, text="Xmin:").grid(row=0, column=0, sticky=tk.E)
entry_xmin = tk.Entry(master=frame_graph_settings,width=6)
entry_xmin.insert(0,"0")
entry_xmin.grid(row=0, column=1, sticky=tk.W)
tk.Label(master=frame_graph_settings, text="Xmax:").grid(row=0, column=2, sticky=tk.E)
entry_xmax = tk.Entry(master=frame_graph_settings,width=6)
entry_xmax.insert(0,"0")
entry_xmax.grid(row=0,column=3, sticky=tk.W)


tk.Label(master=frame_graph_settings, text="Ymin:").grid(row=0, column=4, sticky=tk.E)
entry_ymin = tk.Entry(master=frame_graph_settings,width=6)
entry_ymin.insert(0,"-5")
entry_ymin.grid(row=0, column=5, sticky=tk.W)
tk.Label(master=frame_graph_settings, text="Ymax:").grid(row=0, column=6, sticky=tk.E)
entry_ymax = tk.Entry(master=frame_graph_settings,width=6)
entry_ymax.insert(0,"30")
entry_ymax.grid(row=0, column=7, sticky=tk.W)

tk.Button(master=frame_graph_settings, text="Fc", command=auto_fc).grid(row=0, column=8, padx=4)
tk.Button(master=frame_graph_settings, text="Peak", command=find_peak).grid(row=0, column=9, padx=4)
tk.Button(master=frame_graph_settings, text="Update", command=update_trim_settings).grid(row=0, column=10, padx=4)


radio_extension_button = tk.Radiobutton(master=frame_input_buttons, text="Set Extension", value = "extension", variable=variable_radio_buttons, command=update_trim_entries_ui)
radio_extension_button.grid(row=0, column=1)
radio_extension_button.select()
tk.Radiobutton(master=frame_input_buttons, text="Set Retraction", value= "retraction", variable=variable_radio_buttons, command=update_trim_entries_ui).grid(row=0, column = 2)
tk.Checkbutton(master=frame_input_buttons, text="Set Fit", variable=variable_checkbutton_set_fit, command=radio_button_select).grid(row=0, column=3)

# Buttons to set correction factors
button_1 = Button(master=frame_input_buttons, text="Toggle Time", command=toggle_time)
button_1.grid(row=1,column=1, pady=5)
#button_2 = Button(master=frame_input_buttons, text="Auto Trim", command=auto_trim_data)
#button_2.grid(row=0,column=1)
button_3 = Button(master=frame_input_buttons, text="Trim Data", command=manual_trim_data)
button_3.grid(row=1,column=2, pady=5)

Enlarge_Button = Button(master=frame_input_buttons, text="Expand Graph", command=expand_graph)
Enlarge_Button.grid(row=1, column=3, pady=5)

radio_view_button = tk.Radiobutton(master=frame_input_buttons, text="View Full", value = "full", variable=variable_radio_buttons_view, command=radio_button_select)
radio_view_button.grid(row=2, column=0)
radio_view_button.select()
tk.Radiobutton(master=frame_input_buttons, text="View Extension", value = "extension", variable=variable_radio_buttons_view, command=radio_button_select).grid(row=2, column=1)
tk.Radiobutton(master=frame_input_buttons, text="View Retraction", value= "retraction", variable=variable_radio_buttons_view, command=radio_button_select).grid(row=2, column = 2)
tk.Checkbutton(master=frame_input_buttons, text="View Fit", variable=variable_checkbutton_view_fit, command=radio_button_select).grid(row=2, column=3)


tk.Label(master=frame_input_buttons, text="Lp:").grid(row=3, column=0, sticky=tk.E, pady=(20,0))
Lp_entry = tk.Entry(master=frame_input_buttons,width=6)
Lp_entry.insert(0,"50")
Lp_entry.grid(row=3, column=1, pady=(20,0))
tk.Checkbutton(master=frame_input_buttons, text="Fix", variable=variable_lp_fix).grid(row=3, column=2, sticky=tk.W, pady=(20,0))
Lp_display = tk.Entry(master=frame_input_buttons, text="", width=6, state="readonly")
Lp_display.grid(row=3, column=3, sticky=tk.EW, pady=(20,0))
tk.Label(master=frame_input_buttons, text="Lc:").grid(row=4, column=0, sticky=tk.E)
Lc_entry = tk.Entry(master=frame_input_buttons,width=6)
Lc_entry.insert(0,"16.5")
Lc_entry.grid(row=4, column=1)
tk.Checkbutton(master=frame_input_buttons, text="Fix", variable=variable_lc_fix).grid(row=4, column=2, sticky=tk.W)
Lc_display = tk.Entry(master=frame_input_buttons, text="", width=6, state="readonly")
Lc_display.grid(row=4, column=3, sticky=tk.EW)
tk.Label(master=frame_input_buttons, text="S:").grid(row=5, column=0, sticky=tk.E)
S_entry = tk.Entry(master=frame_input_buttons,width=6)
S_entry.insert(0,"1500")
S_entry.grid(row=5, column=1)
tk.Checkbutton(master=frame_input_buttons, text="Fix", variable=variable_s_fix).grid(row=5, column=2, sticky=tk.W)
S_display = tk.Entry(master=frame_input_buttons, text="", width=6, state="readonly")
S_display.grid(row=5, column=3, sticky=tk.EW)
tk.Label(master=frame_input_buttons, text="F0:").grid(row=6, column=0, sticky=tk.E)
F0_entry = tk.Entry(master=frame_input_buttons,width=6)
F0_entry.insert(0,"0")
F0_entry.grid(row=6, column=1)
tk.Checkbutton(master=frame_input_buttons, text="Fix", variable=variable_f0_fix).grid(row=6, column=2, sticky=tk.W)
F0_display = tk.Entry(master=frame_input_buttons, text="", state="readonly",width=6)
F0_display.grid(row=6, column=3, sticky=tk.EW)

Fit_Button = Button(master=frame_input_buttons, text="Fit", command=fit)
Fit_Button.grid(row=9, column=1, pady=(20,0))

tk.Label(master=frame_input_buttons, text="Sigma:").grid(row=10, column=2, sticky=tk.E)
sigma_display = tk.Entry(master=frame_input_buttons, text="", state="readonly",width=6)
sigma_display.grid(row=10 , column=3, sticky=tk.EW)

canvas_graph_display.bind("<Configure>", window_resize)

on_startup()
window.mainloop()
