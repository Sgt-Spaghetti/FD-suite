# Dependencies
The following python packages should be installed on the host system (or virtual environment):

- tkinter
- numpy
- matplotlib
- scipy
- pandas
- h5py

# Manual

## Downloading

The program is distributed on GitHub as a single python file, making the distribution very simple. The best way to download the program is using a local Git client. Simply open a terminal and navigate to a convenient directory before cloning the git repository. On a Linux distribution, this can be accomplished with:
```
cd ~/path/to/convenient/directory
git clone https://github.com/Sgt-Spaghetti/FD-suite.git
```
This way, whenever there is a new update to the program, it can easily be downloaded:
```
cd ~/path/to/program
git pull
```

Alternatively, the file can be downloaded from the GitHub webpage as a .zip folder, or the contents of the file can be copied and pasted into a local .py file.

## Launching
If the above dependencies are met, the program can be started from the terminal with:
```
python3 main.py
```
Alternatively, an IDE such as Pycharm may be used, which can facilitate the creation of a virtual environment. Simply load in the script from the IDE's user interface, and press the green "run" button. If prompted to install missing dependencies, simply click to accept. If the program runs, a TKinter based GUI will be displayed to the screen, which can be full-screened.

## Loading Data
The program supports loading HDF5 files (.h5) from Lumick's C- and M-Trap optical tweezer systems. The program will attempt to open the "Force LF/Foce 2x" and "Distance/Distance 1" paths in the HDF5 file for data extraction. To load data into the program:

1. Hover the mouse over the "File" tab in the menu bar at the top of the screen
2. Select "Open Folder"
3. Using the folder selection dialogue, navigate to the folder which contains the desired data in HDF5 format
> Note, the files themselves will not be visible. Simply select the folder which contains the files
4. In the folder selection dialogue, press "Ok".
5. The folder will now be loaded, all HDF5 files found in the folder will be listed in the leftmost "all files" panel

## Selecting Data

1. From the leftmost "all files" panel, click through the imported force-distance curves from top to bottom. Curves will be previewed in the rightmost data viewer automatically. Each curve will be selected for future processing if clicked in the leftmost "all files" panel. To deselect a curve, click it once more
2. By default, data will be shown in the distance domain. In order to view the data in the time, which is useful if multiple stretches have been recorded in the same data file, click on the "Toggle Time" button. This button is used to swap between the Distance and Time domains.
3. Once a satisfactory number of good quality data curves have been selected in the "all files" panel, they can be moved into the "selected files" panel for processing
4. Hover the mouse over the "Selection" panel in the topmost menu bar, and then click "Select Highlighted Curves". The curves selected in the "all files" panel will then be copied over to the "selected files" panel, where data processing occurs

> Note, every data acquisition session should include at least 3 good quality reference curves, where torsionally constrained DNA fragments are stretched to the end of overstretching (around $150$pN). These curves should be used later for per-session force correction and may also be used for automatic supercoiling density estimation. A few baselines should also be included, in which bead pairs with no DNA are "stretched", as some instruments may have systematic discrepancies in force measurements which should be corrected for.

## Baseline Subtraction

1. If baseline curves have been included in the data, select one in the "selected files" panel
2. In the top menu, under "Calibration", select "Mark Curve as Baseline"
3. Repeat steps 1 and 2 for all selected baseline curves in the dataset
4. Once all baseline curves have been marked, in the top menu "Calibration" drop-down select "Subtract Baseline"
5. All other curves in the session will now subtract a linearly interpolated average from all the selected baseline curves

## Automatic Session Force Adjustment

When measuring force-distance data of the same DNA sample, under the same buffer conditions and instrumentation setup across multiple temporal sessions, a marked variation in the overstretching force can sometimes be observed. This variation can increase or decrease the observed overstretching force by up to $\approx 10$pN across different sessions, for the exact same experimental sample. Empirical observations have shown that this variation is a multiplicative scale factor in the force domain, two force distance curves from different sessions which appear to differ in their overstretching force can be perfectly overlapped if scaled to each other in the force domain.

Averaging many different measurement sessions of commercial lambda DNA samples combined with cross-referencing with published force-distance data[^1] has shown that lambda DNA should be under $110$pN of force when $4.89\mu$m microspheres are 22.15um apart, or 1.3424x the contour length of lambda. Therefore, to generalise the relationship, force-distance curves should read $110$pN at 1.3424x the contour length of a DNA sample being measured (assuming the elastic properties of the DNA sample is comparable to that of lambda DNA, i.e $\approx 50%$ AT).

1. Select a reference curve from the "selected files" panel. The reference curve should contain at least one extension which reaches the end of the overstretching plateau
2. Mark the selected curve as a reference curve through the "Calibration" drop-down in the top menu (click "Mark as Reference") 
3. If the reference curve is a single good quality extension or one extension-retraction cycle, then set the "ymax" input field to the critical force[^2] before twist-stretch coupling (the peak in the first derivative of the extension curve), and press the "Auto Trim" button
4. Otherwise, if the reference curve is not a single good quality extension or one extension-retraction cycle, see the next section "Trimming Force-Distance Curves" before continuing
5. Repeat steps 1-4 for all desired reference curves, three is a good minimum number
6. In the "Calibration" drop-down in the top menu, press "Auto Force Scale"
7. All curves in the "selected files" panel will now be scaled in the force domain, using an average of the scale factors required to bring the reference curves to be 110pN at 1.3424x their average contour length.

> The reference curves are automatically trimmed to a single extension curve which stretches up to the critical-force. This trimmed data window is used to fit the reference curves with the extensible Odjik Worm-Like Chain Model, and the contour length of each reference curve is extracted. The contour length is averaged, and the resulting mean contour length is used to convert the reference curves from the distance domain to contour length domain (where the distance at the contour length is 1). From here, the force at 1.3424x the contour length is queried for each reference curve, and the scale factor required to bring each curve to 110pN at 1.3424x their average contour length is found. This "force scale factor" is then averaged, and the mean "session force scale factor" is applied to every force-distance curve that is selected. This process normalises each session to a standard value, eliminating per-session systematic variability.

## Trimming Force-Distance Curves

### Manually

1. Select the force-distance curve to be trimmed
2. Inspect the curve in the right-most data visualisation panel. The force-distance data is given in the top panel, the first derivative in the middle panel and the second derivative in the bottom panel. The data can be converted from the distance domain to the time domain using the "Toggle Time" button, which is useful when observing multiple extension/retraction cycles.
3. Use the sliders to select a good region for an extension curve. The slider values will be projected as vertical lines onto the data curves. If the trimmed curve will be used for fitting, the "xmax" value should be at the critical force, which coincides with the peak of the first derivative in the extension curve. If necessary, adjust the "ymax" value to allow trimming up to the critical force.
4. In the check-boxes below the data visualisation panel, check the "Set Extension" checkbox.
5. Press the "Manual Trim" button
> The data will now be trimmed such that only the data points within the box created by "xmin", "xmax", "ymin" and "ymax" are kept. This trimmed data will internally be saved as an extension curve, and can be used for subsequent fitting.
6. You should now see an updated view in the data viewer panel, showing only the points that were within the "trimmed" region. If you are unsatisfied, or to return to view the whole force-distance curve, press the "Display Full" checkbox. 
7. To set a retraction curve, first press the "Display Full" checkbox if it is not already checked.
8. Repeat steps 1-5, except select the "Set Retraction" instead of "Set Extension" Checkbox.
9. Untick the "Display Full" checkbox to view the retraction curve.
10. To swap between viewing the extension or retraction curve, press the corresponding "Set Extension/Retraction" radio button. To view the full curve, press the "Display Full" button which will override the "Set Extension/Retraction" radio buttons

### Automatically
1. If the reference curve is a single good quality extension or one extension-retraction cycle, then set the "ymax" input field to the critical force and press the "Auto Trim" button
2. You should now see an updated view in the data viewer panel, showing only the points that were within the "trimmed" region. If you are unsatisfied, or to return to view the whole force-distance curve, press the "Display Full" checkbox.
3. The extension and the retraction curves will be automatically "trimmed" from the main data, assuming there is a clean peak that separates one single extension and retraction event.

## Fitting Force-Distance Curves

1. If the selected curve has been trimmed well, such that the extension or retraction curves end at the critical force, the trimmed data that is visible in the data-viewing panel can be fit by the extensible Odjik Worm-Like Chain Model[^3]
2. A list of persistence length (Lp), contour length (Lc), Stretch modulus (S) and force offset (F0) initial guess parameters are given. These guess parameters assume standard Lambda DNA is being fit with a contour length of $16.5mu$m. If lambda DNA is not being fit, it is crucial to change the contour length guess parameter to be close to the predicted contour length of the DNA being used[^4]
> It is important to know that the fitting is highly sensitive to a contour length that is too short. It is important to be within $\approx 0.5\mu$m of the predicted contour length if underestimating the contour length, but the fitting is much more lenient when overestimating the contour length.
3. Press the "Fit" button. The display should now update showing the result of the fit as a red line going through the blue data points. The numerical values of the fit parameters will be shown after the corresponding initial guess-parameter entry boxes (rounded to 4 significant figures in this preview).

## Exporting Data

In the top menu, hover over the "File" drop-down and select "Export to CSV". A CSV file for every curve in the "selected files" panel will be exported to a folder called "FD-DATA", which is created within the directory the original HDF5 files are kept. This CSV file will contain all the data generated from the program:

- Raw unmodified force data, distance data and time data (derived from the input sample rate) 
- Processed data: The raw data after it has been baseline subtracted and per-session force scaled
- Extension curve force, distance and time data from the "trimmed" view created in the program
- Retraction curve force, distance and time data from the "trimmed" view created in the program
- The precalculated Extensible Odjik Worm-Like Chain Model fit data, for facilitated plotting of the fit curve in other programs (such as GNUplot)
- The extension curve fit parameters, unrounded: Lp, Lc, S, F0 as well was their errors[^5]
- The retraction curve fit parameters, unrounded: Lp, Lc, S, F0 as well as their errors[^5]
- Descriptive metadata for the curve (eg if it is a reference curve)

Furthermore, a file summarising all the fit parameters of every curve that has been fit in the "selected files" panel is also generated, facilitating downstream data analysis.

The program will also automatically save every graph produced as a png file.

[^1]: Leger et al. 1999 doi: 10.1103/PhysRevLett.83.1066
[^2]: This tends to correspond to the default value of 30pN for good quality torsionally constrained DNA, however using the critical-force is more accurate.
[^3]: Uses the SciPy curve-fit package for nonlinear least-squares regression using the inverted eOdjik model with force as a function of distance, adapted from the FD fit matlab package created by Onno.
[^4]: B-DNA has a rise of $0.34$nm per base-pair.
[^5]: Errors derived from square root of the diagonals of the covariance matrix.
