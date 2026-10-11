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
If the above dependencies are met, once you navigate to the folder containint the program it can be started from the terminal with:
```
python3 main.py
```
Alternatively, an IDE such as Pycharm may be used, which can facilitate the creation of a virtual environment. Simply load in the script from the IDE's user interface, and press the green "run" button. If prompted to install missing dependencies, simply click to accept. If the program runs, a TKinter based GUI will be displayed to the screen, which can be full-screened.

For very regular users, if you are not dependent on an IDE such as Pycharm it may be convenient to "alias" the command to run the program. On a linux distribution, this can usually be achieved by editing the terminal configuration. A common example is shown below:
```
cd ~
vim .bashrc
```
Then simply add the following line:
```
alias fdsuite="python3 ~/path/to/program/main.py"
```

## A Standard Proposed Workflow

- During a session, at least 3 reference curves should be taken per DNA construct investigated. These are torsionally constrained, non-supercoiled curves which reach the end of the overstretching plateau. Save the frame rate and extension speed used in the session (or per curve), as this data is used for smoothing when calculating the derivatives.
- HDF5 files from a Lumick's M- or C-Trap should be exported and saved on the local computer.
- Each session should be baseline subtracted if a systematic baseline drift error is observed when stretching "empty" bead pairs, as has been seen in some Optical Tweezer setups.
- Each reference curve should be split into extension and retraction curves. These curves should then be fit.
- The fit data from each reference curve should then be used to force-calibrate every curve in the session.
- All curves in the session can then be trimmed and fit. Reference curves should be re-trimmed and fit as force scaling can affect Lp and S.
- The supercoiling density estimations for the curves in the session can be calculated.

## First Steps

Once the program is launched, the first screen you are greeted with will ask you to confirm some core settings used throughout the measurement session. It is important to know the framerate and the extension speed used for the C-Trap's data acquisition. We would recommend a framerate of at least 100fps, with an extension speed between 0.3 and 0.15 micrometers per second. Here, you may also change data the program will used when processing the HDF5 files, which will be expanded on in the next section.

## Loading Data
The program supports loading HDF5 files (.h5) from Lumick's C- and M-Trap optical tweezer systems. By default, the program will attempt to open the "Force LF/Trap 2" and "Distance/Distance 1" paths in the HDF5 file for data extraction, although this can be changed in the initial screen shown after the program launches. The non-default options provided will still To load data into the program expecting "standard" Lumick's paths as of October 2026. Should this need to be changed, you may adjust lines 109, 110 and 112 of the program appropriately.

To load data into the program:

1. Hover the mouse over the "File" tab in the menu bar at the top of the screen
2. Select "Open Folder"
3. Using the folder selection dialogue, navigate to the folder which contains the desired data in HDF5 format
> Note, the files themselves will not be visible. Simply select the folder which contains the files
4. In the folder selection dialogue, press "Ok".
5. The folder will now be loaded, all HDF5 files found in the folder will be listed in the leftmost "all files" panel

## Selecting Data

1. From the leftmost "all files" panel, click through the imported force-distance curves from top to bottom. Curves will be previewed in the rightmost data viewer automatically. Each curve will be selected for future processing if clicked in the leftmost "all files" panel. To deselect a curve, click it once more. Shift-selection and Ctrl-selection is supported.
2. By default, data will be shown in the distance domain. In order to view the data in the time, which is useful if multiple stretches have been recorded in the same data file, click on the "Toggle Time" button. This button is used to swap between the Distance and Time domains.
3. Once a satisfactory number of good quality data curves have been selected in the "all files" panel, they can be moved into the "selected files" panel for processing
4. Hover the mouse over the "Selection" panel in the topmost menu bar, and then click "Select Highlighted Curves". The curves selected in the "all files" panel will then be copied over to the "selected files" panel, where data processing occurs

> Note, every data acquisition session should include at least 3 good quality reference curves, where torsionally constrained DNA fragments are stretched to the end of overstretching (around $150$pN). These curves should be used later for per-session force correction and may also be used for automatic supercoiling density estimation. A few baselines should also be included, in which bead pairs with no DNA are "stretched", as some instruments may have systematic discrepancies in force measurements which should be corrected for.

5. For each reference curve selected in the session, select the curve from the "selected files" panel and toggle it as a reference, using the "Selection" menu drop-down followed by the "Mark Curve as Reference". The name of the curve should now be green in the "selected files" panel.
6. For each baseline selected in the session, select the curve from the "selected files" panel and toggle it as a reference, using the "Selection" menu drop-down followed by the "Mark Curve as Baseline". The name of the curve should now be red in the "selected files" panel.

## Baseline Subtraction

The program supports automatic baseline subtraction, as well as using pre-calculated baselines if convenient.

To proceed with automatic baseline subtraction:

1. If baseline curves have been included in the data, select one in the "selected files" panel
2. In the top menu, under "Calibration", select "Mark Curve as Baseline"
3. Repeat steps 1 and 2 for all selected baseline curves in the dataset
4. Once all baseline curves have been marked, in the top menu "Calibration" drop-down select "Calculate Baseline"
5. To view the resulting averaged baseline, in the top menu "Calibration" drop-down select "View Baseline". A graph should pop up in a separate window showing the averaged baseline curves
> The calculate baseline may now be saved as a csv file. In the top menu "Calibration" drop-down select "Save Baseline"
6. To apply the calculated baseline to all the curves in the session, from the top menu "Calibration" drop-down select "Subtract Baseline"
7. All other curves in the session will now subtract a linearly interpolated average from all the selected baseline curves

To use a pre-calculated baseline:

The program will expect a csv file with two columns, with headers "Distance" and "Force"

1. In the top menu "Calibration" drop-down select "Load Baseline"
5. To view the resulting averaged baseline, in the top menu "Calibration" drop-down select "View Baseline". A graph should pop up in a separate window showing the averaged baseline curves
6. To apply the baseline to all the curves in the session, from the top menu "Calibration" drop-down select "Subtract Baseline"
7. All other curves in the session will now subtract a linearly interpolated value derived from the loaded baseline curve

## Trimming Force-Distance Curves

Force distance curves loaded into the program must be processed in two steps. If the data is very clear with only a single extension/retraction curve up to 115pN, the program can attempt to trim and process the curves automatically, however this should not be blindly relied upon. We strongly suggest manually trimming all curves.

## Manually trimming Extension and Retractions

Firstly, the data should be split into designated Extension and Retraction curves, and their associated "fit range" to be used for subsequent data fitting (i.e generally up to 30pN) should be defined. If the data only consists of a single extension or retraction without its counterpart, the other entry will simply be empty and will be ignored in future steps.

1. Once the curves have been selected and baseline subtraction has been carried out, click on the desired curve in the "selected files" panel. If not already done, switch to the Time domain by clicking the "Toggle Time" button. In the time domain, move the xmax slider such that it goes past the force peak and slightly into the retraction section of the data.
2. Press the "Peak" button to automatically move the xmax slider back to the point of largest force.
3. Ensure that the "Set Extension" and "View Full" checkboxes are filled.
4. Press "Trim Data". The "View Full" radio button should now deactivate, and automatically switch to "View Extension". The trimmed extension curve should now be shown. (To go back to viewing the complete curve, simply press the "View Full" button). The "Set Fit" checkbutton should now also have automatically been selected.
> As the data is not yet force calibrated, the following steps should initially only be done on the reference curves in the session.
6. With the "View Extension" and "Set Fit" buttons checked, move the xmax slider just past the sharp rise in the force-distance curve (i.e past the elastic extension stage, or past the peak in the first derivative).
7. Press the "Fc" button. This will automatically bring the xmax slider back to the peak of the first derivative (to the "critical force").
8. At this stage, one may switch back to the Distance domain and set the xmin slider to a set distance, such that all trimmed curves will have an equal amount of "baseline" before the sharp rise with elastic extension; usually only data from 10 micrometers onwards is used.
9. Press "Trim Data". The "View Fit" button should now automatically activate, and only data between 10 micrometers and ~30pN should be displayed. This data can be used for curve fitting. Do note, if the data is not force calibrated, steps 6 onwards should only be done on the reference curves. This is because subsequent force calibration will change the fit parameters S and Lp, but not Lc which is initially used for force calibration.
10. The above steps should be repeated on the retraction curves. The process is very similar, simply the "Set Retraction" checkbutton should be selected instead of "Set Extension". Do not forget to switch to the Time domain before trying to initially trim the retraction curves.


## Fitting Force-Distance Curves

This process should only be done on the reference curves unless the session has been force adjusted, in which case all curves should then be (re-)fit to extract accurate parameters.

1. If the selected curve has been trimmed well, such that the extension or retraction curves end at the critical force, the trimmed data that is visible in the data-viewing panel can be fit by the extensible Odjik Worm-Like Chain Model[^3]
2. A list of persistence length (Lp), contour length (Lc), Stretch modulus (S) and force offset (F0) initial guess parameters are given. These guess parameters assume standard Lambda DNA is being fit with a contour length of $16.5mu$m. If lambda DNA is not being fit, it is crucial to change the contour length guess parameter to be close to the predicted contour length of the DNA being used[^4]
> It is important to know that the fitting is highly sensitive to a contour length that is too short. It is important to be within $\approx 0.5\mu$m of the predicted contour length if underestimating the contour length, but the fitting is much more lenient when overestimating the contour length. Any combination of parameters may be fixed if desired.
3. While viewing the data associated with the "View Fit" checkbutton, press the "Fit" button. The display should now update showing the result of the fit as a red line going through the blue data points. The numerical values of the fit parameters will be shown in the boxes to the right of the "Fix" checkbutton (rounded to 4 significant figures in this preview).


## Automatic Session Force Adjustment

When measuring force-distance data of the same DNA sample, under the same buffer conditions and instrumentation setup across multiple temporal sessions, a marked variation in the overstretching force can sometimes be observed. This variation can increase or decrease the observed overstretching force by up to $\approx 10$pN across different sessions, for the exact same experimental sample. Empirical observations have shown that this variation is a multiplicative scale factor in the force domain, two force distance curves from different sessions which appear to differ in their overstretching force can be perfectly overlapped (when plotted in Lc space) if scaled to each other in the force domain.

Averaging many different measurement sessions of commercial lambda DNA samples combined with cross-referencing with published force-distance data[^1] has shown that lambda DNA should be under $110$pN of force when $4.89\mu$m microspheres are 22.15um apart, or 1.3424x the contour length of lambda. Therefore, to generalise the relationship, we suggest that force-distance curves should read $110$pN at 1.3424x the contour length of a DNA sample being measured (assuming the elastic properties of the DNA sample is comparable to that of lambda DNA, i.e $\approx 50%$ AT).

1. If reference curves for a session have not yet been marked, sequentially select a reference curve from the "selected files" panel. The reference curve should contain at least one extension which reaches the end of the overstretching plateau
2. Mark the selected curve as a reference curve through the "Calibration" drop-down in the top menu (click "Mark as Reference") 
4. If the reference curve is not a single good quality extension or one extension-retraction cycle, see the section "Manually Trimming Extensions and Retractions" before continuing
5. Repeat steps 1-4 for all desired reference curves, three is a good minimum number
6. If the marked reference curves have been well trimmed and fit as described previously, then in the "Calibration" drop-down in the top menu, press "Auto Force Scale"
7. All curves in the "selected files" panel will now be scaled in the force domain, using an average of the scale factors required to bring the reference curves to be 110pN at 1.3424x their average contour length.

> In general, the program can attempt to automatically trim the extension curves to a single extension curve and retraction curve, each of which then stretches up to the critical-force. However, this automatic process should not be relied upon. The well-trimmed data is used to fit the reference curves with the extensible Odjik Worm-Like Chain Model, and the contour length of each reference curve is extracted. The contour length is averaged, and the resulting mean contour length is used to convert the reference curves from the distance domain to contour length domain (where the distance at the contour length is 1). From here, the force at 1.3424x the contour length is queried for each reference curve, and the scale factor required to bring each curve to 110pN at 1.3424x their average contour length is found. This "force scale factor" is then averaged, and the mean "session force scale factor" is applied to every force-distance curve that is selected. This process normalises each session to a standard value, eliminating per-session systematic variability.

Every curve in the session should now be force scaled, such that the overstretching plateau occurs at 110pN. In general, the reference curves can be overlaid to check that they are internally consistent through the "Overlay References" option in the "View" drop down menu.

## After force scaling and Supercoiling Density Estimation

Each curve can now be trimmed and fit. If a session includes reference curves, the supercoiling density of a curve can be estimated by measuring its distance in Lc space at 70pN. Using a calibration plot derived from data presented in Leger et al. 1999 [^1], the distance at 70pN will be converted to supercoiling density and this is the most accurate method for supercoiling density estimation. However, should a curve not reach 70pN, any distance from 5pN to 30pN can be used instead. These distances use a secondary internal calibration plot, which converts lower forces to supercoiling density using data derived from many lambda-DNA curves at different supercoiling densities which were themselves determined using the distance at 70pN method.

> Supercoiling density estimation requires reference curves, as a tiny subtractive correction factor ($\approx 5x10^{-3}$) is applied to all curves in order to zero the estimated supercoiling density of the reference curves. This is because there can be small variations in the twist-stretch coupling region even if all other areas of the curves overlap after force-scaling in Lc space, which may derive from small differences in local heating.

## Exporting Data

In the top menu, hover over the "File" drop-down and select "Export to CSV". A CSV file for every curve in the "selected files" panel will be exported to a folder called "FD-DATA", which is created within the directory the original HDF5 files are kept. This CSV file will contain all the data generated from the program, including:

- Raw unmodified force data, distance data and time data (derived from the timestamp associated with "Distance" value, with the very first timestamp subtracted from all subseqeuent timestamps and converted from nanoseconds to seconds) 
- Processed data: The raw data after it has been baseline subtracted and per-session force scaled
- Extension curve force, distance and time data from the "trimmed" view created in the program
- Retraction curve force, distance and time data from the "trimmed" view created in the program
- The precalculated Extensible Odjik Worm-Like Chain Model fit data for both extension and retraction curves, for facilitated plotting of the fit curve in other programs (such as GNUplot)
- The extension curve fit parameters, unrounded: Lp, Lc, S, F0 as well was their errors[^5]
- The retraction curve fit parameters, unrounded: Lp, Lc, S, F0 as well as their errors[^5]
- Descriptive metadata for the curve (eg if it is a reference curve)

Furthermore, a file summarising all the fit parameters and supercoiling densities of every curve that has been fit in the "selected files" panel is also generated, facilitating downstream data analysis.

The program will also automatically save every graph produced as a png file.

[^1]: Leger et al. 1999 doi: 10.1103/PhysRevLett.83.1066
[^2]: This tends to correspond to the default value of 30pN for good quality torsionally constrained DNA, however using the critical-force is more accurate.
[^3]: Uses the SciPy curve-fit package for nonlinear least-squares regression using the inverted eOdjik model with force as a function of distance, adapted from the FD fit matlab package created by Onno.
[^4]: B-DNA has a rise of $0.34$nm per base-pair.
[^5]: Errors derived from square root of the diagonals of the covariance matrix.
