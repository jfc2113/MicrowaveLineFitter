This file describes the contents of MicrowaveLineFitter, which contains code used in Corby et al. (2015), Monthly Notices of the Royal Astronomical Society.  The code can be applied to any baseline-subtracted broadband spectral line data.  It performs automated line fitting and semi-automated line identification.  We provide instructions for use including a description of the input and output files and system requirements.

The scripts have been tested in IPYTHON version 0.12.1 with PYTHON 2.7.3, NUMPY 1.6.1, and MATPLOTLIB 1.1.1rc. 
PYTHON 2.7.3 includes the MATH and CSV modules that are used in the scripts.  A full description of the code follows below.  

First, we provide a short set of directions:  

1. Edit the parameters at the top of FullCode.py and throughout VelocityEditFile.py based on your source kinematics and data range.

2. Run the code in ipython with: 

         %run FullCode.py

         cd VEL
         
         Optional:  Manually adjust output of FullCode.py to improve fits that are flagged as bad (see bins below).  

         %run VelocityFileMaker.py

The final output files are named velocity_ALL_(r).csv where r is the name of the region or source of interest. 

Second, we provide a description of the code and process:

The MicrowaveLineFitter repository includes eight scripts, five spectra towards positions in Sgr B2 (see Corby et al. (2015)), and csv files holding line data.  Two of the scripts are run directly, namely  'FullCode.py' is run before 'VEL/VelocityFileMaker.py'.  After optional adjustments by the user, 'VelocityFileMaker.py' uses the output of 'FullCode.py' to generate a final csv file that contains Gaussian fits to the full spectrum and line identifications.  Of the remaining six scripts, five contain modules called in 'FullCode.py'.  The modules and their primary purposes include:

1. 'LineFind.py' grabs segments of the spectrum that contain lines, detecting "line-channels" with single-channel intensity values greater than Sigma_{rms} times an input threshold.  The module also selects a set number of channels on either side of  line-channels.  The number of neighboring channels has a default value of 5, but can be adjusted when LineFind is called.

2. 'FTinterp.py' interpolates a segment of the spectrum, outputting a spectrum with a factor of 3 improvement in the spectral resolution.  It performs the interpolation by Fourier-domain zero padding.  This should only be used for data with poor spectral resolution compared to the linewidths.

3. 'Gfit.py' provides the most essential functionality; through 'FullCode.py', the function Gfit.GauFit is passed a segment of the spectrum that contains one or more lines.  Gfit.GauFit returns a 1- or 2-component Gaussian fit to a spectral feature within the segment.  The script selects a specific portion of the line-containing segment and first attempts a one-component fit.  It then follows an algorithm to adjust the section of the segment used for line-fitting and re-fit 1-component Gaussians.  It passes information to 'G2criteria.py', which  evaluates the fit against a set of criteria to determine whether a 2-component Gaussian fit is justified.  The output of G2criteria determines whether Gfit.GauFit fits a 2-component Gaussian, following an algorithm to adjust the section used for line-fitting and attempting 2-component Gaussian fits.

4. 'G2criteria.py' determines if a 2-component Gaussian fit is required.  Gfit.GauFit provides G2criteria with information about the data and the best 1-component Gaussian fit achieved; G2criteria then determines if the data differs significantly from the 1-component fit, requiring a 2-component fit.  G2criteria takes into account the  signal-to-noise ratio of the line in question and the degree to which the data differs from the 1-component Gaussian fit.  The primary factors considered include: (1) the difference between the position of the absolute maximum of the data and the position of the 1-component Gaussian center; (2) the maximum residual between the data and the fit, compared to the rms noise level at that frequency tuning; (3) the mean residual between the data and the fit, compared to the rms noise level; (4) the width of the fit compared to the typical width of features.  The criteria, mostly empirically determined, works very well for the work referenced above.  However, we advise inspecting results obtained on other datasets carefully and likely adjusting the criteria, as they will change with the parameters of the data.  

5. 'LineID.py' provides semi-automated line identification of features fit over the full spectrum.  FullCode compiles the output of Gfit.py and passes it to LineID, along with one or two csv files containing spectral line catalog data.  The input csv files can quickly and easily be output from the ALMA Spectral Line Catalogue - Splatalogue (Available at www.splatalogue.net; Remijan et al. (2007)) after selecting a set of molecules that might appear in your spectrum.  
   
FullCode manages the full data spectrum, compiles the results of Gfit.GauFit, and calls LineID to output the data.  All of the modules except G2criteria are written to 
be easily tunable, and the main parameters that users will want to adjust are specified and explained in the first section of FullCode.  These parameters include:
whether the results should be output; what plots should be produced; the names of input files; whether data interpolation is desired; 
information about the source including guidelines for appropriate linewidths and line velocities; and information about the data structure 
(e.g. how many separate data sections are in a single spectrum; in this work, the line-fitter inspected a spectrum including 11 tunings, each of which had different Sigma_{rms}
noise levels.  If all tunings have the same noise level in the intensity units you are using, they can be treated as a single section).

FullCode requires a spectrum input in standard ascii format with frequency and baseline subtracted intensity.
The file 'K6_fullSpec.txt', provided with the scripts, is in the correct format.  If the user chooses to perform line identification, the script requires files containing recombination line or molecular line data output in csv format from any spectral line catalog.  The line data should be formatted similarly to the files 'recombLines.csv' and 'molecularLines.csv', provided.  If your source has recombination lines, we recommend preserving the file name of 'recombLines.csv', but including lines that are relevant for your frequency range.  All csv files read into the code should be text csv files with unicode UTF-8 characters and a colon (:) as the field delimiter. 
The code outputs csv file data with this format as well.

FullCode plots the full data spectrum with the line fits overlaid, marking the location of spectral lines near the line fits from input line data. It also can output a csv file containing Gaussian line fits and possible line assignments.  This output file is named 'ALLFITS_(r).csv', where (r) is the name of the source or region targeted.  It contains:

1. Gaussian line fit parameters (height, center frequency, and width in frequency units) of all features.

2. two measures of the quality of the fit.  First, the rms residual between the data and the fit; second, for 1-component fits, the file provides the offset between the frequency at which the data has its maximum absolute value and the Gaussian center frequency. Fits with abnormally large values are likely suspect. 

3. the "bin" indicates how the line was handled in Gfit.GauFit. Lines with bin = 0 are 1-component fits that were deemed appropriate.  Lines with bin = 1 are 2-component fits that were found to be reasonable.  For lines with bin = 2, the data significantly differed from the best 1-component fit obtained, however Gfit.GauFit was unable to obtain a reasonable 2-component fit and therefore used the best 1-component fit obtained.  These should be inspected and possibly re-fit manually.  For bin = 3 lines, the residual between the data and the fit is significantly larger than the noise level of the spectrum.  These should be inspected and possibly re-fit manually.  Bin = 5 lines are entirely unreasonable, with very broad line widths or unreasonable line heights.  They are used as placeholders to enable the code to proceed, and they will need to be re-fit.  Of 617 components used to fit the spectrum towards K6, 360 have bin = 0, 242 have bin = 1, one has bin = 2, 12 have bin = 3, and two have bin = 5.

4. possible atomic or molecular carriers of the detected Gaussian features.  The output file lists any transitions from the input 'recombLines.csv' and 'molecularLines.csv' files that are within a specified velocity range of the detected features.  

Upon inspection of the results and manual updates to any bad fits, run VelocityFileMaker to generate a final output csv file containing line identifications and Gaussian fit parameters, with parameters listed in velocity space. Before running VelocityFileMaker, edit parameters in 'VelocityEditFile.py', imported to VelocityFileMaker as 'VE'.  Specify the names of the positions from which you extracted spectra, csv file input names, and velocity ranges that are appropriate for different types of lines, including recombination lines, molecular lines from the primary source, and additional lines observed in absorption by foreground diffuse or translucent clouds.  The latter is only relevant in some lines-of-sight.

The code implements primarily kinematics-based consistency checks for line identifications, and can use three optional input csv files to prioritize line identifications. In order to minimize mis-identification, the final line-identification process operates on a priority system as follows:

1. Priority 1: The code first inspects identified recombination lines for kinematic consistency and consistent line-intensities.  This step requires a sufficiently large number of recombination lines to work well (>10 high signal-to-noise lines or >20 total lines recommended), as it derives the mean Gaussian fit parameters and standard deviations in order to evaluate consistency.  Lines that are self-consistent with the mean and standard deviations are output as firm line identifications, and assigned a Line Type of either "H Recomb" or "He Recomb".
The script is not presently equipped to handle carbon recombination lines, but could be adjusted to do so.
Features that are inconsistent with recombination lines include absorption lines, features that are significantly more broad or narrow than the mean, features that are too far from the mean center frequency,or lines that are significantly too strong or weak (e.g. an H-gamma line is as strong as most H-alpha lines).  These may be output as tentative line identifications, which the user should inspect for line blending or data issues.  If you do not have recombination lines, the code will simply move to Priority 2.  

2. Priority 2: The user has the option of specifying a file entitled 'strongLines.csv'. (Of course, the file name can be changed per the users preference).  The file 'strongLines.csv' contains the strongest lines, which can be firmly identified prior to identifying weaker molecular lines.  Lines that fall within the 
 velocity range set in the function VE.strongVrange (in VelocityEditFile) are output as firmly identified, and assigned the Line Type "strong".  Verify that the first two columns of 
 your 'strongLines.csv' file are formatted the same as the file included in the package, which is different than the 'molecularLines.csv' file.

3. Priority 3: The user can use a file entitled 'SAClines.csv', which indicates which transitions that may be observed in foreground absorption by diffuse and translucent cloud material, often referred to as Spiral Arm Clouds (abbreviated SAC). The foreground absorption components are at velocities that are inconsistent with the targeted cloud, and must be constrained in VelocityEditFile.  Lines that fall within the velocity range set in the function vEDIT.SACvRange are output as firmly identified, and assigned the Line Type "SAC".  Verify that the first two columns of your 'SAClines.csv' file are formatted the same as the file included in the package.

4. Priority 4: The code then assigns remaining lines tentatively assigned in 'ALLFITS_(r)'.  In cases where multiple transitions are near the Gaussian center frequency, the line with the best kinematic match to the velocity, as specified in VelocityEditFile, is output.  While this typically works well 
at centimeter wavelengths, we recommend inspecting lines that have multiple line entries in 'ALLFITS_(r)', especially at higher frequency.  
If the code gets the "wrong" answer at this point, mark the line as blended, and either change the line assignment in velocity_ALL_(r) 
or only keep the preferred line in ALLFITS_(r).

5. Priority 5: Finally, the code outputs all fits that were not firmly identified in Priorities 1-4.  These are given a Line Type of "Unidentified".  If a fit was associated with a known transition in ALLFITS_(r), the transition parameters will be output so that the user can decide whether or not to adopt the transition as identified.

The code does not handle line blending with much sophistication, however the prioritization system significantly improves the likelihood that most of the line radiation can be ascribed to the identified carrier transition.

Within each priority, the VelocityFileMaker also handles wing components.  Some spectral features are non-Gaussian, including recombination lines, strongly masing lines, and optically thick molecular lines.  In FullCode, in which the line fitting is conducted, after subtracting a best-fit 1- or 2-component Gaussian fits from the raw data (called the primary component(s)), high signal-to-noise lines may have wing components with intensities exceeding the detection threshold.  Such wing features, which are typically weak as compared to the primary 1- or 2-component fit, are fit to the residual.  VelocityFileMaker selects these features to prevent them from incorrectly being mis-assigned to other molecular carriers or labelled as unidentified.  Wing features are selected based on (1) their separation from the primary line component as compared to the width of the primary line component and (2) their height as compared to the primary line component.  If a feature meets either of the following criteria, it is identified as a wing of an identified line and the velocity is calculated with respect to the rest frequency of that line. The criteria include:

         1. H_{Feature}/H_{Primary} < 0.22 and v_{Feature} - v_{Primary} < deltaV_{Primary}
         2. H_{Feature}/H_{Primary} < 0.05 and v_{Feature} - v_{Primary} < 1.5 * deltaV_{Primary}
         
for the line heights, velocity centers, and FWHM widths (H, v, and deltaV) of the possible wing component (Feature) and the primary component (Primary).  These are assigned a LineType of "wing". Within each priority, VelocityFileMaker identifies the primary components, which must be consistent with the kinematic information input, and then selects any wings associated with these lines.  As such, the recombination lines and wings of recombination lines are identified, followed by designated "strong lines" and the wings of strong lines, etc.  

VelocityFileMaker generates an output csv file entitled "velocity_ALL_(r).csv".  The file 
contains line rest frequencies, species, transitions, Gaussian fit velocities, heights, widths (in velocity units), 
center frequencies, and line types.
