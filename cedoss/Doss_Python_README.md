# Repository of Christopher Doss

## Thin Plasma Lens Studies

### Introduction


Greetings!  I am writing this in the brief window of June 2023 between finishing my PhD and leaving CU.  If you are reading this then I assume you are here to try and decypher the contents of this folder.  This is unfortunate since all of these scripts were written and rewritten over a period of about 7 years with varying levels of skill and care.  Since I was the only one using these scripts, I adopted bad coding habits that persisted throughout my full PhD.  If you want to be better than me, I recommend adopting good coding practices early (Object-oriented design, lots of comments, input files, output files, etc.)  But anyways, here I review some of the main functionalities in my scripts.  My goal here is to point out which scripts do what, and so you can take the things that work and build a better, more legible code.  Good Luck!
 -Chris

### Table of Contents
-01- Generic Modules - '/modules/'

-02- Complex Q Parameter Propagation - '/UpdatedCodes/QParamProp.py'

-03- Calling Robert's Fourier Refraction - '/UpdatedCodes/FourierRefractionSetup.py'

-04- Plotting results of Fourier Refraction - '/UpdatedCodes/RefractionAnalysisAndPlots.py'

-05- Plotting the Oide Limit - '/UpdatedCodes/SampleOideCalculation_SigVsL.py'

-06- Restricted TPL Matching - '/pyscripts/RestrictedTPLMatching.py'

-07- Plasma Lens Calculations - '/tplscripts/PlasmaLensCalc.py'

-08- Numerical Beam Propagation - '/beamprop_v2/'

-09- Poisson Solver for Electric Fields -/PoissonSolver/'

-10- Generic SLAC Data Analysis Code - '/SLACData/'

-11- Plasma Lens Commissioning Shift Analysis - '/SLACData/'



-----------------------------------------------------------------
### -01- Generic Modules - '/modules/'

Here I will go through some of the contents of this folder, which I use in other scripts.  The ones I don't list are ones I don't recall using too much of.

CalcEmitGrowth.py - Contains functions for calculating the emittance growth in a plasma lens.  I have functions that are either normalized to sqrt(k), or functions that are unnormalized and you just use SI units.  In addition to the thin and thick Chromatic Amplitudes, I also have functions to calculate the projected betafunction in the thin and thick regimes.  For the thick regime, I have a function for either a flattop energy spread or a Gaussian energy spread.

FourierRefraction.py - These functions facilitate running Robert's Fourier propagation code to simulate the effects of refraction.  'GetDefaultParams' initializes a params dictionary with some default values that I typically change as needed.  'RunRefraction' runs Robert's simulation, but we need to have the the Efunc and density defined.  Density is a 3D numpy array, and Efunc is loaded into params dictionary using 'LoadPulseParams' from the Gaussian beam parameters in QParamProp.py (Chap -02-)

GaussianBeam.py - Used mostly with QParamProp (Chap -02-).  Three main components: (1) functions for getting the intensity and electric field profile from the spot size and phase curvature for an ideal Gaussian beam.  (2) functions for initializing the complex Q parameter and propagating it through various optics and interfaces.  (3) functions to grab the necessary parameters to describe the Efield of the Gaussian beam.

OideCalc.py - Functions for calculating various effects for the Oide effect.  'SigMin' is just the beam size from the synchrotron radiation, and 'SigOide' is the specific Oide limit for the particular focusing optic.  SigMin only equals SigOide if the betafunction at the focus is equivalent to 'BetaMin.'

ThreeDimensionAnalysis.py - This module was basically my scratch space for making functions to plot 3D numpy arrays.  Usually I only needed to plot a 2D planes along each of the three axes of x=0, y=0, and z=0.  For some reason I also have some functions to calculate the ionization fraction from an intensity profile here as well.  'VarianceCut' plots 1D lineouts at various offsets for a specific 2D plane, and 'ImageCut' plots the three 2D planes on each axis.  A good example of this in action is in 'RefractionAnalysis.py'.  Lastly, the second half of this module as a ton of functions for fitting 1D and 2D data to various functions.  The generic fitting function is 'FitDataSomething'

TPLFocalLength.py - Similar to CalcEmitGrowth, this module has a ton of functions for calculating various things for plasma lenses.  This module is more focused on the betafunction at the focus and the focus location.


-----------------------------------------------------------------
### -02- Complex Q Parameter Propagation - '/UpdatedCodes/QParamProp.py'

These scripts propagate the complex q parameter of an idealized Gaussian laser pulse through optics and interfaces between regions of different refractive indices.  This q parameter contains the spot size and radius of curvature for a Gaussian laser, and with this information we can quickly model the electric field and intensity profile at the focus or at an arbitrary position.  And since the transverse components of a laser are separable, we can perform this propagation with cylindrical optics to get an asymetric profile of the Gaussian laser pulse.

I have three generations of the main script.  The first iteration, '/tplscripts/OpticalSetup.py', ended up quite messy because I defined both the laser and optics using if/else blocks at the start of the script rather than dedicated input files.  The second iteration, '/tplscripts/OpticalSetup2.py', is much like the first but included an additional feature of adding an M^2 to the beam.  M^2 is a general parameter to represent the non-Gaussian-ness of a beam, and effectively increases the spot size while keeping the Rayleigh length the same.  The third iteration, '/UpdatedCodes/QParamProp.py', is a cleaned up version of OpticalSetup2.py with more comments.

One of the more important parameters to keep track of is 'zi'.  This sets the longitudinal window width over which the zoomed-in analysis takes place.  But also this controls at what position the electric field is calculated when importing into Robert's refraction simulation.  Additionally, I have a quick calculation of the peak intensity at the position of 'zi', so you can quickly change 'zi' to get the script to spit out the intensity at various positions.

BUG NOTE: Importing into Robert's Refraction simulation should only be done if M^2 = 1.  For larger M^2, I do not correctly calculate the electric field.  Probably can and should be fixed.

I import 'modules/GaussianBeam.py', which contains functions for propagating the q parameter through various optics and interfaces.  This module also contains functions for calculating the intensity and electric field of Gaussian beams using the spot size and phase curvature.  Also functions for saving paramters of the electric field.

At the top I have a lot of flags for doing various calculations.  'calcfocal' is very quick and just computes a 1D plasma lens density profile.  'calcdensity' is much slower but computes the full 3D plasma profile.  'save' saves the electric field parameters to the specified filename, and this is used when importing into FourierRefractionSetup.py (Chap -03-)

-----------------------------------------------------------------
### -03- Calling Robert's Fourier Refraction - '/UpdatedCodes/FourierRefractionSetup.py'

This is a close copy of '/tplscripts/FourierRefractionLooper.py', although the I never really used that code for its "loop".  Here I gather all of the parameters required for running Robert's refraction code through FourierRefraction.py.  Can also make your own variable loop close to the bottom with the 'var_loop' array, but I just have it set up for a single simulation here. It's a pretty straightforward script with flags at the top for selecting what you want to run, but there are two important things to keep in mind to use this correctly.

First, make sure that you have the correct files specified.  If we are using a gas density from a 3D numpy file, that is named in 'denfile'.  For an example on how one of these is generated, see '/pyscripts/ApproxGasJet.py'.  The 'pulsefile' must match the 'path'+'filename' from QParamProp to use the ideal Gaussian pulse generated from that script, and only M2=1 is valid currently.  Then, 'outputpath' is just the folder where the outputs of the refraction simulation go

Second, you need to make adjustments to the default params dictionary to properly configure the simulation.  Most importantly, 'X' and 'Y' need to be adjusted to be large enough for the spot sizes, 'Z' needs to match '2*zi' from QParamProp if you want the Fourier propagation to have the laser spot size in the longitudinal center, and 'EI' should be whatever ionization energy we are using.


-----------------------------------------------------------------
### -04- Plotting results of Fourier Refraction - '/UpdatedCodes/RefractionAnalysisAndPlots.py'

This one still rough around the edges, but here you can get the general idea of what can be plotted after running FourierRefractionSetup.py.  Basically just load in output folder from Chapter -03-, specify the density and window size, and a ton of plots will show up.  The current plots are a three-panel of the planes along the origin axes, density lineouts with offsets in various directions, an analysis of how the approximate plasma lens thickness varies with horizontal and vertical offsets, and some fits of the density lineouts to the standard Double-Tanh profile.  At the end there is even an effort to fit the 2D profile in the longitudinal-vertical plane to a 2D double tanh, which I used when simulating plasma lenses in PIC simulations.


-----------------------------------------------------------------
### -05- Plotting the Oide Limit - '/UpdatedCodes/SampleOideCalculation_SigVsL.py'

This one is a copy-and-paste of the original in 'tplscripts/VaryLength_OideCalc.py', but that folder has a ton of different iterations and this one is by far the most interesting.  Just a simple plot of the Oide functions given initial parameters, but by varying the plasma lens thickness L we can get a nice transition from ideal focusing to a synchrotron radiation dominated focus.  The final plot compares the spot size from ideal focusing (orange) to that of just chromatic emittance growth (red) and of just synchrotron radiation (blue).  The 'synchrotron radiation' curve only matches the 'Oide limit' (green) curve when the betafunction at the focus (beta_f^*) is equivalent to the optimal betafunction (beta_opt^*), as calculated from OideCalc.Calc_BetaMin


-----------------------------------------------------------------
### -06- Restricted TPL Matching - '/pyscripts/RestrictedTPLMatching.py'

Just a quick one, but useful for PWFA matching.  If you know the matching condition for the PWFA and the initial betafunction waist for a witness bunch, then the plasma lens location and thickness can be calculated that will achieve matching.  Can see my thesis for more on this, but just a quick calculation.


-----------------------------------------------------------------
### -07- Plasma Lens Calculations - '/tplscripts/PlasmaLensCalc.py'

Another quick but useful script: this one ended up being my plasma lens calculator.  Just a place I could define a bunch of parameters for the ebeam and the plasma lens, and I would just run through and spit out all of the information about the resulting spot size.  For an example of the information you can get from here, this is what an output looks like:
```
Focal length [m]:  0.04093649481481482
sqrt(K)*L:  0.494254456155
Beam Density:  6459236281370099.0

The following from Thick lens equations
Focusing strength K [m^-2]:  2442.8746742873013
Emittance growth:  1.20960456961
Final Emittance [m-rad]:  5.17831716249e-05
Initial rms size [m]:  5.74739734165e-05
Centroid Focus location [m]:  0.0375165281691
Centorid Focus beta [m]:  0.00120395151456
Focus rms size [m]:  1.9633400825e-06

Oide Limit Calculations
SR beam size [m]:  1.62288334579e-06
```

-----------------------------------------------------------------
### -08- Numerical Beam Propagation - '/beamprop_v2/'

Everything in this folder (probably) has to do with the numerical propagation of macroparticle beams through plasma and free space.  Basically a more rough version of WARGSim.  This folder is a successor to '/beampropagation/', and is generally a bit better than that version.  I also won't go through everything here, but I will try to touch on the highlights.

First of all, 'BeamPropFuncs.py' is my general module containing functions for interfacing with the particle tracker code.  There are functions to get default parameters, functions to propagate through long and thin plasmas, and functions to propagate just the CS parameters of a beam.  Electron beams can also be initialized with Gaussian parameters (GaussianBeam) or with the output of a VSim simulation (VorpalBeam)

NOTE:  this was super bad coding practice, but some of the parameters are hard-coded at the top of 'BeamPropFuncs.py'.  I think in particular you need to be sure to change def_nset and def_deltaE to simulate the correct plasma density and beam energy spread.  If you remake these scripts I highly recommend not doing this.

All of the scripts with 'SingleBeamPass' or 'TPL_Chromaticity' are usually just propagating this macroparticle beam through either a plasma lens or a plasma lens into a PWFA.  There are just a billion files because nearly everytime I wanted to simulate something else I would copy-paste the previous version into a new file and modify as needed.  Of these, the nicest one to look at for inspiration is 'SingleBeamPass_TPL_LiOven.py', which uses the restricted TPL matching alogrithm to match into FACET-II's Lithium Oven with a very thin plasma lens far upstream.

All of the scripts that start with 'Vorpal' have to do with analyzing the output of VSim simulations.  Of these, 'VorpalPropagate.py' is a useful one that simply loads in the Witness bunch from a VSim output into the numerical macroparticle propagation to simulate the beam passing through free space.  This one could be upgraded to simulate the beam propagating through a magnetic spectrometer, for instance.  'VorpalCrossSection' codes have to do with plotting the cross section of wakes, particularly for my second paper on linear density gradient plasma wakes.  Of these, 'VorpalCrossSection_Looper2.py' loops over longitudinal slices in a single directory and does a ton of analysis work.  'VCS_Sheath_Jz.py' plots the electron density in the sheath for the central axis and finds the boundaries of the off-center slices.  The OG 'VorpalCrossSection.py' loads in a 2D cross-section from a VSim output, plots the density, and fits the sheath boundary to a 'circular' function.  I say 'circular' because I really just fit the sheath boundary in polar coordinates to a cosine function, which is only really circular if that amplitude is small.  Could be better.

Lastly, the scripts that start with 'Paper2' are for plotting various figures in that second paper.  The rest in this folder are ones I've either forgotten about or are forgettable themself.


-----------------------------------------------------------------
### -09- Poisson Solver for Electric Fields -/PoissonSolver/'

I wrote up a quick algorithm to take a charge distribution on a 2D plane and calculate the electric potential from that distribution.  In the 'Ecalcer_linecharge.py' script, this was modified to tread the 2D charge distribution as a distribution of line charges, which allows for a full 3D contribution from a plasma distribution.  This script is able to calculate out the electric fields from a given ion column in a blowout plasma wake.  The base version does not include any effects from the electron sheath, but could in the future.  The issue is computational resources and resolution.

'Ecalcer_linecharge.py' assumes the wake is circular, and the density can either be uniform within the ion column or have a linear density gradient.

'Ecalcer_linecharge_ellip.py' assumes the wake's shape is elliptical rather than circular.  For a uniform density ion column this produces a slight asymmetry in the focusing in x than y.  Didn't work too much on this subject, but this wake shape does come up if the plasma blowout radius is larger than the plasma dimension for one transverse axis.  This is particularly interesting for the thin plasma lens design using a transversely propagating laser with a small spot size.

'Ecalcer_linecharge_ring.py' removes the ion column and just has a ring of negative charge.  This is what I used to find the electric field for an asymmetric sheath of electrons, as in my linear density gradient paper.  For a linear density gradient in this sheath, you end up with an approximately constant offset in the electric field strength.


-----------------------------------------------------------------
### -10- Generic SLAC Data Analysis Code - '/SLACData/'

The SLAC data part is broken into two chapters.  This first chapter deals with more generic scripts that can be applied more easily to other SLAC analysis efforts.  The second chapter is a brief look at the specific scripts for analyzing the data from the E308 commissioning run.

FACETPlotter.py - This is a module that contains scripts for plotting data from DTOTR2.  In particular, the 'plotDTOTR2_withProjection' function makes the nice log-scale projection plots I used in my dissertation.  Also of note is the 'loadDTOTR2' and 'getZarr' functions, which load necessary data from the .mat matlab file.  These functions demonstrate the basic way that data is extracted from these .mat files, although from my experience they were nearly impossible to read in spyder.  To be able to explore these files on your own, definitely need to download a version of matlab and look at these .mat files through that.  Anyways, at the bottom of this module I have a main function for just quickly plotting a single DTOTR2 image.

MatrixCalc.py - This second module deals with calculating the resulting magnification from the imaging spectrometer.  I originally copied it from a FACET-I matlab file provided by Alex Knetsch.  I modified a lot of this to be used as 1D arrays to calculate the magnification vs the object plane, so some additional work may be needed to make this more general.  'getSpectrometerArrays' loads in the spectrometer values from the .mat file, using the tags 'LI20_LGPS_3141_BACT' 'LI20_LGPS_3261_BACT' and 'LI20_LGPS_3091_BACT' for q0, q1, and q2, respectively.  'Rmat' is a general generator for a 4x4 transfer matrix for either drift space or quadrupoles.  'FACET_Espec_Matrix' calculates the full transfer matrix for the FACET-II imaging spectrometer using the quadrupole values.  'getM11ForAllSteps' is the top-level function to just get the magnification M11 for each object plane location in the scan.


-----------------------------------------------------------------
### -11- Plasma Lens Commissioning Shift Analysis - '/SLACData/'

Lastly, I want to highlight some of the main scripts I used in analyzing my specific Commissioning shift data.  I can forsee some of these scripts requiring significant changes in order to be useful for other applications, but they might provide some ideas for analysis techniques.

AllDirLoop_AnalyzeChargeAndYbeam.py - Loops over all the directories, and calculates both the total amount of light on DTOTR2 for each image and the energy centroid location on the DTOTR2 screen.  Useful for comparing these two basic signals across datasets.

AllDirTopView.py - This one loads in TopView data from different datasets, and plots longitudinal lineouts along the beam axis.  There are flags at lines 83-86 to do various things with these lineouts, such as compare with and without beam, normalize to the backing pressure, smooth the data, and fit the data to a SuperGaussian function.

AllDirToroidLoop.py - Here I plot the toroid signals vs the backing pressure.  Toroids 'TORO_LI20_1988_0_TMIT','TORO_LI20_2040_0_TMIT', and 'TORO_LI20_2452_0_TMIT' are all far upstream in the linac and should show constant beam charge.  Toroid 'TORO_LI20_3163_0_TMIT' is just upstream of the PB, and toroid 'TORO_LI20_3255_0_TMIT' is downstream of the PB.

DTOTR2_ybeamcalibration.py - A simple version of AllDirLoop_AnalyzeChargeAndYbeam.py that only looks at the 0psi dataset and the energy centroid.  I take this y centroid vs object plane and fit a linear function to it, and this calibration is used in FACETPlotter to guess where the 10GeV line actually is.

SingleImage_TopView_NicePlots.py - This is a script to actually plot out TopView, and I overlay it with some details to make it look nice for my thesis.

SingleDirSigmaLoop_SmartCalc_NOGammaSlices.py - One of my main workhorses.  This one takes a single directory, loops over all DTOTR2 data taken, and fits the horizontal beam size projection to a Gaussian function.  During this process, I also filter the data by detecting when either the image is bad or the fit itself is bad.  These filtered images are discarded.  I also include the magnification from the imaging spectrometer to calculate the electron beam sizes at the object plane.  Afterwards, I make a plot of sigma vs object plane location and fit this curve to that of an ideal beam propagating in vacuum.

SingleDirSigmaLoop_SmartCalc_SetGammaSlices - Same as above, but I instead look at energy slices of the beam.  I iterate through bins of electron energies, and within those energy bins I do the same thing as above:  take the projection, fit a Gaussian sigma, plot versus object plane, and fit to vacuum betafunction propagation.  This process is a little tricky since the DTOTR2 screen is not linear with electron energy.  Overall this analysis didn't show too much since our beam wasn't strongly chirped or anything, but this might be useful in the future.

SingleDirSigmaLoop_EnergyContours.py - Another way of plotting the DTOTR2 data.  I take bins of electron energy and find what percentage of the beam is in each bin for every object plane location.  Then, I make a 2D graph of how these percentages change as the object plane is scanned.  It's a very strange type of plot but I like it!

Everything Else - There are a lot here that I didn't mention, and most of the other ones are variations on the 'SigmaLoop' formula.  The ones without 'SmartCalc' are before I started adding in the spectrometer settings, and also probably didn't even try to filter the data.  Then there are variations where I use the energy projection of DTOTR2 to get 'relative energy slices', but this is kind of odd because the energy centroid will change due to imaging conditions so this wasn't too useful.  I also had scripts to loop over all directories and compare, but those honestly take a lot of time to run.  At the end of it all, I mostly just used NOGammaSlices for when I wanted to calculate sigma for the full beam, and SetGammaslices for when I wanted to look at energy slices.
