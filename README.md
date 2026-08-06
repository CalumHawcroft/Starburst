# Starburst
Here is the GitHub repository for pyStarburst99. This is essentially a translation/slight reworking of the Starburst99 fortran code. Please let us know if you have any thoughts or questions on the code at pystarburst99@gmail.com

Last updated 4th August 2026.

## Installation
* Clone the repository to install
```
$ git clone https://github.com/CalumHawcroft/Starburst.git
```

* Or run in browser with Colab

https://colab.research.google.com/drive/1MN4P8Q47jshh3bw2sEODb3eJDNxgAbLc?usp=sharing

## Running the code
The easiest method is by following the link to Colab and running pyStarburst99 in your browser. In this case no installation is required, although the runtime will be longer than a local installation. The code itself will run in the same way locally after cloning the repository to your local machine and installing the prerequisite packages in your python environment.

After that, much like the original Starburst99 code, you need to define input parameters. 
This can either be done through an input file, or instead defined at the start of the pySB99 code file itself. If using an input file you will still have to change the input_option variable at the beginning of the pySB99 file to 'INPUT_FILE'. 
As a disclaimer, currently not all functions and model options from Starburst99 (for example calculating strengths of spectral features such as the CO index) but there also new options available in pyStarburst99 like new metallicities and an increased upper mass limit. We will be adding functionality regularly so let us know if there is something you need that's missing and check for updates!

* **Input from code or file** (input_option) [MAIN_CODE, INPUT_FILE] - The MAIN_CODE option will let you define all variables in the code itself (maybe preferred for running the notebook or using colab), while the INPUT_FILE option lets you define all variables in a separate input file which is read in by the code. An example input file is provided 'pysb_inputs.py', note that this should be a '.py' file so the variables are read as intended, but it is essentially a text file.

* **Star formation** (star_formation_option) [ISF, CSF] - can either run a model of an instantaneous burst of star formation 'ISF' or for continuous star formation 'CSF'. If using instant star formation you have to define the total mass of the burst (in solar masses). If using continuous star formation you have to define the star formation rate (in solar masses per year).

* **Total stellar mass** (M_total) - this is the total stellar mass (spread between the upper and lower mass limits with your chosen IMF) for an instantaneous burst in solar masses.

* **Star formation rate** (SFR) - this is the star formation rate (spread between the upper and lower mass limits with your chosen IMF) in solar masses per year. Variable star formation histories and mixed age populations may be included in a future update, in the meantime the FiCUS code may be useful (https://github.com/asalda/FiCUS).

* **IMF exponent(s)** (IMF_exponents) [KROUPA = 1.3,2.3] - one or more IMF exponent for a power-law can be specified. The exponents refer to the individual power-law intervals, ordered by increasing mass. For instance, 1.3,2.3 specifies an IMF with exponents of 1.3 and 2.3 at low and high masses, respectively, with the boundaries given in the next input field. A single Salpeter-type law would be entered as 2.35. If there is more than one input value, the entries must be comma separated. (note that currently only up to 2 IMF exponents have been robustly verified in pySB99 but in principle more exponents should be fine following the SB99 code)

* **Mass boundaries for the IMF** (IMF_mass_limits) [KROUPA = 0.1,0.5,100 SOLAR MASSES] - the boundaries of the IMF intervals corresponding to the specified exponents. In this specific example we define two intervals from 0.1 to 0.5 and from 0.5 to 100 solar masses, The former would have a slope of 1.3, and the latter 2.3. A single power law between 1 and 100 solar masses would be entered as 1,100. The input values must be comma separated. For many metallicity options now the IMF upper mass limit can be increased to 300 (and up to 500 for a few)

* **Output age** (output_age) - this has no impact on how the model is produced, it is just to set the age of an example population which is taken from the outputs and plotted for visual inspection. The age is defined in Myrs.

* **Code Speed Option** (run_speed_mode) [DEFAULT, FAST, HIGH_RES] - in order to optimise the speed of the code for different purposes we include mulitple options which correspond to setting the resolution of the isochrone interpolation. 
In 'DEFAULT' mode the interpolation is increased in the mass range from 7 to 35 solar masses to ensure smooth outputs for a wide range of ages. 
In 'FAST' mode the interpolation is constant throughout the mass range, this might result in some numerical artifacts beyond 10Myr when lower mass stars are more important for the outputs than at early times. 
In 'HIGH_RES' mode the interpolation is significantly increased across all masses, this ensures smooth outputs at the cost of runtime.
Keen users can also go into the 'interpolate_param' function and set their own custom interpolation resolution to highlight the stellar mass range they are most interested in.

* **Metallicity** (Z) [MWC, MW, LMC, SMC, IZw18, XMP, Z0] - this selects the metallicity of the stellar evolutionary tracks and spectral library. The current options are set to reflect the applicable stellar environments the GENEC evolutionary tracks are tailored to; Milky Way Centre (MWC - Z=0.02), Milky Way (MW - Z=0.014), LMC (Z=0.006), SMC (Z=0.002), IZw18 (Z=0.0004), XMP (Z=1e-5) Z0 (Z=0). More options for other evolutionary tracks are available (for now the MESA models based on Pauli+2026) we plan to include more models in future releases.

* **Spectral libraries** (SED_library [FW, WM, PoWR] and spectral_library [WM, PoWR]) - pySB99 contains two separate OB stellar spectral libraries for different output purposes. First is the low-resolution but wide wavelength coverage SED library for computing quantities such as ionising fluxes. Second is the high-resolution but narrower wavelength coverage spectral line library, used to create composite diagnostic line profiles. 
For the lower resolution OB star stellar library used to generate SEDs (equivalent to .spectrum SB99 output), there is a choice between WMBasic (WM), Fastwind (FW) and PoWR (PoWR) libraries. 
For the higher resolution OB star spectral library used for line synthesis, the WMBasic and PoWR libraries are currently available, with Fastwind libraries planned for a future update. 
The WMBasic libraries are described in Leitherer+2010, the Fastwind library is designed to match the latest GENEC metallicity and mass range for pySB99 (see Hawcroft+2025), and the PoWR libraries (e.g. Hainich+2019) are available at https://www.astro.physik.uni-potsdam.de/PoWR/powrgrid1.php
The WMbasic code is no longer maintained and so may be limited in updates to input physics and atomic data compared to others. However, the WMbasic library is still advantageous in some ways.
For the low resolution SED grid, I would advise to use either the Fastwind or PoWR grids. In addition to code maintenance, they also offer a better match to the metallicities of the evolutionary models, as well as better parameter space coverage in temperature and luminosity/surface gravity for OB stars and the latest mass-loss rate prescriptions. For models at Z=0.02, 0.00001 and 0.0, and/or models which include stars above 100 solar masses I would recommend the Fastwind library as the PoWR library does not yet have coverage for these inputs.
For the high resolution spectra grid, the PoWR grid is likely the better choice, although it does not include X-rays so you might want to check the very high ionization lines with WMbasic spectra. 
For now, all the high resolution libraries are limited for very low metallicities and very high stellar masses, but hopefully this will be remedied soon with the inclusion of Fastwind spectra. At first, the Fastwind spectra will only be available for a few key wind and photospheric features and only as normalized spectra, but with updates to the Fastwind code a future library will be able to include the full stellar continuum at physical flux levels.
We also plan to include new empirical libraries in the future.

* **Stellar evolutionary models** (track_choice) [GENEC, Pauli] - choice of stellar evolutionary model grid used to generate the isochrones and resulting stellar population parameters. The current options are the GENEC grids or the MESA models based on Pauli+2026.
I would suggest consulting the papers for the respective models to make your input choice as both are very good, valid options. However, the Pauli grids do not yet cover the lowest/highest metallicities (Z=0.02, 0.00001 and 0.0), while the GENEC grids do not include stars above 100 solar masses at SMC and IZw18 metallicities (Z=0.002 and 0.0004).
These are both grids of single star models, GENEC includes the option of rotation or no rotation while the Pauli models have a rotation of 100 km/s which is based on the mean rotation for observed OB populations in the LMC and SMC.
Another key factor is the mass-loss prescription. The GENEC models use the Vink+2001 predictions for typical OB stars. The Pauli models use the Vink+2001 prescription scaled down by a factor of 3, in better agreement with empirical results (and closer to the Bjorklund+ or Krticka+ mass-loss theoretical predictions). The Pauli models also include Eddington-limit induced mass ejections which may be more representative of the distribution of stellar types in observed populations. 
Unfortunately no binaries for now.

* **Rotation option** (rot) [True, False] - choice of rotating (v=0.4 break-up) or non-rotating stellar evolutionary tracks, True includes rotation while False is non-rotating. This is only relevant if using the GENEC libraries, the Pauli grid only has one option for rotation as defined above.

* **Output plot options** (*_ion_flux, _wind, _uv_slope, _ew, _colours, _sed, _hires_spectra, _SN_rate, _isochrones) - choice of which output arrays and files to plot with time. Corresponding to ionising fluxes of HI, HeI and HeII, the wind momentum/power, beta uv-slope, equivalent widths of Halpha, Hbeta etc, colours (V, U, I, B, M_V) and an example SED at a specified age. We have recently added the predicted UV spectra from high resolution theoretical spectral libraries, as well as the predicted rates of SN and you can inspect/save the isochrones. These are the current outputs which are fully tested and verified for the latest release. More to come soon including a detailed breakdown of spectral types.

* **End age of the population** (times_spectra_end) [Myr] - choice of how long to output the results for the model population, typically somewhere around 50Myr is used for the instant burst models, you can make outputs for late times but this code is optimised for synthesizing populations of massive stars which makes it maybe not the best choice for older populations.

* **IMF input from file** (input_imf_option) [True, False] - currently testing the option to generate a population using an IMF generated by other programs/codes instead of the IMF defined in pySB99. The easiest outputs we have tested so far generate a long list of masses for each individual star in a population (e.g. from stochastic IMF sampling) and so if this option is set to True that should be the format of the input file. If setting this option to False the IMF will be defined as described above.

* **Output save option** (save_output) [True, False] - choice of whether to save the array files used to generate the outputs described above (note that the plot option must be set to True as well for the save option to work). The majority of the outputs are saved as text files with the SEDs saved as numpy arrays instead.

* **Model designation** (SBmodel_name) [self-defined string) - if save output is set to True you also have to provide a name for the output model, all outputs will be saved to a folder with this designation. 

* **Maximum supernova mass** (maximum_SN_mass) [Msol] - if the SN rate output is selected you can define the maximum mass of star which will be counted as undergoing a supernova, this could be the same as the maximum stellar mass of the input population or it could be somewhere around 30 or 40 solar masses, which is maybe more physically motivated. (this can also be somewhat dependent on the metallicity)
