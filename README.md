# Starburst
Here is the GitHub repository for pyStarburst99. This is essentially a translation/slight reworking of the Starburst99 fortran code. 

## Installation
* Clone the repository to install
```
$ git clone https://github.com/CalumHawcroft/Starburst.git
```

* Or run in browser with Colab

https://colab.research.google.com/drive/1MN4P8Q47jshh3bw2sEODb3eJDNxgAbLc?usp=sharing

## Running the code
The easiest method is by following the link to Colab and running pyStarburst99 in your browser. In this case no installation is required, although the runtime will be longer than a local installation. The code itself will run in the same way locally after cloning the repository to your local machine and installing the prerequisite packages in your python environment.

After that, much like the original Starburst99 code, you need to define input parameters. This is no longer done with an input file, but instead defined at the start of the pySB99 file. As a disclaimer, currently not all functions and model options from Starburst99 (for example calculating strengths of spectral features such as the CO index) but there also new options available in pyStarburst99 like new metallicities and an increased upper mass limit. We will be adding functionality regularly so let us know if there is something you need that's missing and check for updates!

* **Total stellar mass** (M_total) - this is the total stellar mass (spread between the upper and lower mass limits with your chosen IMF) for an instantaneous burst in solar masses. Star formation histories and mixed age populations are planned for a future update, in the meantime the FiCUS code may be useful (https://github.com/asalda/FiCUS).

* **IMF exponent(s)** (IMF_exponents) [KROUPA = 1.3,2.3] - one or more IMF exponent for a power-law can be specified. The exponents refer to the individual power-law intervals, ordered by increasing mass. For instance, 1.3,2.3 specifies an IMF with exponents of 1.3 and 2.3 at low and high masses, respectively, with the boundaries given in the next input field. A single Salpeter-type law would be entered as 2.35. If there is more than one input value, the entries must be comma separated. (note that currently only up to 2 IMF exponents have been robustly verified in pySB99 but in principle more exponents should be fine following the SB99 code)

* **Mass boundaries for the IMF** (IMF_mass_limits) [KROUPA = 0.1,0.5,100 SOLAR MASSES] - the boundaries of the IMF intervals corresponding to the specified exponents. In this specific example we define two intervals from 0.1 to 0.5 and from 0.5 to 100 solar masses, The former would have a slope of 1.3, and the latter 2.3. A single power law between 1 and 100 solar masses would be entered as 1,100. The input values must be comma separated. For many metallicity options now the IMF upper mass limit can be increased to 300 (and up to 500 for a few)

* **Code Speed Option** (run_speed_mode) [DEFAULT, FAST, HIGH_RES] - in order to optimise the speed of the code for different purposes we include mulitple options which correspond to setting the resolution of the isochrone interpolation. 
In 'DEFAULT' mode the interpolation is increased in the mass range from 7 to 35 solar masses to ensure smooth outputs for a wide range of ages. 
In 'FAST' mode the interpolation is constant throughout the mass range, this might result in some numerical artifacts beyond 10Myr when lower mass stars are more important for the outputs than at early times. 
In 'HIGH_RES' mode the interpolation is significantly increased across all masses, and even moreso in the 7 to 35 solar mass range, this ensures smooth outputs at the cost of runtime.
Keen users can also go into the 'interpolate_param' function and set their own custom interpolation resolution to highlight the stellar mass range they are most interested in.



