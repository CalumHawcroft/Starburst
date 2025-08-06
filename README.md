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

* **Total stellar mass (M_total)** - This is the total stellar mass (spread between the upper and lower mass limits with your chosen IMF) for an instantaneous burst in solar masses. Star formation histories and mixed age populations are planned for a future update, in the meantime the FiCUS code may be useful (https://github.com/asalda/FiCUS).




