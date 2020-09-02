# Data Characterization of COVID-19 Pandemic and Lifting Confinement in Mexico
*Cristy Leonor Azanza Ricardo and Esteban A. Hernandez-Vargas*

Raw data is obtined from the [Government] official site. To process this data filenames were standardized and stored in a Google Drive [folder]. 

Some data files are supplied:
- **20200525fits.csv**: Contains all the fits results obtained with the data up to June 25, 2020. 
- **20200525_ResXStates.csv**: The same as above with the total infected individual for each state. This datasheet is used for the map generation
- **20200525_Start1503.csv**: Data up to May 25, 2020.
- **20200815_Start1503.csv**: Data up to August 15, 2020.
- *20200526_LAST*: This is a folder that contains results from simultions in the re-smpling strategy for all states

Some scripts and notebooks are available:
- **dataExtract.py**: Extract required information from raw data and generates Figures 5, 9 and 10. Raw data must be downloaded and unziped in the folder *DirNacEpid*.
- **genStat.py**: Performs the re-sampling strategy based on the distribution functions in Figure 1, for all states generating files as the ones in *20200526_LAST*.
- **genFig02.py**: Generates Figure 2.
- **genFig03.py**: Generates Figure 3.
- **genFig04.py**: Generates Figure 4.
- **genFig06-07.py**: Generates Figures 6 and 7.
- **genFig08.py**: Generates Figure 8.
- **genFig11-12.ipynb**: Generates Figures 11 and 12.
- **importCOVID19.py**: Includes models for COVID19 pandemic characterization in Mexico.
- **PDEparams.py**: Differential Evolution algorithm python library

[Government]: https://www.gob.mx/salud/documentos/datos-abiertos-152127

[folder]: https://drive.google.com/file/d/1gGAdd-YjhH7dD2pGTe5C290VR9nAmc-9/view?usp=sharing
