# GISP2-Remelt-Statistical-Analysis---Honors-Thesis---Rylan-Abel

Data Analysis and Comparison of the New GISP2 Ice Core Data: Continuous High-Resolution Water Isotope Data Vs Discrete 1997 and 2004 Collection 

2026 Honors Thesis for The Department of Geologic Sciences at the University of Colorado at Boulder
By Rylan Abel

Thesis Advisor: 
Dr. Bradley Markle, The Department of Geologic Sciences

Defense Committee:
Dr. Bradley Markle, The Department of Geologic Sciences
Dr. Brian Hynek, The Department of Geologic Sciences
Dr. Tracy Farrell, The Program for Writing and Rhetoric

This repository contains the code used in Rylan Abel's honors thesis "Data Analysis and Comparison of the New GISP2 Ice Core Data: Continuous High-Resolution Water Isotope Data Vs Discrete 1997 and 2004 Collection". It includes a copy of all the code for the interval used within the thesis and any data used. The thesis itself oversees a statistical analysis between the original GISP2 analyses completed on IRMS and the INSTAAR stable isotope lab's new re-analysis of GISP2 using a CRDS-CFA system. The analysis of GISP2 for delta 18O was originally completed by P.M. Grootes and M. Stuvier in 2007 at the University of Washington's Quaternary Isotope Laboratory and the analysis of GISP2 for delta deuterium (D) was completed by J. White in 2004 at the University of Colorado Boulder’s Stable Isotope Laboratory. The majority of the data analysis was created in Python (Spyder (Python 3.13)), some visuals were edited within Fire Alpaca (64bit). 

GISP2_18O_Analysis_1370-1890_Honors_Thesis.py is the script that completes all analysis of delta 18O, including the downsampling of 2025 dataset and the statistical analysis between the 1997 and 2025 datasets.
GISP2_D_Analysis_1370-1890_Honors_Thesis.py is the script that completes all analysis of delta D, including the downsampling of 2025 dataset and the statistical analysis between the 2004 and 2025 datasets.
Statistical Analysis completed:
Lagging Mean - Takes a set interval of points and takes
    the avaerage, done with every point within the arrays. Smooth out
    short-term fluctuations and highlight longer-term trends or cycles.
Difference between datasets, Comparison test of original dataset, Creation of a CTL Mean dataset with and accompanying T-Test, Lag Correlation test with an accompanying re-analyzer that appplies lag to a given interval

"All data used (edited)" contains the version of the 1997/2004 GISP2 datasets that can be directly used within the .py files provided, with metadata clipped and correct column names for direct application into code. Additionally the section of the 2025 data produced at the INSAAR-SIL for this analysis is included within a zip file containing only the interval completed in time for the productiom of this thesis. This file is labeled New_Anlysis.zip while the 1997/2004 datasets are within the file Old_Analysis
"All data used (original)" contains the version of the 1997 and 2004 datasets unedited, thus including all meta data and essnetially appearing as they were originally downloaded. The citations for those datasets are as follows:
- P. M. Grootes, M. Stuiver. (1997). GISP2 Ice Core 110,000 Year Oxygen Isotope Data. NOAA Paleoclimatology Program, National Centers for Environmental Information (NCEI). Retrieved April 27, 2025. https://www.ncei.noaa.gov/access/paleo-search/study/17796, doi:10.25921/jtjy-9030
- White, J. (2004). GISP2 Stable Isotopes (Deuterium, Deuterium Excess, and Oxygen). Boulder, CO: National Center for Atmospheric Research, ARCSS Data Archive. Retrieved October 29, 2025.

Notes for using the files produced for the analysis of the entire 2025 re-analysis of GISP2
The intention of the code created for this project was for it to be quickly addaptable for use when the 2025 re-analysis is completed. There are only three additional changes that
