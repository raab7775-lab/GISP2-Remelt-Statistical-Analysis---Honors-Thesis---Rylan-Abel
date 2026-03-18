# GISP2-Remelt-Statistical-Analysis---Honors-Thesis---Rylan-Abel

## Thesis Information

Data Analysis and Comparison of the New GISP2 Ice Core Data: Continuous High-Resolution Water Isotope Data Vs Discrete 1997 and 2004 Collection 

2026 Honors Thesis for The Department of Geologic Sciences at the University of Colorado at Boulder
By Rylan Abel

Thesis Advisor: 
Dr. Bradley Markle, The Department of Geologic Sciences

Defense Committee:

Dr. Bradley Markle, The Department of Geologic Sciences

Dr. Brian Hynek, The Department of Geologic Sciences

Dr. Tracy Farrell, The Program for Writing and Rhetoric



Thesis Abstract:
The Greenland Ice Sheet Project (GISP2) is an ice core drilled within the Greenland Ice Sheet that has produced several foundational datasets for the study of paleoclimate. The core’s measured stable isotope record of δ18O and δD is one of these foundational datasets, completed in 1997 for the measurement of δ18O and in 2004 for δD, using Isotope Ratio Mass Spectroscopy. Beginning in 2024, there has been an effort to reanalyze the GISP2 ice core for δ18O and δD, completed by the INSTAAR Stable Isotope Lab with Laser Absorption Spectroscopy as a continuous flow analysis. An analysis was completed that aimed to downsample the 2025 dataset to the scale of the older datasets and to complete a statistical analysis between the old and new datasets. This was accomplished with the intent of discovering any apparent difference in the core, and the results of this are the detection of two distinct differences between the old and new datasets. One is the statistical offset between the means of the new and old data. Within the δ18O measurement campaign, where the 2025 dataset is isotopically heavier than the 1997 dataset, and for the δD measurement campaign, the 2025 dataset is isotopically lighter than the 2004 dataset. The second difference is an offset within the depth scale of the core, where the measured values of the 2025 dataset match most closely with the 1997 dataset when every point of the 2025 dataset is moved 0.02 m to a deeper depth. These differences were speculated to be caused in part by minimal core alteration of 30 years of storage, as well as systematic differences introduced between the operations of each melt. Ultimately, these are within the error of the δ18O and δD systems, meaning they are minimal enough to declare that the two datasets can essentially be treated as equivalent in future analysis. Thus, a core that has been in storage for 30 years can successfully be remelted on a completely different LAS system and can still come to the same results as the IRMS system.

## Purpose of Repository and Citations

This repository contains the code used in Rylan Abel's honors thesis, "Data Analysis and Comparison of the New GISP2 Ice Core Data:

Continuous High-Resolution Water Isotope Data Vs Discrete 1997 and 2004 Collection". It includes a copy of all the code for the interval used within the thesis and any data used. The thesis itself oversees a statistical analysis between the original GISP2 analyses completed on IRMS and the INSTAAR stable isotope lab's new re-analysis of GISP2 using a CRDS-CFA system. The analysis of GISP2 for delta 18O was originally completed by P.M. Grootes and M. Stuvier in 2007 at the University of Washington's Quaternary Isotope Laboratory, and the analysis of GISP2 for delta deuterium (D) was completed by J. White in 2004 at the University of Colorado Boulder’s Stable Isotope Laboratory. The majority of the data analysis was created in Python (Spyder (Python 3.13)), and some visuals were edited within Fire Alpaca (64bit). 

GISP2 Original Site: a latitude of 72.5833333 and a longitude of -38.466667.

Core stored within the National Science Foundation Ice Core Facility for 30 years.

The citations for the original 1997 and 2004 datasets:

- P. M. Grootes, M. Stuiver. (1997). GISP2 Ice Core 110,000 Year Oxygen Isotope Data. NOAA Paleoclimatology Program, National Centers for Environmental Information (NCEI). Retrieved April 27, 2025. https://www.ncei.noaa.gov/access/paleo-search/study/17796, doi:10.25921/jtjy-9030
- White, J. (2004). GISP2 Stable Isotopes (Deuterium, Deuterium Excess, and Oxygen). Boulder, CO: National Center for Atmospheric Research, ARCSS Data Archive. Retrieved October 29, 2025.

## Contents of Repository

GISP2_18O_Analysis_1370-1890_Honors_Thesis.py is the script that completes all analyses of delta 18O, including the downsampling of the 2025 dataset and the statistical analysis between the 1997 and 2025 datasets.

GISP2_D_Analysis_1370-1890_Honors_Thesis.py is the script that completes all analyses of delta D, including the downsampling of the 2025 dataset and the statistical analysis between the 2004 and 2025 datasets.

"All data used (edited)" contains the version of the 1997/2004 GISP2 datasets that can be directly used within the .py files provided, with metadata clipped and correct column names for direct application into code. Additionally, the section of the 2025 data produced at the INSAAR-SIL for this analysis is included within a zip file containing only the interval completed in time for the production of this thesis. This file is labeled New_Analysis.zip, while the 1997/2004 datasets are within the file Old_Analysis

"All data used (original)" contains the version of the 1997 and 2004 datasets unedited, thus including all metadata and essentially appearing as they were upon originally being downloaded. 

Downsampling for both delta 18O and D:

- Completed in one part for delta 18O.
- completed in two parts for delta D due to split in dataset resolution, and the combined in the Interval_complete(...) function.

Statistical Analysis completed for both delta 18O and D:

- Lagging Mean: Takes a set interval and takes the mean of it, done with every point within the arrays. Smooth out short-term fluctuations and highlight longer-term trends or cycles.
- Difference between datasets: Subtracts the new Data from the old dataset. Takes the mean and Standard Deviation of the difference.
- Comparison test: tests the dataset's normality with the Shapiro-Wilk test, then applies a parametric or non-parametric test to the new and old datasets, depending on the Shapiro results. Completed to determine if there is any statistically significant difference between the analyses.
- Creation of a CTL Mean dataset with an accompanying T-Test: Creates a normal distribution from both the old and new datasets using the Central Limit Theorem and takes a t-test of those normal datasets. Completed to test of there is any statistically significant difference between the analyses again.
- Lag Correlation test with an accompanying re-analyzer that applies lag to a given interval: moves the 2025 dataset up or down compared to the 1997/2004 dataset and completes a correlation test. Tests for any apparent depth offset in the core, due to alteration or not.

Notes for using the files produced for the analysis of the entire 2025 re-analysis of GISP2
The intention of the code created for this project was for it to be quickly adaptable for use when the 2025 re-analysis is completed. There are only three additional changes to the base code provided to make it ready for use in with the full dataset.
- Change the pd.read_csv(...) to the complete files of the 2025 dataset only
- **Change the Testing_Restabilization (...) function to include the overall interval of the dataset.** The current code contains the 1370-1890 m interval for convenience of analysis, and for analysis with the whole interval, this needs to be adjusted. Done easily by going to the bottom of the function and changing the start and end interval of the bottom four lines.
- Due to several gaps being present at the top of the core, it is advised to make a function within the Interval functions to account for the number of NaNs present within an interval that disregards any intervals with more than 20% NaNs as this might skew the analysis. Not a particular problem for the interval studied, which had more consistent archival pieces.
