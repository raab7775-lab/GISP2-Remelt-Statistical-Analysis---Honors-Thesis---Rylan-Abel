# -*- coding: utf-8 -*-
"""
Created on Thu May 22 13:51:43 2025

@author: Rylan Abel
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
from scipy.stats import shapiro
from scipy.stats import spearmanr
import random

GISP97 = pd.read_csv(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\Old_Analysis\1997_GISP2D_18O.DAT")
GISP97 = pd.read_table(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\Old_Analysis\1997_GISP2D_18O.DAT", sep='\t', skiprows=[0])
depth25 = np.loadtxt(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\New_Melt\2025_depth_model.txt")
O1825 = np.loadtxt(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\New_Melt\2025_data_d18O_SIL_cfa.txt")
ModGISP97 = GISP97[GISP97["(per mil)"] != 999999]  # Gets rid of unvalid data points
# Interval Given is between 1372 and 1890

depth97 = np.array(ModGISP97["Top (m)"])
O1897 = np.array(ModGISP97["(per mil)"])
# Takes pandas column and turns it into a numpy array that can be in function


def Data_Unedited_Figure():
    '''
    Function that outputs a figure of both datasets as it is, without
    any editing. No input, output is a single figure.'

    Input
    None.

    Returns
    None.

    '''
    plt.figure(figsize=(10, 5))
    plt.plot(depth25, O1825,
             label="δ180 2025 data")
    plt.plot(depth97, O1897,
             label="δ180 1997 data", color="orange")
    plt.xlabel("Depth (m)")
    plt.ylabel("δ18O per mL")
    plt.title("Comparison Between 2025 and 1997 GISP2 Data")
    plt.legend()
    plt.show()


def Interval(a, b, Figure=False):
    '''
    Function that performs downsampling of 2025 data to the 1997 resolution

    Parameters
    ----------
    a : int
        Starting point of interval
    b : int
        Ending point of interval
    Figure : bool, optional
        Input True for a figure to be printed of both same resolution samples
        in a single graph in the selected interval.
        The default is False.

    Returns
    -------
    3 Arrays
    DEPTH97Correct : np.float64
        Output is the corrected depth array.
    O1897int : float
        The corresponding δ18O values from the 1997 dataset.
    Average25int : np.float64
        The average of the 2025 δ18O values corresponding to the 1997
        depth resolution

    '''
    DEPTH25int = []
    DEPTH97int = []

    for i in depth25:
        if a <= i <= b:
            DEPTH25int.append(i)
    for i in depth97:
        if a <= i <= b:
            DEPTH97int.append(i)
    start25 = np.where(depth25 == DEPTH25int[0])[0][0]
    end25 = np.where(depth25 == DEPTH25int[-1])[0][0]
    start97 = np.where(depth97 == DEPTH97int[0])[0][0]
    end97 = np.where(depth97 == DEPTH97int[-1])[0][0]
    O1825int = O1825[start25:end25]
    O1897int = O1897[start97:end97]

    def Average(b, c):
        Interval_sum = 0
        start = int(DEPTH97int[b])
        end = int(DEPTH97int[c])
        length = []

        for i in DEPTH25int:
            if start <= i <= end:
                Interval_sum += i
                length.append(i)
        if not length:
            return None
        d = np.where(DEPTH25int == length[0])[0][0]
        e = np.where(DEPTH25int == length[-1])[0][0] + 1
        Interval_O18 = O1825int[d:e]
        Interval_sum25 = 0
        for i in Interval_O18:
            Interval_sum25 += i
        Interval_length25 = len(Interval_O18)
        Average_int25 = Interval_sum25 / Interval_length25
        return (Average_int25)

    Average25int = []

    for a in range(len(DEPTH97int)-1):
        start = a
        end = a + 1
        at = Average(start, end)
        Average25int.append(at)

    DEPTH97Correct = []
    for i in DEPTH97int:
        DEPTH97Correct.append(i + 1)

    if Figure is False:
        print()
    elif Figure is True:
        plt.figure()
        plt.plot(DEPTH97Correct[:-1], O1897int, label="δ18O 1997 data",
                 color='orange')
        plt.plot(DEPTH97Correct[:-1], Average25int,
                 label="Average δ18O from 2025 Data")
        plt.xlabel("Depth (m)")
        plt.ylabel("δ18O‰")
        plt.title("Same-Resolution Comparison of δ18O Between 2025 and 1997 GISP2 Data")
        plt.legend()
        plt.show()

    return (DEPTH97Correct, O1897int, Average25int)


def Moving_Mean(a, b, c, Figure=False):
    '''
    Moving mean function applied to a same resolution interval of the 1997
    and the 2025 data. Takes an interval of points (c) and takes
    the avaerage, done with every point within the arrays. Smooth out
    short-term fluctuations and highlight longer-term trends or cycles.

    Parameters
    ----------
    a : Int
        Input the start interval
    b : Int
        Input the end interval
    c : Int
        Input the width of the smoothing interval
    Figure : bool, optional
        Input true for a figure with a graph of both the 1997 and 2025 datasets
        smoothed by the established smoothing interval.
        The default is False.

    Returns
    -------
    2 Arrays
    Mean97 : float
        An array of the 1997 data smoothed by the smoothing interval
    Mean25 : float
        An array of the 2025 data smoothed by the smoothing interval

    '''
    if c == 1:
        New = 1
    elif c > 1:
        New = int(c/2)

    DEPTH97int, O1897int, Average25int = Interval((a - c), (b + c))

    def Interval_97(a, New):
        start = a
        Depthmeanint = []
        for i in DEPTH97int:
            if start - New <= i <= start + New:
                Depthmeanint.append(i)
        e = np.where(DEPTH97int == Depthmeanint[0])[0][0]
        f = np.where(DEPTH97int == Depthmeanint[-1])[0][0] + 1
        intO1897 = O1897int[e:f]
        MEAN97 = np.mean(intO1897)
        return (MEAN97)

    def Interval_av(a, New):
        start = a
        Depthmeanint = []
        for i in DEPTH97int:
            if start - New <= i <= start + New:
                Depthmeanint.append(i)
        e = np.where(DEPTH97int == Depthmeanint[0])[0][0]
        f = np.where(DEPTH97int == Depthmeanint[-1])[0][0] + 1
        intO18Average = Average25int[e:f]
        MEAN25 = np.mean(intO18Average)
        return (MEAN25)

    MEAN97 = []
    MEAN25 = []
    for i in range(len(DEPTH97int)):
        start = int(DEPTH97int[i])
        at = Interval_97(start, New)
        ab = Interval_av(start, New)
        MEAN97.append(at)
        MEAN25.append(ab)

    title = 'Running Mean', c,
    'Interval, Comparison of δ18O Between 2025 and 1997 GISP2 Data'

    if Figure is False:
        print()
    elif Figure is True:
        plt.figure()
        plt.plot(DEPTH97int[New:-New], MEAN97[New:-New],
                 label="Smooth Mean of δ18O 1997 data", color='orange')
        plt.plot(DEPTH97int[New:-New], MEAN25[New:-New],
                 label="Smoothed Mean of Average δ18O from 2025 Data")
        plt.xlabel("Depth (m)")
        plt.ylabel("Smoothed Mean of δ18O Data")
        plt.title(title)
        plt.legend()
        plt.show()

    return (MEAN97, MEAN25)


def Difference(a, b, Hist=False, Figure=False):
    '''
    Fuction that calculates the difference between the 1997 and 2025 data with
    depth

    Parameters
    ----------
    a : int
        Input starting interval
    b : int
        Input end interval
    Hist : bool, optional
        Outputs a histogram figure of dfference values. The default is False.
    Figure : bool, optional
        Outputs a figure of the difference array ploted with depth.
        The default is False.

    Returns
    -------
    DifferenceARR : float
        Array of difference values.

    '''
    DEPTH97int, O1897int, Average25int = Interval(a, b)
    Difference = Average25int - O1897int
    DifferenceARR = np.array(Difference)

    Difference_mean = np.mean(DifferenceARR)
    Difference_STD = np.std(DifferenceARR)

    print("Mean Difference:", Difference_mean, "Difference STD:",
          Difference_STD)

    if Hist is False:
        print()
    elif Hist is True:
        plt.figure()
        plt.hist(Difference, bins=30)
        plt.title("Distribution of Difference of δ18O with depth between 1997 and 2025 Measurements")
        plt.xlabel("Difference in δ18O‰")
        plt.ylabel("Frequency")
        plt.show()

    if Figure is False:
        print()
    elif Figure is True:
        plt.figure()
        plt.plot(DEPTH97int[:-1], DifferenceARR)
        plt.xlabel("Depth (m)")
        plt.ylabel("Difference in δ18O‰")
        plt.title("Difference between 1997 and 2025 δ18O Measurements with Depth")
        plt.show()

    return (DifferenceARR)


def Comparison_Tests(a, b, Figure=False):
    '''
    Tests for the normality of the 1997 and 2025 data set, then applys a
    subsequent parametric or nonparametric camparison based on the normality.
    Output is two printed statments, the first being the results of the
    shapiro test (>0.05 is normal), and the second being the results of
    the parametric or nonparametric tests (<0.05 is significant difference).

    Parameters
    ----------
    a : int
        Input starting interval.
    b : int
        Input end interval
    Figure : bool, optional
        A figure of the resulting comparison test after testing normality,
        Contains the histogram distribution of the 1997 data, the 2025 data
        and both dataset overlaid. The default is False.

    Returns
    -------
    None

    '''
    DEPTH97int, O1897int, Average25int = Interval(a, b)
    statistic_normal97, p_value_normal97 = shapiro(O1897int)
    statistic_normal25, p_value_normal25 = shapiro(Average25int)

    print("1997 δ18O p-value is:", p_value_normal97,
          "and 2025 p-value is:", p_value_normal25)

    if p_value_normal25 and p_value_normal97 >= 0.05:
        t_statisticPairedT, p_valuePairedT = stats.ttest_rel(O1897int,
                                                             Average25int)
        t_statisticIndT, p_valueiIndT = stats.ttest_ind(O1897int,
                                                        Average25int,
                                                        equal_var=True)
        Caption = ("Parametric test used, paired p-value is", p_valuePairedT,
                   "and independent p-value is", p_valueiIndT)
        print(Caption)
    elif p_value_normal25 or p_value_normal97 < 0.05:
        t_statisticpairW, p_valuepairW = wilcoxon(O1897int, Average25int)
        t_statisticindM, p_valueindM = mannwhitneyu(O1897int,
                                                    Average25int)
        Caption = ("Nonparametric test used, Paired p-value is", p_valuepairW,
                   "and independent p-value is", p_valueindM)
        print(Caption)

    if Figure is False:
        print()
    elif Figure is True:
        fig, axs = plt.subplots(2, 2, figsize=(5, 5))
        plt.subplots_adjust(wspace=0.5, hspace=1)
        plt.title("Distribiution of δ18O Values within 1997 and 2025 dataset")
        axs[0, 0].hist(O1897int, bins=50, color="orange")
        axs[0, 0].set_xlabel("δ18O‰ from 1997 dataset")
        axs[0, 0].set_ylabel("Frequency")
        axs[1, 0].hist(Average25int, bins=50, color="blue")
        axs[1, 0].set_xlabel("δ18O‰ from 2025 dataset")
        axs[1, 0].set_ylabel("Frequency")
        axs[0, 1].hist(O1897int, bins=50, color="orange")
        axs[0, 1].hist(Average25int, bins=50, color="blue", alpha=0.5)
        axs[0, 1].set_xlabel("δ18O‰ from 2025 and 1997 dataset")
        axs[0, 1].set_ylabel("Frequency")


def Testing_Normality(a, b, c, p_value=False, Histograms=False):
    '''
    To compare both data sets with a t-test, central limit theorum is used.
    It forms a normal distribution of the means of both data by taking a
    random sample of an established interval size and takeing the mean. This
    is accomplished with both datasets by using indexing to take from the same
    depth between datasets. By sampling several times, a normal distribution
    is taken.

    Parameters
    ----------
    a : int
        Input the number of trials accomplished.
    b : int
        Input the sample size of each trail.
    c : int
        Set a seed number. Thesis accomplished with seed number 42.
    p_value : bool, optional
         Prints the shapiro test p-value results for the new CLT distribution
         to prove normality of new distribution. (>0.05 is normal)
         The default is False.
    Histograms : bool, optional
        Outputs a figure of the two produced CLT histograms side by side for
        visual comparison of dataset. The default is False.

    Returns
    -------
    2 Arrays
    Normal_distribution25 : float
        The CLT dataset produced for the 2025 dataset.
    Normal_distribution97 : float
        The CLT dataset produced for the 1997 dataset.

    '''
    DEPTH97int, O1897int, Average25int = Interval(1370, 1890)
    O1897int = O1897int[:-1]
    random_index = np.arange(0, len(O1897int))
    Normal_distribution97 = []
    Normal_distribution25 = []
    random.seed(c)
    for i in range(a):
        random_sample = random.choices(random_index, k=b)
        O1897intSAMP = []
        Average25intSAMP = []
        for x in random_sample:
            O1897intSAMP.append(O1897int[x])
            Average25intSAMP.append(Average25int[x])
        Mean97 = np.mean(O1897intSAMP)
        Mean25 = np.mean(Average25intSAMP)
        Normal_distribution25.append(Mean25)
        Normal_distribution97.append(Mean97)

    statistic_normal25, p_value_normal25 = shapiro(Normal_distribution25)
    statistic_normal97, p_value_normal97 = shapiro(Normal_distribution97)

    if p_value is False:
        print()
    elif p_value is True:
        print("1997 p-value is", p_value_normal97, "and 2025 p-value is",
              p_value_normal25)

    if Histograms is False:
        print()
    elif Histograms is True:
        caption = ("1997:", round(p_value_normal97, 5), "2025:",
                   round(p_value_normal25, 5))
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        plt.subplots_adjust(wspace=0.5, hspace=1)
        fig.suptitle("Testing for Normality in Data")
        axs[0].hist(Normal_distribution97, bins=50, color='orange')
        axs[0].set_ylabel("Frequency")
        axs[0].set_xlabel("δ18O Mean Distribution 1997 dataset")
        axs[0].set_title("Mean δ18O Data Distribution 1997 dataset")
        axs[1].hist(Normal_distribution25, bins=50, color='blue')
        axs[1].set_title("Mean δ18O Data Distribution 2025 dataset")
        axs[1].set_ylabel("Frequency")
        axs[1].set_xlabel("Mean δ18O Data Distribution 2025 dataset")
        plt.figtext(0.43, 0, caption)

    return (Normal_distribution25, Normal_distribution97)


def TTest(a, b, c, Histogram=False, Boxplot=False):
    '''
    Function that performs a T test on the normal distributions formed
    by the Testing Normality function. Tests both a paired and independent
    t-test. Output is a printed statement of the t-test results.
    (<0.05 is significant difference)

    Parameters
    ----------
    a : int
        Input the number of trials accomplished.
    b : int
        Input the sample size of each trail.
    c : int
        Set a seed number. Thesis accomplished with seed number 42.
    Histogram : bool, optional
        Outputs a histogram of the 1997 and 2025 CLT dataset histograms
        overlaid. The default is False.
    Boxplot : bool, optional
        Outputs a graph of the 1997 and 2025 CLT datasets as two
        side by side boxplots. The default is False.

    Returns
    -------
    None.

    '''
    Normal_distribution25, Normal_distribution97 = Testing_Normality(a, b, c)
    t_statisticPaired, p_valuePaired = stats.ttest_rel(Normal_distribution97,
                                                       Normal_distribution25)
    t_statisticInd, p_valueiInd = stats.ttest_ind(Normal_distribution97,
                                                  Normal_distribution25,
                                                  equal_var=True)

    Caption = ("Paired p-value is", float(p_valuePaired),
               "and independent p-value is", float(p_valueiInd))
    print(Caption)
    if Histogram is False:
        return None
    elif Histogram is True:
        plt.figure()
        plt.hist(Normal_distribution97, bins=50, color='orange')
        plt.hist(Normal_distribution25, bins=50, color='blue', alpha=0.5)
        plt.ylabel("Frequency")
        plt.xlabel("δ18O‰")
        plt.title("Overlay of 1997 and 2025 Mean Distributions")
        plt.show()

    if Boxplot is False:
        return None
    elif Boxplot is True:
        plt.figure()
        plt.boxplot([Normal_distribution97, Normal_distribution25],
                    labels=['1997 Data', '2025 Data'])
        plt.ylabel("δ18O‰")
        plt.title("1997 and 2025 Mean Distributions Boxplots")
        plt.show()


def TTest_GorL(a, b, c):
    '''
    Function that determines which data set is more positive or negative
    to the other, completed through greater or less than t-tests. Output is 4
    printed statements giving results of t tests, paired being the first two
    results and independent being the last two results.

    Parameters
    ----------
    a : int
        Input the number of trials accomplished.
    b : int
        Input the sample size of each trail.
    c : int
        Set a seed number. Thesis accomplished with seed number 42.

    Returns
    -------
    None.

    '''
    Normal_distribution25, Normal_distribution97 = Testing_Normality(a, b, c)        
    t_statisticPaired_25big, p_valuePaired_25big = stats.ttest_rel(Normal_distribution97, Normal_distribution25, alternative = 'greater')
    t_statisticPaired_25small, p_valuePaired_25small = stats.ttest_rel(Normal_distribution97, Normal_distribution25, alternative = 'less')
    t_statisticInd_25big, p_valueInd_25big = stats.ttest_ind(Normal_distribution97, Normal_distribution25, equal_var=True, alternative = 'greater')
    t_statisticInd_25small, p_valueInd_25small = stats.ttest_ind(Normal_distribution97, Normal_distribution25, equal_var=True, alternative = 'less')

    if p_valuePaired_25big < 0.05:
        print("For a paired TTest, Null hypothesis is rejected. 1997 dataset is consistently greater than the 2025 data set.")
    elif p_valuePaired_25big > 0.05:
        print("For a paired TTest, Null hypothesis cannot be rejected. 1997 dataset is not consistently larger than the 2025 data set.")
    
    if p_valuePaired_25small < 0.05:
        print("For a paired TTest, Null hypothesis is rejected. 1997 dataset is consistently less than the 2025 data set.")
    elif p_valuePaired_25big > 0.05:
        print("For a paired TTest, Null hypothesis cannot be rejected. 1997 dataset is not consistently less than the 2025 data set.")
        
    if p_valueInd_25big < 0.05:
        print("For an Independent TTest, Null hypothesis is rejected. 1997 dataset is consistently greater than the 2025 data set.")
    elif p_valueInd_25big > 0.05:
        print("For an Independent TTest, Null hypothesis cannot be rejected. 1997 dataset is not consistently larger than the 2025 data set.")

    if p_valueInd_25small < 0.05:
        print("For an Independent TTest, Null hypothesis is rejected. 1997 dataset is consistently less than the 2025 data set.")
    elif p_valueInd_25big > 0.05:
        print("For an Independent TTest, Null hypothesis cannot be rejected. 1997 dataset is not consistently less than the 2025 data set.")


def Correlation_Lag(a, b, d, e, f):
    '''
    Function that performs a lag correlation on the 1997 and 2025 δ18O
    data, determining if there is an offset in the data as a result of storage.
    The 2025 dataset is the dataset that is moved up or down core.
    Outputs a figure plotted offset array against correlation coefficent
    values with zero offset marked by a line.

    Parameters
    ----------
    a : int
        Input the starting interval
    b : int
        Input the end interval
    d : float
        Input left end member of interval of depth offset. Ex: -0.2 is moving
        the 2025 dataset up core 0.2 m.
    e : float
        Input right end member of interval of depth offset. Ex: 0.2 is moving
        the 2025 dataset down core 0.2 m.
    f : float
        Input the interval of in which correlation tests occure. Ex: 0.01
        performs a correlation test every 0.01 m from the left end member to
        the right end memeber

    Returns
    -------
    1 array
    correlations : float
        An array of correlation coefficent values with each correlation test
        performed within the interval assigned in f.

    '''
    def Correlation(a, b, c):
        DEPTH25int = []
        DEPTH97int = []

        for i in depth25:
            if a <= i <= b:
                DEPTH25int.append(i)
        for i in depth97:
            if a <= i <= b:
                DEPTH97int.append(i)
        start25 = np.where(depth25 == DEPTH25int[0])[0][0]
        end25 = np.where(depth25 == DEPTH25int[-1])[0][0]
        start97 = np.where(depth97 == DEPTH97int[0])[0][0]
        end97 = np.where(depth97 == DEPTH97int[-1])[0][0]
        O1825int = O1825[start25:(end25+1)]
        O1897int = O1897[start97:(end97+1)]

        Cm1 = []
        for i in range(len(DEPTH25int)):
            Cm1.append(DEPTH25int[i]+c)

        def Average(b, c):
            Interval_sum = 0
            start = DEPTH97int[b]
            end = DEPTH97int[c]
            length = []

            for i in Cm1:
                if start <= i <= end:
                    Interval_sum += i
                    length.append(i)

            d = np.where(Cm1 == length[0])[0][0]
            e = np.where(Cm1 == length[-1])[0][0] + 1
            Interval_O18 = O1825int[d:e]
            Interval_sum25 = 0

            for i in Interval_O18:
                Interval_sum25 += i
            Interval_length25 = len(Interval_O18)
            Average_int25 = Interval_sum25 / Interval_length25
            return (Average_int25)

        Average25int = []

        for a in range(len(DEPTH97int)-1):
            start = a
            end = a + 1
            at = Average(start, end)
            Average25int.append(at)

        Corr_Coef, p_value = spearmanr(O1897int[:-1], Average25int)

        return (Corr_Coef)

    Cm_of_shift = np.arange(d, e, f)
    correlations = []
    for i in Cm_of_shift:
        correlations.append(Correlation(a, b, i))
    plt.figure()
    plt.plot(Cm_of_shift, correlations)
    plt.title("Correlation lag of δ18O, 205 dataset moved")
    plt.xlabel("Meters offset")
    plt.ylabel("Correlation Coefficent")
    plt.axvline(x=0, color='red', linestyle='--')
    plt.show()
    return (correlations)


def Testing_Restablization(a):
    '''
    Function that applies the lag calculated with the Correlation Lag
    function to an interval, and then preforms the Testing Normality and TTest
    functions, to test if this would fix the apparent offset within the data.
    Outputs 3 statements:
        1. The normality of the 1997 and 2025 datasets after applying lag.
        Shapiro test (>0.05 is normal).
        2. Comparison test results of unaltered dataset.
        parametric or nonparametric tests (<0.05 is significant difference).
        3. The shapiro test p-value results for the new CLT distribution
        to prove normality of new distribution after applying lag.
        Shapiro test (>0.05 is normal).
        4. A printed statement of the t-test results of CLT distribution.
        (<0.05 is significant difference)

    Outputs 4 figures:
        1. A figure of the resulting comparison test after testing normality,
        Contains the histogram distribution of the 1997 data, the 2025 data
        and both dataset overlaid.
        2. A figure of the two produced CLT histograms side by side for
        visual comparison of dataset.
        3. A histogram of the 1997 and 2025 CLT dataset histograms
        overlaid.
        4. A graph of the 1997 and 2025 CLT datasets as two
        side by side boxplots. The default is False.

    Parameters
    ----------
    a : float
        Input offset of lag calculated from the correlation lag function

    Returns
    -------
    None.

    '''
    def Lag2025(a, b, c):
        DEPTH25int = []
        DEPTH97int = []

        for i in depth25:
            if a <= i <= b:
                DEPTH25int.append(i)
        for i in depth97:
            if a <= i <= b:
                DEPTH97int.append(i)
        start25 = np.where(depth25 == DEPTH25int[0])[0][0]
        end25 = np.where(depth25 == DEPTH25int[-1])[0][0]
        start97 = np.where(depth97 == DEPTH97int[0])[0][0]
        end97 = np.where(depth97 == DEPTH97int[-1])[0][0]
        O1825int = O1825[start25:(end25+1)]
        O1897int = O1897[start97:(end97+1)]

        Cm1 = []
        for i in range(len(DEPTH25int)):
            Cm1.append(DEPTH25int[i]+c)

        def Average(b, c):
            Interval_sum = 0
            start = DEPTH97int[b]
            end = DEPTH97int[c]
            length = []

            for i in Cm1:
                if start <= i <= end:
                    Interval_sum += i
                    length.append(i)
            d = np.where(Cm1 == length[0])[0][0]
            e = np.where(Cm1 == length[-1])[0][0] + 1
            Interval_O18 = O1825int[d:e]
            Interval_sum25 = 0
            for i in Interval_O18:
                Interval_sum25 += i
            Interval_length25 = len(Interval_O18)
            Average_int25 = Interval_sum25 / Interval_length25
            return (Average_int25)

        Average25int = []

        for a in range(len(DEPTH97int)-1):
            start = a
            end = a + 1
            at = Average(start, end)
            Average25int.append(at)
        return (DEPTH97int, O1897int, Average25int)

    def Testing_Normality(b, c, d, p_value=False, Histograms=False):
        DEPTH97int, O1897int, Average25int = Lag2025(1370, 1890, a)
        O1897int = O1897int[:-1]
        random_index = np.arange(0, len(O1897int))
        Normal_distribution97 = []
        Normal_distribution25 = []
        random.seed(d)
        for i in range(b):
            random_sample = random.choices(random_index, k=c)
            O1897intSAMP = []
            Average25intSAMP = []
            for x in random_sample:
                O1897intSAMP.append(O1897int[x])
                Average25intSAMP.append(Average25int[x])
            Mean97 = np.mean(O1897intSAMP)
            Mean25 = np.mean(Average25intSAMP)
            Normal_distribution25.append(Mean25)
            Normal_distribution97.append(Mean97)

        statistic_normal25, p_value_normal25 = shapiro(Normal_distribution25)
        statistic_normal97, p_value_normal97 = shapiro(Normal_distribution97)

        if p_value is False:
            print()
        elif p_value is True:
            print("1997 p-value is", p_value_normal97, "and 2025 p-value is",
                  p_value_normal25)

        if Histograms is False:
            print()
        elif Histograms is True:
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            plt.subplots_adjust(wspace=0.5, hspace=1)
            fig.suptitle("Testing for Normality in Data")
            axs[0].hist(Normal_distribution97, bins=50, color='orange')
            axs[0].set_ylabel("Frequency")
            axs[0].set_xlabel("δ18O Mean Distribution 1997 dataset")
            axs[0].set_title("Mean δ1O8 Data Distribution 1997 dataset")
            axs[1].hist(Normal_distribution25, bins=50, color='blue')
            axs[1].set_title("Mean δ18O Data Distribution 2025 dataset")
            axs[1].set_ylabel("Frequency")
            axs[1].set_xlabel("Mean δ18O Data Distribution 2025 Dataset")

        return (Normal_distribution25, Normal_distribution97)

    def TTest(a, b, c, Histogram=False, Boxplot=False):
        Normal_distribution25, Normal_distribution97 = Testing_Normality(a, b, c)
        t_statisticPaired, p_valuePaired = stats.ttest_rel(Normal_distribution97,
                                                           Normal_distribution25)
        t_statisticInd, p_valueiInd = stats.ttest_ind(Normal_distribution97,
                                                      Normal_distribution25,
                                                      equal_var=True)

        Caption = ("Paired p-value is", p_valuePaired,
                   "and independent p-value is", p_valueiInd)
        print(Caption)
        if Histogram is False:
            return None
        elif Histogram is True:
            plt.figure()
            plt.hist(Normal_distribution97, bins=50, color='orange')
            plt.hist(Normal_distribution25, bins=50, color='blue', alpha=0.5)
            plt.ylabel("Frequency")
            plt.xlabel("δ18O‰")
            plt.title("Overlay of 1997 and 2025 Mean Distributions")
            plt.show()

        if Boxplot is False:
            return None
        elif Boxplot is True:
            plt.figure()
            plt.boxplot([Normal_distribution97, Normal_distribution25],
                        labels=['1997', '2025'])
            plt.ylabel("δ18O‰")
            plt.title("1997 and 2025 Mean Distributions Boxplots")
            plt.show()

    def Comparison_Tests(a, b, c, Figure=False):
        DEPTH97int, O1897int, Average25int = Lag2025(a, b, c)

        statistic_normal97, p_value_normal97 = shapiro(O1897int)
        statistic_normal25, p_value_normal25 = shapiro(Average25int)

        print("1997 δ18O p-value is:", p_value_normal97,
              "and 2025 p-value is:", p_value_normal25)

        if p_value_normal25 and p_value_normal97 >= 0.05:
            t_statisticPairedT, p_valuePairedT = stats.ttest_rel(O1897int[:-1],
                                                                 Average25int)
            t_statisticIndT, p_valueiIndT = stats.ttest_ind(O1897int[:-1],
                                                            Average25int,
                                                            equal_var=True)
            Caption = ("Parametric test used, paired p-value is",
                       p_valuePairedT, "and independent p-value is",
                       p_valueiIndT)
            print(Caption)
        elif p_value_normal25 or p_value_normal97 < 0.05:
            t_statisticpairW, p_valuepairW = wilcoxon(O1897int[:-1],
                                                      Average25int)
            t_statisticindM, p_valueindM = mannwhitneyu(O1897int[:-1],
                                                        Average25int)
            Caption = ("Nonparametric test used, Paired p-value is",
                       p_valuepairW, "and independent p-value is",
                       p_valueindM)
            print(Caption)

        if Figure is False:
            print()
        elif Figure is True:
            fig, axs = plt.subplots(2, 2, figsize=(5, 5))
            plt.subplots_adjust(wspace=0.5, hspace=1)
            plt.title("Distribiution of δ18O Values within 1997 and 2025 measurement campaigns")
            axs[0, 0].hist(O1897int, bins=50, color="orange")
            axs[0, 0].set_xlabel("δ18O‰ from 1997 data set")
            axs[0, 0].set_ylabel("Frequency")
            axs[1, 0].hist(Average25int, bins=50, color="blue")
            axs[1, 0].set_xlabel("δ18O‰ from 2025 data set")
            axs[1, 0].set_ylabel("Frequency")
            axs[0, 1].hist(O1897int, bins=50, color="orange")
            axs[0, 1].hist(Average25int, bins=50, color="blue", alpha=0.5)
            axs[0, 1].set_xlabel("δ18O‰ from 2025 and 1997 data set")
            axs[0, 1].set_ylabel("Frequency")

    DEPTH97int, O1897int, Average25int = Lag2025(1370, 1890, a)
    Comparison_Tests(1370, 1890, a, True)
    Testing_Normality(2000, 200, 24, True, True)
    TTest(2000, 200, 24, True, True)
