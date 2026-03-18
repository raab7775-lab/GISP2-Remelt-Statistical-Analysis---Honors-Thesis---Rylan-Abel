# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 21:12:08 2025

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

GISP04 = pd.read_csv(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\Old_Analysis\2004_GISP2_Holocene_deuterium.txt")
GISP04 = pd.read_table(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\Old_Analysis\2004_GISP2_Holocene_deuterium.txt", sep='\t', skiprows=[0])
depth25 = np.loadtxt(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\New_Melt\2025_depth_model.txt")
D225 = np.loadtxt(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\New_Melt\2025_data_dD_SIL_cfa.txt")
ModGISP04 = GISP04[GISP04["deuterium"] != 999999]
GISP04pt2 = pd.read_csv(r"C:\Users\omgit\OneDrive\Desktop\GISP2 Statistical Analysis - Rylan Abel\Original Analysis\All data used (edited)\Old_Analysis\2004_GISP2_Below_1600_deuterium.csv")
ModGISP04pt2 = GISP04pt2[GISP04pt2["delta D (per mil)"] != 999999]

# Interval Given is between 1372 and 1598
# Overlap between dataset is 1372 to 1598

top04 = np.array(ModGISP04["Top"])
bottom04 = np.array(ModGISP04["Bottom"])
D204 = np.array(ModGISP04["deuterium"])
top04pt2 = np.array(ModGISP04pt2["top depth (m)"])
bottom04pt2 = np.array(ModGISP04pt2["bottom depth (m)"])
D204pt2 = np.array(ModGISP04pt2["delta D (per mil)"])
# Takes pandas column and turns it into a numpy array that can be in function

# Function that removes the anomalous -300 per mil values from the Holocene
# dataset


def take_above_300(D204):
    indicies = np.where(D204 < -300)[0].tolist()
    mask = np.ones(len(D204), dtype=bool)
    mask[indicies] = False
    D204filt = D204[mask]
    top04filt = top04[mask]
    bottom04filt = bottom04[mask]
    return (D204filt, top04filt, bottom04filt)


D204, top04, bottom04 = take_above_300(D204)


def Data_Unedited_figure():
    '''
    Function that outputs a figure of both datasets as it is, without
    any editing. No input, output is a single figure.'

    Input
    None.

    Returns
    None.

    '''
    plt.figure(figsize=(10, 5))
    plt.plot(depth25, D225,
             label="δD 2025 data", color="red")
    plt.plot(top04, D204,
             label="δD 2004 data sporadic interval", color="#ffd042")
    plt.plot(top04pt2, D204pt2,
             label="δD 2004 data 0.5 interval", color="orange")
    plt.xlabel("Depth (m)")
    plt.ylabel("δD per mL")
    plt.title("Comparison Between 2025 and 2004 GISP2 Data")
    plt.legend()
    plt.show()


def IntervalHolo(a, b, Figure=False):
    '''
    Function that performs downsampling of 2025 data to the 2004 resolution but
    only for the 2004 Holocene data

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
    DEPTH04Correct : float
        Output is the corrected depth array.
    D204int : float
        The corresponding δD values from the 2004 dataset.
    Average25int : float
        The average of the 2025 δD values corresponding to the 2004
        depth resolution

    '''
    DEPTH25int = []
    DEPTH04top = []

    for i in depth25:
        if a <= i <= b:
            DEPTH25int.append(i)
    for i in top04:
        if a <= i <= b:
            DEPTH04top.append(i)
    start25 = np.where(depth25 == DEPTH25int[0])[0][0]
    end25 = np.where(depth25 == DEPTH25int[-1])[0][0] + 1
    start04 = np.where(top04 == DEPTH04top[0])[0][0]
    end04 = np.where(top04 == DEPTH04top[-1])[0][0] + 1
    D225int = D225[start25:end25]
    D204int = D204[start04:end04]
    DEPTH04bottom = bottom04[start04:end04]

    DEPTH25int = np.array(DEPTH25int)

    def Average(b, c):
        Interval_sum = 0
        start = DEPTH04top[b]
        end = DEPTH04bottom[c]
        length = []
        for i in DEPTH25int:
            if start <= i <= end:
                Interval_sum += i
                length.append(i)
        if not length:
            return None
        d = np.where(DEPTH25int == length[0])[0][0]
        e = np.where(DEPTH25int == length[-1])[0][0] + 1
        Interval_D2 = D225int[d:e]
        Interval_sum2 = 0
        for i in Interval_D2:
            Interval_sum2 += i
        Interval_length2 = len(Interval_D2)
        Average_int2 = Interval_sum2 / Interval_length2
        return (Average_int2)

    Average25int = []

    for a in range(len(DEPTH04top)):
        start = a
        end = a
        at = Average(start, end)
        if at is None:
            continue
        Average25int.append(at)

    half = np.array(ModGISP04["length"])/2

    DEPTH04correct = []
    for i in range(len(DEPTH04top)):
        DEPTH04correct.append(DEPTH04top[i] + half[i])

    if Figure is False:
        print()
    elif Figure is True:
        plt.plot(DEPTH04correct, D204int, label="δD 2004 data",
                 color='orange')
        plt.plot(DEPTH04correct, Average25int,
                 label="Average δD from 2025 Data", color="red")
        plt.xlabel("Depth (m)")
        plt.ylabel("δD per mL")
        plt.title("Comparison of δD Between 2025 and 2004 GISP2 Data")
        plt.legend()
        plt.show()

    return (DEPTH04correct, D204int, Average25int)


def Intervalpt2(a, b, Figure=False):
    '''
    Function that performs downsampling of 2025 data to the 2004 resolution
    but only for the 2004 dataset that is below 1600 meters at a resolution
    of 0.5 m.

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
    DEPTH04Correct : float
        Output is the corrected depth array.
    D204int : float
        The corresponding δD values from the 2004 dataset.
    Average25int : float
        The average of the 2025 δD values corresponding to the 2004
        depth resolution

    '''
    DEPTH25int = []
    DEPTH04top = []

    for i in depth25:
        if a <= i <= b:
            DEPTH25int.append(i)
    for i in top04pt2:
        if a <= i <= b:
            DEPTH04top.append(i)

    start25 = np.where(depth25 == DEPTH25int[0])[0][0]
    end25 = np.where(depth25 == DEPTH25int[-1])[0][0] + 1
    start04 = np.where(top04pt2 == DEPTH04top[0])[0][0]
    end04 = np.where(top04pt2 == DEPTH04top[-1])[0][0] + 1
    D225int = D225[start25:end25]
    D204int = D204pt2[start04:end04]
    DEPTH04bottom = bottom04pt2[start04:end04]

    DEPTH25int = np.array(DEPTH25int)

    def Average(b, c):
        Interval_sum = 0
        start = DEPTH04top[b]
        end = DEPTH04bottom[c]
        length = []
        for i in DEPTH25int:
            if start <= i <= end:
                Interval_sum += i
                length.append(i)
        if not length:
            return None
        d = np.where(DEPTH25int == length[0])[0][0]
        e = np.where(DEPTH25int == length[-1])[0][0] + 1
        Interval_D2 = D225int[d:e]
        Interval_sum2 = 0
        for i in Interval_D2:
            Interval_sum2 += i
        Interval_length2 = len(Interval_D2)
        Average_int2 = Interval_sum2 / Interval_length2
        return (Average_int2)

    Average25int = []

    for a in range(len(DEPTH04top)):
        start = a
        end = a
        at = Average(start, end)
        if at is None:
            continue
        Average25int.append(at)

    DEPTH04Correct = []
    for i in DEPTH04top:
        DEPTH04Correct.append(i + 0.25)

    if Figure is False:
        print()
    elif Figure is True:
        plt.plot(DEPTH04Correct, D204int, label="dD 2004 data",
                 color='orange')
        plt.plot(DEPTH04Correct, Average25int,
                 label="Average δD 2025 Data", color="red")
        plt.xlabel("Depth (m)")
        plt.ylabel("δD‰")
        plt.title("Comparison of δD Between 2025 and 2004 GISP2 Data")
        plt.legend()
        plt.show()

    return (DEPTH04Correct, D204int, Average25int)


def Interval_complete(a, b, Figure=False):
    '''
    Function that performs downsampling of 2025 data to the 2004 resolution.
    Full dataset, combines the funtions IntervalHolo and Intervalpart2

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
    DEPTH04Correct : float
        Output is the corrected depth array.
    D204int : float
        The corresponding δD values from the 2004 dataset.
    Average25int : float
        The average of the 2025 δD values corresponding to the 2004
        depth resolution

    '''
    if b < 1600:
        DEPTH04int, D204int, Average25int = IntervalHolo(a, b)
    elif a > 1600:
        DEPTH04int, D204int, Average25int = Intervalpt2(a, b)
    elif b >= 1600 and a <= 1600:
        DEPTH04intpt1, D204intpt1, Average25intpt1 = IntervalHolo(a, 1600)
        DEPTH04intpt2, D204intpt2, Average25intpt2 = Intervalpt2(1600, b)
        DEPTH04int = np.concatenate((DEPTH04intpt1[:-1], DEPTH04intpt2[:-1]))
        D204int = np.concatenate((D204intpt1[:-1], D204intpt2[:-1]))
        Average25int = np.concatenate((Average25intpt1[:-1],
                                       Average25intpt2[:-1]))

    title = "Same-Resolution Comparison of δD Between 2025 and 2004 GISP2 Data"

    if Figure is False:
        print()
    elif Figure is True:
        if len(DEPTH04int) == len(Average25int):

            plt.figure()
            plt.plot(DEPTH04int, D204int, label="dD 2004 data",
                     color='orange')
            plt.plot(DEPTH04int, Average25int,
                     label="Average δD 2025 Data", color="red")
            plt.xlabel("Depth (m)")
            plt.ylabel("δD‰")
            plt.title(title)
            plt.legend()
            plt.show()
        else:
            plt.figure()
            plt.plot(DEPTH04int, D204int, label="dD 2004 data",
                     color='orange')
            plt.plot(DEPTH04int[:-1], Average25int,
                     label="Average δD 2025 Data", color="red")
            plt.xlabel("Depth (m)")
            plt.ylabel("δD‰")
            plt.title(title)
            plt.legend()
            plt.show()

    return (DEPTH04int, D204int, Average25int)


def Moving_Mean(a, b, c, Figure=False):
    '''
    Moving mean function applied to a same resolution interval of the 2004
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
        Input true for a figure with a graph of both the 2004 and 2025 datasets
        smoothed by the established smoothing interval.
        The default is False.

    Returns
    -------
    2 Arrays
    Mean04 : float
        An array of the 2004 data smoothed by the smoothing interval
    Mean25 : float
        An array of the 2025 data smoothed by the smoothing interval

    '''
    if c == 1:
        New = 1
    elif c > 1:
        New = int(c/2)

    if b < 1600:
        DEPTH04int, D204int, Average25int = IntervalHolo((a - c), (b + c))
    elif a > 1600:
        DEPTH04int, D204int, Average25int = Intervalpt2((a - c), (b + c))
    elif b >= 1600 and a <= 1600:
        DEPTH04intpt1, D204intpt1, Average25int1 = IntervalHolo((a - c), 1600)
        DEPTH04intpt2, D204intpt2, Average25int2 = Intervalpt2(1600, (b + c))
        DEPTH04int = np.concatenate((DEPTH04intpt1, DEPTH04intpt2))
        D204int = np.concatenate((D204intpt1, D204intpt2))
        Average25int = np.concatenate((Average25int1, Average25int2))

    def Interval_04(a, New):

        start = a
        Depthmeanint = []
        for i in DEPTH04int:
            if start - New <= i <= start + New:
                Depthmeanint.append(i)
        e = np.where(DEPTH04int == Depthmeanint[0])[0][0]
        f = np.where(DEPTH04int == Depthmeanint[-1])[0][0] + 1
        intD204 = D204int[e:f]
        MEAN04 = np.mean(intD204)
        return (MEAN04)

    def Interval_av(a, New):
        start = a
        Depthmeanint = []
        for i in DEPTH04int:
            if start - New <= i <= start + New:
                Depthmeanint.append(i)
        e = np.where(DEPTH04int == Depthmeanint[0])[0][0]
        f = np.where(DEPTH04int == Depthmeanint[-1])[0][0] + 1
        intD2Average = Average25int[e:f]
        MEAN25 = np.mean(intD2Average)
        return (MEAN25)

    MEAN04 = []
    MEAN25 = []
    for i in range(len(DEPTH04int)):
        start = int(DEPTH04int[i])
        at = Interval_04(start, New)
        ab = Interval_av(start, New)
        MEAN04.append(at)
        MEAN25.append(ab)

    title = 'Running Mean', c,
    'Interval, Comparison of δD Between 2025 and 2004 GISP2 Data'

    if Figure is False:
        print()
    elif Figure is True:
        plt.figure()
        plt.plot(DEPTH04int[c:-c], MEAN04[c:-c],
                 label="Smoothed Mean of δD 2004 data", color='orange')
        plt.plot(DEPTH04int[c:-c], MEAN25[c:-c],
                 label="Smoothed Mean of Average δD from 2025 Data",
                 color='red')
        plt.xlabel("Depth (m)")
        plt.ylabel("Smoothed Mean of δD Data")
        plt.title(title)
        plt.legend()
        plt.show()

    return (MEAN04, MEAN25)


def Difference(a, b, Figure=False, Figure2=False):
    '''
    Fuction that calculates the difference between the 2004 and 2025 data with
    depth.

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
    DEPTH04int, D204int, Average25int = Interval_complete(a, b)

    Difference = Average25int - D204int[:-1]
    DifferenceARR = np.array(Difference)

    Difference_mean = np.mean(DifferenceARR)
    Difference_STD = np.std(DifferenceARR)

    print(Difference_mean, Difference_STD)

    if Figure is False:
        print()
    elif Figure is True:
        plt.figure()
        plt.hist(Difference, bins=50, color="red")
        plt.title("Distribution of Difference in δD in the 2004 and 2025 Data")
        plt.xlabel("Difference in δD‰")
        plt.ylabel("Frequency")
        plt.show()

    if Figure2 is False:
        print()
    elif Figure2 is True:
        plt.figure()
        plt.plot(DEPTH04int[:-1], DifferenceARR, color="red")
        plt.xlabel("Depth (m)")
        plt.ylabel("Difference in δD‰")
        plt.title("Difference between 2004 and 2025 δD datasets with Depth")
        plt.show()

    return (DifferenceARR)


def Comparison_Tests(a, b, Figure=False):
    '''
    Tests for the normality of the 2004 and 2025 data set, then applys a
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
        Contains the histogram distribution of the 2004 data, the 2025 data
        and both dataset overlaid. The default is False.

    Returns
    -------
    None

    '''

    DEPTH04int, D204int, Average25int = Interval_complete(a, b)

    statistic_normal04, p_value_normal04 = shapiro(D204int)
    statistic_normal25, p_value_normal25 = shapiro(Average25int)

    print("2004 dD p-value is:", p_value_normal04,
          "and 2025 p-value is:", p_value_normal25)

    if p_value_normal25 and p_value_normal04 >= 0.05:
        t_statisticPairedT, p_valuePairedT = stats.ttest_rel(D204int,
                                                             Average25int)
        t_statisticIndT, p_valueiIndT = stats.ttest_ind(D204int,
                                                        Average25int,
                                                        equal_var=True)
        Caption = ("Parametric test used, paired p-value is", p_valuePairedT,
                   "and independent p-value is", p_valueiIndT)
        print(Caption)
    elif p_value_normal25 or p_value_normal04 < 0.05:
        t_statisticpairW, p_valuepairW = wilcoxon(D204int, Average25int)
        t_statisticindM, p_valueindM = mannwhitneyu(D204int,
                                                    Average25int)
        Caption = ("Nonparametric test used, Paired p-value is", p_valuepairW,
                   "and independent p-value is", p_valueindM)
        print(Caption)

    if Figure is False:
        print()
    elif Figure is True:
        fig, axs = plt.subplots(2, 2, figsize=(5, 5))
        plt.subplots_adjust(wspace=0.5, hspace=1)
        plt.title("Distribiution of δD Values within 2004 and 2025 datasets")
        axs[0, 0].hist(D204int, bins=50, color="orange")
        axs[0, 0].set_xlabel("δD‰ from 2004 data set")
        axs[0, 0].set_ylabel("Frequency")
        axs[1, 0].hist(Average25int, bins=50, color="red")
        axs[1, 0].set_xlabel("δD‰ from 2025 data set")
        axs[1, 0].set_ylabel("Frequency")
        axs[0, 1].hist(D204int, bins=50, color="orange")
        axs[0, 1].hist(Average25int, bins=50, color="darkred", alpha=0.5)
        axs[0, 1].set_xlabel("δD‰ from 2004 and 2025 data set")
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
    Normal_distribution04 : float
        The CLT dataset produced for the 2004 dataset.

    '''

    DEPTH04int, D204int, Average25int = Interval_complete(1370, 1890)

    D204int = D204int
    random_index = np.arange(0, len(D204int))
    Normal_distribution04 = []
    Normal_distribution25 = []
    random.seed(c)
    for i in range(a):
        random_sample = random.choices(random_index, k=b)
        D204intSAMP = []
        Average25intSAMP = []
        for x in random_sample:
            D204intSAMP.append(D204int[x])
            Average25intSAMP.append(Average25int[x])
        Mean04 = np.mean(D204intSAMP)
        Mean25 = np.mean(Average25intSAMP)
        Normal_distribution25.append(Mean25)
        Normal_distribution04.append(Mean04)

    statistic_normal25, p_value_normal25 = shapiro(Normal_distribution25)
    statistic_normal04, p_value_normal04 = shapiro(Normal_distribution04)

    if p_value is False:
        print()
    elif p_value is True:
        print("2004 p-value is", p_value_normal04,
              "and 2025 p-value is", p_value_normal25)

    if Histograms is False:
        print()
    elif Histograms is True:
        caption = ("2004:", round(p_value_normal04, 5), "2025:",
                   round(p_value_normal25, 5))
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        plt.subplots_adjust(wspace=0.5, hspace=1)
        fig.suptitle("Testing for Normality in Data")
        axs[0].hist(Normal_distribution04, bins=50, color='orange')
        axs[0].set_title("Mean δD Data Distribution 2004 dataset")
        axs[0].set_xlabel("Mean δD Distribution 2004 dataset")
        axs[0].set_title("Mean δD Data Distribution 2004 dataset")
        axs[1].hist(Normal_distribution25, bins=50, color='red')
        axs[1].set_title("Mean δD Data Distribution 2025 dataset")
        axs[1].set_ylabel("Frequency")
        axs[1].set_xlabel("Mean δD Data Distribution 2025 dataset")
        plt.figtext(0.43, 0, caption)

    return (Normal_distribution25, Normal_distribution04)


def Wilcoxon(a, b, c, Histogram=False, Boxplot=False):
    '''
    Function that performs a Wilcoxon on the normal distributions formed
    by the Testing Normality function. Tests both a paired and independent
    nonparametric test. Output is a printed statement of the test results.
    (<0.05 is significant difference). Provided because the CLT method has
    had some problems with the 2004 data specifically with appearing normal
    depending on seed chosen. 42 is seed used in analysis, use the t-test
    variation with that seed.

    Parameters
    ----------
    a : int
        Input the number of trials accomplished.
    b : int
        Input the sample size of each trail.
    c : int
        Set a seed number.
    Histogram : bool, optional
        Outputs a histogram of the 2004 and 2025 CLT dataset histograms
        overlaid. The default is False.
    Boxplot : bool, optional
        Outputs a graph of the 2004 and 2025 CLT datasets as two
        side by side boxplots. The default is False.

    Returns
    -------
    None.

    '''
    Normal_distribution25, Normal_distribution04 = Testing_Normality(a, b, c)
    t_statisticPaired, p_valuePaired = wilcoxon(Normal_distribution04,
                                                Normal_distribution25)

    t_statisticInd, p_valueiInd = mannwhitneyu(Normal_distribution04,
                                               Normal_distribution25)

    Caption = ("Paired p-value is", p_valuePaired,
               "and independent p-value is", p_valueiInd)
    print(Caption)

    if Histogram is False:
        return None
    elif Histogram is True:
        plt.figure()
        plt.hist(Normal_distribution04, bins=50, color='orange')
        plt.hist(Normal_distribution25, bins=50, color='red', alpha=0.5)
        plt.ylabel("Frequency")
        plt.xlabel("δD per ml")
        plt.title("Overlay of 2004 and 2025 Mean Distributions")
        plt.show()

    if Boxplot is False:
        return None
    elif Boxplot is True:
        plt.figure()
        plt.boxplot([Normal_distribution04, Normal_distribution25],
                    labels=['2004', '2025'])
        plt.ylabel("δD per ml")
        plt.title("2004 and 2025 Mean Distributions Boxplots")
        plt.show()


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
        Outputs a histogram of the 2004 and 2025 CLT dataset histograms
        overlaid. The default is False.
    Boxplot : bool, optional
        Outputs a graph of the 2004 and 2025 CLT datasets as two
        side by side boxplots. The default is False.

    Returns
    -------
    None.

    '''
    Normal_distribution25, Normal_distribution04 = Testing_Normality(a, b, c)
    t_statisticPaired, p_valuePaired = stats.ttest_rel(Normal_distribution04,
                                                       Normal_distribution25)
    t_statisticInd, p_valueiInd = stats.ttest_ind(Normal_distribution04,
                                                  Normal_distribution25,
                                                  equal_var=True)

    Caption = ("Paired p-value is", p_valuePaired,
               "and independent p-value is", p_valueiInd)
    print(Caption)

    if Histogram is False:
        return None
    elif Histogram is True:
        plt.figure()
        plt.hist(Normal_distribution04, bins=50, color='orange')
        plt.hist(Normal_distribution25, bins=50, color='red', alpha=0.5)
        plt.ylabel("Frequency")
        plt.xlabel("δD per ml")
        plt.title("Overlay of 2004 and 2025 Mean Distributions")
        plt.figtext(0.01, -0.1, Caption)
        plt.show()

    if Boxplot is False:
        return None
    elif Boxplot is True:
        plt.figure()
        plt.boxplot([Normal_distribution04, Normal_distribution25],
                    labels=['2004', '2025'])
        plt.ylabel("δD per ml")
        plt.title("2004 and 2025 Mean Distributions Boxplots")
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
    Normal_distribution25, Normal_distribution04 = Testing_Normality(a, b, c)        
    t_statisticPaired_25big, p_valuePaired_25big = stats.ttest_rel(Normal_distribution04, Normal_distribution25, alternative = 'greater')
    t_statisticPaired_25small, p_valuePaired_25small = stats.ttest_rel(Normal_distribution04, Normal_distribution25, alternative = 'less')
    t_statisticInd_25big, p_valueInd_25big = stats.ttest_ind(Normal_distribution04, Normal_distribution25, equal_var=True, alternative = 'greater')
    t_statisticInd_25small, p_valueInd_25small = stats.ttest_ind(Normal_distribution04, Normal_distribution25, equal_var=True, alternative = 'less')

    if p_valuePaired_25big < 0.05:
        print("For a paired TTest, Null hypothesis is rejected. 2004 dataset is consistently greater than the 2025 data set.")
    elif p_valuePaired_25big > 0.05:
        print("For a paired TTest, Null hypothesis cannot be rejected. 2004 dataset is not consistently larger than the 2025 data set.")
    
    if p_valuePaired_25small < 0.05:
        print("For a paired TTest, Null hypothesis is rejected. 2004 dataset is consistently less than the 2025 data set.")
    elif p_valuePaired_25big > 0.05:
        print("For a paired TTest, Null hypothesis cannot be rejected. 2004 dataset is not consistently less than the 2025 data set.")
        
    if p_valueInd_25big < 0.05:
        print("For an Independent TTest, Null hypothesis is rejected. 2004 dataset is consistently greater than the 2025 data set.")
    elif p_valueInd_25big > 0.05:
        print("For an Independent TTest, Null hypothesis cannot be rejected. 2004 dataset is not consistently larger than the 2025 data set.")

    if p_valueInd_25small < 0.05:
        print("For an Independent TTest, Null hypothesis is rejected. 2004 dataset is consistently less than the 2025 data set.")
    elif p_valueInd_25big > 0.05:
        print("For an Independent TTest, Null hypothesis cannot be rejected. 2004 dataset is not consistently less than the 2025 data set.")


def Correlation_Lag(a, b, d, e, f):
    '''
    Function that performs a lag correlation on the 2004 and 2025 δD
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

        for i in depth25:
            if a <= i <= b:
                DEPTH25int.append(i)
        start25 = np.where(depth25 == DEPTH25int[0])[0][0]
        end25 = np.where(depth25 == DEPTH25int[-1])[0][0] + 1
        D225int = D225[start25:end25]

        if b < 1600:
            DEPTH04top = []
            for i in top04:
                if a <= i <= b:
                    DEPTH04top.append(i)
            start04 = np.where(top04 == DEPTH04top[0])[0][0]
            end04 = np.where(top04 == DEPTH04top[-1])[0][0] + 1
            D204int = D204[start04:end04]
            DEPTH04bottom = bottom04[start04:end04]
        elif a > 1600:
            DEPTH04top = []
            for i in top04pt2:
                if a <= i <= b:
                    DEPTH04top.append(i)
            start04 = np.where(top04pt2 == DEPTH04top[0])[0][0]
            end04 = np.where(top04pt2 == DEPTH04top[-1])[0][0] + 1
            D204int = D204pt2[start04:end04]
            DEPTH04bottom = bottom04pt2[start04:end04]
        elif b >= 1600 and a <= 1600:
            DEPTH04toppt1 = []
            DEPTH04toppt2 = []
            for i in top04:
                if a <= i <= b:
                    DEPTH04toppt1.append(i)
            start04pt1 = np.where(top04 == DEPTH04toppt1[0])[0][0]
            end04pt1 = np.where(top04 == DEPTH04toppt1[-1])[0][0] + 1
            D204intpt1 = D204[start04pt1:end04pt1]
            DEPTH04bottompt1 = bottom04[start04pt1:end04pt1]

            for i in top04pt2:
                if a <= i <= b:
                    DEPTH04toppt2.append(i)
            start04pt2 = np.where(top04pt2 == DEPTH04toppt2[0])[0][0]
            end04pt2 = np.where(top04pt2 == DEPTH04toppt2[-1])[0][0] + 1
            D204intpt2 = D204pt2[start04pt2:end04pt2]
            DEPTH04bottompt2 = bottom04pt2[start04pt2:end04pt2]

            DEPTH04top = np.concatenate((DEPTH04toppt1[:-1], DEPTH04toppt2))
            DEPTH04bottom = np.concatenate((DEPTH04bottompt1[:-1],
                                            DEPTH04bottompt2))
            D204int = np.concatenate((D204intpt1[:-1], D204intpt2))

        Cm1 = []
        for i in range(len(DEPTH25int)):
            Cm1.append(DEPTH25int[i]+c)

        def Average(b, c):
            Interval_sum = 0
            start = DEPTH04top[b]
            end = DEPTH04bottom[c]
            length = []

            for i in Cm1:
                if start <= i <= end:
                    Interval_sum += i
                    length.append(i)

            Cm1Arr = np.array(Cm1)
            d = np.where(Cm1Arr == length[0])[0][0]
            e = np.where(Cm1Arr == length[-1])[0][0] + 1
            Interval_D2 = D225int[d:e]
            Interval_sum25 = 0

            for i in Interval_D2:
                Interval_sum25 += i
            Interval_length25 = len(Interval_D2)
            Average_int25 = Interval_sum25 / Interval_length25
            return (Average_int25)

        Average25int = []

        for a in range(len(DEPTH04top)-1):
            start = a
            end = a
            at = Average(start, end)
            Average25int.append(at)

        Corr_Coef, p_value = spearmanr(D204int[:-1], Average25int)

        return (Corr_Coef)

    Cm_of_shift = np.arange(d, e, f)
    correlations = []
    for i in Cm_of_shift:
        correlations.append(Correlation(a, b, i))
    plt.figure()
    plt.plot(Cm_of_shift, correlations, color="red")
    plt.axvline(x=0, color='blue', linestyle='--')
    plt.show()
    return (correlations)


def Testing_Restablization(a):
    ''' Function that applies the lag calculated with the Correlation Lag
    function to an interval, and then preforms the Testing Normality and TTest
    functions, to test if this would fix the apparent offset within the data.
    Input the starting interval (a), the end interval (b), the amount of
    offset (c), the seed of the mean distribution produced (d), set p_value to
    true to print shapiro test results, set Histograms to true to print a
    figure of the mean distributions for each data, set Histograms2 to True to
    output a histogram of the mean distributions overlapped, and set
    Boxplot to output a boxplot of both normal distributions next to eachother.
    Outputs Several figures and returns the shifted Lag2025 depth array,
    deuterium 2004 data array, and the deuterium 2025 data array.
    '''
    def Lag2025(a, b, c):
        DEPTH25int = []
        for i in depth25:
            if a <= i <= b:
                DEPTH25int.append(i)
        start25 = np.where(depth25 == DEPTH25int[0])[0][0]
        end25 = np.where(depth25 == DEPTH25int[-1])[0][0] + 1
        D225int = D225[start25:end25]

        DEPTH04toppt1 = []
        DEPTH04toppt2 = []
        for i in top04:
            if a <= i <= b:
                DEPTH04toppt1.append(i)
        start04pt1 = np.where(top04 == DEPTH04toppt1[0])[0][0]
        end04pt1 = np.where(top04 == DEPTH04toppt1[-1])[0][0] + 1
        D204intpt1 = D204[start04pt1:end04pt1]
        DEPTH04bottompt1 = bottom04[start04pt1:end04pt1]

        for i in top04pt2:
            if a <= i <= b:
                DEPTH04toppt2.append(i)
        start04pt2 = np.where(top04pt2 == DEPTH04toppt2[0])[0][0]
        end04pt2 = np.where(top04pt2 == DEPTH04toppt2[-1])[0][0] + 1
        D204intpt2 = D204pt2[start04pt2:end04pt2]
        DEPTH04bottompt2 = bottom04pt2[start04pt2:end04pt2]

        DEPTH04top = np.concatenate((DEPTH04toppt1[:-1], DEPTH04toppt2))
        DEPTH04bottom = np.concatenate((DEPTH04bottompt1[:-1],
                                        DEPTH04bottompt2))
        D204int = np.concatenate((D204intpt1[:-1], D204intpt2))

        Cm1 = []
        for i in range(len(DEPTH25int)):
            Cm1.append(DEPTH25int[i]+c)

        def Average(b, c):
            Interval_sum = 0
            start = DEPTH04top[b]
            end = DEPTH04bottom[c]
            length = []

            for i in Cm1:
                if start <= i <= end:
                    Interval_sum += i
                    length.append(i)

            Cm1Arr = np.array(Cm1)
            d = np.where(Cm1Arr == length[0])[0][0]
            e = np.where(Cm1Arr == length[-1])[0][0] + 1
            Interval_D2 = D225int[d:e]
            Interval_sum25 = 0

            for i in Interval_D2:
                Interval_sum25 += i
            Interval_length25 = len(Interval_D2)
            Average_int25 = Interval_sum25 / Interval_length25
            return (Average_int25)

        Average25int = []

        for a in range(len(DEPTH04top)-1):
            start = a
            end = a
            at = Average(start, end)
            Average25int.append(at)

        half = np.array(ModGISP04["length"])/2

        DEPTH04correct = []
        for i in range(len(DEPTH04top)):
            DEPTH04correct.append(DEPTH04top[i] + half[i])

        return (DEPTH04correct, D204int[:-1], Average25int)

    def Comparison_Tests(a, b, c, Figure=False):
        DEPTH04int, D204int, Average25int = Lag2025(a, b, c)

        statistic_normal04, p_value_normal04 = shapiro(D204int)
        statistic_normal25, p_value_normal25 = shapiro(Average25int)

        print("2004 dD p-value is:", p_value_normal04,
              "and 2025 p-value is:", p_value_normal25)

        if p_value_normal25 and p_value_normal04 >= 0.05:
            t_statisticPairedT, p_valuePairedT = stats.ttest_rel(D204int,
                                                                 Average25int)
            t_statisticIndT, p_valueiIndT = stats.ttest_ind(D204int,
                                                            Average25int,
                                                            equal_var=True)
            Caption = ("Parametric test used, paired p-value is",
                       p_valuePairedT, "and independent p-value is",
                       p_valueiIndT)
            print(Caption)
        elif p_value_normal25 or p_value_normal04 < 0.05:
            t_statisticpairW, p_valuepairW = wilcoxon(D204int, Average25int)
            t_statisticindM, p_valueindM = mannwhitneyu(D204int,
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
            plt.title("Distribiution of δD Values within 2004 and 2025 Measurment Campaigns")
            axs[0, 0].hist(D204int, bins=50, color="orange")
            axs[0, 0].set_xlabel("δD‰ from 2004 data set")
            axs[0, 0].set_ylabel("Frequency")
            axs[1, 0].hist(Average25int, bins=50, color="red")
            axs[1, 0].set_xlabel("δD‰ from 2025 data set")
            axs[1, 0].set_ylabel("Frequency")
            axs[0, 1].hist(D204int, bins=50, color="orange")
            axs[0, 1].hist(Average25int, bins=50, color="darkred", alpha=0.5)
            axs[0, 1].set_xlabel("δD‰ from 2004 and 2025 data set")
            axs[0, 1].set_ylabel("Frequency")

    def Testing_Normality_SHIFT(d, b, c, p_value=False, Histograms=False):
        DEPTH04int, D204int, Average25int = Lag2025(1370, 1890, a)
        random_index = np.arange(0, len(D204int))
        Normal_distribution04 = []
        Normal_distribution25 = []
        random.seed(c)
        for i in range(d):
            random_sample = random.choices(random_index, k=b)
            D204intSAMP = []
            Average25intSAMP = []
            for x in random_sample:
                D204intSAMP.append(D204int[x])
                Average25intSAMP.append(Average25int[x])
            Mean04 = np.mean(D204intSAMP)
            Mean25 = np.mean(Average25intSAMP)
            Normal_distribution25.append(Mean25)
            Normal_distribution04.append(Mean04)

        statistic_normal25, p_value_normal25 = shapiro(Normal_distribution25)
        statistic_normal04, p_value_normal04 = shapiro(Normal_distribution04)

        if p_value is False:
            print()
        elif p_value is True:
            print("2004 p-value is", p_value_normal04,
                  "and 2025 p-value is", p_value_normal25)

        if Histograms is False:
            print()
        elif Histograms is True:
            caption = ("2004:", round(p_value_normal04, 5), "2025:",
                       round(p_value_normal25, 5))
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            plt.subplots_adjust(wspace=0.5, hspace=1)
            fig.suptitle("Testing for Normality in Data")
            axs[0].hist(Normal_distribution04, bins=50, color='orange')
            axs[0].set_title("Mean δD Data Distribution 2004 dataset")
            axs[0].set_xlabel("δD Mean Distribution 2004 dataset")
            axs[0].set_title("Mean δD Data Distribution 2004 dataset")
            axs[1].hist(Normal_distribution25, bins=50, color='red')
            axs[1].set_title("Mean δD Data Distribution 2025 dataset")
            axs[1].set_ylabel("Frequency")
            axs[1].set_xlabel("Mean δD Data Distribution 2025 datast")
            plt.figtext(0.43, 0, caption)

        return (Normal_distribution25, Normal_distribution04)

    def TTest_SHIFT(a, b, c, Histogram=False, Boxplot=False):
        Normal_distribution25, Normal_distribution04 = Testing_Normality_SHIFT(a, b, c)
        t_statisticPaired, p_valuePaired = stats.ttest_rel(Normal_distribution04,
                                                           Normal_distribution25)
        t_statisticInd, p_valueiInd = stats.ttest_ind(Normal_distribution04,
                                                      Normal_distribution25,
                                                      equal_var=True)

        Caption = ("Paired p-value is", p_valuePaired,
                   "and independent p-value is", p_valueiInd)
        print(Caption)

        if Histogram is False:
            return None
        elif Histogram is True:
            plt.figure()
            plt.hist(Normal_distribution04, bins=50, color='orange')
            plt.hist(Normal_distribution25, bins=50, color='red', alpha=0.5)
            plt.ylabel("Frequency")
            plt.xlabel("δD per ml")
            plt.title("Overlay of 2004 and 2025 Mean Distributions")
            plt.figtext(0.01, -0.1, Caption)
            plt.show()

        if Boxplot is False:
            return None
        elif Boxplot is True:
            plt.figure()
            plt.boxplot([Normal_distribution04, Normal_distribution25],
                        labels=['2004', '2025'])
            plt.ylabel("δD per ml")
            plt.title("2004 and 2025 Mean Distributions Boxplots")
            plt.show()

    DEPTH04int, D204int, Average25int = Lag2025(1370, 1890, a)
    Comparison_Tests(1370, 1890, a, True)
    Testing_Normality_SHIFT(2000, 200, 42, True, True)
    TTest_SHIFT(2000, 200, 42, True, True)
