#
# read MB csv file and generate report
#
import csv
import os
import sys
import time
from datetime import datetime
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tabulate import tabulate
from scipy import stats
#
# 12/30/25 remove CD90/CD117 from CD11,CD14,CD16
#
def classify_doublet(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    if total < 10:
        return False
    r_cd8 = cd8_value / total
    r_cd11b = cd11b_value / total
    r_cd14 = cd14_value / total
    r_cd16 = cd16_value / total
    r_cd20 = cd20_value / total
    r_cd90 = cd90_value / total
    r_cd117 = cd117_value / total
    r_cd34 = cd34_value / total
    r_cd4 = cd4_value / total

    if r_cd4 > 0.05 and ((r_cd14 > 0.05) or
                         (r_cd16 > 0.05) or
                         (r_cd20 > 0.05) or
                         (r_cd90 > 0.05) or
                         (r_cd117 > 0.05) or
                         (r_cd34 > 0.05)):
        return True

    if r_cd8 > 0.05 and ((r_cd11b > 0.5) or
                         (r_cd14 > 0.05) or
                         (r_cd16 > 0.05) or
                         (r_cd20 > 0.05) or
                         (r_cd90 > 0.05) or
                         (r_cd117 > 0.05) or
                         (r_cd34 > 0.05)):
        return True

    if r_cd11b > 0.05 and ((r_cd20 > 0.05) or
                         #(r_cd90 > 0.05) or
                         (r_cd117 > 0.05) or
                         (r_cd34 > 0.05)):
        return True

    if r_cd14 > 0.05 and ((r_cd16 > 0.5) or
                          (r_cd20 > 0.05) or
                         #(r_cd90 > 0.05) or
                         #(r_cd117 > 0.05) or
                         (r_cd34 > 0.05)):
        return True

    if r_cd16 > 0.05 and ((r_cd20 > 0.05) or
                         #(r_cd90 > 0.05) or
                         #(r_cd117 > 0.05) or
                         (r_cd34 > 0.05)):
        return True

    if r_cd20 > 0.05 and ((r_cd90 > 0.05) or
                         (r_cd117 > 0.05) or
                         (r_cd34 > 0.05)):
        return True

    return False





def classify_CD4_1(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd4_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r = cd4_value / total
    if r >= 0.95:
        return "CD4"
    else:
        return ""

def classify_CD4_2(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd4_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd4_value / total
    r2 = cd8_value / total
    max_noise = max(cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value) / total

    if r1 >= 0.90 and r2 < 0.10 and max_noise < 0.05:
        return "CD4"
    else:
        return ""

def classify_CD4_3(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd4_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd4_value / total
    max_noise = max(cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value) / total

    if r1 >= 0.80 and max_noise < 0.05:
        return "CD4"
    else:
        return ""
#
#----------------------------------------------------------------------------------------------------------------
#
def classify_CD8_1(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd8_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r = cd8_value / total
    if r >= 0.95:
        #print("1    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "CD8"
    else:
        return ""

def classify_CD8_2(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd8_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd8_value / total
    r2 = cd4_value / total
    max_noise = max(cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value) / total
    if r1 >= 0.90 and r2 < 0.10 and max_noise < 0.05:
        #print("2    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "CD8"
    else:
        return ""

def classify_CD8_3(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd8_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd8_value / total
    max_noise = max(cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value) / total

    if r1 >= 0.80 and max_noise < 0.05:
        #print("3    ",r1,max_noise)
        #print("3    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "CD8"
    else:
        return ""
#
#----------------------------------------------------------------------------------------------------------------
#  B Cells
#
def classify_B_1(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd20_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd20_value / total
    if r1 >= 0.95:
        #print("1    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "B"
    else:
        return ""

def classify_B_2(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd20_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd20_value / total
    max_noise = max(cd8_value, cd11b_value, cd14_value, cd16_value, cd90_value, cd117_value, cd34_value, cd4_value) / total

    if r1 >= 0.80 and max_noise < 0.05:
        #print("2    ",r1,max_noise)
        #print("2    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "B"
    else:
        return ""
#
#----------------------------------------------------------------------------------------------------------------
# NK Cells
def classify_NK_1(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd16_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd16_value / total
    if r1 >= 0.95:
        #print("1    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "NK"
    else:
        return ""

def classify_NK_2(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd16_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd16_value / total
    r2 = cd11b_value / total

    if r1  + r2 >= 0.90:
        #print("2    ",r1,max_noise)
        #print("2    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "NK"
    else:
        return ""

def classify_NK_3(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd16_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd16_value / total
    r2 = cd11b_value / total
    r3 = cd90_value / total

    if r1 + r2 + r3 >= 0.90:
        #print("3    ",r1,max_noise)
        #print("3    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "NK"
    else:
        return ""
#
#----------------------------------------------------------------------------------------------------------------
# Granulocytes
def classify_Granulocytes_1(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd90_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd90_value / total
    if r1 >= 0.95:
        #print("1    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Granulocytes"
    else:
        return ""

def classify_Granulocytes_2(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd11b_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd11b_value / total

    if r1 >= 0.90:
        #print("2    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Granulocytes"
    else:
        return ""

def classify_Granulocytes_3(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd11b_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd11b_value / total
    r2 = cd90_value / total

    if r1 + r2 >= 0.90:
        #print("3    ",r1,r2)
        #print("3    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Granulocytes"
    else:
        return ""
def classify_Granulocytes_4(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd11b_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd90_value / total
    r2 = cd11b_value / total
    r3 = cd117_value / total

    if r1 + r2 + r3 >= 0.90:
        #print("4    ",r1,max_noise)
        #print("4    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Granulocytes"
    else:
        return ""
def classify_Granulocytes_5(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd11b_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd90_value / total
    r2 = cd11b_value / total
    r3 = cd117_value / total
    r4 = cd4_value / total

    if r1 + r2 + r3 + r4 >= 0.90:
        #print("5    ",r1,max_noise)
        #print("5    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Granulocytes"
    else:
        return ""
#
#----------------------------------------------------------------------------------------------------------------
# Monocytes
def classify_Monocytes_1(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd14_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd14_value / total
    if r1 >= 0.95:
        #print("1    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Monocytes"
    else:
        return ""
def classify_Monocytes_2(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd14_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd14_value / total
    r2 = cd11b_value / total

    if r1 + r2 >= 0.90:
        #print("2    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Monocytes"
    else:
        return ""
def classify_Monocytes_3(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd14_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd14_value / total
    r2 = cd90_value / total

    if r1 + r2 >= 0.90:
        #print("3    ",r1,r2)
        #print("3    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Monocytes"
    else:
        return ""

def classify_Monocytes_4(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd14_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd14_value / total
    r2 = cd90_value / total
    r3 = cd117_value / total

    if r1 + r2 + r3 >= 0.90:
        #print("3    ",r1,r2)
        #print("3    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Monocytes"
    else:
        return ""
def classify_Monocytes_5(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd14_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value

    r1 = cd14_value / total
    r2 = cd90_value / total
    r3 = cd117_value / total
    r4 = cd11b_value / total

    if r1 + r2 + r3 + r4 >= 0.90:
        #print("3    ",r1,r2)
        #print("3    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Monocytes"
    else:
        return ""
#
#----------------------------------------------------------------------------------------------------------------
# HSC
def classify_HSC_1(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd34_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd90_value / total
    r2 = cd34_value / total
    if (r1 + r2 >= 0.9) and (r1 > 0.1):
        #print("HSC1    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "HSC"
    else:
        return ""

def classify_HSC_2(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd90_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd90_value / total
    r2 = cd117_value / total
    if (r1 + r2 >= 0.9) and (r1 > 0.1) and (r2 > 0.1):
        #print("HSC2    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "HSC"
    else:
        return ""
def classify_HSC_3(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if cd34_value < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd90_value / total
    r2 = cd34_value / total
    r3 = cd117_value / total
    if (r1 + r2 + r3 >= 0.9) and (r1 > 0.1) and (r2 > 0.1):
        #print("HSC3   ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "HSC"
    else:
        return ""
#
#----------------------------------------------------------------------------------------------------------------
# Progrenitor Cells

def classify_Progenitor(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
    if (cd117_value + cd34_value) < 5:
        return ""
    total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
    r1 = cd117_value / total
    r2 = cd34_value / total
    if (r1 + r2 >= 0.95):
        #print("Progenitor    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
        return "Progenitor"
    else:
        return ""

# def classify_Late_Progenitor(cd8_value, cd11b_value, cd14_value, cd16_value, cd20_value, cd90_value, cd117_value, cd34_value, cd4_value):
#     if cd34_value < 5:
#         return ""
#     total = cd8_value + cd11b_value + cd14_value + cd16_value + cd20_value + cd90_value + cd117_value + cd34_value + cd4_value
#     r1 = cd34_value / total
#     if r1 >= 0.95:
#         #print("Late Progenitor    ",cd4_value,cd8_value,cd11b_value,cd14_value,cd16_value,cd20_value,cd90_value,cd117_value,cd34_value)
#         return "Late Progenitor"
#     else:
#         return ""


#
# classify cell types based on marker values
#
def ClassifyCellTypes(cd8_values, cd11b_values, cd14_values, cd16_values,
                      cd20_values, cd90_values, cd117_values, cd34_values, cd4_values):
    pCellType_values = []
    unknown_values = []
    doublte_values = []
    for cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4 in zip(
            cd8_values, cd11b_values, cd14_values, cd16_values,
            cd20_values, cd90_values, cd117_values, cd34_values, cd4_values):
        if classify_doublet(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4):
            pCellType_values.append("Doublet")
            doublte_values.append([cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4])
            continue
        # first check if CD4 cell
        if classify_CD4_1(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "CD4":
            pCellType_values.append("CD4")
            continue
        if classify_CD4_2(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "CD4":
            pCellType_values.append("CD4")
            continue
        if classify_CD4_3(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "CD4":
            pCellType_values.append("CD4")
            continue
        if classify_CD8_1(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "CD8":
            pCellType_values.append("CD8")
            continue
        if classify_CD8_2(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "CD8":
            pCellType_values.append("CD8")
            continue
        if classify_CD8_3(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "CD8":
            pCellType_values.append("CD8")
            continue
        if classify_B_1(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "B":
            pCellType_values.append("B")
            continue
        if classify_B_2(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "B":
            pCellType_values.append("B")
            continue
        if classify_NK_1(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "NK":
            pCellType_values.append("NK")
            continue
        if classify_NK_2(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "NK":
            pCellType_values.append("NK")
            continue
        if classify_NK_3(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "NK":
            pCellType_values.append("NK")
            continue
        if classify_Granulocytes_1(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Granulocytes":
            pCellType_values.append("Granulocytes")
            continue
        if classify_Granulocytes_2(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Granulocytes":
            pCellType_values.append("Granulocytes")
            continue
        if classify_Granulocytes_3(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Granulocytes":
            pCellType_values.append("Granulocytes")
            continue
        if classify_Granulocytes_4(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Granulocytes":
            pCellType_values.append("Granulocytes")
            continue
        if classify_Granulocytes_5(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Granulocytes":
            pCellType_values.append("Granulocytes")
            continue
        if classify_Monocytes_1(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Monocytes":
            pCellType_values.append("Monocytes")
            continue
        if classify_Monocytes_2(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Monocytes":
            pCellType_values.append("Monocytes")
            continue
        if classify_Monocytes_3(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Monocytes":
            pCellType_values.append("Monocytes")
            continue
        if classify_Monocytes_4(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Monocytes":
            pCellType_values.append("Monocytes")
            continue
        if classify_Monocytes_5(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Monocytes":
            pCellType_values.append("Monocytes")
            continue
        if classify_HSC_1(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "HSC":
            pCellType_values.append("HSC")
            continue
        if classify_HSC_2(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "HSC":
            pCellType_values.append("HSC")
            continue
        if classify_HSC_3(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "HSC":
            pCellType_values.append("HSC")
            continue
        if classify_Progenitor(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Progenitor":
            pCellType_values.append("Progenitor")
            continue
        # if classify_Late_Progenitor(cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4) == "Late Progenitor":
        #     pCellType_values.append("Late Progenitor")
        #     continue
        # if no classification, assign "Unknown"
        s = sum([cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4])
        if s >= 2:
            unknown_values.append((cd8, cd11b, cd14, cd16, cd20, cd90, cd117, cd34, cd4))
        pCellType_values.append("Unknown")

    return (pCellType_values, unknown_values,doublte_values)
#
# classify edit based on number of edited -v- unedited sequences
#
#  - is defined as knockout...meaning sequence was edited
#  + means WT, no editing
#
#
def classifyEdit(nEdit,nTotal):
    # if not enouth data, return "Unknown"
    if nTotal < 5:
        return "Unknown"
    # claculate edit fraction
    r = nEdit / nTotal
    # classify based on values
    if r < 0.05:
        return "+/+"
    elif r >= 0.05 and r < 0.95:
        return "+/-"
    elif r >= 0.95:
        return "-/-"
    else:
        return "Unknown"
#
#
#   build a table for data summary
#
def buildEditTable1(label1,column_data_1, pCellType_values):
    #
    # create data structure to hold counts
    #
    output_counts = defaultdict(lambda: defaultdict(int))
    #
    # initialize output_counts for all combinations of input gene
    #
    for l1 in ['-/-','+/-','+/+','Unknown']:
        editClassd = f"{l1}"
        for pCellType in ['B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Myeloid Progenitor','Late Progenitor','Unknown','Doublet']:
            output_counts[editClassd][pCellType] = 0
    #
    # count occurrences of each combination of gene_edit and pCellType
    #
    for index in range(len(column_data_1)):
        c1_edit = column_data_1[index]
        editClassd = f"{c1_edit}"
        pCellType = pCellType_values[index]
        output_counts[editClassd][pCellType] += 1
    #
    #
    #
    table = []
    pCellType_list = ['B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Progenitor','Unknown','Doublet']
    header = [label1,'B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Progenitor','Unknown','Doublet']
    table.append(header)
    #
    # create table rows
    #

    for kClass,outputDict in output_counts.items():
        row = kClass.split(',')
        for pCellType in pCellType_list:
            row.append(outputDict[pCellType])
        table.append(row)

    return table
#
# 2 targets
#
def buildEditTable2(label1,label2,column_data_1, column_data_2, pCellType_values):
    #
    # create data structure to hold counts
    #
    output_counts = defaultdict(lambda: defaultdict(int))
    #
    # initialize output_counts for all combinations of CD33_edit and CD33P_edit and pCellType
    #
    for l1 in ['-/-','+/-','+/+','Unknown']:
        for l2 in ['-/-','+/-','+/+','Unknown']:
            editClassd = f"{l1},{l2}"
            for pCellType in ['B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Myeloid Progenitor','Late Progenitor','Unknown','Doublet']:
                output_counts[editClassd][pCellType] = 0
    #
    # count occurrences of each combination of CD33_edit, CD33P_edit and pCellType
    #
    for index in range(len(column_data_1)):
        c1_edit = column_data_1[index]
        c2_edit = column_data_2[index]
        editClassd = f"{c1_edit},{c2_edit}"
        pCellType = pCellType_values[index]
        output_counts[editClassd][pCellType] += 1
    #
    #
    #
    table = []
    pCellType_list = ['B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Progenitor','Unknown','Doublet']
    header = [label1,label2,'B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Progenitor','Unknown','Doublet']
    table.append(header)
    #
    # create table rows
    #

    for kClass,outputDict in output_counts.items():
        row = kClass.split(',')
        for pCellType in pCellType_list:
            row.append(outputDict[pCellType])
        table.append(row)

    return table
#
#
#   build a table for three data types
#
#
def buildEditTable3(label1,label2,label3,column_data_1, column_data_2, column_data_3,pCellType_values):
    #
    # create data structure to hold counts
    #
    output_counts_3 = defaultdict(lambda: defaultdict(int))
    #
    # initialize output_counts for all combinations of l1,l2,l3 and pCellType
    #
    for l1 in ['-/-','+/-','+/+','Unknown']:
        for l2 in ['-/-','+/-','+/+','Unknown']:
          for l3 in ['-/-','+/-','+/+','Unknown']:
            editClassd = f"{l1},{l2},{l3}"
            for pCellType in ['B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Myeloid Progenitor','Late Progenitor','Unknown','Doublet']:
                output_counts_3[editClassd][pCellType] = 0
    #
    # count occurrences of each combination of CD33_edit, CD33P_edit and pCellType
    #
    for index in range(len(column_data_1)):
        c1_edit = column_data_1[index]
        c2_edit = column_data_2[index]
        c3_edit = column_data_3[index]
        editClassd = f"{c1_edit},{c2_edit},{c3_edit}"
        pCellType = pCellType_values[index]
        output_counts_3[editClassd][pCellType] += 1
    #
    #
    #
    table = []
    pCellType_list = ['B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Progenitor','Unknown','Doublet']
    header = [label1,label2,label3,'B','CD4','CD8','NK','Monocytes','Granulocytes','HSC','Progenitor','Unknown','Doublet']
    table.append(header)
    #
    # create table rows
    #
    for kClass,outputDict in output_counts_3.items():
        row = kClass.split(',')
        for pCellType in pCellType_list:
            row.append(outputDict[pCellType])
        table.append(row)
    return table
#-----------------------------------------------------------------------------------------------------------------
#
#  Main program
#
#
#
#
#
#-----------------------------------------------------------------------------------------------------------------
#
if len(sys.argv) < 2:
    print("Usage: python report.py <csv_file>")
    sys.exit(1)

csv_file = sys.argv[1]
if not os.path.isfile(csv_file):
    print(f"File {csv_file} does not exist.")
    sys.exit(1)

# read tsv file
data = pd.read_csv(csv_file, sep='\t',error_bad_lines=False)
init_Lenght = len(data)
# print first 5 rows
#print(data.head())
# print columns
#print(data.columns)
# print number of rows
print(f"Number of rows: {len(data)}")
# check if required columns exist
required_columns = ['CellBarcode','nCap','CD33_e', 'CD33_t',
                    'CD33P_e','CD33P_t',
                    'HBG1_e','HBG1_t',
                    'HBG2_e','HBG2_t',
                    'CD8', 'CD11b', 'CD14', 'CD16', 'CD20', 'CD90', 'CD117', 'CD34', 'CD4']
if not all(col in data.columns for col in required_columns):
    print(f"Missing columns in input file: {set(required_columns) - set(data.columns)}")
    sys.exit(1)
#
# for each row, get CD33_e, CD33_ne, CD33P_e, CD33P_ne ,HBG1_e, HBG1_ne, HBG2_e, HBG2_ne
# filter out rows where each pair is less than 5

#
data = data[data['nCap'] >= 1000]
#
print(f"length of data after nCap 1000 filter: {len(data)}")

data = data[(data['CD33_t'] >= 5) |
            (data['CD33P_t'] >= 5) |
            (data['HBG1_t'] >= 5) |
            (data['HBG2_t'] >= 5)]
# print number of rows after filtering
print(f"Number of rows after DNA edit sequence filtering: {len(data)}")
#
# column sum of CD33_total, CD33P_total, HBG1_total, HBG2_total
#
s33_e = sum(data['CD33_e'])
s33_t = sum(data['CD33_t'])
s33p_e = sum(data['CD33P_e'])
s33p_t = sum(data['CD33P_t'])
sHBG1_e = sum(data['HBG1_e'])
sHBG1_t = sum(data['HBG1_t'])
sHBG2_e = sum(data['HBG2_e'])
sHBG2_t = sum(data['HBG2_t'])
print(f"Total CD33 sequences:  {s33_t:7d}, Edited: {s33_e:7d} % = {s33_e/(s33_t)*100:5.2f}%")
print(f"Total CD33P sequences: {s33p_t:7d}, Edited: {s33p_e:7d} % = {s33p_e/(s33p_t)*100:5.2f}%")
print(f"Total HBG1 sequences:  {sHBG1_t:7d}, Edited: {sHBG1_e:7d} % = {sHBG1_e/(sHBG1_t)*100:5.2f}%")
print(f"Total HBG2 sequences:  {sHBG2_t:7d}, Edited: {sHBG2_e:7d} % = {sHBG2_e/(sHBG2_t)*100:5.2f}%")
#
# get data columans as lists
#
cd8_values = data['CD8'].tolist()
cd11b_values = data['CD11b'].tolist()
cd14_values = data['CD14'].tolist()
cd16_values = data['CD16'].tolist()
cd20_values = data['CD20'].tolist()
cd90_values = data['CD90'].tolist()
cd117_values = data['CD117'].tolist()
cd34_values = data['CD34'].tolist()
cd4_values = data['CD4'].tolist()

pCellType_values,unknown_values,doublet_values = ClassifyCellTypes(cd8_values, cd11b_values, cd14_values, cd16_values,
                                    cd20_values, cd90_values, cd117_values, cd34_values, cd4_values)

#
# summarize editing
#
cd33_e_values = data['CD33_e'].tolist()
cd33_t_values = data['CD33_t'].tolist()
cd33p_e_values = data['CD33P_e'].tolist()
cd33p_t_values = data['CD33P_t'].tolist()
hbg1_e_values = data['HBG1_e'].tolist()
hbg1_t_values = data['HBG1_t'].tolist()
hbg2_e_values = data['HBG2_e'].tolist()
hbg2_t_values = data['HBG2_t'].tolist()
#
# + = WT
# - = edited
#
# for each CD33 AND CD33P value, define CD33_Edit as follows:
# CD33P < 0.05 -> "+/+
# CD33P >= 0.05 and < 0.95 -> "+/-"
# CD33P >= 0.95 -> -/-
#
CD33_edit = []
for ce,cne in zip(cd33_e_values, cd33_t_values):
    e = classifyEdit(ce, cne)
    CD33_edit.append(e)
#
#
#
CD33P_edit = []
for ce,cne in zip(cd33p_e_values, cd33p_t_values):
    e = classifyEdit(ce, cne)
    CD33P_edit.append(e)
#
#
#
HBG1_edit = []
for ce,cne in zip(hbg1_e_values, hbg1_t_values):
    e = classifyEdit(ce, cne)
    HBG1_edit.append(e)

HBG2_edit = []
for ce,cne in zip(hbg2_e_values, hbg2_t_values):
    e = classifyEdit(ce, cne)
    HBG2_edit.append(e)
#
# combined tables
#
cd33_single_table = buildEditTable1("CD33",CD33_edit,pCellType_values)
#
# print fancy table
#
print()
print(tabulate(cd33_single_table, headers="firstrow", tablefmt="grid"))
#
#
#
cd33p_single_table = buildEditTable1("CD33P",CD33P_edit,pCellType_values)
#
# print fancy table
#
print()
print(tabulate(cd33p_single_table, headers="firstrow", tablefmt="grid"))
print("\n\n")
#
#
#
hbg1_single_table = buildEditTable1("HBG1",HBG1_edit,pCellType_values)
#
# print fancy table
#
print('---hbg1---')
print(tabulate(hbg1_single_table, headers="firstrow", tablefmt="grid"))
print("\n\n")
#
#
#
hbg2_single_table = buildEditTable1("HBG2",HBG2_edit,pCellType_values)
#
# print fancy table
#
print('---hbg2---')
print(tabulate(hbg2_single_table, headers="firstrow", tablefmt="grid"))
print("\n\n")
#
# cd33/cd33p table
#
cd33_table = buildEditTable2("CD33","CD33P",CD33_edit, CD33P_edit, pCellType_values)
#
# print fancy table
#
print(tabulate(cd33_table, headers="firstrow", tablefmt="grid"))
print("\n\n")
#
# hbg1/hbg2 table
#
hbg_table = buildEditTable2("HBG1","HBG2",HBG1_edit, HBG2_edit, pCellType_values)
#
# print fancy table
#
print(tabulate(hbg_table, headers="firstrow", tablefmt="grid"))
print("\n\n")
#
# 3-way
#
cd33_hbg1_hbg2 = buildEditTable3("CD33","HBG1","HBG2",CD33_edit,HBG1_edit,HBG2_edit,pCellType_values)
print(tabulate(cd33_hbg1_hbg2, headers="firstrow", tablefmt="grid"))

cd33p_hbg1_hbg2 = buildEditTable3("CD33P","HBG1","HBG2",CD33P_edit,HBG1_edit,HBG2_edit,pCellType_values)
print(tabulate(cd33p_hbg1_hbg2 , headers="firstrow", tablefmt="grid"))
#
# generate output csv file with a row for each output value and columns for each pCellType value
#
output_csv_file = os.path.splitext(csv_file)[0] + "_report.csv"
with open(output_csv_file, 'w', newline='') as f:
    # header
    f.write(f"Report for file: {csv_file}\n")
    f.write(f"Total rows in input file: {init_Lenght}\n")
    f.write(f"Number of rows after DNA edit sequence filtering: {len(data)}\n")
    f.write(f"Total CD33 sequences:  {s33_t:7d}, Edited: {s33_e:7d} % = {s33_e/(s33_t)*100:5.2f}\n")
    f.write(f"Total CD33P sequences: {s33p_t:7d}, Edited: {s33p_e:7d} % = {s33p_e/(s33p_t)*100:5.2f}\n")
    f.write(f"Total HBG1 sequences:  {sHBG1_t:7d}, Edited: {sHBG1_e:7d} % = {sHBG1_e/(sHBG1_t)*100:5.2f}\n")
    f.write(f"Total HBG2 sequences:  {sHBG2_t:7d}, Edited: {sHBG2_e:7d} % = {sHBG2_e/(sHBG2_t)*100:5.2f}\n\n")
    #
    # summarize CD33
    #
    s_cd33dict = defaultdict(int)
    totalKnwon = 0
    for e in CD33_edit:
        if e != "Unknown":
            s_cd33dict[e] += 1
            totalKnwon += 1
    f.write("CD33 Summary\n")
    for k in ['-/-','+/-','+/+']:
        f.write(f"  {k:7s} : {s_cd33dict[k]:7d} cells, {s_cd33dict[k]/totalKnwon*100:5.2f}%\n")
    f.write("\n")
    #
    # summarize CD33P
    #
    s_cd33pdict = defaultdict(int)
    totalKnown = 0
    for e in CD33P_edit:
        if e != "Unknown":
            s_cd33pdict[e] += 1
            totalKnown += 1
    f.write("CD33P Summary\n")
    for k in ['-/-','+/-','+/+']:
        f.write(f"  {k:7s} : {s_cd33pdict[k]:7d} cells, {s_cd33pdict[k]/totalKnown*100:5.2f}%\n")
    f.write("\n")
    #
    # summarize HBG1
    #
    s_hbg1dict = defaultdict(int)
    totalKnown = 0
    for e in HBG1_edit:
        if e != "Unknown":
            s_hbg1dict[e] += 1
            totalKnown += 1
    f.write("HBG1 Summary\n")
    for k in ['-/-','+/-','+/+']:
        f.write(f"  {k:7s} : {s_hbg1dict[k]:7d} cells, {s_hbg1dict[k]/totalKnown*100:5.2f}%\n")
    f.write("\n")
    #
    # summarize HBG2
    #
    s_hbg2dict = defaultdict(int)
    totalKnown = 0
    for e in HBG2_edit:
        if e != "Unknown":
            s_hbg2dict[e] += 1
            totalKnown += 1
    f.write("HBG2 Summary\n")
    for k in ['-/-','+/-','+/+']:
        f.write(f"  {k:7s} : {s_hbg2dict[k]:7d} cells, {s_hbg2dict[k]/totalKnown*100:5.2f}%\n")
    f.write("\n")
    #
    # write CD33 table to csv file
    #
    writer = csv.writer(f)
    writer.writerows(cd33_single_table)
    f.write("\n\n\n")
    #
    # write CD33P table to csv file
    #
    writer = csv.writer(f)
    writer.writerows(cd33p_single_table)
    f.write("\n\n\n")
    #
    # write HBG1 table to csv file
    #
    writer = csv.writer(f)
    writer.writerows(hbg1_single_table)
    f.write("\n\n\n")
    #
    # write HBG2 table to csv file
    #
    writer = csv.writer(f)
    writer.writerows(hbg2_single_table)
    f.write("\n\n\n")
    #
    # write CD33/CD33P table to csv file
    #
    writer = csv.writer(f)
    writer.writerows(cd33_table)
    f.write("\n\n\n")
    #
    # write HBG1/HBG2 table to csv file
    #
    writer = csv.writer(f)
    writer.writerows(hbg_table)
    f.write("\n\n\n")
    #
    # write cd33_bhg1_hbg2
    #
    writer = csv.writer(f)
    writer.writerows(cd33_hbg1_hbg2)
    f.write("\n\n\n")
    #
    # write cd33p_bhg1_hbg2
    #
    writer = csv.writer(f)
    writer.writerows(cd33p_hbg1_hbg2)
print(f"Report saved to {output_csv_file}")




#
# save unknown values
#
unknown_csv_file = os.path.splitext(csv_file)[0] + "_unknown.csv"
with open(unknown_csv_file, 'w', newline='') as f:
    # header
    writer = csv.writer(f)
    writer.writerow(['CD8', 'CD11b', 'CD14', 'CD16', 'CD20', 'CD90', 'CD117', 'CD34', 'CD4'])
    writer.writerows(unknown_values)
print(f"Unknown values saved to {unknown_csv_file}")
#
# save doublet values
#
doublet_csv_file = os.path.splitext(csv_file)[0] + "_doublet.csv"
with open(doublet_csv_file, 'w', newline='') as f:
    # header
    writer = csv.writer(f)
    writer.writerow(['CD8', 'CD11b', 'CD14', 'CD16', 'CD20', 'CD90', 'CD117', 'CD34', 'CD4'])
    writer.writerows(doublet_values)
print(f"Doublet values saved to {doublet_csv_file}")
#
# save barcodes
#
bc_file = os.path.splitext(csv_file)[0] + "_bc.csv"
print(type(data))
with open(bc_file,'w') as f:
  for bc in data['CellBarcode']:
    f.write(f"{bc}\n")

#
# look at editing
#
# for index in range(0,len(data)):
# #     print('CD33  `',CD33_edit[index], cd33_e_values[index], cd33_ne_values[index])
# #     print('CD33P `',CD33P_edit[index], cd33p_e_values[index], cd33p_ne_values[index])
# #     #      cd33p_e_values[index], cd33p_ne_values[index],
# #     #      HBG1_edit[index], hbg1_e_values[index], hbg1_ne_values[index],
# #     #      HBG2_edit[index], hbg2_e_values[index], hbg2_ne_values[index],
# #     #      pCellType_values[index])
# #     #HBG1_edit[index], HBG2_edit[index], pCellType_values[index])
#         if CD33P_edit[index] == "+/+":
#             print("--------------------------------------------------")
#             print(f"{' '*54}  CD4\tCD8\tCD11b\tCD14\tCD16\tCD20\tCD90\tCD117\tCD34")
#             print(f"Row {index} is {CD33_edit[index]} {CD33P_edit[index]} Cell Type {pCellType_values[index]:20s}\t{cd4_values[index]}\t{cd8_values[index]}\t{cd11b_values[index]}\t{cd14_values[index]}\t{cd16_values[index]}\t{cd20_values[index]}\t{cd90_values[index]}\t{cd117_values[index]}\t{cd34_values[index]}")
#             print(f"  CD33  : Edited {cd33_e_values[index]}, Total {cd33_t_values[index]}")
#             print(f"  CD33P : Edited {cd33p_e_values[index]}, Total {cd33p_t_values[index]}")
