import pandas as pd
import numpy as np
import configparser
import math
import sys
import random
import gzip

from argparse import ArgumentParser
from columnInfo import get_columns

parser = ArgumentParser(description="%prog name")
parser.add_argument("-i", "--input", dest="input", type=str,
                    default=None,
                    help="File location of variant data in tsv format. If not specified use stdin")
parser.add_argument("-o", "--output", dest="output", type=str,
                    default=None,
                    help="Location were the generated file is stored. If not specified use stdout")
parser.add_argument("-c", "--config", dest="config", type=str,
                    default='./impute_GRCh37_v1.6.cfg',
                    help="Config file that specifies used tracks")
parser.add_argument("-b", "--cat2bool", dest="cat2bool", action='store_true',
                    help="Specify whether categories are split into multiple boolean classifier")
parser.add_argument("--noheader", dest="noheader", default=False,
                    action="store_true",
                    help="Do not print header line")

args = parser.parse_args()

# define input and output sources
if args.input is None:
    stdin = sys.stdin
else:
    stdin = open(args.input, 'r')

if args.output is None:
    stdout = sys.stdout
else:
    stdout = open(args.output, 'w')

ImpFile = "/h/xutong/project_RV/21_nonSigGWAS/10_regBaseM/ScoreMetrix_ALL.tsv"
ImpData = pd.read_csv(ImpFile, sep="\t", header=0, index_col=0)

CA = ["DNase-seq", "ATAC-seq"]
HM = ["H2AK5ac", "H2AK9ac", "H2BK120ac", "H2BK12ac", "H2BK15ac", "H2BK20ac", 
      "H2BK5ac", "H3F3A", "H3K14ac", "H3K18ac", "H3K23ac", "H3K23me2", "H3K27ac", 
      "H3K27me3", "H3K36me3", "H3K4ac", "H3K4me1", "H3K4me2", "H3K4me3", "H3K56ac", 
      "H3K79me1", "H3K79me2", "H3K9ac", "H3K9me1", "H3K9me2","H3K9me3", "H3T11ph", 
      "H4K12ac", "H4K20me1", "H4K5ac", "H4K8ac", "H4K91ac"]
TF = ["CTCF", "POLR2A", "RAD21", "SMC3", "EP300", "H2AFZ"]

RM = ["DNase", "H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3", "H3K9me3", 
      "H3K27me3", "H3K9ac", "H3K4me2", "H2AFZ", "H3K79me2", "H4K20me1"]

CA_h = list(map(lambda x: "EpiMap_" + x, CA))
HM_h = list(map(lambda x: "EpiMap_" + x, HM))
TF_h = list(map(lambda x: "EpiMap_" + x, TF))
RM_h = list(map(lambda x: "RoadMap_" + x, RM))

# MissZero (0 or mean)
def score_change(score_list, MissValue = 0):
    for i in range(len(score_list)):
        try:
            score_list[i] = float(score_list[i])
        except:
            score_list[i] = MissValue

    return(score_list)


def RMark_change(mark_list):
    for i in range(len(mark_list)):
        try:
            mark_list[i] = mark_list[i].split("-")[1]
        except:
            mark_list[i] = "."
            
    return(mark_list)


def ColValue_split(Value_list):
    for i in range(len(Value_list)):
        try:
            Value_list[i] = Value_list[i].split(",")[0]
            if Value_list[i] == ".":
                Value_list[i] = "NA"
        except:
            Value_list[i] = "NA"
            
    return(Value_list)


def regBase_Pro(regLine_in, colnms_in, ImpData):
    
    for ColNum in range(len(colnms_in)):
        colnm = colnms_in[ColNum]
        MissValue = ImpData.loc[colnm][0]

        regLine_in[ColNum] = max(score_change(regLine_in[ColNum].split(","), MissValue))
    
    return(regLine_in)


def EpiData_Pro(EpiLine_in, CA, HM, TF, RM):
    
    HumanCAge_score = score_change(EpiLine_in[0].split(","))
    HumanCAge_mark = EpiLine_in[1].split(",")
    HumanHMge_score = score_change(EpiLine_in[2].split(","))
    HumanHMge_mark = EpiLine_in[3].split(",")
    HumanTFge_score = score_change(EpiLine_in[4].split(","))
    HumanTFge_mark = EpiLine_in[5].split(",")

    RoadMap_score = score_change(EpiLine_in[6].split(","))
    RoadMap_cellMark = RMark_change(EpiLine_in[7].split(","))

    # CA
    CA_num = [0] * len(CA)
    for m in range(len(CA)):
        for i in range(len(HumanCAge_mark)):
            if HumanCAge_mark[i] == CA[m]:
                CA_num[m] = max(CA_num[m], HumanCAge_score[i])
            else:
                continue

    # HM
    HM_num = [0] * len(HM)
    for m in range(len(HM)):
        for i in range(len(HumanHMge_mark)):
            if HumanHMge_mark[i] == HM[m]:
                HM_num[m] = max(HM_num[m], HumanHMge_score[i])
            else:
                continue

    # TF
    TF_num = [0] * len(TF)
    for m in range(len(TF)):
        for i in range(len(HumanTFge_mark)):
            if HumanTFge_mark[i] == TF[m]:
                TF_num[m] = max(TF_num[m], HumanTFge_score[i])
            else:
                continue

    # RM            
    RM_num = [0] * len(RM)
    for m in range(len(RM)):
        for i in range(len(RoadMap_cellMark)):
            if RoadMap_cellMark[i] == RM[m]:
                RM_num[m] = max(RM_num[m], RoadMap_score[i])
            else:
                continue

    CA_num = list(map(lambda x: math.log10(x+1), CA_num))
    HM_num = list(map(lambda x: math.log10(x+1), HM_num))
    TF_num = list(map(lambda x: math.log10(x+1), TF_num))
    RM_num = list(map(lambda x: math.log10(x+1), RM_num))

    EpiLine_out = CA_num + HM_num + TF_num + RM_num
    
    return(EpiLine_out)


def CADD_Pro(CADDL_in, colnms_in):
    
    columnNames = list(map(str.lower, colnms_in))

    # associate fields with column names
    CADDL_in = ColValue_split(CADDL_in)
    fieldsDict = dict(zip(columnNames, CADDL_in))
    outFields = []
    indicatorFields = []

    for trackName, status in config['Tracks']:
        if 'colname' in trackData[trackName].keys():
            trackName = trackData[trackName]['colname']
        track = trackData[trackName]

        if status == 'Ignore':
            continue
        if status != 'True':
            continue

        if track['type'] == 'combined':
            if args.cat2bool:
                i = trackData[track['base']]['id']
                baseArray = np.array(outFields[i:i+len(trackData[track['base']]['categories'])])
            else:
                baseValue = outFields[trackData[track['base']]['id']]
                baseArray = (np.array(trackData[track['base']]['categories']) == baseValue).astype(int)
            values = []
            for child in track['child']:
                values.extend(baseArray * outFields[trackData[child]['id']])
            outFields.extend([value for value in values])
        else:  
            try:
                if 'derive' in track.keys():
                    value = track['derive'](fieldsDict)
                else:
                    value = fieldsDict[trackName]
                if track['type'] in [float, int]:
                    value = track['type'](value)
                if 'transformation' in track.keys(): # transform is slightly redundant to derive
                    value = track['transformation'](value)
                if track['type'] is list:
                    assert(value in track['categories'])        
                if 'indicator' in track.keys():
                    indicatorFields.append('0')
            except:
                value = track['na_value']
                if 'indicator' in track.keys():
                    indicatorFields.append('1')

            if args.cat2bool and track['type'] is list:
                values = (np.array(track['categories']) == value).astype(int)
                outFields.extend(values)
            else:
                outFields.append(value)

    # minimize zeros and stringify
    outFields = ['0' if f == 0 else str(f) for f in outFields]
    outFields.extend(indicatorFields)
    
    return(outFields)


### TRACK PREPARATION ###
newColumnNames, _, trackData, config = \
    get_columns(args.config, None, args.cat2bool, False)

### DATA PROCESSING ###
for line in stdin:
    line = line.strip("\n")
    
    # detect header line
    if line.startswith('#') or line.startswith('CHROM'):
        columnNames = line.strip('#').split('\t')
        CADD_colnm = columnNames[0:114]
        reg_colnm = columnNames[125:164]
        newColumnNames = newColumnNames + reg_colnm + CA_h + HM_h + TF_h + RM_h
        if not args.noheader:
            stdout.write('\t'.join(newColumnNames) + '\n')

        continue
        
    # associate fields with column names
    line = line.split('\t')
    
    CADDLine = line[0:114]
    EpiLine = line[114:125]
    regLine = line[125:164]
    
    Epi_outl = EpiData_Pro(EpiLine, CA, HM, TF, RM)
    reg_outl = regBase_Pro(regLine, reg_colnm, ImpData)
    CADD_outl = CADD_Pro(CADDLine, CADD_colnm)

    out_line = CADD_outl + Epi_outl + reg_outl
    out_line = [str(x) for x in out_line]
    stdout.write('\t'.join(out_line) + '\n')
