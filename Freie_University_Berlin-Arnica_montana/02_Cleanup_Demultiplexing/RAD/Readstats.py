#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 25 16:20:54 2025

@author: katja
"""
import pandas as pd
#import numpy as np
import re as re

with open("/Users/katja/Desktop/ArnicaSNP/03_Preprocessed/Testseq_FlexbarSheetWithCounts.tsv") as f:
    dat = f.read().split("\n")

#print(dat[1:20])

heads=[]
for i in range(len(dat)):
    if re.match(".*individualID.*$",dat[i]):
        #print(i)
        heads.append(i)
heads.append(len(dat)+11-3)

runstats = []
readstats = []
textout = str()

for i in range(len(heads)-1):
    name = dat[heads[i]-7]
    if name[0] == "#":
        total = dat[heads[i]-5].split(":")[1].split("(")[0].strip()
        processed = dat[heads[i]-3].split(":")[1].split("(")[0].strip()
        propercent = dat[heads[i]-3].split("(")[2].split("%")[0].strip()
        
        textout += name + "\t" + total + "\t" + processed + "\t" + propercent + "\n"
        runstats.append([name, int(total), int(processed), float(propercent)])
        
        start = heads[i]+2
        stop = heads[i+1]-11
        for j in range(start, stop):
            adapters = dat[j].split("|")[1].strip()
            if adapters[0] == "X":
                samples = dat[j].split("|")[2].strip()
                counts = dat[j].split("|")[3].strip()
                percent = dat[j].split("|")[4].split("%")[0].strip()
            else:
                adapters = dat[j+2].split("|")[1].strip()
                samples = dat[j+2].split("|")[2].strip()
                counts = dat[j].split("|")[3].strip()
                percent = dat[j].split("|")[4].split("%")[0].strip()
            
            textout += adapters + "\t" + samples + "\t" + counts + "\t" + percent + "\n"
            readstats.append([name, adapters, samples, int(counts), float(percent)])
        
runstats = pd.DataFrame(runstats, columns=["Library", "Total Reads", "Processed Reads", "Percent Processed"])

readstats = pd.DataFrame(readstats, columns=["Library", "Adapter Combination", "Sample", "Reads", "Percent Reads"])

#samples = pd.unique(readstats[1])

samplestats = readstats.groupby("Sample")["Reads"].sum()

samplestats.sort_values(ascending=False)

        
        
        