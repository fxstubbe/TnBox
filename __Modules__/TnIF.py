import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from __Modules__.Library import *


##### -------- ---------------- ---------------- ---------------- -------- A TnIF/R100 Algorithm  -------- ---------------- ---------------- ---------------- -------- #####


def vectorizedCoverageMean(x):
    x = x + 1
    x = x.astype(float)
    x = np.log10(x)
    return np.mean(x)


def TnIF(gff, tnif , trim=10):
    means = []
    for index, row in gff.iterrows():
        correction_factor = (row.end-row.start)/trim
        start = round(row.start + correction_factor)
        end = round(row.end - correction_factor)
        data = tnif[tnif.Position.between(start, end)]
        data = data[data.Chromosome == row.seqname]
        if data.empty:
            means.append(0)
            continue
    mean = vectorizedCoverageMean(data.Coverage.to_numpy())
    means.append(mean)
    return means


def vectorizedCoverageMean(x):
    x = x + 1
    x = x.astype(float)
    x = np.log10(x)
    return np.mean(x)


def compile_table(df, gff_p, output_p, graph_p, method, trim, Rwindow, Slide):
    
    '''
    Function that computes a TnIF. It is wrapped in launch compile tr
    sposon in tkin2
    
    method can take either 1 (TnIF) or 2 (R100) as argument
    
    '''

    gff = pd.read_csv(gff_p)
    
    #Loop over the libraries      
    for my_library in df:

    #Read the file
        my_read = f"/Users/stubbf02/Desktop/Py_project/Toolbox/data/{my_library}.txt"
        print(my_library)
        tnif = pd.read_csv(my_read, names = ['Chromosome', 'Position', 'Coverage'] ,delimiter = "\t")
        
        #Get the metric of interest
        
        metric = [] 
        

        # R slide metrics
        #Rwindow = 100 #Still needs to implement the Entry box
        #Slide = 5 #Still needs to implement the Entry box
        #trim = 10 #Should I

        
        for index, row in gff.iterrows():
            correction_factor = (row.end-row.start)/trim
            start = round(row.start + correction_factor)
            end = round(row.end - correction_factor)
            data = tnif[tnif.Position.between(start, end)]
            data = data[data.Chromosome == row.seqname]
            
            #Compute the TnIF (method = 1)
            if method == 1 :
                if data.empty:
                    metric.append(0)
                    continue
                mean = vectorizedCoverageMean(data.Coverage.to_numpy())
                metric.append(mean)
                
            #Compute the R100 (method = 2)    
            else : 
            	#What about defining the window size and the slide ?
                rlist = data.Coverage.to_numpy()
                my_insertion_list = [list(data.Coverage[i : i + Rwindow]).count(0) for i in range(0, len(data.Coverage), Slide)]
                my_Rwindow = my_insertion_list.count(Rwindow)
                metric.append(my_Rwindow)
                

        #append library column    
        gff[my_library] = metric
        
       



    
    gff.to_csv(output_p)
    print("Non Indexed Summary table is done !")
