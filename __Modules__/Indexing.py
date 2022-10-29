import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from __Modules__.Library import *



def FindCorr(array, index):
    if float(index + 6) > np.size(array):
        return (np.size(array))
    else : 
        return (index+6)

def FindMode(cdt): 

	Mode = None
	binwidth = 0.1

	if(pd.isna(cdt).any()):

		print ("Array is empty, return 0")

		return 0

	else :

		values = plt.hist(cdt, bins = np.arange(min(cdt), max(cdt) + binwidth, binwidth))
		Y_axis, X_axis, patches = values

		

		#Flip the axis
		Y_axis_flip = np.flip(Y_axis)
		X_axis_flip = np.flip(X_axis)

		#Find the mode
		for count, value in enumerate(Y_axis_flip[1:], start = 1):
			if Y_axis_flip[count] > Y_axis_flip[count - 1] and Y_axis_flip[count] > Y_axis_flip[count + 1]:
				if Y_axis_flip[count] > np.mean(Y_axis_flip[count : FindCorr(Y_axis_flip, count)]): 
					true_mode = round(X_axis_flip[count+1],2)
					Mode = true_mode
					break
				

	if Mode is None : 
		return 0

	else : 
		return Mode

def Indexing(Reference, Condition):

	'''
	There is an issue if the given array is full of 0 or NAs
	Maybe try passing a value error ? 
	'''

	Index_coef =  FindMode(Reference) / FindMode(Condition) ### that's the issue
	Condition_indexed = Condition*Index_coef
	return(Condition_indexed)


def compile_index(no_index_p,  output_p ,control):

	#Read the non indexed file
	no_indexed = pd.read_csv(no_index_p)


	#Get the name out
	colnames = list(no_indexed.columns.values)

	#Core file
	core_indexed = no_indexed.iloc[:, :12]

	#Libraries
	libraries = colnames[12:]


	#Get the chromosomes
	chromo = pd.unique(no_indexed.seqname)

	for my_library in libraries : 
		print(my_library)

		ctd_indexed = []
		for count, chromosome_name in enumerate(chromo):

			#Get the data for the right chromosome
			data_chromosome = no_indexed.loc[no_indexed['seqname'] == chromosome_name]

			#Get the data to work with
			ctrl_data = data_chromosome[control].to_numpy()
			cdt_data = data_chromosome[my_library].to_numpy()

			#Index

			#I need to think of a foolproof solution
			indexed = Indexing(ctrl_data,cdt_data) ###### here
			indexed_list = indexed.tolist()
			for element in indexed_list : ctd_indexed.append(element)


		core_indexed[my_library] = ctd_indexed

	core_indexed.to_csv(output_p)	
	print("Indexed table is printed")
