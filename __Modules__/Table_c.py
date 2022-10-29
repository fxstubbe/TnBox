
import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from __Modules__.Indexing import *

class Table_c(object):

	def __init__(self, input):


		if Path(input).is_file():

			self.path = Path(input)
			#self.data = pd.read_csv(Path(input))

		else :
			print("File doesn't exist !")

	def col_names(self):
		data = pd.read_csv(self.path)
		colnames = list(data)[12:]
		print(colnames)
		return(colnames)


	def indexing(self, control):

		no_indexed = pd.read_csv(self.path) 

		#Get the chromosomes
		chromo = pd.unique(no_indexed.seqname)

		libraries = list(no_indexed)[12:]

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


			#core_indexed[my_library] = ctd_indexed

		#core_indexed.to_csv(output_p)	
		print("Indexed table is printed")
