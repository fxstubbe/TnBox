import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from __Modules__.Library import *


def Artemis(input_p, output_p):

	#
	# Transform a file into an Artemis readable file
	#
	#

	# Read in the table
	my_table = pd.read_table(input_p, names=["Chromosome", "Position", "Count"])

	#Get the pseudocount
	my_table["pseudoCount"] = np.log10(my_table["Count"] + 1)

	#Get the chromosomes to loop over
	Chro = my_table.Chromosome.unique()

	#Write out Artemis Files
	for ch in Chro :
		my_chro =  my_table[my_table.Chromosome.eq(ch)].loc[:, ["Position","pseudoCount"]].to_csv(Path(output_p.replace(".txt", f"_{ch}.txt")), sep = "\t", index = False)
