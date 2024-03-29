import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from __Modules__.Library import *


##### -------- ---------------- ---------------- ---------------- -------- A Tnseq pipeline  -------- ---------------- ---------------- ---------------- -------- #####


def make_reference(filename):
    Ref = Reference(filename, "Ref")
    Ref.bwa_index()


def tn_seq(filename, reference, method, trim_trans,thread=4):

	mini_tn5 = "GGTTGAGATGTGTATAAGAGACAG"

    # Create a library object
	Exp = Library(filename, "Exp")

    #Trim adapter
	Exp.trim_adapters()

	#If needed, fish the transposon
	if trim_trans >  0 :

		if trim_trans == 1 :
			Exp.trim_transposon(mini_tn5)
		else :
			with open(f"{os.getcwd()}/custom_sequence.txt") as f:custom_sequence = f.read()

			Exp.trim_transposon(custom_sequence.strip())

	else : Exp.quality_trim()

	#Align to the reference
	Exp.bwa_aln(str(reference),thread)

	#Extract both the TnIF and TA site
	Exp.get_TnIF()
	Exp.get_TA()

	#Clear intermediate files 
	Exp.clean_temp()