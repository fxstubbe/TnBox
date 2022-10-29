
######
####
##
# Author : François-Xavier STUBBE 
# Institution : Université de Namur
#
##
####
######

##### -------- ---------------- ---------------- ---------------- -------- Libraries and packages needed  -------- ---------------- ---------------- ---------------- -------- #####

#Standard library packages

import os
import sys
import shutil 
import time 
import subprocess 
import threading 
from pathlib import Path
from tkinter import * 
from tkinter import filedialog as fd
import importlib.util

#Import my packages

from __Modules__.tkin2 import *
from __Modules__.Library import *
from __Modules__.Function import *
from __Modules__.Indexing import *
from __Modules__.TnSeq import *
from __Modules__.TnIF import *
from __Modules__.GFF_parser import *

##### -------- ---------------- ---------------- ---------------- -------- Libraries and packages needed  -------- ---------------- ---------------- ---------------- -------- #####


# I think the easiest is to ask the user to download anaconda


# -------- -------- Non standard libray packages -------- -------- #

# Doens't check if 
#packages = ['numpy', 'pandas', 'seaborn', 'pillow', 'seqio']


#for pckg in packages : 

#	if pckg in sys.modules:

#		print(f"{pckg!r} already in sys.modules")
#		pass

#	elif (spec := importlib.util.find_spec(pckg)) is not None:

#	# Install using conda
#		subprocess.run(f"conda install - anaconda {pckg}" , shell = True)
#
#		print(f"{pckg!r} has been imported")

#	else:
#		print(f"can't find the {pckg!r} module")


# -------- -------- Biological Tools & Git -------- -------- # Using ANACONDA


#How to use the API ? 

#The aligner
#subprocess.run("conda install -c bioconda bwa" , shell = True)

#The genome coverage
#subprocess.run("conda install -c bioconda bedtools" , shell = True)

#Git (to downaload BBMAP)
#subprocess.run("conda install -c anaconda git" , shell = True)


# -------- -------- BBAMP -------- -------- # USING GIT or TAR

###Clone BBMAP in the working directory

if os.path.exists(f"{os.getcwd()}/__Modules__/bbmap/") is False : 
	subprocess.run("git clone https://github.com/BioInfoTools/BBMap.git ./__Modules__/bbmap" , shell = True)

#Is than an ok way ? I could put a zip file and unzip it
#Extract tar.gz fil
#tar -xf ./ToolBox/Modules/BBMap_38.86.tar.gz -C ./ToolBox/Modules/







##### -------- ---------------- ---------------- ---------------- -------- Create the TnBox instance  -------- ---------------- ---------------- ---------------- -------- #####


root = Tk()
e = TnBox(root)
root.mainloop()

