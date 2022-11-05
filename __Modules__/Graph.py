import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path



def graph_library(cdt, library) :
	
	#cdt is a numpy array 
	binwidth = 0.1
	plt.hist(cdt, bins = np.arange(min(cdt), max(cdt) + binwidth, binwidth))
