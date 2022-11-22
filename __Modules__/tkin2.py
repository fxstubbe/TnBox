
# import tkinter and all its functions
from tkinter import *
from tkinter import filedialog as fd
from PIL import ImageTk, Image
from Bio import SeqIO
import time
import subprocess
import threading
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
from __Modules__.Library import *
from __Modules__.Function import *
from __Modules__.Indexing import *
from __Modules__.TnSeq import *
from __Modules__.TnIF import *
from __Modules__.GFF_parser import *
from __Modules__.Table_c import *
from __Modules__.Artemis import *


        
class TnBox():

    def __init__(self, master):

    ##### -------- ---------------- ---------------- ---------------- -------- Class Methods -------- ---------------- ---------------- ---------------- -------- #####
        #Catch the data
        self.my_data = {}

        #--------------The list of references available
        my_databases = list(Path(f"{os.getcwd()}/Database/").glob('**/*'))
        self.my_data.update( { "Databases" : list(set([file.stem for file in my_databases]))} )

        #The list of files already processed and present in data
        my_ann = list(Path(f"{os.getcwd()}/Annotation/").glob('**/*'))
        self.my_data.update( { "Annotation" : list(set([file.stem for file in my_ann])) } )

        #The list of files already processed and present in data
        my_pics = list(Path(f"{os.getcwd()}/Graph/").glob('**/*'))
        self.my_data.update( { "Graph" : list(set([file.stem for file in my_pics])) } )

        self.my_data["Experiment"] = list()

        #Initiate an empty index key
        self.my_data["Index"] = ["--"]
        self.my_data["Index_path"] = list()

        #Initiate an empty index key
        self.my_data["Delta"] = ["--"]
        self.my_data["Delta_path"] = list()

        #Initiate an empty graph key
        self.my_data["Graph"] = [""]
        self.my_data["Graph_path"] = list()

        #Empty file to compare
        self.my_data["File_compare"] = ["", ""]

        #Empty artemis path
        self.my_data["Artemis_path"] = list()

        #
        self.my_data["File"] = [[""], [""]]

        #The list of files already processed and present in data
        self.update_mydata()

    ##### -------- ---------------- ---------------- ---------------- -------- GUI -------- ---------------- ---------------- ---------------- -------- #####

    # -------- -------- Define the working canvas -------- -------- #

        master.title("A TN-seq Toolbox")
        # Create left and right frames
        left_frame = Frame(master, width=200, height=620, bg='lightgrey')
        left_frame.grid(row=0, column=0, padx=10, pady=5, rowspan = 2)
        left_frame.grid_propagate(False)

        right_frame = Frame(master, width=400, height=450, bg='lightgrey')
        right_frame.grid(row=0, column=1, padx=10, pady=5)
        right_frame.grid_propagate(False)

        lower_frame = Frame(master, width = 400, height = 250, bg='lightgrey')
        lower_frame.grid(row=1, column=1, padx=10, pady=5)
        lower_frame.grid_propagate(False)

        image_frame =  Frame(master, width=510, height=430, bg='lightgrey')
        image_frame.grid(row=0, column=2, padx=10, pady=5, rowspan = 2)
        image_frame.grid_propagate(False)

    # -------- -------- Define sub panels -------- -------- #

        #------- The Library Processing panel
        # Create tool bar frame
        library_bar = Frame(left_frame, width=190, height=70)
        library_bar.grid(row=1, column=0, padx=5, pady=5)
        library_bar.grid_propagate(False)

        # Create tool bar frame
        tool_bar = Frame(left_frame, width=190, height=220)
        tool_bar.grid(row=2, column=0, padx=5, pady=5)
        tool_bar.grid_propagate(False)

        # Create tool bar frame
        file_bar = Frame(left_frame, width=190, height=220)
        file_bar.grid(row=3, column=0, padx=5, pady=5)
        file_bar.grid_propagate(False)

        # Create tool bar frame
        file_lower = Frame(left_frame, width=190, height=30)
        file_lower.grid(row=4, column=0, padx=5, pady=5)
        file_lower.grid_propagate(False)

        #------- The tools and metrics panel
 
        #GFF and metrisc
        gff_bar = Frame(right_frame, width = 300, height = 185)
        gff_bar.grid(row = 1, column = 0, padx=5, pady=5)

        #algorith choice
        gff_bar_list = Frame(right_frame, width = 390, height = 250)
        gff_bar_list.grid(row=3,column=0, padx=5, pady=5)
        gff_bar_list.grid_propagate(False)

        #------- The indexing panel

        index_bar = Frame(lower_frame, width = 390, height = 100)
        index_bar.grid(row = 1, column = 0, padx=5, pady=5)
        index_bar.grid_propagate(False)

        #------- The delta panel

        delta_bar = Frame(lower_frame, width = 390, height = 100)
        delta_bar.grid(row = 2, column = 0, padx=5, pady=5)
        delta_bar.grid_propagate(False)

        #------- An image panel

        image_artemis = Frame(image_frame, width = 500, height = 100)
        image_artemis.grid(row = 1, column = 0, padx=5, pady=5)
        image_artemis.grid_propagate(False)
        image_top = Frame(image_frame, width = 500, height = 130)
        image_top.grid(row = 2, column = 0, padx=5, pady=5)
        image_top.grid_propagate(False)
        image_bottom = Frame(image_frame, width = 500, height = 130)
        image_bottom.grid(row = 3, column = 0, padx=5, pady=5)
        image_bottom.grid_propagate(False)



    # -------- -------- Set up the Alignment panel -------- -------- #

        # Name the left frames
        Label(left_frame, text = "1. Process Libraries",  font=("Arial Bold", 15)).grid(row=0, column=0, padx=5, pady=5)

        # Add a checkbutton 
        Label(library_bar, text="1.1 Transposon fishing", font=("Arial Bold", 13)).grid(row=0, column=0, padx=5, pady=0)
        self.CheckVar= IntVar()
        self.check_library = Checkbutton(library_bar, text='miniTn5',variable=self.CheckVar, onvalue=1, offvalue=0)
        self.check_library.grid(column=0, row=1, sticky = W)
        self.check_library_mariner = Checkbutton(library_bar, text='mariner',variable=self.CheckVar, onvalue=2, offvalue=0)
        self.check_library_mariner.grid(column=0, row=2, sticky = W)


        # Add a listbox within the tool_bar
        Label(tool_bar, text="1.2 Select Reference Genome", font=("Arial Bold", 13)).grid(row=1, column=1, padx=5, pady=0)
        self.list_box = Listbox(tool_bar)
        self.list_box.grid(column=0, row=2, columnspan = 2, padx = 5)
        self.refresh_list(self.list_box, "Databases")
        self.btn_list_ref = Button(tool_bar, text="Add", command = lambda : self.get_Reference(self.btn_list_ref))
        self.btn_list_ref.grid(column=1, row=3)

         # Add a listbox within the file_bar
        Label(file_bar, text="1.3 Select Sequencing files", font=("Arial Bold", 13)).grid(row=0, column=0, columnspan = 2 , padx=5, pady=0)
        self.fastq_box = Listbox(file_bar)
        self.btn_fastq_add = Button(file_bar, text="Add", command = lambda : self.get_fastq(self.btn_fastq_rem))
        self.btn_fastq_rem = Button(file_bar, text="Remove", command = lambda: self.delete_item_listbox(self.fastq_box, "Experiment"))
        self.fastq_box.grid(column=0, row=1, columnspan = 2, padx = 5)
        self.btn_fastq_add.grid(column=0, row=2, columnspan = 1)
        self.btn_fastq_rem.grid(column=1, row=2)

        # Add a button on the bottom left
        self.button_get_tnseq = Button(file_lower, text = "Start Analysis", command = lambda : self.launch_tnseq(self.button_get_tnseq))
        self.button_get_tnseq.grid(row=0, column=0, padx=35)


    # -------- -------- Set up the GFF panel & Metrics  -------- -------- #

        # Labelling the r
        Label(right_frame, text = "2. Get transposon  insertion sites",  font=("Arial Bold", 15)).grid(row=0, column=0, padx=5, pady=0, columnspan = 3)

        
        #Subpanel labelling
        Label(gff_bar, text="2.1. Select Reference", font = ("Arial Bold", 13)).grid(row=0, column=0, columnspan = 2 , padx=5, pady=0)
        
        #Option box and adding GFF
        self.ann_set = StringVar()
        self.ann_set.set("    ")

        self.ann_optionlist = OptionMenu(gff_bar, self.ann_set, *self.my_data["Annotation"])
        self.ann_btn_add = Button(gff_bar, text = "Add .gff", command = lambda : self.get_gff(self.ann_btn_add))
        self.ann_optionlist.grid(row=1, column=0)
        self.ann_btn_add.grid(row=1, column=1)


        Label(gff_bar, text="2.2. Define metrics", font = ("Arial Bold", 13) ).grid(row=2, column=0, columnspan = 2 , padx=5, pady=0)
        #Metrics widget
        Label(gff_bar, text = "Trim end : ", font=("Arial", 13)).grid(row=3, column=0, sticky=W )
        self.trim_entry = Entry( gff_bar, width = 6)
        self.trim_entry.grid(row=3, column = 1)
        self.trim_entry.insert(0, 10)

        Label(gff_bar, text = "R window : ", font=("Arial", 13)).grid(row=4, column=0, sticky=W )
        self.Rwindow_entry = Entry( gff_bar,  width = 6)
        self.Rwindow_entry.grid(row=4, column = 1)
        self.Rwindow_entry.insert(0, 200)

        Label(gff_bar, text = "Slide : ", font=("Arial", 13)).grid(row=5, column=0, sticky=W )
        self.RSlide_entry = Entry( gff_bar, width = 6)
        self.RSlide_entry.grid(row=5, column = 1)
        self.RSlide_entry.insert(0, 5)

    # -------- -------- Set up the algorithm panel -------- -------- #

        #Subpanel labelling
        Label(gff_bar_list, text="2.3. Select files and algorithm", font = ("Arial Bold", 13)).grid(row=0, column=0 ,columnspan = 4, padx=5, pady=0)

        self.TA_box = Listbox(gff_bar_list, selectmode = "multiple")
        self.TnIF_box = Listbox(gff_bar_list, selectmode = "multiple")
        self.btn_tnif = Button(gff_bar_list, text = "TnIF", command = lambda : self.launch_compile_transposons(self.btn_tnif, self.TnIF_box, 1))
        self.btn_ta = Button(gff_bar_list, text = "TnIF" , command = lambda : self.launch_compile_transposons(self.btn_ta, self.TA_box, 1))
        self.btn_tnif_R100 = Button(gff_bar_list, text = "R Slide", command = lambda : self.launch_compile_transposons(self.btn_tnif_R100, self.TnIF_box, 2)) 
        self.btn_ta_R100 = Button(gff_bar_list, text = "R Slide",  command = lambda : self.launch_compile_transposons(self.btn_ta_R100, self.TA_box, 2))

        self.refresh_list(self.TA_box, "TA")
        self.refresh_list(self.TnIF_box, "TnIF")

        Label(gff_bar_list, text = "TA (Anchor site)").grid(row=1, column=0, columnspan=2 ,padx=0, pady=0)
        self.TA_box.grid(row=2, column = 0, columnspan=2 ,padx=5, pady=0)
        self.btn_ta.grid(row=3, column = 0, padx=5, pady=0 )
        self.btn_ta_R100.grid(row=3, column = 1, padx=5, pady=0 )

        Label(gff_bar_list, text = "TnIF (Coverage)").grid(row=1, column=2, columnspan=2 ,padx=0, pady=0)
        self.TnIF_box.grid(row=2,column=2, columnspan=2 ,padx=5, pady=0)
        self.btn_tnif.grid(row=3, column = 2, padx=5, pady=0 )
        self.btn_tnif_R100.grid(row=3, column = 3, padx=5, pady=0 )


    # -------- -------- Set up the Indexing panel -------- -------- #

        Label(lower_frame, text = "3. Indexing & Delta (TnIF)",  font=("Arial Bold", 15)).grid(row=0, column=0, padx=5, pady=0, columnspan = 2)


        Label(index_bar, text = "3.1 Index on a reference ",  font=("Arial Bold", 14)).grid(row=0, column=0, sticky=W)

        # Open file of interest
        Label(index_bar, text = "Choose file : ",  font=("Arial ", 13)).grid(row=1, column=0, sticky=W)
        self.file_to_ref_add = Button(index_bar, text="Load file", command = lambda : self.Load_OptionMenu(self.file_to_ref_add, self.my_data["Index"], self.my_data["Index_path"], [(self.index_optionlist,self.index_set )], [self.file_to_ref_del, self.file_to_ref_add, self.index_btn]))
        self.file_to_ref_add.grid(column=1, row=1, padx=5, pady=0)

        # remove file of interest
        self.file_to_ref_del = Button(index_bar, text="Unload", state =DISABLED , command = lambda: self.Unload_OptionMenu(self.file_to_ref_add, self.my_data["Index"], self.my_data["Index_path"], [(self.index_optionlist,self.index_set )], [self.file_to_ref_del, self.file_to_ref_add, self.index_btn]))  
        self.file_to_ref_del.grid(column=2, row=1, padx=5, pady=0)

        # OptionMenu to select the library you wanna use as index
        Label(index_bar, text = "Select control : ",  font=("Arial ", 13)).grid(row=2, column=0, sticky=W)
        self.index_set = StringVar()
        #self.index_set.set("--")
        self.index_optionlist = OptionMenu(index_bar, self.index_set, *self.my_data["Index"])
        self.index_optionlist.grid(row=2,column=1,columnspan=2)

        #Button to start the indexing
        Label(index_bar, text = "Index : ",  font=("Arial ", 13)).grid(row=3, column=0, sticky=W)
        self.index_btn = Button(index_bar, text = "Index file", state=DISABLED ,command = lambda : self.launch_indexing(self.index_btn))
        self.index_btn.grid(row=3, column=1, columnspan = 2)


    # -------- -------- Set up the Delta panel -------- -------- #

        Label(delta_bar, text = "3.2 Compute a delta",  font=("Arial Bold", 14)).grid(row=0, column=0, sticky=W)

        # Open file of interest
        Label(delta_bar, text = "Choose file : ",  font=("Arial ", 13)).grid(row=1, column=0, sticky=W)
        self.file_to_delta_add = Button(delta_bar, text="Load file", command = lambda : self.Load_OptionMenu(self.file_to_delta_add, self.my_data["Delta"], self.my_data["Delta_path"], [(self.delta_optionlist, self.delta_set)], [self.file_to_delta_del,self.file_to_delta_add, self.delta_btn]))
        self.file_to_delta_add.grid(column=1, row=1, padx=5, pady=0)

        # remove file of interest
        self.file_to_delta_del = Button(delta_bar, text="Unload", state =DISABLED , command = lambda: self.Unload_OptionMenu(self.file_to_delta_add, self.my_data["Delta"], self.my_data["Delta_path"], [(self.delta_optionlist, self.delta_set)], [self.file_to_delta_del,self.file_to_delta_add, self.delta_btn]))  
        self.file_to_delta_del.grid(column=2, row=1, padx=5, pady=0)

        # OptionMenu to select the library you wanna use as index
        Label(delta_bar, text = "Select control : ",  font=("Arial ", 13)).grid(row=2, column=0, sticky=W)
        self.delta_set = StringVar()
        #self.index_set.set("--")
        self.delta_optionlist = OptionMenu(delta_bar, self.delta_set, *self.my_data["Delta"])
        self.delta_optionlist.grid(row=2,column=1,columnspan=2)

        #Button to start the indexing
        Label(delta_bar, text = "Delta : ",  font=("Arial ", 13)).grid(row=3, column=0, sticky=W)
        self.delta_btn = Button(delta_bar, text = "Compute Delta", state=DISABLED ,command = lambda : self.launch_delta(self.delta_btn)) # Not good algorithm
        self.delta_btn.grid(row=3, column=1, columnspan = 2)

    # -------- -------- Set up the Graph panel -------- -------- #

        Label(image_frame, text = "4. Explore (Graphical options)",  font=("Arial Bold", 15)).grid(row=0, column=0, padx=5, pady=5, columnspan = 2)
        

        #-------- 2 Libraries from a single file

        Label(image_top, text = "Explore up to 2 libraries from a single file",  font=("Arial Bold", 14)).grid(row=0, column=0, columnspan = 2 ,sticky=W)

        #Make a button to load file to graph
        Label(image_top, text = "Load file : ",  font=("Arial", 13)).grid(row=1, column=0,sticky=W)
        self.file_to_graph_add = Button(image_top, text="Load file", command = lambda : self.Load_OptionMenu(self.file_to_graph_add, self.my_data["Graph"], self.my_data["Graph_path"], [(self.pic_optionlist_1, self.pic_set_1), (self.pic_optionlist_2, self.pic_set_2)], [self.file_to_graph_remove,self.draw_top, self.file_to_graph_add]))
        self.file_to_graph_add.grid(row=1, column = 1, padx=5, pady=0)

        #Make a button to load file to graph
        self.file_to_graph_remove = Button(image_top, text="Unload", state =DISABLED ,command = lambda : self.Unload_OptionMenu(self.file_to_graph_add, self.my_data["Graph"], self.my_data["Graph_path"], [(self.pic_optionlist_1, self.pic_set_1), (self.pic_optionlist_2, self.pic_set_2)], [self.file_to_graph_add], [self.file_to_graph_remove,self.draw_top]))
        self.file_to_graph_remove.grid(row=1, column = 2, padx=5, pady=0)

        #Load a file to graph from
        Label(image_top, text = "Select library (1): ",  font=("Arial", 13)).grid(row=2, column=0, sticky=W)
        self.pic_set_1 = StringVar()
        self.pic_optionlist_1 = OptionMenu(image_top, self.pic_set_1, *self.my_data["Graph"])
        self.pic_optionlist_1.grid(row=2,column=1,padx=5, pady=0)

        #Load a file to graph from
        Label(image_top, text = "Select library (2): ",  font=("Arial", 13)).grid(row=3, column=0, sticky=W)
        self.pic_set_2 = StringVar()
        self.pic_optionlist_2 = OptionMenu(image_top, self.pic_set_2, *self.my_data["Graph"])
        self.pic_optionlist_2.grid(row=3,column=1,padx=5, pady=0)

        #Draw the graph
        self.draw_top = Button(image_top, text="Draw Graph", state =DISABLED ,command = lambda : self.Graph_library(self.draw_top,self.pic_set_1, self.pic_set_2 ))
        self.draw_top.grid(row=4, column = 0, padx=5, pady=0)


        #-------- 2 libraries from 2 different files

        Label(image_bottom, text = "Compare 2 libraries from 2 different files (e.g. RSlide vs TnIF)",  font=("Arial Bold", 14)).grid(row=0, column=0, columnspan = 3,sticky=W)
        Label(image_bottom, text = "Files must use the same reference (e.g. Melitensis 16MM) ",  font=("Arial", 13)).grid(row=1, column=0, columnspan = 3,sticky=W)

        #Make a button to load file to graph
        Label(image_bottom, text = "Select File (1): ",  font=("Arial", 13)).grid(row=2, column=0, sticky=W)
        self.file_low_1_add = Button(image_bottom, text="Load file", command = lambda : self.Load_file_to_compare(self.file_low_1_add, 0, self.file_optionlist_1, self.file_set_1 ))
        self.file_low_1_add.grid(row=2, column = 1, padx=5, pady=0)
        self.file_set_1 = StringVar()
        self.file_optionlist_1 = OptionMenu(image_bottom, self.file_set_1, *self.my_data["File"][0])
        self.file_optionlist_1.grid(row=2,column=3,padx=5, pady=0)


        Label(image_bottom, text = "Select File (2): ",  font=("Arial", 13)).grid(row=3, column=0, sticky=W)
        self.file_low_2_add = Button(image_bottom, text="Load file", command = lambda : self.Load_file_to_compare(self.file_low_2_add, 1, self.file_optionlist_2, self.file_set_2 ))
        self.file_low_2_add.grid(row=3, column = 1, padx=5, pady=0)
        self.file_set_2 = StringVar()
        self.file_optionlist_2 = OptionMenu(image_bottom, self.file_set_2, *self.my_data["File"][1])
        self.file_optionlist_2.grid(row=3,column=3,padx=5, pady=0)

        #Draw the graph
        self.draw_bottom = Button(image_bottom, text="Draw Graph", state =NORMAL ,command = lambda : self.Graph_comparison(self.draw_top,self.file_set_1, self.file_set_2 ))
        self.draw_bottom.grid(row=4, column = 0, padx=5, pady=0)


        #-------- An artemis Panel
        Label(image_artemis, text = "Make an Artemis File",  font=("Arial Bold", 14)).grid(row=0, column=0, columnspan = 3,sticky=W)
        Label(image_artemis, text = "Convert a coverage file into an artemis compatible file) ",  font=("Arial", 13)).grid(row=1, column=0, columnspan = 3,sticky=W)

        Label(image_artemis, text = "Select covetage file : ",  font=("Arial", 13)).grid(row=2, column=0, sticky=W)
        self.file_artemis = Button(image_artemis, text="Load file", command = lambda : self.Load_artemis(self.file_artemis))
        self.file_artemis.grid(row=2, column = 1, padx=5, pady=0)

        self.draw_artemis = Button(image_artemis, text="Convert", state =DISABLED ,command = lambda : self.launch_artemis(self.draw_artemis))
        self.draw_artemis.grid(row=4, column = 0, padx=5, pady=0)

    ##### -------- ---------------- ---------------- ---------------- -------- Class Methods -------- ---------------- ---------------- ---------------- -------- #####

    # -------- -------- -------- -------- A few general class functions -------- -------- -------- -------- #

    def update_mydata(self):
        my_output = [file.stem for file in Path(f"{os.getcwd()}/data/").glob('**/*')]
        #print(my_output)
        self.my_data.update( {"TA" : list(filter(lambda x: x.endswith("TA"), my_output)) } )
        self.my_data.update( {"TnIF" : list(filter(lambda x: x.endswith("TnIF"), my_output)) } )
        #print(self.my_data)

    # ----- ----- ------ # LISTBOX functions

    def selected_item_listbox(self,lbox):

        #""" Class function returning the selected item in a listbox """

        items = [lbox.get(i) for i in lbox.curselection()]
        return items

    def selected_reference(self):

         #""" Single use Class function returning the selected fasta reference contained in listbox  """

        for i in self.list_box.curselection():
            return self.list_box.get(i)


    def delete_item_listbox(self, lb, key):
        for item in lb.curselection():
            temp = list(self.my_data[str(key)])
            del temp[item]
            self.my_data.update({str(key) : temp })
            self.refresh_list(lb, key)

    def refresh_list(self, lbox, key):

        #""" Class function refreshing the displayed list in the listbox"""

        lbox.delete(0, END)
        for item in self.my_data[str(key)]:
            lbox.insert("end", Path(item).stem)  

    def Change_color(self,lbox, position, color):

        #""" Class function Changing the background color of a given element in a listbox  """

        lbox.itemconfig(position, {'bg':color})    

    # ----- ----- ------ # OPTIOMENU functions

    def update_option_menu(self, om, val, val_set):
        menu = om["menu"]
        menu.delete(0, "end")
        for string in val:
            menu.add_command(label=string, command=lambda value=string: val_set.set(value))

    # -------- -------- -------- -------- Class function referring to path (get_file) -------- -------- -------- -------- #

    def get_Reference(self, button):

        self.files=fd.askopenfilename(initialdir="./Desktop/", title="Select file")
        if check_ext(Path(self.files), [".fasta" , ".fna"]) is True :

            #Update the reference
            self.my_data.update({"Reference" : self.files})
            make_reference(self.my_data["Reference"])

            #Update the database list
            created_database = Path(self.files).stem
            print(created_database)
            self.my_data["Databases"].append(created_database)

            #Refresh list_box
            self.refresh_list(self.list_box, "Databases")
        else : 

            print("Not a valid file ! Database was not created !")

    def get_fastq(self, button):

        self.files=list(fd.askopenfilenames(initialdir="/", title="Select file"))
 
        if len(self.files) == 0 :
            return
        else :

            #Check that all files are fastq
            my_paths = list(map(is_fastq, [self.files, [".fastq", ".fastq.gz"]]))
            print(my_paths)
            if all(my_paths) :

                #What if there were files in experiment already ?

                if len(list(self.my_data["Experiment"])) > 0 :
                    self.files.extend(list(self.my_data["Experiment"]))

                #Catch files
                self.my_data.update({"Experiment": list(set(self.files))})
                print(self.my_data["Experiment"])

                #Refresh the listbox
                self.refresh_list(self.fastq_box, "Experiment")
            else :
                return

    # -------- -------- -------- -------- Class function referring to path (get_file) -------- -------- -------- -------- #

    def get_gff(self, button):
        self.files=fd.askopenfilename(initialdir="/", title="Select file")
        self.my_data["Annotation"].append(Path(self.files).stem)
        parse_gff(self.files)
        self.update_option_menu(self.ann_optionlist, self.my_data["Annotation"], self.ann_set )

    # -------- -------- -------- -------- Class function referring to path (get_file) -------- -------- -------- -------- #


    # -------- -------- Load/unload files -------- -------- #

    #Options Menus
    def Load_OptionMenu(self, button, key_file, key_path, tuple_option_list, buttons_switch) : 

        #Make sure the catching keys/arguments are empty
        key_file.clear()
        key_path.clear()

        #Get the file to read
        my_data = fd.askopenfilename(initialdir=os.getcwd(), title="Select file")
        data = pd.read_csv(my_data)

        #Store the Index path
        key_path.append(my_data)

        #Get the columns
        my_files = list(data.columns.values)[12:]

        if len(my_files) < 1 : 
            print("Not a valid file")

        else : 

            #Update the key_file argument
            key_file.extend(my_files) 

            print(key_file)
            print(key_path)
     

            #Update the option_menu(s) 
            for tp in tuple_option_list : 
                option_list, option_set = tp
                self.update_option_menu(option_list, key_file, option_set)

            #Change button state
            for btn in buttons_switch : switchButtonstate(btn)

    def Unload_OptionMenu(self, button, key_file, key_path, tuple_option_list, buttons_switch):

        #Make sure to empty the dictionnary from previous input files
        key_file.clear()
        key_path.clear()

        #Update the option Menu
        for tp in tuple_option_list : 
            option_list, option_set = tp
            self.update_option_menu(option_list, key_file, option_set)
            option_set.set("   ")

        #Change button state
        for btn in buttons_switch : switchButtonstate(btn)

    # Two files for the comparison grahs    
    def Load_file_to_compare(self, button, position, optionlist, set_var):

        #Make sure to empty the dictionnary from previous input files
        my_file = fd.askopenfilename(initialdir="/", title="Select file")
    
        #Read the file
        data = pd.read_csv(my_file)

        #Get the colnames of interest
        my_files = list(data.columns.values)[12:]

        self.my_data["File_compare"][position] = my_file
        self.my_data["File"][position] = my_files
        
        self.update_option_menu(optionlist, self.my_data["File"][position], set_var)

    # an artemis    
    def Load_artemis(self, button):

        self.my_data["Artemis_path"].clear()
        my_data = fd.askopenfilename(initialdir=f"{os.getcwd()}/data/", title="Select file")
        self.my_data["Artemis_path"].append(my_data)
        button.config(state=DISABLED, text = "Loaded")
        self.draw_artemis.config(state=NORMAL)

    # -------- -------- -------- -------- Class function referring to the aligment -------- -------- -------- -------- #

    def tn(self, button, reference):

    #"""
    #Class function calling the alignment functions
    #"""

        #This computes the Tnseq
        Reference = self.selected_reference()
        trim_trans = self.CheckVar.get()

        for id,vale in enumerate(self.my_data["Experiment"]):
            self.Change_color(self.fastq_box, id, "orange")
            tn_seq(self.my_data["Experiment"][id], reference , 8, trim_trans)
            self.Change_color(self.fastq_box, id, "green")
            self.update_mydata()
            self.refresh_list(self.TA_box, "TA")
            self.refresh_list(self.TnIF_box, "TnIF")
            #self.update_option_menu(self.reftnif_optionlist, self.my_data["TnIF"], self.reftnif_set )
            #self.update_option_menu(self.refta_optionlist, self.my_data["TnIF"], self.refta_set )
        button.config(state =NORMAL)

    def launch_tnseq(self,button):

    #"""
    #Function wrapping tn. It gets the necessary argument and launch it on a different thread
    #"""

        button.config(text = "Running")

        #Get the selected reference
        ref = self.selected_reference()

        #Launch tn within its own thread
        threading.Thread(target = self.tn, args=[button, ref]).start()     

        #Inactivate button   
        self.btn_fastq_add.config(state=DISABLED)
        self.btn_fastq_rem.config(state=DISABLED)
        button.config(state=DISABLED)


    # -------- -------- -------- -------- Class function referring to TnIF processing -------- -------- -------- -------- #

    def compile_transposons(self, button, files, annotation, output_path, graph_path, method, trim, Rwindow, Slide):

       # """ Class function calling the TnIF/R100 algorithm """

        # Compute files
        button.config(state=DISABLED)  

        #Makes a non indexed table
        compile_table(files, annotation ,output_path, graph_path, method, trim, Rwindow, Slide)

        #Makes an indexed table
        #compile_index(files, annotation, no_index_p, index_p , CONTROL) #output_pth = no index file  
 
        #Enable the use to start again 
        button.config(state=NORMAL)



    def launch_compile_transposons(self, button, lbox, method): #Needs to add a selection for the reference (used to index)

        # """ Class function wrapping compile_transposons """


        # Let's check that the metrics have been properly set
        # If not, break the function

        my_input_metrics = self.trim_entry.get(), self.Rwindow_entry.get(), self.RSlide_entry.get()
        metrics = ConvertTuple(my_input_metrics)
        
        if metrics is None : 
            print("Metrics are unvalid")

        else : 

            #Check the method
            if method == 1 : 
                my_method = "TnIF"

            else : 
                my_method = "R100"


            # Init variables 
            gff_n = self.ann_set.get()
            my_files = self.selected_item_listbox(lbox) 
            gff = f"{os.getcwd()}/Annotation/{gff_n}.csv"
            #no_index_p = f"{os.getcwd()}/Sum_{gff_n}_{my_method}_noIndex.csv"
            index_p = f"{os.getcwd()}/Sum_{gff_n}_{my_method}_Index.csv"
            graph_p = f"{os.getcwd()}/Graph/"
            output_p = fd.asksaveasfilename(confirmoverwrite=False, title = "Select file",filetypes = (("CSV Files","*.csv"),))

            print(output_p)
            

        #Check a save file has been selected and start computing
        if len(output_p) == 0:  # user selected file
            return
        else : 
            #Launch the TnIF computation
            threading.Thread(target = self.compile_transposons, args=[button, my_files, gff, output_p, graph_p ,method, metrics[0], metrics[1], metrics[2]]).start() 


    # -------- -------- -------- -------- Class function referring to TnIF processing -------- -------- -------- -------- #

    def launch_indexing(self, button):

        #Deactivate the button while indexing
        button.config(state=DISABLED)

        control = self.index_set.get()
        no_index_p = self.my_data["Index_path"][0] #that's no good
        #output_p = no_index_p.replace("noIndex", "Index")
        #output_p = no_index_p.replace(".csv", f"_{control}.csv")

        output_p = fd.asksaveasfilename(confirmoverwrite=True, title = "Select file",filetypes = (("CSV Files","*.csv"),))

        #Check a save file has been selected and start computing
        if output_p is None:  # user selected file
            button.config(state=NORMAL)
            return
        else : compile_index(no_index_p,  output_p ,control)

        #Reactivate the button
        button.config(state=NORMAL)

    def launch_delta(self, button):

        #Deactivate the button while indexing
        button.config(state=DISABLED)

        control = self.delta_set.get()
        input_p = self.my_data["Delta_path"][0] #that's no good

        output_p = fd.asksaveasfilename(confirmoverwrite=True, title = "Select file",filetypes = (("CSV Files","*.csv"),))

        #Check a save file has been selected and start computing
        if output_p is None:  # user selected file
            button.config(state=NORMAL)
            return
        else : Delta(input_p, output_p ,control)

        #Reactivate the button
        button.config(state=NORMAL)


    # -------- -------- -------- -------- Class function for the artemis processing -------- -------- -------- -------- #


    def artemis(self, input_p, output_p):

        #Launch function
        Artemis(input_p, output_p)

        #Empty the Artemis path
        self.my_data["Artemis_path"].clear()

        #Change the button state
        self.file_artemis.config(state=NORMAL, text = "Load")


    def launch_artemis(self, button):

        #Get the file to index
        input_p = self.my_data["Artemis_path"][0]

        #Check the output
        output_p = fd.asksaveasfilename(confirmoverwrite=True, title = "Select file",filetypes = (("TXT Files","*.txt"),))


        #Check a save file has been selected and start computing
        if output_p is None:  # user selected file
            button.config(state=NORMAL)
            return
        else : 
            #Start here
            button.config(state=DISABLED)
            lala = [input_p, output_p]
            print(lala)
            threading.Thread(target = self.artemis, args=lala).start() 



     # -------- -------- -------- -------- Make some Graph -------- -------- -------- -------- #

   # def Create_figure(self, button):

   #     return Figure

    def Graph_library(self, button, set_var_1, set_var_2):

        #Deactivate the button while indexing
        #button.config(state=DISABLED)

        my_data = pd.read_csv(self.my_data["Graph_path"][0])

        library_1 = set_var_1.get()
        library_2 = set_var_2.get()


        if len(library_1) and len(library_2) > 5 : 

            g = sns.FacetGrid(my_data, col = "seqname")
            g.map(sns.scatterplot, library_1, library_2, size=0.5 , alpha=.3)
            plt.show()
            self.pic_set_2.set("     ")

        else : 

            f = sns.FacetGrid(my_data, col = "seqname")
            f.map(sns.histplot,  library_1,binwidth=0.1, kde=True)
            plt.show()

    def Graph_comparison(self, button ,set_var_1, set_var_2) : 

        #Get the data to graph
        data_1 = pd.read_csv(self.my_data["File_compare"][0])
        data_2 = pd.read_csv(self.my_data["File_compare"][1])

        #Select the column to graph
        library_1 = set_var_1.get()
        library_2 = set_var_2.get()


        #Substet those columns
        data_1_sub = data_1[["seqname", "locus_tag" ,library_1]].rename(columns = {library_1:'my_x'})
        data_2_sub = data_2[["seqname", "locus_tag" ,library_2]].rename(columns = {library_2:'my_y'})

        #Merge into a new dataframe
        data = data_1_sub.join(data_2_sub.set_index(['locus_tag', 'seqname']), on=['locus_tag', 'seqname'])

        #Grpah it
        g = sns.FacetGrid(data, col = "seqname")
        g.map(sns.scatterplot, "my_x" , "my_y" ,size=0.5 , alpha=.3)
        g.set(xlabel=library_1, ylabel=library_2)
        plt.show()

    ##### -------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- -------- #####


