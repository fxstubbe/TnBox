# TnBox
## _A Tn-seq Toolbox_


TnBox is a user friendly toolbox for the analysis of high throughput mutant libraries (Tn-seq). 
It provides a set of tools allowing the user to parse and align sequencing reads on a prokaryotic genome. It then implements two algorithms (Rslide & TnIF) for the study of gene essentiality.  

For support, questions or requests, please contact: francois-xavier.stubbe@unamur.be

## Installation & quick start

TnBox comes as a simple python program that works on both macOS and Linux systems. However, a few dependencies need to be installed and set in your PATH. Similarly, a couple python packages need to be downloaded. There are many ways to install these but the eaisest way might be to use the python package manager [Anaconda][conda]. Since Anaconda comes along a version of python, it is not necessary to have python pre-installed on your machine.

#### Install dependencies with conda on linux/macOS systems

- [Burrow Wheeler Aligner][bwa] 
- [Samtools][samtools]
- [Bedtools][bedtools]
- [Git][git]

Paste the following commands in a terminal

```sh
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c anaconda git
```
#### Install python packages with conda on linux/macOS systems

- [Numpy][python]
- [Pandas][bwa] 
- [Matplotlib][bwa] 
- [Searborn][samtools]

Paste the following commands in a terminal

```sh
conda install -c anaconda numpy
conda install -c anaconda pandas
conda install -c anaconda matplotlib
conda install -c anaconda searborn
```

#### Setting up and testing TnBox

##### Clone TnBox

To install TnBox, clone this repository using git : 
```sh
git clone https://github.com/fxstubbe/TnBox
```
Alternatively, you can also download TnBox using the green button at page's top.

##### Setting up TnBox

If you correctly installed all the required repositories and packages, TnBox is ready to use as it is. To start the program, you need to cd into TnBox folder and TnBox.py.

Assuming you've installed TnBox on the desktop; paste the following commands :

```sh
cd Desktop/TnBox/
python TnBox.py
```

On the first use, TnBox will automatically clone [BBMAP][bbmap] into the TnBox folder. This steps is crucial and if not executed properly will lead to failure. If an interface opens up : congrats ! TnBox is now woriking on your computer

## Analysis

![](https://github.com/fxstubbe/TnBox/blob/main/Images/TnBox.png)

### 1. Process libraries

The first step int any high throughput sequencing project is to map the sequencing reads onto a genome of interest. To do so, reads must first be trimmed of adapter and quality trimmed. Then, the processed reads are aliged onto an indexed genome using an aligner. TnBox uses [BBMAP][bbmap] for the filtering and [bwa][bwa] for the alignment. 

#### Transposon fishing

We recommend sequencing of miniTn5 libraries with the [in-house transposon sequencing service][fasteris] provided by fasteris.

In the past, Tn-seq libraries where often sequenced on Illumina HiSeq intruments. With those machines a nested strategy using a specific primer matching the 3'end of the transposon. Therefore, the first read base was the transposon insertion site. On newer machines (e.g. NextSeq, NovaSeq), sequencing reactions with specific primers aren't reliable. Libraires are now sequenced from the illumina adapters, which then re-sequence the end of the transposon. To separate specific from aspecific reads, it is necessary to separate reads containing the transposon from reads that do not (up to 35% of the library). By checking the tickbox "miniTn5", TnBox will filter transposon containing reads and trim out the transposon sequence prior to mapping (bottom panel, image below). If omitted, lots of non-specific reads will map and weakened downstrean analysis. 

Both Mariner and miniTn5 are available to parse. However, this methode has only been used and confirmed on miniTn5 libraires.

![](https://github.com/fxstubbe/TnBox/blob/main/Images/aspecific_mapping.png)

#### Select or Add a reference

As previoulsy stated, sequencing reads need to be mapped onto a reference genome. Genomes (.fasta , .fna) can be dowaloaded from NCBI (e.g. [*Brucella abortus*][abortus]). A few references are provided with TnBox. Once added, your reference genome will appear in the list. Select the genome of interest (highlighted in blue) and proceed to the enxt step. 


#### Add sequencing files (.fastq.gz)

Add your sequencing reads. The files must use the extension .fastq.gz otherwise TnBox will fail. Tn-seq are usually sequence in single-end but, if your reads are paired-end, only provide the forward reads (R1) as this is where the transposon lies. Similarly, sometimes, sequencing reads will be split over mutiple file. You can simply concatenate them together as follow : 

```sh
cat file_1 file_2 > concatenated_file
```

#### Start the anaysis

Click on start. That's it. If nothing starts, make sure all the previous steps are correctly fulfilled.

This process is time consuming. It largely depends on the amount of reads, the genome size ... but also your machine performance. Be patient, make yourself a cofee and come back later. TnBox provides a visual cue of where it is in the process. When the library name is orange, it's under process. When it turns green, TnBox is done with this library.

*FYI : TnBox will produced both an aligmnet file (.bam) and two coverage file (more details below). Those are respectively stored in ./TnBox/bam and ./TnBox/data*

### 2. Get transposons insertion sites 

#### Choose or Add a reference (.gff)
Now that the reads have been aligned to the genome of interest, it's time extract where the transposons are. The first step is to select the appropirate reference or the add your reference of choice in gff format. For example, you can dowload Brucella abortus reference [here][abortus]. 

#### Define metrics for the algorithms

`Trim end` 

Since many genes can tolerate insertions in their 5' or 3', it is wise to remove those extremities while looking for insertion sites. By default, TnBox trims 10% of the transcript length on both 5' and 3' ends. 
This parameter is used for both the Rslide and TnIF algorithms (see below)

`R window`, `Slide`

While using the RSlide algorithm, a sliding window of size n (`R window`) is moved with an increment p (`Slide`) over the genome. By default, the R window is set to be 200 as described in [Sternon et al.,][sternon] and the increment is set to 5. 


#### Select files and algorithm

When processing the libraries, TnBox generated 2 types of files :
i. TA files, where only 5' end of each mapped read has been counted. This is equivalent to the exact insertion site. This is the recommended method when looking for essential genes.
ii. TnIF files, where the whole mapped read has been counted. Even if this method is less sensitive to detect essential genes, it performs better when using the TnIF algorithm which has for goal to detect essentiality variations across several conditions as described in [Potemberg et al., ][potemberg] 



`Rslide`

Briefly; the coverage file is split into windows matching the input metrics. For example a 3 278 307 bp genome is split into 655,662 windows. For each window, the coverage sum is computed . Each gene is then assigned an *Essentiality score* equivalent to the number of windows having a sum of 0 (aka, no transposons jumped in that genomic region). 


`TnIF`



### 3. Indexing

### 4. Explore

![](https://github.com/fxstubbe/TnBox/blob/main/Images/com.png)

## License

MIT

**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[python]: https://www.python.org/
[conda]: https://www.anaconda.com/
[git]: https://git-scm.com/
[pandas]: https://pandas.pydata.org/
[numpy]: https://numpy.org/
[plt]: https://matplotlib.org/
[sns]: https://seaborn.pydata.org/
[bbmap]: https://github.com/BioInfoTools/BBMap
[abortus]:https://www.ncbi.nlm.nih.gov/genome/520
[fasteris]:https://www.fasteris.com/en-us/NGS/DNA-sequencing/Tn-Seq
[sternon]:https://journals.asm.org/doi/10.1128/IAI.00312-18
[potemberg]: https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1010621
