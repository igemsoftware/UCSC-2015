UCSC iGEM 2015
======

# THIS PROJECT IS FROZEN AS OF 9-18-2015

Currently only works on Linux or Mac OS X. Some configuration is needed for Windows compatability

#### Dependencies

* [Python 3.4 or higher](https://www.python.org/downloads/)
* an sh-compatible shell (bash works just fine)
* IDLE (if working on widnows)

#### Description
We are submitting three unique programs that simplify and improve the process of protein analysis for organisms across various domains of life. These tools are useful not only for our wet lab research, but for any projects incorporating proteins from different organisms. Firstly, the ```sequenceAnalysis.py``` program allows for quick, large scale protein analysis to provide an efficient filtering process for identifying proteins of interest. Secondly, the ```CodonBiasGenerator.py``` produces a GCG codon bias table of an organism for the process of codon optimization. Finally, the Frequency Optimized Codon Usage Strategy (```FOCUS.py```) tool is a working prototype that allows for the detection of rare codons based on an organism’s codon bias. This tool provides a novel method for codon optimization with respect to rare codons, which may improve protein production without loss of enzymatic function.


#### Usage
User guides are located in the ```User_Guides``` folder

##### Setup
There is no setup required to run any of the python programs. 
The only setup needed is to install Python 3.4 on your machine



##### Pitfalls
```FOCUS.py``` is a working prototype and is not finished. A new scoring method needs to be created to reduce the amount of noise for detecting rare codons. 

