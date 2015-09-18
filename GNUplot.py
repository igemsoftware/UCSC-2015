'''
Author: Jairo Navarro
Date: July 23, 2015

Description:
This program will create a file that can be used in a 
GNUplot script. Currently, the program only contains the class that is
used in the FOCUS python program, and does not work on its own.

Formatting:
The first column will be the amino acid position.
    second column will be the amino acid
    third column will be the number 1, 2 or 4
        -Coil: 1
        -Helix: 2
        -Sheet: 4
    fourth column will be the percent per thousand for the codon

'''

class GNUmaker:

    def __init__(self,proteinSeq,freqSeq,SSpred,Thousand):
        #Protein sequence string
        self.proteinSeq = proteinSeq

        #FOCUS score string
        self.freqSeq = freqSeq
        #List with a string for the Secondary structure prediction as 
        # the first entry
        self.SSpred = SSpred
        #List with entries correspoding to the percent per thousand
        #The first character in the self.proteinSeq corresponds to the
        # first entry in the self.Thousand list. 
        self.Thousand = Thousand
        #Dictionary to convert between Secondary Structure character
        # to a number which can be used in GNUplot. 
        self.reverseTrans = {'C':1,'H':2,'E':4}

    def makeTable(self):
        '''
           Prints a table that can be used in GNUplot
        '''
        for count,b  in enumerate(zip(self.proteinSeq,self.freqSeq)):
            AA = b[0]
            freq = b[1]
            thousand = self.Thousand[count]
            
            try:
                SS = self.reverseTrans[self.SSpred[0][count]]
            except IndexError:
                continue
            print("%d  %s  %s  %d  %.2f" % (count+1,AA,freq, SS, thousand))


def main():
    pass

if __name__ == "__main__":
    main()
