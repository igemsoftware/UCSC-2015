#!/usr/bin/env python3
#Name: Raymond Bryan (rvbryan)
#Created: 6/30/2015 Last Updated: 9/17/2015
'''
The following program will take a FASTA file given as an input
by the user. The file is then taken and seperated by the FastAreader
program given in the accompanied SequenceAnalysis.py file, to
header + sequence. Each sequence is then taken where each codon is
analysed and acounted for. From the gathered data information Codon
Frequency and Percent/Thousand is calculated. The gathered and calculated
data is then displayed and written to a seperate output(.dic) file
for a Genome of seperate genes.  

user inputs: Genome FASTA File(.fa) in a gene format
program outputs: A formated print table and .dic file as

print: AA   Codon   #count  Percent/Thousand   CodonFrequncy 
'''
################Start: psuedo code#################################
# 1. Take the FastaAreader program to read each gene sequence     #
#                                                                 #
# 2. Take the sequence and document/count each codon              #
#                                                                 #
# 3. Access the count                                             #
#                                                                 #
# 4. Use count to calculate Frequncy -> store in similar fashion  #
#                                                                 #
# 5. Use count to calculte the total time an AA is called         #
#                                                                 #
# 6. Use total to calculate codon frequency                       #
#                                                                 #
# 7. Use total the create a Grade Total                           #
#                                                                 #
# 8. Use count and Grand Total to calculate Percent per Thousand  #
#                                                                 #
# 9. Print Results in a output file as                            #
#     -> AA   Codon #count  Percent/Thousand   CodonFrequncy      #
#                                                                 #
################End: psuedo code###################################

class CodonFreq:
    """
    This class will give the frequentcies of each condon in the form a 
    dictionary. The dictionary will be coded from a check system that is provided 
    """
    rnaCodonTable = {
	# Second Base
	     # U             C             A             G
	#U
	'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
	'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
	'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
	'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
	#C 
	'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
	'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
	'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
	'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
	#A
	'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
	'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
	'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
	'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
	#G
	'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
	'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
	'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
	'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
                }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    def __init__ (self):
        '''In this method we will be able to analyse the sequences and then
        calculating the frequencies of each codon. '''
        #The following creates a three clear AA list of all possible
        #codon dictionaries as [AA[{codon:#}]]

        
        self.CodonCount = {}
        #Dictionary used to store the # of times the codon is called for in the Genome

        for codon, AA in CodonFreq.dnaCodonTable.items():
            try:
                self.CodonCount[AA].append( {codon : 0})
            except KeyError:
                self.CodonCount[AA] = [{codon: 0}]

        
        self.CodonFrequen = {}
        #Dictionary used to store the Frequency of which the codon appears in the Genome
        for codon, AA in CodonFreq.dnaCodonTable.items():
            try:
                self.CodonFrequen[AA].append( {codon : 0})
            except KeyError:
                self.CodonFrequen[AA] = [{codon: 0}]

        self.CodonPerThou = {}
        #Dictionary used to store the %/Thousand calculated for the codon in the Genome
        for codon, AA in CodonFreq.dnaCodonTable.items():
            try:
                self.CodonPerThou[AA].append( {codon : 0})
            except KeyError:
                self.CodonPerThou[AA] = [{codon: 0}]
                

    def readSeq(self, OriginalSeq):
        '''The following function allows the forward read to be	analysed of
        codons and append their count to a master list dictionary by matching
        the codon present to the possible codon in the dictionary. '''
        #Go over each nucleotide in the sequence 
        for i in range (0, len(OriginalSeq), 3):
            codon = OriginalSeq[i : i+3]

            #Find what AA the codon codes for
            try:
                Amino = CodonFreq.dnaCodonTable[codon]
            except KeyError:
                continue
            #Find how many codons code for that AA
            possibleCodons = len(self.CodonCount[Amino])

            #search possible Codons for the existing count
            for position in range(0,possibleCodons):
                #increment count of that codon by 1
                try:
                    #access the existing count of the codon
                    #by serching and matching the codon in the list
                    exCount = self.CodonCount[Amino][position][codon]
                    exCount += 1
                    #sets count as the new added one
                    self.CodonCount[Amino][position][codon] = exCount
                    #print ('total'+str(self.CodonCount[Amino][position][codon]))
                except KeyError:
                    #try the next one
                    continue
					
        return (self.CodonCount)
			
    def TableMaker(self):
        '''This program calculates the final Frequency of the
        codon as it appears in the genome.'''
        Total = []
        AAcount = {}
        i=0
        for AA in self.CodonCount.keys():

            #get total times the AA is coded
            for i in range(0, len(self.CodonCount[AA])):
                for codon in self.CodonCount[AA][i]:
                    Total.append(self.CodonCount[AA][i][codon])
                TotalPoss = sum(Total)
                AAcount[AA] = TotalPoss
            Total = []
            #print (AAcount)
        #for each codon within the AA group calc frequency
        #and replace the count for the frequency
        for AA in self.CodonCount.keys():
            for i in range(0,len(self.CodonCount[AA])):
                for codon in self.CodonCount[AA][i]:
                    number = self.CodonCount[AA][i][codon]
                    
                    try:
                        Freq = number/(AAcount[AA])
                        #print (Freq)
                        FreqPercent = float("{0:.2f}".format(Freq))
                        self.CodonFrequen[AA][i][codon] = FreqPercent
                    except ZeroDivisionError:
                        FreqPercent = 0


                
			
        return (self.CodonFrequen)

    def PercentPerThou(self):
        '''This portion withh take all the counted codons and place them in
    a sum that will be used to calculate percent per 1000 '''
        #add all counted codons together in a master Total
        Total = []
        AAcount = {}
        i=0
        MasterCount = 0
        for AA in self.CodonCount.keys():

            #get total times the AA is coded
            for i in range(0, len(self.CodonCount[AA])):
                for codon in self.CodonCount[AA][i]:
                    Total.append(self.CodonCount[AA][i][codon])
                TotalPoss = sum(Total)
                AAcount[AA] = TotalPoss
            Total = []
        for TotalCounts in AAcount.values():
            MasterCount += TotalCounts
        #from the Mastercount the Percent per thousand will be calc'd and stored
            
        for AA in self.CodonCount.keys():
            for i in range(0,len(self.CodonCount[AA])):
                for codon in self.CodonCount[AA][i]:
                    number = self.CodonCount[AA][i][codon]
                    
                    try:
                        PerThou = (number/ MasterCount)*1000
                        #print (PerThou)
                        PercentPerThou = float("{0:.2f}".format(PerThou))
                        self.CodonPerThou[AA][i][codon] = PercentPerThou
                    except ZeroDivisionError:
                        PercentPerThou = 0
        

import sequenceAnalysis
def main():
    '''Implements the Usage exception handler that can be raised from anywhere in process.'''
    import sys
    usage = "\nUsage: \npython3 %s FASTA.fa\n" % sys.argv[0]
    if len(sys.argv) != 2:
        print("Please enter a single FASTA file.")
        print(usage)
        sys.exit(3)
    OUTFILE = sys.argv[1].split(".")[0] + ".dic"
    FASTA = sys.argv[1]
    Class = CodonFreq()
    FileOutput = open(OUTFILE, 'w+')

    #IDEL TESTing code w/o cmd line
##    FASTA = input('file:')
##    FileOutput = open("ThisFileRightHere.fa", 'w+')
##    Class = CodonFreq()

    
    seqRead = sequenceAnalysis.FastAreader(FASTA)
    print('wait for it......')
 

    for header, Sequence in seqRead.readFasta():

        #print (header)
        Class.readSeq(Sequence)
        Class.TableMaker()
        Class.PercentPerThou()

    print('wait for it......')

    for AA in Class.CodonFrequen.keys():
            for i in range(0,len(Class.CodonFrequen[AA])):
                for codon in Class.CodonFrequen[AA][i]:
                    #print: codon, AA, Frequency, number

                    print('{:s} {:s} {:.0f} {:.3f} {:.3f}'.format(AA, codon, Class.CodonCount[AA][i][codon], Class.CodonPerThou[AA][i][codon], Class.CodonFrequen[AA][i][codon]))
                    print('{:s} {:s} {:.0f} {:.3f} {:.3f}'.format(AA, codon, Class.CodonCount[AA][i][codon], Class.CodonPerThou[AA][i][codon], Class.CodonFrequen[AA][i][codon]), file = FileOutput)
    #print (Class.CodonCount)		
 
if __name__ == "__main__":
	main()



		
		
	
	
	
