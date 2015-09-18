#!/usr/bin/env python3
#Name: Cristian Camacho (cacamach@ucsc.edu)
#Group Members: David Bernick (dbernick@soe.ucsc.edu)

"""
The MIT License (MIT)

Copyright (c) <2015> <David L. Bernick, Cristian Camacho>

This module, called sequenceAnalysis, includes the following classes: ProteinParam,
and FastAreader. The ProteinParam class calculates and outputs various characteristics of
a given amino acid sequence. The FastAreader class, designed by Professor David Bernick of
UC Santa Cruz, provides the ability to read through multiple given FASTA sequences, and parse
the header from the sequence. It also allows for command line usability. 

"""

class ProteinParam:
    """Creates a ProteinParam object that which has access to the below dictionaries and methods.
    Attributes:
     attr1 (dict): Dictionary of molecular weights of amino acids
     attr2 (float): Molecular weight of H20
     attr3 (dict): Dictionary of absorbance values
     attr4 (dict): Dictionary of positive charges
     attr5 (dict): Dictionary of negative charges
     attr6 (float): Charge of the N terminus
     attr7 (float): Charge of the C terminus
     attr8 (dict): Dictionary of valid aa's
    
    Methods:
    aaCount(): Returns aa count
    pI(): Returns float pH value
    aaComposition(): Returns validAA dictionary
    charge(): Returns float of net charge of a protein
    molarExtinction(): Returns float molar extinction coefficient
    massExtinction(): Returns float mass extinction
    molecularWeight(): Returns float molecular weight of protein
    """
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
 
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
    validAA = {} 

    def __init__(self, protein):
        """Takes in sequence from main input, creates a list of all strings and
        concatnates. Then creates and saves the aaComposition dictionary.
        Args:
         param1: String of amino acids 
        """
        l = ''.join(protein).split()
        protString = ''.join(l).upper() #removes spaces and makes uppercase

        for aa in self.aa2mw.keys(): #only counts aa, not unexpected characters
            self.validAA[aa] = protString.count(aa)

    def aaCount(self): 
        """Iterates through every key in validAA dictionary, adds its value to aaNum and returns aaNum
        Return:
         Float amino acid count value
        """
        aaNum = 0
        for aa in self.validAA.keys():
            aaNum += self.validAA[aa]
        return aaNum
    
    def pI (self): 
        """Estimates the theoretical isoelectric point by finding the pH that yeilds
        a neutral net charge (as close to zero as possible, accurate to two decimal places)
        Returns:
         Float of best pH
        """
        bestCharge = 100000000 #utilizes a large number so the conditional statement is always passed
        bestPH = 0
        PH = 0
        while PH < 14.01: 
            thisCharge = abs(self.charge(PH))
            if thisCharge < bestCharge :
                bestCharge = thisCharge
                bestPH = PH
            PH += 0.01 #interates through every pH in decimal places. Range() only accounts for integers
        return (bestPH)
    
    def aaComposition(self):
        """Returns:
        ValidAA dictionary created in the init method """
        return self.validAA
    
    def charge(self, pH):
        """Calculates the net charge on the protein at specific pH using pKa of each
        charged amino acid, Nterminus and Cterminus.
        Args:
         param1: pH values from pI
        Returns:
         Float charge of the protein
        """
        posChrg = 0
        negChrg = 0
        for aa in self.validAA.keys():
            if aa in self.aa2chargePos.keys():
                posChrg += self.validAA.get(aa)*((10**self.aa2chargePos.get(aa))/(10**self.aa2chargePos.get(aa)+ 10**pH)) 
        posChrg += ((10**self.aaNterm)/(10**self.aaNterm + 10**pH)) #takes Nterm into account
            
        for aa in self.validAA.keys():
            if aa in self.aa2chargeNeg.keys():
                negChrg += self.validAA.get(aa)*((10**pH)/(10**self.aa2chargeNeg.get(aa)+10**pH))
        negChrg += ((10**pH)/(10**self.aaCterm +10**pH)) #takes Cterm into account

        return posChrg - negChrg
    
    def molarExtinction(self):
        """Estimates molar extinction coefficient based on number and extinction coefficients of
        tyrosines, tryptophans, and cysteines. Assume all Cysteine residues are represented as Cystine.
        Return:
         float molar extinction value
        """
        molEx = 0
        for aa in self.validAA.keys():
            if aa in self.aa2abs280.keys():
                molEx += self.validAA.get(aa)* self.aa2abs280.get(aa)
        return molEx
                
    
    def massExtinction(self):
        """Calculated by dividing molar extinction by the molecular weight of corresponding protein.
        Returns:
         float mass Extinction value
        """
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0
    
    def molecularWeight(self):
        """Iterate through every aa in validAA, multiply every mw for the appropriate aa by its count
        in validAA, sum all products, and subtract total mw of H2O from the summation.
        Return:
         Float molecular weight
        """
        mw = 0
        totmwH2O = self.mwH2O * (self.aaCount() - 1)
        for aa in self.validAA.keys():
            mw += (self.validAA.get(aa)*self.aa2mw.get(aa))
        return mw - totmwH2O
    
    
import sys
class FastAreader :
    '''
    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):
 
    object attributes:
    fname: the initial file name
 
    methods:
    doOpen () : Checks if a file name is given, else takes in file from system.in
    readFasta() : returns header and sequence as strings.
    Author: David Bernick
    Date: April 19, 2013
    '''
    def __init__ (self, fname = ''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):           
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence. If filename was not provided, std.in is
        used. Whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:
        #with open(self.fname) as fileH:     #original
            # initialize return containers
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
 
            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence
 
# presumed object instantiation and example usage
# myReader = FastAreader ('testTiny.fa');
# for head, seq in myReader.readFasta() :
#     print (head,seq)

