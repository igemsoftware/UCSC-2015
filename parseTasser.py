'''
Author: Jairo Navarro
Date: July 23,2015

The MIT License (MIT)

Copyright (c) 2015 Jairo Navarro

Python program to parse I-TASSER secondary structure format files
The program as currently written returns the 3rd and 4th column from
the I-TASSER file, which is used in the FOCUS python program.
'''

class parseTASSER:

    def __init__(self,TASSER):

        #Saves file name given to the program
        self.TASSER = TASSER
        self.trans = {1:'C',2:'H',4:'E'}

    def OpenFile (self,thisFile):
        '''Opens the file and prints an error message if the file is not found.'''
        try:
            return open(thisFile,'r')
        except FileNotFoundError:
            print("\n%s does not exist.\n" % thisFile)

    def returnList (self):
        '''Returns a list with the first entry being the secondary
           structure string and the second entry being the 
           confidence for the secondary structure prediction.
        '''
        returnList = []
        #Open the secondary structure file
        tassFile = self.OpenFile(self.TASSER)
        #Secondary Structure string container
        SecondaryStructure = ''
        #Confidence for the secondary structure string container
        Confidence = ''
        # Read line by line in the file
        for line in tassFile:
            #Create a list of the current line
            line = line.split()
            #Translate the number in the file to a Secondary Structure character
            letter = self.trans[int(line[2])]
            #Confidence of secondary structure
            number = line[3]
            
            #Add the secondary structure character to the container
            SecondaryStructure += letter
            #Add the confidence to its container
            Confidence += number

        #Return the list
        return [SecondaryStructure,Confidence]

def main():
    import sys
    #File given as the argument after the program name in the command line
    FILE = sys.argv[1]
    parser = parseTASSER(FILE)
    LIST = parser.returnList()
    #Print both strings in the list
    for item in LIST:
        print(item)

if __name__ == "__main__":
    main()
