## Author: Ethan Moran
## Date: 16th March 2023

import tkinter as tk
from tkinter import filedialog
import sys


def main():
    print("creating file explorer window minimise spkder to see it. :....")
    root = tk.Tk()
    root.wm_withdraw() # this completely hides the root window
    # use windows explorer to input the file name
    FileName = filedialog.askopenfilename(filetypes = [('Fasta files','.fasta')])
    root.destroy()
    dnaSeq = "";    
    CodonList = [];
    CodonList2 = []
    CodonList3 = [];
    CodonListCompliment1 = []
    Comp1 = ''
    complimentList = []
    revCompList = []
    revCompList2 = []
    revCompList3 = []
    limit = 300
    
    
    
    CodonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    
    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #complement dictionary
   
    try:
        Fp1 = open(FileName,'r')

        #gets the first line
        first_line = Fp1.readline();
        for line in Fp1:
                dnaSeq += line; #adds all other lines to the dna seq
        

        ListSeq = dnaSeq.split('\n');
        compSeq = dnaSeq.replace('\n', '');
        newdnaSeq = ('').join(ListSeq);

        

        complimentDnaSeq = ('').join(compSeq)

        
        #gets the dna sequence used for the complement
        for i in range(0, len(complimentDnaSeq)):
            complimentList.append(complimentDnaSeq[i])
       
        
        #first reading frame
        for i in range(0,len(newdnaSeq), 3):
            CodonList.append(newdnaSeq[i:i+3])
            
        #second reading frame

        for i in range(1, len(newdnaSeq), 3):
                CodonList2.append(newdnaSeq[i:i+3])

        #third reading frame

        for i in range(2, len(newdnaSeq), 3):
            CodonList3.append(newdnaSeq[i: i+3])
                

        print("Description: " + first_line)

        
        
        i = 0
        
        #changes the list to a complement list 
        while (i < len(complimentList)):
            if complimentList[i] in complement:  #checks if the character is in the dictionary for complements
                complement1 = complement[complimentList[i]]   
                CodonListCompliment1.append(complement1)     #adds to the a new complement list if necessary
                Comp1 += complement1 #adds to a string to be reversed
                i += 1       
            else:
                complimentList.remove(complimentList[i])
                continue


        complimentRev = Comp1[::-1] #reverses the compliment


        #4th reading frame -1 
        for i in range(0, len(complimentRev),3):
            revCompList.append(complimentRev[i:i+3])    

        #5th reading frame -2
        for i in range(1, len(complimentRev),3):
            revCompList2.append(complimentRev[i:i+3])  
        
        #6th reading frame -3
        for i in range(2, len(complimentRev),3):
            revCompList3.append(complimentRev[i:i+3])  



        #translates the codon from a dna sequence to an amino acid sequence
        def CodonTranslation(ListCodon):
            CodonStrToList = []
            CodonStrToList = [ListCodon[i:i+3] for i in range(0, len(ListCodon), 3)] #ensures that the ListCodon entered is a list.
            AminoAcidList = []
            i = 0
            AminoAcidSeq = ''
            while (i < len(CodonStrToList)): #loops through the list
                if CodonStrToList[i] in CodonTable: #checks if codon in table
                    AminoAcid = CodonTable[CodonStrToList[i]]     
                    AminoAcidList.append(AminoAcid) #adds the amino acid to the list
                    AminoAcidSeq += AminoAcid #adds it to a string
                    i += 1         
                else:
                    CodonStrToList.remove(CodonStrToList[i]) #ensures no codons that are less than 3 get translated
                    continue
            return AminoAcidSeq #returns the string

        
        
        







        def orf(AminoList, RFNumber, ORFList1=None):
                    
                    #start and stop codons declared
                    startCodon = ['ATG']
                    stopCodon = ["TAA", "TAG", "TGA"]

                    #ensures that the Open reading frame list is empty as ORFs from other reading frames would appear here otherwise
                    if ORFList1 is None: 
                        ORFList1 = []

                    i = 0
                    #this is a flag used to ensure that the code does not add any characters between the start and stop codon in an ORF.If there is an M in an ORF this is ignored until a stop codon is found.
                    startFound = False

                    ORFEnd = 0
                    #loops through the AminoList provided
                    for i in range(0,len(AminoList)):
                        #checks if the character is an M
                        if AminoList[i] in startCodon:
                            #this checks if startFound is false which means a new orf has been started. This is because startFound would be false if no previous start codon has been found
                            if not startFound:
                                ORFStart = i
                                ORF = AminoList[i]
                                startFound = True
                            elif startFound:
                                 ORF += AminoList[i]
                        elif startFound:
                            
                                #this checks any character thats not a start codon to see if it is a stop codon. If it is, the characters from ORF are added to the list
                                if AminoList[i] in stopCodon:
                                    stop = AminoList[i]
                                    ORFEnd = i
                                    if (((ORFEnd - ORFStart)*3)+3) >= limit:
                                        if(RFNumber > 0):
                                            if(RFNumber == 1):
                                                ORFList1.append(("RF" + str(RFNumber), "NucleotideStart:  " + str((ORFStart*3)+1), "NucleotideEnd:  " + str((ORFEnd+1)*3),"aa:  " + str(ORFEnd - ORFStart), "nt: " + str(((ORFEnd - ORFStart) *3)+3), "ORF:  " + ORF+str(stop), "    ORF Amino Acids:  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                            elif(RFNumber == 2):
                                                ORFList1.append(("RF" + str(RFNumber), "NucleotideStart:  " + str((ORFStart*3)+2), "NucleotideEnd:  " + str(((ORFEnd)*3)+4), "aa:  " + str(ORFEnd - ORFStart), "nt: " + str(((ORFEnd - ORFStart) *3)+3), "ORF:  " + ORF + str(stop), "    ORF Amino Acids::  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                            elif(RFNumber == 3):
                                                ORFList1.append(("RF" + str(RFNumber), "NucleotideStart:  " + str((ORFStart*3)+3), "NucleotideEnd:  " + str(((ORFEnd)*3)+5), "aa:  " + str(ORFEnd - ORFStart), "nt: " + str(((ORFEnd - ORFStart) *3)+3), "ORF:  " + ORF + str(stop), "    ORF Amino Acids::  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                        else: 
                                            ORFStart = len(AminoList) - ORFStart
                                            ORFEnd = len(AminoList) - ORFEnd
                                            if(RFNumber == -1):
                                                ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)-1),  "NucleotideEnd:  " + str((ORFEnd*3)-3), "aa:  " + str(ORFStart-ORFEnd), "nt: " + str(((ORFStart-ORFEnd) *3)+3), "ORF:  " + ORF + str(stop), "    ORF Amino Acids::  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                            elif(RFNumber == -2):
                                                ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)-2),  "NucleotideEnd:  " + str(((ORFEnd)*3)-4), "aa:  " + str(ORFStart-ORFEnd), "nt: " + str(((ORFStart-ORFEnd) *3)+3), "ORF:  " + ORF +str(stop), "    ORF Amino Acids::  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                            elif(RFNumber == -3):
                                                ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)),  "NucleotideEnd:  " + str(((ORFEnd)*3)-2) ,"aa:  " + str(ORFStart-ORFEnd), "nt: " + str(((ORFStart-ORFEnd) *3)+3), "ORF:  " + ORF +str(stop), "    ORF Amino Acids::  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.  
                                    startFound = False
                              
                                
                                else:
                                    ORF += AminoList[i] # this fills up ORF with the characters between the start and stop codons.
                                        
                    
                                if startFound:
                                        if (((ORFEnd - ORFStart)*3)+3) >= limit:
                                                if(RFNumber > 0):
                                                    ORFEnd = len(AminoList)
                                                    if(RFNumber == 1):
                                                        ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)+1),  "NucleotideEnd:  " + str((ORFEnd+1)*3), "aa:  " + str(ORFEnd - ORFStart), "nt: " + str(((ORFEnd - ORFStart) *3)+3), "ORF:  " + ORF + str(stop), "    ORF Amino Acids::  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                                    elif(RFNumber == 2):
                                                        ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)+2),  "NucleotideEnd:  " + str(((ORFEnd)*3)+4), "aa:  " + str(ORFEnd - ORFStart), "nt: " + str(((ORFEnd - ORFStart) *3)+3), "ORF:  " + ORF+str(stop), "    ORF Amino Acids::  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                                    elif(RFNumber == 3):
                                                        ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)+3),  "NucleotideEnd:  " + str(((ORFEnd)*3)+5), "aa:  " + str(ORFEnd - ORFStart), "nt: " + str(((ORFEnd - ORFStart) *3)+3), "ORF:  " + ORF+str(stop), "    ORF Amino Acids:  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                                else: 
                                                    ORFStart = len(AminoList) - ORFStart
                                                    ORFEnd = len(AminoList) - ORFEnd
                                                    if(RFNumber == -1):
                                                        ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)-1),  "NucleotideEnd:  " + str(((ORFEnd)*3)-3), "aa:  " + str(ORFStart - ORFEnd), "nt: " + str(((ORFStart-ORFEnd) *3)+3), "ORF:  " + ORF+str(stop), "    ORF Amino Acids:  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                                    elif(RFNumber == -2):
                                                        ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)-2),  "NucleotideEnd:  " + str(((ORFEnd)*3)-4), "aa:  " + str(ORFStart - ORFEnd), "nt: " + str(((ORFStart-ORFEnd) *3)+3), "ORF:  " + ORF+str(stop), "    ORF Amino Acids:  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.
                                                    elif(RFNumber == -3):
                                                        ORFList1.append(("RF" + str(RFNumber),  "NucleotideStart:  " + str((ORFStart*3)),  "NucleotideEnd:  " + str(((ORFEnd)*3)-2), "aa:  " + str(ORFStart - ORFEnd), "nt: " + str(((ORFStart-ORFEnd) *3)+3), "ORF:  " + ORF+str(stop), "    ORF Amino Acids:  " + CodonTranslation(ORF))) # this adds the reading frame number, the start of the reading frame, end of the reading frame and the ORF to a list.  

                        
                        
                    return ORFList1

        #prints the output of the program in the correct format, printing each tuple returned on a separate line
        def printOutput(rf):
            #f = open("orf.txt", "a") #this can be uncommented to add to a file in the same directory
            if rf != []:
                for tpl in rf:
                    print(tpl)
                    #f.write(str(tpl))
                    print("\n")
        
        

        printOutput(orf(CodonList,1))
        printOutput(orf(CodonList2,2))
        printOutput(orf(CodonList3,3))
        printOutput(orf(revCompList, -1))
        printOutput(orf(revCompList2, -2))
        printOutput(orf(revCompList3, -3))
        
                    

    
        # closing the file
        Fp1.close();

    except IOError:
        print("error unable to read file or file does not exist!!!")
        print("Exiting the program")
        Fp1.close()
        sys.exit(1)


#****************  executing the program ***************

main()
              