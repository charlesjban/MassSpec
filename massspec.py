
import argparse
import sys
from pathlib import Path


## dictionaries for the monoisotopic and average isotopic masses
monoMassesDict = {'A' :  71.0371, 'C' : 103.0092, 'D' : 115.0269,
'E' : 129.0426, 'F' : 147.0684, 'G' :  57.0215, 
'H' : 137.0589, 'I' : 113.0841,'K' : 128.0950, 
'L' : 113.0841, 'M' : 131.0405, 'N' : 114.0429,
'P' :  97.0528, 'Q' : 128.0586, 'R' : 156.1011, 
'S' :  87.0320, 'T' : 101.0477, 'V' :  99.0684, 
'W' : 186.0793, 'Y' : 163.0633, '*' : 0.0, 
'H2O' : 18.0106, 'proton' : 1, 'PO3': 79.9663}

averageMassesDict = {'A' :  71.08, 'C' : 103.14, 'D' : 115.09, 
'E' : 129.12,'F' : 147.18, 'G' :  57.05, 
'H' : 137.14, 'I' : 113.16,'K' : 128.17, 
'L' : 113.16, 'M' : 131.19, 'N' : 114.10,
'P' :  97.12, 'Q' : 128.13, 'R' : 156.19, 
'S' :  87.08,'T' : 101.10, 'V' :  99.13, 
'W' : 186.21, 'Y' : 163.18, '*' : 0.0, 
'X' : 0, 'H2O' : 18.0153, 'proton' : 1, 'PO3': 79.98}


## create and parse arguments for command line use
options = argparse.ArgumentParser(description='Determine isotopic masses for each digested peptide')
options.add_argument("-f", help="choose .fasta file with which to  perform task 3, mass spec analysis")
options.add_argument("-i", help="choose either 'monoisotopic' or 'averageisotopic' mass values ['m','a']", default='a')
options.add_argument("-c", help="choose a charge for peptides[1,2,3]", default='1')
options.add_argument("-t", help="choose to report back only N or C terminal peptides ['n','c','a' ]", default='a')
options.add_argument("-s", help="choose to have basic stats files output. y/n. ['y','n']", default='n')
options.add_argument("-p", help="choose to modify Ser, Try and Thr by adding a phosphorus group. y/n. ['y','n']", default='n') 
args = options.parse_args()


### Argument error catching - ensures all arguments are correct, and assign to corresponsing variables if required 
if (args.s is not 'y') and (args.s is not 'n'):  
	sys.exit("Error: Stats argument ('-s') must be 'y' or 'n'. Default is 'n'.")

if (args.i is not 'm') and (args.i is not 'a'): 
	sys.exit("Error: Isotopic mass (-i) argument must be 'm' (monoisotopic masses) or 'a' (average isotopic). Default is 'a'.")  # halts the program
if args.i == 'm':
	massDictionary = monoMassesDict
else:
	massDictionary = averageMassesDict 

if (args.c is not '1') and (args.c is not '2') and (args.c is not '3'): 
	sys.exit("Error: Charge (-c) can only take value of either 1, 2, or 3. Default is '1'.")
if args.c == '3':
	charge = 3
elif args.c == '2':
	charge = 2
else:
	charge = 1

if (args.t is not 'n') and (args.t is not 'c') and (args.t is not 'a'):
	sys.exit("Error: terminal (-t) argument must be 'n' (n-terminal), 'c' (c-terminal) or 'a' (all peptides). Default is 'a'.")
if args.t == 'n': 
	terminal = 'n'
elif args.t == 'c':
	terminal = 'c'
else:
	terminal = 'a'

if args.f == None: 
	sys.exit("Error: Please choose fasta file (using -f tag) to perform analysis eg. 'massspec.py -f digested_a.fasta'")
file = args.f 
if '.fasta' not in file[-6:]: 
	sys.exit("Error: please make sure file is named with '.fasta' suffix")

if (args.p is not 'y') and (args.p is not 'n'):  #error checking if the argument is not 'y' or 'n'
	sys.exit("Error: phosphorus '-p' argument must must be 'y' for yes or 'n' for no")  # halts the program

## check that the file can be found at location
fileCheck = Path(file) 
if not fileCheck.is_file(): 
	sys.exit("Error: Please check chosen file (-f) exists and is in working directory.")  # halts the program


fileObj = open(file, 'r')
outFile = open(file[:-6]+'.masses', 'w') 


# dictionary to store N- or C-terminal peptides
peptideDictionary = {} 
# terminal dictionary stores 
terminalDictionary = {} 

## read each line, treating as groups of 2
lines = fileObj.readlines() 
for index in range(0, len(lines),2):
	## create variable for each key component 
	heading = lines[index] 
	splitHeading = heading.split() 
	peptideName = splitHeading[0][1:]  
	peptideNumber = splitHeading[1]
	missedCleaves = splitHeading[2][7]
	enzyme = splitHeading[3]
	## next line is amino acid sequence
	nextLine = lines[index + 1] 
	aaSeq = nextLine.replace("\n","")

	## fetch the values for each amino acid from the dictionary 
	residueValue = [] 
	for count in aaSeq:
		residueValue.append(massDictionary.get(count, 0))
		if massDictionary.get(count, 0) == 0:
			print("Warning: unknown amino acid: "+count+"in "+ peptideName +". M/Z value for this AA: 0")
		## additional value for phosphate if this argument is chosen
		if args.p is 'y':
			if count is "S":
					residueValue.append(massDictionary['PO3'])
			elif count is "Y":
				residueValue.append(massDictionary['PO3'])
			elif count is "T":
				residueValue.append(massDictionary['PO3'])
		residueValueList = (residueValue) 
	
	## sum to get total  weight and then divide by the charge to get M/Z ratio
	peptideValue = sum(residueValueList)
	peptideValueFull = (float(peptideValue) + charge + massDictionary['H2O'])/charge 
	peptideValue4sf = format(peptideValueFull, '.4f') 
	
	## function tp print in correct format to the output file
	def outputPrint(): 
		print(peptideName.ljust(20), peptideNumber.rjust(2), 
			str(peptideValue4sf).rjust(10), missedCleaves.rjust(1), 
			repr(charge).rjust(1), enzyme.rjust(1), aaSeq, file=outFile) 
	
	## terminal == 'n' - print only the first instance of a peptide name, which will be the n-terminal peptide
	if terminal =='n':  
		if peptideName not in terminalDictionary:
			terminalDictionary[peptideName] = peptideName 
			outputPrint()
	## each instance of peptide name will replace the previous seq, so only the c-terminal is saved
	elif terminal == 'c':
		terminalDictionary[peptideName] = {'number': peptideNumber, 'mass': peptideValue4sf, 'missed': missedCleaves, 'charge': charge, 'enzyme' :enzyme, 'sequence': aaSeq}
	else:
		outputPrint()
	## counts the intances of each peptide
	if peptideName not in peptideDictionary: # if the peptide is not already in dictionary
		peptideDictionary[peptideName] =  1 # creates a key of peptide and assigns it a value of 1
	else: 									# else, if the peptide has been encountered already, add 1
		peptideDictionary[peptideName] += 1

# print c-terminal peptides from the dictionary to output format
if terminal == 'c':  
	for keys in terminalDictionary:  
		print(keys.ljust(20), terminalDictionary[keys]['number'].rjust(2), 
			terminalDictionary[keys]['mass'].rjust(10), 
			terminalDictionary[keys]['missed'].rjust(1), 
			repr(terminalDictionary[keys]['charge']).rjust(1),
			terminalDictionary[keys]['enzyme'].rjust(1),
			terminalDictionary[keys]['sequence'], file=outFile)

fileObj.close() #close files
outFile.close()

## additional stats output if requested
if args.s == 'y':
	statsFile = open(file[:-6]+'.csv', 'w')
	statsOverview = open(file[:-6]+'.stats', 'w')
	peptideNumberList= [] 
	totalProteins = len(peptideDictionary) 
	print('There are',totalProteins, 'digested proteins which cleave to make peptides in the required range', file=statsOverview) #prints the number of proteins
	for keys in peptideDictionary: 
		peptideNumberList.append(peptideDictionary[keys])  
		print(keys, ',',peptideDictionary[keys], file=statsFile)
	totalPeptides = sum(peptideNumberList) 
	print('There are',totalPeptides,'of these peptides in total', file=statsOverview)  
	averagePeptides = format((totalPeptides / totalProteins),'.4f')  
	print('The average number of peptides (in given range) per ("useful") protein for',file ,'is', averagePeptides, file=statsOverview)  #prints the average
