import re
import datetime
import math

def writeRowToFile(row, outFile):
	InteractionTypeMIs = ["psi-mi:\"MI:2232\"(molecular interaction)", "psi-mi:\"MI:2235\"(up-regulate)", "psi-mi:\"MI:2240\"(down-regulate)"]
	InteractorTypeMIs = ["", "", "psi-mi:\"MI:0326\"(protein)"]
	EvidenceMIs = ["psi-mi:\"MI:2231\"(coexpression)", "psi-mi:\"MI:0037\"(domain profile pairs)", "", "psi-mi:\"MI:1055\"(internally-curated)", "psi-mi:\"MI:0085\"(phylogenetic profile)"]
	line = ""
	for x, item in enumerate(row):
		value = item
		if str(value) == "None" or str(value) == "NULL" or value == None or value == '':
			value = "-"
		elif x == 0 or x == 1:
			taxVal = row[x + 9]
			if taxVal == "4530":
				value = "msu:" + value
			elif taxVal == "3702":
				value = "agi:" + value
			elif taxVal == "3964":
				value = "jgi:" + value
			elif taxVal == "39947":
				value = "stringdb:" + value
		elif x == 6:
			value = EvidenceMIs[value]
		elif x == 8:
			value = "pubmed:" + str(value)
		elif x == 9 or x == 10:
			value = "taxid:" + str(value)
		elif x == 11:
			value = InteractionTypeMIs[value]
		elif x == 13:
			value = "jaiswal:" + str(value)
		elif x == 14:
			if str(value) == "-":
				value = "intact-miscore:0.0"
			else:
				value = "intact-miscore:" + str(value)
		line += (str(value) + "\t")
	line = line + "\n"
	outFile.write(line)

def main():
	counts = [0 for x in range(5)]
	combinations = [[0 for y in range(5)] for x in range(5)]
	dataFile = open("LLSmod.tab", "r")
	rowCount = 0
	outFile = open("data/psimi_" + datetime.datetime.now().strftime("%y-%m-%d-%H-%M") + '.tab' , "w")
	field_names = ["ID(s) interactor A", "ID(s) interactor B", "Alt. ID(s) interactor A", "Alt. ID(s) interactor B", "Alias(es) interactor A", "Alias(es) interactor B", "Interaction detection method(s)", "Publication 1st author(s)",
		"Publication Identifier(s)", "Taxid Interactor A", "Taxid Interactor B", "Interaction Types", "Source Databases", "Interaction Identifier(s)", "Confidence value(s)"]
	outFile.write("\t".join(field_names) + "\n")
	outputRow = ['' for x in range(15)]
	for fileRow in dataFile:
		row = re.split(r'\t+', fileRow)[:8]
		for y in range(3, 8):
			if (y != 6 and row[y] != "NA"):
				include = True
				for z in range(3, 8):
					if (z != y and row[z] != "NA"):
						include = False
						break
				if (include == True):
					outputRow[0] = str("AGI:" + row[1])
					outputRow[1] = str("AGI:" + row[2])
					outputRow[6] = y - 3
					outputRow[7] = 'Lee'
					outputRow[8] = 20118918
					outputRow[9] = 3702
					outputRow[10] = 3702
					outputRow[11] = 0
					outputRow[12] = 'jaiswal'
					outputRow[13] = row[0]
					confScore = round(math.log10(float(row[y])), 2)
					outputRow[14] = str(confScore)
					writeRowToFile(outputRow, outFile)
					rowCount = rowCount + 1
					# print(outputRow)
					if (rowCount % 1000 == 0):
						print(outputRow)
		
main()

