import re

def main():
	counts = [0 for x in range(5)]
	dataFile = open("LLS.tab", "r")
	outFile = open("LLSmod.tab", "w")
	rowNum = 0
	for row in dataFile:
		row = re.split(r'\t+', row)[:7]
		rowNum = rowNum + 1
		for x in range(5):
			if (row[x+2] != "NA"):
				row.insert(0, str(rowNum))
				modRow = '\t'.join(row)
				outFile.write(modRow + "\t\n")
				break
		
main()