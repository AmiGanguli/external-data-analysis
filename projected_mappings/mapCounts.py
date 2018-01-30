import sys


def mapCounts():
	fnames = [line.rstrip('\n') for line in open(sys.argv[2])]
	path_to_data = sys.argv[3]
	with open(str((path_to_data.split('/')[-2:])[0]) + "_geneCounts.tab", "w") as oFile:
		oFile.write("file_name\t#os_genes\t#total_maps\tavg_maps/gene")
		for fileName in fnames:
			mappedGeneCount = 0
			fileLength = 0
			with open(path_to_data + fileName, "r") as file:
				for line in file:
					mappedGenes=(line.split('\t')[1]).split(' ')
					mappedGeneCount = mappedGeneCount + len(mappedGenes)
					fileLength = fileLength + 1
				oFile.write(fileName + "\t" + str(fileLength) + "\t" + str(mappedGeneCount) + "\t" + str(mappedGeneCount/fileLength) + "\n")

def mapCountDiff():
	path_to_data1 = sys.argv[2]
	path_to_data2 = sys.argv[3]
	oFile = open("mapCountDiff.tab", "w")
	oFile.writelines("file_name\t\u0394os_genes\t\u0394total_maps\t\u0394avg_maps/gene\n")
	with open(path_to_data1, "r") as iFile1, open(path_to_data2, "r") as iFile2:
		iFile1.readline()
		iFile2.readline()
		print ("thing")
		for line1, line2 in zip(iFile1, iFile2):
			line1Vals = line1.split('\t')
			line2Vals = line2.split('\t')
			name = line1Vals[0]
			osDiff = int(line2Vals[1]) - int(line1Vals[1])
			mapsDiff = int(line2Vals[2]) - int(line1Vals[2])
			avgDiff = float(line2Vals[3].rstrip('\n')) - float(line1Vals[3].rstrip('\n'))
			oFile.write(name + "\t" + str(osDiff) + "\t" + str(mapsDiff) + '\t' + str(('+' if (avgDiff > 0) else '') + str(avgDiff)) + '\n')
			print (name + "\t" + str(osDiff) + "\t" + str(mapsDiff) + '\t' + str(avgDiff))
				
def main():
	command = sys.argv[1]
	if (command == "mapCounts"):
		mapCounts()
	elif (command == "mapCountDiff"):
		mapCountDiff()
main()