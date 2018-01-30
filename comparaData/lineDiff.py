lc14 = open("lineTotals14", "r")
lc15 = open("lineTotals15", "r")

orthCountFile = open("orthCountDiff", "w")

orthCountFile.writelines("Species_Name\tCompara_EG29-32\tCompara_EG38\tchange\n")

for line14, line15 in zip(lc14, lc15):
	row14 = (line14.split(" "))
	row15 = (line15.split(" "))
	lenR14 = len(row14)
	lenR15 = len(row15)
	
	orthCount14 = int(row14[lenR14-2])
	specName14 = ' '.join(row14[lenR14-1].split('_')[:2])
	orthCount15 = int(row15[lenR15-2])
	specName15 = ' '.join(row15[lenR15-1].split('_')[:2])
	orthCountDiff = orthCount15 - orthCount14
	orthCountFile.write(specName14 + "\t" + str(orthCount14) + "\t" + str(orthCount15) + "\t" + ('+' if orthCountDiff > 0 else '') + str(orthCountDiff) + "\n")
orthCountFile.close()