import os
import sys

slice15ls = os.listdir("./slice_15")
slice14ls = os.listdir("./slice_14")

index = 0

fnames14 = [None] * len(slice14ls)
fnames15 = [None] * len(slice15ls)

for x in slice14ls:
	splitStr = x.split("_")
	splitLen = len(splitStr)
	fnames14[index] = splitStr[splitLen - 2] + "_" + splitStr[splitLen - 1]
	os.rename("./slice_14/" + slice14ls[index], "./slice_14/" + fnames14[index])
	index = index + 1
	
index = 0

for x in slice15ls:
	splitStr = x.split("_")
	splitLen = len(splitStr)
	fnames15[index] = splitStr[splitLen - 2] + "_" + splitStr[splitLen - 1]
	os.rename("./slice_15/" + slice15ls[index], "./slice_15/" + fnames15[index])
	index = index + 1

#print ("length of fnames14: " + str(len(fnames14)))
#print ("length of fnames15: " + str(len(fnames15)))

fnameIntersect = set(fnames14) & set(fnames15)
#print (fnameIntersect)
#print ("length of fnameIntersect: " + str(len(fnameIntersect)))

