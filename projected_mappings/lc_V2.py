import os
import sys

slice15ls = os.listdir("./slice_15")
slice14ls = os.listdir("./slice_14")



fnameIntersect = set(slice15ls) & set(slice14ls)

fnameFile = open ("fnames", "w")
for x in fnameIntersect:
	fnameFile.writelines(x + "\n")