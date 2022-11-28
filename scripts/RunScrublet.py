import scrublet as scr
import scipy.io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np

import sys

args=sys.argv

matrix=args[1]
savename=args[2]

print("Load Matrix!")
counts_matrix = scipy.io.mmread(matrix).T.tocsc()


print("Make Scrublet Object")
scrub = scr.Scrublet(counts_matrix)
print("Run Scrublet")
doub_scores, predicted_doublets = scrub.scrub_doublets()

print("Save!")
np.savetxt(savename,doub_scores)

