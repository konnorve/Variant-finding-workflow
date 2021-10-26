from Bio import pairwise2
from Bio.Seq import Seq
a = Seq("AACTGGC")
b = Seq("GAGGC")

alignment = pairwise2.align.globalms(a, b, 2, 0, -1, 0)[0]

aligned_ref = alignment.seqA
aligned_alt = alignment.seqB

mismatches = []
insertions = []
deletions = []

for i in range(len(aligned_ref)):
    if aligned_ref[i] != aligned_alt[i]:
        if aligned_ref[i] == "-":
            insertions.append(i)
        elif aligned_alt[i] == "-":
            deletions.append(i)
        else:
            mismatches.append(i)


i = 0
while i < len(insertions):
    start = insertions[i]
    while i + 1 < len(insertions):
        if insertions[i] + 1 == insertions[i+1]:
            i += 1
    end = insertions[i]
    i += 1
    print("insertion from {} to {}".format(start, end))

i = 0
while i < len(deletions):
    start = deletions[i]
    while i + 1 < len(deletions):
        if deletions[i] + 1 == deletions[i+1]:
            i += 1
    end = deletions[i]
    i += 1
    print("deletion from {} to {}".format(start, end))

import deep_variant_calling as dvc

print(dvc.determineEffect(a, b))

