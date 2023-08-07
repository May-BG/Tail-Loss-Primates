# import packages
import sys
import csv
from Bio import AlignIO

# read input fasta file and name output csv file 
name=str(sys.argv[1]) + "_new_8species.fasta"
output=str(sys.argv[1]) + "_8species_03.csv"
alignment = AlignIO.read(name,"fasta")

# variants shared by six hominoid species ("hg38", "gorGor5", "panTro5", "panPan2", "ponAbe2", "nomLeu3"), while different in two non-hominoid species("macNem1","calJac3")
ninser=0
y=0
ngap=0
a=[]
for r in range(0,len(alignment[0].seq)):
    if alignment[0,r] == "-":
        ngap=ngap+1
        if alignment[6,r] not in ["-","N","*"] and alignment[7,r] not in ["-","N","*"] and alignment[1,r] in ["-","N","*"] and alignment[2,r] in ["-","N","*"] and alignment[3,r] in ["-","N","*"] and alignment[4,r] in ["-","N","*"] and alignment[5,r] in ["-","N","*"]:
            a.append([r-ngap, alignment[0,r],alignment[5,r], alignment[6,r], alignment[7,r], y, ngap,ninser])
    else:
        ngap=ngap
        if alignment[0,r] == alignment[5,r] and alignment[1,r] == alignment[5,r] and alignment[2,r] == alignment[5,r] and alignment[3,r] == alignment[5,r] and alignment[4,r] == alignment[5,r] and alignment[5,r] != alignment[6,r] and alignment[5,r] != alignment[7,r] and alignment[0,r] not in ["-","N","*"] and alignment[6,r] not in ["-","N","*"] and alignment[7,r] not in ["-","N","*"]:
            y=y+1
            a.append([r-ngap, alignment[0,r],alignment[5,r], alignment[6,r],alignment[7,r], y, ngap,ninser])
        else:
            if alignment[6,r] in ["-","N","*"] and alignment[7,r] in ["-","N","*"] and alignment[0,r] not in ["-","N","*"] and alignment[1,r] not in ["-","N","*"] and alignment[2,r] not in ["-","N","*"] and alignment[3,r] not in ["-","N","*"] and alignment[4,r] not in ["-","N","*"] and alignment[5,r] not in ["-","N","*"]:
                ninser=ninser+1
                a.append([r-ngap, alignment[0,r],alignment[5,r], alignment[6,r],alignment[7,r], y, ngap,ninser])
with open(output, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(a)
