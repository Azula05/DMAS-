#!/usr/bin/env python3

# temperary script to test tm prediction

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq


myseq = Seq('TCTTTACTCACCTGTAGATgtaTC')
template = Seq('AGAAATGAGTGGACATCTACAAAG')
template_correct= Seq("AGAAATGAGTGGACATCTACATAG")

print("match")
match = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
print('%0.2f' % match)
print("mismatch")
mismatch = mt.Tm_NN(myseq, c_seq=template , nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
print('%0.2f' % mismatch)
#mt.Tm_NN Arguments(nn_table: 4, dnac1: 250.00, dnac2: 250.00, Na: 50.00, K: 0.00, Tris: 75.00, Mg: 3, dNTPs: 1.2, saltcorr: 7)
print("Delta Tm")
print('%0.2f' % (match - mismatch))