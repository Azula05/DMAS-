#!/usr/bin/env python3

# temperary script to test tm prediction

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

dnac = 250
Na_conc = 50
K_conc = 0
Tris_conc = 75
Mg_conc = 3
dNTPs_conc = 1.2


myseq = Seq('TCTCTGCTCTGCAGCAactTC')
template = Seq('TCTCTGCTCTGCAGAAactTC')
template = template.complement()
print(template)

print("match")
match = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
print('%0.2f' % match)
print("mismatch")
mismatch = mt.Tm_NN(myseq, c_seq=template , nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
print('%0.2f' % mismatch)
#mt.Tm_NN Arguments(nn_table: 4, dnac1: 250.00, dnac2: 250.00, Na: 50.00, K: 0.00, Tris: 75.00, Mg: 3, dNTPs: 1.2, saltcorr: 7)
print("Delta Tm")
print('%0.2f' % (match - mismatch))