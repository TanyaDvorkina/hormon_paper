from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

import edlib

def edist(lst):
    if len(str(lst[0])) == 0:
        return 100500
    if len(str(lst[1])) == 0:
        return 100500
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW")
    if result["editDistance"] == -1:
        return 100500
    return result["editDistance"]

def load_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    for i in range(len(records)):
        records[i] = records[i].upper()
    return records


alphasat_seq = load_fasta("./AlphaSat.fa")[0]
mons = load_fasta("./monomers.fa")

new_mons = []
reverted = 0
for m in mons:
    ed, ed_rev = edist([m.seq, alphasat_seq.seq]), edist([m.seq.reverse_complement(), alphasat_seq.seq])
    if ed > ed_rev :
        print("Revert: ", m.id, ed, ed_rev)
        reverted += 1
        new_record = SeqRecord(m.seq.reverse_complement(), id=m.id, name=m.name, description = m.description)
        new_mons.append(new_record)
    else:
        new_mons.append(m)

print("Total:", len(new_mons), "Reverted: ", reverted)

with open("./onestrand_monomers.fa", "w") as fw:
        SeqIO.write(new_mons, fw, "fasta")
print("Saved to ./converted_monomers.fa")

