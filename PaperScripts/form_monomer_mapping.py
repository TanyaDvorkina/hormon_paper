from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

import sys
import os


import edlib

#print("Add -B HOR manually!")
#exit(-1)

def edist(lst):
    if len(str(lst[0])) == 0:
        return 100500
    if len(str(lst[1])) == 0:
        return 100500
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW")
    if result["editDistance"] == -1:
        return 100500
    return result["editDistance"]


def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            records[r] = records[r].upper()
    else:
        records = list(SeqIO.parse(filename, "fasta"))
        for i in range(len(records)):
            records[i] = records[i].upper()
    return records

def get_mono_name(num):
    if num < 26:
        return chr(ord("A") + num)
    else:
        return chr(ord("A") + num//26 - 1) + chr(ord("A") + num%26)

def get_most_similar(mon, lst):
    best_mono, best_mono_ed = None, 100500
    for m in lst:
        cur_ed = edist([mon.seq, m.seq])
        if best_mono_ed > cur_ed:
            best_mono, best_mono_ed = m.id, cur_ed
    return best_mono, best_mono_ed

cens = [str(i) for i in range(1, 23)] + ["X"]
cens_mono = {}
for cen in cens:
    cens_mono[cen] = set()

mono_cens = {}
with open("./final_decomposition.tsv", "r") as fin:
    for ln in fin.readlines():
        ref, mon, start, end, idnt = ln.strip().split("\t")[:5]
        cen = ref.split("_")[0][3:]
        if float(idnt) > 95:
            if mon.endswith("'"):
                mon = mon[:-1]
            cens_mono[cen].add(mon)
            if mon not in mono_cens:
                mono_cens[mon] = set()
            mono_cens[mon].add(cen)

hormon_mono = load_fasta("./onestrand_monomers.fa", "map")
new_naming_map = {}
MAX_ED = 100
for cen in cens:
    if not os.path.exists("./ManualNamesMapping/cen" + cen):
        os.makedirs("./ManualNamesMapping/cen" + cen)
    manual_mono = load_fasta("../MonomersT2T/cen" + cen + "Monomers.fa")
    mono_map = {}
    cen_mono_lst = []
    used = set()
    used_manual = set()
    cen_mono_lst_old = [hormon_mono[m] for m in sorted(cens_mono[cen])]
    min_ind = 1
    for m in cen_mono_lst_old:
        best_mono, best_ed = get_most_similar(m, manual_mono)
        for mm in manual_mono:
            if mm.id == best_mono:
                best_ca_mono, best_ca_ed = get_most_similar(mm, cen_mono_lst_old)
                break
        #print(m.id, best_mono, best_ed, best_ca_mono, best_ca_ed, used)
        if m.id == best_ca_mono and best_ed < MAX_ED:
            let = best_mono.split(".")[1]
            if let.isnumeric():
                if "-B" in best_mono:
                    add = 16
                else:
                    add = 0
                new_name = get_mono_name(int(let) - 1 + add)
                m_name = new_name + "/".join(sorted(list(mono_cens[m.id])))
                mono_map[m_name] = [best_mono, str(best_ed)]
                new_record = SeqRecord(m.seq, id=m_name, name=m_name, description = best_mono)
                cen_mono_lst.append(new_record)
                new_naming_map[m.id] = m_name
                used.add(m.id)
                used_manual.add(best_mono)
                print(" ", m.id, best_mono, m_name, best_ed)
                min_ind = max(min_ind, int(let) + add, len(cen_mono_lst))

    for m in cen_mono_lst_old:
        if m.id not in used:
            m_name = get_mono_name(min_ind) + "/".join(sorted(list(mono_cens[m.id])))
            best_mono, best_ed = get_most_similar(m, manual_mono)
            print(m_name, best_mono, best_ed, min_ind, len(cen_mono_lst))
            mono_map[m_name] = [best_mono, str(best_ed)]
            if best_mono not in used_manual:
                print("Not used: ", m_name, best_mono, best_ed)
            new_record = SeqRecord(m.seq, id=m_name, name=m_name, description = best_mono)
            cen_mono_lst.append(new_record)
            new_naming_map[m.id] = m_name
            min_ind += 1

    cen_mono_lst = sorted(cen_mono_lst, key = lambda x: x.id)
    with open("./ManualNamesMapping/cen" + cen + "/cen" + cen + "_monomers.fa" , "w") as fw:
        SeqIO.write(cen_mono_lst, fw, "fasta")
    print("Saved to ./ManualNamesMapping/cen" + cen + "/cen" + cen + "_monomers.fa")
    with open("./ManualNamesMapping/cen" + cen + "/cen" + cen + "_manualmap.tsv" , "w") as fw:
        for m in mono_map:
            fw.write("\t".join([m, mono_map[m][0], mono_map[m][1]]) + "\n")

difstrand_mons = load_fasta("monomers.fa", "map")
strand_mon = {}
for m in hormon_mono:
    if hormon_mono[m].seq != difstrand_mons[m].seq:
        strand_mon[m] = True
    else:
        strand_mon[m] = False

with open("./final_decomposition.tsv", "r") as fin:
    with open("./final_decomposition_paper.tsv", "w") as fout:
        for ln in fin.readlines():
            lst = ln.split("\t")
            mon, rev = lst[1], False
            if mon.endswith("'"):
                mon = mon[:-1]
                rev = True
            rev = not rev if strand_mon[mon] else rev
            rev_str = "'" if rev else ""
            new_ln = "\t".join([lst[0], new_naming_map[mon] + rev_str] + lst[2:])
            fout.write(new_ln)






