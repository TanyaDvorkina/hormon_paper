import os

def load_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    for i in range(len(records)):
        records[i] = records[i].upper()
    return records

def get_manual_mapping(cen):
    filename = "../CA_monomers/ManualNamesMapping/cen" + cen +"/cen" + cen + "_manualmap.tsv"
    res = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            m_ca, m_manual = ln.split("\t")[:2]
            res[m_ca] = m_manual
    return res

def remove_cen(m):
    new_m = list(m)[:-1]
    while new_m[-1].isdigit() or new_m[-1] == "/":
        new_m = new_m[:-1]
    return "".join(new_m)

res = []
cens = [str(i) for i in range(1, 23)] + ["X"]
for cen in cens:
    if not os.path.exists(os.path.join("./cen" + cen, "MonoRunRaw", "L.csv")):
        continue
    mp = get_manual_mapping(cen)
    with open(os.path.join("./cen" + cen, "MonoRunRaw", "L.csv"), "r") as fin:
        isfirst = True
        for ln in fin.readlines():
            cycle_name, mons, weight = ln.strip().split()
            mons = mons.replace("(", "").split(")")[:-1]
            new_mons, manual_mons = [], []
            hor_name = ""
            for m in mons:
                new_m = []
                for mm in m.split("+"):
                    if len(mm.split(".")) > 1:
                        new_m.append(remove_cen(mm.split(".")[0]) + "." + mm.split(".")[1])
                    else:
                        new_m.append(remove_cen(mm.split(".")[0]))
                    if len(new_m[-1]) > 1 and len(m.split("+")) == 1:
                        new_m[-1] = "(" + new_m[-1] + ")"
                if len(new_m) > 1:
                    new_mons.append("(" + "+".join(new_m) + ")")
                else:
                    new_mons.append(new_m[0])
                hor_name = mp[m.split("+")[0].split(".")[0]].split('.')[0]
                manual_mons.append(mp[m.split("+")[0].split(".")[0]][len(hor_name) + 1:])

            if isfirst:
                res.append("\t".join([cen, cycle_name, "".join(new_mons), hor_name + "." + ",".join(manual_mons), weight]))
            else:
                res.append("\t".join(["", cycle_name, "".join(new_mons), hor_name + "." + ",".join(manual_mons), weight]))
            isfirst = False

print("\n".join(res))


