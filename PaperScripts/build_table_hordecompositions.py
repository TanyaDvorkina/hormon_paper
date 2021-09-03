
MONOIDNT = 95


cens = [str(i) for i in range(1, 23)] + ["X"]
lens = [6,4,17,19,6,18,6,11,7,8,5,8,11,8,11,10,16,12,2,16,11,8,12]
canonical_len = {}

for i in range(len(cens)):
    canonical_len[cens[i]] = {"H1": lens[i]}
canonical_len["17"]["H2"] = 14

def load_monodec(filename):
    dec = []
    monomers = set()
    rc_num = 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
           if len(ln.strip().split("\t")) < 5:
              continue
           ref, mon, start, end, idnt = ln.strip().split("\t")[:5]
           if float(idnt) > MONOIDNT:
               dec.append([ref, mon, start, end, idnt])
               monomers.add(mon)
               if mon.endswith("'"):
                   rc_num += 1
    if rc_num > 0.5*len(dec):
         revert = True
    else:
         revert = False
    return dec, monomers, revert

def load_hordec(filename):
    dec = []
    monomers = set()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.strip().split("\t")) < 5:
                continue
            ref, hor, start, end, idnt, hor_name = ln.strip().split("\t")[:6]
            add = ""
            if hor.endswith("'"):
                add = "'"
            hor_ini = hor.split("<sup>")[0]
            if len(hor.split("<sup>")) > 1:
                hor_run = hor.split("<sup>")[1].split("<")[0]
            else:
                hor_run = 1
            dec.append([ref, hor_ini, int(hor_run), start, end, idnt, hor_name])
    return dec

def get_len(h, cen, hor_name):
    if h.startswith("c"):
        return canonical_len[cen][hor_name]
    elif h.startswith("p"):
        can_len = canonical_len[cen][hor_name]
        if "<sub>-" in h:
            l, r = h.split("<sub>-")[1].split("--")
            r = r.split("</sub>")[0]
            h = "p" + "<sub>" + r + "-" + l + "</sub>"
        l, r = map(int, h.split("<sub>")[1].split("<")[0].split("-"))
        if l < r:
            return r - l + 1
        else:
            return can_len - l + 1 + r
    else:
        return 1

def calc_stats(hor_dec, cen):
    hors = {}
    for h in hor_dec:
        if h[1].startswith("c"):
            if h[1].split("<sub>")[0] not in hors:
                hors[h[1].split("<sub>")[0]] = {"runs": 0, "cnt": 0, "len": get_len(h[1], cen, h[6])}
            hors[h[1].split("<sub>")[0]]["runs"] += 1
            hors[h[1].split("<sub>")[0]]["cnt"] += h[2]
        else:
            if h[1] not in hors:
                hors[h[1]] = {"runs": 0, "cnt": 0, "len": get_len(h[1], cen, h[6])}
            hors[h[1]]["runs"] += 1
            hors[h[1]]["cnt"] += h[2]
    return hors


res = []
for cen in cens:
    print(cen)
    mono_dec, mons, revert = load_monodec("./MonomerInference/cen" + cen + "/final_decomposition_paper.tsv")
    hor_dec = load_hordec("./HORInference/cen" + cen + "_reverted_collapsed.tsv")
    hors_stat = calc_stats(hor_dec, cen)
    can_nums = []
    for i in range(2):
        can = "c" + str(i)
        if can in hors_stat:
            can_nums.append(str(hors_stat[can]["cnt"]) + "/" + str(hors_stat[can]["runs"]))
    mon_total, mon_partial = 0, 0
    for h in hors_stat:
        mon_total += hors_stat[h]["len"]*hors_stat[h]["cnt"]
        if not h.startswith("c"):
            mon_partial += hors_stat[h]["len"]*hors_stat[h]["cnt"]
    hors_lst = sorted([[h.replace("<sub>", "").replace("</sub>", ""), hors_stat[h]["cnt"]] for h in hors_stat if not h.startswith("c")], key = lambda x: -x[1])
    res.append("\t".join([cen, str(len(mons))+ "/" + str(len(mono_dec)), " ".join(can_nums),\
           str(mon_partial) + "/" + str(int(mon_partial*100/mon_total)),\
           " ".join(["/".join(map(str, x)) for x in hors_lst[:3]]), str(len(hor_dec))]))
print("\n".join(res))
