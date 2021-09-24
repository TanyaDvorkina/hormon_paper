import os

res = []
cens = [str(i) for i in range(1, 23)] + ["X"]
for cen in cens:
    if not os.path.exists(os.path.join("./cen" + cen, "merge_split", "MergeSplitStat.csv")):
        continue
    with open(os.path.join("./cen" + cen, "merge_split", "MergeSplitStat.csv"), "r") as fin:
        ln1 = fin.readline().strip().split(",")
        ln2 = fin.readline().strip().split(",")
        mp = {}
        for i in range(len(ln1)):
            mp[ln1[i]] = ln2[i]
    mn, mn_plus = mp["#Init mon"], mp["#Final mon"]
    merge, split = mp["#Merge"], mp["#Split"]
    sqmean, dbindex = "{:.2f}".format(float(mp["Final SqMean"])), "{:.2f}".format(float(mp["Final DBIndex"]))

    with open(os.path.join("./cen" + cen, "MonoRunRaw", "Monorunscnt.csv"), "r") as fin:
        edges = fin.readline().strip().split()[-1]
    res.append("\t".join([cen, mn + "/" + mn_plus, edges, merge + "/" + split, sqmean, dbindex]))

print("\n".join(res))
