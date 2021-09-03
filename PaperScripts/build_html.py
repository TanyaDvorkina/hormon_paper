import os

cens = ["X"] #[str(i) for i in range(1, 23)] + ["X"]
for cen in cens:
    if os.path.exists(os.path.join("./cen" + cen + "_reverted_collapsed.tsv")):
        print("<h3>cen" + cen + "</h3>")
        hors = {}
        best_hors = {}
        dec = []
        rc_num = 0
        with open(os.path.join("./cen" + cen + "_reverted_collapsed.tsv"), "r") as fin:
            for ln in fin.readlines():
                hor_collapsed = ln.split("\t")[1]
                hor_collapsed1 = hor_collapsed.split("<sup>")[0]
                if hor_collapsed1 not in hors:
                    hors[hor_collapsed1] = 0
                if "sup" in hor_collapsed:
                    hors[hor_collapsed1] += int(hor_collapsed.split("<sup>")[1][:-len("</sup>")])
                else:
                    hors[hor_collapsed1] += 1
                if hor_collapsed1.startswith("c"):
                    best_hors[hor_collapsed1] = "red"
                if hor_collapsed.startswith("c0"):
                    dec.append([hor_collapsed.replace("c0", "c"), hor_collapsed1])
                else:
                    dec.append([hor_collapsed, hor_collapsed1])
                if "'" in hor_collapsed:
                    rc_num += 1
        if rc_num > 0.5*len(dec):
            revert = True
        else:
            revert = False
        hors_lst = sorted([[h, hors[h]] for h in hors if not h.startswith("c")], key = lambda x: -x[1])
        best_hors[hors_lst[0][0]] = "blue"
        if len(hors_lst) > 1:
           best_hors[hors_lst[1][0]] = "green"
        if len(hors_lst) > 2:
           best_hors[hors_lst[2][0]] = "brown"
        final_str = "<p>"
        if revert:
            for d in dec[::-1]:
                if "'" in d[0]:
                    horname = d[0].replace("'", "")
                else:
                    if "</sub>" in d[0]:
                        horname = d[0].split("</sub>")[0] + "</sub>'" + d[0].split("</sub>")[1]
                    elif "<sup>" in d[0]:
                        horname =  d[0].split("<sup>")[0] + "'<sup>" + d[0].split("<sup>")[1]
                    else:
                        horname = d[0]+"'"
                if d[1] in best_hors:
                    final_str += '<span style="color:' +  best_hors[d[1]] + ';">' + horname + '</span>'
                else:
                    final_str += horname
        else:
            for d in dec:
                if d[1] in best_hors:
                    final_str += '<span style="color:' +  best_hors[d[1]] + ';">' + d[0] + '</span>'
                else:
                    final_str += d[0]
        final_str += "</p>"
        print(final_str)
        print("")
    else:
        print("No file", os.path.join("./cen" + cen + "_collapsed.tsv"))