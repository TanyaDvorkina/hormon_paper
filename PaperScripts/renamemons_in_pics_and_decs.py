#!/usr/bin/env python3
import math

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO

import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from networkx.drawing.nx_agraph import read_dot

import os

from subprocess import check_call

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


def draw_monomer_graph(G, outdir):
    write_dot(G, os.path.join(outdir, "graph_paper.dot"))
    try:
        check_call(['dot', '-Tpng', os.path.join(outdir, "graph_paper.dot"), '-o', os.path.join(outdir, "graph_paper.png")])
    except Exception:
        return

def read_monomer_graph(filename):
    G = nx.DiGraph(read_dot(filename))
    return G

def fix_nodes_names(G, nodes_mp):
    nx.relabel_nodes(G, mapping)
    return G

# def get_manual_mapping(cen):
#     filename = "../CA_monomers/ManualNamesMapping/cen" + cen +"/cen" + cen + "_manualmap.tsv"
#     res = {}
#     with open(filename, "r") as fin:
#         for ln in fin.readlines():
#             m_ca, m_manual = ln.split("\t")[:2]
#             res[m_ca] = m_manual
#     return res

def get_manual_mapping(cen, mon_lst):
    manual_mons = load_fasta("../MonomersT2T/cen" + cen + "Monomers.fa")
    res = {}
    for m in mon_lst:
        best_m, best_ed = None, 1000500
        for mm in manual_mons:
            ed = edist([m.seq, mm.seq])
            if ed < best_ed:
                best_m, best_ed = mm.id, ed
        res[m.id] = best_m
    return res

def get_ca_mapping(cen, mon_lst):
    ca_mons = load_fasta("../CA_monomers/ManualNamesMapping/cen" + cen +"/cen" + cen + "_monomers.fa")
    res = {}
    for m in mon_lst:
        best_m, best_ed = None, 1000500
        for mm in ca_mons:
            ed = edist([m.seq, mm.seq])
            if ed < best_ed:
                best_m, best_ed = mm.id, ed
        res[m.id] = best_m
    return res

def build_new_names_from_final_mons(mon_lst, labels_mp, cen):
    res = {}
    manual_mp = get_manual_mapping(cen, mon_lst)
    print(manual_mp)
    for m in mon_lst:
        print(m.description)
        if m.id in labels_mp:
            res[m.id] = m.description.split(" ")[1]+ "(" + manual_mp[m.id] + ")[" + labels_mp[m.id].split('[')[1][:-1]
    return res

def build_new_names_from_ca_mons(mon_lst, labels_mp, cen):
    res = {}
    manual_mp = get_manual_mapping(cen, mon_lst)
    for m in mon_lst:
        print(m.description)
        if m.id in labels_mp:
            res[m.id] = m.id+ "(" + manual_mp[m.id] + ")[" + labels_mp[m.id].split('[')[1]
    return res

def update_dec(filename, mapping):
    with open(filename, "r") as fin:
        with open(filename[:-len(".tsv")] + "_paper.tsv", "w") as fout:
            for ln in fin.readlines():
                mon = ln.split("\t")[1]
                rev = ""
                if mon.endswith("'"):
                    rev = "'"
                    mon = mon[:-1]
                new_ln = "\t".join([ln.split("\t")[0], mapping[mon].split("(")[0] + rev] + ln.split("\t")[2:])
                fout.write(new_ln)

def convert2old(mon_str, mapping):
    res = []
    for m in mon_str.split(","):
        res.append(mapping[m].split("(")[0])
    return ",".join(res)

def update_chor(filename, mapping):
    with open(filename, "r") as fin:
        with open(filename[:-len(".tsv")] + "_paper.tsv", "w") as fout:
            for ln in fin.readlines():
                lst = ln.strip().split("\t")
                fout.write(lst[0] + "\t" + convert2old(lst[1], mapping) + "\n")
                print(lst[0] + "\t" + convert2old(lst[1], mapping) + "\n")

def load_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    for i in range(len(records)):
        records[i] = records[i].upper()
    return records

cens = [str(i) for i in range(1, 23)] + ["X"]
for cen in cens:
    if os.path.exists(os.path.join("./cen" + cen, "graph.dot")):
        mon_lst = load_fasta(os.path.join("./cen" + cen, "mn.fa"))
        G = read_monomer_graph(os.path.join("./cen" + cen, "graph.dot"))
        labels = G.nodes.data("label")
        labels_mp = {}
        for l in labels:
            labels_mp[l[0]] = l[1]
        mapping = build_new_names_from_final_mons(mon_lst, labels_mp, cen)
        print(labels_mp)
        print(mapping)
        nx.set_node_attributes(G, mapping, name="label")

        draw_monomer_graph(G, "./cen" + cen)
        update_dec(os.path.join("./cen" + cen, "final_decomposition.tsv"), mapping)
        update_chor(os.path.join("./cen" + cen, "HORs.tsv"), mapping)

        G = read_monomer_graph(os.path.join("./cen" + cen, "final_simplified_graph", "graph.dot"))
        nx.set_node_attributes(G, mapping, name="label")
        draw_monomer_graph(G, os.path.join("./cen" + cen, "final_simplified_graph") )

        mon_lst = load_fasta(os.path.join("./cen" + cen, "merge_split", "mn.fa"))
        G = read_monomer_graph(os.path.join("./cen" + cen, "simplified_graph", "graph.dot"))
        labels = G.nodes.data("label")
        labels_mp = {}
        for l in labels:
            labels_mp[l[0]] = l[1]
        mapping = build_new_names_from_ca_mons(mon_lst, labels_mp, cen)
        print(labels_mp)
        print(mapping)
        nx.set_node_attributes(G, mapping, name="label")
        draw_monomer_graph(G, os.path.join("./cen" + cen, "simplified_graph"))

        G = read_monomer_graph(os.path.join("./cen" + cen, "merge_split", "graph.dot"))
        nx.set_node_attributes(G, mapping, name="label")
        draw_monomer_graph(G, os.path.join("./cen" + cen, "merge_split") )

    else:
        print("No monomers", cen)







