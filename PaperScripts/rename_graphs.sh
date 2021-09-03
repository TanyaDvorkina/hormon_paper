#!/usr/bin/env bash

array2=(./cen*)
for f in "${array2[@]}"
do
    echo "ls "${f}"/graph_paper.png"
    echo "$(basename ${f})"
    cp "${f}"/graph_paper.png "${f}"/"$(basename ${f})"_graph_paper.png
    cp "${f}"/merge_split/graph_paper.png "${f}"/"$(basename ${f})"_mergesplit_graph_paper.png
    cp "${f}"/simplified_graph/graph_paper.png "${f}"/"$(basename ${f})"_simplified_graph_paper.png
    cp "${f}"/final_simplified_graph/graph_paper.png "${f}"/"$(basename ${f})"_final_simplified_graph_paper.png
    cp "${f}"/MonoRunRaw/monorun.png "${f}"/"$(basename ${f})"_monorunraw.png
    cp "${f}"/MonoRunRaw/monorun_splv.png "${f}"/"$(basename ${f})"_monorunraw_splv.png
    cp "${f}"/MonoRun/monorun.png "${f}"/"$(basename ${f})"_monorun.png
    cp "${f}"/MonoRun/monorun_splv.png "${f}"/"$(basename ${f})"_monorun_splv.png
    cp "${f}"/MonoRun/IterMnrun.png "${f}"/"$(basename ${f})"_IterMnrun.png
    cp "${f}"/MonoRun/SplitIterMnrun.png "${f}"/"$(basename ${f})"_SplitIterMnrun.png
    cp "${f}"/MonoRunRaw/IterMnrun.png "${f}"/"$(basename ${f})"_IterMnrunRaw.png
    cp "${f}"/MonoRunRaw/SplitIterMnrun.png "${f}"/"$(basename ${f})"_SplitIterMnrunRaw.png
done
