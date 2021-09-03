import os
from shutil import copyfile

dirs = os.listdir("./FinalConsensus")

for d in dirs:
  if d.endswith("L"):
      if os.path.exists(os.path.join("./FinalConsensus", d, "monomer_consensus.fasta")):
          cen = "cen" + d.split("_")[1]
          if cen not in {"cen3", "cen4", "cen5", "cen8", "cen11"} \
              or d.startswith("hor_3_2") or d.startswith("hor_4_2") \
              or d.startswith("hor_5_5") or d.startswith("hor_8_2") or d.startswith("hor_11_2"):
              print("Copy",os.path.join("./FinalConsensus", d, "monomer_consensus.fasta"), "to",  os.path.join("./MonomersT2T", cen + "Monomers.fa"))
              copyfile(os.path.join("./FinalConsensus", d, "monomer_consensus.fasta"), \
                       os.path.join("./MonomersT2T", cen + "Monomers.fa" ))
              print("Copy",os.path.join("./FinalConsensus", d, "hor_consensus.fasta"), "to",  os.path.join("./HORsT2T", cen + "HOR.fa"))
              copyfile(os.path.join("./FinalConsensus", d, "hor_consensus.fasta"), \
                       os.path.join("./HORsT2T", cen + "HOR.fa" ))
