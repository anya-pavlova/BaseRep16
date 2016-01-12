#!/bin/bash

rm -fr feat_rar_4000
multiple_rarefactions.py -i otu_rarified/otu_table.biom -m 4000 -x 4000 -s 20 -n 1 -o feat_rar_4000/ --lineages_included
alpha_diversity.py -i feat_rar_4000/ -o alpha_rare_4000/ -t "/home/qiime/qiime_software/gg_13_5_otus/trees/97_otus.tree" -m chao1
collate_alpha.py -i alpha_rare_4000/ -o alpha_collated_4000/ 
mv feat_rar_4000/rarefaction_4000_0.biom feat_rar_4000/otu_table.biom
convert_biom.py -i feat_rar_4000/otu_table.biom -o feat_rar_4000/otu_table.txt -b --header_key taxonomy
summarize_taxa.py -i feat_rar_4000/otu_table.biom -o feat_rar_4000 -L 2,3,4,5,6,7
summarize_taxa.py -a -i feat_rar_4000/otu_table.biom -o feat_rar_4000_counts -L 2,3,4,5,6,7

