# Genome V1 annotations

Annotate The *M. capitata* predicted proteins so that we know what their functions are.



Run `DIAMOND` against nr (2022_07) to annotate the proteins.

```bash
./run_diamond_nr.sh
```



Filter `DIAMOND` results so that we have just the top hits.

```bash
./blast_top_hits.py -n 1 -s 1 \
  -i Mcap.protein.fa.blastp_nr.outfmt6.gz \
  -o Mcap.protein.fa.blastp_nr.outfmt6.tophit
```



Assign KEGG IDs using `eggnog-mapper`.

```bash
sed -e 's/*/X/g' Mcap.protein.fa > Mcap.protein.cleaned.fa

./run_eggnog-mapper.sh
./run_InterProScan.sh
```



Extract annotations and add them to proteins (i.e., top hits + KEGG annotations)

```bash
cat Mcap.protein.fa \
| awk 'BEGIN{print "Name"}$1~"^>"{gsub(">","",$1); print $1}' \
| ~/scripts/add_value_to_table.py \
   -d $'NA\tNA\tNA' \
   -a <(cat Mcap.protein.fa.blastp_nr.outfmt6.tophit \
          | awk -F'\t' 'BEGIN{print "Name\ttop_hit_id\ttop_hit_evalue\ttop_hit_description"}
                             {print $1"\t"$2"\t"$11"\t"$15}'\
          | sed -e 's/ >.*//') \
| ~/scripts/add_value_to_table.py \
   -d '-' \
   -a <(cat Mcap.protein.cleaned.fa.emapper.annotations \
          | awk -F'\t' 'BEGIN{print "Name\tKO_numbers"} NR>1{print $1"\t"$12}' \
          | sed -e 's/ko://g') \
   > Mcap.protein.fa.annotations
```



Extract KEGG KO numbers from "target" metabolism pathways.

```bash
~/scripts/grepf_column.py \
  --keep_header \
  -i /scratch/timothy/databases/KEGG/ko00001.keg.unpacked \
  -c 5 \
  -f <(cut -f2 target_KEGG_pathways.tsv | sed -e 's/ko//') \
  -o target_KEGG_pathways.keg

awk -F'\t' 'NR>1{print $7}' target_KEGG_pathways.keg | sort | uniq > target_KEGG_pathways.keg.KOlist
```



