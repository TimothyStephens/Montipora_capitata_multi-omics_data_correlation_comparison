# Re-analysis of genome V1 protein-vs-transcriptome results

## Download original data

Download original data from Zenodo: https://zenodo.org/record/6861688#.Y9xZhOzMJTa

Extract Data S1 and S3 into tab delimited tables for separate analysis. 

Upload to Coral server.

```bash
scp Supplemental_Data-S*.txt timothy@coral.rutgers.edu:/scratch/timothy/projects/0042_Mcapitata_TP3_PRJNA694677/03_Analysis/2023-02-02/01_genomeV1_Prot-vs-Trans/
```

Convert files to unix format and remove unused information.

```bash
dos2unix Supplemental_Data-S1.txt
dos2unix Supplemental_Data-S3.txt
```

**Proteomic data**

- Remove table name
- Remove columns other then "Accession" and the count data
- Rename "Accession" to "Name"
- Remove non-Mcap proteins (i.e., background contaminants database entries)
- Replace blank values with zeros ("0")

```bash
cat Supplemental_Data-S1.txt \
  | awk -F'\t' 'NR>1 {print $4"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32}' | sed -e 's/^Accession/Name/' \
  | awk -F'\t' 'NR==1 || $1~"^Montipora"' \
  | awk 'BEGIN{OFS=FS="\t"}{ for(i=1; i<=NF; i++){if($i==""){$i=0}}; print }' \
  > Mcapitata_V1_proteomic_data.tsv
```



**Transcriptomic data**

- Remove table name

```bash
cat Supplemental_Data-S3.txt \
  | awk -F'\t' 'NR>1' \
  > Mcapitata_V1_transcriptomic_data.tsv
```





## Differential Abundance Analysis

Run differential expression/abundance analysis.

```bash
./run_diffAccum_Analysis.sh
./run_diffExpr_Analysis.sh
```

Partition results based on timepoint.

```bash
PROT="Mcapitata_V1_proteomic_data.tsv_DiffAccumResults.txt"
awk -F'\t' 'NR==1 || ($6=="MC-289_T1-Amb" && $7=="MC-289_T1-HiT")' "${PROT}" > "${PROT}.TP1"
awk -F'\t' 'NR==1 || ($6=="MC-289_T3-Amb" && $7=="MC-289_T3-HiT")' "${PROT}" > "${PROT}.TP3"
awk -F'\t' 'NR==1 || ($6=="MC-289_T5-Amb" && $7=="MC-289_T5-HiT")' "${PROT}" > "${PROT}.TP5"

TRAN="Mcapitata_V1_transcriptomic_data.tsv_DiffExprResults.txt"
awk -F'\t' 'NR==1 || ($8=="MC-289_T1-Amb" && $9=="MC-289_T1-HiT")' "${TRAN}" > "${TRAN}.TP1"
awk -F'\t' 'NR==1 || ($8=="MC-289_T3-Amb" && $9=="MC-289_T3-HiT")' "${TRAN}" > "${TRAN}.TP3"
awk -F'\t' 'NR==1 || ($8=="MC-289_T5-Amb" && $9=="MC-289_T5-HiT")' "${TRAN}" > "${TRAN}.TP5"
```



## Create big data table

Create a master data table with all abundance (transcriptome+proteome), differential expression (transcriptome+proteome), and functional annotation information.

Convert transcript expression counts to TPM.

```bash
./run_Normalize_Expression_data.sh
```

Rename headers. R removes `-` from column names. Add this back in to prevent issues downstream. 

```bash
sed -e '1 s/MC.289_Field_289./MC-289_Field_289-/g' \
    -e '1 s/MC.289_T1./MC-289_T1-/g' \
    -e '1 s/MC.289_T3./MC-289_T3-/g' \
    -e '1 s/MC.289_T5./MC-289_T5-/g' \
    -i Mcapitata_V1_transcriptomic_data.tsv.normalized_counts.txt
```

```bash
sed -e '1 s/MC.289_Field_289./MC-289_Field_289-/g' \
    -e '1 s/MC.289_T1./MC-289_T1-/g' \
    -e '1 s/MC.289_T3./MC-289_T3-/g' \
    -e '1 s/MC.289_T5./MC-289_T5-/g' \
    -i Mcapitata_V1_transcriptomic_data.tsv.TPM.txt
```

Create master table.

Protein file header

> 1	Name
> 2	FC
> 3	VIP
> 4	pvalue
> 5	padj
> 6	Treatment1
> 7	Treatment2

Transcript file header

> 1	seqName
> 2	baseMean
> 3	log2FoldChange
> 4	lfcSE
> 5	stat
> 6	pvalue
> 7	padj
> 8	Treatment1
> 9	Treatment2

- Print fold change [FC] and if significant *p*-value.
    - Use "*p*-value" for **proteome** (un-adjusted)
    - Use "adjusted *p*-value" for **transcriptome**
- Also summarize significance of both data layers (i.e., DEP, DEG, or Both)

```bash
cat ../01_genomeV1-01_Annotations/Mcap.protein.fa.annotations \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(sed '1 s/MC/Protein-MC/g' Mcapitata_V1_proteomic_data.tsv) \
      -d $'0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0' \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(sed '1 s/MC/Transcript-MC/g' Mcapitata_V1_transcriptomic_data.tsv.TPM.txt) \
      -d $'0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0' \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(awk -F'\t' 'BEGIN{print "Name\tTP1-Protein-FC\tTP1-Protein-pvalue"}NR>1{
                             print $1"\t"$2"\t"$4}' "${PROT}.TP1") \
      -d $'0\tNA' \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(awk -F'\t' 'BEGIN{print "Name\tTP1-Transcript-FC\tTP1-Transcript-adj_pvalue"}NR>1{
                             print $1"\t"$3"\t"$7}' "${TRAN}.TP1") \
      -d $'0\tNA' \
  | awk -F'\t' '{ if(NR==1){
                    print $0"\tTP1-Regulation"
                  }else{ 
                         if($(NF-2)<0.05  && $NF>=0.05){print $0"\tDEP"} 
                    else if($(NF-2)>=0.05 && $NF<0.05 ){print $0"\tDEG"} 
                    else if($(NF-2)<0.05  && $NF<0.05 ){print $0"\tDEP&DEG"} 
                    else{print $0"\tnone"}
                  } }' \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(awk -F'\t' 'BEGIN{print "Name\tTP3-Protein-FC\tTP3-Protein-pvalue"}NR>1{
                             print $1"\t"$2"\t"$4}' "${PROT}.TP3") \
      -d $'0\tNA' \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(awk -F'\t' 'BEGIN{print "Name\tTP3-Transcript-FC\tTP3-Transcript-adj_pvalue"}NR>1{
                             print $1"\t"$3"\t"$7}' "${TRAN}.TP3") \
      -d $'0\tNA' \
  | awk -F'\t' '{ if(NR==1){
                    print $0"\tTP3-Regulation"
                  }else{ 
                         if($(NF-2)<0.05  && $NF>=0.05){print $0"\tDEP"} 
                    else if($(NF-2)>=0.05 && $NF<0.05 ){print $0"\tDEG"} 
                    else if($(NF-2)<0.05  && $NF<0.05 ){print $0"\tDEP&DEG"} 
                    else{print $0"\tnone"}
                  } }' \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(awk -F'\t' 'BEGIN{print "Name\tTP5-Protein-FC\tTP5-Protein-pvalue"}NR>1{
                             print $1"\t"$2"\t"$4}' "${PROT}.TP5") \
      -d $'0\tNA' \
  | ~/scripts/add_value_to_table_SQLite3.py \
      -a <(awk -F'\t' 'BEGIN{print "Name\tTP5-Transcript-FC\tTP5-Transcript-adj_pvalue"}NR>1{
                             print $1"\t"$3"\t"$7}' "${TRAN}.TP5") \
      -d $'0\tNA' \
  | awk -F'\t' '{ if(NR==1){
                    print $0"\tTP5-Regulation"
                  }else{ 
                         if($(NF-2)<0.05  && $NF>=0.05){print $0"\tDEP"} 
                    else if($(NF-2)>=0.05 && $NF<0.05 ){print $0"\tDEG"} 
                    else if($(NF-2)<0.05  && $NF<0.05 ){print $0"\tDEP&DEG"} 
                    else{print $0"\tnone"}
                  } }' \
  > Mcapitata_V1_multiomics_results.tsv
```



## Stress response genes

Check the activity (up vs. down) for Amanda's selected stress response genes. 

Extract here annotation from the original Table S2 and use it to annotate the genes in the new dataset.

```bash
~/scripts/add_value_to_table.py \
  -i Mcapitata_V1_multiomics_results.tsv \
  -a stress_gene_list.txt \
  -o stress_gene_list_multiomics_results.tsv \
  -d $'NA\tNA'
```

All appears to be consistent. 



## Plot prot-vs-trans data

Download results and plot values.

```bash
rsync -avzP timothy@coral.rutgers.edu:/scratch/timothy/projects/0042_Mcapitata_TP3_PRJNA694677/03_Analysis/2023-02-02/01_genomeV1-02_Prot-vs-Trans/*.tsv .
```







## Count stats

Count number of proteins in proteomic data.

```bash
awk 'NR>1' Mcapitata_V1_proteomic_data.tsv | wc -l
#4036
```



**Stress Response Genes**

```bash
# total "stress" genes
awk -F'\t' 'NR>1 && $57!="NA"' stress_gene_list_multiomics_results.tsv | wc -l
#138
```

Number of genes that are DEGs or DEPs at any of the timepoints.

```bash
awk -F'\t' 'NR>1 && $57!="NA" {print $45";"$50";"$55}' stress_gene_list_multiomics_results.tsv | grep -v 'none;none;none' | wc -l
#55
```

Number of genes which are DEPs&DEGs at any time point.

```bash
awk -F'\t' 'NR>1 && $57!="NA" {print $45";"$50";"$55}' stress_gene_list_multiomics_results.tsv | grep 'DEP&DEG' | wc -l
#9
```

Number of DEPs and DEGs at each time point.

```bash
## 138 Stree-response genes
#TP1
awk -F'\t' 'NR>1 && $57!="NA" {print $45}' stress_gene_list_multiomics_results.tsv | sort | uniq -c
      5 DEG
      8 DEP
    125 none

#TP3
awk -F'\t' 'NR>1 && $57!="NA" {print $50}' stress_gene_list_multiomics_results.tsv | sort | uniq -c
     11 DEG
     20 DEP
      9 DEP&DEG
     98 none

#TP5
awk -F'\t' 'NR>1 && $57!="NA" {print $55}' stress_gene_list_multiomics_results.tsv | sort | uniq -c
     16 DEG
      6 DEP
    116 none


## Total genes
#TP1
awk -F'\t' 'NR>1 {print $45}' stress_gene_list_multiomics_results.tsv | sort | uniq -c
    209 DEG
    292 DEP
      2 DEP&DEG
  62724 none

#TP3
awk -F'\t' 'NR>1 {print $50}' stress_gene_list_multiomics_results.tsv | sort | uniq -c
    695 DEG
    485 DEP
     50 DEP&DEG
  61997 none

#TP5
awk -F'\t' 'NR>1 {print $55}' stress_gene_list_multiomics_results.tsv | sort | uniq -c
    786 DEG
    110 DEP
      8 DEP&DEG
  62323 none
```







## KO numbers

Count the number of total genes with assigned KO numbers.

```bash
awk -F'\t' 'NR>1 && $5!="-"' stress_gene_list_multiomics_results.tsv | wc -l
#18684
```

Count number of proteins with assigned KO numbers that are part of major metabolic pathways.

```bash
awk -F'\t' 'NR>1 && $5!="-"' stress_gene_list_multiomics_results.tsv \
  | awk -F'\t' '{ split($5,a,","); for(i=1; i<=length(a); i++){print $1"\t"a[i]} }' \
  | ~/scripts/grepf_column.py \
    -f ../01_genomeV1-01_Annotations/target_KEGG_pathways.keg.KOlist \
    -c 2 \
  | cut -f1 | sort | uniq \
  | wc -l
#1925
```



Count number of proteome-identified proteins with assigned KO numbers.

```bash
~/scripts/grepf_column.py \
    -i stress_gene_list_multiomics_results.tsv \
    -f <(awk -F'\t' 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$5!="-"' \
  | wc -l
#2760
```

Count number of proteome-identified proteins with assigned KO numbers that are part of major metabolic pathways.

```bash
~/scripts/grepf_column.py \
    -i stress_gene_list_multiomics_results.tsv \
    -f <(awk -F'\t' 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$5!="-" { split($5,a,","); for(i=1; i<=length(a); i++){print $1"\t"a[i]} }' \
  | ~/scripts/grepf_column.py \
    -f ../01_genomeV1-01_Annotations/target_KEGG_pathways.keg.KOlist \
    -c 2 \
  | cut -f1 | sort | uniq \
  | wc -l
#414
```







## Direction of stress-response genes

Shared direction of stress response genes (i.e., do both proteins and trans go up or down together).

**TP1**

```bash
# Stress-response genes
awk -F'\t' 'NR>1 && $57!="NA" {print $41"\t"$43"\t"$45}' stress_gene_list_multiomics_results.tsv | awk '($1>0 && $2>0) || ($1<=0 && $2<=0)' | wc -l
#80


# All genes
~/scripts/grepf_column.py \
    -i Mcapitata_V1_multiomics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '{print $41"\t"$43"\t"$45}' \
  | awk '($1>0 && $2>0) || ($1<=0 && $2<=0)' \
  | wc -l
#2129
```

**TP3**

```bash
# Stress-response genes
awk -F'\t' 'NR>1 && $57!="NA" {print $46"\t"$48"\t"$50}' stress_gene_list_multiomics_results.tsv | awk '($1>0 && $2>0) || ($1<=0 && $2<=0)' | wc -l
#94


# All genes
~/scripts/grepf_column.py \
    -i Mcapitata_V1_multiomics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '{print $46"\t"$48"\t"$50}' \
  | awk '($1>0 && $2>0) || ($1<=0 && $2<=0)' \
  | wc -l
#2393
```











































Direction of regulation port-vs-trans for all proteins with proteomic evidence.

**TP1**

```bash
~/scripts/grepf_column.py \
    -i stress_gene_list_multiomics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$41>0 && $43>0' | wc -l
#898

~/scripts/grepf_column.py \
    -i stress_gene_list_multiomics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$41>0 && $43<0' | wc -l
#878

~/scripts/grepf_column.py \
    -i stress_gene_list_multiomics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$41<0 && $43<0' | wc -l
#1115

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$41<0 && $43>0' | wc -l
#925

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$41==0 || $43==0' | wc -l
#220
```



**TP3**

```bash
~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$46>0 && $48>0' | wc -l
#886

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$46>0 && $48<0' | wc -l
#626

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$46<0 && $48<0' | wc -l
#1364

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$46<0 && $48>0' | wc -l
#906

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$46==0 || $48==0' | wc -l
#254
```



**TP5**

```bash
~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$51>0 && $53>0' | wc -l
#965

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$51>0 && $53<0' | wc -l
#842

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$51<0 && $53<0' | wc -l
#1073

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$51<0 && $53>0' | wc -l
#936

~/scripts/grepf_column.py \
    -i Mcapitata_V1.omics_results.tsv \
    -f <(awk 'NR>1{print $1}' Mcapitata_V1_proteomic_data.tsv) \
  | awk -F'\t' '$51==0 || $53==0' | wc -l
#220
```







