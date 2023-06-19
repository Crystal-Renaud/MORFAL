# MORFAL (Molecular Open Reading Frame AnaLyzer)
MORFAL is a pipeline build for the prediction of TIS and annotation of uORF (upstream open reading frame). 

## Requirements
- Samtools (tested on version 1.14) with hg19 fasta file
- Netstart (tested on version 1.0c) in singularity container
- VEP (tested on release-109.3)
- UTRannotator (<https://github.com/ImperialCardioGenetics/UTRannotator>) in singularity conntainer

## Usage

```
Usage: python3 /path/to/MORFAL.py [-gtf gtf_file] [-vcf vcf_file] [-O outdir]

Option:
-gtf,    --gtf  TEXT      /Path/to/gtf file
-vcf,    --vcf  TEXT      /Path/to/vcf files
-O,      --outdir  TEXT    /Path/to/output directory
-hg,     --hg19fasta  TEXT      /Path/to/hg19 fasta file
-net,    --netstart_container  TEXT      /Path/to/netstart singularity container
-UTR,    --UTRannotator  TEXT      /Path/to/UTRannotator singularity container
```

## Annotation Output
The output annotation is named output_annotation.txt . The other files are outputs from tools or from the concatenation and filtering stage.
The annotations can be read in 4 sections: VCF and GTF information, Netstart prediction and delta scores, UTRannotator annotations and GnomAD population frequencies. 
  
**VCF and GTF information** : 
  - *identifiant* : unique identifier created for each variant/transcript pair
  - *chr* : chromosome carrying of the variant
  - *position* : genomic position of the variant
  - *c.* : position of the variant relative to the start of the CDS
  - *ref* : affected nucleotide of wild-type sequence
  - *alt* : mutated nucleotide of the mutated sequence
  - *name_gene* : name of the gene affected by the variant
  - *gene_id* : Ensembl identifier of the gene
  - *transcript_id* : Ensembl identifier of the transcript
  - *strand* : strand direction -> + if strand is sense, - if strand is antisense

**Netstart prediction and delta scores**
  - *score_ATG_CDS(pos_ns)* : score predicted by Netstart for TIS and their position in the transcript
  - *delta_ATG_CDS* : delta score for TIS
  - *delta_min* : smallest delta score for transcript ATGs
  - *score_delta_min(pos)* : score predicted by Netstart for the smallest delta score and their position in the transcript
  - *delta_max* : highest delta score for transcript ATGs
  - *score_delta_max(pos)* : score predicted by Netstart for the highest delta score and their position in the transcript
  - *abolition(pos)* : score and position of lost ATG if available
  - *new(pos)* : score and position of created ATG if available
  - *diff_phy_new* : difference between CDS score and new ATG score

**UTRannotator annotations**
  - *Distance_CDS* : distance between start of uORF and CDS
  - *Distance_stop* : distance between start and stop of uORF
  - *Distance_Cap_start* : distance between CAP and start of uORF
  - *Kozak_context* : kozak sequence
  - *Kozak_strength* : strength of sequence kozak
  - *type* :  type of uORF: uORF, OutOfFrame_oORF, inFrame_oORF
  - *type_ref* : uORF reference type if this change
  - *Distance_stop_ref* : distance between start and stop of uORF reference if it change
  - *evidence* : whether the uORF disrupted has any translation evidence
  - *UTR_consequence* : output the variant consequences of a given 5 prime UTR variant: uAUG_gained, uAUG_lost, uSTOP_gained, uSTOP_lost, uFrameshift
  - *nb oORF_inframe* : number of existing inframe overlapping ORFs (inFrame_oORF) already within the 5 prime UTR
  - *nb oORF_outframe* : number of existing out-of-frame overlapping ORFs (inFrame_oORF) already within the 5 prime UTR
  - *nb uORF* : number of existing uORFs with a stop codon within the 5 prime UTR
  - *consequence* : variant consequence on the transcript
  - *impact* : strength of impact of variant on the transcript
   
**Frequencies GnomAD**
  - *Clin_sig* : ClinVar classification
  - *gnomad_combi* : frequency of existing variant in gnomAD exomes combined population
  - *gnomad_afr* : frequency of existing variant in gnomAD exomes African/American population
  - *gnomad_ame* : frequency of existing variant in gnomAD exomes American population
  - *gnomad_jew* : frequency of existing variant in gnomAD exomes Ashkenazi Jewish population
  - *gnomad_east_asi* : frequency of existing variant in gnomAD exomes East Asian population
  - *gnomad_fin* : frequency of existing variant in gnomAD exomes Finnish population
  - *gnomad_eur* : frequency of existing variant in gnomAD exomes Non-Finnish European population
  - *gnomad_other_combi* : frequency of existing variant in gnomAD exomes other combined population
  - *gnomad_south_asi* : frequency of existing variant in gnomAD exomes South Asian population
