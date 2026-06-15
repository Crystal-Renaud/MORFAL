import os
import subprocess
import argparse
import time
import random
import string


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-gtf', '--gtf', dest='gtf', type=str,
                        help='/path/to/gtf')
    parser.add_argument('-vcf', '--vcf', dest='vcf', type=str,
                        help='/path/to/vcf')
    parser.add_argument('-gen_fa', '--genome_fasta', dest='genome_fasta', type=str,
                        help='/path/to/fasta')
    parser.add_argument('-vep_cache', '--vep_cache', dest='vep_cache', type=str,
                        help='/path/to/vep_cache')
    parser.add_argument('-vep_file', '--UTRannotator_file', dest='UTRannotator_file', type=str,
                        help='/path/to/UTRannotator_file')
    parser.add_argument('-genome', '--genome', dest='genome', type=str,
                        help='GRCh37 or 38?')
    parser.add_argument('-MORFEE_ind', '--MORFEEdb_index', dest='MORFEEdb_index', type=str,
                        nargs='?', const =None, default=None, help='/path/to/MORFEEdb_index.tai')
    parser.add_argument('-MORFEE', '--MORFEEdb', dest='MORFEEdb', type=str,
                        nargs='?', const =None, default=None, help='/path/to/MORFEEdb')
    parser.add_argument('-O', '--outdir', dest='outdir', type=str,
                        help='/path/to/output directory')
    return parser.parse_args()


#function for inverting the antisense strand 
def invers_seq(seq):
    "creates the sense strand from the antisense strand"
    seq_sens = []
    for i in seq:
        seq_sens.append(i)
    seq_sens = seq_sens [::-1]
    comp = []
    dict_comp = { 'A':'T' , 'T':'A' , 'G':'C' , 'C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N'}  
    for base in seq_sens:  
        comp.append(dict_comp[base])
    comp = "".join(comp)
    return comp

#function for creating mutated sequence
def mut_seq(seq_wt, start_exon, pos_variant, alt, ref):
    "creates the mutated sequence from WT sequence"
    seq = "".join(seq_wt)
    start_del_ins = int(pos_variant) - int(start_exon)
    len_del = len(ref)
    if alt == "*" or alt == "-":
        seq_mut = seq[:start_del_ins].split() + seq[(start_del_ins + len_del):].split()
    else:
        seq_mut = seq[:start_del_ins].split() + alt.split() + seq[(start_del_ins + len_del):].split()
    seq_mut = "".join(seq_mut)
    return seq_mut

#function for integrating UTRannotator into the final table
def recup_vep(line, final_table, id):
    "Retrieves columns of interest from the VEP output file"
    nb_uORF = line[39].split()
    nb_uORF = "".join(nb_uORF)
    if line[35] != "-" and line[35] != "":
        ligne = line[35].replace(",", ":")
        UTR_annotation = ligne.split(":")
        #creation of dictionnary which takes all the information in the last colunm
        dico_utr_annot = {}
        if line[36] != "5_prime_UTR_uORF_frameshift_variant":
            for info in UTR_annotation:
                info = info.split("=")
                dico_utr_annot[info[0]] = info[1]
            if line[36] == "5_prime_UTR_premature_start_codon_gain_variant":
                final_table[id].update({"Distance_CDS":dico_utr_annot["DistanceToCDS"], "Distance_stop": dico_utr_annot["DistanceToStop"], "Distance_Cap": dico_utr_annot["CapDistanceToStart"],
                                        "Kozak_context": dico_utr_annot["KozakContext"], "Kozak_strength": dico_utr_annot["KozakStrength"], "type": dico_utr_annot["type"], 
                                        "type_ref": "-", "Distance_stop_ref": "-", "evidence": "-", "UTR_consequence": line[36], "oORF_inframe": line[37], 
                                        "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], "impact": line[13], "Clin_sig": line[32], 
                                        "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], "gnomad_jew": line[26], "gnomad_east_asi": line[27], 
                                        "gnomad_fin": line[28], "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
            elif line[36] == "5_prime_UTR_stop_codon_gain_variant":
                final_table[id].update({"Distance_CDS":dico_utr_annot["ref_StartDistanceToCDS"], "Distance_stop": int(dico_utr_annot["ref_StartDistanceToCDS"]) -  int(dico_utr_annot["newSTOPDistanceToCDS"]), 
                                        "Distance_Cap": "-", "Kozak_context": dico_utr_annot["KozakContext"], "Kozak_strength": dico_utr_annot["KozakStrength"], "type": "-", 
                                        "type_ref": dico_utr_annot["ref_type"], "Distance_stop_ref": "-", "evidence": "-", "UTR_consequence": line[36], "oORF_inframe": line[37], 
                                        "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], "impact": line[13], "Clin_sig": line[32], 
                                        "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], "gnomad_jew": line[26], "gnomad_east_asi": line[27], 
                                        "gnomad_fin": line[28], "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
            elif line[36] == "5_prime_UTR_stop_codon_loss_variant":
                final_table[id].update({"Distance_CDS": "-", "Distance_stop": "-", "Distance_Cap": "-", "Kozak_context": dico_utr_annot["KozakContext"], 
                                        "Kozak_strength": dico_utr_annot["KozakStrength"], "type": dico_utr_annot["FrameWithCDS"], "type_ref": "-", 
                                        "Distance_stop_ref": "-", "evidence": dico_utr_annot["Evidence"], "UTR_consequence": line[36], "oORF_inframe": line[37], 
                                        "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], "impact": line[13], "Clin_sig": line[32], 
                                        "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], "gnomad_jew": line[26], "gnomad_east_asi": line[27], 
                                        "gnomad_fin": line[28], "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
            else:
                final_table[id].update({"Distance_CDS":dico_utr_annot["DistanceToCDS"], "Distance_stop": dico_utr_annot["DistanceToStop"], "Distance_Cap": dico_utr_annot["CapDistanceToStart"],
                                        "Kozak_context": dico_utr_annot["KozakContext"], "Kozak_strength": dico_utr_annot["KozakStrength"], "type": dico_utr_annot["type"], 
                                        "type_ref": "-", "Distance_stop_ref": "-", "evidence": dico_utr_annot["Evidence"], "UTR_consequence": line[36], 
                                        "oORF_inframe": line[37], "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], "impact": line[13], 
                                        "Clin_sig": line[32], "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], "gnomad_jew": line[26], "gnomad_east_asi": line[27],
                                        "gnomad_fin": line[28], "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
        else:
            dico_utr_annot = {'ref_StartDistanceToCDS':[], 'alt_type':[], 'ref_type':[], 'alt_type_length':[], 'ref_type_length':[], 'KozakContext':[], 'KozakStrength':[], 'Evidence':[]}
            for info in UTR_annotation:
                info = info.split("=")
                dico_utr_annot[info[0]].append(info[1])
            final_table[id].update({"Distance_CDS":dico_utr_annot["ref_StartDistanceToCDS"], "Distance_stop": dico_utr_annot["alt_type_length"], "Distance_Cap": "-",
                                    "Kozak_context": dico_utr_annot["KozakContext"], "Kozak_strength": dico_utr_annot["KozakStrength"], "type": dico_utr_annot["alt_type"],
                                    "type_ref": dico_utr_annot["ref_type"], "Distance_stop_ref": dico_utr_annot["ref_type_length"], "evidence": dico_utr_annot["Evidence"], 
                                    "UTR_consequence": line[36], "oORF_inframe": line[37], "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], 
                                    "impact": line[13], "Clin_sig": line[32],  "gnomad_combi": line[23],  "gnomad_afr": line[24], "gnomad_ame": line[25], "gnomad_jew": line[26], 
                                    "gnomad_east_asi": line[27], "gnomad_fin": line[28], "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
    else:
        final_table[id].update({"Distance_CDS":"-", "Distance_stop":"-", "Distance_Cap": "-", "Kozak_context": "-", "Kozak_strength": "-", "type": "-", "evidence": "-", 
                                "UTR_consequence": line[36], "oORF_inframe": line[37], "oORF_outframe": line[38], "uORF": nb_uORF, "type_ref": "-", "Distance_stop_ref": "-", 
                                "consequence": line[6], "impact": line[13], "Clin_sig": line[32],  "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], 
                                "gnomad_jew": line[26], "gnomad_east_asi": line[27], "gnomad_fin": line[28],  "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
    return final_table[id]

#function for obtaining MORFEE from line number
def found_my_line(chr, pos, MORFEE_ind):
    with open(MORFEE_ind+chr+".txt.tai", 'r', encoding="utf-8") as file:
        search_line=[]
        for i, line in enumerate(file, start=0):
            if not line.startswith('#'):
                s_line = line.split("\t")
                start = int(s_line[1])
                if chr in line and pos == start :
                        my_line = i
                        search_line.append(my_line)
        return search_line

#function for adding MORFEE information  into the final table   
def add_MORFEE(line_morf, id, final_table, snv_type):
    clnsig = line_morf[42].split()
    clnsig = "".join(clnsig)
    if snv_type != "uTIS":
        final_table[id].update({"RefSeq_transcript": line_morf[5], "avsnp150":line_morf[7], "Mutation":line_morf[8], "Func.ensGene": line_morf[9], "orfSNVs_type": [line_morf[10]],
                                        "orfSNVs_frame": [line_morf[11]], "type_of_generated_ORF": [line_morf[12]], "Ratio_length_pred_obs": [line_morf[13]], "NewAALength": [line_morf[14]],
                                        "ORF_size":[line_morf[15]], "MORFEE_uTIS": [line_morf[16]], "MORFEE_dTIS": [line_morf[17]], "MORFEE_intTIS":[line_morf[18]], "MORFEE_uSTOP": [line_morf[19]],
                                        "MORFEE_dSTOP": [line_morf[20]], "MORFEE_intSTOP": [line_morf[21]], "MORFEE_uSTOP_creation": [line_morf[22]], "MORFEE_dSTOP_creation": [line_morf[23]],
                                        "MORFEE_intSTOP_creation": [line_morf[24]], "TIS_sequence": [line_morf[25]], "TIS_type": [line_morf[26]], "modification_type": [line_morf[27]], 
                                        "STOP_sequence": [line_morf[28]], "CDS_start_position": [line_morf[29]], "CDS_end_position": [line_morf[30]], "TIS_position": [line_morf[31]],
                                        "STOP_position": [line_morf[32]], "Kozak_sequence": [line_morf[33]], "Kozak_score": [line_morf[34]], "Kozak_strength_MORFEE": [line_morf[35]], 
                                        "TIS_c_pos": [line_morf[36]], "STOP_c_pos": [line_morf[37]], "gnomad312_AF": [line_morf[38]], "gnomad312_AF_popmax": [line_morf[39]], "CLNDN": [line_morf[40]], 
                                        "CLNDISB": [line_morf[41]], "CLNSIG": [clnsig]})
    else:
        final_table[id].update({"RefSeq_transcript": line_morf[5], "avsnp150":line_morf[7], "Mutation":line_morf[8], "Func.ensGene": line_morf[9], "orfSNVs_type": line_morf[10],
                                        "orfSNVs_frame": line_morf[11], "type_of_generated_ORF": line_morf[12], "Ratio_length_pred_obs": line_morf[13], "NewAALength": line_morf[14],
                                        "ORF_size": line_morf[15], "MORFEE_uTIS": line_morf[16], "MORFEE_dTIS": line_morf[17], "MORFEE_intTIS":line_morf[18], "MORFEE_uSTOP": line_morf[19],
                                        "MORFEE_dSTOP": line_morf[20], "MORFEE_intSTOP": line_morf[21], "MORFEE_uSTOP_creation": line_morf[22], "MORFEE_dSTOP_creation": line_morf[23],
                                        "MORFEE_intSTOP_creation": line_morf[24], "TIS_sequence": line_morf[25], "TIS_type": line_morf[26], "modification_type": line_morf[27], 
                                        "STOP_sequence": line_morf[28], "CDS_start_position": line_morf[29], "CDS_end_position": line_morf[30], "TIS_position": line_morf[31],
                                        "STOP_position": line_morf[32], "Kozak_sequence": line_morf[33], "Kozak_score": line_morf[34], "Kozak_strength_MORFEE": line_morf[35], 
                                        "TIS_c_pos": line_morf[36], "STOP_c_pos": line_morf[37], "gnomad312_AF": line_morf[38], "gnomad312_AF_popmax": line_morf[39], "CLNDN": line_morf[40], 
                                        "CLNDISB": line_morf[41], "CLNSIG": clnsig})
    return(final_table[id])
        
#function for adding MORFEE information when there are two lines for 1 variant
def append_MORFEE(line_morf, id, final_table):
    clnsig = line_morf[42].split()
    clnsig = "".join(clnsig)
    final_table[id]["orfSNVs_type"].append(line_morf[10])
    final_table[id]["orfSNVs_frame"].append(line_morf[11])
    final_table[id]["type_of_generated_ORF"].append(line_morf[12])
    final_table[id]["Ratio_length_pred_obs"].append(line_morf[13])
    final_table[id]["NewAALength"].append(line_morf[14])
    final_table[id]["ORF_size"].append(line_morf[15])
    final_table[id]["MORFEE_uTIS"].append(line_morf[16])
    final_table[id]["MORFEE_dTIS"].append(line_morf[17])
    final_table[id]["MORFEE_intTIS"].append(line_morf[18])
    final_table[id]["MORFEE_uSTOP"].append(line_morf[19])
    final_table[id]["MORFEE_dSTOP"].append(line_morf[20])
    final_table[id]["MORFEE_intSTOP"].append(line_morf[21])
    final_table[id]["MORFEE_uSTOP_creation"].append(line_morf[22])
    final_table[id]["MORFEE_dSTOP_creation"].append(line_morf[23])
    final_table[id]["MORFEE_intSTOP_creation"].append(line_morf[24])
    final_table[id]["TIS_sequence"].append(line_morf[25])
    final_table[id]["TIS_type"].append(line_morf[26])
    final_table[id]["modification_type"].append(line_morf[27])
    final_table[id]["STOP_sequence"].append(line_morf[28])
    final_table[id]["CDS_start_position"].append(line_morf[29])
    final_table[id]["CDS_end_position"].append(line_morf[30])
    final_table[id]["TIS_position"].append(line_morf[31])
    final_table[id]["STOP_position"].append(line_morf[32])
    final_table[id]["Kozak_sequence"].append(line_morf[33])
    final_table[id]["Kozak_score"].append(line_morf[34])
    final_table[id]["Kozak_strength_MORFEE"].append(line_morf[35])
    final_table[id]["TIS_c_pos"].append(line_morf[36])
    final_table[id]["STOP_c_pos"].append(line_morf[37])
    final_table[id]["gnomad312_AF"].append(line_morf[38])
    final_table[id]["gnomad312_AF_popmax"].append(line_morf[39])
    final_table[id]["CLNDN"].append(line_morf[40])
    final_table[id]["CLNDISB"].append(line_morf[41])
    final_table[id]["CLNSIG"].append(clnsig)
    return(final_table[id])

#principal function
def MORFAL(gtf, vcf, gen_fasta, outdir, vep_cache, UTRAnnotator_file, MORFEE, MORFEE_ind, genome):
    if os.path.exists(outdir+'/seq2netstart.fasta'):
        os.remove(outdir+'/seq2netstart.fasta')
    
    alphabet = list(string.ascii_letters)
    num = list(string.digits)
    punc = list(string.punctuation)
    alphanum = []

    for i in range(len(alphabet)):
        alphanum.append(alphabet[i])
    for i in range(len(num)):
        alphanum.append(num[i])
    for i in range(len(punc)):
        alphanum.append(punc[i])
    alphanum.remove("\\")
    
    #################### GTF storage in dictionary form
    list_chr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
    "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    dico_gtf = {chr: [] for chr in list_chr}

    with open (gtf) as f:
        gene_id = None
        transcript_id = None
        start_codon = None
        exon = []

        for line_GTF in f:
            if not line_GTF.startswith('#'): 
                line_GTF = line_GTF.split("\t")
                colonne = line_GTF[8].split("; ")
                #creation of a dictionary which takes all the information in the last column
                dico_colonne8 = {}
                for col in colonne:
                    col = col.split(" ")
                    dico_colonne8[col[0]] = col[1].replace('"', '')
                
                if line_GTF[2] == 'transcript':
                    #before changing transcript, create dico temp which will be added to dico gtf with the important info 
                    if (transcript_id) and (dico_colonne8["transcript_id"] != transcript_id):
                        dico_temp = {"gene_id": gene_id, "transcript_id": transcript_id, "gene_name": gene_name, "strand": strand,
                        "start_end_transcript": start_end_transcrit, "start_codon": start_codon, "exons": exon}
                        for chr in dico_gtf.keys():
                            if chrom == chr:
                                if dico_temp["transcript_id"] not in dico_gtf.values():
                                    dico_gtf[chr].append(dico_temp)
                        exon = []
                        start_codon = None
                        transcript_id = None
                        start_end_transcrit = ()
                    
                    #for the new transcript, retrieve the important info
                    if (dico_colonne8["transcript_id"] != transcript_id):
                        if (dico_colonne8["gene_id"] != gene_id): 
                            chrom = line_GTF[0]
                            strand = line_GTF[6]
                            gene_name = dico_colonne8["gene_name"]
                            gene_id = dico_colonne8["gene_id"]
                        transcript_id = dico_colonne8["transcript_id"]
                        start_end_transcrit = (int(line_GTF[3]), int(line_GTF[4]))


                #recovers the position of the start codon
                elif (line_GTF[2] == "start_codon"):
                    if (dico_colonne8["transcript_id"] == transcript_id):
                        if strand == '+':
                            start_codon = int(line_GTF[3])
                        elif strand == '-':
                            start_codon = int(line_GTF[4])
                
                #retrieve the start and end of each exon 
                elif (line_GTF[2] == "exon"):
                    if (dico_colonne8["transcript_id"] == transcript_id):
                        if len(exon) < 15:
                            exon.append((int(line_GTF[3]), int(line_GTF[4])))
        
        if transcript_id:
            dico_temp = {"gene_id": gene_id, "transcript_id": transcript_id, "gene_name": gene_name, "strand": strand,
                        "start_end_transcript": start_end_transcrit, "start_codon": start_codon, "exons": exon}
            for chr in dico_gtf.keys():
                if chrom == chr:
                    if dico_temp["transcript_id"] not in dico_gtf.values():
                        dico_gtf[chr].append(dico_temp)
            exon = []
            start_codon = None
            transcript_id = None
            start_end_transcrit = ()


    #################### VCF file processing
    with open(f'{vcf}', 'r') as file:
        canonical_start_ns = {}
        variant_ns ={}
        final_table = {} 
        for line_VCF in file:
            #recover columns of interest 
            if not line_VCF.startswith('#'):
                line_VCF = line_VCF.split()
                chrom = line_VCF[0]
                pos = int(line_VCF[1])
                ref = line_VCF[3]
                alt = line_VCF[4]
                
                #search in the dico_gtf according to the chr and if the variant is in the transcript
                for dico in dico_gtf[chrom]:
                    if dico["start_end_transcript"][0] <= pos <= dico["start_end_transcript"][1]:
                        #checks whether a start codon exists, otherwise the transcript is not kept
                        if dico["start_codon"]:
                            compt = 0
                            #checks whether the variant is in an exon 
                            for i in range(len(dico["exons"])):
                                if dico["exons"][i][0] <= pos <= dico["exons"][i][1]:
                                    compt += 1
                            
                            if compt != 0:
                                strand = (dico["strand"])
                                seq_merge_mut = []
                                seq_merge_wt = []
                                if strand == "+":
                                    start_transcript = dico["start_end_transcript"][0]
                                elif strand == "-":
                                    start_transcript = dico["start_end_transcript"][1]

                                #creation of a unique identifier for each transcript
                                id_var = []
                                while len(id_var) < 15:
                                    id_var.append(random.choice(alphanum))
                                id_var = "".join(id_var)
                                if id_var.startswith('"'):
                                    id_var = id_var.replace('"', "'")
                                
                                #calculate the size of the variant (positive if ins, negative if del)
                                if alt == "*" or alt == "-":
                                    for i in range(len(dico["exons"])):
                                        if dico["exons"][i][0] <= pos <= dico["exons"][i][1]:
                                            if pos + len(ref) < dico["exons"][i][1]:
                                                len_var = - len(ref)
                                            else:
                                                len_var = - (dico["exons"][i][1] - pos +1)
                                else:
                                    len_var = len(alt) - len(ref)

                                #gives the position of the start codon relative to the beginning of the transcript (useful for netstart)
                                pos_phy_ns = 0
                                found = 0
                                pos_var_ns = 0
                                find = 0
                                for start_seq, end_seq in dico["exons"]:
                                    if start_seq <= dico["start_codon"] <= end_seq:
                                        if strand == "+":
                                            pos_phy_ns += int(dico["start_codon"] - start_seq +1)
                                        if strand == "-":
                                            pos_phy_ns += int(abs(dico["start_codon"] - end_seq)+1)
                                        found += 1
                                    else:
                                        if found == 0:
                                            pos_phy_ns += int(end_seq - start_seq +1)

                                    #gives the position of the variant relative to the start of the transcript
                                    if start_seq <= pos <= end_seq:
                                        if strand == "+":
                                            pos_var_ns += int(pos - start_seq +1)
                                        if strand == "-":
                                            pos_var_ns += int(abs(pos - end_seq)+1)
                                        find += 1
                                    else:
                                        if find == 0:
                                            pos_var_ns += int(end_seq - start_seq +1)
                                    
                                    c = int(pos_var_ns - pos_phy_ns)
                                    if c >= 0:
                                        c += 1
                                       
                                    if c < 200 : 
                                        #uses SAMtools to create sequences 
                                        output = subprocess.run([f"samtools faidx {gen_fasta} {chrom}:{start_seq}-{end_seq}"], shell=True, capture_output=True, text=True)
                                        samtools_output = output.stdout
                                        samtools_output = samtools_output.split("\n")
                                        exon_seq_wt = "".join(samtools_output[1:])

                                        #uses the function to create the mutated sequence
                                        if (start_seq <= pos) and (pos<= end_seq):
                                            if (len(ref) == 1) and (len(alt) == 1) and (alt != "*") and (alt != "-"):
                                                exon_mut = mut_seq(exon_seq_wt, start_seq, pos, alt, ref)
                                                variant = "subs"
                                            elif len(ref) < len(alt):
                                                exon_mut = mut_seq(exon_seq_wt, start_seq, pos, alt, ref)
                                                variant = "ins"
                                            elif len(ref) > len(alt) or (alt == "*") or (alt == "-"):
                                                exon_mut = mut_seq(exon_seq_wt, start_seq, pos, alt, ref)
                                                variant = "del"
                                                c -= len_var
                                                if alt == "*":
                                                    alt = "-"
                                            elif len(ref) == len(alt) and len(ref) > 1:
                                                exon_mut = mut_seq(exon_seq_wt, start_seq, pos, alt, ref)
                                                variant = "subs"
                                        else: 
                                            exon_mut = exon_seq_wt
                                        
                                        
                                        #if antisense strand, use the sequence inversion function
                                        if strand == "-":
                                            exon_mut = invers_seq(exon_mut)
                                            exon_seq_wt = invers_seq(exon_seq_wt)
                                        seq_merge_mut.append(exon_mut)
                                        seq_merge_wt.append(exon_seq_wt)
                                if c < 200:         
                                    seq_merge_mut = "".join(seq_merge_mut)
                                    seq_merge_wt = "".join(seq_merge_wt)

                                    #creation of the fasta file to be used by Netstart
                                    samtools_output_wt = f">WT:{id_var}"
                                    samtools_output_mut = f">mute:{id_var}"
                                    
                                    output_test = samtools_output_wt.split() + seq_merge_wt.split() + samtools_output_mut.split() + seq_merge_mut.split()
                                    
                                    with open(outdir+'/seq2netstart.fasta', 'a') as output_file:
                                        for item in output_test:
                                            output_file.write(f"{item}\n")
                                    

                                    #creation dictionary for final file 
                                    dico_t = {"chrom": chrom, "pos": pos, "c.":c, "ref": ref, "alt": alt, "gene_name": dico['gene_name'], "gene_id": dico['gene_id'], "transcript_id": dico['transcript_id'], 
                                    "start_transcript": start_transcript, "start_codon": dico['start_codon'], "strand": strand, "variant": variant, "len_variant": len_var}
                                    final_table[id_var] = dico_t
                                    canonical_start_ns[id_var] = pos_phy_ns
                                    variant_ns[id_var] = pos_var_ns

                                    #if MORFEE analysis is requested, add the keys to the final dictionary
                                    if MORFEE != None:
                                        index_MORFEE = found_my_line(chrom, pos, MORFEE_ind)
                                        final_table[id_var].update({"index_MORFEE": index_MORFEE, "RefSeq_transcript": "-", "avsnp150":"-", "Mutation":"-", "Func.ensGene": "-", "orfSNVs_type": "-",
                                                 "orfSNVs_frame": "-", "type_of_generated_ORF": "-", "Ratio_length_pred_obs": "-", "NewAALength": "-", "ORF_size": "-", 
                                                 "MORFEE_uTIS": "-", "MORFEE_dTIS": "-", "MORFEE_intTIS": "-", "MORFEE_uSTOP": "-", "MORFEE_dSTOP": "-", "MORFEE_intSTOP": "-",
                                                 "MORFEE_uSTOP_creation": "-", "MORFEE_dSTOP_creation": "-","MORFEE_intSTOP_creation": "-", "TIS_sequence": "-", "TIS_type": "-", 
                                                 "modification_type": "-", "STOP_sequence": "-", "CDS_start_position": "-", "CDS_end_position": "-", "TIS_position": "-",
                                                 "STOP_position": "-", "Kozak_sequence": "-", "Kozak_score": "-", "Kozak_strength_MORFEE": "-", "TIS_c_pos": "-", "STOP_c_pos": "-", 
                                                 "gnomad312_AF": "-", "gnomad312_AF_popmax": "-", "CLNDN": "-", "CLNDISB": "-","CLNSIG": "-"})
                                  
    
    print(f"Netstart part")
    #################### Launch Netstart and put its output in a dictionary
    #Netstart only reads files < 490 Kb, so the fasta file must be divided into several smaller fasta files
    os.chdir(f"{outdir}")
    subprocess.run([f"split -l 68 {outdir}/seq2netstart.fasta split_seq2netstart_"], shell=True, text=True, check=False)

    #go to the Netstart container
    os.chdir("/opt/netstart-1.0c")
    #command to use Netstart and create an output file
    netstart_result = []
    file = []
    for filename in os.listdir(outdir):
        if filename.startswith("split_seq2netstart_"):
            file.append(filename)

    for var_fasta in file:
        fasta_path = "{}/{}".format(outdir, var_fasta)
        netstart = subprocess.run(
            ["/usr/bin/tcsh", "/opt/netstart-1.0c/netstart", "-vert", fasta_path],
            capture_output=True,
            text=True,
        )
        netstart_output = netstart.stdout
        netstart_output = netstart_output.split("\n")
        netstart_result.append(netstart_output)
    
    #store the positions and scores of each ATG for each transcript 
    dico_ns = {} 
    id = None
    name = None
    pos_ns_wt = []
    score_wt = []
    pos_ns_mut = []
    score_mut = []
    tout = []
    for ns_output in netstart_result:
        for line in ns_output:
            if line[:7] == '  Name:':
                if id and (id != line[-15:].strip()):
                    dico_temp = {"pos_ns_wt": pos_ns_wt, "score_wt": score_wt, "pos_ns_mut": pos_ns_mut, "score_mut": score_mut}
                    dico_ns[id] = dico_temp
                    pos_ns_wt = []
                    score_wt = []
                    pos_ns_mut = []
                    score_mut = []
                    id = None
                if (id != line[-15:].strip()):
                    id = line[-15:].strip()
                if name != line[8:]:
                    name = line[8:].strip()
            elif line[:6] == '   Pos':
                pass
            elif line[:1] == ' ':
                if name.startswith("WT"):
                    pos_ns_wt.append(int(line[1:6].strip()))
                    score_wt.append(float(line[10:15].strip()))
                    tout.append(["WT", id, f'{int(line[1:6].strip())}', f'{float(line[10:15].strip())}', f'{line[20:].strip()}'])
                elif name.startswith("mute"):
                    pos_ns_mut.append(int(line[1:6].strip()))
                    score_mut.append(float(line[10:15].strip()))
                    tout.append(["Mut", id, f'{int(line[1:6].strip())}', f'{float(line[10:15].strip())}', f'{line[20:].strip()}'])

        #store the last
        if id : 
            dico_temp = {"pos_ns_wt": pos_ns_wt, "score_wt": score_wt, "pos_ns_mut": pos_ns_mut, "score_mut": score_mut}
            dico_ns[id] = dico_temp
            pos_ns_wt = []
            score_wt = []
            pos_ns_mut = []
            score_mut = []
            id = None

    with open(outdir+"/positive_prediction.txt", "w") as f:
        f.write("info\tID\tpos\tscore\tpred\n")
        for liste in tout:
            if liste[-1] == "Yes":
                lign = "\t".join(liste)
                f.write(f"{lign}\n")

    #################### Delta score and Netstart output file
    result = []
    score_pos = []
    abo = None
    new = None
    for id_v, dico in dico_ns.items():
        #save start codon score
        for key, val in canonical_start_ns.items():
            if key == id_v:
                for wt in dico["pos_ns_wt"]:
                    if wt == val:
                        index_wt = dico["pos_ns_wt"].index(wt)
                        result.append([id_v, str(val), f"canonical TIS {(dico['score_wt'][index_wt])}"])
                        for id_variant, dico_t in final_table.items():
                            if id_variant == id_v:
                                final_table[id_v].update({"canonical_pos": val, "canonical_score": (dico['score_wt'][index_wt])})
                        if final_table[id_v]["variant"] == "subs":
                            #calculate delta score and add to result if non-zero
                            if wt in dico["pos_ns_mut"]:
                                index_mut = dico["pos_ns_mut"].index(wt)
                                delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                                delta = float("{:.3f}".format(delta))
                                final_table[id_v].update({"delta_start_cds": delta})
                            #check for loss of canonical ATG
                            elif wt not in dico["pos_ns_mut"]:
                                delta = "loss"
                                final_table[id_v].update({"delta_start_cds": delta})
                        if final_table[id_v]["variant"] == "del" or final_table[id_v]["variant"] == "ins":
                            if final_table[id_v]["c."] > 3:
                                if wt in dico["pos_ns_mut"]:
                                    index_mut = dico["pos_ns_mut"].index(wt)
                                    delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                                    delta = float("{:.3f}".format(delta))
                                    final_table[id_v].update({"delta_start_cds": delta})
                                #check for loss of canonical ATG
                                elif wt not in dico["pos_ns_mut"]:
                                    delta = "loss"
                                    final_table[id_v].update({"delta_start_cds": delta})
                            elif final_table[id_v]["c."] <= 3:
                                wt = int(wt) + final_table[id_v]["len_variant"]
                                if wt in dico["pos_ns_mut"]:
                                    index_mut = dico["pos_ns_mut"].index(wt)
                                    delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                                    delta = float("{:.3f}".format(delta))
                                    final_table[id_v].update({"delta_start_cds": delta})
                                #check for loss of canonical ATG
                                elif wt not in dico["pos_ns_mut"]:
                                    delta = "loss"
                                    final_table[id_v].update({"delta_start_cds": delta})


        if final_table[id_v]["variant"] == "subs":
            #delta score calculation and addition to result if not zero
            for atg_wt in dico["pos_ns_wt"]:
                if atg_wt in dico["pos_ns_mut"]:
                    index_wt = dico["pos_ns_wt"].index(atg_wt)
                    index_mut = dico["pos_ns_mut"].index(atg_wt)
                    delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                    delta = float("{:.3f}".format(delta))
                    if delta != 0.000:
                        result.append([id_v, str(atg_wt), f"{delta}% ({(dico['score_wt'][index_wt])})"])
                        score_pos.append([delta, ((dico['score_wt'][index_wt]), atg_wt)])
                #Check ATG loss
                elif atg_wt not in dico["pos_ns_mut"]:
                    index_wt = dico["pos_ns_wt"].index(atg_wt)
                    delta = "abolition"
                    result.append([id_v, str(atg_wt), f"{delta} ({(dico['score_wt'][index_wt])})"])
                    abo =((dico['score_wt'][index_wt]), atg_wt)
            #Check for creation
            for atg_mut in dico["pos_ns_mut"]:
                if atg_mut not in (dico["pos_ns_wt"]):
                    index_mut = dico["pos_ns_mut"].index(atg_mut)
                    delta = "new"
                    dif_phy_new = final_table[id_v]["canonical_score"] - (dico['score_mut'][index_mut])
                    dif_phy_new = float("{:.3f}".format(dif_phy_new))
                    result.append([id_v, str(atg_mut), f"{delta} ({(dico['score_mut'][index_mut])})"])
                    new = ((dico['score_mut'][index_mut]), atg_mut)
            
            for id_variant, dico_t in final_table.items():
                rien = ["-", ("-", "-")]
                if id_variant == id_v:
                    if score_pos:
                        final_table[id_v].update({"delta_min(pos)": min(score_pos), "delta_max(pos)": max(score_pos)})
                    else:
                        final_table[id_v].update({"delta_min(pos)": rien, "delta_max(pos)": rien})
                    if abo: 
                        final_table[id_v].update({"abolition(pos)": abo})
                    else:
                        final_table[id_v].update({"abolition(pos)": "--"})
                    if new:
                        final_table[id_v].update({"new(pos)": new, "diff_phy_new": dif_phy_new})
                    else:
                        final_table[id_v].update({"new(pos)": "--", "diff_phy_new": "-"})

                    final_table[id_v].update({"Distance_CDS":"-", "Distance_stop":"-", "Distance_Cap": "-", "Kozak_context": "-", "Kozak_strength": "-", "type": "-", "evidence": "-", 
                            "UTR_consequence": "-", "oORF_inframe": "-", "oORF_outframe": "-", "uORF": "-", "type_ref": "-", "Distance_stop_ref": "-", 
                            "consequence": "-", "impact":"-", "Clin_sig": "-", "gnomad_combi": "-", "gnomad_afr": "-", "gnomad_ame": "-", "gnomad_jew": "-",
                            "gnomad_east_asi": "-", "gnomad_fin": "-", "gnomad_eur": "-", "gnomad_other_combi": "-", "gnomad_south_asi": "-"})
                            
                    
                    abo = None
                    new = None
                    score_pos = []
        
        elif final_table[id_v]["variant"] == "del" or final_table[id_v]["variant"] == "ins":
            for i in range(len(dico["pos_ns_wt"])):
                #For ATG before variant
                if dico["pos_ns_wt"][i] < variant_ns[id_v]:
                    #calculate delta score and add to result if not zero
                    wt = dico["pos_ns_wt"][i]
                    if wt in dico["pos_ns_mut"]:
                        index_wt = dico["pos_ns_wt"].index(wt)
                        index_mut = dico["pos_ns_mut"].index(wt)
                        delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                        delta = float("{:.3f}".format(delta))
                        if delta != 0.000:
                            result.append([id_v, str(wt), f"{delta}% ({(dico['score_wt'][index_wt])})"])
                            score_pos.append([delta, ((dico['score_wt'][index_wt]), wt)])
                    #Check for ATG loss
                    elif wt not in dico["pos_ns_mut"]:
                        index_wt = dico["pos_ns_wt"].index(wt)
                        delta = "abolition"
                        result.append([id_v, str(wt), f"{delta} ({(dico['score_wt'][index_wt])})"])
                        abo =((dico['score_wt'][index_wt]), wt)
                #the same for ATG after variant
                else:
                    wt = dico["pos_ns_wt"][i]
                    wt_del_ins = dico["pos_ns_wt"][i] + final_table[id_v]["len_variant"]
                    if wt_del_ins in dico["pos_ns_mut"]:
                            index_wt = dico["pos_ns_wt"].index(wt)
                            index_mut = dico["pos_ns_mut"].index(wt_del_ins)
                            delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                            delta = float("{:.3f}".format(delta))
                            if delta != 0.000:
                                result.append([id_v, str(wt), f"{delta}% ({(dico['score_wt'][index_wt])})"])
                                score_pos.append([delta, ((dico['score_wt'][index_wt]), wt)])
                    elif wt_del_ins not in dico["pos_ns_mut"]:
                            index_wt = dico["pos_ns_wt"].index(wt)
                            delta = "abolition"
                            result.append([id_v, str(wt), f"{delta} ({(dico['score_wt'][index_wt])})"])
                            abo =((dico['score_wt'][index_wt]), wt)

            #Check if creation
            for i in range(len(dico["pos_ns_mut"])):
                if dico["pos_ns_mut"][i] < variant_ns[id_v]:
                    mut = dico["pos_ns_mut"][i]
                    if mut not in (dico["pos_ns_wt"]):
                        index_mut = dico["pos_ns_mut"].index(mut)
                        delta = "new"
                        dif_phy_new = final_table[id_v]["canonical_score"] - (dico['score_mut'][index_mut])
                        dif_phy_new = float("{:.3f}".format(dif_phy_new))
                        result.append([id_v, str(mut), f"{delta} ({(dico['score_mut'][index_mut])})"])
                        new = ((dico['score_mut'][index_mut]), mut)
                else:
                    mut = dico["pos_ns_mut"][i]
                    mut_del_ins = dico["pos_ns_mut"][i] - final_table[id_v]["len_variant"]
                    if mut_del_ins not in (dico["pos_ns_wt"]):
                        index_mut = dico["pos_ns_mut"].index(mut)
                        delta = "new"
                        dif_phy_new = final_table[id_v]["canonical_score"] - (dico['score_mut'][index_mut])
                        dif_phy_new = float("{:.3f}".format(dif_phy_new))
                        result.append([id_v, str(mut), f"{delta} ({(dico['score_mut'][index_mut])})"])
                        new = ((dico['score_mut'][index_mut]), mut)
            
            for id_variant, dico_t in final_table.items():
                rien = ["-", ("-", "-")]
                if id_variant == id_v:
                    if score_pos:
                        final_table[id_v].update({"delta_min(pos)": min(score_pos), "delta_max(pos)": max(score_pos)})
                    else:
                        final_table[id_v].update({"delta_min(pos)": rien, "delta_max(pos)": rien})
                    if abo: 
                        final_table[id_v].update({"abolition(pos)": abo})
                    else:
                        final_table[id_v].update({"abolition(pos)": "--"})
                    if new:
                        final_table[id_v].update({"new(pos)": new, "diff_phy_new": dif_phy_new})
                    else:
                        final_table[id_v].update({"new(pos)": "--", "diff_phy_new": "-"})
                        
                    #Add keys used for UTRAnnotator for all IDs
                    final_table[id_v].update({"Distance_CDS":"-", "Distance_stop":"-", "Distance_Cap": "-", "Kozak_context": "-", "Kozak_strength": "-", "type": "-", "evidence": "-", 
                            "UTR_consequence": "-", "oORF_inframe": "-", "oORF_outframe": "-", "uORF": "-", "type_ref": "-", "Distance_stop_ref": "-", 
                            "consequence": "-", "impact":"-", "Clin_sig": "-", "gnomad_combi": "-", "gnomad_afr": "-", "gnomad_ame": "-", "gnomad_jew": "-",
                            "gnomad_east_asi": "-", "gnomad_fin": "-", "gnomad_eur": "-", "gnomad_other_combi": "-", "gnomad_south_asi": "-"})
                                  
                    abo = None
                    new = None
                    score_pos = []
                   

    with open(outdir+"/delta.txt", 'w') as file:
        file.write("ID des transcrits var\tpos\tdelta(score WT ou new)\n")
        for r in result:
            line="\t".join(r)
            file.write(f"{line}\n")
    
    with open(outdir+"/UTR_FirstExons.vcf", "w") as file_out:
        file_out.write("#File containing the variants in 5' UTR and near the start codon\n")
        ligne = []
        liste =[]
        for id, dico in final_table.items():
            ligne.append(f"{dico['chrom']}\t{dico['pos']}\t.\t{dico['ref']}\t{dico['alt']}\n")
            for element in ligne:
                if element not in liste:
                    liste.append(element)
        print(len(liste))
        file_out.write(''.join(liste))

    print("UTRAnnotator")
    #################### Use UTRannotator on our variants in the UTR and near the start codon
    vep_res = subprocess.run([
        "/opt/vep/src/ensembl-vep/vep",
        "--fork", "4",
        "--fasta", gen_fasta,
        "--cache",
        "--offline",
        "--merged",
        "--format", "vcf",
        "--tab",
        "--force_overwrite",
        "--dir", vep_cache,
        "--af_gnomade",
        "--assembly", genome,
        "--use_transcript_ref",
        "--plugin", "UTRAnnotator,file={}".format(UTRAnnotator_file),
        "-i", "{}/UTR_FirstExons.vcf".format(outdir),
        "-o", "{}/UTR.annotated.vcf".format(outdir),], capture_output=True, text=True)
    
    
    with open(f"{outdir}/UTR.annotated.vcf") as file:
        for line in file:
            if not line.startswith('#'): 
                line = line.split("\t")
                if line[18] == "Ensembl":
                    if line[14] == "-":
                        uploaded = line[0].split("_")
                        for id, dico in final_table.items():
                            if uploaded[0] == dico["chrom"]:
                                if line[3] == dico["gene_id"][:15]:
                                    if (len(dico["alt"]) == len(dico["ref"])) and (int(uploaded[1]) == dico["pos"]):
                                        if line[4] == dico["transcript_id"][:15]:
                                            recup_vep(line, final_table, id)
                                    elif (len(dico["alt"]) != len(dico["ref"])):
                                        if int(uploaded[1]) == dico["pos"]:
                                            if line[4] == dico["transcript_id"][:15]:
                                                nuc = uploaded[2].split("/")
                                                if len(dico["ref"]) > 1 and len(dico["alt"]) == 1:
                                                    if dico["alt"] == "-":
                                                        if nuc[0] == dico["ref"]:
                                                            recup_vep(line, final_table, id)
                                        if (int(uploaded[1])-1 == dico["pos"]):
                                            if line[4] == dico["transcript_id"][:15]:
                                                nuc = uploaded[2].split("/")
                                                if len(dico["ref"]) > 1 and len(dico["alt"]) == 1 and (dico["alt"] != "-"):
                                                    if nuc[0] == dico["ref"][1:]:
                                                        recup_vep(line, final_table, id)
                                                elif len(dico["ref"]) > 1 and len(dico["alt"]) > 1:
                                                    if (nuc[0] == dico["ref"][1:]) and (nuc[1] == dico["alt"][1:]):
                                                        recup_vep(line, final_table, id)
                                                elif len(dico["ref"]) == 1 and len(dico["alt"]) > 1:
                                                    if nuc[1] == dico["alt"][1:]:
                                                        recup_vep(line, final_table, id)

    
    ###############Add MORFEE if asked
    if MORFEE != None:
        print(f"MORFEE part")
        
        list_chr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        index_to_ids = {chr: {} for chr in list_chr}
        for id, dico in final_table.items():
            chr = dico["chrom"]
            for ind in dico["index_MORFEE"]:
                index_to_ids[chr].setdefault(ind, []).append(id)

        for chr, dico_ind in index_to_ids.items():
            if dico_ind != {}:
                with open(MORFEE+chr+".txt", "r", encoding="utf-8") as f:
                    for i, line in enumerate(f, start=0):
                        if i not in index_to_ids[chr]:
                            continue  # Skip the lines that are not relevant

                        line_morf = line.split("\t")
                        for id in index_to_ids[chr][i]:
                            dico = final_table[id]
                            transcript = dico["transcript_id"][:15]
                            alt = dico["alt"]

                            if line_morf[4][:15] == transcript and line_morf[3] == alt:
                                snv_type = line_morf[10]
                                if snv_type == "new_uSTOP":
                                    add_MORFEE(line_morf, id, final_table, snv_type)
                                elif snv_type == "uSTOP":
                                    if dico["Func.ensGene"] == '-':
                                        add_MORFEE(line_morf, id, final_table, snv_type)
                                    else:
                                        append_MORFEE(line_morf, id, final_table)
                                elif snv_type == "uTIS":
                                    if dico["Func.ensGene"] == '-':
                                        add_MORFEE(line_morf, id, final_table, snv_type)
                                    else:
                                        append_MORFEE(line_morf, id, final_table)


        

    print(f"Final_output")                                        
    ##########Write final fichier
    if MORFEE == None: 
        with open(outdir+"/output_annotation.txt", "w") as out_file:
            out_file.write("##MORFAL output\n"\
            f"##Analysis in {genome}\n"\
            "##Column descriptions:\n"\
            "#start_transcript : beginning of the transcript on the chr\n"\
            "#start_CDS : beginning of the protein on the chr\n"\
            "#score_start_CDS(pos_ns) : score predicted by Netstart and position in the transcript\n"\
            "#delta_min : smallest percentage of variation between the reference sequence scores and those of the mutated sequence\n"\
            "#score_delta_min : reference sequence score associated with delta\n"\
            "#delta_max : highest percentage of variation between the reference sequence scores and those of the mutated sequence\n"\
            "#score_delta_max : reference sequence score associated with delta\n"\
            "#abolition(pos): score and position of the lost ATG\n"\
            "#new(pos) : score and position of the gain ATG\n"\
            "#difference_canonical_vs_new_ATG : difference between the score of the canonical atg and the one created\n"\
            "#evidence : whether the uORF disrupted has any translation evidence\n"\
            "#identifiant : indentifier that links the mutation to its sequence in FASTA file\n"\
            "chr\tposition\tc.\tref\talt\tname_gene\tgene_id\ttranscript_id\tstrand\tstart_transcript\tstart_CDS\tscore_ATG_CDS(pos_ns)\t"\
                "delta_ATG_CDS\tdelta_min\tscore_delta_min(pos)\tdelta_max\tscore_delta_max(pos)\tabolition(pos)\tnew(pos)\tdifference_canonical_vs_new_ATG\tDistance_CDS\t"\
                    "Distance_stop\tDistance_Cap_start\tKozak_context\tKozak_strength\ttype\ttype_ref(if type change)\tDistance_stop_ref(if uFrameshift)\t"\
                        "evidence\tUTR_consequence\toORF_inframe\toORF_outframe\tuORF\tconsequence\timpact\tClin_sig\tgnomad_combi\tgnomad_afr\t"\
                            "gnomad_ame\tgnomad_jew\tgnomad_east_asi\tgnomad_fin\tgnomad_eur\tgnomad_other_combi\tgnomad_south_asi\tidentifier\n""")

            for id, dico in final_table.items():
                #print(f"{id} donne {dico}")
                out_file.write(f"{dico['chrom']}\t{dico['pos']}\t{dico['c.']}\t{dico['ref']}\t{dico['alt']}\t{dico['gene_name']}\t{dico['gene_id']}\t"\
                    f"{dico['transcript_id']}\t{dico['strand']}\t{dico['start_transcript']}\t{dico['start_codon']}\t{dico['canonical_score']} ({dico['canonical_pos']})\t"\
                        f"{dico['delta_start_cds']}%\t{dico['delta_min(pos)'][0]}%\t{dico['delta_min(pos)'][1][0]} ({dico['delta_min(pos)'][1][1]})\t{dico['delta_max(pos)'][0]}%\t"\
                            f"{dico['delta_max(pos)'][1][0]} ({dico['delta_max(pos)'][1][1]})\t{dico['abolition(pos)'][0]} ({dico['abolition(pos)'][1]})\t"\
                                f"{dico['new(pos)'][0]} ({dico['new(pos)'][1]})\t{dico['diff_phy_new']}\t{dico['Distance_CDS']}\t{dico['Distance_stop']}\t"\
                                    f"{dico['Distance_Cap']}\t{dico['Kozak_context']}\t{dico['Kozak_strength']}\t{dico['type']}\t{dico['type_ref']}\t{dico['Distance_stop_ref']}\t"\
                                        f"{dico['evidence']}\t{dico['UTR_consequence']}\t{dico['oORF_inframe']}\t{dico['oORF_outframe']}\t{dico['uORF']}\t{dico['consequence']}\t{dico['impact']}\t"\
                                            f"{dico['Clin_sig']}\t{dico['gnomad_combi']}\t{dico['gnomad_afr']}\t{dico['gnomad_ame']}\t{dico['gnomad_jew']}\t"\
                                                f"{dico['gnomad_east_asi']}\t{dico['gnomad_fin']}\t{dico['gnomad_eur']}\t{dico['gnomad_other_combi']}\t{dico['gnomad_south_asi']}\t{id}\n")
    
    else:
        with open(outdir+"/output_annotation.txt", "w") as out_file:
            out_file.write("##MORFAL output\n"\
            f"##Analysis in {genome}\n"\
            "##Column descriptions:\n"\
            "#start_transcript : beginning of the transcript on the chr\n"\
            "#start_CDS : beginning of the protein on the chr\n"\
            "#score_start_CDS(pos_ns) : score predicted by Netstart and position in the transcript\n"\
            "#delta_min : smallest percentage of variation between the reference sequence scores and those of the mutated sequence\n"\
            "#score_delta_min : reference sequence score associated with delta\n"\
            "#delta_max : highest percentage of variation between the reference sequence scores and those of the mutated sequence\n"\
            "#score_delta_max : reference sequence score associated with delta\n"\
            "#abolition(pos): score and position of the lost ATG\n"\
            "#new(pos) : score and position of the gain ATG\n"\
            "#difference_canonical_vs_new_ATG : difference between the score of the canonical atg and the one created\n"\
            "#evidence : whether the uORF disrupted has any translation evidence\n"\
            "#identifiant : indentifier that links the mutation to its sequence in FASTA file\n"\
            "chr\tposition\tc.\tref\talt\tname_gene\tgene_id\ttranscript_id\tstrand\tstart_transcript\tstart_CDS\tscore_ATG_CDS(pos_ns)\t"\
                "delta_ATG_CDS\tdelta_min\tscore_delta_min(pos)\tdelta_max\tscore_delta_max(pos)\tabolition(pos)\tnew(pos)\tdifference_canonical_vs_new_ATG\tDistance_CDS\t"\
                    "Distance_stop\tDistance_Cap_start\tKozak_context\tKozak_strength\ttype\ttype_ref(if type change)\tDistance_stop_ref(if uFrameshift)\t"\
                        "evidence\tUTR_consequence\toORF_inframe\toORF_outframe\tuORF\tconsequence\timpact\tClin_sig\tgnomad_combi\tgnomad_afr\t"\
                            "gnomad_ame\tgnomad_jew\tgnomad_east_asi\tgnomad_fin\tgnomad_eur\tgnomad_other_combi\tgnomad_south_asi\tRefSeq_transcript\tavsnp150\t"\
                                "Mutation\tFunc.ensGene\torfSNVs_type\torfSNVs_frame\ttype_of_generated_ORF\tRatio_length_pred_obs\tNewAALength\tORF_size\tMORFEE_uTIS\t"\
                                     "MORFEE_dTIS\tMORFEE_intTIS\tMORFEE_uSTOP\tMORFEE_dSTOP\tMORFEE_intSTOP\tMORFEE_uSTOP_creation\tMORFEE_dSTOP_creation\tMORFEE_intSTOP_creation\t"\
                                        "TIS_sequence\tTIS_type\tmodification_type\tSTOP_sequence\tCDS_start_position\tCDS_end_position\tTIS_position\tSTOP_position\t"\
                                            "Kozak_sequence\tKozak_score\tKozak_strength_MORFEE\tTIS_c_pos\tSTOP_c_pos\tgnomad312_AF\tgnomad312_AF_popmax\tCLNDN\tCLNDISB\tCLNSIG\tidentifier\n""")

            for id, dico in final_table.items():
                out_file.write(f"{dico['chrom']}\t{dico['pos']}\t{dico['c.']}\t{dico['ref']}\t{dico['alt']}\t{dico['gene_name']}\t{dico['gene_id']}\t"\
                    f"{dico['transcript_id']}\t{dico['strand']}\t{dico['start_transcript']}\t{dico['start_codon']}\t{dico['canonical_score']} ({dico['canonical_pos']})\t"\
                        f"{dico['delta_start_cds']}%\t{dico['delta_min(pos)'][0]}%\t{dico['delta_min(pos)'][1][0]} ({dico['delta_min(pos)'][1][1]})\t{dico['delta_max(pos)'][0]}%\t"\
                            f"{dico['delta_max(pos)'][1][0]} ({dico['delta_max(pos)'][1][1]})\t{dico['abolition(pos)'][0]} ({dico['abolition(pos)'][1]})\t"\
                                f"{dico['new(pos)'][0]} ({dico['new(pos)'][1]})\t{dico['diff_phy_new']}\t{dico['Distance_CDS']}\t{dico['Distance_stop']}\t"\
                                    f"{dico['Distance_Cap']}\t{dico['Kozak_context']}\t{dico['Kozak_strength']}\t{dico['type']}\t{dico['type_ref']}\t{dico['Distance_stop_ref']}\t"\
                                        f"{dico['evidence']}\t{dico['UTR_consequence']}\t{dico['oORF_inframe']}\t{dico['oORF_outframe']}\t{dico['uORF']}\t{dico['consequence']}\t{dico['impact']}\t"\
                                            f"{dico['Clin_sig']}\t{dico['gnomad_combi']}\t{dico['gnomad_afr']}\t{dico['gnomad_ame']}\t{dico['gnomad_jew']}\t"\
                                                f"{dico['gnomad_east_asi']}\t{dico['gnomad_fin']}\t{dico['gnomad_eur']}\t{dico['gnomad_other_combi']}\t{dico['gnomad_south_asi']}\t{dico['RefSeq_transcript']}\t"\
                                                    f"{dico['avsnp150']}\t{dico['Mutation']}\t{dico['Func.ensGene']}\t{dico['orfSNVs_type']}\t{dico['orfSNVs_frame']}\t{dico['type_of_generated_ORF']}\t"\
                                                        f"{dico['Ratio_length_pred_obs']}\t{dico['NewAALength']}\t{dico['ORF_size']}\t{dico['MORFEE_uTIS']}\t{dico['MORFEE_dTIS']}\t"\
                                                            f"{dico['MORFEE_intTIS']}\t{dico['MORFEE_uSTOP']}\t{dico['MORFEE_dSTOP']}\t{dico['MORFEE_intSTOP']}\t{dico['MORFEE_uSTOP_creation']}\t"\
                                                                f"{dico['MORFEE_dSTOP_creation']}\t{dico['MORFEE_intSTOP_creation']}\t{dico['TIS_sequence']}\t{dico['TIS_type']}\t{dico['modification_type']}\t"\
                                                                    f"{dico['STOP_sequence']}\t{dico['CDS_start_position']}\t{dico['CDS_end_position']}\t{dico['TIS_position']}\t{dico['STOP_position']}\t"\
                                                                        f"{dico['Kozak_sequence']}\t{dico['Kozak_score']}\t{dico['Kozak_strength_MORFEE']}\t{dico['TIS_c_pos']}\t{dico['STOP_c_pos']}\t"\
                                                                            f"{dico['gnomad312_AF']}\t{dico['gnomad312_AF_popmax']}\t{dico['CLNDN']}\t{dico['CLNDISB']}\t{dico['CLNSIG']}\t{id}\n")
    

    subprocess.run([f"rm {outdir}/split_seq2netstart_*"], shell=True, text=True)



def main():
    start_time = time.time()

    # get arguments
    args = get_arguments()
    gtf = os.path.realpath(args.gtf)
    vcf = os.path.realpath(args.vcf)
    genome_fasta = os.path.realpath(args.genome_fasta)
    vep_cache = os.path.realpath(args.vep_cache)
    UTRannotator_file = os.path.realpath(args.UTRannotator_file)
    genome = args.genome
    MORFEEdb = args.MORFEEdb
    MORFEE_ind = args.MORFEEdb_index
    outdir = os.path.realpath(args.outdir)


    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # MORFAL
    MORFAL(gtf=gtf, vcf=vcf, gen_fasta=genome_fasta, vep_cache=vep_cache, UTRAnnotator_file=UTRannotator_file, genome=genome, MORFEE=MORFEEdb, MORFEE_ind=MORFEE_ind, outdir=outdir)

    step = f"--- Program {os.path.basename(__file__)} executed in {round((float(time.time() - start_time) / 60), 2)} minutes ---"; print(step)

    # MORFEE=MORFEEdb,

if __name__ == "__main__":
    main()
