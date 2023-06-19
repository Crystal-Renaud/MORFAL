import os
import subprocess
import pandas as pd
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
    parser.add_argument('-O', '--outdir', dest='outdir', type=str,
                        help='/path/to/output directory')
    parser.add_argument('-s', '--hg19fasta', dest='hg19fasta', type=str,
                        help='/path/to/hg19 fasta')
    parser.add_argument('-net', '--netstart container', dest='netstart_container', type=str,
                        help='/path/to/netstart container')
    parser.add_argument('-UTR', '--UTRannotator', dest='UTRannotator', type=str,
                        help='/path/to/UTRannotator')
    return parser.parse_args()


#inverse le brin antisens en brin sens
def invers_seq(seq):
    seq_sens = []
    for i in seq:
        seq_sens.append(i)
    seq_sens = seq_sens [::-1]
    comp = []
    dict_comp = { 'A':'T' , 'T':'A' , 'G':'C' , 'C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c' }  
    for base in seq_sens:  
        comp.append(dict_comp[base])
    comp = "".join(comp)
    return comp

#fonction pour creer sequence mute
def mut_seq(seq_wt, start_exon, pos_variant, alt, ref):
    seq = "".join(seq_wt)
    start_del_ins = int(pos_variant) - int(start_exon)
    len_del = len(ref)
    if alt == "*" or alt == "-":
        seq_mut = seq[:start_del_ins].split() + seq[(start_del_ins + len_del):].split()
    else:
        seq_mut = seq[:start_del_ins].split() + alt.split() + seq[(start_del_ins + len_del):].split()
    seq_mut = "".join(seq_mut)
    return seq_mut

def recup_vep(line, final_table, id):
    nb_uORF = line[39].split()
    nb_uORF = "".join(nb_uORF)
    if line[35] != "-" and line[35] != "":
        ligne = line[35].replace(",", ":")
        UTR_annotation = ligne.split(":")
        #crea dico qui prend toutes les info de la derniere colonne
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
                                        "somatic": line[33], "pheno": line[34], "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], 
                                        "gnomad_jew": line[26], "gnomad_east_asi": line[27], "gnomad_fin": line[28], "gnomad_eur": line[29], 
                                        "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
            elif line[36] == "5_prime_UTR_stop_codon_gain_variant":
                final_table[id].update({"Distance_CDS":dico_utr_annot["ref_StartDistanceToCDS"], "Distance_stop": int(dico_utr_annot["ref_StartDistanceToCDS"]) -  int(dico_utr_annot["newSTOPDistanceToCDS"]), 
                                        "Distance_Cap": "-", "Kozak_context": dico_utr_annot["KozakContext"], "Kozak_strength": dico_utr_annot["KozakStrength"], "type": "-", 
                                        "type_ref": dico_utr_annot["ref_type"], "Distance_stop_ref": "-", "evidence": "-", "UTR_consequence": line[36], "oORF_inframe": line[37], 
                                        "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], "impact": line[13], "Clin_sig": line[32], 
                                        "somatic": line[33], "pheno": line[34], "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], 
                                        "gnomad_jew": line[26], "gnomad_east_asi": line[27], "gnomad_fin": line[28], "gnomad_eur": line[29], 
                                        "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
            elif line[36] == "5_prime_UTR_stop_codon_loss_variant":
                final_table[id].update({"Distance_CDS": "-", "Distance_stop": "-", "Distance_Cap": "-", "Kozak_context": dico_utr_annot["KozakContext"], 
                                        "Kozak_strength": dico_utr_annot["KozakStrength"], "type": dico_utr_annot["FrameWithCDS"], "type_ref": "-", 
                                        "Distance_stop_ref": "-", "evidence": dico_utr_annot["Evidence"], "UTR_consequence": line[36], "oORF_inframe": line[37], 
                                        "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], "impact": line[13], "Clin_sig": line[32], 
                                        "somatic": line[33], "pheno": line[34], "gnomad_combi": line[23], "gnomad_afr": line[24], "gnomad_ame": line[25], 
                                        "gnomad_jew": line[26], "gnomad_east_asi": line[27], "gnomad_fin": line[28], "gnomad_eur": line[29], 
                                        "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
            else:
                final_table[id].update({"Distance_CDS":dico_utr_annot["DistanceToCDS"], "Distance_stop": dico_utr_annot["DistanceToStop"], "Distance_Cap": dico_utr_annot["CapDistanceToStart"],
                                        "Kozak_context": dico_utr_annot["KozakContext"], "Kozak_strength": dico_utr_annot["KozakStrength"], "type": dico_utr_annot["type"], 
                                        "type_ref": "-", "Distance_stop_ref": "-", "evidence": dico_utr_annot["Evidence"], "UTR_consequence": line[36], 
                                        "oORF_inframe": line[37], "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], "impact": line[13], 
                                        "Clin_sig": line[32], "somatic": line[33], "pheno": line[34], "gnomad_combi": line[23], "gnomad_afr": line[24], 
                                        "gnomad_ame": line[25], "gnomad_jew": line[26], "gnomad_east_asi": line[27], "gnomad_fin": line[28], "gnomad_eur": line[29], 
                                        "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
        else:
            dico_utr_annot = {'ref_StartDistanceToCDS':[], 'alt_type':[], 'ref_type':[], 'alt_type_length':[], 'ref_type_length':[], 'KozakContext':[], 'KozakStrength':[], 'Evidence':[]}
            for info in UTR_annotation:
                info = info.split("=")
                dico_utr_annot[info[0]].append(info[1])
            final_table[id].update({"Distance_CDS":dico_utr_annot["ref_StartDistanceToCDS"], "Distance_stop": dico_utr_annot["alt_type_length"], "Distance_Cap": "-",
                                    "Kozak_context": dico_utr_annot["KozakContext"], "Kozak_strength": dico_utr_annot["KozakStrength"], "type": dico_utr_annot["alt_type"],
                                    "type_ref": dico_utr_annot["ref_type"], "Distance_stop_ref": dico_utr_annot["ref_type_length"], "evidence": dico_utr_annot["Evidence"], 
                                    "UTR_consequence": line[36], "oORF_inframe": line[37], "oORF_outframe": line[38], "uORF": nb_uORF, "consequence": line[6], 
                                    "impact": line[13], "Clin_sig": line[32], "somatic": line[33], "pheno": line[34], "gnomad_combi": line[23], 
                                    "gnomad_afr": line[24], "gnomad_ame": line[25], "gnomad_jew": line[26], "gnomad_east_asi": line[27], "gnomad_fin": line[28], 
                                    "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
    else:
        final_table[id].update({"Distance_CDS":"-", "Distance_stop":"-", "Distance_Cap": "-", "Kozak_context": "-", "Kozak_strength": "-", "type": "-", "evidence": "-", 
                                "UTR_consequence": line[36], "oORF_inframe": line[37], "oORF_outframe": line[38], "uORF": nb_uORF, "type_ref": "-", "Distance_stop_ref": "-", 
                                "consequence": line[6], "impact": line[13], "Clin_sig": line[32], "somatic": line[33], "pheno": line[34], "gnomad_combi": line[23], 
                                "gnomad_afr": line[24], "gnomad_ame": line[25], "gnomad_jew": line[26], "gnomad_east_asi": line[27], "gnomad_fin": line[28], 
                                "gnomad_eur": line[29], "gnomad_other_combi": line[30], "gnomad_south_asi": line[31]})
    return final_table[id]


def MORFAL(gtf, vcf, outdir, hg19, netstart, UTRannotator):
    if os.path.exists(outdir+'/seq2netstart.fasta'):
        os.remove(outdir+'/seq2netstart.fasta')
    
    alphabet = list(string.ascii_letters)

    #################### Creation fichier VCF concatener
    subprocess.run([f"grep -h -v '#' {vcf}/*.vcf | cut -f1-5 | sort -k1,1 -k2,2n | uniq > {outdir}/recup.vcf"], shell=True, capture_output=True, text=True)
    
    with open(f"{outdir}/recup.vcf", "r") as file_in:
        new = []
        for line in file_in:
            line = line.split()
            alt = line[4]
            if "," in alt: 
                alt = alt.split(",")
                for i in range(len(alt)):
                    new.append(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{alt[i]}\n")
            else:
                new.append(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{alt}\n")

    with open(f"{outdir}/recup.vcf", "w") as file_out:
        file_out.write("".join(new))

    subprocess.run([f"sort -k1,1 -k2,2n {outdir}/recup.vcf | uniq > {outdir}/all.vcf"], shell=True, capture_output=True, text=True)

    
    #################### Stockage GTF sous forme de dico
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
                #crea dico qui prend toutes les info de la derniere colonne
                dico_colonne8 = {}
                for col in colonne:
                    col = col.split(" ")
                    dico_colonne8[col[0]] = col[1].replace('"', '')
                
                if line_GTF[2] == 'transcript':
                    #avant de changer de transcrit, creation dico temp qui sera ajoute dans dico gtf avec les infos importantes 
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
                    
                    #pour le nouveau transcrit, recupere les elements importants
                    if (dico_colonne8["transcript_id"] != transcript_id):
                        if (dico_colonne8["gene_id"] != gene_id): 
                            chrom = line_GTF[0]
                            strand = line_GTF[6]
                            gene_name = dico_colonne8["gene_name"]
                            gene_id = dico_colonne8["gene_id"]
                        transcript_id = dico_colonne8["transcript_id"]
                        start_end_transcrit = (int(line_GTF[3]), int(line_GTF[4]))


                #recupere position du codon start 
                elif (line_GTF[2] == "start_codon"):
                    if (dico_colonne8["transcript_id"] == transcript_id):
                        if strand == '+':
                            start_codon = int(line_GTF[3])
                        elif strand == '-':
                            start_codon = int(line_GTF[4])
                
                #recupere debut et fin de chaque exon 
                elif (line_GTF[2] == "exon"):
                    if (dico_colonne8["transcript_id"] == transcript_id):
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

    #################### Correspondance VCF
    with open(f'{outdir}/all.vcf', 'r') as file:
        start_physio_ns = {}
        variant_ns ={}
        final_table = {} 
        for line_VCF in file:
            #recuperation des colonnes d'interet 
            if not line_VCF.startswith('#'):
                line_VCF = line_VCF.split()
                chrom = line_VCF[0]
                pos = int(line_VCF[1])
                ref = line_VCF[3]
                alt = line_VCF[4]

                #recherche dans le dico_gtf en fonction du chr et si le variant est dans le transcrit
                for dico in dico_gtf[chrom]:
                    if dico["start_end_transcript"][0] <= pos <= dico["start_end_transcript"][1]:
                        #verifie si un codon start existe, sinon on ne garde pas le transcrit
                        if dico["start_codon"]:
                            compt = 0
                            #verifie si le variant est dans un exon
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

                                #creation d'un identifiant aleatoire pour chaque transcrit
                                id_var = []
                                while len(id_var) < 15:
                                    id_var.append(random.choice(alphabet))
                                id_var = "".join(id_var)
                                
                                #calcul taille variant (positif si ins et neg si del)
                                if alt == "*" or alt == "-":
                                    len_var = - len(ref)
                                else:
                                    len_var = len(alt) - len(ref)

                                #donne la position du codon start par rapport au debut du transcrit (utile pour netstart)
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

                                    #donne la position du transcrit par rapport au debut du transcrit
                                    if start_seq <= pos < end_seq:
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
                                       
                                    if c < 300 : 
                                        #utilise samtools pour avoir les séquences 
                                        output = subprocess.run([f"samtools faidx {hg19} {chrom}:{start_seq}-{end_seq}"], shell=True, capture_output=True, text=True)
                                        samtools_output = output.stdout
                                        samtools_output = samtools_output.split("\n")
                                        exon_seq_wt = "".join(samtools_output[1:])

                                        #utilise la fonction pour creer la sequence mut
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
                                        else: 
                                            exon_mut = exon_seq_wt
                                        
                                        #si brin antisens, on va utiliser la fonction d'inversion sur seq, alt et ref
                                        if strand == "-":
                                            exon_mut = invers_seq(exon_mut)
                                            exon_seq_wt = invers_seq(exon_seq_wt)
                                        seq_merge_mut.append(exon_mut)
                                        seq_merge_wt.append(exon_seq_wt)
                                if c < 300:         
                                    seq_merge_mut = "".join(seq_merge_mut)
                                    seq_merge_wt = "".join(seq_merge_wt)

                                    #creation du fichier fasta qui servira pour netstart
                                    samtools_output_wt = f">WT:{id_var}"
                                    samtools_output_mut = f">mute:{id_var}"
                                    
                                    output_test = samtools_output_wt.split() + seq_merge_wt.split() + samtools_output_mut.split() + seq_merge_mut.split()
                                    
                                    with open(outdir+'/seq2netstart.fasta', 'a') as output_file:
                                        for item in output_test:
                                            output_file.write(f"{item}\n")

                                    #dictionnaire pour le fichier final 
                                    dico_t = {"chrom": chrom, "pos": pos, "c.":c, "ref": ref, "alt": alt, "gene_name": dico['gene_name'], "gene_id": dico['gene_id'], "transcript_id": dico['transcript_id'], 
                                    "start_transcript": start_transcript, "start_codon": dico['start_codon'], "strand": strand, "variant": variant, "len_variant": len_var}
                                    final_table[id_var] = dico_t
                                    start_physio_ns[id_var] = pos_phy_ns
                                    variant_ns[id_var] = pos_var_ns
    

    #################### Lancement de Netstart et mise de sa sortie dans un dico
    #Netstart ne lit que les fichiers < 490 Ko , il faut donc diviser le fichier fasta en plusieurs petits fasta
    os.chdir(f"{outdir}")
    subprocess.run([f"split -l 112 {outdir}/seq2netstart.fasta split_seq2netstart_"], shell=True, capture_output=True, text=True)

    #pour aller dans le container ou se trouve Netstart
    os.chdir(f"{netstart}")
    #commande pour utiliser Netstart et creer un fichier de sortie
    netstart_result = []
    file = []
    for filename in os.listdir(outdir):
        if filename.startswith("split_seq2netstart_"):
            file.append(filename)

    for fasta in file:
        netstart = subprocess.run([f"singularity run netstart_tcsh.simg -vert {outdir}/{fasta}"], shell=True, capture_output=True, text=True)
        netstart_output = netstart.stdout
        netstart_output = netstart_output.split("\n")
        netstart_result.append(netstart_output)
    
    #dico pour stocker les positions et les scores de chaque atg pour chaque transcrit 
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

        #pour stocker aussi le dernier
        if id : 
            dico_temp = {"pos_ns_wt": pos_ns_wt, "score_wt": score_wt, "pos_ns_mut": pos_ns_mut, "score_mut": score_mut}
            dico_ns[id] = dico_temp
            pos_ns_wt = []
            score_wt = []
            pos_ns_mut = []
            score_mut = []
            id = None

    with open(outdir+"/prediction_positive.txt", "w") as f:
        f.write("info\tID\tpos\tscore\tpred\n")
        for liste in tout:
            if liste[-1] == "Yes":
                lign = "\t".join(liste)
                f.write(f"{lign}\n")

    #################### Delta score et fichier sortie de netstart
    #se servir dico pour faire delta score
    result = []
    score_pos = []
    abo = None
    new = None
    for id_v, dico in dico_ns.items():
        #pour faire toujours apparaitre le codon start dans le fichier final
        for key, val in start_physio_ns.items():
            if key == id_v:
                for wt in dico["pos_ns_wt"]:
                    if wt == val:
                        index_wt = dico["pos_ns_wt"].index(wt)
                        result.append([id_v, str(val), f"site physio {(dico['score_wt'][index_wt])}"])
                        for id_variant, dico_t in final_table.items():
                            if id_variant == id_v:
                                final_table[id_v].update({"pos_physio": val, "score_physio": (dico['score_wt'][index_wt])})
                        if final_table[id_v]["variant"] == "subs":
                            #calcul delta score et ajout au resultat si non nul
                            if wt in dico["pos_ns_mut"]:
                                index_mut = dico["pos_ns_mut"].index(wt)
                                delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                                delta = float("{:.3f}".format(delta))
                                final_table[id_v].update({"delta_start_cds": delta})
                            #verifie si perte d'ATG physio
                            elif wt not in dico["pos_ns_mut"]:
                                delta = "loss"
                                final_table[id_v].update({"delta_start_cds": delta})
                        if final_table[id_v]["variant"] == "del" or final_table[id_v]["variant"] == "ins":
                            if final_table[id_v]["c."] > 0:
                                if wt in dico["pos_ns_mut"]:
                                    index_mut = dico["pos_ns_mut"].index(wt)
                                    delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                                    delta = float("{:.3f}".format(delta))
                                    final_table[id_v].update({"delta_start_cds": delta})
                                #verifie si perte d'ATG physio
                                elif wt not in dico["pos_ns_mut"]:
                                    delta = "loss"
                                    final_table[id_v].update({"delta_start_cds": delta})
                            elif final_table[id_v]["c."] < 0:
                                wt = int(wt) + final_table[id_v]["len_variant"]
                                if wt in dico["pos_ns_mut"]:
                                    index_mut = dico["pos_ns_mut"].index(wt)
                                    delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                                    delta = float("{:.3f}".format(delta))
                                    final_table[id_v].update({"delta_start_cds": delta})
                                #verifie si perte d'ATG physio
                                elif wt not in dico["pos_ns_mut"]:
                                    delta = "loss"
                                    final_table[id_v].update({"delta_start_cds": delta})


        if final_table[id_v]["variant"] == "subs":
            #calcul delta score et ajout au resultat si non nul
            for atg_wt in dico["pos_ns_wt"]:
                if atg_wt in dico["pos_ns_mut"]:
                    index_wt = dico["pos_ns_wt"].index(atg_wt)
                    index_mut = dico["pos_ns_mut"].index(atg_wt)
                    delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                    delta = float("{:.3f}".format(delta))
                    if delta != 0.000:
                        result.append([id_v, str(atg_wt), f"{delta}% ({(dico['score_wt'][index_wt])})"])
                        score_pos.append([delta, ((dico['score_wt'][index_wt]), atg_wt)])
                #verifie si perte d'ATG
                elif atg_wt not in dico["pos_ns_mut"]:
                    index_wt = dico["pos_ns_wt"].index(atg_wt)
                    delta = "abolition"
                    result.append([id_v, str(atg_wt), f"{delta} ({(dico['score_wt'][index_wt])})"])
                    abo =((dico['score_wt'][index_wt]), atg_wt)
            #verifie si creation
            for atg_mut in dico["pos_ns_mut"]:
                if atg_mut not in (dico["pos_ns_wt"]):
                    index_mut = dico["pos_ns_mut"].index(atg_mut)
                    delta = "new"
                    dif_phy_new = final_table[id_v]["score_physio"] - (dico['score_mut'][index_mut])
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
                            "consequence": "-", "impact":"-", "Clin_sig": "-", "somatic": "-", "pheno": "-", "gnomad_combi": "-", 
                            "gnomad_afr": "-", "gnomad_ame": "-", "gnomad_jew": "-", "gnomad_east_asi": "-", "gnomad_fin": "-", 
                            "gnomad_eur": "-", "gnomad_other_combi": "-", "gnomad_south_asi": "-"})
                    
                    abo = None
                    new = None
                    score_pos = []
        
        elif final_table[id_v]["variant"] == "del" or final_table[id_v]["variant"] == "ins":
            for i in range(len(dico["pos_ns_wt"])):
                #pour atg avant variant
                if dico["pos_ns_wt"][i] < variant_ns[id_v]:
                    #calcul delta score et ajout au resultat si non nul
                    wt = dico["pos_ns_wt"][i]
                    if wt in dico["pos_ns_mut"]:
                        index_wt = dico["pos_ns_wt"].index(wt)
                        index_mut = dico["pos_ns_mut"].index(wt)
                        delta = float((((dico["score_mut"][index_mut]) - (dico["score_wt"][index_wt])) / (dico["score_wt"][index_wt])) * 100)
                        delta = float("{:.3f}".format(delta))
                        if delta != 0.000:
                            result.append([id_v, str(wt), f"{delta}% ({(dico['score_wt'][index_wt])})"])
                            score_pos.append([delta, ((dico['score_wt'][index_wt]), wt)])
                    #verifie si perte d'ATG
                    elif wt not in dico["pos_ns_mut"]:
                        index_wt = dico["pos_ns_wt"].index(wt)
                        delta = "abolition"
                        result.append([id_v, str(wt), f"{delta} ({(dico['score_wt'][index_wt])})"])
                        abo =((dico['score_wt'][index_wt]), wt)
                #la meme chose pour atg apres variant
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

            #verifie si creation
            for i in range(len(dico["pos_ns_mut"])):
                if dico["pos_ns_mut"][i] < variant_ns[id_v]:
                    mut = dico["pos_ns_mut"][i]
                    if mut not in (dico["pos_ns_wt"]):
                        index_mut = dico["pos_ns_mut"].index(mut)
                        delta = "new"
                        dif_phy_new = final_table[id_v]["score_physio"] - (dico['score_mut'][index_mut])
                        dif_phy_new = float("{:.3f}".format(dif_phy_new))
                        result.append([id_v, str(mut), f"{delta} ({(dico['score_mut'][index_mut])})"])
                        new = ((dico['score_mut'][index_mut]), mut)
                else:
                    mut = dico["pos_ns_mut"][i]
                    mut_del_ins = dico["pos_ns_mut"][i] - final_table[id_v]["len_variant"]
                    if mut_del_ins not in (dico["pos_ns_wt"]):
                        index_mut = dico["pos_ns_mut"].index(mut)
                        delta = "new"
                        dif_phy_new = final_table[id_v]["score_physio"] - (dico['score_mut'][index_mut])
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
                        
                    #ajoute pour tous les id les keys servant pour UTRannotator
                    final_table[id_v].update({"Distance_CDS":"-", "Distance_stop":"-", "Distance_Cap": "-", "Kozak_context": "-", "Kozak_strength": "-", "type": "-", "evidence": "-", 
                            "UTR_consequence": "-", "oORF_inframe": "-", "oORF_outframe": "-", "uORF": "-", "type_ref": "-", "Distance_stop_ref": "-", 
                            "consequence": "-", "impact":"-", "Clin_sig": "-", "somatic": "-", "pheno": "-", "gnomad_combi": "-", 
                            "gnomad_afr": "-", "gnomad_ame": "-", "gnomad_jew": "-", "gnomad_east_asi": "-", "gnomad_fin": "-", 
                            "gnomad_eur": "-", "gnomad_other_combi": "-", "gnomad_south_asi": "-"})
                    
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


    #################### Utilisation de UTRannotator sur nos variants dans l'UTR et proche du codon start
    os.chdir(f"{UTRannotator}")
    subprocess.run([f"""singularity run VEP_refseq.simg --fork 4 --dir [{UTRannotator}] --fasta /mnt/diag_references/hg19_full_2022_03_24/hg19_order.fa \
                    --cache --offline --merged --format vcf --tab --force_overwrite --af_gnomade --use_transcript_ref --plugin UTRAnnotator,file=/mnt/recherche/uORF/UTRAnnotator/uORF_5UTR_GRCh37_PUBLIC.txt \
                    -i {outdir}/UTR_FirstExons.vcf -o {outdir}/UTR.annotated.vcf"""], shell=True, capture_output=True, text=True)

    

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

                                            

    with open(outdir+"/output_annotation.txt", "w") as out_file:
        out_file.write("##MORFAL output\n"\
        "##GENCODE version gencode 19\n"\
        "##Column descriptions:\n"\
        "#identifiant : indentifier that links the mutation to its sequence in FASTA file\n"\
        "#start_transcript : beginning of the transcript on the chr\n"\
        "#start_CDS : beginning of the protein on the chr\n"\
        "#score_start_CDS(pos_ns) : score predicted by Netstart and position in the transcript\n"\
        "#delta_min : smallest percentage of variation between the reference sequence scores and those of the mutated sequence\n"\
        "#score_delta_min : reference sequence score associated with delta\n"\
        "#delta_max : highest percentage of variation between the reference sequence scores and those of the mutated sequence\n"\
        "#score_delta_max : reference sequence score associated with delta\n"\
        "#abolition(pos): score and position of the lost ATG\n"\
        "#new(pos) : score and position of the gain ATG\n"\
        "#diff_phy_new : difference between the score of the physiological atg and the one created\n"\
        "#evidence : whether the uORF disrupted has any translation evidence\n"\
        "identifiant\tchr\tposition\tc.\tref\talt\tname_gene\tgene_id\ttranscript_id\tstrand\tstart_transcript\tstart_CDS\tscore_ATG_CDS(pos_ns)\t"\
            "delta_ATG_CDS\tdelta_min\tscore_delta_min(pos)\tdelta_max\tscore_delta_max(pos)\tabolition(pos)\tnew(pos)\tdiff_phy_new\tDistance_CDS\t"\
                "Distance_stop\tDistance_Cap_start\tKozak_context\tKozak_strength\ttype\ttype_ref(si changement type)\tDistance_stop_ref(si uFrameshift)\t"\
                    "evidence\tUTR_consequence\tnb oORF_inframe\tnb oORF_outframe\tnb uORF\tconsequence\timpact\tClin_sig\tsomatic\tpheno\tgnomad_combi\tgnomad_afr\t"\
                        "gnomad_ame\tgnomad_jew\tgnomad_east_asi\tgnomad_fin\tgnomad_eur\tgnomad_other_combi\tgnomad_south_asi\n""")

        for id, dico in final_table.items():
            out_file.write(f"{id}\t{dico['chrom']}\t{dico['pos']}\t{dico['c.']}\t{dico['ref']}\t{dico['alt']}\t{dico['gene_name']}\t{dico['gene_id']}\t"\
                f"{dico['transcript_id']}\t{dico['strand']}\t{dico['start_transcript']}\t{dico['start_codon']}\t{dico['score_physio']} ({dico['pos_physio']})\t"\
                    f"{dico['delta_start_cds']}%\t{dico['delta_min(pos)'][0]}%\t{dico['delta_min(pos)'][1][0]} ({dico['delta_min(pos)'][1][1]})\t{dico['delta_max(pos)'][0]}%\t"\
                        f"{dico['delta_max(pos)'][1][0]} ({dico['delta_max(pos)'][1][1]})\t{dico['abolition(pos)'][0]} ({dico['abolition(pos)'][1]})\t"\
                            f"{dico['new(pos)'][0]} ({dico['new(pos)'][1]})\t{dico['diff_phy_new']}\t{dico['Distance_CDS']}\t{dico['Distance_stop']}\t"\
                                f"{dico['Distance_Cap']}\t{dico['Kozak_context']}\t{dico['Kozak_strength']}\t{dico['type']}\t{dico['type_ref']}\t{dico['Distance_stop_ref']}\t"\
                                    f"{dico['evidence']}\t{dico['UTR_consequence']}\t{dico['oORF_inframe']}\t{dico['oORF_outframe']}\t{dico['uORF']}\t{dico['consequence']}\t{dico['impact']}\t"\
                                        f"{dico['Clin_sig']}\t{dico['somatic']}\t{dico['pheno']}\t{dico['gnomad_combi']}\t{dico['gnomad_afr']}\t{dico['gnomad_ame']}\t{dico['gnomad_jew']}\t"\
                                            f"{dico['gnomad_east_asi']}\t{dico['gnomad_fin']}\t{dico['gnomad_eur']}\t{dico['gnomad_other_combi']}\t{dico['gnomad_south_asi']}\n")
    
    subprocess.run([f"rm {outdir}/split_seq2netstart_*"], shell=True, text=True)



def main():
    start_time = time.time()

    # get arguments
    args = get_arguments()
    gtf = os.path.realpath(args.gtf)
    vcf = os.path.realpath(args.vcf)
    outdir = os.path.realpath(args.outdir)
    hg19 = os.path.realpath(args.hg19fasta)
    netstart = os.path.realpath(args.netstart_container)
    UTRannotator = os.path.realpath(args.UTRannotator)


    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # MORFAL
    MORFAL(gtf=gtf, vcf=vcf, outdir=outdir)

    step = f"--- Program {os.path.basename(__file__)} executed in {round((float(time.time() - start_time) / 60), 2)} minutes ---"; print(step)



if __name__ == "__main__":
    main()