#!/usr/bin/python
import subprocess
import argparse
import vcf
import os
import shutil

def parse_svaba(input_vcf, SDID, output, vcftype):
    """
    vawk '{print $1, $2, $2+1, "P-00356971_svaba", "BND", $5, S$*$AD, S$*$DP}'
    """
    header = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_normal\tSUPPORT_tumor\tDPnormal\tDPtumor NOTES\tGENES\"" + " > "+ output +  "_" + vcftype + "_svaba.mut"
    svaba_cmd = "vawk '{print $1, $2, $2+1, \"" + SDID + "_svaba\", \"BND\", $5, S$*$AD, S$*$DP}' " + \
                " " + input_vcf + " >> " + output + "_" + vcftype + "_svaba.mut"

    subprocess.call(" && ".join([header, svaba_cmd]), shell=True)


def parse_svict(input_vcf, SDID, output):
    """
    CHROM START   END SDID    SVTYPE  ALT SUPPORTING_READS    NOTES   GENES
    vawk '{if (I$SUPPORT>8 && I$SVTYPE ~ "INS" ) print $1, $2, I$END, "$SDID_svict", I$SVTYPE, "<INS>", I$SUPPORT 
        else if (I$SUPPORT>8 && I$SVTYPE ~ "INV" ) print $1, $2, I$END, "P-00356971_svict", I$SVTYPE, "<INV>", I$SUPPORT
        else if (I$SUPPORT>8 && I$SVTYPE ~ "DEL" ) print $1, $2, I$END, "P-00356971_svict", I$SVTYPE, "<DEL>", I$SUPPORT
        else if (I$SUPPORT>8 && I$SVTYPE ~ "BND" ) print $1, $2, $2+1, "$SDID_svict", I$SVTYPE, $5, I$SUPPORT}' $input.vcf
    """
    header8 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\tNOTES\tGENES\"" + " > "+ output + "_svict_SR8.mut"
    header12 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\tNOTES\tGENES\"" + " > "+ output + "_svict_SR12.mut"

    sup_8 = "vawk '{if (I$SUPPORT>8 && I$SVTYPE ~ \"INS\" ) print $1, $2, I$END, \"" + SDID + "_svict\", I$SVTYPE, \"<INS>\", I$SUPPORT ;" + \
        " else if (I$SUPPORT>8 && I$SVTYPE ~ \"INV\" ) print $1, $2, I$END, \"" + SDID + "_svict\", I$SVTYPE, \"<INV>\", I$SUPPORT ; " + \
        " else if (I$SUPPORT>8 && I$SVTYPE ~ \"DEL\" ) print $1, $2, I$END, \"" + SDID + "_svict\", I$SVTYPE, \"<DEL>\", I$SUPPORT ; " + \
        " else if (I$SUPPORT>8 && I$SVTYPE ~ \"BND\" ) print $1, $2, $2+1, \"" + SDID + "_svict\", I$SVTYPE, $5, I$SUPPORT}' " + input_vcf + \
        " >> " + output + "_svict_SR8.mut"

    sup_12 = "vawk '{if (I$SUPPORT>12 && I$SVTYPE ~ \"INS\" ) print $1, $2, I$END, \"" + SDID + "_svict\", I$SVTYPE, \"<INS>\", I$SUPPORT ; " + \
        " else if (I$SUPPORT>12 && I$SVTYPE ~ \"INV\" ) print $1, $2, I$END, \"" + SDID + "_svict\", I$SVTYPE, \"<INV>\", I$SUPPORT ; " + \
        " else if (I$SUPPORT>12 && I$SVTYPE ~ \"DEL\" ) print $1, $2, I$END, \"" + SDID + "_svict\", I$SVTYPE, \"<DEL>\", I$SUPPORT ;" + \
        " else if (I$SUPPORT>12 && I$SVTYPE ~ \"BND\" ) print $1, $2, $2+1, \"" + SDID + "_svict\", I$SVTYPE, $5, I$SUPPORT}' " + input_vcf + \
        " >> " + output + "_svict_SR12.mut"
    
    subprocess.call(" && ".join([header8, header12, sup_8, sup_12]), shell=True)

def parse_lumpy(input_vcf, SDID, output):
    """
    ##CHROM START   END SDID    SVTYPE  ALT SVLENGTH    SUPPORT_READS    NOTES   GENES
    vawk '{if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != "hs37d5" && $1 !~ "GL" && $5 !~ "hs37d5" && I$SU>50 && I$SVTYPE ~ "BND") print $1, $2, $2+1, "P-00356971_lumpy", I$SVTYPE, $5, "NA", I$SU
         else if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != "hs37d5" && $1 !~ "GL" && $5 !~ "hs37d5" && I$SU>50 && I$SVTYPE !~ "BND") print $1, $2, I$END, "P-00356971_lumpy", I$SVTYPE, $5, I$SVLEN, I$SU}'
    """
    header1k_sup_50 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSVLENGTH\tSUPPORT_READS\tNOTES\tGENES\"" + " > "+ output + "_lumpy_len1k_SU50.mut"
    header500_sup_24 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSVLENGTH\tSUPPORT_READS\tNOTES\tGENES\"" + " > "+ output + "_lumpy_len500_SU24.mut"

    len1k_sup_50 = "vawk '{if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>50 && I$SVTYPE ~ \"BND\") print $1, $2, $2+1, \""+ SDID +"_lumpy\", I$SVTYPE, $5, \"NA\", I$SU ;"  + \
                " else if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>50 && I$SVTYPE !~ \"BND\") print $1, $2, I$END, \""+ SDID +"_lumpy\", I$SVTYPE, $5, I$SVLEN, I$SU}' " + input_vcf + \
                " >> " + output + "_lumpy_len1k_SU50.mut"

    len500_sup_24 = "vawk '{if ((I$SVLEN>500 || I$SVLEN<-500) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>24 && I$SVTYPE ~ \"BND\") print $1, $2, $2+1, \"" + SDID + "_lumpy\", I$SVTYPE, $5, I$SVLEN, I$SU ;" + \
                " else if ((I$SVLEN>500 || I$SVLEN<-500) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>24 && I$SVTYPE !~ \"BND\") print $1, $2, I$END, \"" + SDID + "_lumpy\", I$SVTYPE, $5, I$SVLEN, I$SU }' " + input_vcf + \
                " >> " + output + "_lumpy_len500_SU24.mut"

    subprocess.call(" && ".join([header1k_sup_50, header500_sup_24, len1k_sup_50, len500_sup_24]), shell=True)


def create_symlink(travers_dir_name, src_dir, igvnav_dirname_dst, suffix):
    "Recursively Traverse through the directory and create symlink"
    for root, dirs, files in os.walk(travers_dir_name):
        for each_file in files:
            if each_file.endswith(suffix) and not os.path.exists(os.path.join(igvnav_dirname_dst,each_file)):
                os.symlink(os.path.join(root,each_file), os.path.join(igvnav_dirname_dst,each_file))
    return



parser = argparse.ArgumentParser(description=
    'A MUT file (.mut) is a tab-delimited text file that lists mutations. \
    The first row contains column headings and each subsequent row identifies a mutation. \
    IGV ignores the column headings.It reads the first five columns as shown below and \
    ignores all subsequent columns:  \
    1. chromosome \
    2. start location (location of the first base pair in the mutated region) \
    3. end location (location of the last base pair in the mutated region) \
    4. sample or patient ID \
    5. mutation type (for example, Synonymous, Missense, Nonsense, Indel, etc.)')
parser.add_argument('--vcf', required=True, help="Input VCF file")
parser.add_argument('--sdid', required=True, help="SDID from analysis")
parser.add_argument('--vcftype', help="somatic (or) germline vcf (only for svaba)")
parser.add_argument('--tool', required=True, help="Tool name - Variant callers")
parser.add_argument('--output', required=True,
                    help="output tab demilited file for IGVNav, format=output.mut")
args = parser.parse_args()

vcftype = args.vcftype
vcf = args.vcf
sv_caller = args.tool
sdid = args.sdid
output = args.output

if sv_caller == 'svict':
    parse_svict(vcf, sdid, output)
elif sv_caller == 'lumpy':
    parse_lumpy(vcf, sdid, output)
elif sv_caller == 'svaba':
    parse_svaba(vcf, sdid, output, vcftype)

igvnav_dirname_dst = os.path.join(os.path.dirname(os.path.abspath(args.output)), 'IGVnav')
src_dir = os.path.abspath(os.path.dirname(os.path.abspath(args.output)))

try:
    if not os.path.exists(igvnav_dirname_dst): os.mkdir(igvnav_dirname_dst)
    for each_input in [('bams','-nodups.bam'), ('bams','.overlapped.bam'), ('variants','.vep.vcf'), ('svs/igv','.mut'), ('svs','.gtf'),('svs','.bam'),('svs/svaba', '.contigs.bam')]:
        dir_name = os.path.join(src_dir,each_input[0])
        create_symlink(dir_name, src_dir, igvnav_dirname_dst, each_input[1])
except Exception as e:
    print(e)


