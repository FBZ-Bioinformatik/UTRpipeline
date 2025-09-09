###### COMBINED PIPELINE ##########
# 1. I have two different pipelines: new_vcf_referece_38.py (hg38) and expanded_pipeline_37.py (hg19)
# 2. They are both nearly identical but one processes hg19 files and the other hg38 files
# GOAL: One pipeline that can
# 1. Use either hg19 or hg38
# 2. Load the correct reference genome (hs37d5.fa or GRCh38.fa)
# 3. Parse the matching MANE annotation file (reference_MANE_hg19.txt or reference_MANE_hg38.txt)
# 4. Use the right VCF input (vcf_37.vcf.gz or vcf_38.vcf.gz)
# 5. Everything else remains the same (run all the same logic)
# 5.1 Identify variants in the 5' untranslated regions of MANE transcripts
# 5.2 Check if any of these variants introduce a new start codon within the 5'UTR
# 5.3 Extract the surrounding sequence (Kozak sequence) and evaluate the Kozak consensus strength of the introduced startcodon


import pandas as pd
import pysam 
import os
import argparse

col_names = [
    "bin", "name", "chrom", "strand", "txStart", "txEnd",
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
    "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"
]

parser = argparse.ArgumentParser()
#parser.add_argument("--build", choices=["hg19","hg38"], required =True)
parser.add_argument("--vcf", required=True, help="Optional: Path to a costum VCF file")
parser.add_argument("--transcript", required=True, help="Optional: Path to a custom MANE annotation file")
parser.add_argument("--ref", required = True, help="Optional: Path to referece genome file")

args = parser.parse_args()

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.abspath(os.path.join(script_dir, ".."))

vcf_path = os.path.join(base_dir, args.vcf)
mane_path = os.path.join(base_dir, args.transcript)

ref_genome = os.path.join(base_dir, args.ref)


#print(f"Using build: {genome_build}")
print(f"Reference genome: {ref_genome}")
print(f"MANE file: {mane_path}")
print(f"VCF file: {vcf_path}")

######

print("Loading reference genome...")
fasta = pysam.FastaFile(ref_genome)


mane_df = pd.read_csv(mane_path, sep="\t", header=None, names=col_names, skiprows=1)
#print(mane_df.head())

#mane_df.to_excel("/mnt/c/Users/maria/OneDrive/Documents/Master Arbeit Brust- und Eierstockkrebs/MANE_38.xlsx", index =False)

############## 
# Translation table 
tab = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement(seq):
     return seq.translate(tab)[::-1]     

# Constructs the full transcript sequence (mRNA) from exon coordinates and
# returns the full transcript sequence (reverse complemented if the exon is on the - strand)
def build_transcript(fasta, chromosome, exon_starts, exon_ends, strand):
    
    # Function concatenates all exonic sequences inorder
    transcript_seq = ""
    for start, end in zip(exon_starts, exon_ends):
        seq = fasta.fetch(chromosome, start, end)
        transcript_seq += seq.upper()

    # The sequence is reversed complemented for variants in the negatice strand
    #if strand == "-":
        #transcript_seq = reverse_complement(transcript_seq)
    
    return transcript_seq

# if the genomic position falls inside an exon, it calculates where it lands within the transcript
# at which position in the transcript does the genomic coordinate fall?

def map_genom_to_transcript(genomic_pos, exon_starts, exon_ends, strand):
    transcript_index = 0 # counts the positions within the transcript (not genome)

    
    for start, end in zip(exon_starts, exon_ends):
        # Check if the genomic position lies within this exon
        if start <= genomic_pos < end:
            return transcript_index + (genomic_pos - start) # position in transcript wo die variante ist ausgeben!!
        
        transcript_index += (end - start)

    return None 

# transcript_index = position in the mRNA sequence (concatenated transcript), z.B. 197
# function maps it back to the corresponding genomic coordinate based on exon boundaries
# genomic_index = tracks how far along the transcript
# exon_length= how many transcript bases come from this exon
# does the exon contain the transcript index?
# yes? compute offset --> substract how many transcript bases we already walked from the target index
# genomic positon = start + offset


# Simulate translation of an ORF starting from the position where the startcodon is found in the transcript
# Takes the full transcript sequence (5' to 3') and the startcodon position where the ORF starts
# Returns a dictionary with the ORF sequences and it's length

def analyze_orf_transcript(transcript_seq, startcodon_index):
    stop_codons = {"TAA", "TAG", "TGA"}
    orf_seq = ""
    
    # Scans the transcript 3 bases at a time starting from a start codon
    for i in range(startcodon_index, len(transcript_seq), 3):
        codon = transcript_seq[i:i+3]
        if len(codon) < 3:
            break
        orf_seq += codon
        if codon in stop_codons:
            break
    return {
        "orf_seq": orf_seq,
        "orf_length": len(orf_seq)
    }

# Determines which exon a the startcodon falls into
def find_exon_number(genomic_pos, exon_starts, exon_ends, strand):
    if strand == "+":
        for i, (start, end) in enumerate(zip(exon_starts, exon_ends), start =1):
            if start <= genomic_pos < end:
                return i
    if strand =="-":
        for i, (start, end) in enumerate(zip(exon_starts[::-1], exon_ends[::-1]), start =1):
            if start <= genomic_pos < end:
                return i
            
        return None

#### dynamic gene_info dictionary!
gene_info = {}

for index, row in mane_df.iterrows():
        
        transcript_id= row['name']
        #print('Transcript ID:', row['name'])
        chromosome = row['chrom'].replace("chr","")
        #print('Chromosome:', row['chrom'])
        strand = row['strand']
        #print('Strand:', row['strand'])
        tx_start = row['txStart']
        tx_end = row['txEnd']
        gene_name = row['name2']
        #print('Gene name:', row['name2'])
        cds_start = row['cdsStart']
        cds_end = row['cdsEnd']
        exon_starts = [int (x) for x in row['exonStarts'].strip(',').split(',')]
        exon_ends = [int(x) for x in row['exonEnds'].strip(',').split(',')]

        if strand == "+":
            zipped_exons = zip(exon_starts, exon_ends)
            utr_start = tx_start
            utr_end = cds_start
        else:
            zipped_exons = zip(exon_starts[::-1], exon_ends[::-1])
            utr_start = cds_end
            utr_end = tx_end

        exonic_regions = []    

        for i, (start, end) in enumerate(zipped_exons,start=1):
            if end <= utr_start or start >= utr_end:
                continue 

            label = "unknown"
            if start >= utr_start and end <= utr_end:
                label = "full"
                #print(f"Exon {i} is part of the 5'UTR")
            elif start < utr_start and end > utr_start:
                label = "left_overlap"
                #print(f"Exon {i} is partially part of the 5'UTR (overlap)")
            elif start < utr_end and end > utr_end:
                label = "right_overlap"
                #print(f"Exon {i} is partially part of the 5'UTR (overlap)")
            else:
                label = "internal split"
                #print(f"Exon {i} is partially part of the 5'UTR (split")
            
            region = {
            "start": max(start, utr_start),
            "end": min(end, utr_end),
            "exon": i,
            "label": label
            }
            exonic_regions.append(region)

        if exonic_regions:
            if chromosome not in gene_info:
                gene_info[chromosome] = []

            gene_info[chromosome].append({
                "exonic_regions": exonic_regions,
                "strand": strand,
                "gene_name": gene_name,
                "transcript_id": transcript_id
            })

def correct_region(vcf_record, exonic_regions):
    for region in exonic_regions:
        if region["start"] <= vcf_record.pos <= region["end"]:
            return True
    return False 

# 1. 
################# QUALITY FILTER 1: DOES THE REFERENCE MATCH? 5'UTR REGIO? ######################
def compare_vcf_to_ref(vcf_path, fasta, gene_info):
    mismatch_found = False 

    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            chrom = record.chrom
            pos = record.pos
            #print(pos)
            ref = record.ref
            #print(ref)
            alt = record.alts
            #print(alt)

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            if chrom not in gene_info:
                continue

            ref_seq = fasta.fetch(chrom, pos-1, pos-1 + len(ref))
        
            if str(ref_seq) != ref:
                mismatch_found = True
                print(f"Mismatch at {chrom}:{pos} | VCF REF: {ref} | Genome REF: {ref_seq}")
            else:
                if not mismatch_found:
                    pass
            
            
            for entry in gene_info[chrom]:
                exonic_regions = entry['exonic_regions']
                strand = entry['strand']
                gene = entry['gene_name']
                transcript = entry['transcript_id']

            # Process only the variants in the 5'UTR
                print(f"Checking variant {chrom}:{pos} against {gene} ({transcript})....")
                if correct_region(record, exonic_regions):
                    transcript_row = mane_df[mane_df['name'] == transcript].iloc[0]
                    exon_starts = [int(x) for x in transcript_row['exonStarts'].strip(',').split(',')]
                    exon_ends = [int(x) for x in transcript_row['exonEnds'].strip(',').split(',')]
                    cds_start = int(transcript_row["cdsStart"])  # make sure it's an int
                    cds_end = int(transcript_row["cdsEnd"])
                    transcript_seq = build_transcript(fasta, chrom, exon_starts, exon_ends, strand)

        # Now call with the correct data
                    process_vcf(fasta, chrom, pos, ref, alt, strand, exon_starts, exon_ends, transcript_seq, cds_start, cds_end, gene)
                    #print(f"Variant at {chrom}:{pos} in {gene} is in the 5'UTR!")
                    
                
                else:
                #print(f"Variant at {chrom}:{pos} is not in the 5'UTR region")
                    pass
    


# 2. 

################# VCF FILE PROCESSING: ADD FLANKING BASES TO FIND START CODON ################

# Generates the reverse complement of the DNA sequence
# Reverses the input sequence
# Adjust the sequences based on the strand orientation
# to match the strand where the variant resides
# For reverse strand the sequence are processsed differently


def process_vcf(fasta, chrom, pos, ref, alt, strand, exon_starts, exon_ends, transcript_seq, cds_start, cds_end, gene):

    # Extended sequences for REF and ALT

    extended_ref = extraq_flanking_bases(fasta, chrom, pos, ref, strand)
    extended_alts = alt_with_flanking_bases(fasta, chrom, pos, alt, ref, strand)
    print(f"Processing variant at {pos} in chromosome {chrom} in {gene}")
    print(f"The variant is in the strand: {strand}")
    print(f"Extended REF sequence: {extended_ref}")

    for extended_alt in extended_alts:
        print(f"Extended ALT sequence: {extended_alt}")

    start_codon = "ATG" 
    
    for extended_alt in extended_alts:

        if start_codon in extended_alt:
            print(f"Start codon {start_codon} created at position {pos} in chromosome {chrom}")
            #print(extended_alt)
            startcodon_index = extended_alt.find(start_codon)
            #print(startcodon_index)
            
            genomic_start_codon_pos = pos - 3 + startcodon_index 
            
            if start_codon in extended_alt:
                    for alt_allele in alt:
                        if strand == "+":
                            wider_ref = fasta.fetch(chrom, pos -6, pos + len(ref)+2)
                            wider_alt = wider_ref[:5] + alt_allele + wider_ref[-3:]

                        if strand =="-":
                            wider_ref = fasta.fetch(chrom, pos-4, pos+len(ref)+4)
                            wider_alt = wider_ref[:3] + alt_allele + wider_ref[-5:]
                            wider_alt = reverse_complement(wider_alt)
                        
                        index = wider_alt[3:].find('ATG') + 3
                        #print(wider_alt)
                        
                        print("KOZAK1", wider_alt[index-3:index], "KOZAK2", wider_alt[index+3:index+4])
                        kozak_seq = wider_alt[index-3:index+4]
                        print(f"This is the Kozak sequence: {kozak_seq}")
                        strength = evaluate_kozak_strength(kozak_seq, strand)
                        print(f"This is the strength: {strength}")
                        #kozak_index = transcript_seq.find(kozak_seq)
                        #print(kozak_index)

                        transcript_seq_plus = build_transcript(fasta, chrom, exon_starts, exon_ends, "+") 

                        tx_index = map_genom_to_transcript(pos-1, exon_starts, exon_ends, strand)
                        print(f"this is the tx_index {tx_index} in strand {strand}")
                        mutant_tx = transcript_seq_plus[:tx_index] + alt_allele + transcript_seq_plus[tx_index+len(ref):]
                        
                        if strand =="-":
                            mutant_tx = reverse_complement(mutant_tx)
                            atg_index = mutant_tx[3:].find("ATG")
                            print(f"this is the ATG index {atg_index}")
                            orf_start_index = atg_index + 3

                        else:
                            search_start = max(tx_index - 6, 0)
                            atg_index = mutant_tx[3:].find("ATG", search_start)
                            orf_start_index = atg_index 
                
                        print(f"ORF start {orf_start_index} in strand {strand}")
            
                        orf_result = analyze_orf_transcript(mutant_tx, orf_start_index)
                        print(f"ORF starts at this index {orf_start_index}")
                        print(f"ORF length: {orf_result['orf_length']}")
                        print(f"ORF sequence: {orf_result['orf_seq']}")
                    
                        orf_end_index = orf_start_index + orf_result["orf_length"]
                        print(f"This is the index of the ORF´s end: {orf_end_index}")

                        # get the genomic coordinate of the CDS in the sequence with the start codon
                        cds_start_genomic = cds_start if strand == "+" else cds_end -1
                        # map the transcript space using exon_starts/exon_ends
                        cds_start_index = map_genom_to_transcript(cds_start_genomic, exon_starts, exon_ends, strand)
                        print(f"This is the CDS start index of the UCSC file {cds_start_index}")

                        # if the variant is in the negative strand, reverse-complement the transcript and adjust the index
                        if strand =="-":
                            mutant_tx = reverse_complement(mutant_tx)
                            cds_start_index = len(mutant_tx) - cds_start_index -1
                            print(f"This is the reverse-complemented CDS start index for the negative variants {cds_start_index}")
                        
                        if (cds_start_index - orf_start_index) % 3 == 0:
                            print("New start codon is in-frame with the CDS")
                        else:
                            print("New start codon is out-of-frame with the CDS")

                        # if (len(alt_allele) - len(ref)) % 3 != 0:
                        #     print("Resulting uORF is out-of-frame")
                        # else:
                        #     print("Resulting uORF is in-frame")

                        if orf_end_index >= cds_start_index:
                            print("ORF overlaps with the CDS")
                        else:
                            print("ORF remains in the 5'UTR")
                        
                        overlaps_cds = (orf_end_index > cds_start_index)
                        out_of_frame_vs_cds = ((cds_start_index - orf_start_index) % 3 != 0)
                        if overlaps_cds and out_of_frame_vs_cds:
                            print("Out-of-Frame oORF likely to affect native protein translation")
                        elif overlaps_cds and not out_of_frame_vs_cds:
                            print("In-frame oORF, N-terminal extension (Not CDS frameshift)")
                        else:
                            print("Non-overlapping uORF, does not affect CDS")


                        exon_number = find_exon_number(genomic_start_codon_pos, exon_starts, exon_ends, strand)
                        print(f"Start codon falls in exon {exon_number}")
            
def extraq_flanking_bases(fasta, chromosome, pos, ref, strand): 
    ref_length = len(ref)
    extended_ref = fasta.fetch(chromosome, pos-3, pos+ref_length+1)

    return extended_ref

def alt_with_flanking_bases(fasta, chrom, pos, alt, ref, strand):
    # All the alt alleles are stored in the list
    alt_sequences = []

    ref_seq = fasta.fetch(chrom, pos-3, pos+len(ref)+1)

    for alt_allele in alt:
        alt_seq = ref_seq[:2] + alt_allele + ref_seq[-2:] 
        #print(f"Printing extended ALT: {alt_seq}")
        if strand == "-":
            alt_seq = reverse_complement(alt_seq)
            #print(f"Printing reversed extended ALT:{alt_seq}")
        alt_sequences.append(alt_seq)


    return alt_sequences

####### 3. KOZAK SEQUENCE ####################
# ANN AUG G --> strong
# GNNAUGG --> strong
# ANNAUGN --> moderate
# GNNAUGN --> moderate
# NNNAUGG --> moderate
# NNNAUGN --> weak  

def evaluate_kozak_strength(kozak_seq, strand):
    
    if len(kozak_seq) != 7:
        return "invalid"  # Ensure the sequence is exactly 7 bases

    position_1 = kozak_seq[0]  # First position (1st base)
    start_codon = kozak_seq[3:6]
    position_7 = kozak_seq[6]  # Last position (7th base)

    if start_codon != "ATG":
        return "invalid"
    
    # Check for strong Kozak sequence
    if position_1 in ('A', 'G') and position_7 == 'G':
        return "strong"

    # Check for moderate Kozak sequence
    elif position_1 in ('A', 'G') or position_7 == 'G':
        return "moderate"

    # If neither strong nor moderate, it's weak
    else:
        return "weak"


print(f"Processing VCF file ... {vcf_path}")
vcf_7 = compare_vcf_to_ref(vcf_path, fasta, gene_info)
