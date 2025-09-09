# UTRpipeline

Pipeline was created to detect and characterize variants in the 5' untranslated regions (5'UTRs) of MANE transcripts.
The goal is to identify variants that create novel start codons in 5'UTRs, extract the Kozak sequence, simulate translation of the resulting uORF, and evaluate whether the new ORF overlaps with or disrupts the coding sequence (CDS), and it differentiates between in-frame relative to the CDS (possible N-terminal extension), out-of-frame (potential frameshift effect), or non-overlapping (remains in UTR without affecting CDS).
The pipeline works with hg38 (GRCh38)

The pipeline dynamically loads:
- The matching transcript annotation file (`GRCh38.fa`)
- The correct reference genome (`reference_MANE_hg38.txt`)
- The appropriate VCF input (`vcf_38.vcf.gz`)

# Requirements

- Python 3.8+  
- [pandas](https://pandas.pydata.org/)  
- [pysam](https://pysam.readthedocs.io/)

The pipeline can be ran using the commands: `python UTRpipeline.py \
    --vcf path/to/vcf_38.vcf.gz \
    --transcript path/to/reference_MANE_hg38.txt \
    --ref path/to/GRCh38.fa` 
