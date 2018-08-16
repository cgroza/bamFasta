import pysam
import sys

fasta_path = sys.argv[1]
bam_path = sys.argv[2]
new_fasta_path = sys.argv[3]

bam = pysam.AlignmentFile(bamPath, "rb")
fasta = pysam.FastaFile(fastaPath)

#NOTE: BAM, VCF coordinates are 1 based. Make sure to adjust when substituting into 0 based strings
#NOTE: Maybe integrate some collision detection between reads?

def log(qname, rname, s, e, o):
    """Logging function to simplify output to stderr"""
    sys.stderr.write(qname + " in " + rname + ": " + str(s) + " - " + str(e) + " offset = " + str(o) + "\n")

def individualize_contig(contig_name):
    ref_seq = fasta.fetch(reference=contig_name)
    # initially, no offset because no substitutions yet
    offset = 0
    for read in bam.fetch(contig=contig_name):
        # not needed
        #chrom = read.reference_name
        query_seq = read.query_alignment_sequence
        # length spanned by the alignment in the reference sequence
        aligned_length = read.reference_length
        # length of read
        read_length = read.infer_query_length()
        start = read.reference_start + offset
        end = read.reference_end + offset
        # make the substitution in the reference sequence
        ref_seq = ref_seq[:start] + query_seq + ref_seq[end:]
        log(read.query_name, contig_name, start, end, offset)
        # update the offset of the contig coordinate
        offset = (read_length - aligned_length) + offset
    return ref_seq

out = open(newFastaPath 'w')
for contig in fasta.references:
    new_seq = individualize_contig(contig)
    out.write(">" + contig + "\n")
    out.write(out)

bam.close()
fasta.close()
