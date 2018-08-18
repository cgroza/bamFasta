import pysam
import sys

if len(sys.argv) != 4:
    print "Usage: bamFasta in.fa in.bam out.fa"
    exit(1)

fasta_path = sys.argv[1]
bam_path = sys.argv[2]
new_fasta_path = sys.argv[3]

bam = pysam.AlignmentFile(bam_path, "rb")
fasta = pysam.FastaFile(fasta_path)

n_seqs = 0
n_bases = 0

#NOTE: SAM, VCF coordinates are 1 based. But pysam converts all coordinates to 0-based.
#NOTE: Maybe integrate some collision detection between reads?

def log(qname, rname, s, e, o, d):
    """Logging function to simplify output to stderr"""
    sys.stderr.write(qname + " in " + rname + ": " + str(s) + " - " + str(e) + " offset = " + str(o) + " delta_offset= " +  str(d) + "\n")

def individualize_contig(contig_name, overlapped = None):
    global n_seqs, n_bases
    ref_seq = fasta.fetch(reference=contig_name)
    # initially, no offset because no substitutions yet
    offset = 0
    for read in bam.fetch(contig=contig_name):
        if read in overlapped:
            continue
        # not needed
        #chrom = read.reference_name
        query_seq = read.query_alignment_sequence
        # print query_seq
        # length spanned by the alignment in the reference sequence
        aligned_length = read.reference_length
        # length of read
        read_length = read.infer_query_length()
        start = read.reference_start + offset
        end = read.reference_end + offset
        # make the substitution in the reference sequence.
        # offset the end by 1 base since it belongs to the read. start is not included in the prefix by default
        ref_seq = ref_seq[:start] + query_seq + ref_seq[end:]
        # update the offset of the contig coordinate
        delta_offset = read_length - aligned_length
        offset = delta_offset + offset
        log(read.query_name, contig_name, start, end, offset, delta_offset)
        n_seqs += 1
        n_bases += aligned_length
    return ref_seq + "\n"

def report_overlaps():
    """This function is intended to detect overlapping reads for later exclusion."""
    overlaps = []
    for read1 in bam.fetch():
        for read2 in bam.fetch(contig=read1.reference_name):
            # Do not compare read with self
            if read1 == read2:
                continue
            # Check if read2 start or end is within read1 interval
            if (read1.reference_start < read2.reference_start and read2.reference_start < read1.reference_end) or (read1.reference_start < read2.reference_end  and read2.reference_end < read1.reference_end):
                overlaps.append(read2)
    return set(overlaps)


overlapped_reads = report_overlaps()
print "Total reads: " +  str(bam.count())
print "Overlapping reads: " + str(len(overlapped_reads))


out = open(new_fasta_path, 'w')
for contig in fasta.references:
    new_seq = individualize_contig(contig, overlapped_reads)
    out.write(">" + contig + "\n")
    out.write(new_seq)

print "Substituted sequences: " + str(n_seqs)
print "Substituted reference bases: " + str(n_bases)

bam.close()
fasta.close()
out.close()

