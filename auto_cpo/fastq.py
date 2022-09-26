import gzip
import io



def get_first_n_reads(fastq_path, num_reads):
    num_lines = num_reads * 4
    lines_read = 0
    num_reads_read = 0
    reads = []
    if fastq_path.endswith('.gz'):
        f = io.TextIOWrapper(io.BufferedReader(gzip.open(fastq_path, 'rb')))
    else:
        f = open(fastq_path, 'r')

    while num_reads_read < num_reads:
        header = next(f).strip()
        seq = next(f).strip()
        plus = next(f).strip()
        quality = next(f).strip()
        read = {
            'header': header,
            'seq': seq,
            'quality': quality,
        }
        reads.append(read)
        num_reads_read += 1 
    f.close()

    return reads


def estimate_read_length(reads):
    num_reads = len(reads)
    total_len_seqs = 0
    avg_read_len = 0
    estimated_read_len = 150

    for read in reads:
        total_len_seqs += len(read['seq'])
    if num_reads > 0:
        avg_read_len = total_len_seqs / num_reads

    if avg_read_len > 50 and avg_read_len < 100:
        estimated_read_len = 100
    elif avg_read_len > 100 and avg_read_len < 150:
        estimated_read_len = 150
    elif avg_read_len > 150 and avg_read_len < 200:
        estimated_read_len = 200
    elif avg_read_len > 150 and avg_read_len < 200:
        estimated_read_len = 250

    return estimated_read_len
        
