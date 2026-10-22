import os
import csv
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def gen_file_path(input_file):
    """
    """
    result_dir = os.path.dirname(input_file)
    basename = os.path.splitext(os.path.basename(input_file))[0]
    oligominer_fastq = os.path.join(result_dir, basename) + '.fastq'
    oligominer_fasta = os.path.join(result_dir, basename) + '_oligominer.fasta'
    oligo_sets = os.path.join(result_dir, basename) + '_oligosets.csv'
    param_file = os.path.join(result_dir, basename) + '_parameters.csv'
    blast_result = os.path.join(result_dir, basename) + '_blast.csv'
    return oligominer_fastq, oligominer_fasta, oligo_sets, param_file, blast_result


def record_param(mFISH3D_param, oligominer_param, param_file):
    """
    """
    f = open(param_file, 'w')
    w = csv.writer(f)
    for i in sorted(mFISH3D_param, key=str.lower):
        w.writerow([i, mFISH3D_param[i]])
    for i in sorted(oligominer_param, key=str.lower):
        w.writerow([i, oligominer_param[i]])
    f.close()


def read_blast(blast_result):
    """
    """
    df_blast = pd.read_csv(blast_result,
                           names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                  'sstart', 'send', 'evalue', 'bitscore'])
    return df_blast


def convert_fastq2fasta(fastq_file, fasta_file):
    """
    convert Fastq file to fasta format.
    """
    count = SeqIO.convert(fastq_file, 'fastq', fasta_file, 'fasta')


def add_id(fasta_file):
    """
    add ID to each oligo in a fasta file.
    """
    newidsequence = []
    for counter, record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        newid = record.id + '_' + str(counter)
        newidsequence.append(SeqRecord(record.seq, id=newid, description=record.description))
    SeqIO.write(newidsequence, fasta_file, 'fasta')


def dataframe_from_oligominer(oligominer_fastq):
    """
    convert fastq of oligominer to pandas dataframe.
    """
    gene_id_list = []
    oligo_id_list = []
    seq_list = []
    hyb_start_list = []
    hyb_end_list = []
    with open(oligominer_fastq) as handle:
        for count, record in enumerate(SeqIO.parse(handle, 'fastq')):
            gene_id_list.append(record.id)
            oligo_id_list.append(record.id + '_' + str(count))
            seq_list.append(record.seq)
            hyb_start_end = record.description.split(':')[-1].split('-')
            hyb_start_list.append(int(hyb_start_end[0]))
            hyb_end_list.append(int(hyb_start_end[1]))
    fastq_df = pd.DataFrame(
        {'gene': gene_id_list, 'oligo_id': oligo_id_list, 'seq': seq_list, 'start': hyb_start_list,
         'end': hyb_end_list},
        columns=['gene', 'oligo_id', 'seq', 'start', 'end']
    )
    return fastq_df
