import mygene
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import pandas as pd


def refseqid_from_fasta(fasta):
    record = list(SeqIO.parse(fasta, 'fasta'))

    return record[0].id


def fetch_mygene(genename, spiecies):
    genename = genename.upper()
    mg = mygene.MyGeneInfo()
    geneinfo = mg.querymany(genename, scopes='symbol', fields='ensembl.gene', species=spiecies,
                            verbose=False)
    return geneinfo


def fetch_mygene_from_refseqid(refseqid, species=None):
    mg = mygene.MyGeneInfo()
    geneinfo = mg.querymany(refseqid, scopes='refseq', fields='ensembl.gene', species=species,
                            verbose=False)
    return geneinfo


def check_refseqid_exist(fasta, return_entrezid=True):
    refseqid = refseqid_from_fasta(fasta)
    gene_info = fetch_mygene_from_refseqid(refseqid)

    if '_id' in gene_info[0].keys():
        print(f'refseq id {refseqid} was found')
        entrezid = gene_info[0]['_id']
    else:
        refseqid = None
        entrezid = None
        print(f'no refseq id was found')

    if return_entrezid:
        return entrezid
    else:
        return refseqid


def run_blast(input_file, output_file, database, task='blastn-short', strand='plus', num_threads=1):
    """
    run blast against transcriptome database. negative_seqid option does not work for an unknown reason.
    """
    blastn_cline = NcbiblastnCommandline(
        task=task,
        query=input_file,
        db=database,
        num_threads=num_threads,
        strand=strand,
        outfmt='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        out=output_file
    )
    stdout, stderr = blastn_cline()


def run_blast_df(input_file, database, task='blastn-short', strand='plus', num_threads=1):
    blastn_cline = NcbiblastnCommandline(
        task=task,
        query=input_file,
        db=database,
        num_threads=num_threads,
        strand=strand,
        outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        out='-'
    )
    output = blastn_cline()[0].strip()
    rows = [line.split() for line in output.splitlines()]
    cols = ['qseqid', 'sseqid', 'pident', 'length',
            'mismatch', 'gapopen', 'qstart', 'qend',
            'sstart', 'send', 'evalue', 'bitscore']

    data_types = {'pident': float, 'length': int, 'mismatch': int,
                  'gapopen': int, 'qstart': int, 'qend': int,
                  'sstart': int, 'send': int, 'evalue': float,
                  'bitscore': float}

    df = pd.DataFrame(rows, columns=cols).astype(data_types)

    return df


def get_homology_in_database(fasta, database, scopes='ensembl.transcript', num_threads=1, similarity_threshold=0.15):
    """

    """
    df = run_blast_df(fasta, database, task='blastn', strand='plus', num_threads=num_threads)
    transcript_series = df['sseqid'].str.split('.').str[0]

    mg = mygene.MyGeneInfo()
    df_transcript_id = pd.DataFrame(
        mg.querymany(transcript_series, scopes=scopes, fields='id', verbose=False))['_id']


    gp = pd.Series((df['pident'] * df['length']).values, index=df_transcript_id.values)
    gp = gp.groupby(level=0).sum()
    gp = (gp / gp.sum()).sort_values(ascending=False)

    top_hit = gp.index[0]
    seqid = top_hit


    res = mg.querymany(top_hit, scopes='entrezgene', fields='symbol', verbose=False)[0]['symbol']
    print(f'The sequence is the most similar to {res} with the score of {gp[0]}')

    if gp[gp > similarity_threshold].size > 1:
        print(f'Homologous transcript(s) was identified. This may impact on the selection of the oligos')
        for hit in gp[gp > similarity_threshold].index:
            if hit != top_hit:
                res = mg.querymany(hit, scopes='entrezgene', fields='symbol', verbose=False)[0]['symbol']
                print(f'Homologous transcript {res} was found with the score of {gp.loc[hit]}')

    return seqid


def self_hit_mask(df_blast, gene_id):
    """
    return a boolean Series aligned to df_blast. True marks a hit on the gene of interest,
    judged by referring the gene information from mygene.
    """
    mg = mygene.MyGeneInfo()
    transcript_series = df_blast['sseqid'].str.split('.').str[0]
    df_transcript_id = pd.DataFrame(mg.querymany(transcript_series.drop_duplicates(keep='first'),
                                                 scopes='ensembl.transcript', fields='id', verbose=False))

    if not '_id' in df_transcript_id.columns:
        return pd.Series(False, index=df_blast.index, name='self_hit')
    else:
        self_transcripts = df_transcript_id[df_transcript_id['_id'] == gene_id]['query']
        return transcript_series.isin(self_transcripts).rename('self_hit')


def exclude_self(df_blast, gene_id):
    """
    omit the gene of interest from the blast result
    by referring the gene information from mygene.
    """
    return df_blast[~self_hit_mask(df_blast, gene_id)].reset_index(drop=True)


def get_off_targeting_oligo(minimum_offtarget_gap, df_blast):
    """

    """
    oligo_to_be_removed = []
    df_nonspecific_oligo = df_blast[df_blast['sseqid'].duplicated(keep=False)]
    # check each oligo whether the oligo potentially causes a off-target signal
    for nonspecific_oligo in list(df_nonspecific_oligo['qseqid'].unique()):
        nonspecific_gene_set = (df_nonspecific_oligo[df_nonspecific_oligo['qseqid'] == nonspecific_oligo])['sseqid']
        # check each nonspecific binding site.
        for nonspecific_gene in list(nonspecific_gene_set.unique()):
            # calculate distance to the nearest off-target binding site.
            sstart_send = df_nonspecific_oligo[df_nonspecific_oligo['sseqid'] == nonspecific_gene]
            sstart_send = (sstart_send.set_index('qseqid', drop=True))[['sstart', 'send']]
            sstart_send = sstart_send.sort_values(by=['sstart'], ascending=True)
            dist_before = (sstart_send['send'] - sstart_send['sstart'].shift(periods=1)).rename('dist_before')
            sstart_send = pd.concat([sstart_send, dist_before, (dist_before.shift(periods=-1)).rename('dist_after')],
                                    axis=1)
            # see whether the binding site is within the intolerable range.
            leaky_set = sstart_send[(sstart_send['dist_before'].abs() < minimum_offtarget_gap) | (
                        sstart_send['dist_after'].abs() < minimum_offtarget_gap)]

            # if oligo is causing off-target signal, remove.
            if nonspecific_oligo in list(leaky_set.index):
                df_nonspecific_oligo = df_nonspecific_oligo[~(df_nonspecific_oligo['qseqid'] == nonspecific_oligo)]
                # update the oligo list
                oligo_to_be_removed.append(nonspecific_oligo)
    return oligo_to_be_removed


def remove_oligo(oligominer_df, oligo_to_be_removed):
    """
    """
    selected_oligo_df = oligominer_df[~(oligominer_df['oligo_id'].isin(oligo_to_be_removed))]
    selected_oligo_df = selected_oligo_df.reset_index(drop=True)
    return selected_oligo_df


def reverse_complement_dataframe(df):
    """
    """
    df['seq'] = df['seq'].apply(lambda x: x.reverse_complement())

    return df

def calc_oligo_interval(oligo_df, mean_scan_range=48):
    """
    """
    dist_before = (oligo_df['start'] - oligo_df['end'].shift(periods=1)).rename('interval_before')
    dist_after = (oligo_df['start'].shift(periods=-1) - oligo_df['end']).rename('interval_after')

    mean_interval = dist_before.rolling(mean_scan_range, min_periods=mean_scan_range).mean()
    mean_interval = mean_interval.rename('mean_interval')
    interval_df = pd.concat([dist_after, mean_interval], axis=1)

    return interval_df


def add_hcr_seq_v3(oligo_df, hcr_initiator_set):
    """
    add HCR sequences to both sides of split probes.
    """
    seq_even_l = hcr_initiator_set['seq_even_l']
    seq_odd_l = hcr_initiator_set['seq_odd_l']
    seq_even_r = hcr_initiator_set['seq_even_r']
    seq_odd_r = hcr_initiator_set['seq_odd_r']

    # anchorseq_ll = hcr_initiator_set['anchorseq_ll']
    # anchorseq_rr = hcr_initiator_set['anchorseq_rr']
    # anchorseq_lr = hcr_initiator_set['anchorseq_lr']
    # anchorseq_rl = hcr_initiator_set['anchorseq_rl']

    oligo_df = oligo_df.sort_values(by='start')

    even_oligo = oligo_df[oligo_df.index % 2 == 0]
    even_hcr_oligo = (seq_even_l + even_oligo['seq'] + seq_even_r).rename('hcr_seq')
    even_oligo = pd.concat([even_oligo, even_hcr_oligo], axis=1)

    odd_oligo = oligo_df[oligo_df.index % 2 == 1]
    odd_hcr_oligo = (seq_odd_l + odd_oligo['seq'] + seq_odd_r).rename('hcr_seq')
    odd_oligo = pd.concat([odd_oligo, odd_hcr_oligo], axis=1)

    hcr_oligo_set = (pd.concat([even_oligo, odd_oligo], axis=0)).sort_values(by='start')

    return hcr_oligo_set