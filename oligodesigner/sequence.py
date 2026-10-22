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


def get_off_targeting_oligo_detail(minimum_offtarget_gap, df_blast):
    """
    Identify oligos that potentially cause an off-target signal and record why.

    An oligo is flagged when it and another candidate hit the same (non-self) transcript
    within minimum_offtarget_gap of each other, because such a pair of split probes can
    initiate HCR on the off-target transcript.

    Returns:
        oligo_to_be_removed (list): qseqid of the flagged oligos, in the order they were flagged.
        detail_df (DataFrame): one row per flagged oligo with columns
            'oligo_id', 'exclusion_sseqid' (transcript on which the pair was found),
            'exclusion_sstart', 'exclusion_send' (hit coordinates on that transcript),
            'exclusion_partner' (neighbouring oligo(s) within the gap, ';'-joined),
            'exclusion_gap' (distance to the nearest partner, measured from the start of one binding
            site to the end of the other; this is the value compared against minimum_offtarget_gap).
    """
    oligo_to_be_removed = []
    records = []
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
                records.append(_exclusion_record(nonspecific_oligo, nonspecific_gene, sstart_send,
                                                 minimum_offtarget_gap))
    detail_df = pd.DataFrame(records, columns=['oligo_id', 'exclusion_sseqid', 'exclusion_sstart',
                                               'exclusion_send', 'exclusion_partner', 'exclusion_gap'])
    return oligo_to_be_removed, detail_df


def _exclusion_record(oligo_id, sseqid, sstart_send, minimum_offtarget_gap):
    """
    Build the reason record for one flagged oligo from the sorted hit table of one transcript.
    sstart_send is indexed by qseqid and has columns sstart, send, dist_before, dist_after.
    """
    ids = list(sstart_send.index)
    # the same oligo may hit the same transcript more than once; take the first flagged occurrence.
    positions = [i for i, q in enumerate(ids) if q == oligo_id]
    partners = []
    gaps = []
    for pos in positions:
        row = sstart_send.iloc[pos]
        if pos > 0 and abs(row['dist_before']) < minimum_offtarget_gap:
            partners.append(ids[pos - 1])
            gaps.append(abs(row['dist_before']))
        if pos < len(ids) - 1 and abs(row['dist_after']) < minimum_offtarget_gap:
            partners.append(ids[pos + 1])
            gaps.append(abs(row['dist_after']))
        if partners:
            break
    row = sstart_send.iloc[positions[0]]
    return {'oligo_id': oligo_id,
            'exclusion_sseqid': sseqid,
            'exclusion_sstart': int(row['sstart']),
            'exclusion_send': int(row['send']),
            'exclusion_partner': ';'.join(dict.fromkeys(partners)),
            'exclusion_gap': int(min(gaps)) if gaps else pd.NA}


def get_off_targeting_oligo(minimum_offtarget_gap, df_blast):
    """
    Return the list of oligos that potentially cause an off-target signal.
    See get_off_targeting_oligo_detail for the criterion and for the per-oligo reason.
    """
    oligo_to_be_removed, _ = get_off_targeting_oligo_detail(minimum_offtarget_gap, df_blast)
    return oligo_to_be_removed


def candidate_report(oligominer_df, df_blast, self_hit, detail_df):
    """
    Build a debugging table with one row per OligoMiner candidate, in the same orientation and
    with the same leading columns as the oligo set csv, followed by the BLAST verdict.

    Arguments:
        oligominer_df: dataframe from io.dataframe_from_oligominer (target orientation).
        df_blast: the full BLAST result (before self hits are removed).
        self_hit: boolean Series aligned to df_blast, True for hits on the gene of interest.
        detail_df: second return value of get_off_targeting_oligo_detail.

    Columns:
        gene, oligo_id, seq (probe orientation), start, end,
        excluded, exclusion_sseqid, exclusion_sstart, exclusion_send, exclusion_partner, exclusion_gap,
        n_offtarget_hits, n_offtarget_transcripts, offtarget_sseqids (';'-joined), n_self_hits
    """
    report = oligominer_df.copy()
    report['seq'] = report['seq'].apply(lambda x: str(x.reverse_complement()))

    report['excluded'] = report['oligo_id'].isin(detail_df['oligo_id'])
    report = report.merge(detail_df, on='oligo_id', how='left')

    offtarget = df_blast[~self_hit]
    grouped = offtarget.groupby('qseqid')['sseqid']
    n_hits = grouped.size().rename('n_offtarget_hits')
    n_transcripts = grouped.nunique().rename('n_offtarget_transcripts')
    sseqids = grouped.apply(lambda s: ';'.join(s.drop_duplicates())).rename('offtarget_sseqids')
    n_self = df_blast[self_hit].groupby('qseqid').size().rename('n_self_hits')

    summary = pd.concat([n_hits, n_transcripts, sseqids, n_self], axis=1)
    report = report.merge(summary, left_on='oligo_id', right_index=True, how='left')

    for col in ['n_offtarget_hits', 'n_offtarget_transcripts', 'n_self_hits']:
        report[col] = report[col].fillna(0).astype(int)
    report['offtarget_sseqids'] = report['offtarget_sseqids'].fillna('')
    report['exclusion_partner'] = report['exclusion_partner'].fillna('')
    for col in ['exclusion_sstart', 'exclusion_send', 'exclusion_gap']:
        report[col] = report[col].astype('Int64')

    return report


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