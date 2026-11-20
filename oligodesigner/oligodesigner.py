from oligodesigner import blockParse_py3
from oligodesigner import io
from oligodesigner import parse
from oligodesigner import sequence
import pandas as pd


def generator(mFISH3D_param, oligominer_param, num_threads=1):
    """
    Arguments:
        mFISH3D_param (dict)
            'fasta': fasta file of the template sequence
            'database': database for blast
            'minimum_offtarget_gap': If the gap between two non-specific binding is more than minimum_offtarget_gap,
                the pair is not regarded to cause a off-target signal. Recommended value: 100.
            'hcr_seqs': sequences of HCR fragments. The example design can be found in oligodesigner/parameters.py
            'self_remove': Recommended: True. set self_remove True if you want to remove the template sequence from
                off-target analysis. Otherwise, the gene of your interest could be regarded as an off-target product.
                Turn this to False when your template is not found in database (e.g. GFP).
            'export_blast': Optional, default False. Set True to write two debugging files next to the
                input fasta. <fasta basename>_blast.csv is the raw BLAST result; its extra column 'self_hit'
                is True for hits on the template gene, which are excluded from the off-target analysis
                (all False when self_remove is False). <fasta basename>_offtarget_screening.csv lists every
                OligoMiner candidate in the oligo set format, with 'excluded' and the reason: the
                transcript ('exclusion_sseqid') and partner oligo ('exclusion_partner') that form a
                potential off-target pair, plus all non-self transcripts the oligo hits ('offtarget_sseqids').
            'minimum_offtarget_homology': Optional, default 15. BLAST hits whose alignment length is
                shorter than this many bases are ignored in the off-target analysis. Set 0 to disable;
                then the shortest hit is decided by the e-value cutoff in 'blast_options' and by the
                database size (about 14 nt for a 20-mer against a transcriptome). The alignment length
                counts mismatches and gaps, so a 16 nt alignment with 2 mismatches passes 15.
            'blast_options': Optional. Dict of blastn options for the off-target screening, merged over
                sequence.BLASTN_SHORT_OPTIONS (word_size 7, evalue 10, reward 1, penalty -3, gapopen 5,
                gapextend 2; the implicit defaults of '-task blastn-short'). Keys are blastn option names
                without the leading dash, e.g. {'evalue': 1000}. Lowering evalue or raising word_size
                reports fewer and longer hits; raising evalue reports shorter hits (down to ~13 nt at 1000).
            The values actually used, including the defaults, are written to <fasta basename>_parameters.csv.
    """
    fasta = mFISH3D_param['fasta']
    database = mFISH3D_param['database']
    minimum_offtarget_gap = mFISH3D_param['minimum_offtarget_gap']
    hcr_initiator_set = mFISH3D_param['hcr_seqs']
    self_remove = mFISH3D_param['self_remove']
    export_blast = mFISH3D_param.get('export_blast', False)
    minimum_offtarget_homology = mFISH3D_param.get('minimum_offtarget_homology', 15)
    blast_options = {**sequence.BLASTN_SHORT_OPTIONS, **mFISH3D_param.get('blast_options', {})}

    # Make path to temp files and result file
    (oligominer_fastq, oligominer_fasta, oligo_sets_path, param_path, blast_path, screening_path) = io.gen_file_path(fasta)

    # Save parameters to csv file, including the resolved values of the optional parameters.
    io.record_param({**mFISH3D_param,
                     'minimum_offtarget_homology': minimum_offtarget_homology,
                     'blast_options': blast_options},
                    oligominer_param, param_path)

    # Run OligoMiner
    blockParse_py3.runSequenceCrawler(fasta, *parse.oligominer_param_parser(oligominer_param))

    # read oligominer fastq as dataframe
    oligominer_df = io.dataframe_from_oligominer(oligominer_fastq)

    # converting fastq to fasta. writing fasta is a must for blast.
    io.convert_fastq2fasta(oligominer_fastq, oligominer_fasta)
    io.add_id(oligominer_fasta)

    # run blast
    df_blast = sequence.run_blast_df(oligominer_fasta, database, task='blastn-short', strand='plus',
                                     num_threads=num_threads, **blast_options)

    # Flag hits that are too short to count as off-target binding.
    short_hit = (df_blast['length'] < minimum_offtarget_homology).rename('short_hit')

    # Flag hits on the gene of interest. Nothing is flagged when self_remove is False.
    if self_remove:
        seqid = sequence.check_refseqid_exist(fasta, return_entrezid=True)
        if seqid is None:
            seqid = sequence.get_homology_in_database(fasta, database, num_threads=num_threads)
        else:
            pass
        self_hit = sequence.self_hit_mask(df_blast, seqid)
    else:
        self_hit = pd.Series(False, index=df_blast.index, name='self_hit')

    # Exclude hits on the gene of interest and hits shorter than minimum_offtarget_homology from the result
    df_blast_all = df_blast
    df_blast = df_blast[~self_hit & ~short_hit].reset_index(drop=True)

    # Remove combination if sstart and send are close enough and if they are rective combination.
    oligo_to_be_removed, exclusion_detail = sequence.get_off_targeting_oligo_detail(minimum_offtarget_gap, df_blast)
    print("These sequences were removed because BLAST detected potential off-target binding", oligo_to_be_removed)

    if export_blast:
        df_blast_all.assign(self_hit=self_hit, short_hit=short_hit).to_csv(blast_path, index=False)
        sequence.candidate_report(oligominer_df, df_blast_all, self_hit, exclusion_detail, short_hit).to_csv(
            screening_path, index=False)

    oligo_df = sequence.remove_oligo(oligominer_df, oligo_to_be_removed)

    oligo_df = sequence.reverse_complement_dataframe(oligo_df)

    interval = sequence.calc_oligo_interval(oligo_df)
    oligo_df = pd.concat([oligo_df, interval], axis=1)

    # Make oligos with HCR sequences.
    res = sequence.add_hcr_seq_v3(oligo_df, hcr_initiator_set)
    res.to_csv(oligo_sets_path, sep=',', index=False, header=True)

    return res

