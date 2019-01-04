#!/usr/bin/env python3

# ---------------------------------------------------------------------------- #
# covcurv CLI entrypoint for use on a single node with hyperthreading.
# ---------------------------------------------------------------------------- #

import sys
from covcurv.reads import *
from covcurv.coverage import *
from covcurv.gene_processing import *
from collections import OrderedDict


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments, display welcome message, create output directory
    # ---------------------------------------------------------------------------- #
    args = parse_args()
    n_jobs = args.proc_per_node

    output_dir = create_output_dir(args.output_dir)

    configure_logger(output_dir)
    welcome()
    logging.info('covcurv output directory -- {0}'.format(output_dir))

    # ---------------------------------------------------------------------------- #
    # If any Bam index (.bai) files need to be created, do that now.
    # ---------------------------------------------------------------------------- #
    if args.create_bai_files:
        for file_idx in range(len(args.create_bai_files)):
            bam_file = args.create_bai_files[file_idx]
            logging.info('creating index file for {0} -- {1} / {2}'
                         .format(bam_file, file_idx + 1, len(args.create_bai_files)))
            out = create_index_file(bam_file)

    # ---------------------------------------------------------------------------- #
    # .bam file preprocessing:
    # Determine intersection of chromosomes across samples from .bam files
    # ---------------------------------------------------------------------------- #
    sample_ids = list()
    chroms = list()
    cov_files = OrderedDict()
    read_count_dict = dict()
    n_samples = len(args.bam_files)

    # load each .bam file's header, find joint intersection of read chromosomes.
    for idx in range(n_samples):
        header_dat = BamReadsProcessor(args.bam_files[idx]
                                        , index_file=args.bai_files[idx]).header
        new_chroms = header_dat.chr.values.tolist()

        if not chroms:
            chroms = new_chroms

        else:
            chroms = np.intersect1d(chroms, new_chroms).tolist()

    # ---------------------------------------------------------------------------- #
    # Load .gtf or .gff files and run processing pipeline.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #
    logging.info('Begin genome annotation file processing...')
    gap = GeneAnnotationProcessor(args.genome_annotation
                                  , verbose=True
                                  , chroms=chroms)
    exon_df = gap.run()
    genes_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

    # take intersection of chromosomes available in genome annotation file and those in the reads data,
    # if for some reason annotation file only contains subset.
    chroms = np.intersect1d(chroms, genes_df.chr.unique()).tolist()
    logging.info('Found {0} chromosomes in intersection of all experiments and gene annotation data:\n'
                 '\t{1}'.format(len(chroms), ', '.join(chroms)))

    # ---------------------------------------------------------------------------- #
    # Load .bam files while simultaneously parsing into coverage arrays. Store chromosome
    # coverage arrays and and compute gene read counts.
    # ---------------------------------------------------------------------------- #

    # iterate over .bam files; compute each sample's chromosomes' coverage arrays
    # and save them to .npz files.
    for idx in range(n_samples):
        logging.info('Loading RNA-seq data file {0} / {1}'.format(idx + 1, n_samples))

        reader = BamReadsProcessor(bam_file=args.bam_files[idx]
                                   , index_file=args.bai_files[idx]
                                   , chroms=chroms
                                   , n_jobs=n_jobs
                                   , output_dir=output_dir
                                   , verbose=True)

        sample_id = reader.sample_id
        sample_ids.append(sample_id)
        cov_files[sample_id], read_count_dict[sample_id] = reader.coverage_read_counts(genes_df)

    logging.info('Successfully processed chromosome read coverage and gene read counts for all {0} experiments'
                 .format(len(sample_ids)))

    del reader
    gc.collect()

    # ---------------------------------------------------------------------------- #
    # Merge per-sample gene read count matrices:
    # obtain read count DataFrame containing X, an n (genes) x p (samples) matrix.
    # ---------------------------------------------------------------------------- #
    read_count_df = read_count_dict[sample_ids[0]]
    read_count_df.rename(columns={'read_count': sample_ids[0]}, inplace=True)

    for sample_id in sample_ids[1:]:
        read_count_df = pd.merge(read_count_df
                                 , read_count_dict[sample_id].rename(columns={'read_count': sample_id})
                                 , how='inner'
                                 , on=['chr', 'gene'])

    logging.info('Successfully merged sample read counts -- shape: {0}'.format(read_count_df.shape))

    del read_count_dict
    gc.collect()

    # ---------------------------------------------------------------------------- #
    # Slice up genome coverage matrix for each gene according to exon positioning.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #
    gene_cov_dicts = Parallel(n_jobs=min(n_jobs, len(chroms))
                              , verbose=0
                              , backend='threading')(delayed(gene_coverage)(
        exon_df=exon_df,
        chrom=chrom,
        coverage_files=cov_files,
        output_dir=output_dir,
        verbose=True) for chrom in chroms)

    # convert list of tuples into 2-d dictionary: {chrom: {gene: coverage matrix}}
    chrom_gene_cov_dict = {chroms[idx]: gene_cov_dicts[idx] for idx in range(len(chroms))}

    del gene_cov_dicts
    gc.collect()

    # ---------------------------------------------------------------------------- #
    # Remove traces of genes with null coverage.
    # ---------------------------------------------------------------------------- #

    # subset genes, exons to genes in intersection of experiments.
    genes_df = genes_df[genes_df.gene.isin(read_count_df.gene.unique())]
    exon_df = exon_df[exon_df.gene.isin(read_count_df.gene.unique())]

    # re-order genes, read counts so that they're parsimonious in chromosome + gene ordering.
    genes_df.sort_values(['chr', 'gene'], inplace=True)
    genes_df.reset_index(inplace=True, drop=True)
    read_count_df.sort_values(['chr', 'gene'], inplace=True)
    read_count_df.reset_index(inplace=True, drop=True)

    # storage for row indices of genes with null coverage.
    delete_idx = list()

    for i in range(genes_df.shape[0]):
        chrom = genes_df.chr.iloc[i]
        gene = genes_df.gene.iloc[i]

        # extract gene's p x Li coverage matrix.
        cov_mat = chrom_gene_cov_dict[chrom][gene]

        # do not add gene if gene has no reads covering it.
        if cov_mat.max() == 0:
            delete_idx.append(i)
            del chrom_gene_cov_dict[chrom][gene]

    # remove trace of genes with no reads.
    if delete_idx:
        read_count_df = read_count_df.drop(delete_idx
                                           , axis=0).reset_index(drop=True)
        genes_df = genes_df.drop(delete_idx
                                 , axis=0).reset_index(drop=True)

    # ---------------------------------------------------------------------------- #
    # Run quality control checks and
    # save gene annotation metadata and original read counts.
    # ---------------------------------------------------------------------------- #

    # check that read counts and coverage matrices contain data for same number of genes.
    n_cov_mats = sum([len(chrom_gene_cov_dict[chrom]) for chrom in chrom_gene_cov_dict])
    if n_cov_mats != read_count_df.shape[0]:
        raise ValueError('Number of coverage matrices not equal to number of genes in read count DataFrame!')

    # quality control.
    if genes_df.shape[0] != read_count_df.shape[0]:
        raise ValueError('Genes DataFrame and read counts DataFrame do not have same number of rows!')

    # save gene annotation metadata.
    exon_output_file = os.path.join(output_dir, 'gene_exon_metadata.csv')
    logging.info('Saving gene-exon metadata to {0}'.format(exon_output_file))
    exon_df.to_csv(exon_output_file
                   , index=False)

    # save read counts.
    read_count_file = os.path.join(output_dir, 'read_counts.csv')
    logging.info('Saving read counts to {0}'.format(read_count_file))
    read_count_df.to_csv(read_count_file
                         , index=False)

    del exon_df, read_count_df
    gc.collect()

    # ---------------------------------------------------------------------------- #
    # Save coverage matrices into per-chromosome directories.
    # ---------------------------------------------------------------------------- #
    for chrom in chrom_gene_cov_dict:
        n_genes = len(chrom_gene_cov_dict[chrom])

        if n_genes > 0:
            logging.info('Chromosome {0} -- saving {1} coverage matrices'.format(chrom, n_genes))
            chrom_dir = os.path.join(output_dir, chrom)

            if not os.path.exists(chrom_dir):
                os.makedirs(chrom_dir)

            with open(os.path.join(chrom_dir, 'coverage_matrices_{0}.pkl'.format(chrom)), 'wb') as f:
                pkl.dump(chrom_gene_cov_dict[chrom], f)

    del chrom_gene_cov_dict
    gc.collect()

    logging.info('covcurv pipeline complete! Exiting...')
    sys.exit(0)


if __name__ == "__main__":
    main()