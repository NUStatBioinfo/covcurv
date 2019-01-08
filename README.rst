=================================================================
covcurv: per-gene coverage arrays pipeline and visualization
=================================================================

.. image:: img/app_live.png
   :height: 150px
   :width: 500px
   :align: center


While there already exist many useful tools for computing genome-wide reads coverage from individual next-generation sequencing experiments,
`covcurv` is the first **per-gene**, **multi-sample** reads coverage parsing and visualization tool. `covcurv` is a command line interface (CLI) tool
that takes aligned read (.bam) files, along with a gene annotation (.gtf) file, and gives you per-gene reads coverage arrays and visualization capabilities.

+++++++++++++++++++
Why covcurv?
+++++++++++++++++++

It's difficult to access gene transcript coverage data, let alone from multiple sequencing experiments. For example, the `genomecov` function (as part of the `bedtools suite <https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html>`_.) will
report coverage for all genome positions in convenient bedgraph format, but it cannot process multiple .bam files at once. Its counting logic for paired reads is undesirable:
genomecov either double counts paired read overlap regions (using `-split`), or it considers all positions between two paired reads, even when there's a gap between them (with `-pc`). Further,
there is no option to break up coverage into the coding regions that comprise a gene.
`htseq-count`'s accounting principles are more preferable, but it still does not provide coverage for multiple samples at once, and it only counts reads per feature, not per base position.

If you want gene transcript coverage based on desirable CIGAR string-based read accounting logic along with an easy way to view and export this coverage data, `covcurv` is for you.

+++++++++++++++++++
How does it work?
+++++++++++++++++++


1. **Parse .gtf file** into (gene, exon)-tuple metadata. Exons at the intersection of more than a single gene are removed.
    - Coverage arrays are a concatenation of a gene's exons.

2. **Aligned reads files** are loaded (optionally in parallel, over chromosomes) and transformed into chromosome-wide coverage arrays.
    - Coverage is computed according to CIGAR score. Only "M" (match) segments contribute to coverage.
    - Paired read overlap only counts for "+1" coverage (it is not double counted).

3. **Chromosome-wide coverage arrays are diced into gene arrays** in a memory efficient manner.

.. image:: img/chrom_coverage_dicing.png
   :height: 200px
   :width: 50px
   :align: center

+++++++++++++++++++
Getting started
+++++++++++++++++++

`covcurv` requires aligned reads (.bam) files. Those files must be sorted in genome order with the `samtools index` function. If you haven't created them, required
bam index (.bai) files will also be created for you if you have `samtools` in your `$PATH`. Otherwise, you will need to supply bam index files. Finally, you'll need
a genome annotation (.gtf) file. `covcurv` works with paired and single reads.

1. **Sort and index your alignments** with `samtools`
    # sort and index alignment files in data directory
    for FILE in ezh2_data/*.bam
    do
        samtools sort $FILE -o ${FILE/.bam/}'_sorted.bam'
        samtools index ${FILE/.bam/}'_sorted.bam' ${FILE/.bam/.bai}
    done

2. **Run the `covcurv` pipeline**. Use `-p` to run with the pipeline in parallel with hyperthreading.

    # run pipeline on all alignments in ezh2_data directory
    covcurv --bam-dir ezh2_data \
     -o output_dir/covcurv_ezh2 \
     -g ezh2_data/genes.gtf \
     -p 4

3. **Start the web app visualization tool with the `covcurv_app` command`**. Just point `covcurv_app` to a `covcurv` output directory using the `-d` flag (for "data directory") and navigate to the web app's URL in your browser.

    covcurv_app -d output_dir/covcurv_ezh2

.. image:: img/app_server.png
   :height: 150px
   :width: 500px
   :align: center

--------------
Installation
--------------

