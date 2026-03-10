rule deeptools_bamcoverage:
    input:
        # Required input.
        bam = 'YOUR_BAM_FILE',
    output:
        # Required output.
        # Output file format should be one of ['bw', 'bigwig', 'bigWig', 'bedgraph', 'bedGraph'].
        output = 'YOUR_OUTPUT_FILE'
    params:
        # Optional parameters.
        extra = '',
        # The computed scaling factor (or 1, if not applicable) will be multiplied by this.
        scale_factor = 1,
        # Bin size (default: 50).
        # The smaller the bin size, the bigger the output file will be.
        bin_size = 1,
        # Region of the genome to limit the operation to.
        # e.g. region = 'chr2:10000000-10100000'
        region = '',
        # The effective genome size is the portion of the genome that is mappable.
        # A table of values is available:
        # GRCh37 	2864785220
        # GRCh38 	2913022398
        # GRCm37 	2620345972
        # GRCm38 	2652783500
        # dm3 	162367812
        # dm6 	142573017
        # GRCz10 	1369631918
        # WBcel235 	100286401
        effective_genome_size = 2913022398,  # For GRCh38.
        # Use one of the entered methods to normalize the number
        # of reads per bin. By default, no normalization is perfomred.
        # Available options are:
        # RPKM: Reads Per Kilobase per Million mapped reads.
        # CPM: Counts Per Million mapped reads.
        # BPM: Bins Per Million mapped reads.
        # RPGC: Reads Per Genomic Context (1x normalization).
        normalize_using = 'BPM',
        # This parameter determines if non-covered regions
        # (regions without overlapping reads) in a BAM file
        # should be skipped.
        skip_non_covered_regions = False,
        # The smooth length defines a window, larger than the binSize,
        # to average the number of reads. For example, if the --binSize is set to 20
        # and the --smoothLength is set to 60, then for each bin, the average of the bin
        # and its left and right neighbors is considered.
        # Any value smaller than --binSize will be ignored and no smoothing will be applied.
        smooth_length = 5,
        # This parameter allows the extension of reads to fragment size. If set, each read is
        # extended, without exception.
        # Single end: Requires a user specified value for the final fragment length.
        # Reads that already exceed this fragment length will not be extended.
        # Paired end: Reads with mates are always extended to match the fragment
        # size defined by the two read mates.
        extend_reads = 200,
    threads: 1
    log: 'logs/deeptools_bamcoverage/{sample}.log'
    benchmark: 'benchmark/deeptools_bamcoverage/{sample}.txt'
    wrapper: 'http://dohlee-bio.info:9193/deeptools/bamcoverage'