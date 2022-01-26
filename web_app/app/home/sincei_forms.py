from flask_wtf import FlaskForm as Form
from wtforms import validators, BooleanField, StringField, DecimalField, FileField, MultipleFileField, FloatField, IntegerField, SelectField, SelectMultipleField, SubmitField, StringField


class form_scFilterStats(Form):
    tagName=StringField(label='Tag  Name', validators=[], default='BC', description='Name of the BAM tag from which to extract barcodes.')
    duplicateFilter=SelectField(label='Duplicate  Filter', validators=[],
                     choices=['start_bc', 'start_bc_umi', 'start_end_bc', 'start_end_bc_umi'], default=None, description='How to filter for duplicates? Different combinations (using start/end/umi) are possible. Read start position and read barcode are always considered. Default (None) would consider all reads. Note that in case of paired end data, both reads in the fragment are considered (and kept). So if you wish to keep only read1, combine this option with samFlagInclude ')
    motifFilter=StringField(label='Motif  Filter', validators=[], default=None, description='Check whether a given motif is present in the read and the corresponding reference genome. This option checks for the motif at the 5-end of the read and at the 5-overhang in the genome, which is useful in identifying reads properly cut by a restriction-enzyme or MNAse. For example, if you want to search for an "A" at the 5\'-end of the read and "TA" at 5\'-overhang, use "-m \'A,TA\'". Reads not containing the given motif are discarded. ')
#    genome2bit=StringField(label='Genome2Bit', validators=[], default=None, description='If --motifFilter is provided, please also provide the genome sequence (in 2bit format). ')
    GCcontentFilter=StringField(label=' G Ccontent  Filter', validators=[], default=None, description='Check whether the GC content of the read falls within the provided range. If the GC content of the reads fall outside the range, they are discarded. ')
    minAlignedFraction=FloatField(label='Min  Aligned  Fraction', validators=[], default=None, description='Minimum fraction of the reads which should be aligned to be counted. This includes mismatches tolerated by the aligners, but excludes InDels/Clippings')
#    outFile=FileField(label='Out  File', validators=[], default=None, description='The file to write results to. By default, results are printed to the console')
    binSize=IntegerField(label='Bin  Size', validators=[], default='1000000', description='Length in bases of the window used to sample the genome.')
    distanceBetweenBins=IntegerField(label='Distance  Between  Bins', validators=[], default='1000000', description='To reduce the computation time, not every possible genomic bin is sampled. This option allows you to set the distance between bins actually sampled from. Larger numbers are sufficient for high coverage samples, while smaller values are useful for lower coverage samples. Note that if you specify a value that results in too few (<1000) reads sampled, the value will be decreased.')
    numberOfProcessors=IntegerField(label='Number  Of  Processors', validators=[], default=1, description='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors.')
    filterRNAstrand=SelectField(label='Filter RNA strand', validators=[],
                     choices=['forward', 'reverse'], default=None, description='Selects RNA-seq reads (single-end or paired-end) in the given strand.')
    minMappingQuality=IntegerField(label='Min  Mapping  Quality', validators=[], default=None, description='If set, only reads that have a mapping quality score of at least this are considered.')
    samFlagInclude=IntegerField(label='Sam  Flag  Include', validators=[], default=None, description='Include reads based on the SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage.')
    samFlagExclude=IntegerField(label='Sam  Flag  Exclude', validators=[], default=None, description='Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use --samFlagExclude 16, where 16 is the SAM flag for reads that map to the reverse strand.')
#    blackListFileName=StringField(label='Black  List  File  Name', validators=[], default=None, description='A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.')
    submit = SubmitField('Submit')


class form_scPlotUMAP(Form):
    geneName=StringField(label='Gene', validators=[], default=None, description='Provide a gene or region to explore')
    submit=SubmitField('Submit')
