from flask_wtf import FlaskForm as Form
from wtforms import validators, BooleanField, StringField, DecimalField, FileField, MultipleFileField, FloatField, IntegerField, SelectField, SelectMultipleField, SubmitField, StringField


class testForm(Form):
    outFile=StringField(label='Out  File', validators=[], description='The file to write results to. By default, results are printed to the console')
    labels=StringField(label='Labels', validators=[], description='Labels for the samples. The default is to use the file name of the sample. The sample labels should be separated by spaces and quoted if a label itselfcontains a space E.g. -labels label-1 "label 2"  ')
    smartLabels=StringField(label='Smart  Labels', validators=[], description='Instead of manually specifying labels for the input BAM files, this causes sincei to use the file name after removing the path and extension.')
    binSize=IntegerField(label='Bin  Size', validators=[], description='Length in bases of the window used to sample the genome. ')
    distanceBetweenBins=IntegerField(label='Distance  Between  Bins', validators=[], description='To reduce the computation time, not every possible genomic bin is sampled. This option allows you to set the distance between bins actually sampled from. Larger numbers are sufficient for high coverage samples, while smaller values are useful for lower coverage samples. Note that if you specify a value that results in too few (<1000) reads sampled, the value will be decreased. ')
    numberOfProcessors=IntegerField(label='Number  Of  Processors', validators=[], description='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. ')
    filterRNAstrand=SelectField(label='Filter Rn Astrand', validators=[], description='Selects RNA-seq reads (single-end or paired-end) in the given strand. ')
    minMappingQuality=IntegerField(label='Min  Mapping  Quality', validators=[], description='If set, only reads that have a mapping quality score of at least this are considered.')
    samFlagInclude=IntegerField(label='Sam  Flag  Include', validators=[], description='Include reads based on the SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage. ')
    samFlagExclude=IntegerField(label='Sam  Flag  Exclude', validators=[], description='Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use -samFlagExclude 16, where 16 is the SAM flag for reads that map to the reverse strand. ')
    blackListFileName=MultipleFileField(label='Black  List  File  Name', validators=[], description='A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.')
    submit = SubmitField('Submit')
