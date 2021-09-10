import os
import sys
import shutil
import numpy as np
import pandas as pd
import pyBigWig
import math

# deeptools modules
from deeptools import mapReduce
from deeptools.utilities import getCommonChrNames
from deeptools import bamHandler
from deeptools import utilities
from deeptools.writeBedGraph import bedGraphToBigWig, getGenomeChunkLength

# own modules
import ReadCounter as cr

debug = 0


def scaleCoverage(tile_coverage, args):
    """
#    Return coverage per cluster as sum of cells.
#    tileCoverage should be an list with only one element
#    """
    return args['scaleFactor'] * tile_coverage

def writeBedGraph_wrapper(args):
    """
    Passes the arguments to writeBedGraph_worker.
    This is a step required given
    the constrains from the multiprocessing module.
    The args var, contains as first element the 'self' value
    from the WriteBedGraph object

    """
    return WriteBedGraph.writeBedGraph_worker(*args)



class WriteBedGraph(cr.CountReadsPerBin):

    r"""Reads bam files coverages and writes a bedgraph or bigwig file

    Extends the CountReadsPerBin object such that the coverage
    of bam files is writen to multiple bedgraph files at once.

    The bedgraph files are later merge into one and converted
    into a bigwig file if necessary.

    The constructor arguments are the same as for CountReadsPerBin. However,
    when calling the `run` method, the following parameters have
    to be passed

    Examples
    --------

    Given the following distribution of reads that cover 200 on
    a chromosome named '3R'::

          0                              100                           200
          |------------------------------------------------------------|
        A                                ===============
                                                        ===============


        B                 ===============               ===============
                                         ===============
                                                        ===============

    >>> import tempfile
    >>> test_path = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"

    >>> outFile = tempfile.NamedTemporaryFile()
    >>> bam_file = test_path +  "testA.bam"

    For the example a simple scaling function is going to be used. This function
    takes the coverage found at each region and multiplies it to the scaling factor.
    In this case the scaling factor is 1.5

    >>> function_to_call = scaleCoverage
    >>> funcArgs = {'scaleFactor': 1.5}

    Restrict process to a region between positions 0 and 200 of chromosome 3R

    >>> region = '3R:0:200'

    Set up such that coverage is computed for consecutive bins of length 25 bp
    >>> bin_length = 25
    >>> step_size = 25

    >>> num_sample_sites = 0 #overruled by step_size
    >>> c = WriteBedGraph([bam_file], binLength=bin_length, region=region, stepSize=step_size)
    >>> c.run(function_to_call, funcArgs, outFile.name)
    >>> f = open(outFile.name, 'r')
    >>> f.readlines()
    ['3R\t0\t100\t0\n', '3R\t100\t200\t1.5\n']
    >>> f.close()
    >>> outFile.close()


    """

    def run(self, func_to_call, func_args, out_file_prefix, blackListFileName=None, format="bedgraph", smoothLength=0, normUsing=None):
        r"""
        Given a list of bamfiles, a function and a function arguments,
        this method writes a bedgraph file (or bigwig) file
        for a partition of the genome into tiles of given size
        and a value for each tile that corresponds to the given function
        and that is related to the coverage underlying the tile.

        Parameters
        ----------
        func_to_call : str
            function name to be called to convert the list of coverages computed
            for each bam file at each position into a single value.
        func_args : dict
            dict of arguments to pass to `func`. E.g. {'scaleFactor':1.0}

        out_file_prefix : str
            name of the file to save the resulting data.


        """
        #self.__dict__["smoothLength"] = smoothLength
        getStats = len(self.mappedList) < len(self.bamFilesList)
        bam_handles = []
        for x in self.bamFilesList:
            if getStats:
                bam, mapped, unmapped, stats = bamHandler.openBam(x, returnStats=True, nThreads=self.numberOfProcessors)
                self.mappedList.append(mapped)
                self.statsList.append(stats)
            else:
                bam = bamHandler.openBam(x)
            bam_handles.append(bam)

        genome_chunk_length = getGenomeChunkLength(bam_handles, self.binLength, self.mappedList)
        # check if both bam files correspond to the same species
        # by comparing the chromosome names:
        chrom_names_and_size, non_common = getCommonChrNames(bam_handles, verbose=False)

        if self.region:
            # in case a region is used, append the tilesize
            self.region += ":{}".format(self.binLength)

        for x in list(self.__dict__.keys()):
            if x in ["mappedList", "statsList", "barcodes", "clusterInfo"]:
                continue
            sys.stderr.write("{}: {}\n".format(x, self.__getattribute__(x)))

        # below we get the same ouput as in deeptools, except that the 3rd list
        # element contains multiple tmp file names, one tmp file per cluster
        res = mapReduce.mapReduce([func_to_call, func_args],
                                  writeBedGraph_wrapper,
                                  chrom_names_and_size,
                                  self_=self,
                                  genomeChunkLength=genome_chunk_length,
                                  region=self.region,
                                  blackListFileName=blackListFileName,
                                  numberOfProcessors=self.numberOfProcessors)

        # Determine the sorted order of the temp files
        chrom_order = dict()
        for i, _ in enumerate(chrom_names_and_size):
            chrom_order[_[0]] = i
        res = [[chrom_order[x[0]], x[1], x[2], x[3]] for x in res]
        res.sort()

        # write output for each cluster
        cluster_info = self.clusterInfo
        clusters = cluster_info.cluster.unique().tolist()
        prefix = os.path.splitext(os.path.basename(out_file_prefix))[0]
        for cl in clusters:
            print('Writing output for cluster: {}'.format(cl))
            if pd.isna(cl):
                continue
            # concatenate the coverages
            tmp_out = "/tmp/{}_{}.tmp".format(prefix, cl)
            out_file = open(tmp_out, 'wb')
            for r in res:
                if r[3][cl]:
                    _foo = open(r[3][cl], 'rb')
                    shutil.copyfileobj(_foo, out_file)
                    _foo.close()
                    os.remove(r[3][cl])
            out_file.close()

            ## read back and normalize
            cl_idx = cluster_info.index[pd.Series(cluster_info.cluster == cl)].tolist()
            nCells = float(len(cl_idx))
            out_file = pd.read_csv(tmp_out, sep = "\t", index_col=None, header = None)
            # CPM norm
            if normUsing=='CPM':
                mil_reads_mapped = float(np.sum(out_file[3])) / 1e6
                # per mil counts
                out_file[3] *= 1.0 / (mil_reads_mapped)
            elif normUsing=='Mean':
                # divided by nCells
                out_file[3] *= 1.0/(nCells)

            # out
            bg_out = "{}_{}.bedgraph".format(out_file_prefix, cl)

            out_file.to_csv(bg_out, sep ="\t", index = False, header = False)
            os.remove(tmp_out)
            if format == 'bigwig':
                bedGraphToBigWig(chrom_names_and_size, [bg_out],
                                         "{}_{}.bw".format(out_file_prefix, cl))

    def writeBedGraph_worker(self, chrom, start, end, func_to_call, func_args, bed_regions_list=None):

        r"""Writes a bedgraph based on the read coverage per group of cells, indicated by cluster_info data frame.

        The given func is called to compute the desired bedgraph value
        using the funcArgs

        Parameters
        ----------
        chrom : str
            Chrom name
        start : int
            start coordinate
        end : int
            end coordinate
        func_to_call : str
            function name to be called to convert the list of coverages computed
            for each bam file at each position into a single value. An example
            is a function that takes the ratio between the coverage of two
            bam files.
        func_args : dict
            dict of arguments to pass to `func`.
        smoothLength : int
            Distance in bp for smoothing the coverage per tile.
        bed_regions_list: list
            List of tuples of the form (chrom, start, end)
            corresponding to bed regions to be processed.
            If not bed file was passed to the object constructor
            then this list is empty.

        Returns
        -------
        A list of [chromosome, start, end, temporary file], where the temporary file contains the bedgraph results for the region queried.

        Examples
        --------
        >>> test_path = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        >>> bamFile1 = test_path +  "testA.bam"
        >>> bin_length = 50
        >>> number_of_samples = 0 # overruled by step_size
        >>> func_to_call = scaleCoverage
        >>> funcArgs = {'scaleFactor': 1.0}

        >>> c = WriteBedGraph([bamFile1], bin_length, number_of_samples, stepSize=50)
        >>> tempFile = c.writeBedGraph_worker( '3R', 0, 200, func_to_call, funcArgs)
        >>> f = open(tempFile[3], 'r')
        >>> f.readlines()
        ['3R\t0\t100\t0\n', '3R\t100\t200\t1\n']
        >>> f.close()
        >>> os.remove(tempFile[3])
        """
        if start > end:
                raise NameError("start position ({0}) bigger "
                                "than end position ({1})".format(start, end))
        coverage, _, r = self.count_reads_in_region(chrom, start, end)

        ## get groups (clusters)
        cluster_info = self.clusterInfo
        clusters = cluster_info.cluster.unique().tolist()
        tempfilenames = dict.fromkeys(clusters)
        ## sum up tilecoverage group-wise
        for cl in clusters:
            if pd.isna(cl):
                continue
            _file = open(utilities.getTempFileName(suffix='.bg'), 'w')
            previous_value = None
            line_string = "{}\t{}\t{}\t{:g}\n"
            cl_idx = cluster_info.index[pd.Series(cluster_info.cluster == cl)].tolist()
#            nCells = len(cl_idx)

            for tileIndex in range(coverage.shape[0]):
                ## smoothing disabled for now
                tileCoverage = coverage[tileIndex, :]
                if self.skipZeroOverZero and np.sum(tileCoverage) == 0:
                    continue

                value = func_to_call(np.sum(tileCoverage[cl_idx]), func_args)

                if previous_value is None:
                    writeStart = start + tileIndex * self.binLength
                    writeEnd = min(writeStart + self.binLength, end)
                    previous_value = value

                elif previous_value == value:
                    writeEnd = min(writeEnd + self.binLength, end)

                elif previous_value != value:
                    if not np.isnan(previous_value):
                        _file.write(
                            line_string.format(chrom, writeStart, writeEnd, previous_value))
                    previous_value = value
                    writeStart = writeEnd
                    writeEnd = min(writeStart + self.binLength, end)

            # write remaining value if not a nan
            if previous_value is not None and writeStart != end and not np.isnan(previous_value):
                _file.write(line_string.format(chrom, writeStart,end, previous_value))

            tempfilenames[cl] = _file.name
            _file.close()

        return chrom, start, end, tempfilenames
