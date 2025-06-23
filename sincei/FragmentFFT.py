import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scipy as sp
import numpy as np
import pandas as pd


## fast fouriour transform
def ffttable(selected):
    r"""Computes the FFT of the fragment length distribution

    Parameters
    ----------
    selected : numpy array
        Array of fragment lengths

    Returns
    -------
    pandas dataframe
        Dataframe with two columns. The first column contains the frequencies
        and the second column the power spectrum.

    Examples
    --------

    >>> test = Tester()
    >>> c = CountReadsPerBin([test.bamFile1, test.bamFile2], 50, 4)
    >>> num_reads_per_bin, regionList = c.run()
    >>> selected = num_reads_per_bin[:,0]
    >>> ffttable(selected)
       freq         value
    0   0.0  7.812500e+06
    1   1.0  1.812500e+06
    2   2.0  1.812500e+06
    3   3.0  1.812500e+06
    4   4.0  1.812500e+06
    5   5.0  1.812500
    """
    selected = np.log2(selected + 1)
    selected2 = [y - x for x, y in zip(selected, selected[1:])]
    fragment_fft = sp.fftpack.fft(selected2)
    fragment_psd = np.abs(fragment_fft) ** 2
    fftfreq = sp.fftpack.fftfreq(len(fragment_psd), 10.0)
    d = pd.DataFrame({"freq": fftfreq, "value": fragment_psd})
    d2 = d[d.freq > 0]
    return d2


# get fragment size disribution and fft value from dict{barcode:fragment_size_list}
def fragment_distribution(fragment_len_dict, length_plot):
    r"""Plot fragment length distribution

    Parameters
    ----------
    fragment_len_dict : dict
        Dictionary containing the fragment length for each barcode.

    length_plot : str
        File name to save the plot.

    Returns
    -------
    dict
        Dictionary containing the top 20 frequencies for each barcode.


    Examples
    --------

    >>> test = Tester()
    >>> fragment_len_dict = {'AAACCTGAGAGGTTCT': [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
    ...                      'AAACCTGAGAGGTTCT': [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]}
    >>> fragment_distribution(fragment_len_dict, test.length_plot)
    """
    plt.style.use("classic")
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)

    outdict = dict.fromkeys(fragment_len_dict.keys())
    for barcode in outdict.keys():
        fragment_len = fragment_len_dict[barcode]

        n, bins, patches = ax1.hist(fragment_len, bins=100, color="orange", range=(0, 1000), alpha=0.3)
        # calculating fft
        dflist = []
        for num in range(80, 101, 1):
            dflist.append(ffttable(n[10:num]))
        d2 = pd.concat(dflist, ignore_index=True)
        d2 = d2.sort_values(by=["freq"], ascending=False)
        outdict[barcode] = d2

    # plot fragment sizes
    ax1.xaxis.set_ticks_position("bottom")
    ax1.yaxis.set_ticks_position("left")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    plt.xlabel("Fragment length", fontsize=20)
    plt.ylabel("Fragment count", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax1.set_title("Fragment distribution", fontsize=20)
    plt.savefig(length_plot, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return outdict


## plot the fragment periodicity per barcode using the output of
def fftplot(outdict, plot):
    r"""Computes the periodicity of the fragment distribution

    Parameters
    ----------
    outdict : dict
        Dictionary containing the fragment distribution for each barcode

    plot : str
        File name to save the plot

    Returns
    -------
    dict
        Dictionary containing the periodicity of the fragment distribution for each barcode

    Examples
    --------

    >>> test = Tester()
    >>> d = test.get_fragment_distribution()
    >>> fftplot(d, 'test.png')
    """

    plt.style.use("classic")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    periodicity = dict.fromkeys(outdict.keys())
    for barcode in outdict.keys():
        d2 = outdict[barcode]
        max_y_pos = 1 / d2.loc[d2.value.idxmax(), "freq"]
        mean_value = np.mean(d2.value)

        ## color barcodes where max periodicity is not between 120 and 160, or 240 and 320
        p1 = range(120, 160, 1)
        p2 = range(240, 320, 1)

        mononuc_periodicity = np.mean(d2.loc[[int(x) in p1 for x in 1 / d2.freq], "value"])
        dinuc_periodicity = np.mean(d2.loc[[int(x) in p2 for x in 1 / d2.freq], "value"])

        # if int(mean_value) in p1:
        #    col='red'
        # elif int(mean_value) in p2:
        #    col='blue'
        # else:
        col = "grey"
        ax.plot(1 / d2.freq, 10 * np.log10(d2.value + 1), color=col, alpha=0.2)

        periodicity[barcode] = [mononuc_periodicity, dinuc_periodicity]

    # ax.axvline(max_y_pos,color = 'r',linestyle= '--')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.set_xlim(0, 400)
    ax.set_xlabel("Period (bp)", fontsize=20)
    ax.set_ylabel("Power", fontsize=20)
    ax.set_title("Period of fragment distribution", fontsize=20)

    plt.savefig(plot, dpi=300, bbox_inches="tight")
    return periodicity
