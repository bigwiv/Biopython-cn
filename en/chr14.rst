Chapter 14  Sequence motif analysis using Bio.motifs
====================================================

This chapter gives an overview of the functionality of the
``Bio.motifs`` package included in Biopython. It is intended for people
who are involved in the analysis of sequence motifs, so I’ll assume that
you are familiar with basic notions of motif analysis. In case something
is unclear, please look at Section \ `14.8 <#sec:links>`__ for some
relevant links.

Most of this chapter describes the new ``Bio.motifs`` package included
in Biopython 1.61 onwards, which is replacing the older ``Bio.Motif``
package introduced with Biopython 1.50, which was in turn based on two
older former Biopython modules, ``Bio.AlignAce`` and ``Bio.MEME``. It
provides most of their functionality with a unified motif object
implementation.

Speaking of other libraries, if you are reading this you might be
interested in `TAMO <http://fraenkel.mit.edu/TAMO/>`__, another python
library designed to deal with sequence motifs. It supports more
*de-novo* motif finders, but it is not a part of Biopython and has some
restrictions on commercial use.

14.1  Motif objects
-------------------

Since we are interested in motif analysis, we need to take a look at
``Motif`` objects in the first place. For that we need to import the
Bio.motifs library:

.. code:: verbatim

    >>> from Bio import motifs

and we can start creating our first motif objects. We can either create
a ``Motif`` object from a list of instances of the motif, or we can
obtain a ``Motif`` object by parsing a file from a motif database or
motif finding software.

14.1.1  Creating a motif from instances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we have these instances of a DNA motif:

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> instances = [Seq("TACAA"),
    ...              Seq("TACGC"),
    ...              Seq("TACAC"),
    ...              Seq("TACCC"),
    ...              Seq("AACCC"),
    ...              Seq("AATGC"),
    ...              Seq("AATGC"),
    ...             ]

then we can create a Motif object as follows:

.. code:: verbatim

    >>> m = motifs.create(instances)

The instances are saved in an attribute ``m.instances``, which is
essentially a Python list with some added functionality, as described
below. Printing out the Motif object shows the instances from which it
was constructed:

.. code:: verbatim

    >>> print m
    TACAA
    TACGC
    TACAC
    TACCC
    AACCC
    AATGC
    AATGC
    <BLANKLINE>

The length of the motif defined as the sequence length, which should be
the same for all instances:

.. code:: verbatim

    >>> len(m)
    5

The Motif object has an attribute ``.counts`` containing the counts of
each nucleotide at each position. Printing this counts matrix shows it
in an easily readable format:

.. code:: verbatim

    >>> print m.counts
            0      1      2      3      4
    A:   3.00   7.00   0.00   2.00   1.00
    C:   0.00   0.00   5.00   2.00   6.00
    G:   0.00   0.00   0.00   3.00   0.00
    T:   4.00   0.00   2.00   0.00   0.00
    <BLANKLINE>

You can access these counts as a dictionary:

.. code:: verbatim

    >>> m.counts['A']
    [3, 7, 0, 2, 1]

but you can also think of it as a 2D array with the nucleotide as the
first dimension and the position as the second dimension:

.. code:: verbatim

    >>> m.counts['T',0]
    4
    >>> m.counts['T',2]
    2
    >>> m.counts['T',3]
    0

You can also directly access columns of the counts matrix

.. code:: verbatim

    >>> m.counts[:,3]
    {'A': 2, 'C': 2, 'T': 0, 'G': 3}

Instead of the nucleotide itself, you can also use the index of the
nucleotide in the sorted letters in the alphabet of the motif:

.. code:: verbatim

    >>> m.alphabet
    IUPACUnambiguousDNA()
    >>> m.alphabet.letters
    'GATC'
    >>> sorted(m.alphabet.letters)
    ['A', 'C', 'G', 'T']
    >>> m.counts['A',:]
    (3, 7, 0, 2, 1)
    >>> m.counts[0,:]
    (3, 7, 0, 2, 1)

The motif has an associated consensus sequence, defined as the sequence
of letters along the positions of the motif for which the largest value
in the corresponding columns of the ``.counts`` matrix is obtained:

.. code:: verbatim

    >>> m.consensus
    Seq('TACGC', IUPACUnambiguousDNA())

as well as an anticonsensus sequence, corresponding to the smallest
values in the columns of the ``.counts`` matrix:

.. code:: verbatim

    >>> m.anticonsensus
    Seq('GGGTG', IUPACUnambiguousDNA())

You can also ask for a degenerate consensus sequence, in which ambiguous
nucleotides are used for positions where there are multiple nucleotides
with high counts:

.. code:: verbatim

    >>> m.degenerate_consensus
    Seq('WACVC', IUPACAmbiguousDNA())

Here, W and R follow the IUPAC nucleotide ambiguity codes: W is either A
or T, and V is A, C, or G [`10 <#cornish1985>`__\ ]. The degenerate
consensus sequence is constructed following the rules specified by
Cavener [`11 <#cavener1987>`__\ ].

We can also get the reverse complement of a motif:

.. code:: verbatim

    >>> r = m.reverse_complement()
    >>> r.consensus
    Seq('GCGTA', IUPACUnambiguousDNA())
    >>> r.degenerate_consensus
    Seq('GBGTW', IUPACAmbiguousDNA())
    >>> print r
    TTGTA
    GCGTA
    GTGTA
    GGGTA
    GGGTT
    GCATT
    GCATT
    <BLANKLINE>

The reverse complement and the degenerate consensus sequence are only
defined for DNA motifs.

14.1.2  Reading motifs
~~~~~~~~~~~~~~~~~~~~~~

Creating motifs from instances by hand is a bit boring, so it’s useful
to have some I/O functions for reading and writing motifs. There are no
really well established standards for storing motifs, but there’s a
couple of formats which are more used than others. The most important
distinction is whether the motif representation is based on instances or
on some version of PWM matrix.

JASPAR
^^^^^^

One of the most popular motif databases
`JASPAR <http://jaspar.genereg.net>`__ stores motifs either as a list of
instances, or as a frequency matrix. As an example, these are the
beginning and ending lines of the JASPAR ``Arnt.sites`` file showing
known binding sites of the mouse helix-loop-helix transcription factor
Arnt:

.. code:: verbatim

    >MA0004 ARNT    1
    CACGTGatgtcctc
    >MA0004 ARNT    2
    CACGTGggaggtac
    >MA0004 ARNT    3
    CACGTGccgcgcgc
    ...
    >MA0004 ARNT    18
    AACGTGacagccctcc
    >MA0004 ARNT    19
    AACGTGcacatcgtcc
    >MA0004 ARNT    20
    aggaatCGCGTGc

The parts of the sequence in capital letters are the motif instances
that were found to align to each other.

We can create a ``Motif`` object from these instances as follows:

.. code:: verbatim

    >>> from Bio import motifs
    >>> arnt = motifs.read(open("Arnt.sites"), "sites")

The instances from which this motif was created is stored in the
``.instances`` property:

.. code:: verbatim

    >>> print arnt.instances[:3]
    [Seq('CACGTG', IUPACUnambiguousDNA()), Seq('CACGTG', IUPACUnambiguousDNA()), Seq('CACGTG', IUPACUnambiguousDNA())]
    >>> for instance in arnt.instances:
    ...     print instance
    ... 
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    CACGTG
    AACGTG
    AACGTG
    AACGTG
    AACGTG
    CGCGTG

The counts matrix of this motif is automatically calculated from the
instances:

.. code:: verbatim

    >>> print arnt.counts
            0      1      2      3      4      5
    A:   4.00  19.00   0.00   0.00   0.00   0.00
    C:  16.00   0.00  20.00   0.00   0.00   0.00
    G:   0.00   1.00   0.00  20.00   0.00  20.00
    T:   0.00   0.00   0.00   0.00  20.00   0.00
    <BLANKLINE>

The JASPAR database also makes motifs available directly as a count
matrix, without the instances from which it was created. For example,
this is the JASPAR file ``SRF.pfm`` containing the count matrix for the
human SRF transcription factor:

.. code:: verbatim

     2  9  0  1 32  3 46  1 43 15  2  2
     1 33 45 45  1  1  0  0  0  1  0  1
    39  2  1  0  0  0  0  0  0  0 44 43
     4  2  0  0 13 42  0 45  3 30  0  0

We can create a motif for this count matrix as follows:

.. code:: verbatim

    >>> srf = motifs.read(open("SRF.pfm"),"pfm")
    >>> print srf.counts
            0      1      2      3      4      5      6      7      8      9     10     11
    A:   2.00   9.00   0.00   1.00  32.00   3.00  46.00   1.00  43.00  15.00   2.00   2.00
    C:   1.00  33.00  45.00  45.00   1.00   1.00   0.00   0.00   0.00   1.00   0.00   1.00
    G:  39.00   2.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  44.00  43.00
    T:   4.00   2.00   0.00   0.00  13.00  42.00   0.00  45.00   3.00  30.00   0.00   0.00
    <BLANKLINE>

As this motif was created from the counts matrix directly, it has no
instances associated with it:

.. code:: verbatim

    >>> print srf.instances
    None

We can now ask for the consensus sequence of these two motifs:

.. code:: verbatim

    >>> print arnt.counts.consensus
    CACGTG
    >>> print srf.counts.consensus
    GCCCATATATGG

MEME
^^^^

MEME [`12 <#bailey1994>`__\ ] is a tool for discovering motifs in a
group of related DNA or protein sequences. It takes as input a group of
DNA or protein sequences and outputs as many motifs as requested.
Therefore, in contrast to JASPAR files, MEME output files typically
contain multiple motifs. This is an example.

At the top of an output file generated by MEME shows some background
information about the MEME and the version of MEME used:

.. code:: verbatim

    ********************************************************************************
    MEME - Motif discovery tool
    ********************************************************************************
    MEME version 3.0 (Release date: 2004/08/18 09:07:01)
    ...

Further down, the input set of training sequences is recapitulated:

.. code:: verbatim

    ********************************************************************************
    TRAINING SET
    ********************************************************************************
    DATAFILE= INO_up800.s
    ALPHABET= ACGT
    Sequence name            Weight Length  Sequence name            Weight Length
    -------------            ------ ------  -------------            ------ ------
    CHO1                     1.0000    800  CHO2                     1.0000    800
    FAS1                     1.0000    800  FAS2                     1.0000    800
    ACC1                     1.0000    800  INO1                     1.0000    800
    OPI3                     1.0000    800
    ********************************************************************************

and the exact command line that was used:

.. code:: verbatim

    ********************************************************************************
    COMMAND LINE SUMMARY
    ********************************************************************************
    This information can also be useful in the event you wish to report a
    problem with the MEME software.

    command: meme -mod oops -dna -revcomp -nmotifs 2 -bfile yeast.nc.6.freq INO_up800.s
    ...

Next is detailed information on each motif that was found:

.. code:: verbatim

    ********************************************************************************
    MOTIF  1        width =   12   sites =   7   llr = 95   E-value = 2.0e-001
    ********************************************************************************
    --------------------------------------------------------------------------------
            Motif 1 Description
    --------------------------------------------------------------------------------
    Simplified        A  :::9:a::::3:
    pos.-specific     C  ::a:9:11691a
    probability       G  ::::1::94:4:
    matrix            T  aa:1::9::11:

To parse this file (stored as ``meme.dna.oops.txt``), use

.. code:: verbatim

    >>> handle = open("meme.dna.oops.txt")
    >>> record = motifs.parse(handle, "meme")
    >>> handle.close()

The ``motifs.parse`` command reads the complete file directly, so you
can close the file after calling ``motifs.parse``. The header
information is stored in attributes:

.. code:: verbatim

    >>> record.version
    '3.0'
    >>> record.datafile
    'INO_up800.s'
    >>> record.command
    'meme -mod oops -dna -revcomp -nmotifs 2 -bfile yeast.nc.6.freq INO_up800.s'
    >>> record.alphabet
    IUPACUnambiguousDNA()
    >>> record.sequences
    ['CHO1', 'CHO2', 'FAS1', 'FAS2', 'ACC1', 'INO1', 'OPI3']

The record is an object of the ``Bio.motifs.meme.Record`` class. The
class inherits from list, and you can think of ``record`` as a list of
Motif objects:

.. code:: verbatim

    >>> len(record)
    2
    >>> motif = record[0]
    >>> print motif.consensus
    TTCACATGCCGC
    >>> print motif.degenerate_consensus
    TTCACATGSCNC

In addition to these generic motif attributes, each motif also stores
its specific information as calculated by MEME. For example,

.. code:: verbatim

    >>> motif.num_occurrences
    7
    >>> motif.length
    12
    >>> evalue = motif.evalue
    >>> print "%3.1g" % evalue
    0.2
    >>> motif.name
    'Motif 1'

In addition to using an index into the record, as we did above, you can
also find it by its name:

.. code:: verbatim

    >>> motif = record['Motif 1']

Each motif has an attribute ``.instances`` with the sequence instances
in which the motif was found, providing some information on each
instance:

.. code:: verbatim

    >>> len(motif.instances)
    7
    >>> motif.instances[0]
    Instance('TTCACATGCCGC', IUPACUnambiguousDNA())
    >>> motif.instances[0].motif_name
    'Motif 1'
    >>> motif.instances[0].sequence_name
    'INO1'
    >>> motif.instances[0].start
    620
    >>> motif.instances[0].strand
    '-'
    >>> motif.instances[0].length
    12
    >>> pvalue = motif.instances[0].pvalue

.. code:: verbatim

    >>> print "%5.3g" % pvalue
    1.85e-08

MAST
^^^^

TRANSFAC
^^^^^^^^

TRANSFAC is a manually curated database of transcription factors,
together with their genomic binding sites and DNA binding profiles
[`27 <#matys2003>`__\ ]. While the file format used in the TRANSFAC
database is nowadays also used by others, we will refer to it as the
TRANSFAC file format.

A minimal file in the TRANSFAC format looks as follows:

.. code:: verbatim

    ID  motif1
    P0      A      C      G      T
    01      1      2      2      0      S
    02      2      1      2      0      R
    03      3      0      1      1      A
    04      0      5      0      0      C
    05      5      0      0      0      A
    06      0      0      4      1      G
    07      0      1      4      0      G
    08      0      0      0      5      T
    09      0      0      5      0      G
    10      0      1      2      2      K
    11      0      2      0      3      Y
    12      1      0      3      1      G
    //

This file shows the frequency matrix of motif ``motif1`` of 12
nucleotides. In general, one file in the TRANSFAC format can contain
multiple motifs. For example, this is the contents of the example
TRANSFAC file ``transfac.dat``:

.. code:: verbatim

    VV  EXAMPLE January 15, 2013
    XX
    //
    ID  motif1
    P0      A      C      G      T
    01      1      2      2      0      S
    02      2      1      2      0      R
    03      3      0      1      1      A
    ...
    11      0      2      0      3      Y
    12      1      0      3      1      G
    //
    ID  motif2
    P0      A      C      G      T
    01      2      1      2      0      R
    02      1      2      2      0      S
    ...
    09      0      0      0      5      T
    10      0      2      0      3      Y
    //

To parse a TRANSFAC file, use

.. code:: verbatim

    >>> handle = open("transfac.dat")
    >>> record = motifs.parse(handle, "TRANSFAC")
    >>> handle.close()

The overall version number, if available, is stored as
``record.version``:

.. code:: verbatim

    >>> record.version
    'EXAMPLE January 15, 2013'

Each motif in ``record`` is in instance of the
``Bio.motifs.transfac.Motif`` class, which inherits both from the
``Bio.motifs.Motif`` class and from a Python dictionary. The dictionary
uses the two-letter keys to store any additional information about the
motif:

.. code:: verbatim

    >>> motif = record[0]
    >>> motif.degenerate_consensus # Using the Bio.motifs.Motif method
    Seq('SRACAGGTGKYG', IUPACAmbiguousDNA())
    >>> motif['ID'] # Using motif as a dictionary
    'motif1'

TRANSFAC files are typically much more elaborate than this example,
containing lots of additional information about the motif. Table
`14.1.2 <#table:transfaccodes>`__ lists the two-letter field codes that
are commonly found in TRANSFAC files:

    --------------

    +-------------------------------------------------------+
    | Table 14.1: Fields commonly found in TRANSFAC files   |
    +-------------------------------------------------------+

    +----------+---------------------------------------------------+
    | ``AC``   | Accession number                                  |
    +----------+---------------------------------------------------+
    | ``AS``   | Accession numbers, secondary                      |
    +----------+---------------------------------------------------+
    | ``BA``   | Statistical basis                                 |
    +----------+---------------------------------------------------+
    | ``BF``   | Binding factors                                   |
    +----------+---------------------------------------------------+
    | ``BS``   | Factor binding sites underlying the matrix        |
    +----------+---------------------------------------------------+
    | ``CC``   | Comments                                          |
    +----------+---------------------------------------------------+
    | ``CO``   | Copyright notice                                  |
    +----------+---------------------------------------------------+
    | ``DE``   | Short factor description                          |
    +----------+---------------------------------------------------+
    | ``DR``   | External databases                                |
    +----------+---------------------------------------------------+
    | ``DT``   | Date created/updated                              |
    +----------+---------------------------------------------------+
    | ``HC``   | Subfamilies                                       |
    +----------+---------------------------------------------------+
    | ``HP``   | Superfamilies                                     |
    +----------+---------------------------------------------------+
    | ``ID``   | Identifier                                        |
    +----------+---------------------------------------------------+
    | ``NA``   | Name of the binding factor                        |
    +----------+---------------------------------------------------+
    | ``OC``   | Taxonomic classification                          |
    +----------+---------------------------------------------------+
    | ``OS``   | Species/Taxon                                     |
    +----------+---------------------------------------------------+
    | ``OV``   | Older version                                     |
    +----------+---------------------------------------------------+
    | ``PV``   | Preferred version                                 |
    +----------+---------------------------------------------------+
    | ``TY``   | Type                                              |
    +----------+---------------------------------------------------+
    | ``XX``   | Empty line; these are not stored in the Record.   |
    +----------+---------------------------------------------------+

    --------------

Each motif also has an attribute ``.references`` containing the
references associated with the motif, using these two-letter keys:

    --------------

    +-----------------------------------------------------------------+
    | Table 14.2: Fields used to store references in TRANSFAC files   |
    +-----------------------------------------------------------------+

    +----------+---------------------+
    | ``RN``   | Reference number    |
    +----------+---------------------+
    | ``RA``   | Reference authors   |
    +----------+---------------------+
    | ``RL``   | Reference data      |
    +----------+---------------------+
    | ``RT``   | Reference title     |
    +----------+---------------------+
    | ``RX``   | PubMed ID           |
    +----------+---------------------+

    --------------

Printing the motifs writes them out in their native TRANSFAC format:

.. code:: verbatim

    >>> print record
    VV  EXAMPLE January 15, 2013
    XX
    //
    ID  motif1
    XX
    P0      A      C      G      T
    01      1      2      2      0      S
    02      2      1      2      0      R
    03      3      0      1      1      A
    04      0      5      0      0      C
    05      5      0      0      0      A
    06      0      0      4      1      G
    07      0      1      4      0      G
    08      0      0      0      5      T
    09      0      0      5      0      G
    10      0      1      2      2      K
    11      0      2      0      3      Y
    12      1      0      3      1      G
    XX
    //
    ID  motif2
    XX
    P0      A      C      G      T
    01      2      1      2      0      R
    02      1      2      2      0      S
    03      0      5      0      0      C
    04      3      0      1      1      A
    05      0      0      4      1      G
    06      5      0      0      0      A
    07      0      1      4      0      G
    08      0      0      5      0      G
    09      0      0      0      5      T
    10      0      2      0      3      Y
    XX
    //
    <BLANKLINE>

You can export the motifs in the TRANSFAC format by capturing this
output in a string and saving it in a file:

.. code:: verbatim

    >>> text = str(record)
    >>> handle = open("mytransfacfile.dat", 'w')
    >>> handle.write(text)
    >>> handle.close()

14.1.3  Writing motifs
~~~~~~~~~~~~~~~~~~~~~~

Speaking of exporting, let’s look at export functions in general. To
export a motif in the JASPAR ``.pfm`` format, use

.. code:: verbatim

    >>> print m.format("pfm")
    3       7       0       2       1
    0       0       5       2       6
    0       0       0       3       0
    4       0       2       0       0
    <BLANKLINE>

To write the motif in a TRANSFAC-like matrix format, use

.. code:: verbatim

    >>> print m.format("transfac")
    P0      A      C      G      T
    01      3      0      0      4      W
    02      7      0      0      0      A
    03      0      5      0      2      C
    04      2      2      3      0      V
    05      1      6      0      0      C
    XX
    //
    <BLANKLINE>

To write out multiple motifs, you can use ``motifs.write``. This
function can be used regardless of whether the motifs originated from a
TRANSFAC file. For example,

.. code:: verbatim

    >>> two_motifs = [arnt, srf]
    >>> print motifs.write(two_motifs, 'transfac')
    P0      A      C      G      T
    01      4     16      0      0      C
    02     19      0      1      0      A
    03      0     20      0      0      C
    04      0      0     20      0      G
    05      0      0      0     20      T
    06      0      0     20      0      G
    XX
    //
    P0      A      C      G      T
    01      2      1     39      4      G
    02      9     33      2      2      C
    03      0     45      1      0      C
    04      1     45      0      0      C
    05     32      1      0     13      A
    06      3      1      0     42      T
    07     46      0      0      0      A
    08      1      0      0     45      T
    09     43      0      0      3      A
    10     15      1      0     30      T
    11      2      0     44      0      G
    12      2      1     43      0      G
    XX
    //
    <BLANKLINE>

14.1.4  Creating a sequence logo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If we have internet access, we can create a
`weblogo <http://weblogo.berkeley.edu>`__:

.. code:: verbatim

    >>> arnt.weblogo("Arnt.png")

We should get our logo saved as a PNG in the specified file.

14.2  Position-Weight Matrices
------------------------------

The ``.counts`` attribute of a Motif object shows how often each
nucleotide appeared at each position along the alignment. We can
normalize this matrix by dividing by the number of instances in the
alignment, resulting in the probability of each nucleotide at each
position along the alignment. We refer to these probabilities as the
position-weight matrix. However, beware that in the literature this term
may also be used to refer to the position-specific scoring matrix, which
we discuss below.

Usually, pseudocounts are added to each position before normalizing.
This avoids overfitting of the position-weight matrix to the limited
number of motif instances in the alignment, and can also prevent
probabilities from becoming zero. To add a fixed pseudocount to all
nucleotides at all positions, specify a number for the ``pseudocounts``
argument:

.. code:: verbatim

    >>> pwm = m.counts.normalize(pseudocounts=0.5)
    >>> print pwm
            0      1      2      3      4
    A:   0.39   0.83   0.06   0.28   0.17
    C:   0.06   0.06   0.61   0.28   0.72
    G:   0.06   0.06   0.06   0.39   0.06
    T:   0.50   0.06   0.28   0.06   0.06
    <BLANKLINE>

Alternatively, ``pseudocounts`` can be a dictionary specifying the
pseudocounts for each nucleotide. For example, as the GC content of the
human genome is about 40%, you may want to choose the pseudocounts
accordingly:

.. code:: verbatim

    >>> pwm = m.counts.normalize(pseudocounts={'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6})
    >>> print pwm
            0      1      2      3      4
    A:   0.40   0.84   0.07   0.29   0.18
    C:   0.04   0.04   0.60   0.27   0.71
    G:   0.04   0.04   0.04   0.38   0.04
    T:   0.51   0.07   0.29   0.07   0.07
    <BLANKLINE>

The position-weight matrix has its own methods to calculate the
consensus, anticonsensus, and degenerate consensus sequences:

.. code:: verbatim

    >>> pwm.consensus
    Seq('TACGC', IUPACUnambiguousDNA())
    >>> pwm.anticonsensus
    Seq('GGGTG', IUPACUnambiguousDNA())
    >>> pwm.degenerate_consensus
    Seq('WACNC', IUPACAmbiguousDNA())

Note that due to the pseudocounts, the degenerate consensus sequence
calculated from the position-weight matrix is slightly different from
the degenerate consensus sequence calculated from the instances in the
motif:

.. code:: verbatim

    >>> m.degenerate_consensus
    Seq('WACVC', IUPACAmbiguousDNA())

The reverse complement of the position-weight matrix can be calculated
directly from the ``pwm``:

.. code:: verbatim

    >>> rpwm = pwm.reverse_complement()
    >>> print rpwm
            0      1      2      3      4
    A:   0.07   0.07   0.29   0.07   0.51
    C:   0.04   0.38   0.04   0.04   0.04
    G:   0.71   0.27   0.60   0.04   0.04
    T:   0.18   0.29   0.07   0.84   0.40
    <BLANKLINE>

14.3  Position-Specific Scoring Matrices
----------------------------------------

Using the background distribution and PWM with pseudo-counts added, it’s
easy to compute the log-odds ratios, telling us what are the log odds of
a particular symbol to be coming from a motif against the background. We
can use the ``.log_odds()`` method on the position-weight matrix:

.. code:: verbatim

    >>> pssm = pwm.log_odds()
    >>> print pssm
            0      1      2      3      4
    A:   0.68   1.76  -1.91   0.21  -0.49
    C:  -2.49  -2.49   1.26   0.09   1.51
    G:  -2.49  -2.49  -2.49   0.60  -2.49
    T:   1.03  -1.91   0.21  -1.91  -1.91
    <BLANKLINE>

Here we can see positive values for symbols more frequent in the motif
than in the background and negative for symbols more frequent in the
background. 0.0 means that it’s equally likely to see a symbol in the
background and in the motif.

This assumes that A, C, G, and T are equally likely in the background.
To calculate the position-specific scoring matrix against a background
with unequal probabilities for A, C, G, T, use the ``background``
argument. For example, against a background with a 40% GC content, use

.. code:: verbatim

    >>> background = {'A':0.3,'C':0.2,'G':0.2,'T':0.3}
    >>> pssm = pwm.log_odds(background)
    >>> print pssm
            0      1      2      3      4
    A:   0.42   1.49  -2.17  -0.05  -0.75
    C:  -2.17  -2.17   1.58   0.42   1.83
    G:  -2.17  -2.17  -2.17   0.92  -2.17
    T:   0.77  -2.17  -0.05  -2.17  -2.17
    <BLANKLINE>

The maximum and minimum score obtainable from the PSSM are stored in the
``.max`` and ``.min`` properties:

.. code:: verbatim

    >>> print "%4.2f" % pssm.max
    6.59
    >>> print "%4.2f" % pssm.min
    -10.85

The mean and standard deviation of the PSSM scores with respect to a
specific background are calculated by the ``.mean`` and ``.std``
methods.

.. code:: verbatim

    >>> mean = pssm.mean(background)
    >>> std = pssm.std(background)
    >>> print "mean = %0.2f, standard deviation = %0.2f" % (mean, std)
    mean = 3.21, standard deviation = 2.59

A uniform background is used if ``background`` is not specified. The
mean is particularly important, as its value is equal to the
Kullback-Leibler divergence or relative entropy, and is a measure for
the information content of the motif compared to the background. As in
Biopython the base-2 logarithm is used in the calculation of the
log-odds scores, the information content has units of bits.

The ``.reverse_complement``, ``.consensus``, ``.anticonsensus``, and
``.degenerate_consensus`` methods can be applied directly to PSSM
objects.

14.4  Searching for instances
-----------------------------

The most frequent use for a motif is to find its instances in some
sequence. For the sake of this section, we will use an artificial
sequence like this:

.. code:: verbatim

    >>> test_seq=Seq("TACACTGCATTACAACCCAAGCATTA",m.alphabet)
    >>> len(test_seq)
    26

14.4.1  Searching for exact matches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to find instances, is to look for exact matches of the
true instances of the motif:

.. code:: verbatim

    >>> for pos,seq in m.instances.search(test_seq):
    ...     print pos, seq
    ... 
    0 TACAC
    10 TACAA
    13 AACCC

We can do the same with the reverse complement (to find instances on the
complementary strand):

.. code:: verbatim

    >>> for pos,seq in r.instances.search(test_seq):
    ...     print pos, seq
    ... 
    6 GCATT
    20 GCATT

14.4.2  Searching for matches using the PSSM score
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It’s just as easy to look for positions, giving rise to high log-odds
scores against our motif:

.. code:: verbatim

    >>> for position, score in pssm.search(test_seq, threshold=3.0):
    ...     print "Position %d: score = %5.3f" % (position, score)
    ... 
    Position 0: score = 5.622
    Position -20: score = 4.601
    Position 10: score = 3.037
    Position 13: score = 5.738
    Position -6: score = 4.601

The negative positions refer to instances of the motif found on the
reverse strand of the test sequence, and follow the Python convention on
negative indices. Therefore, the instance of the motif at ``pos`` is
located at ``test_seq[pos:pos+len(m)]`` both for positive and for
negative values of ``pos``.

You may notice the threshold parameter, here set arbitrarily to 3.0.
This is in *log*\ :sub:`2`, so we are now looking only for words, which
are eight times more likely to occur under the motif model than in the
background. The default threshold is 0.0, which selects everything that
looks more like the motif than the background.

You can also calculate the scores at all positions along the sequence:

.. code:: verbatim

    >>> pssm.calculate(test_seq)
    array([  5.62230396,  -5.6796999 ,  -3.43177247,   0.93827754,
            -6.84962511,  -2.04066086, -10.84962463,  -3.65614533,
            -0.03370807,  -3.91102552,   3.03734159,  -2.14918518,
            -0.6016975 ,   5.7381525 ,  -0.50977498,  -3.56422281,
            -8.73414803,  -0.09919716,  -0.6016975 ,  -2.39429784,
           -10.84962463,  -3.65614533], dtype=float32)

In general, this is the fastest way to calculate PSSM scores. The scores
returned by ``pssm.calculate`` are for the forward strand only. To
obtain the scores on the reverse strand, you can take the reverse
complement of the PSSM:

.. code:: verbatim

    >>> rpssm = pssm.reverse_complement()
    >>> rpssm.calculate(test_seq)
    array([ -9.43458748,  -3.06172252,  -7.18665981,  -7.76216221,
            -2.04066086,  -4.26466274,   4.60124254,  -4.2480607 ,
            -8.73414803,  -2.26503372,  -6.49598789,  -5.64668512,
            -8.73414803, -10.84962463,  -4.82356262,  -4.82356262,
            -5.64668512,  -8.73414803,  -4.15613794,  -5.6796999 ,
             4.60124254,  -4.2480607 ], dtype=float32)

14.4.3  Selecting a score threshold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to use a less arbitrary way of selecting thresholds, you can
explore the distribution of PSSM scores. Since the space for a score
distribution grows exponentially with motif length, we are using an
approximation with a given precision to keep computation cost
manageable:

.. code:: verbatim

    >>> distribution = pssm.distribution(background=background, precision=10**4)

The ``distribution`` object can be used to determine a number of
different thresholds. We can specify the requested false-positive rate
(probability of “finding” a motif instance in background generated
sequence):

.. code:: verbatim

    >>> threshold = distribution.threshold_fpr(0.01)
    >>> print "%5.3f" % threshold
    4.009

or the false-negative rate (probability of “not finding” an instance
generated from the motif):

.. code:: verbatim

    >>> threshold = distribution.threshold_fnr(0.1)
    >>> print "%5.3f" % threshold
    -0.510

or a threshold (approximately) satisfying some relation between the
false-positive rate and the false-negative rate (fnr/fpr≃ *t*):

.. code:: verbatim

    >>> threshold = distribution.threshold_balanced(1000)
    >>> print "%5.3f" % threshold
    6.241

or a threshold satisfying (roughly) the equality between the
false-positive rate and the −\ *log* of the information content (as used
in patser software by Hertz and Stormo):

.. code:: verbatim

    >>> threshold = distribution.threshold_patser()
    >>> print "%5.3f" % threshold
    0.346

For example, in case of our motif, you can get the threshold giving you
exactly the same results (for this sequence) as searching for instances
with balanced threshold with rate of 1000.

.. code:: verbatim

    >>> threshold = distribution.threshold_fpr(0.01)
    >>> print "%5.3f" % threshold
    4.009
    >>> for position, score in pssm.search(test_seq,threshold=threshold):
    ...     print "Position %d: score = %5.3f" % (position, score)
    ... 
    Position 0: score = 5.622
    Position -20: score = 4.601
    Position 13: score = 5.738
    Position -6: score = 4.601

14.5  Each motif object has an associated Position-Specific Scoring Matrix
--------------------------------------------------------------------------

To facilitate searching for potential TFBSs using PSSMs, both the
position-weight matrix and the position-specific scoring matrix are
associated with each motif. Using the Arnt motif as an example:

.. code:: verbatim

    >>> from Bio import motifs
    >>> handle = open("Arnt.sites")
    >>> motif = motifs.read(handle, 'sites')
    >>> print motif.counts
            0      1      2      3      4      5
    A:   4.00  19.00   0.00   0.00   0.00   0.00
    C:  16.00   0.00  20.00   0.00   0.00   0.00
    G:   0.00   1.00   0.00  20.00   0.00  20.00
    T:   0.00   0.00   0.00   0.00  20.00   0.00
    <BLANKLINE>
    >>> print motif.pwm
            0      1      2      3      4      5
    A:   0.20   0.95   0.00   0.00   0.00   0.00
    C:   0.80   0.00   1.00   0.00   0.00   0.00
    G:   0.00   0.05   0.00   1.00   0.00   1.00
    T:   0.00   0.00   0.00   0.00   1.00   0.00
    <BLANKLINE>

.. code:: verbatim

    >>> print motif.pssm
            0      1      2      3      4      5
    A:  -0.32   1.93   -inf   -inf   -inf   -inf
    C:   1.68   -inf   2.00   -inf   -inf   -inf
    G:   -inf  -2.32   -inf   2.00   -inf   2.00
    T:   -inf   -inf   -inf   -inf   2.00   -inf
    <BLANKLINE>

The negative infinities appear here because the corresponding entry in
the frequency matrix is 0, and we are using zero pseudocounts by
default:

.. code:: verbatim

    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.pseudocounts[letter])
    ...
    A: 0.00
    C: 0.00
    G: 0.00
    T: 0.00

If you change the ``.pseudocounts`` attribute, the position-frequency
matrix and the position-specific scoring matrix are recalculated
automatically:

.. code:: verbatim

    >>> motif.pseudocounts = 3.0
    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.pseudocounts[letter])
    ...
    A: 3.00
    C: 3.00
    G: 3.00
    T: 3.00

.. code:: verbatim

    >>> print motif.pwm
            0      1      2      3      4      5
    A:   0.22   0.69   0.09   0.09   0.09   0.09
    C:   0.59   0.09   0.72   0.09   0.09   0.09
    G:   0.09   0.12   0.09   0.72   0.09   0.72
    T:   0.09   0.09   0.09   0.09   0.72   0.09
    <BLANKLINE>

.. code:: verbatim

    >>> print motif.pssm
            0      1      2      3      4      5
    A:  -0.19   1.46  -1.42  -1.42  -1.42  -1.42
    C:   1.25  -1.42   1.52  -1.42  -1.42  -1.42
    G:  -1.42  -1.00  -1.42   1.52  -1.42   1.52
    T:  -1.42  -1.42  -1.42  -1.42   1.52  -1.42
    <BLANKLINE>

You can also set the ``.pseudocounts`` to a dictionary over the four
nucleotides if you want to use different pseudocounts for them. Setting
``motif.pseudocounts`` to ``None`` resets it to its default value of
zero.

The position-specific scoring matrix depends on the background
distribution, which is uniform by default:

.. code:: verbatim

    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.background[letter])
    ...
    A: 0.25
    C: 0.25
    G: 0.25
    T: 0.25

Again, if you modify the background distribution, the position-specific
scoring matrix is recalculated:

.. code:: verbatim

    >>> motif.background = {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2}
    >>> print motif.pssm
            0      1      2      3      4      5
    A:   0.13   1.78  -1.09  -1.09  -1.09  -1.09
    C:   0.98  -1.68   1.26  -1.68  -1.68  -1.68
    G:  -1.68  -1.26  -1.68   1.26  -1.68   1.26
    T:  -1.09  -1.09  -1.09  -1.09   1.85  -1.09
    <BLANKLINE>

Setting ``motif.background`` to ``None`` resets it to a uniform
distribution:

.. code:: verbatim

    >>> motif.background = None
    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.background[letter])
    ...
    A: 0.25
    C: 0.25
    G: 0.25
    T: 0.25

If you set ``motif.background`` equal to a single value, it will be
interpreted as the GC content:

.. code:: verbatim

    >>> motif.background = 0.8
    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.background[letter])
    ...
    A: 0.10
    C: 0.40
    G: 0.40
    T: 0.10

Note that you can now calculate the mean of the PSSM scores over the
background against which it was computed:

.. code:: verbatim

    >>> print "%f" % motif.pssm.mean(motif.background)
    4.703928

as well as its standard deviation:

.. code:: verbatim

    >>> print "%f" % motif.pssm.std(motif.background)
    3.290900

and its distribution:

.. code:: verbatim

    >>> distribution = motif.pssm.distribution(background=motif.background)
    >>> threshold = distribution.threshold_fpr(0.01)
    >>> print "%f" % threshold
    3.854375

Note that the position-weight matrix and the position-specific scoring
matrix are recalculated each time you call ``motif.pwm`` or
``motif.pssm``, respectively. If speed is an issue and you want to use
the PWM or PSSM repeatedly, you can save them as a variable, as in

.. code:: verbatim

    >>> pssm = motif.pssm

14.6  Comparing motifs
----------------------

Once we have more than one motif, we might want to compare them.

Before we start comparing motifs, I should point out that motif
boundaries are usually quite arbitrary. This means we often need to
compare motifs of different lengths, so comparison needs to involve some
kind of alignment. This means we have to take into account two things:

-  alignment of motifs
-  some function to compare aligned motifs

To align the motifs, we use ungapped alignment of PSSMs and substitute
zeros for any missing columns at the beginning and end of the matrices.
This means that effectively we are using the background distribution for
columns missing from the PSSM. The distance function then returns the
minimal distance between motifs, as well as the corresponding offset in
their alignment.

To give an exmaple, let us first load another motif, which is similar to
our test motif ``m``:

.. code:: verbatim

    >>> m_reb1 = motifs.read(open("REB1.pfm"), "pfm")
    >>> m_reb1.consensus
    Seq('GTTACCCGG', IUPACUnambiguousDNA())
    >>> print m_reb1.counts
            0      1      2      3      4      5      6      7      8
    A:  30.00   0.00   0.00 100.00   0.00   0.00   0.00   0.00  15.00
    C:  10.00   0.00   0.00   0.00 100.00 100.00 100.00   0.00  15.00
    G:  50.00   0.00   0.00   0.00   0.00   0.00   0.00  60.00  55.00
    T:  10.00 100.00 100.00   0.00   0.00   0.00   0.00  40.00  15.00
    <BLANKLINE>

To make the motifs comparable, we choose the same values for the
pseudocounts and the background distribution as our motif ``m``:

.. code:: verbatim

    >>> m_reb1.pseudocounts = {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
    >>> m_reb1.background = {'A':0.3,'C':0.2,'G':0.2,'T':0.3}
    >>> pssm_reb1 = m_reb1.pssm
    >>> print pssm_reb1
            0      1      2      3      4      5      6      7      8
    A:   0.00  -5.67  -5.67   1.72  -5.67  -5.67  -5.67  -5.67  -0.97
    C:  -0.97  -5.67  -5.67  -5.67   2.30   2.30   2.30  -5.67  -0.41
    G:   1.30  -5.67  -5.67  -5.67  -5.67  -5.67  -5.67   1.57   1.44
    T:  -1.53   1.72   1.72  -5.67  -5.67  -5.67  -5.67   0.41  -0.97
    <BLANKLINE>

We’ll compare these motifs using the Pearson correlation. Since we want
it to resemble a distance measure, we actually take 1−\ *r*, where *r*
is the Pearson correlation coefficient (PCC):

.. code:: verbatim

    >>> distance, offset = pssm.dist_pearson(pssm_reb1)
    >>> print "distance = %5.3g" % distance
    distance = 0.239
    >>> print offset
    -2

This means that the best PCC between motif ``m`` and ``m_reb1`` is
obtained with the following alignment:

.. code:: verbatim

    m:      bbTACGCbb
    m_reb1: GTTACCCGG

where ``b`` stands for background distribution. The PCC itself is
roughly 1−0.239=0.761.

14.7  *De novo* motif finding
-----------------------------

Currently, Biopython has only limited support for *de novo* motif
finding. Namely, we support running and parsing of AlignAce and MEME.
Since the number of motif finding tools is growing rapidly,
contributions of new parsers are welcome.

14.7.1  MEME
~~~~~~~~~~~~

Let’s assume, you have run MEME on sequences of your choice with your
favorite parameters and saved the output in the file ``meme.out``. You
can retrieve the motifs reported by MEME by running the following piece
of code:

.. code:: verbatim

    >>> from Bio import motifs
    >>> motifsM = motifs.parse(open("meme.out"), "meme")

.. code:: verbatim

    >>> motifsM
    [<Bio.motifs.meme.Motif object at 0xc356b0>]

Besides the most wanted list of motifs, the result object contains more
useful information, accessible through properties with self-explanatory
names:

-  ``.alphabet``
-  ``.datafile``
-  ``.sequence_names``
-  ``.version``
-  ``.command``

The motifs returned by the MEME Parser can be treated exactly like
regular Motif objects (with instances), they also provide some extra
functionality, by adding additional information about the instances.

.. code:: verbatim

    >>> motifsM[0].consensus
    Seq('CTCAATCGTA', IUPACUnambiguousDNA())
    >>> motifsM[0].instances[0].sequence_name
    'SEQ10;'
    >>> motifsM[0].instances[0].start
    3
    >>> motifsM[0].instances[0].strand
    '+'

.. code:: verbatim

    >>> motifsM[0].instances[0].pvalue
    8.71e-07

14.7.2  AlignAce
~~~~~~~~~~~~~~~~

We can do very similar things with the AlignACE program. Assume, you
have your output in the file ``alignace.out``. You can parse your output
with the following code:

.. code:: verbatim

    >>> from Bio import motifs
    >>> motifsA = motifs.parse(open("alignace.out"),"alignace")

Again, your motifs behave as they should:

.. code:: verbatim

    >>> motifsA[0].consensus
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

In fact you can even see, that AlignAce found a very similar motif as
MEME. It is just a longer version of a reverse complement of the MEME
motif:

.. code:: verbatim

    >>> motifsM[0].reverse_complement().consensus
    Seq('TACGATTGAG', IUPACUnambiguousDNA())

If you have AlignAce installed on the same machine, you can also run it
directly from Biopython. A short example of how this can be done is
shown below (other parameters can be specified as keyword parameters):

.. code:: verbatim

    >>> command="/opt/bin/AlignACE"
    >>> input_file="test.fa"
    >>> from Bio.motifs.applications import AlignAceCommandline
    >>> cmd = AlignAceCommandline(cmd=command,input=input_file,gcback=0.6,numcols=10)
    >>> stdout,stderr= cmd()

Since AlignAce prints all of its output to standard output, you can get
to your motifs by parsing the first part of the result:

.. code:: verbatim

    >>> motifs = motifs.parse(stdout,"alignace")

14.8  Useful links
------------------

-  `Sequence motif <http://en.wikipedia.org/wiki/Sequence_motif>`__ in
   wikipedia
-  `PWM <http://en.wikipedia.org/wiki/Position_weight_matrix>`__ in
   wikipedia
-  `Consensus
   sequence <http://en.wikipedia.org/wiki/Consensus_sequence>`__ in
   wikipedia
-  `Comparison of different motif finding
   programs <http://bio.cs.washington.edu/assessment/>`__

14.9  Obsolete Bio.Motif module
-------------------------------

The rest of this chapter above describes the ``Bio.motifs`` package
included in Biopython 1.61 onwards, which is replacing the older
``Bio.Motif`` package introduced with Biopython 1.50, which was in turn
based on two older former Biopython modules, ``Bio.AlignAce`` and
``Bio.MEME``.

To allow for a smooth transition, the older ``Bio.Motif`` package will
be maintained in parallel with its replacement ``Bio.motifs`` at least
two more releases, and at least one year.

14.9.1  Motif objects
~~~~~~~~~~~~~~~~~~~~~

Since we are interested in motif analysis, we need to take a look at
``Motif`` objects in the first place. For that we need to import the
Motif library:

.. code:: verbatim

    >>> from Bio import Motif

and we can start creating our first motif objects. Let’s create a DNA
motif:

.. code:: verbatim

    >>> from Bio.Alphabet import IUPAC
    >>> m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)

This is for now just an empty container, so let’s add some sequences to
our newly created motif:

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> m.add_instance(Seq("TATAA",m.alphabet))
    >>> m.add_instance(Seq("TATTA",m.alphabet))
    >>> m.add_instance(Seq("TATAA",m.alphabet))
    >>> m.add_instance(Seq("TATAA",m.alphabet))

Now we have a full ``Motif`` instance, so we can try to get some basic
information about it. Let’s start with length and consensus sequence:

.. code:: verbatim

    >>> len(m)
    5
    >>> m.consensus()
    Seq('TATAA', IUPACUnambiguousDNA())

In case of DNA motifs, we can also get a reverse complement of a motif:

.. code:: verbatim

    >>> m.reverse_complement().consensus()
    Seq('TTATA', IUPACUnambiguousDNA())
    >>> for i in m.reverse_complement().instances:
    ...     print i
    TTATA
    TAATA
    TTATA
    TTATA

We can also calculate the information content of a motif with a simple
call:

.. code:: verbatim

    >>> print "%0.2f" % m.ic()
    5.27

This gives us a number of bits of information provided by the motif,
which tells us how much differs from background.

The most common representation of a motif is a PWM (Position Weight
Matrix). It summarizes the probabilities of finding any symbol (in this
case nucleotide) in any position of a motif. It can be computed by
calling the ``.pwm()`` method:

.. code:: verbatim

    >>> m.pwm()
    [{'A': 0.05, 'C': 0.05, 'T': 0.85, 'G': 0.05}, 
     {'A': 0.85, 'C': 0.05, 'T': 0.05, 'G': 0.05}, 
     {'A': 0.05, 'C': 0.05, 'T': 0.85, 'G': 0.05}, 
     {'A': 0.65, 'C': 0.05, 'T': 0.25, 'G': 0.05}, 
     {'A': 0.85, 'C': 0.05, 'T': 0.05, 'G': 0.05}]

The probabilities in the motif’s PWM are based on the counts in the
instances, but we can see, that even though there were no Gs and no Cs
in the instances, we still have non-zero probabilities assigned to them.
These come from pseudo-counts which are, roughly speaking, a commonly
used way to acknowledge the incompleteness of our knowledge and avoid
technical problems with calculating logarithms of 0.

We can control the way that pseudo-counts are added with two properties
of Motif objects ``.background`` is the probability distribution over
all symbols in the alphabet that we assume represents background,
non-motif sequences (usually based on the GC content of the respective
genome). It is by default set to a uniform distribution upon creation of
a motif:

.. code:: verbatim

    >>> m.background  
    {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25}

The other parameter is ``.beta``, which states the amount of
pseudo-counts we should add to the PWM. By default it is set to 1.0,

.. code:: verbatim

    >>> m.beta
    1.0

so that the total input of pseudo-counts is equal to that of one
instance.

Using the background distribution and pwm with pseudo-counts added, it’s
easy to compute the log-odds ratios, telling us what are the log odds of
a particular symbol to be coming from a motif against the background. We
can use the ``.log_odds()`` method:

.. code:: verbatim

     >>> m.log_odds() 
    [{'A': -2.3219280948873622, 
      'C': -2.3219280948873622, 
      'T': 1.7655347463629771, 
      'G': -2.3219280948873622}, 
     {'A': 1.7655347463629771, 
      'C': -2.3219280948873622, 
      'T': -2.3219280948873622, 
      'G': -2.3219280948873622}, 
     {'A': -2.3219280948873622, 
      'C': -2.3219280948873622, 
      'T': 1.7655347463629771, 
      'G': -2.3219280948873622}, 
     {'A': 1.3785116232537298, 
      'C': -2.3219280948873622, 
      'T': 0.0, 
      'G': -2.3219280948873622}, 
     {'A': 1.7655347463629771, 
      'C': -2.3219280948873622, 
      'T': -2.3219280948873622, 
      'G': -2.3219280948873622}
    ]

Here we can see positive values for symbols more frequent in the motif
than in the background and negative for symbols more frequent in the
background. 0.0 means that it’s equally likely to see a symbol in
background and in the motif (e.g. ‘T’ in the second-last position).

14.9.1.1  Reading and writing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating motifs from instances by hand is a bit boring, so it’s useful
to have some I/O functions for reading and writing motifs. There are no
really well established standards for storing motifs, but there’s a
couple of formats which are more used than others. The most important
distinction is whether the motif representation is based on instances or
on some version of PWM matrix. On of the most popular motif databases
`JASPAR <http://jaspar.genereg.net>`__ stores motifs in both formats, so
let’s look at how we can import JASPAR motifs from instances:

.. code:: verbatim

    >>> from Bio import Motif
    >>> arnt = Motif.read(open("Arnt.sites"),"jaspar-sites")

and from a count matrix:

.. code:: verbatim

    >>> srf = Motif.read(open("SRF.pfm"),"jaspar-pfm")

The ``arnt`` and ``srf`` motifs can both do the same things for us, but
they use different internal representations of the motif. We can tell
that by inspecting the ``has_counts`` and ``has_instances`` properties:

.. code:: verbatim

    >>> arnt.has_instances
    True
    >>> srf.has_instances
    False
    >>> srf.has_counts
    True

.. code:: verbatim

    >>> srf.counts
    {'A': [2, 9, 0, 1, 32, 3, 46, 1, 43, 15, 2, 2],
     'C': [1, 33, 45, 45, 1, 1, 0, 0, 0, 1, 0, 1],
     'G': [39, 2, 1, 0, 0, 0, 0, 0, 0, 0, 44, 43],
     'T': [4, 2, 0, 0, 13, 42, 0, 45, 3, 30, 0, 0]}

There are conversion functions, which can help us convert between
different representations:

.. code:: verbatim

    >>> arnt.make_counts_from_instances()
    {'A': [8, 38, 0, 0, 0, 0],
     'C': [32, 0, 40, 0, 0, 0],
     'G': [0, 2, 0, 40, 0, 40],
     'T': [0, 0, 0, 0, 40, 0]}

    >>> srf.make_instances_from_counts()
    [Seq('GGGAAAAAAAGG', IUPACUnambiguousDNA()),
     Seq('GGCCAAATAAGG', IUPACUnambiguousDNA()),
     Seq('GACCAAATAAGG', IUPACUnambiguousDNA()),
    ....

The important thing to remember here is that the method
``make_instances_from_counts()`` creates fake instances, because usually
there are very many possible sets of instances which give rise to the
same pwm, and if we have only the count matrix, we cannot reconstruct
the original one. This does not make any difference if we are using the
PWM as the representation of the motif, but one should be careful with
exporting instances from count-based motifs.

Speaking of exporting, let’s look at export functions. We can export to
fasta:

.. code:: verbatim

    >>> print m.format("fasta")
    >instance0
    TATAA
    >instance1
    TATTA
    >instance2
    TATAA
    >instance3
    TATAA

or to TRANSFAC-like matrix format (used by some motif processing
software)

.. code:: verbatim

    >>> print m.format("transfac")
    XX
    TY Motif
    ID 
    BF undef
    P0 G A T C
    01 0 0 4 0
    02 0 4 0 0
    03 0 0 4 0
    04 0 3 1 0
    05 0 4 0 0
    XX

Finally, if we have internet access, we can create a
`weblogo <http://weblogo.berkeley.edu>`__:

.. code:: verbatim

    >>> arnt.weblogo("Arnt.png")

We should get our logo saved as a png in the specified file.

14.9.2  Searching for instances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most frequent use for a motif is to find its instances in some
sequence. For the sake of this section, we will use an artificial
sequence like this:

.. code:: verbatim

    test_seq=Seq("TATGATGTAGTATAATATAATTATAA",m.alphabet)

The simplest way to find instances, is to look for exact matches of the
true instances of the motif:

.. code:: verbatim

    >>> for pos,seq in m.search_instances(test_seq):
    ...     print pos,seq.tostring()
    ... 
    10 TATAA
    15 TATAA
    21 TATAA

We can do the same with the reverse complement (to find instances on the
complementary strand):

.. code:: verbatim

    >>> for pos,seq in m.reverse_complement().search_instances(test_seq):
    ...     print pos,seq.tostring()
    ... 
    12 TAATA
    20 TTATA

It’s just as easy to look for positions, giving rise to high log-odds
scores against our motif:

.. code:: verbatim

    >>> for pos,score in m.search_pwm(test_seq,threshold=5.0):
    ...     print pos,score
    ... 
    10 8.44065060871
    -12 7.06213898545
    15 8.44065060871
    -20 8.44065060871
    21 8.44065060871

You may notice the threshold parameter, here set arbitrarily to 5.0.
This is in *log*\ :sub:`2`, so we are now looking only for words, which
are 32 times more likely to occur under the motif model than in the
background. The default threshold is 0.0, which selects everything that
looks more like the motif than the background.

If you want to use a less arbitrary way of selecting thresholds, you can
explore the ``Motif.score_distribution`` class implementing an
distribution of scores for a given motif. Since the space for a score
distribution grows exponentially with motif length, we are using an
approximation with a given precision to keep computation cost
manageable:

.. code:: verbatim

    >>> sd = Motif.score_distribution(m,precision=10**4)

The sd object can be used to determine a number of different thresholds.

We can specify the requested false-positive rate (probability of
“finding” a motif instance in background generated sequence):

.. code:: verbatim

    >>> sd.threshold_fpr(0.01)
    4.3535838726139886

or the false-negative rate (probability of “not finding” an instance
generated from the motif):

.. code:: verbatim

    >>> sd.threshold_fnr(0.1)
    0.26651713652234044

or a threshold (approximately) satisfying some relation between fpr and
fnr *fnr*/*fpr*\ ≃ *t*:

.. code:: verbatim

    >>> sd.threshold_balanced(1000)
    8.4406506087056368

or a threshold satisfying (roughly) the equality between the
false-positive rate and the −\ *log* of the information content (as used
in patser software by Hertz and Stormo).

For example, in case of our motif, you can get the threshold giving you
exactly the same results (for this sequence) as searching for instances
with balanced threshold with rate of 1000.

.. code:: verbatim

    >>> for pos,score in m.search_pwm(test_seq,threshold=sd.threshold_balanced(1000)):
    ...     print pos,score
    ... 
    10 8.44065060871
    15 8.44065060871
    -20 8.44065060871
    21 8.44065060871

14.9.3  Comparing motifs
~~~~~~~~~~~~~~~~~~~~~~~~

Once we have more than one motif, we might want to compare them. For
that, we have currently three different methods of ``Bio.Motif``
objects.

Before we start comparing motifs, I should point out that motif
boundaries are usually quite arbitrary. This means, that we often need
to compare motifs of different lengths, so comparison needs to involve
some kind of alignment. This means, that we have to take into account
two things:

-  alignment of motifs
-  some function to compare aligned motifs

In ``Bio.Motif`` we have 3 different functions for motif comparison,
which are based on the same idea behind motif alignment, but use
different functions to compare aligned motifs. Briefly speaking, we are
using ungapped alignment of PWMs and substitute the missing columns at
the beginning and end of the matrices with background distribution. All
three comparison functions are written in such a way, that they can be
interpreted as distance measures, however only one (``dist_dpq``)
satisfies the triangle inequality. All of them return the minimal
distance and the corresponding offset between motifs.

To show how these functions work, let us first load another motif, which
is similar to our test motif ``m``:

.. code:: verbatim

    >>> ubx=Motif.read(open("Ubx.pfm"),"jaspar-pfm")
    <Bio.Motif.Motif.Motif object at 0xc29b90>
    >>> ubx.consensus()
    Seq('TAAT', IUPACUnambiguousDNA())

The first function we’ll use to compare these motifs is based on Pearson
correlation. Since we want it to resemble a distance measure, we
actually take 1−\ *r*, where *r* is the Pearson correlation coefficient
(PCC):

.. code:: verbatim

    >>> m.dist_pearson(ubx)
    (0.41740393308237722, 2)

This means, that the best PCC between motif ``m`` and ``Ubx`` is
obtained with the following alignment:

.. code:: verbatim

    bbTAAT
    TATAAb

where ``b`` stands for background distribution. The PCC itself is
roughly 1−0.42=0.58. If we try the reverse complement of the Ubx motif:

.. code:: verbatim

    >>> m.dist_pearson(ubx.reverse_complement())
    (0.25784180151584823, 1)

We can see that the PCC is better (almost 0.75), and the alignment is
also different:

.. code:: verbatim

    bATTA
    TATAA

There are two other functions: ``dist_dpq``, which is a true metric
(satisfying traingle inequality) based on the Kullback-Leibler
divergence

.. code:: verbatim

    >>> m.dist_dpq(ubx.reverse_complement())
    (0.49292358382899853, 1)

and the ``dist_product`` method, which is based on the product of
probabilities which can be interpreted as the probability of
independently generating the same instance by both motifs.

.. code:: verbatim

    >>> m.dist_product(ubx.reverse_complement())
    (0.16224587301064275, 1)

14.9.4  *De novo* motif finding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, Biopython has only limited support for *de novo* motif
finding. Namely, we support running and parsing of AlignAce and MEME.
Since the number of motif finding tools is growing rapidly,
contributions of new parsers are welcome.

14.9.4.1  MEME
^^^^^^^^^^^^^^

Let’s assume, you have run MEME on sequences of your choice with your
favorite parameters and saved the output in the file ``meme.out``. You
can retrieve the motifs reported by MEME by running the following piece
of code:

.. code:: verbatim

    >>> motifsM = list(Motif.parse(open("meme.out"),"MEME"))
    >>> motifsM
    [<Bio.Motif.MEMEMotif.MEMEMotif object at 0xc356b0>]

Besides the most wanted list of motifs, the result object contains more
useful information, accessible through properties with self-explanatory
names:

-  ``.alphabet``
-  ``.datafile``
-  ``.sequence_names``
-  ``.version``
-  ``.command``

The motifs returned by MEMEParser can be treated exactly like regular
Motif objects (with instances), they also provide some extra
functionality, by adding additional information about the instances.

.. code:: verbatim

    >>> motifsM[0].consensus()
    Seq('CTCAATCGTA', IUPACUnambiguousDNA())

    >>> motifsM[0].instances[0].pvalue
    8.71e-07
    >>> motifsM[0].instances[0].sequence_name
    'SEQ10;'
    >>> motifsM[0].instances[0].start
    3
    >>> motifsM[0].instances[0].strand
    '+'

14.9.4.2  AlignAce
^^^^^^^^^^^^^^^^^^

We can do very similar things with AlignACE program. Assume, you have
your output in the file ``alignace.out``. You can parse your output with
the following code:

.. code:: verbatim

    >>> motifsA=list(Motif.parse(open("alignace.out"),"AlignAce"))

Again, your motifs behave as they should:

.. code:: verbatim

    >>> motifsA[0].consensus()
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

In fact you can even see, that AlignAce found a very similar motif as
MEME, it is just a longer version of a reverse complement of MEME motif:

.. code:: verbatim

    >>> motifsM[0].reverse_complement().consensus()
    Seq('TACGATTGAG', IUPACUnambiguousDNA())

If you have AlignAce installed on the same machine, you can also run it
directly from Biopython. Short example of how this can be done is shown
below (other parameters can be specified as keyword parameters):

.. code:: verbatim

    >>> command="/opt/bin/AlignACE"
    >>> input_file="test.fa"
    >>> from Bio.Motif.Applications import AlignAceCommandline
    >>> cmd = AlignAceCommandline(cmd=command,input=input_file,gcback=0.6,numcols=10)
    >>> stdout,stderr= cmd()

Since AlignAce prints all its output to standard output, you can get to
your motifs by parsing the first part of the result:

.. code:: verbatim

    motifs=list(Motif.parse(stdout,"AlignAce"))


