Chapter 2  Quick Start – What can you do with Biopython?
========================================================

This section is designed to get you started quickly with Biopython, and
to give a general overview of what is available and how to use it. All
of the examples in this section assume that you have some general
working knowledge of Python, and that you have successfully installed
Biopython on your system. If you think you need to brush up on your
Python, the main Python web site provides quite a bit of free
documentation to get started with
(```http://www.python.org/doc/`` <http://www.python.org/doc/>`__).

Since much biological work on the computer involves connecting with
databases on the internet, some of the examples will also require a
working internet connection in order to run.

Now that that is all out of the way, let’s get into what we can do with
Biopython.

2.1  General overview of what Biopython provides
------------------------------------------------

As mentioned in the introduction, Biopython is a set of libraries to
provide the ability to deal with “things” of interest to biologists
working on the computer. In general this means that you will need to
have at least some programming experience (in Python, of course!) or at
least an interest in learning to program. Biopython’s job is to make
your job easier as a programmer by supplying reusable libraries so that
you can focus on answering your specific question of interest, instead
of focusing on the internals of parsing a particular file format (of
course, if you want to help by writing a parser that doesn’t exist and
contributing it to Biopython, please go ahead!). So Biopython’s job is
to make you happy!

One thing to note about Biopython is that it often provides multiple
ways of “doing the same thing.” Things have improved in recent releases,
but this can still be frustrating as in Python there should ideally be
one right way to do something. However, this can also be a real benefit
because it gives you lots of flexibility and control over the libraries.
The tutorial helps to show you the common or easy ways to do things so
that you can just make things work. To learn more about the alternative
possibilities, look in the Cookbook (Chapter `18 <#chapter:cookbook>`__,
this has some cools tricks and tips), the Advanced section
(Chapter `20 <#chapter:advanced>`__), the built in “docstrings” (via the
Python help command, or the `API
documentation <http://biopython.org/DIST/docs/api/>`__) or ultimately
the code itself.

2.2  Working with sequences
---------------------------

Disputably (of course!), the central object in bioinformatics is the
sequence. Thus, we’ll start with a quick introduction to the Biopython
mechanisms for dealing with sequences, the ``Seq`` object, which we’ll
discuss in more detail in Chapter \ `3 <#chapter:Bio.Seq>`__.

Most of the time when we think about sequences we have in my mind a
string of letters like ‘\ ``AGTACACTGGT``\ ’. You can create such
``Seq`` object with this sequence as follows - the “\ ``>>>``\ ”
represents the Python prompt followed by what you would type in:

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> my_seq = Seq("AGTACACTGGT")
    >>> my_seq
    Seq('AGTACACTGGT', Alphabet())
    >>> print my_seq
    AGTACACTGGT
    >>> my_seq.alphabet
    Alphabet()

What we have here is a sequence object with a *generic* alphabet -
reflecting the fact we have *not* specified if this is a DNA or protein
sequence (okay, a protein with a lot of Alanines, Glycines, Cysteines
and Threonines!). We’ll talk more about alphabets in
Chapter \ `3 <#chapter:Bio.Seq>`__.

In addition to having an alphabet, the ``Seq`` object differs from the
Python string in the methods it supports. You can’t do this with a plain
string:

.. code:: verbatim

    >>> my_seq
    Seq('AGTACACTGGT', Alphabet())
    >>> my_seq.complement()
    Seq('TCATGTGACCA', Alphabet())
    >>> my_seq.reverse_complement()
    Seq('ACCAGTGTACT', Alphabet())

The next most important class is the ``SeqRecord`` or Sequence Record.
This holds a sequence (as a ``Seq`` object) with additional annotation
including an identifier, name and description. The ``Bio.SeqIO`` module
for reading and writing sequence file formats works with ``SeqRecord``
objects, which will be introduced below and covered in more detail by
Chapter \ `5 <#chapter:Bio.SeqIO>`__.

This covers the basic features and uses of the Biopython sequence class.
Now that you’ve got some idea of what it is like to interact with the
Biopython libraries, it’s time to delve into the fun, fun world of
dealing with biological file formats!

2.3  A usage example
--------------------

Before we jump right into parsers and everything else to do with
Biopython, let’s set up an example to motivate everything we do and make
life more interesting. After all, if there wasn’t any biology in this
tutorial, why would you want you read it?

Since I love plants, I think we’re just going to have to have a plant
based example (sorry to all the fans of other organisms out there!).
Having just completed a recent trip to our local greenhouse, we’ve
suddenly developed an incredible obsession with Lady Slipper Orchids (if
you wonder why, have a look at some `Lady Slipper Orchids photos on
Flickr <http://www.flickr.com/search/?q=lady+slipper+orchid&s=int&z=t>`__,
or try a `Google Image
Search <http://images.google.com/images?q=lady%20slipper%20orchid>`__).

Of course, orchids are not only beautiful to look at, they are also
extremely interesting for people studying evolution and systematics. So
let’s suppose we’re thinking about writing a funding proposal to do a
molecular study of Lady Slipper evolution, and would like to see what
kind of research has already been done and how we can add to that.

After a little bit of reading up we discover that the Lady Slipper
Orchids are in the Orchidaceae family and the Cypripedioideae sub-family
and are made up of 5 genera: *Cypripedium*, *Paphiopedilum*,
*Phragmipedium*, *Selenipedium* and *Mexipedium*.

That gives us enough to get started delving for more information. So,
let’s look at how the Biopython tools can help us. We’ll start with
sequence parsing in Section \ `2.4 <#sec:sequence-parsing>`__, but the
orchids will be back later on as well - for example we’ll search PubMed
for papers about orchids and extract sequence data from GenBank in
Chapter \ `9 <#chapter:entrez>`__, extract data from Swiss-Prot from
certain orchid proteins in Chapter \ `10 <#chapter:swiss_prot>`__, and
work with ClustalW multiple sequence alignments of orchid proteins in
Section \ `6.4.1 <#sec:align_clustal>`__.

2.4  Parsing sequence file formats
----------------------------------

A large part of much bioinformatics work involves dealing with the many
types of file formats designed to hold biological data. These files are
loaded with interesting biological data, and a special challenge is
parsing these files into a format so that you can manipulate them with
some kind of programming language. However the task of parsing these
files can be frustrated by the fact that the formats can change quite
regularly, and that formats may contain small subtleties which can break
even the most well designed parsers.

We are now going to briefly introduce the ``Bio.SeqIO`` module – you can
find out more in Chapter \ `5 <#chapter:Bio.SeqIO>`__. We’ll start with
an online search for our friends, the lady slipper orchids. To keep this
introduction simple, we’re just using the NCBI website by hand. Let’s
just take a look through the nucleotide databases at NCBI, using an
Entrez online search
(```http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?db=Nucleotide`` <http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?db=Nucleotide>`__)
for everything mentioning the text Cypripedioideae (this is the
subfamily of lady slipper orchids).

When this tutorial was originally written, this search gave us only 94
hits, which we saved as a FASTA formatted text file and as a GenBank
formatted text file (files
```ls_orchid.fasta`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.fasta>`__
and
```ls_orchid.gbk`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.gbk>`__,
also included with the Biopython source code under
``docs/tutorial/examples/``).

If you run the search today, you’ll get hundreds of results! When
following the tutorial, if you want to see the same list of genes, just
download the two files above or copy them from ``docs/examples/`` in the
Biopython source code. In
Section \ `2.5 <#sec:connecting-with-biological-databases>`__ we will
look at how to do a search like this from within Python.

2.4.1  Simple FASTA parsing example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you open the lady slipper orchids FASTA file
```ls_orchid.fasta`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.fasta>`__
in your favourite text editor, you’ll see that the file starts like
this:

.. code:: verbatim

    >gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTG
    AATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGTGACCCTGATTTGTTGTTGGG
    ...

It contains 94 records, each has a line starting with “\ ``>``\ ”
(greater-than symbol) followed by the sequence on one or more lines. Now
try this in Python:

.. code:: verbatim

    from Bio import SeqIO
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print seq_record.id
        print repr(seq_record.seq)
        print len(seq_record)

You should get something like this on your screen:

.. code:: verbatim

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC', SingleLetterAlphabet())
    740
    ...
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC', SingleLetterAlphabet())
    592

Notice that the FASTA format does not specify the alphabet, so
``Bio.SeqIO`` has defaulted to the rather generic
``SingleLetterAlphabet()`` rather than something DNA specific.

2.4.2  Simple GenBank parsing example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now let’s load the GenBank file
```ls_orchid.gbk`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.gbk>`__
instead - notice that the code to do this is almost identical to the
snippet used above for the FASTA file - the only difference is we change
the filename and the format string:

.. code:: verbatim

    from Bio import SeqIO
    for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
        print seq_record.id
        print repr(seq_record.seq)
        print len(seq_record)

This should give:

.. code:: verbatim

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC', IUPACAmbiguousDNA())
    740
    ...
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC', IUPACAmbiguousDNA())
    592

This time ``Bio.SeqIO`` has been able to choose a sensible alphabet,
IUPAC Ambiguous DNA. You’ll also notice that a shorter string has been
used as the ``seq_record.id`` in this case.

2.4.3  I love parsing – please don’t stop talking about it!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Biopython has a lot of parsers, and each has its own little special
niches based on the sequence format it is parsing and all of that.
Chapter \ `5 <#chapter:Bio.SeqIO>`__ covers ``Bio.SeqIO`` in more
detail, while Chapter \ `6 <#chapter:Bio.AlignIO>`__ introduces
``Bio.AlignIO`` for sequence alignments.

While the most popular file formats have parsers integrated into
``Bio.SeqIO`` and/or ``Bio.AlignIO``, for some of the rarer and unloved
file formats there is either no parser at all, or an old parser which
has not been linked in yet. Please also check the wiki pages
```http://biopython.org/wiki/SeqIO`` <http://biopython.org/wiki/SeqIO>`__
and
```http://biopython.org/wiki/AlignIO`` <http://biopython.org/wiki/AlignIO>`__
for the latest information, or ask on the mailing list. The wiki pages
should include an up to date list of supported file types, and some
additional examples.

The next place to look for information about specific parsers and how to
do cool things with them is in the Cookbook
(Chapter `18 <#chapter:cookbook>`__ of this Tutorial). If you don’t find
the information you are looking for, please consider helping out your
poor overworked documentors and submitting a cookbook entry about it!
(once you figure out how to do it, that is!)

2.5  Connecting with biological databases
-----------------------------------------

One of the very common things that you need to do in bioinformatics is
extract information from biological databases. It can be quite tedious
to access these databases manually, especially if you have a lot of
repetitive work to do. Biopython attempts to save you time and energy by
making some on-line databases available from Python scripts. Currently,
Biopython has code to extract information from the following databases:

-  `Entrez <http://www.ncbi.nlm.nih.gov/Entrez/>`__ (and
   `PubMed <http://www.ncbi.nlm.nih.gov/PubMed/>`__) from the NCBI – See
   Chapter \ `9 <#chapter:entrez>`__.
-  `ExPASy <http://www.expasy.org/>`__ – See
   Chapter \ `10 <#chapter:swiss_prot>`__.
-  `SCOP <http://scop.mrc-lmb.cam.ac.uk/scop/>`__ – See the
   ``Bio.SCOP.search()`` function.

The code in these modules basically makes it easy to write Python code
that interact with the CGI scripts on these pages, so that you can get
results in an easy to deal with format. In some cases, the results can
be tightly integrated with the Biopython parsers to make it even easier
to extract information.

2.6  What to do next
--------------------

Now that you’ve made it this far, you hopefully have a good
understanding of the basics of Biopython and are ready to start using it
for doing useful work. The best thing to do now is finish reading this
tutorial, and then if you want start snooping around in the source code,
and looking at the automatically generated documentation.

Once you get a picture of what you want to do, and what libraries in
Biopython will do it, you should take a peak at the Cookbook
(Chapter `18 <#chapter:cookbook>`__), which may have example code to do
something similar to what you want to do.

If you know what you want to do, but can’t figure out how to do it,
please feel free to post questions to the main Biopython list (see
```http://biopython.org/wiki/Mailing_lists`` <http://biopython.org/wiki/Mailing_lists>`__).
This will not only help us answer your question, it will also allow us
to improve the documentation so it can help the next person do what you
want to do.

Enjoy the code!
