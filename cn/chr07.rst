第7章  BLAST
================

嘿，每个人都爱BLAST对吧?我的意思是，还有比用BLAST来比较你的一条序列和世
界上其他序列更简单的吗？但是，当然这部分不是来讲BLAST有多酷，因为这是众
所周知的。而是关于一些BLAST的问题 - 处理BLAST产生的大量的数据以及自动
BLAsT非常困难。

幸运的是，Biopython的开发者们对这非常了解，所以他们已经开发了很多工具来
处理BLAST，让事情变得更简单。这部分详细介绍怎么使用这些工具，并且用它们
做一些有用的工作。

处理BLAST可以分成两步，都可以用Biopython来完成。首先，用你的查询序列来作
BLAST获得输出结果。其次，Python来解析结果用于进一步分析。

你第一次作BLAST可能是通过NCBI的网页服务。实际上，有很多的途径可以作NCBI，
大致可以分成几类。最重要的区别是本地BLAST（在自己的机器上），和远程BLAST
（在其他机器上，一般是在NCBI服务器上）。我们将以一个Python脚本调用NCBI的
在线BLAST服务开始本章节。

*NOTE*: 接下来第 \ `8 <#chapter:searchio>`__ 章介绍``Bio.SearchIO``，一个
Biopython的*实验*模块。我们打算用这个模块最终代替老的``Bio.Blast``模块，因
为它提供一个可以处理其他相关序列搜索工具的通用框架。但是，在发布稳定版前，
产品代码中还是请用``Bio.Blast``模块处理NCBI的BLAST。

7.1  在线BLAST
------------------------------------

我们用 ``Bio.Blast.NCBIWWW`` 模块的 ``qblast()`` 方法调用在线版BLAST。有三
个固定的参数：

-  第一个参数是用于搜索的blast程序，小写字母。程序的参数和描述可以找到在
   ```http://www.ncbi.nlm.nih.gov/BLAST/blast_program.shtml`` <http://www.ncbi.nlm.nih.gov/BLAST/blast_program.shtml>`__.
   当前的 ``qblast`` 只支持blastn, blastp, blastx, tblast and tblastx.
-  第二个参数指定搜索的数据库，同样，参数选项可以再NCBI的网站找到
   ```http://www.ncbi.nlm.nih.gov/BLAST/blast_databases.shtml`` <http://www.ncbi.nlm.nih.gov/BLAST/blast_databases.shtml>`__.
-  第三个参数是一个包含你查询序列的字符串，可以是fasta格式的序列或是像GI一
样的标识符。

``qblast`` 方法同样接受很多其他的参数选项，和你在在线BLAST网页设置的参数
相似。我们将展示少数参数：

-  ``qblast`` 方法能以多种格式返回BLAST结果，你能选择可选格式的关键字：
   ``"HTML"``, ``"Text"``, ``"ASN.1"``, or ``"XML"``.默认是XML格式，能够被
   后面\ `7.3 <#sec:parsing-blast>`__部分提到的解析器解析。

-  ``expect``参数设置e-value的阈值。

想了解跟多的BLAST可选参数，我们推荐你参考NCBI的文档，或者嵌入到Biopython中
的文档：

.. code:: verbatim

    >>> from Bio.Blast import NCBIWWW
    >>> help(NCBIWWW.qblast)
    ...

注意，NCBI在线BLAST的默认设置和QBLAST的不太一样。如果获得不同的结果，你需要
检查参数（如e-value的阈值和空位罚分值）。

例如，如果你有一条核酸序列，想用BLASTN搜索核酸数据库，并且你知道查询序列的GI
号，你可以用：

.. code:: verbatim

    >>> from Bio.Blast import NCBIWWW
    >>> result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")

或者，如果我们有FASTA格式的查询序列，我们只需要打开并读取文件得到一个字符串
并用它作为一个查询参数：

.. code:: verbatim

    >>> from Bio.Blast import NCBIWWW
    >>> fasta_string = open("m_cold.fasta").read()
    >>> result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

我们同样可以读取FASTA文件得到一个``SeqRecord``，只用序列本身用于参数：

.. code:: verbatim

    >>> from Bio.Blast import NCBIWWW
    >>> from Bio import SeqIO
    >>> record = SeqIO.read("m_cold.fasta", format="fasta")
    >>> result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

Supplying just the sequence means that BLAST will assign an identifier
for your sequence automatically. You might prefer to use the
``SeqRecord`` object’s format method to make a fasta string (which will
include the existing identifier):

.. code:: verbatim

    >>> from Bio.Blast import NCBIWWW
    >>> from Bio import SeqIO
    >>> record = SeqIO.read("m_cold.fasta", format="fasta")
    >>> result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

This approach makes more sense if you have your sequence(s) in a
non-FASTA file format which you can extract using ``Bio.SeqIO`` (see
Chapter \ `5 <#chapter:Bio.SeqIO>`__).

Whatever arguments you give the ``qblast()`` function, you should get
back your results in a handle object (by default in XML format). The
next step would be to parse the XML output into Python objects
representing the search results (Section `7.3 <#sec:parsing-blast>`__),
but you might want to save a local copy of the output file first. I find
this especially useful when debugging my code that extracts info from
the BLAST results (because re-running the online search is slow and
wastes the NCBI computer time).

We need to be a bit careful since we can use ``result_handle.read()`` to
read the BLAST output only once – calling ``result_handle.read()`` again
returns an empty string.

.. code:: verbatim

    >>> save_file = open("my_blast.xml", "w")
    >>> save_file.write(result_handle.read())
    >>> save_file.close()
    >>> result_handle.close()

After doing this, the results are in the file ``my_blast.xml`` and the
original handle has had all its data extracted (so we closed it).
However, the ``parse`` function of the BLAST parser (described
in \ `7.3 <#sec:parsing-blast>`__) takes a file-handle-like object, so
we can just open the saved file for input:

.. code:: verbatim

    >>> result_handle = open("my_blast.xml")

Now that we’ve got the BLAST results back into a handle again, we are
ready to do something with them, so this leads us right into the parsing
section (see Section \ `7.3 <#sec:parsing-blast>`__ below). You may want
to jump ahead to that now ….

7.2  Running BLAST locally
--------------------------

7.2.1  Introduction
~~~~~~~~~~~~~~~~~~~

Running BLAST locally (as opposed to over the internet, see
Section \ `7.1 <#sec:running-www-blast>`__) has at least major two
advantages:

-  Local BLAST may be faster than BLAST over the internet;
-  Local BLAST allows you to make your own database to search for
   sequences against.

Dealing with proprietary or unpublished sequence data can be another
reason to run BLAST locally. You may not be allowed to redistribute the
sequences, so submitting them to the NCBI as a BLAST query would not be
an option.

Unfortunately, there are some major drawbacks too – installing all the
bits and getting it setup right takes some effort:

-  Local BLAST requires command line tools to be installed.
-  Local BLAST requires (large) BLAST databases to be setup (and
   potentially kept up to date).

To further confuse matters there are at least four different standalone
BLAST packages, and there are also other tools which can produce
imitation BLAST output files, such as BLAT.

7.2.2  Standalone NCBI “legacy” BLAST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`NCBI “legacy”
BLAST <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__
included command line tools ``blastall``, ``blastpgp`` and ``rpsblast``.
This was the most widely used standalone BLAST tool up until its
replacement BLAST+ was released by the NCBI.

The ``Bio.Blast.Applications`` module has wrappers for the “legacy” NCBI
BLAST tools like ``blastall``, ``blastpgp`` and ``rpsblast``, and there
are also helper functions in ``Bio.Blast.NCBIStandalone``. These are now
considered obsolete, and will be deprecated and eventually removed from
Biopython as people move over to the replacement BLAST+ suite.

To try and avoid confusion, we will not cover calling these old tools
from Biopython in this tutorial. Have a look at the older edition of
this tutorial included with Biopython 1.52 if you are curious (look at
the Tutorial PDF or HTML file in the Doc directory within
``biopython-1.52.tar.gz`` or ``biopython-1.52.zip``).

7.2.3  Standalone NCBI BLAST+
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`NCBI “new”
BLAST+ <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__
was released in 2009. This replaces the old NCBI “legacy” BLAST package.
The ``Bio.Blast.Applications`` module has wrappers for these “new” tools
like ``blastn``, ``blastp``, ``blastx``, ``tblastn``, ``tblastx`` (which
all used to be handled by ``blastall``), ``psiblast`` (replacing
``blastpgp``) and ``rpsblast`` and ``rpstblastn`` (which replace the old
``rpsblast``). We don’t include a wrapper for the ``makeblastdb`` used
in BLAST+ to build a local BLAST database from FASTA file, nor the
equivalent tool ``formatdb`` in “legacy” BLAST.

This section will show briefly how to use these tools from within
Python. If you have already read or tried the alignment tool examples in
Section \ `6.4 <#sec:alignment-tools>`__ this should all seem quite
straightforward. First, we construct a command line string (as you would
type in at the command line prompt if running standalone BLAST by hand).
Then we can execute this command from within Python.

For example, taking a FASTA file of gene nucleotide sequences, you might
want to run a BLASTX (translation) search against the non-redundant (NR)
protein database. Assuming you (or your systems administrator) has
downloaded and installed the NR database, you might run:

.. code:: verbatim

    blastx -query opuntia.fasta -db nr -out opuntia.xml -evalue 0.001 -outfmt 5

This should run BLASTX against the NR database, using an expectation
cut-off value of 0.001 and produce XML output to the specified file
(which we can then parse). On my computer this takes about six minutes -
a good reason to save the output to a file so you and repeat any
analysis as needed.

From within Biopython we can use the NCBI BLASTX wrapper from the
``Bio.Blast.Applications`` module to build the command line string, and
run it:

.. code:: verbatim

    >>> from Bio.Blast.Applications import NcbiblastxCommandline
    >>> help(NcbiblastxCommandline)
    ...
    >>> blastx_cline = NcbiblastxCommandline(query="opuntia.fasta", db="nr", evalue=0.001,
    ...                                      outfmt=5, out="opuntia.xml")
    >>> blastx_cline
    NcbiblastxCommandline(cmd='blastx', out='opuntia.xml', outfmt=5, query='opuntia.fasta',
    db='nr', evalue=0.001)
    >>> print blastx_cline
    blastx -out opuntia.xml -outfmt 5 -query opuntia.fasta -db nr -evalue 0.001
    >>> stdout, stderr = blastx_cline()

In this example there shouldn’t be any output from BLASTX to the
terminal, so stdout and stderr should be empty. You may want to check
the output file ``opuntia.xml`` has been created.

As you may recall from earlier examples in the tutorial, the
``opuntia.fasta`` contains seven sequences, so the BLAST XML output
should contain multiple results. Therefore use
``Bio.Blast.NCBIXML.parse()`` to parse it as described below in
Section \ `7.3 <#sec:parsing-blast>`__.

7.2.4  WU-BLAST and AB-BLAST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may also come across `Washington University
BLAST <http://blast.wustl.edu/>`__ (WU-BLAST), and its successor,
`Advanced Biocomputing BLAST <http://blast.advbiocomp.com>`__ (AB-BLAST,
released in 2009, not free/open source). These packages include the
command line tools ``wu-blastall`` and ``ab-blastall``.

Biopython does not currently provide wrappers for calling these tools,
but should be able to parse any NCBI compatible output from them.

7.3  Parsing BLAST output
-------------------------

As mentioned above, BLAST can generate output in various formats, such
as XML, HTML, and plain text. Originally, Biopython had parsers for
BLAST plain text and HTML output, as these were the only output formats
offered at the time. Unfortunately, the BLAST output in these formats
kept changing, each time breaking the Biopython parsers. Our HTML BLAST
parser has been removed, but the plain text BLAST parser is still
available (see Section \ `7.5 <#sec:parsing-blast-deprecated>`__). Use
it at your own risk, it may or may not work, depending on which BLAST
version you’re using.

As keeping up with changes in BLAST became a hopeless endeavor,
especially with users running different BLAST versions, we now recommend
to parse the output in XML format, which can be generated by recent
versions of BLAST. Not only is the XML output more stable than the plain
text and HTML output, it is also much easier to parse automatically,
making Biopython a whole lot more stable.

You can get BLAST output in XML format in various ways. For the parser,
it doesn’t matter how the output was generated, as long as it is in the
XML format.

-  You can use Biopython to run BLAST over the internet, as described in
   section \ `7.1 <#sec:running-www-blast>`__.
-  You can use Biopython to run BLAST locally, as described in
   section \ `7.2 <#sec:running-local-blast>`__.
-  You can do the BLAST search yourself on the NCBI site through your
   web browser, and then save the results. You need to choose XML as the
   format in which to receive the results, and save the final BLAST page
   you get (you know, the one with all of the interesting results!) to a
   file.
-  You can also run BLAST locally without using Biopython, and save the
   output in a file. Again, you need to choose XML as the format in
   which to receive the results.

The important point is that you do not have to use Biopython scripts to
fetch the data in order to be able to parse it. Doing things in one of
these ways, you then need to get a handle to the results. In Python, a
handle is just a nice general way of describing input to any info source
so that the info can be retrieved using ``read()`` and ``readline()``
functions (see Section sec:appendix-handles).

If you followed the code above for interacting with BLAST through a
script, then you already have ``result_handle``, the handle to the BLAST
results. For example, using a GI number to do an online search:

.. code:: verbatim

    >>> from Bio.Blast import NCBIWWW
    >>> result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")

If instead you ran BLAST some other way, and have the BLAST output (in
XML format) in the file ``my_blast.xml``, all you need to do is to open
the file for reading:

.. code:: verbatim

    >>> result_handle = open("my_blast.xml")

Now that we’ve got a handle, we are ready to parse the output. The code
to parse it is really quite small. If you expect a single BLAST result
(i.e. you used a single query):

.. code:: verbatim

    >>> from Bio.Blast import NCBIXML
    >>> blast_record = NCBIXML.read(result_handle)

or, if you have lots of results (i.e. multiple query sequences):

.. code:: verbatim

    >>> from Bio.Blast import NCBIXML
    >>> blast_records = NCBIXML.parse(result_handle)

Just like ``Bio.SeqIO`` and ``Bio.AlignIO`` (see
Chapters \ `5 <#chapter:Bio.SeqIO>`__
and \ `6 <#chapter:Bio.AlignIO>`__), we have a pair of input functions,
``read`` and ``parse``, where ``read`` is for when you have exactly one
object, and ``parse`` is an iterator for when you can have lots of
objects – but instead of getting ``SeqRecord`` or
``MultipleSeqAlignment`` objects, we get BLAST record objects.

To be able to handle the situation where the BLAST file may be huge,
containing thousands of results, ``NCBIXML.parse()`` returns an
iterator. In plain English, an iterator allows you to step through the
BLAST output, retrieving BLAST records one by one for each BLAST search
result:

.. code:: verbatim

    >>> from Bio.Blast import NCBIXML
    >>> blast_records = NCBIXML.parse(result_handle)
    >>> blast_record = blast_records.next()
    # ... do something with blast_record
    >>> blast_record = blast_records.next()
    # ... do something with blast_record
    >>> blast_record = blast_records.next()
    # ... do something with blast_record
    >>> blast_record = blast_records.next()
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    StopIteration
    # No further records

Or, you can use a ``for``-loop:

.. code:: verbatim

    >>> for blast_record in blast_records:
    ...     # Do something with blast_record

Note though that you can step through the BLAST records only once.
Usually, from each BLAST record you would save the information that you
are interested in. If you want to save all returned BLAST records, you
can convert the iterator into a list:

.. code:: verbatim

    >>> blast_records = list(blast_records)

Now you can access each BLAST record in the list with an index as usual.
If your BLAST file is huge though, you may run into memory problems
trying to save them all in a list.

Usually, you’ll be running one BLAST search at a time. Then, all you
need to do is to pick up the first (and only) BLAST record in
``blast_records``:

.. code:: verbatim

    >>> from Bio.Blast import NCBIXML
    >>> blast_records = NCBIXML.parse(result_handle)
    >>> blast_record = blast_records.next()

or more elegantly:

.. code:: verbatim

    >>> from Bio.Blast import NCBIXML
    >>> blast_record = NCBIXML.read(result_handle)

I guess by now you’re wondering what is in a BLAST record.

7.4  The BLAST record class
---------------------------

A BLAST Record contains everything you might ever want to extract from
the BLAST output. Right now we’ll just show an example of how to get
some info out of the BLAST report, but if you want something in
particular that is not described here, look at the info on the record
class in detail, and take a gander into the code or automatically
generated documentation – the docstrings have lots of good info about
what is stored in each piece of information.

To continue with our example, let’s just print out some summary info
about all hits in our blast report greater than a particular threshold.
The following code does this:

.. code:: verbatim

    >>> E_VALUE_THRESH = 0.04

    >>> for alignment in blast_record.alignments:
    ...     for hsp in alignment.hsps:
    ...         if hsp.expect < E_VALUE_THRESH:
    ...             print '****Alignment****'
    ...             print 'sequence:', alignment.title
    ...             print 'length:', alignment.length
    ...             print 'e value:', hsp.expect
    ...             print hsp.query[0:75] + '...'
    ...             print hsp.match[0:75] + '...'
    ...             print hsp.sbjct[0:75] + '...'

This will print out summary reports like the following:

.. code:: verbatim

    ****Alignment****
    sequence: >gb|AF283004.1|AF283004 Arabidopsis thaliana cold acclimation protein WCOR413-like protein
    alpha form mRNA, complete cds
    length: 783
    e value: 0.034
    tacttgttgatattggatcgaacaaactggagaaccaacatgctcacgtcacttttagtcccttacatattcctc...
    ||||||||| | ||||||||||| || ||||  || || |||||||| |||||| |  | |||||||| ||| ||...
    tacttgttggtgttggatcgaaccaattggaagacgaatatgctcacatcacttctcattccttacatcttcttc...

Basically, you can do anything you want to with the info in the BLAST
report once you have parsed it. This will, of course, depend on what you
want to use it for, but hopefully this helps you get started on doing
what you need to do!

An important consideration for extracting information from a BLAST
report is the type of objects that the information is stored in. In
Biopython, the parsers return ``Record`` objects, either ``Blast`` or
``PSIBlast`` depending on what you are parsing. These objects are
defined in ``Bio.Blast.Record`` and are quite complete.

Here are my attempts at UML class diagrams for the ``Blast`` and
``PSIBlast`` record classes. If you are good at UML and see
mistakes/improvements that can be made, please let me know. The Blast
class diagram is shown in Figure \ `7.4 <#fig:blastrecord>`__.

|image1|

The PSIBlast record object is similar, but has support for the rounds
that are used in the iteration steps of PSIBlast. The class diagram for
PSIBlast is shown in Figure \ `7.4 <#fig:psiblastrecord>`__.

|image2|

7.5  Deprecated BLAST parsers
-----------------------------

Older versions of Biopython had parsers for BLAST output in plain text
or HTML format. Over the years, we discovered that it is very hard to
maintain these parsers in working order. Basically, any small change to
the BLAST output in newly released BLAST versions tends to cause the
plain text and HTML parsers to break. We therefore recommend parsing
BLAST output in XML format, as described in
section \ `7.3 <#sec:parsing-blast>`__.

Depending on which BLAST versions or programs you’re using, our plain
text BLAST parser may or may not work. Use it at your own risk!

7.5.1  Parsing plain-text BLAST output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The plain text BLAST parser is located in ``Bio.Blast.NCBIStandalone``.

As with the XML parser, we need to have a handle object that we can pass
to the parser. The handle must implement the ``readline()`` method and
do this properly. The common ways to get such a handle are to either use
the provided ``blastall`` or ``blastpgp`` functions to run the local
blast, or to run a local blast via the command line, and then do
something like the following:

.. code:: verbatim

    >>> result_handle = open("my_file_of_blast_output.txt")

Well, now that we’ve got a handle (which we’ll call ``result_handle``),
we are ready to parse it. This can be done with the following code:

.. code:: verbatim

    >>> from Bio.Blast import NCBIStandalone
    >>> blast_parser = NCBIStandalone.BlastParser()
    >>> blast_record = blast_parser.parse(result_handle)

This will parse the BLAST report into a Blast Record class (either a
Blast or a PSIBlast record, depending on what you are parsing) so that
you can extract the information from it. In our case, let’s just use
print out a quick summary of all of the alignments greater than some
threshold value.

.. code:: verbatim

    >>> E_VALUE_THRESH = 0.04
    >>> for alignment in blast_record.alignments:
    ...     for hsp in alignment.hsps:
    ...         if hsp.expect < E_VALUE_THRESH:
    ...             print '****Alignment****'
    ...             print 'sequence:', alignment.title
    ...             print 'length:', alignment.length
    ...             print 'e value:', hsp.expect
    ...             print hsp.query[0:75] + '...'
    ...             print hsp.match[0:75] + '...'
    ...             print hsp.sbjct[0:75] + '...'

If you also read the section \ `7.3 <#sec:parsing-blast>`__ on parsing
BLAST XML output, you’ll notice that the above code is identical to what
is found in that section. Once you parse something into a record class
you can deal with it independent of the format of the original BLAST
info you were parsing. Pretty snazzy!

Sure, parsing one record is great, but I’ve got a BLAST file with tons
of records – how can I parse them all? Well, fear not, the answer lies
in the very next section.

7.5.2  Parsing a plain-text BLAST file full of BLAST runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can do this using the blast iterator. To set up an iterator, we first
set up a parser, to parse our blast reports in Blast Record objects:

.. code:: verbatim

    >>> from Bio.Blast import NCBIStandalone
    >>> blast_parser = NCBIStandalone.BlastParser()

Then we will assume we have a handle to a bunch of blast records, which
we’ll call ``result_handle``. Getting a handle is described in full
detail above in the blast parsing sections.

Now that we’ve got a parser and a handle, we are ready to set up the
iterator with the following command:

.. code:: verbatim

    >>> blast_iterator = NCBIStandalone.Iterator(result_handle, blast_parser)

The second option, the parser, is optional. If we don’t supply a parser,
then the iterator will just return the raw BLAST reports one at a time.

Now that we’ve got an iterator, we start retrieving blast records
(generated by our parser) using ``next()``:

.. code:: verbatim

    >>> blast_record = blast_iterator.next()

Each call to next will return a new record that we can deal with. Now we
can iterate through this records and generate our old favorite, a nice
little blast report:

.. code:: verbatim

    >>> for blast_record in blast_iterator:
    ...     E_VALUE_THRESH = 0.04
    ...     for alignment in blast_record.alignments:
    ...         for hsp in alignment.hsps:
    ...             if hsp.expect < E_VALUE_THRESH:
    ...                 print '****Alignment****'
    ...                 print 'sequence:', alignment.title
    ...                 print 'length:', alignment.length
    ...                 print 'e value:', hsp.expect
    ...                 if len(hsp.query) > 75:
    ...                     dots = '...'
    ...                 else:
    ...                     dots = ''
    ...                 print hsp.query[0:75] + dots
    ...                 print hsp.match[0:75] + dots
    ...                 print hsp.sbjct[0:75] + dots

The iterator allows you to deal with huge blast records without any
memory problems, since things are read in one at a time. I have parsed
tremendously huge files without any problems using this.

7.5.3  Finding a bad record somewhere in a huge plain-text BLAST file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One really ugly problem that happens to me is that I’ll be parsing a
huge blast file for a while, and the parser will bomb out with a
ValueError. This is a serious problem, since you can’t tell if the
ValueError is due to a parser problem, or a problem with the BLAST. To
make it even worse, you have no idea where the parse failed, so you
can’t just ignore the error, since this could be ignoring an important
data point.

We used to have to make a little script to get around this problem, but
the ``Bio.Blast`` module now includes a ``BlastErrorParser`` which
really helps make this easier. The ``BlastErrorParser`` works very
similar to the regular ``BlastParser``, but it adds an extra layer of
work by catching ValueErrors that are generated by the parser, and
attempting to diagnose the errors.

Let’s take a look at using this parser – first we define the file we are
going to parse and the file to write the problem reports to:

.. code:: verbatim

    >>> import os
    >>> blast_file = os.path.join(os.getcwd(), "blast_out", "big_blast.out")
    >>> error_file = os.path.join(os.getcwd(), "blast_out", "big_blast.problems")

Now we want to get a ``BlastErrorParser``:

.. code:: verbatim

    >>> from Bio.Blast import NCBIStandalone
    >>> error_handle = open(error_file, "w")
    >>> blast_error_parser = NCBIStandalone.BlastErrorParser(error_handle)

Notice that the parser take an optional argument of a handle. If a
handle is passed, then the parser will write any blast records which
generate a ValueError to this handle. Otherwise, these records will not
be recorded.

Now we can use the ``BlastErrorParser`` just like a regular blast
parser. Specifically, we might want to make an iterator that goes
through our blast records one at a time and parses them with the error
parser:

.. code:: verbatim

    >>> result_handle = open(blast_file)
    >>> iterator = NCBIStandalone.Iterator(result_handle, blast_error_parser)

We can read these records one a time, but now we can catch and deal with
errors that are due to problems with Blast (and not with the parser
itself):

.. code:: verbatim

    >>> try:
    ...     next_record = iterator.next()
    ... except NCBIStandalone.LowQualityBlastError, info:
    ...     print "LowQualityBlastError detected in id %s" % info[1]

The ``.next()`` method is normally called indirectly via a ``for``-loop.
Right now the ``BlastErrorParser`` can generate the following errors:

-  ``ValueError`` – This is the same error generated by the regular
   BlastParser, and is due to the parser not being able to parse a
   specific file. This is normally either due to a bug in the parser, or
   some kind of discrepancy between the version of BLAST you are using
   and the versions the parser is able to handle.
-  ``LowQualityBlastError`` – When BLASTing a sequence that is of really
   bad quality (for example, a short sequence that is basically a
   stretch of one nucleotide), it seems that Blast ends up masking out
   the entire sequence and ending up with nothing to parse. In this case
   it will produce a truncated report that causes the parser to generate
   a ValueError. ``LowQualityBlastError`` is reported in these cases.
   This error returns an info item with the following information:

   -  ``item[0]`` – The error message
   -  ``item[1]`` – The id of the input record that caused the error.
      This is really useful if you want to record all of the records
      that are causing problems.

As mentioned, with each error generated, the BlastErrorParser will write
the offending record to the specified ``error_handle``. You can then go
ahead and look and these and deal with them as you see fit. Either you
will be able to debug the parser with a single blast report, or will
find out problems in your blast runs. Either way, it will definitely be
a useful experience!

Hopefully the ``BlastErrorParser`` will make it much easier to debug and
deal with large Blast files.

7.6  Dealing with PSI-BLAST
---------------------------

You can run the standalone version of PSI-BLAST (the legacy NCBI command
line tool ``blastpgp``, or its replacement ``psiblast``) using the
wrappers in ``Bio.Blast.Applications`` module.

At the time of writing, the NCBI do not appear to support tools running
a PSI-BLAST search via the internet.

Note that the ``Bio.Blast.NCBIXML`` parser can read the XML output from
current versions of PSI-BLAST, but information like which sequences in
each iteration is new or reused isn’t present in the XML file. If you
care about this information you may have more joy with the plain text
output and the ``PSIBlastParser`` in ``Bio.Blast.NCBIStandalone``.

7.7  Dealing with RPS-BLAST
---------------------------

You can run the standalone version of RPS-BLAST (either the legacy NCBI
command line tool ``rpsblast``, or its replacement with the same name)
using the wrappers in ``Bio.Blast.Applications`` module.

At the time of writing, the NCBI do not appear to support tools running
an RPS-BLAST search via the internet.

You can use the ``Bio.Blast.NCBIXML`` parser to read the XML output from
current versions of RPS-BLAST.
