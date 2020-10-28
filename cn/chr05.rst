.. _chapter-Bio.SeqIO:

第5章  序列输入和输出
================================

本章将详细讨论 ``Bio.SeqIO`` 模块，该模块在第 :ref:`2 <chapter-quick-start>` 章已经做过简单的介绍并在第 :ref:`4 <chapter-SeqRecord>` 章使用过，它旨在提供一个简单的接口，实现对各种不同格式序列文件进行统一的处理。详细信息请查阅 ``Bio.SeqIO`` 维基页面（ `http://biopython.org/wiki/SeqIO <http://biopython.org/wiki/SeqIO>`__ ）和内置文档（ `SeqIO <http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html>`__ ）:

.. code:: python

    >>> from Bio import SeqIO
    >>> help(SeqIO)
    ...

学习本章的要领是学会使用 ``SeqRecord`` 对象（请见第 :ref:`4 <chapter-SeqRecord>` 章），该对象包含一个 ``Seq`` 对象（请见第 :ref:`3 <chapter-Bio.Seq>` 章）和注释信息（如序列ID和描述信息）。

5.1 解析/读取序列
---------------------------------

该模块的主要函数是 ``Bio.SeqIO.parse()`` ，它用于读取序列文件生成 ``SeqRecord`` 对象，包含两个参数：

#. 第一个参数是一个文件名或者一个句柄（ *handle* ）。句柄可以是打开的文件，命令行程序的输出，或者来自下载的数据(请见第 :ref:`5.3 <sec-SeqIO_Online>` 节)。更多关于句柄的信息请见第 :ref:`22.1 <sec-appendix-handles>` 节。
#. 第二个参数是一个小写字母字符串，用于指定序列格式（我们并不推测文件格式！），支持的文件格式请见 `http://biopython.org/wiki/SeqIO <http://biopython.org/wiki/SeqIO>`__ 。

``Bio.SeqIO.parse()`` 函数返回一个 ``SeqRecord`` 对象迭代器（ *iterator* ），迭代器通常用在循环中。

有时你需要处理只包含一个序列条目的文件，此时请使用函数 ``Bio.SeqIO.read()`` 。它使用与函数 ``Bio.SeqIO.parse()`` 相同的参数，当文件有且仅有一个序列条目时返回一个 ``SeqRecord`` 对象，否则触发异常。

5.1.1 读取序列文件
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

总的来说， ``Bio.SeqIO.parse()`` 用于读取序列文件并返回 ``SeqRecord`` 对象，并通常用在循环中，如：

.. code:: python

    from Bio import SeqIO
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

上面的示例来自第 :ref:`2.4 <sec-sequence-parsing>` 节，它将读取来自FASTA格式文件 `ls\_orchid.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk>`__ 的兰花DNA序列。如果你想读取GenBank格式文件，如 `ls\_orchid.gbk <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk>`__ ，只需要更改文件名和格式字符串：

.. code:: python

    from Bio import SeqIO
    for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
        print(seq_record.id)
        print(seq_record.seq)
        print(len(seq_record))

同样地，如果需要读取其他格式文件，并且 ``Bio.SeqIO.parse()`` 支持该文件格式，你只需要修改到相应的格式字符串，如“swiss”为SwissProt格式文件，“embl”为EMBL格式文本文件。详细的清单请见维基页面（ `http://biopython.org/wiki/SeqIO <http://biopython.org/wiki/SeqIO>`__ ）和内置文档（ `在线文档 <http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html>`__ ）。

另外一个非常常见的使用Python迭代器的地方是在列表解析（list comprehension，或者生成器表达式generator expression）。例如，如果需要从文件中提取序列ID列表，我们可以通过以下的列表推导很容易地实现：

.. code:: python

    >>> from Bio import SeqIO
    >>> identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank")]
    >>> identifiers
    ['Z78533.1', 'Z78532.1', 'Z78531.1', 'Z78530.1', 'Z78529.1', 'Z78527.1', ..., 'Z78439.1']

更多关于 ``SeqIO.parse()`` 在列表推导中运用的示例请见第 :ref:`18.2 <sec-sequence-parsing-plus-pylab>` 节（e.g. 对序列长度或GC%作图）。

5.1.2 遍历序列文件
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

在上述示例中，我们通常使用for循环遍历所有的序列条目（records）。你可以对for循环使用所有类型的支持迭代接口的Python对象（包括列表，元组（tuple）和字符串）。

``Bio.SeqIO`` 返回的对象实际上是一个返回 ``SeqRecord`` 对象的迭代器。你将顺序地获得每个条目，但是有且仅有一次；优势是，当处理大文件时，迭代器可以有效地节约内存空间。

除了使用for循环，还可以使用迭代器的 ``.next()`` 方法遍历序列条目，如：

.. code:: python

    from Bio import SeqIO
    record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")

    first_record = record_iterator.next()
    print(first_record.id)
    print(first_record.description)

    second_record = record_iterator.next()
    print(second_record.id)
    print(second_record.description)

注意：如果使用 ``.next()`` 方法，当没有序列条目时，将抛出 ``StopIteration`` 异常。

一种特殊情形是，序列文件包含多个序列条目，而你只需要第一个条目。在这种情况下，可使用以下代码，非常简洁：

.. code:: python

    from Bio import SeqIO
    first_record  = next(SeqIO.parse("ls_orchid.gbk", "genbank"))

注意：像上述示例中使用 ``.next()`` 方法将忽略文件中其余的序列。如果序列文件“有且仅有”一条序列条目，如本章后面的某些在线示例、包含单条染色体序列的GenBank文件，请使用 ``Bio.SeqIO.read()`` 函数。该函数会检查文件是否包含额外的序列条目。

5.1.3  获得序列文件中序列条目列表
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

在上一节中，我们讨论了如何使用 ``Bio.SeqIO.parse()`` 返回一个 ``SeqRecord`` 迭代器，然后顺序地获取序列条目。往往我们需要以任意顺序获取序列条目，Python列表数据类型便可以达到这个目的。使用Python内置函数 ``list()`` ，我们可以将序列条目迭代器转变成 ``SeqRecord`` 对象列表，如下：

.. code:: python

    from Bio import SeqIO
    records = list(SeqIO.parse("ls_orchid.gbk", "genbank"))

    print("Found %i records" % len(records))

    print("The last record")
    last_record = records[-1] #using Python's list tricks
    print(last_record.id)
    print(repr(last_record.seq))
    print(len(last_record))

    print("The first record")
    first_record = records[0] #remember, Python counts from zero
    print(first_record.id)
    print(repr(first_record.seq))
    print(len(first_record))

运行结果:

.. code:: python

    Found 94 records
    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592
    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740

当然，你仍然可以对 ``SeqRecord`` 对象列表使用for循环。使用列表比使用迭代器灵活得多（例如，可以根据列表大小知道序列条目数量），但缺点是for循环要同时读取所有的内容，需要更多的内存空间。

5.1.4 提取数据
~~~~~~~~~~~~~~~~~~~~~~

``SeqRecord`` 对象及其注释信息在第 :ref:`4 <chapter-SeqRecord>` 章中有更详细的介绍。为了解释注释信息是如何存储的，我们从GenBank文件 `ls\_orchid.gbk <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk>`__ 中解析出第一个序列条目，并将其输出：

.. code:: python

    from Bio import SeqIO
    record_iterator = SeqIO.parse("ls_orchid.gbk", "genbank")
    first_record = next(record_iterator.next())
    print(first_record)

输出结果:

.. code:: python

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA.
    Number of features: 5
    /sequence_version=1
    /source=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', ..., 'Cypripedium']
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', ..., 'ITS1', 'ITS2']
    /references=[...]
    /accessions=['Z78533']
    /data_file_division=PLN
    /date=30-NOV-2006
    /organism=Cypripedium irapeanum
    /gi=2765658
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')

这可以得到 ``SeqRecord`` 大部分的易读的注释汇总信息。在此例中，我们将使用 ``.annotations`` 属性-即Python字典（dictionary）。该注释字典的内容如上述示例结果，你也可以直接输出：

.. code:: python

    print(first_record.annotations)

与其他Python字典一样，你可以轻松地获得键列表：

.. code:: python

    print(first_record.annotations.keys())

或者值列表:

.. code:: python

    print(first_record.annotations.values())

通常，注释值是字符串或者字符串列表。一个特例是，文件中的所有参考文献（references）都以引用（reference）对象方式存储。

例如你想从GenBank文件 `ls\_orchid.gbk <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk>`__ 中提取出物种列表。我们需要的信息 *Cypripedium irapeanum* 被保存在这个注释字典的‘source’和‘organism’键中，我们可以用下面的方式获取：

.. code:: python

    >>> print(first_record.annotations["source"])
    Cypripedium irapeanum

或:

.. code:: python

    >>> print(first_record.annotations["organism"])
    Cypripedium irapeanum

通常，‘organism’ 用于学名（拉丁名，e.g. *Arabidopsis thaliana* ），而 ‘source’ 用于俗名（common name）（e.g. thale cress）。在此例中，以及在通常情况下，这两个字段是相同的。

现在，让我们遍历所有的序列条目， 创建一个包含所有兰花序列的物种列表：

.. code:: python

    from Bio import SeqIO
    all_species = []
    for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
        all_species.append(seq_record.annotations["organism"])
    print(all_species)

另外一种方式是使用列表解析：

.. code:: python

    from Bio import SeqIO
    all_species = [seq_record.annotations["organism"] for seq_record in \
                   SeqIO.parse("ls_orchid.gbk", "genbank")]
    print(all_species)

两种方式的输出结果相同：

.. code:: python

    ['Cypripedium irapeanum', 'Cypripedium californicum', ..., 'Paphiopedilum barbatum']

因为GenBank文件注释是以标准方式注释，所以相当简单。

现在，假设你需要从一个FASTA文件而不是GenBank文件提取出物种列表，那么你不得不多写一些代码，用以从序列条目的描述行提取需要的数据。使用的示例FASTA文件 `ls\_orchid.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk>`__ 格式如下：

.. code:: python

    >gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTG
    AATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGTGACCCTGATTTGTTGTTGGG
    ...

你可以手动检查，对于每一个序列条目，物种名都是描述行的第二个单词。这意味着如果我们以空白分割序列条目的 ``.description`` ，物种名将会是第1个元素（第0个元素是序列ID），我们可以这样做：

.. code:: python

    from Bio import SeqIO
    all_species = []
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        all_species.append(seq_record.description.split()[1])
    print(all_species)

将得到:

.. code:: python

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', ..., 'P.barbatum']

使用更简洁的列表解析：

.. code:: python

    from Bio import SeqIO
    all_species == [seq_record.description.split()[1] for seq_record in \
                    SeqIO.parse("ls_orchid.fasta", "fasta")]
    print(all_species)

通常，对FASTA描述行提取信息不是那么方便。如果你能获得对目标序列注释很好的文件格式如GenBank或者EMBL，那么这类注释信息就很容易处理。

5.1.5 修改数据
~~~~~~~~~~~~~~~~~~~~~~

在上一节中，我们演示了如何提取SeqRecord。另一个常见任务是更改此数据。一个SeqRecord的属性可以直接修改，例如：

.. code:: python

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
    >>> first_record = next(record_iterator)
    >>> first_record.id
    'gi|2765658|emb|Z78533.1|CIZ78533'
    >>> first_record.id = "new_id"
    >>> first_record.id
    'new_id'

注意，如果您想要改变写入文件时输出FASTA格式(见第 :ref:`5.5 <writing-sequence-files>` 章），应该同时修改 ``id`` 和 ``description`` 属性。为了确保正确的行为，最好在所需的 ``description`` 的开头添加一个 ``id`` 和空格：

.. code:: python

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
    >>> first_record = next(record_iterator)
    >>> first_record.id = "new_id"
    >>> first_record.description = first_record.id + " " + "desired new description"
    >>> print(first_record.format("fasta")[:200])
    >new_id desired new description
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGT

.. _sec-SeqIO_compressed:

5.2 从压缩文档读取解析序列信息
--------------------------------------------

在上一节中，我们研究了从文件中解析序列信息。除了使用文件名，你可以让 ``Bio.SeqIO`` 使用文件句柄（请见第 :ref:`22.1 <sec-appendix-handles>` 节）。在这一节，我们将使用文件句柄从压缩文件中解析序列信息。

正如你上面看到的，我们可以使用文件名作为 ``Bio.SeqIO.read()`` 或 ``Bio.SeqIO.parse()`` 的参数 - 例如在这个例子中，我们利用生成器表达式计算GenBank文件中多条序列条目的总长：

.. code:: python

    >>> from Bio import SeqIO
    >>> print(sum(len(r) for r in SeqIO.parse("ls_orchid.gbk", "gb")))
    67518

此处，我们使用文件句柄，并使用 ``with`` 语句自动关闭句柄：

.. code:: python

    >>> from __future__ import with_statement #Needed on Python 2.5
    >>> from Bio import SeqIO
    >>> with open("ls_orchid.gbk") as handle:
    ...     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))
    67518

或者，用旧版本的方式，手动关闭句柄：

.. code:: python

    >>> from Bio import SeqIO
    >>> handle = open("ls_orchid.gbk")
    >>> print(sum(len(r) for r in SeqIO.parse(handle, "gb")))
    67518
    >>> handle.close()

现在，如果我们有一个gzip压缩的文件呢？这种类型的文件在Linux系统中被普遍使用。我们可以使用Python的 ``gzip`` 模块打开压缩文档以读取数据 - 返回一个句柄对象：

.. code:: python

    >>> import gzip
    >>> from Bio import SeqIO
    >>> with gzip.open("ls_orchid.gbk.gz", "rt") as handle:
    ...     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))
    ...
    67518

相同地，如果我们有一个bzip2压缩文件：

.. code:: python

    >>> import bz2
    >>> from Bio import SeqIO
    >>> with bz2.open("ls_orchid.gbk.bz2", "rt") as handle:
    ...     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))
    ...
    67518

有一种gzip（GNU zip）变种称为BGZF（Blocked GNU Zip Format），它可以作为普通gzip文件被读取，但具有随机读取的优点，我们将在稍后的第 :ref:`5.4.4 <sec-SeqIO-index-bgzf>` 节讨论。

.. _sec-SeqIO_Online:

5.3 解析来自网络的序列
-----------------------------------

在上一节中，我们研究了从文件（使用文件名或者文件句柄）和压缩文件（使用文件句柄）解析序列数据。这里我们将使用 ``Bio.SeqIO`` 的另一种类型句柄，网络连接，从网络下载和解析序列。

请注意，你可以一气呵成地下载序列并解析成为 ``SeqRecord`` 对象，这并不意味这是一个好主意。通常，你可能需要下载序列并存入文件以重复使用。

.. _sec-SeqIO_GenBank_Online:

5.3.1 解析来自网络的GenBank序列条目
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

第 :ref:`9.6 <sec-efetch>` 节将更详细地讨论Entrez EFetch接口，但是现在我们将通过它连接到NCBI，通过GI号从GenBank获得 *Opuntia* （刺梨）序列。

首先，我们只获取一条序列条目。如果你不关注注释和相关信息，下载FASTA文件是个不错的选择，因为他们相对紧凑。请记住，当你希望处理的对象包含有且仅有一条序列条目时，使用 ``Bio.SeqIO.read()`` 函数：

.. code:: python

    from Bio import Entrez
    from Bio import SeqIO

    Entrez.email = "A.N.Other@example.com"
    with Entrez.efetch(
        db="nucleotide", rettype="fasta", retmode="text", id="6273291"
    ) as handle:
        seq_record = SeqIO.read(handle, "fasta")
    print("%s with %i features" % (seq_record.id, len(seq_record.features)))

输出结果为:

.. code:: python

    gi|6273291|gb|AF191665.1|AF191665 with 0 features

NCBI也允许你获取其它格式文件，尤其是GenBank文件。直到2009年复活节，Entrez EFetch API使用“genbank”作为返回类型。然而NCBI现在坚持使用“gb” （蛋白使用“gp”）作为官方返回类型，具体描述参见 `EFetch for Sequence and other Molecular Biology Databases <http://www.ncbi.nlm.nih.gov/entrez/query/static/efetchseq_help.html>`__ 。因此，Biopython1.50及以后版本的 ``Bio.SeqIO`` 中，我们支持“gb”作为“genbank”的别名。

.. code:: python

    from Bio import Entrez
    from Bio import SeqIO

    Entrez.email = "A.N.Other@example.com"
    with Entrez.efetch(
        db="nucleotide", rettype="gb", retmode="text", id="6273291"
    ) as handle:
        seq_record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"
    print("%s with %i features" % (seq_record.id, len(seq_record.features)))

输出结果为：

.. code:: python

    AF191665.1 with 3 features

请注意，这次我们获得3个特征。

现在，让我们获取多个序列条目。这次句柄包含多条序列条目，因此我们必须使用 ``Bio.SeqIO.parse()`` 函数：

.. code:: python

    from Bio import Entrez
    from Bio import SeqIO

    Entrez.email = "A.N.Other@example.com"
    with Entrez.efetch(
        db="nucleotide", rettype="gb", retmode="text", id="6273291,6273290,6273289"
    ) as handle:
        for seq_record in SeqIO.parse(handle, "gb"):
            print("%s %s..." % (seq_record.id, seq_record.description[:50]))
            print(
                "Sequence length %i, %i features, from: %s"
                % (
                    len(seq_record),
                    len(seq_record.features),
                    seq_record.annotations["source"],
                )
            )

输出结果为：

.. code:: python

    AF191665.1 Opuntia marenae rpl16 gene; chloroplast gene for c...
    Sequence length 902, 3 features, from: chloroplast Opuntia marenae
    AF191664.1 Opuntia clavata rpl16 gene; chloroplast gene for c...
    Sequence length 899, 3 features, from: chloroplast Grusonia clavata
    AF191663.1 Opuntia bradtiana rpl16 gene; chloroplast gene for...
    Sequence length 899, 3 features, from: chloroplast Opuntia bradtianaa

更多关于 ``Bio.Entrez`` 模块的信息请见第 :ref:`9 <chapter-entrez>` 章，并阅读NCBI Entrez使用指南（第 :ref:`9.1 <sec-entrez-guidelines>` 节）。

.. _sec-SeqIO_ExPASy_and_SwissProt:

5.3.2 解析来自网络的SwissProt序列条目
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

现在我们使用句柄下载来自ExPASy的SwissProt文件，更深入的信息请见第 :ref:`10 <chapter-swiss_prot>` 章。如上面提到的，当你希望处理的对象包含有且仅有一条序列条目时，使用 ``Bio.SeqIO.read()`` 函数：

.. code:: python

    from Bio import ExPASy
    from Bio import SeqIO

    with ExPASy.get_sprot_raw("O23729") as handle:
        seq_record = SeqIO.read(handle, "swiss")
    print(seq_record.id)
    print(seq_record.name)
    print(seq_record.description)
    print(repr(seq_record.seq))
    print("Length %i" % len(seq_record))
    print(seq_record.annotations["keywords"])

如果网络连接正常，你将会得到：

.. code:: python

    O23729
    CHS3_BROFI
    RecName: Full=Chalcone synthase 3; EC=2.3.1.74; AltName: Full=Naringenin-chalcone synthase 3;
    Seq('MAPAMEEIRQAQRAEGPAAVLAIGTSTPPNALYQADYPDYYFRITKSEHLTELK...GAE')
    Length 394
    ['Acyltransferase', 'Flavonoid biosynthesis', 'Transferase']

5.4 序列文件作为字典
-----------------------------------

我们将介绍 ``Bio.SeqIO`` 模块中3个相关函数，用于随机读取多序列文件。这里需要权衡灵活性和内存使用。总之：

-   ``Bio.SeqIO.to_dict()`` 最灵活但内存占用最大 （请见第 :ref:`5.4.1 <sec-SeqIO-to-dict>` 节）。这基本上是一个辅助函数，用于建立Python ``字典`` ，每个条目以 ``SeqRecord`` 对象形式存储在内存中，允许你修改这些条目。
-   ``Bio.SeqIO.index()`` 处于中间水平，类似于只读字典，当需要时解析序列到 ``SeqRecord`` 对象（请见第 :ref:`5.4.2 <sec-SeqIO-index>` 节）。
-   ``Bio.SeqIO.index_db()`` 也类似于只读字典，但是将文件中的ID和文件偏移值存储到硬盘（SQLite3数据库），这意味着它对内存需求很低（请见第 :ref:`5.4.3 <sec-SeqIO-index-db>` 节），但会慢一点。

全面的概述请见讨论部分（第 :ref:`5.4.5 <sec-SeqIO-indexing-discussion>` 节）。

.. _sec-SeqIO-to-dict:

5.4.1 序列文件作为字典-在内存中
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

我们对兰花数据文件接下来的处理将用于展示如何对他们建立索引，以及使用Python的 ``dictionary``  数量类型（与Perl中hash类似）以类似于数据库的方式读取数据。这常用于从中等大小的文件中读取某些特定元素，形成一个很好的快速数据库。如果处理较大的文件，内存将是个问题，请见下面第 :ref:`5.4.2 <sec-SeqIO-index>` 节。

你可以使用 ``Bio.SeqIO.to_dict()`` 函数创建一个 ``SeqRecord`` 字典（在内存中）。默认会使用每条序列条目的ID（i.e.  ``.id`` 属性）作为键。让我们用GenBank文件试一试：

.. code:: python

    >>> from Bio import SeqIO
    >>> orchid_dict = SeqIO.to_dict(SeqIO.parse("ls_orchid.gbk", "genbank"))

``Bio.SeqIO.to_dict()`` 仅需一个参数，即能够得到 ``SeqRecord`` 对象的列表或生成器，这里我们使用 ``SeqIO.parse`` 函数输出。顾名思义， ``Bio.SeqIO.to_dict()`` 返回一个Python字典。

因为变量 ``orchid_dict``  是一个普通的Python字典，我们可以查看所有的键：

.. code:: python

    >>> len(orchid_dict)
    94

.. code:: python

    >>> list(orchid_dict.keys())
    ['Z78484.1', 'Z78464.1', 'Z78455.1', 'Z78442.1', 'Z78532.1', 'Z78453.1', ..., 'Z78471.1']

在python3中，字典方法（如“.keys()“ and “.values()“）是迭代器而不是列表。

如果你确实需要，你甚至可以一次性查看所有的序列条目：

.. code:: python

    >>> orchid_dict.values() #lots of output!
    ...

我们可以通过键读取单个 ``SeqRecord``  对象并操作改对象：

.. code:: python

    >>> seq_record = orchid_dict["Z78475.1"]
    >>> print(seq_record.description)
    P.supardii 5.8S rRNA gene and ITS1 and ITS2 DNA.
    >>> seq_record.seq
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')

因此，可以用我们的GenBank序列条目轻松地在内存中创建一个数据库（in memory “database”）。接下来我们将尝试使用FASTA文件。

值得注意的是，对有Python使用经验的人来说，可以轻松地创建一个类似的字典。然而，典型的字典构建方法不能很好地处理重复键的情况。使用 ``Bio.SeqIO.to_dict()`` 函数将明确检查重复键，如果发现任何重复键将引发异常并退出。

.. _seq-seqio-todict-functionkey:

5.4.1.1 指定字典键
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

使用上述相同的代码，仅将文件改为FASTA文件：

.. code:: python

    from Bio import SeqIO
    orchid_dict = SeqIO.to_dict(SeqIO.parse("ls_orchid.fasta", "fasta"))
    print(orchid_dict.keys())

这次键为：

.. code:: python

    ['gi|2765596|emb|Z78471.1|PDZ78471', 'gi|2765646|emb|Z78521.1|CCZ78521', ...
     ..., 'gi|2765613|emb|Z78488.1|PTZ78488', 'gi|2765583|emb|Z78458.1|PHZ78458']

这结果是之前在第 :ref:`2.4.1 <sec-fasta-parsing>` 节中我们解析的FASTA文件结果。如果你需要别的作为键，如登录号（Accession Number），可使用 ``SeqIO.to_dict()`` 的可选参数 ``key_function`` ，它允许你根据你的序列条目特点，自定义字典键。

首先，你必须写一个函数，当使用 ``SeqRecord`` 对象作为参数时，可以返回你需要的键（字符串）。通常，函数的细节依赖于你要处理的序列条目的特点。但是对于我们的兰花数据，我们只需要使用“管道”符号（|）切分ID并返回第四个条目（第三个元素）：

.. code:: python

    def get_accession(record):
        """"Given a SeqRecord, return the accession number as a string.
      
        e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "Z78488.1"
        """
        parts = record.id.split("|")
        assert len(parts) == 5 and parts[0] == "gi" and parts[2] == "emb"
        return parts[3]

然后我们可以将此函数赋与 ``SeqIO.to_dict()`` 函数用于构建字典：

.. code:: python

    from Bio import SeqIO
    orchid_dict = SeqIO.to_dict(SeqIO.parse("ls_orchid.fasta", "fasta"), key_function=get_accession)
    print(orchid_dict.keys())

最终可到到新的字典键：

.. code:: python

    >>> print(orchid_dict.keys())
    ['Z78484.1', 'Z78464.1', 'Z78455.1', 'Z78442.1', 'Z78532.1', 'Z78453.1', ..., 'Z78471.1']

不是太困难！

5.4.1.2 使用SEGUID校验和对字典建立索引
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

为了介绍另外一个 ``SeqRecord`` 对象字典的示例，我们将使用SEGUID校验和函数。这是一个相对较新的校验和，冲突非常罕见（i.e. 两条不同序列具有相同的校验和），相对CRC64校验和有所提升。

让我们再一次处理兰花GenBank文件：

.. code:: python

    from Bio import SeqIO
    from Bio.SeqUtils.CheckSum import seguid
    for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
        print(record.id, seguid(record.seq))

将得到：

.. code:: python

    Z78533.1 JUEoWn6DPhgZ9nAyowsgtoD9TTo
    Z78532.1 MN/s0q9zDoCVEEc+k/IFwCNF2pY
    ...
    Z78439.1 H+JfaShya/4yyAj7IbMqgNkxdxQ

现在，再次调用 ``Bio.SeqIO.to_dict()`` 函数 ``key_function`` 参数， ``key_function`` 参数需要一个函数将 ``SeqRecord`` 转变为字符串。我们不能直接使用`seguid() `` 函数，因为它需要 ``Seq`` 对象（或字符串）作为参数。不过，我们可以使用Python的 ``lambda`` 特性创建一个一次性（“one off”）函数，然后传递给 ``Bio.SeqIO.to_dict()`` ：

.. code:: python

    >>> from Bio import SeqIO
    >>> from Bio.SeqUtils.CheckSum import seguid
    >>> seguid_dict = SeqIO.to_dict(SeqIO.parse("ls_orchid.gbk", "genbank"),
    ...                             lambda rec : seguid(rec.seq))
    >>> record = seguid_dict["MN/s0q9zDoCVEEc+k/IFwCNF2pY"]
    >>> print(record.id)
    Z78532.1
    >>> print(record.description)
    C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA.

将会返回文件中第二个序列条目 ``Z78532.1`` 。

.. _sec-SeqIO-index:

5.4.2 序列文件作为字典 - 索引文件
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

之前众多示例试图解释的是使用 ``Bio.SeqIO.to_dict()`` 的灵活性。然而，因为它将所有的信息都存储在内存中，你能处理的文件大小受限于电脑的RAM。通常，这仅能处理一些小文件或中等大小文件。

对于更大的文件，应该考虑使用 ``Bio.SeqIO.index()`` ，工作原理上略有不同。尽管仍然是返回一个类似于字典的对象，它并不将所有的信息存储在内存中。相反，它仅仅记录每条序列条目在文件中的位置 - 当你需要读取某条特定序列条目时，它才进行解析。

让我们使用之前相同的GenBank文件作为示例：

.. code:: python

    >>> from Bio import SeqIO
    >>> orchid_dict = SeqIO.index("ls_orchid.gbk", "genbank")
    >>> len(orchid_dict)
    94

.. code:: python

    >>> orchid_dict.keys()
    ['Z78484.1', 'Z78464.1', 'Z78455.1', 'Z78442.1', 'Z78532.1', 'Z78453.1', ..., 'Z78471.1']

.. code:: python

    >>> seq_record = orchid_dict["Z78475.1"]
    >>> print(seq_record.description)
    P.supardii 5.8S rRNA gene and ITS1 and ITS2 DNA.
    >>> seq_record.seq
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    >>> orchid_dict.close()

注意： ``Bio.SeqIO.index()`` 不接受句柄参数，仅仅接受文件名。这有充分的理由，但是过于技术性。第二个参数是文件格式（与其它 ``Bio.SeqIO`` 函数一样的小写字符串）。你可以使用许多其他的简单的文件格式，包括FASTA和FASTQ文件（示例参见第 :ref:`18.1.11 <sec-fastq-indexing>` 节），但不支持比对文件格式，如PHYLIP或Clustal。最后有个可选参数，你可以指定字符集或者键函数。

下面是使用FASTA文件做的相同的示例 - 仅改变了文件名和格式：

.. code:: python

    >>> from Bio import SeqIO
    >>> orchid_dict = SeqIO.index("ls_orchid.fasta", "fasta")
    >>> len(orchid_dict)
    94
    >>> orchid_dict.keys()
    ['gi|2765596|emb|Z78471.1|PDZ78471', 'gi|2765646|emb|Z78521.1|CCZ78521', ...
     ..., 'gi|2765613|emb|Z78488.1|PTZ78488', 'gi|2765583|emb|Z78458.1|PHZ78458']

5.4.2.1 指定字典键
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

如果想使用与之前一样的键，像第 :ref:`5.4.1.1 <seq-seqio-todict-functionkey>` 节 ``Bio.SeqIO.to_dict()`` 示例，你需要写一个小函数，从FASTA ID（字符串）中匹配你想要的键：

.. code:: python

    def get_acc(identifier):
        """"Given a SeqRecord identifier string, return the accession number as a string.
      
        e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "Z78488.1"
        """
        parts = identifier.split("|")
        assert len(parts) == 5 and parts[0] == "gi" and parts[2] == "emb"
        return parts[3]

然后我们将此函数赋与 ``Bio.SeqIO.index()`` 函数用于构建字典：

.. code:: python

    >>> from Bio import SeqIO
    >>> orchid_dict = SeqIO.index("ls_orchid.fasta", "fasta", key_function=get_acc)
    >>> print(orchid_dict.keys())
    ['Z78484.1', 'Z78464.1', 'Z78455.1', 'Z78442.1', 'Z78532.1', 'Z78453.1', ..., 'Z78471.1']

当你知道怎样实现就变得很简单了。

.. _sec-seqio-index-getraw:

5.4.2.2 获取序列条目原始数据
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

来自 ``Bio.SeqIO.index()`` 的字典样对象以 ``SeqRecord`` 对象形式返回序列条目。但是，有时候从文件中直接获取原始数据非常有用。对于此种情况，使用 ``get_raw()`` 方法，它仅需要一个参数（序列ID），然后返回一个字符串（提取自文件的未处理数据）。

一个重要的例子就是从大文件中提取出一个序列子集，特别是当 ``Bio.SeqIO.write()`` 还不支持这种输出格式（e.g. SwissProt文件格式的文本文件 ） 或者需要完整地保留源文件（Biopython的GenBank和EMBL格式输出并不会保留每一点注释信息）。

假如你已经从UniProt FTP站点下载了整个数据库的SwissPort格式文本文件（ `ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz>`__ ），并也已经解压为文件 ``uniprot_sprot.dat`` ，你需要从中提取一部分序列条目：

.. code:: python

    >>> from Bio import SeqIO
    >>> uniprot = SeqIO.index("uniprot_sprot.dat", "swiss")
    >>> with open("selected.dat", "wb") as out_handle:
    ...     for acc in ["P33487", "P19801", "P13689", "Q8JZQ5", "Q9TRC7"]:
    ...         out_handle.write(uniprot.get_raw(acc))
    ...

请注意，从Python 3开始，我们必须打开文件以二进制模式进行写入，因为 ``get_raw()`` 方法返回字节字符串。

在第 :ref:`18.1.5 <sec-SeqIO-sort>` 节有更多关于使用 ``SeqIO.index()`` 函数对大文件序列排序的示例（不需要一次加载所有信息到内存）。

.. _sec-SeqIO-index-db:

5.4.3 序列文件作为字典 - 数据库索引文件
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Biopython 1.57引入一个替代的函数， ``Bio.SeqIO.index_db()`` 。由于它将序列信息以文件方式存储在硬盘上（使用SQLite3数据库）而不是内存中，因此它可以处理超大文件。同时，你可以同时对多个文件建立索引（前提是所有序列条目的ID是唯一的）。

``Bio.SeqIO.index()`` 函数有三个参数：

-  索引文件名，我们建议使用以 ``.idx`` 结尾的字符，改索引文件实质上是SQLite3数据库；
-  要建立索引的文件列表（或者单个文件名）；
-  文件格式（与 ``SeqIO`` 模块中其它函数一样的小写字符串）。

将以NCBI FTP站点 `ftp://ftp.ncbi.nih.gov/genbank/ <ftp://ftp.ncbi.nih.gov/genbank/>`__ 的GenBank文本文件为例，这些文件为gzip压缩文件。对于GenBank版本210，病毒序列共包含38个文件， ``gbvrl1.seq``  -  ``gbvrl138.seq`` ，解压缩后大约占用8GB磁盘空间，总共包含近200万条记录。

如果您对病毒感兴趣，则可以使用 ``rsync`` 命令非常轻松地从命令行下载所有病毒文件，然后使用 ``gunzip`` 解压缩它们：

.. code:: python

    # For illustration only, see reduced example below
    $ rsync -avP "ftp.ncbi.nih.gov::genbank/gbvrl*.seq.gz" .
    $ gunzip gbvrl*.seq.gz

除非您关心病毒，否则仅此示例需要下载大量数据-因此，让我们仅下载前四个块（每个压缩块约25MB），然后解压缩（占用全部1GB的空间）：

.. code:: python

    # Reduced example, download only the first four chunks
    $ curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl1.seq.gz
    $ curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl2.seq.gz
    $ curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl3.seq.gz
    $ curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl4.seq.gz
    $ gunzip gbvrl*.seq.gz

现在，在Python中，按如下所示索引这些GenBank文件：

.. code:: python

    >>> import glob
    >>> from Bio import SeqIO
    >>> files = glob.glob("gbvrl*.seq")
    >>> print("%i files to index" % len(files))
    4
    >>> gb_vrl = SeqIO.index_db("gbvrl.idx", files, "genbank")
    >>> print("%i sequences indexed" % len(gb_vrl))
    272960 sequences indexed

在我的个人电脑上，对全套病毒GenBank文件进行索引大约需要十分钟，而仅前四个文件大约需要一分钟左右。但一旦完成，重复此操作将在不到一秒钟的时间内重新加载索引文件gbvrl.idx

您可以将索引用作只读的Python字典-不必担心序列来自哪个文件，例如：

.. code:: python

    >>> print(gb_vrl["GQ333173.1"].description)
    Equine encephalosis virus NS3 gene, complete cds, isolate: Kimron1.

5.4.3.1 获取序列条目原始数据
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

与第 :ref:`5.4.2.2 <sec-seqio-index-getraw>` 节讨论的 ``Bio.SeqIO.index()`` 函数一样，该字典样对象同样允许你获取每个序列条目的原始文件：

.. code:: python

    >>> print(gb_vrl.get_raw("AB811634.1"))
    LOCUS       AB811634                 723 bp    RNA     linear   VRL 17-JUN-2015
    DEFINITION  Equine encephalosis virus NS3 gene, complete cds, isolate: Kimron1.
    ACCESSION   AB811634
    ...
    //

.. _sec-SeqIO-index-bgzf:

5.4.4 对压缩文件建立索引
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

经常你要建立索引的文件可能非常大，因此你想对它进行压缩。不幸的是，对常规的文件格式如gzip和bzip2高效的随机读取通常很困难。在这种情况下，BGZF (Blocked GNU Zip Format)非常有用。它是gzip变体（也可以使用标准的gzip工具解压），因BAM文件格式得到推广，`samtools <http://samtools.sourceforge.net/>`__ 和 `tabix <http://samtools.sourceforge.net/tabix.shtml>`__ ；

你可以使用samtools的命令行工具 ``bgzip`` 创建BGZF格式压缩文件。在我们的示例中，使用文件扩展名 ``*.bgz`` ，以区分于普通的压缩文件（命名为 ``*.gz`` ）。你也可以在Python中使用 ``Bio.bgzf`` 模块读写BGZF文件。

``Bio.SeqIO.index()`` 和 ``Bio.SeqIO.index_db()`` 函数均可以用于BGZF压缩文件。例如，如果使用过未压缩的GenBank文件：

.. code:: python

    >>> from Bio import SeqIO
    >>> orchid_dict = SeqIO.index("ls_orchid.gbk", "genbank")
    >>> len(orchid_dict)
    94
    >>> orchid_dict.close()

你可以使用如下的命令行命令压缩该文件（同时保留源文件） - 不需要担心，压缩文件和别的示例及已经包含：

.. code:: python

    $ bgzip -c ls_orchid.gbk > ls_orchid.gbk.bgz

你可以用相同的方式使用压缩文件：

.. code:: python

    >>> from Bio import SeqIO
    >>> orchid_dict = SeqIO.index("ls_orchid.gbk.bgz", "genbank")
    >>> len(orchid_dict)
    94
    >>> orchid_dict.close()

或：

.. code:: python

    >>> from Bio import SeqIO
    >>> orchid_dict = SeqIO.index_db("ls_orchid.gbk.bgz.idx", "ls_orchid.gbk.bgz", "genbank")
    >>> len(orchid_dict)
    94
    >>> orchid_dict.close()

``SeqIO`` 建立索引时自动检测是否为BGZF压缩格式。注意：压缩文件和未压缩文件不能使用相同的索引文件。

.. _sec-SeqIO-indexing-discussion:

5.4.5 讨论
~~~~~~~~~~~~~~~~~

这些方法你该使用哪种及其原因，取决于你要做什么（以及你要处理的数据有多大）。然而，通常 ``Bio.SeqIO.index()`` 是个不错的选择。如果你正在处理上百万条序列条目，多个文件，或者重复性分析，那么看看 ``Bio.SeqIO.index_db()`` 。

选择 ``Bio.SeqIO.to_dict()`` 而不选择 ``Bio.SeqIO.index()`` 或 ``Bio.SeqIO.index_db()`` 的原因主要是它的灵活性，尽管会占用更多内存。存储 ``SeqRecord`` 对象到内存的优势在于可以随意被改变，添加或者删除。除了高内存消耗这个缺点外，建立索引也可能花费更长的时间，因为所有的条目都需要被完全解析。

``Bio.SeqIO.index()`` 和 ``Bio.SeqIO.index_db()`` 都是在需要时才解析序列条目。当建立索引时，他们扫描文件，寻找每个序列条目的起始，并做尽可能少的工作提取出ID信息。

选择 ``Bio.SeqIO.index()`` 而不选择 ``Bio.SeqIO.index_db()`` 的原因包括以下：

-  建立索引更快（需要注意的是简单文件格式）
-  读取 ``SeqRecord`` 对象稍快（但是这种差异只有在解析简单格式文件是才可见）
-  可以使用不可变的Python对象作为字典键而不仅仅是字符串（e.g. 如字符串元组、不可变容器（frozen set））
-  如果被建立索引的序列文件改变，不需要担心索引数据库过期。

选择 ``Bio.SeqIO.index_db()`` 而不选择 ``Bio.SeqIO.index()`` 的原因包括以下：

-  没有内存限制 - 这对通常多达10亿的二代测序文件来说非常重要，如果使用 ``Bio.SeqIO.index()`` 可能需要超过4G的RAM和64位Python
-  索引数据量保存在硬盘上，可重复使用。尽管建立索引数据库需要花费更多的时间，但是从长远看来。如果你有个脚本重新运行这个相同的数据库，可以节约时间
-  可以同时对多个文件建立索引
-  `get_raw() `` 方法可以快得多，因为对于大多数文件格式只需要存储序列条目的长度和偏移量（offset）

5.5 写入序列文件
---------------------------

我们已经讨论了使用 ``Bio.SeqIO.parse()`` 输入序列（读取文件），现在我们将研究使用 ``Bio.SeqIO.write()`` 输出序列（写入文件）。该函数需要三个参数：某些 ``SeqRecord`` 对象，要写入的句柄或文件名，和序列格式。

我们先用硬编码方式（手动创建而不是从文件中加载）创建一个些新的 ``SeqRecord`` 对象，示例如下：

.. code:: python

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rec1 = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD" \
                        +"GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK" \
                        +"NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM" \
                        +"SSAC"),
                     id="gi|14150838|gb|AAK54648.1|AF376133_1",
                     description="chalcone synthase [Cucumis sativus]")

    rec2 = SeqRecord(Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ" \
                        +"DMVVVEIPKLGKEAAVKAIKEWGQ"),
                     id="gi|13919613|gb|AAK33142.1|",
                     description="chalcone synthase [Fragaria vesca subsp. bracteata]")

    rec3 = SeqRecord(Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC" \
                        +"EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP" \
                        +"KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN" \
                        +"NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV" \
                        +"SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW" \
                        +"IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT" \
                        +"TGEGLEWGVLFGFGPGLTVETVVLHSVAT"),
                     id="gi|13925890|gb|AAK49457.1|",
                     description="chalcone synthase [Nicotiana tabacum]")
                   
    my_records = [rec1, rec2, rec3]

现在我们得到一个 ``SeqRecord`` 对象列表，将它写入一个FASTA格式文件：

.. code:: python

    from Bio import SeqIO
    SeqIO.write(my_records, "my_example.faa", "fasta")

如果用你喜欢的文本编辑软件打开，可得到：

.. code:: python

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    >gi|13919613|gb|AAK33142.1| chalcone synthase [Fragaria vesca subsp. bracteata]
    YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ
    DMVVVEIPKLGKEAAVKAIKEWGQ
    >gi|13925890|gb|AAK49457.1| chalcone synthase [Nicotiana tabacum]
    MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC
    EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP
    KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN
    NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV
    SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW
    IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT
    TGEGLEWGVLFGFGPGLTVETVVLHSVAT

怎样才能知道 ``Bio.SeqIO.write()`` 函数写入了多少条序列条目到句柄呢？如果你的序列条目保存在一个列表中，只需要使用 ``len(my_records)`` ，但是你不能对来自生成器/迭代器的序列条目。 ``Bio.SeqIO.write()`` 函数本身就返回写入文件的 ``SeqRecord`` 对象个数。

* 注意 - 如果你 ``Bio.SeqIO.write()`` 函数要写入的文件已经存在，旧文件将会被覆写，并且不会得到任何警告信息。

5.5.1 可逆读写（Round trips）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

某些人需要他们的解析器是“可逆”的，即当你读入某个文件后和可以按原样写回。这需要解析器提取足够多的信息用于 * 精确 * 还原原始文件， ``Bio.SeqIO`` 不打算这么做。

一个简单的例子是，FASTA文件中，允许序列以任意字符数换行。解析以下两条序列得到一个相同的 ``SeqRecord`` 对象，这两条序列仅在换行上不同：

.. code:: python

    >YAL068C-7235.2170 Putative promoter sequence
    TACGAGAATAATTTCTCATCATCCAGCTTTAACACAAAATTCGCACAGTTTTCGTTAAGA
    GAACTTAACATTTTCTTATGACGTAAATGAAGTTTATATATAAATTTCCTTTTTATTGGA

    >YAL068C-7235.2170 Putative promoter sequence
    TACGAGAATAATTTCTCATCATCCAGCTTTAACACAAAATTCGCA
    CAGTTTTCGTTAAGAGAACTTAACATTTTCTTATGACGTAAATGA
    AGTTTATATATAAATTTCCTTTTTATTGGA

为了创建一个可逆读写的FASTA解析器，需要记录序列换行发生的位置，而这些额外的信息通常毫无意义。因此，Biopython在输出时使用默认的60字符换行。空白字符在许多其他文件格式中运用也存在相同的问题。另外一个问题是，在某些情况下，Biopython并不能保存每一点注释信息（e.g. GenBank和EMBL）。

少数时候，重要的是保留原来的布局（这可能有点怪异），第 :ref:`5.4.2.2 <sec-seqio-index-getraw>` 节关于 ``Bio.SeqIO.index()`` 字典样对象的 ``get_raw()`` 方法提供了可能的解决方案。

.. _sec-SeqIO-conversion:

5.5.2 序列格式间的转换
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

在之前的例子中我们使用 ``SeqRecord`` 对象列表作为 ``Bio.SeqIO.write()`` 函数的输入，但是它也接受如来自于 ``Bio.SeqIO.parse()`` 的 ``SeqRecord`` 迭代器 - 这允许我们通过结合使用这两个函数实现文件转换。

在这个例子中，我们将读取GenBank格式文件 `ls\_orchid.gbk <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk>`__ ，然后输出为FASTA格式文件：

.. code:: python

    from Bio import SeqIO
    records = SeqIO.parse("ls_orchid.gbk", "genbank")
    count = SeqIO.write(records, "my_example.fasta", "fasta")
    print("Converted %i records" % count)

这仍然有点复杂，因为文件格式转换是比较常见的任务，有一个辅助函数可以替代上述代码：

.. code:: python

    from Bio import SeqIO
    count = SeqIO.convert("ls_orchid.gbk", "genbank", "my_example.fasta", "fasta")
    print("Converted %i records" % count)

``Bio.SeqIO.convert()`` 函数可以使用句柄或文件名。然而需要注意的是，如果输出文件已存在，将覆写该文件。想了解更多信息，请使用内置帮助文档：

.. code:: python

    >>> from Bio import SeqIO
    >>> help(SeqIO.convert)
    ...

原理上讲，只需要改变文件名和格式字符串，该代码即可实现Biopython支持的文件格式间的转换。然而，写入某种格式时需要某些特定的信息（e.g. 质量值），而其他格式文件不包含此信息。例如，你可以将FASTQ转化为FASTA文件，却不能进行逆操作。不同FASTQ格式间的相互转变请见cookbook章第 :ref:`18.1.9 <sec-SeqIO-fastq-conversion>` 节和第 :ref:`18.1.10 <sec-SeqIO-fasta-qual-conversion>` 节。

最后，使用 ``Bio.SeqIO.convert()`` 函数额外的好处是更快，（最大的好处是代码会更短）原因是该转换函数可以利用几个文件格式特殊的优化条件和技巧。

.. _sec-SeqIO-reverse-complement:

5.5.3 转化序列到反向互补序列
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

假设你有一个核苷酸序列文件，需要转换成一个包含其反向互补的文件。这时，需要做些工作，将从文件得到的 ``SeqRecord`` 对象转化为适合存储到输出文件的信息。

首先，我们将使用 ``Bio.SeqIO.parse()`` 加载文件中的核酸序列，然后使用 ``Seq`` 对象的内置方法 ``.reverse_complement()`` 输出其反向互补序列（请见第 :ref:`3.6 <sec-seq-reverse-complement>` 节）。

.. code:: python

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    ...     print(record.id)
    ...     print(record.seq.reverse_complement())

现在，如果我们想保存这些反向互补序列到某个文件，需要创建 ``SeqRecord`` 对象。我们可以使用 ``SeqRecord`` 对象的内置方法 ``.reverse_complement()`` （请见第 :ref:`4.9 <sec-SeqRecord-reverse-complement>` 节），但是我们必须决定新的序列条目怎么命名。

这是一个绝好的展示列表解析效率地方，列表解析通过在内存中创建一个列表实现：

.. code:: python

    >>> from Bio import SeqIO
    >>> records = [rec.reverse_complement(id="rc_"+rec.id, description = "reverse complement") \
    ...            for rec in SeqIO.parse("ls_orchid.fasta", "fasta")]
    >>> len(records)

这时就用到了列表解析的绝妙之处，在其中添加一个条件语句：

.. code:: python

    >>> records = [rec.reverse_complement(id="rc_"+rec.id, description = "reverse complement") \
    ...            for rec in SeqIO.parse("ls_orchid.fasta", "fasta") if len(rec)<700]
    >>> len(records)
    18

这将在内存中创建一个序列小于700bp的反向互补序列列表。我们可以以相同的方式使用生成器表达式 - 但是更有优势的是，它不需要同时在内存中创建所有序列条目的列表：

.. code:: python

    >>> records = (rec.reverse_complement(id="rc_"+rec.id, description = "reverse complement") \
    ...           for rec in SeqIO.parse("ls_orchid.fasta", "fasta") if len(rec)<700)

完整的示例如下：

.. code:: python

    >>> from Bio import SeqIO
    >>> records = (rec.reverse_complement(id="rc_"+rec.id, description = "reverse complement") \
    ...            for rec in SeqIO.parse("ls_orchid.fasta", "fasta") if len(rec)<700)
    >>> SeqIO.write(records, "rev_comp.fasta", "fasta")
    18

在第 :ref:`18.1.3 <sec-SeqIO-translate>` 节有一个相关的示例，将FASTA文件中核酸序列翻译为氨基酸。

.. _sec-Bio.SeqIO-and-StringIO:

5.5.4 获得格式化为字符串的 ``SeqRecord`` 对象
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

有时你不需要将序列条目写入文件或者句柄，而是想获得包含特定格式序列条目的字符串。 ``Bio.SeqIO`` 接口基于句柄，但是Python有一个有用的内置模块，提供基于字符串的句柄。

举个例子来说明如果使用这个功能，我们先从兰花GenBank文件加载一系列 ``SeqRecord`` 对象，然后创建一个包含FASTA格式序列条目的字符串：

.. code:: python

    from Bio import SeqIO
    from StringIO import StringIO

    records = SeqIO.parse("ls_orchid.gbk", "genbank")
    out_handle = StringIO()
    SeqIO.write(records, out_handle, "fasta")
    fasta_data = out_handle.getvalue()
    print(fasta_data)

当你第一次看到，会觉得这并不够简单明了。在特殊情况下，你希望得到一个只包含特定格式的单条序列条目的字符串，可以使用 ``SeqRecord`` 类的 ``format()`` （请见第 :ref:`4.5 <sec-SeqRecord-format>` 节）。

注意：尽管我们不鼓励这么做，你可以使用 ``format()`` 方法写入文件，示例如下：

.. code:: python

    from Bio import SeqIO
    with open("ls_orchid_long.tab", "w") as out_handle:
        for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
            if len(record) > 100:
                out_handle.write(record.format("tab"))

这类代码可以处理顺序文件格式如FASTA或者此处使用的简单的制表符分割文件，但不能处理更复杂的或是交错式文件格式。这就是为什么我们仍然强调使用 ``Bio.SeqIO.write()`` 的原因，如下面的示例：

.. code:: python

    from Bio import SeqIO
    records = (rec for rec in SeqIO.parse("ls_orchid.gbk", "genbank") if len(rec) > 100)
    SeqIO.write(records, "ls_orchid.tab", "tab")

同时 ，单次调用 ``SeqIO.write(...)`` 也比多次调用 ``SeqRecord.format(...)`` 方法更快。

5.6 低级FASTA和FASTQ解析器
---------------------------

在处理速度很重要的大型高通量FASTA或FASTQ排序文件时，与``Bio.SeqIO.parse``相比，使用低级``SimpleFastaParser``或``FastqGeneralIterator``通常更实用。 如本章引言中所述，文件格式中立的``Bio.SeqIO``接口具有创建许多对象的开销，即使对于FASTA这样的简单格式也是如此。

解析FASTA文件时，在内部``Bio.SeqIO.parse()``使用文件句柄调用低级SimpleFastaParser。 您可以直接使用它-遍历文件句柄，以两个字符串的元组，标题行（>字符之后的所有字符）和序列（以纯字符串形式）返回每个记录：

.. code:: python

    >>> from Bio.SeqIO.FastaIO import SimpleFastaParser
    >>> count = 0
    >>> total_len = 0
    >>> with open("ls_orchid.fasta") as in_handle:
    ...     for title, seq in SimpleFastaParser(in_handle):
    ...         count += 1
    ...         total_len += len(seq)
    ...
    >>> print("%i records with total sequence length %i" % (count, total_len))
    94 records with total sequence length 67518

只要您不关心换行（并且您可能不期待短时间读取高吞吐量数据），那么从这些字符串输出FASTA格式也非常快：

.. code:: python

    ...
    out_handle.write(">%s\n%s\n" % (title, seq))
    ...

同样，在解析FASTQ文件时，在内部``Bio.SeqIO.parse()`` 会使用文件句柄调用低级``FastqGeneralIterator``。 如果您不需要将质量得分转换为整数，或者可以将其作为ASCII字符串使用，则可以选择：

.. code:: python

    >>> from Bio.SeqIO.QualityIO import FastqGeneralIterator
    >>> count = 0
    >>> total_len = 0
    >>> with open("example.fastq") as in_handle:
    ...     for title, seq, qual in FastqGeneralIterator(in_handle):
    ...         count += 1
    ...         total_len += len(seq)
    ...
    >>> print("%i records with total sequence length %i" % (count, total_len))
    3 records with total sequence length 75

Cookbook（第 :ref:`20 <sec-cool-things> 章）中有更多示例，包括如何使用以下代码片段从字符串有效地输出FASTQ：

.. code:: python

    ...
    out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    ...

