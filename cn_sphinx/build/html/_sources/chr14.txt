第14章   使用Bio.motifs进行模体序列分析
====================================================

这章主要的介绍Biopython中的 ``Bio.motifs`` 包。这个包是为了方便那些需要进行模体序列分析的人们而特意提供的，所以我想你们在使用时肯定对模体序列分析的一些相关要点都很熟悉。假如在使用中遇到不清楚的地方，请您查阅 :ref:`14.8 <sec-links>` 相关章节以获得有关的信息。

这章的大部分内容是介绍Biopython 1.61 之前版本中新加入的 ``Bio.motifs`` 包，该包替代了Biopython 1.50版本中的 ``Bio.Motif`` 包，而 ``Bio.Motif`` 包是基于较早版本的Biopython 中的两个模块 ``Bio.AlignAce`` 和 ``Bio.MEME`` 。``Bio.motifs`` 包较好地综合了上述的几个模块的功能，做为一个统一模块工具。

说到其他库，看到这里，你或许会对 `TAMO <http://fraenkel.mit.edu/TAMO/>`__ 感兴趣，这是另一个分析模体序列的Python库。它能提供更多关于 *de-novo* 模体的查找方式，不过它并没有纳入到Biopython中，而且在商业用途上还有一些限制。

14.1  模体对象
-------------------

由于我们感兴趣的是模体分析，所以我们需要先看看 ``Motif`` 对象。对此我们需要先导入Bio.motifs包：

.. code:: python

    >>> from Bio import motifs

然后我们可以开始创建我们第一个模体对象。我们可以从模体的实例列表中创建一个 ``Motif`` 对象，也可以通过读取模体数据库中或模体查找软件产生的文件来获得一个 ``Motif`` 对象。

.. code:: python

    >>> from Bio import motifs


14.1.1  从实例中创建一个模体
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

假设我们有一些DNA模体的实例：

.. code:: python

    >>> from Bio.Seq import Seq
    >>> instances = [Seq("TACAA"),
    ...              Seq("TACGC"),
    ...              Seq("TACAC"),
    ...              Seq("TACCC"),
    ...              Seq("AACCC"),
    ...              Seq("AATGC"),
    ...              Seq("AATGC"),
    ...             ]

然后我们可以如下创建一个模体对象：

.. code:: python

    >>> m = motifs.create(instances)

这些实例被存储在一个名为 ``m.instances`` 的属性中，这个其实也就是一个Python的列表，只不过附加了一些功能，这些功能将在之后介绍。将这些模体对象打印出来后就可以看出这些实例是从哪构建出来的。

.. code:: python

    >>> print m
    TACAA
    TACGC
    TACAC
    TACCC
    AACCC
    AATGC
    AATGC
    <BLANKLINE>

模体的长度像其他一些实例一些被定义为序列的长度：

.. code:: python

    >>> len(m)
    5

模体对象有一个 ``.counts`` 属性，可以用来查看碱基在每个位置的数目。可以把这个统计表用易读的格式打印出来：

.. code:: python

    >>> print m.counts
            0      1      2      3      4
    A:   3.00   7.00   0.00   2.00   1.00
    C:   0.00   0.00   5.00   2.00   6.00
    G:   0.00   0.00   0.00   3.00   0.00
    T:   4.00   0.00   2.00   0.00   0.00
    <BLANKLINE>

你也可以像使用字典一样获取这些数目：

.. code:: python

    >>> m.counts['A']
    [3, 7, 0, 2, 1]

但是你也可以把它看成一个二维数列，核苷酸作为列，位置作为行：

.. code:: python

    >>> m.counts['T',0]
    4
    >>> m.counts['T',2]
    2
    >>> m.counts['T',3]
    0

你还可以直接获得核苷酸数目矩阵中的列

.. code:: python

    >>> m.counts[:,3]
    {'A': 2, 'C': 2, 'T': 0, 'G': 3}

除了使用核苷酸本身，你还可以使用模体碱基序列按字符排序后的核苷酸索引：

.. code:: python

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

模体有一个相关联的一致序列，这个序列被定义为由 ``.counts`` 矩阵相应列中具有最大值的碱基，这些碱基是按模体序列排列的：

.. code:: python

    >>> m.consensus
    Seq('TACGC', IUPACUnambiguousDNA())

反一致序列也一样，只不过是由 ``.counts`` 矩阵中相应列的最小值来选：

.. code:: python

    >>> m.anticonsensus
    Seq('GGGTG', IUPACUnambiguousDNA())

你也可以利用简并一致序列，用不确定核苷酸来表示序列某一位置的所有核苷酸：

.. code:: python

    >>> m.degenerate_consensus
    Seq('WACVC', IUPACAmbiguousDNA())

此处，W和R都是按照IUPAC不确定核苷酸表规定的：W代表A或T，V代表A，C或G [:ref:`10 <cornish1985>`] 。这些简并一致序列是按照Cavener指定的规则 [:ref:`11 <cavener1987>`] 来建立的。

.. code:: python

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

反向互补序列和简并一致序列都只在DNA模体中有。

14.1.2  读取模体
~~~~~~~~~~~~~~~~~~~~~~

从实例手动创建一个模体确实有点无趣，所以用一些I/O函数来读写模体是很有用的。目前对于如何存储模体还没有一些真正的标准，不过有一些格式用得比其他更经常。这其中最重要的区别在于模体表示是基于实例还是某种PWM矩阵。

JASPAR
^^^^^^

作为一个最流行的模体数据库 `JASPAR <http://jaspar.genereg.net>`__ 它不是以一系列的实例就是频率矩阵。比如，下面就是JASPAR ``Arnt.sites`` 文件的开头和结尾行显示了老鼠螺旋-环-螺旋转录因子Arnt的结合位点：


.. code:: python

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

那些用大字字母表示的序列的一部分就是被用来相互比对的模体实例。

我们可以从下面的实例创建一个 ``Motif`` 对象：

.. code:: python

    >>> from Bio import motifs
    >>> arnt = motifs.read(open("Arnt.sites"), "sites")

从这个模体创建的实例存储在该模体的 ``.instances`` 属性：

.. code:: python

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

这个模体的计数矩阵可以从这些实例中自动计算出来：

.. code:: python

    >>> print arnt.counts
            0      1      2      3      4      5
    A:   4.00  19.00   0.00   0.00   0.00   0.00
    C:  16.00   0.00  20.00   0.00   0.00   0.00
    G:   0.00   1.00   0.00  20.00   0.00  20.00
    T:   0.00   0.00   0.00   0.00  20.00   0.00
    <BLANKLINE>

JASPAR数据库也可以让模体像计数矩阵一样获得，不需要那些创建它们的实例。比如，下面这个JASPAR文件 ``SRF.pfm`` 包含了人类SRF转录因子的计数矩阵：

.. code:: python

     2  9  0  1 32  3 46  1 43 15  2  2
     1 33 45 45  1  1  0  0  0  1  0  1
    39  2  1  0  0  0  0  0  0  0 44 43
     4  2  0  0 13 42  0 45  3 30  0  0

我们可以如下为计数矩阵创建一个模体：

.. code:: python

    >>> srf = motifs.read(open("SRF.pfm"),"pfm")
    >>> print srf.counts
            0      1      2      3      4      5      6      7      8      9     10     11
    A:   2.00   9.00   0.00   1.00  32.00   3.00  46.00   1.00  43.00  15.00   2.00   2.00
    C:   1.00  33.00  45.00  45.00   1.00   1.00   0.00   0.00   0.00   1.00   0.00   1.00
    G:  39.00   2.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  44.00  43.00
    T:   4.00   2.00   0.00   0.00  13.00  42.00   0.00  45.00   3.00  30.00   0.00   0.00
    <BLANKLINE>

由于这个模体是由计数矩阵直接创建的，所以它没有相关的实例：

.. code:: python

    >>> print srf.instances
    None

我们可以获得这两个模体的一致序列：

.. code:: python

    >>> print arnt.counts.consensus
    CACGTG
    >>> print srf.counts.consensus
    GCCCATATATGG

MEME
^^^^

MEME [:ref:`12 <bailey1994>`] 是一个用来在一堆相关DNA或蛋白质序列中发现模体的工具。它输入一组相关DNA或蛋白质序列，输出所要求的模体。因此和JASPAR文件相比，MEME输出文件里面一般是含有多个模体。例子如下。

在输出文件的开头，有一些MEME生成的关于MEME和所用MEME版本的背景信息：

.. code:: python

    ********************************************************************************
    MEME - Motif discovery tool
    ********************************************************************************
    MEME version 3.0 (Release date: 2004/08/18 09:07:01)
    ...

再往下，简要概括了输入的训练序列集：

.. code:: python

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

以及所使用到的命令：

.. code:: python

    ********************************************************************************
    COMMAND LINE SUMMARY
    ********************************************************************************
    This information can also be useful in the event you wish to report a
    problem with the MEME software.

    command: meme -mod oops -dna -revcomp -nmotifs 2 -bfile yeast.nc.6.freq INO_up800.s
    ...

接下来就是每个被发现模体的详细信息：

.. code:: python

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

使用下面的方法来读取这个文件（以 ``meme.dna.oops.txt`` 存储）：

.. code:: python

    >>> handle = open("meme.dna.oops.txt")
    >>> record = motifs.parse(handle, "meme")
    >>> handle.close()

``motifs.parse`` 命令直接读取整个文件，所以在使用后可以关闭这个文件。其中头文件信息被存储于属性中

.. code:: python

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

这个数据记录是 ``Bio.motifs.meme.Record`` 类的一个对象。这个类继承于列表（list），所以你可以把这个 ``record`` 看成模体对象的一个列表：

.. code:: python

    >>> len(record)
    2
    >>> motif = record[0]
    >>> print motif.consensus
    TTCACATGCCGC
    >>> print motif.degenerate_consensus
    TTCACATGSCNC

除了一般的模体属性外，每个模体还同时保存着它们由MEME计算的各自特异信息。例如：

.. code:: python

    >>> motif.num_occurrences
    7
    >>> motif.length
    12
    >>> evalue = motif.evalue
    >>> print "%3.1g" % evalue
    0.2
    >>> motif.name
    'Motif 1'

除了像上面所做的用索引来获得相关记录，你也可以用它的名称来找到这个记录：

.. code:: python

    >>> motif = record['Motif 1']

每个模体都有一个 ``.instances`` 属性与在这个被发现模体中的序列实例，能够为每个实例提供一些信息：

.. code:: python

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

.. code:: python

    >>> print "%5.3g" % pvalue
    1.85e-08

MAST
^^^^

TRANSFAC
^^^^^^^^

TRANSFAC是一个为转录因子手动创建的一个专业数据库，同时还包括染色体结合位点和DNA结合的描述 [:ref:`27 <matys2003>`] 。TRANSFAC数据库中所用的文件格式至今还被其他工具所使用，我们下面将介绍TRANSFAC文件格式。

TRANSFAC文件格式简单概括如下：

.. code:: python

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

这个文件显示了模体 ``motif1`` 中12个核苷酸的频率矩阵。总的来说，一个TRANSFAC文件里面可以包含多个模体。以下是示例文件 ``transfac.dat`` 的内容：

.. code:: python

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

可用如下方法读取TRANSFAC文件：

.. code:: python

    >>> handle = open("transfac.dat")
    >>> record = motifs.parse(handle, "TRANSFAC")
    >>> handle.close()

如果有总版本号的话，它是存储在 ``record.version`` 中：

.. code:: python

    >>> record.version
    'EXAMPLE January 15, 2013'

每个在 ``record`` 中的模体都是 ``Bio.motifs.transfac.Motif`` 类的实例，这些实例同时继承 ``Bio.motifs.Motif`` 类和Python字典的属性。这些字典用双字母的键来存储关于这个模体的其他附加信息：

.. code:: python

    >>> motif = record[0]
    >>> motif.degenerate_consensus # Using the Bio.motifs.Motif method
    Seq('SRACAGGTGKYG', IUPACAmbiguousDNA())
    >>> motif['ID'] # Using motif as a dictionary
    'motif1'

TRANSFAC文件一般比这些例子更详细，包含了许多关于模体的附加信息。表格 :ref:`14.1.2 <table-transfaccodes>` 列出了在TRANSFAC文件常见的双字母含义：

--------------

.. _table-transfaccodes:

+-------------------------------------------------------+
| Table 14.1: TRANSFAC文件中常见的字段                  |
+-------------------------------------------------------+

+----------+---------------------------------------------------+
| ``AC``   | Accession numbers 序列号                          |
+----------+---------------------------------------------------+
| ``AS``   | Accession numbers, secondary 第二序列号           |
+----------+---------------------------------------------------+
| ``BA``   | Statistical basis 统计依据                        |
+----------+---------------------------------------------------+
| ``BF``   | Binding factors 结合因子                          |
+----------+---------------------------------------------------+
| ``BS``   | Factor binding sites underlying the matrix        |
|          | 基于矩阵的转录结合位点                            | 
+----------+---------------------------------------------------+
| ``CC``   | Comments 注解                                     |
+----------+---------------------------------------------------+
| ``CO``   | Copyright notice 版权事项                         |
+----------+---------------------------------------------------+
| ``DE``   | Short factor description 短因子说明               |
+----------+---------------------------------------------------+
| ``DR``   | External databases 外部数据库                     |
+----------+---------------------------------------------------+
| ``DT``   | Date created/updated 创建或更新日期               |
+----------+---------------------------------------------------+
| ``HC``   | Subfamilies 亚家庭名称                            |
+----------+---------------------------------------------------+
| ``HP``   | Superfamilies 超家庭名称                          |
+----------+---------------------------------------------------+
| ``ID``   | Identifier 身份证                                 |
+----------+---------------------------------------------------+
| ``NA``   | Name of the binding factor 结合因子的名称         |
+----------+---------------------------------------------------+
| ``OC``   | Taxonomic classification 分类                     |
+----------+---------------------------------------------------+
| ``OS``   | Species/Taxon 种类或分类                          |
+----------+---------------------------------------------------+
| ``OV``   | Older version 旧版本                              |
+----------+---------------------------------------------------+
| ``PV``   | Preferred version 首选版本                        |
+----------+---------------------------------------------------+
| ``TY``   | Type 类型                                         |
+----------+---------------------------------------------------+
| ``XX``   | Empty line; these are not stored in the Record.   |
|          | 空白行;没在记录中存储的数据                       | 
+----------+---------------------------------------------------+

--------------

每个模体同时也有一个包含与这个模体相关参考资料的 ``references`` 属性，用下面的双字母键来获得：

--------------

+-----------------------------------------------------------------+
| Table 14.2: TRANSFAC文件中用来存储参考资料的字段                |
+-----------------------------------------------------------------+

+----------+-------------------------------+
| ``RN``   | Reference number 参考数目     |
+----------+-------------------------------+
| ``RA``   | Reference authors 参考资料作者|
+----------+-------------------------------+
| ``RL``   | Reference data 参考数据       |
+----------+-------------------------------+
| ``RT``   | Reference title 参考标题      |
+----------+-------------------------------+
| ``RX``   | PubMed ID                     |
+----------+-------------------------------+

--------------

将TRANSFAC文件按原来格式打印出来：

.. code:: python

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

通过用字符串形式来截取输出并且保存在文件中，你可以按TRANSFAC的格式导出这些模体：

.. code:: python

    >>> text = str(record)
    >>> handle = open("mytransfacfile.dat", 'w')
    >>> handle.write(text)
    >>> handle.close()

14.1.3  模体写出
~~~~~~~~~~~~~~~~~~~~~~

说到导出，我们可以先看看导出函数。以JASPAR ``.pfm`` 格式导出模体文件，可以用：

.. code:: python

    >>> print m.format("pfm")
    3       7       0       2       1
    0       0       5       2       6
    0       0       0       3       0
    4       0       2       0       0
    <BLANKLINE>

用类似TRANSFAC的格式导出一个模体：

.. code:: python

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

你可以用 ``motifs.write`` 来写出多个模体。这个函数在使用的时候不必担心这些模体来自于TRANSFAC文件。比如：

.. code:: python

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

14.1.4  绘制序列标识图
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

如果能够联网，我们可以创建一个 `weblogo <http://weblogo.berkeley.edu>`__ ：

.. code:: python

    >>> arnt.weblogo("Arnt.png")

将得到的标识图存储成PNG格式。

14.2  位置权重矩阵
------------------------------

模体对象的 ``.counts`` 属性能够显示在序列上每个位置核苷酸出现的次数。我们可以把这矩阵除以序列中的实例数目来标准化这矩阵，得到每个核苷酸在序列位置上出现概率。我们把这概率看作位置权重矩阵。不过，要知道在字面上，这个术语也可以用来说明位置特异性得分矩阵，这个我们将会在下面讨论。

通常来说，伪计数（pseudocounts）在归一化之前都已经加到每个位置中。这样可以避免在这序列上过度拟合位置权重矩阵以至趋向于模体的实例的有限数量，还可以避免概率为0。向每个位置的核苷酸添加一个固定的伪计数，可以为 ``pseudocounts`` 参数指定一个数值：

.. code:: python

    >>> pwm = m.counts.normalize(pseudocounts=0.5)
    >>> print pwm
            0      1      2      3      4
    A:   0.39   0.83   0.06   0.28   0.17
    C:   0.06   0.06   0.61   0.28   0.72
    G:   0.06   0.06   0.06   0.39   0.06
    T:   0.50   0.06   0.28   0.06   0.06
    <BLANKLINE>

另外， ``pseudocounts`` 可以利用字典为每个核苷酸指定一个伪计数值。例如，由于在人类基因组中GC含量大概为40%,因此可以选择下面这些伪计数值：

.. code:: python

    >>> pwm = m.counts.normalize(pseudocounts={'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6})
    >>> print pwm
            0      1      2      3      4
    A:   0.40   0.84   0.07   0.29   0.18
    C:   0.04   0.04   0.60   0.27   0.71
    G:   0.04   0.04   0.04   0.38   0.04
    T:   0.51   0.07   0.29   0.07   0.07
    <BLANKLINE>

位置权重矩阵有它自己的方法计算一致序列、反向一致序列和简并一致序列：

.. code:: python

    >>> pwm.consensus
    Seq('TACGC', IUPACUnambiguousDNA())
    >>> pwm.anticonsensus
    Seq('GGGTG', IUPACUnambiguousDNA())
    >>> pwm.degenerate_consensus
    Seq('WACNC', IUPACAmbiguousDNA())

应当注意到由于伪计数的原因，由位置仅重矩阵计算得到的简并一致序列和由模体中实例计算得到的简并一致序列有一点不同：

.. code:: python

    >>> m.degenerate_consensus
    Seq('WACVC', IUPACAmbiguousDNA())

位置权重矩阵的反向互补矩阵可以直接用 ``pwm`` 计算出来：

.. code:: python

    >>> rpwm = pwm.reverse_complement()
    >>> print rpwm
            0      1      2      3      4
    A:   0.07   0.07   0.29   0.07   0.51
    C:   0.04   0.38   0.04   0.04   0.04
    G:   0.71   0.27   0.60   0.04   0.04
    T:   0.18   0.29   0.07   0.84   0.40
    <BLANKLINE>

14.3  位置特异性得分矩阵
----------------------------------------

使用背景分布和加入伪计数的PWM，很容易就能计算出log-odds比率，提供特定标记的log odds值，这值来自于在这个背景的模体。我们可以用在位置仅重矩阵中 ``.log-odds()`` 方法：

.. code:: python

    >>> pssm = pwm.log_odds()
    >>> print pssm
            0      1      2      3      4
    A:   0.68   1.76  -1.91   0.21  -0.49
    C:  -2.49  -2.49   1.26   0.09   1.51
    G:  -2.49  -2.49  -2.49   0.60  -2.49
    T:   1.03  -1.91   0.21  -1.91  -1.91
    <BLANKLINE>

这时我们可以更经常看到特定标记和背景下的正值和负值。0.0意味着在模体和背景中观察到一个标记有相等的可能性。

上面是假设A,C,G和T在背景中出现的概率是相同的。那在A,C,G和T出现概率不同的情况下，为了计算特定背景下的位置特异性得分矩阵，可以使用 ``background`` 参数。例如，在40%GC含量的背景下，可以用：

.. code:: python

    >>> background = {'A':0.3,'C':0.2,'G':0.2,'T':0.3}
    >>> pssm = pwm.log_odds(background)
    >>> print pssm
            0      1      2      3      4
    A:   0.42   1.49  -2.17  -0.05  -0.75
    C:  -2.17  -2.17   1.58   0.42   1.83
    G:  -2.17  -2.17  -2.17   0.92  -2.17
    T:   0.77  -2.17  -0.05  -2.17  -2.17
    <BLANKLINE>

从PSSM中得到的最大和最小值被存储在 ``.max`` 和 ``.min`` 属性中：

.. code:: python

    >>> print "%4.2f" % pssm.max
    6.59
    >>> print "%4.2f" % pssm.min
    -10.85

在特定背景下计算平均值和标准方差使用 ``.mean`` 和 ``.std`` 方法。

.. code:: python

    >>> mean = pssm.mean(background)
    >>> std = pssm.std(background)
    >>> print "mean = %0.2f, standard deviation = %0.2f" % (mean, std)
    mean = 3.21, standard deviation = 2.59

如果没有指定特定的背景，就会使用一个统一的背景。因为同KL散度或相对熵的值相同，所以平均值就显得特别重要，并且它也是同背景相比的模体信息含量的测量方法。由于在Biopython中用以2为底的对数来计算log-odds值，信息含量的的单位是bit。

``.reverse_complement``, ``.consensus``, ``.anticonsensus`` 和 ``.degenerate_consensus`` 方法可以直接对PSSM使用。

14.4  搜索实例
-----------------------------

模体最常用的功能就是在序列中的查找它的实例。在这节，我们会用如下的序列作为例子：

.. code:: python

    >>> test_seq=Seq("TACACTGCATTACAACCCAAGCATTA",m.alphabet)
    >>> len(test_seq)
    26

14.4.1  搜索准确匹配实例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

查找实例最简单的方法就是查找模体实例的准确匹配：

.. code:: python

    >>> for pos,seq in m.instances.search(test_seq):
    ...     print pos, seq
    ... 
    0 TACAC
    10 TACAA
    13 AACCC

我们可获得反向互补序列（找到互补链的实例）：

.. code:: python

    >>> for pos,seq in r.instances.search(test_seq):
    ...     print pos, seq
    ... 
    6 GCATT
    20 GCATT

14.4.2  用PSSM得分搜索匹配实例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

在模体中很容易找出相应的位置,引起对模体的高log-odds值：

.. code:: python

    >>> for position, score in pssm.search(test_seq, threshold=3.0):
    ...     print "Position %d: score = %5.3f" % (position, score)
    ... 
    Position 0: score = 5.622
    Position -20: score = 4.601
    Position 10: score = 3.037
    Position 13: score = 5.738
    Position -6: score = 4.601

负值的位置是指在测试序列的反向链中找到的模体的实例，而且得力于Python的索引。在 ``pos`` 的模体实例可以用 ``test_seq[pos:pos+len(m)]`` 来定位，不管 ``pos`` 值是正还是负。

你可能注意到阀值参数，在这里随意地设为3.0。这里是 *log*\ :sub:`2` ，所以我们现在开始寻找那些在模体中出现概率为背景中出现概率8倍序列。默认的阀值是0.0,在此阀值下，会把所有比背景中出现概率大的模体实例都找出来。

.. code:: python

    >>> pssm.calculate(test_seq)
    array([  5.62230396,  -5.6796999 ,  -3.43177247,   0.93827754,
            -6.84962511,  -2.04066086, -10.84962463,  -3.65614533,
            -0.03370807,  -3.91102552,   3.03734159,  -2.14918518,
            -0.6016975 ,   5.7381525 ,  -0.50977498,  -3.56422281,
            -8.73414803,  -0.09919716,  -0.6016975 ,  -2.39429784,
           -10.84962463,  -3.65614533], dtype=float32)

通常来说，上述是计算PSSM得分的最快方法。这些得分只能由前导链用 ``pssm.calculate`` 计算得到。为了得到互补链的PSSM值，你可以利用PSSM的互补矩阵：

.. code:: python

    >>> rpssm = pssm.reverse_complement()
    >>> rpssm.calculate(test_seq)
    array([ -9.43458748,  -3.06172252,  -7.18665981,  -7.76216221,
            -2.04066086,  -4.26466274,   4.60124254,  -4.2480607 ,
            -8.73414803,  -2.26503372,  -6.49598789,  -5.64668512,
            -8.73414803, -10.84962463,  -4.82356262,  -4.82356262,
            -5.64668512,  -8.73414803,  -4.15613794,  -5.6796999 ,
             4.60124254,  -4.2480607 ], dtype=float32)

14.4.3  选择得分阀值
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

如果不想刚才那么随意设定一个阀值，你可以探究一下PSSM得分的分布。由于得分的空间分布随着模体长度而成倍增长，我们用一个近似于给定精度值来计算，如此可使计算成本更容易控制：

.. code:: python

    >>> distribution = pssm.distribution(background=background, precision=10**4)

``distribution`` 对象可以用来决定许多不同的阀值。我们可以指定一个需要的的假阳性率（找到一个由序列在此背景下产生的模体实例的概率）：

.. code:: python

    >>> threshold = distribution.threshold_fpr(0.01)
    >>> print "%5.3f" % threshold
    4.009

或者假阴性率（找不到模体产生的实例概率）：

.. code:: python

    >>> threshold = distribution.threshold_fnr(0.1)
    >>> print "%5.3f" % threshold
    -0.510

或者一个阀值（近似），满足假阳性率和假阴性率之间的关系（fnr/fpr≃ *t*)：

.. code:: python

    >>> threshold = distribution.threshold_balanced(1000)
    >>> print "%5.3f" % threshold
    6.241

或者一个阀值能够大体满足假阳性率和信息含量的 −\ *log* 值之间的相等关系（与Hertz和Stormo的Patser软件所用的一样）：

.. code:: python

    >>> threshold = distribution.threshold_patser()
    >>> print "%5.3f" % threshold
    0.346

比如在我们这个模体中，当以1000比率的平衡阀值查找实例，你可以得到一个让你获得相同结果的阀值（对这个序列来说）。

.. code:: python

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

14.5  模体对象自身相关的位置特异性得分矩阵
--------------------------------------------------------------------------

为了更好的利用PSSMs来查找潜在的TFBSs，每个模体都同位置权重矩阵和位置特异性得分矩阵相关联。用Arnt模体来举个例子：

.. code:: python

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

.. code:: python

    >>> print motif.pssm
            0      1      2      3      4      5
    A:  -0.32   1.93   -inf   -inf   -inf   -inf
    C:   1.68   -inf   2.00   -inf   -inf   -inf
    G:   -inf  -2.32   -inf   2.00   -inf   2.00
    T:   -inf   -inf   -inf   -inf   2.00   -inf
    <BLANKLINE>

在这出现的负无穷大是由于在频率矩阵中相关项的值为0,并且我们默认使用0作为伪计数：

.. code:: python

    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.pseudocounts[letter])
    ...
    A: 0.00
    C: 0.00
    G: 0.00
    T: 0.00

如果你更改了 ``.pseudocouts`` 属性，那么位置频率矩阵和位置特异性得分矩阵就都会自动重新计算：

.. code:: python

    >>> motif.pseudocounts = 3.0
    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.pseudocounts[letter])
    ...
    A: 3.00
    C: 3.00
    G: 3.00
    T: 3.00

.. code:: python

    >>> print motif.pwm
            0      1      2      3      4      5
    A:   0.22   0.69   0.09   0.09   0.09   0.09
    C:   0.59   0.09   0.72   0.09   0.09   0.09
    G:   0.09   0.12   0.09   0.72   0.09   0.72
    T:   0.09   0.09   0.09   0.09   0.72   0.09
    <BLANKLINE>

.. code:: python

    >>> print motif.pssm
            0      1      2      3      4      5
    A:  -0.19   1.46  -1.42  -1.42  -1.42  -1.42
    C:   1.25  -1.42   1.52  -1.42  -1.42  -1.42
    G:  -1.42  -1.00  -1.42   1.52  -1.42   1.52
    T:  -1.42  -1.42  -1.42  -1.42   1.52  -1.42
    <BLANKLINE>

如果你想对4个核苷酸使用不同的伪计数，可以使用字典来设定4个核苷酸的 ``pseudocounts`` 。把 ``motif.pseudocounts`` 设为 ``None`` 会让伪计数重置为0的默认值。

位置特异性得分矩阵依赖于一个默认均一的背景分布：

.. code:: python

    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.background[letter])
    ...
    A: 0.25
    C: 0.25
    G: 0.25
    T: 0.25

同样，如果你更改了背景分布，位置特异性得分矩阵也会重新计算：

.. code:: python

    >>> motif.background = {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2}
    >>> print motif.pssm
            0      1      2      3      4      5
    A:   0.13   1.78  -1.09  -1.09  -1.09  -1.09
    C:   0.98  -1.68   1.26  -1.68  -1.68  -1.68
    G:  -1.68  -1.26  -1.68   1.26  -1.68   1.26
    T:  -1.09  -1.09  -1.09  -1.09   1.85  -1.09
    <BLANKLINE>

把 ``motif.backgroud`` 设为 ``None`` 后会将其重置为均一的分布。

.. code:: python

    >>> motif.background = None
    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.background[letter])
    ...
    A: 0.25
    C: 0.25
    G: 0.25
    T: 0.25

如果你把 ``motif.background`` 设为一个单一值，这个值将会被看成是GC含量：

.. code:: python

    >>> motif.background = 0.8
    >>> for letter in "ACGT":
    ...     print "%s: %4.2f" % (letter, motif.background[letter])
    ...
    A: 0.10
    C: 0.40
    G: 0.40
    T: 0.10

应当注意到你能够在当前计算背景下计算PSSM的平均值：

.. code:: python

    >>> print "%f" % motif.pssm.mean(motif.background)
    4.703928

它的标准方差也是一样：

.. code:: python

    >>> print "%f" % motif.pssm.std(motif.background)
    3.290900

和它的分布：

.. code:: python

    >>> distribution = motif.pssm.distribution(background=motif.background)
    >>> threshold = distribution.threshold_fpr(0.01)
    >>> print "%f" % threshold
    3.854375

请注意，每当你调用 ``motif.pwm`` 或 ``motif.pssm`` ，位置仅重矩阵和位置特异性得分矩阵都会重新计算。如果看重速度并且需要重复用到PWM或PSSM时，你可以把他们保存成变量，如下所示：

.. code:: python

    >>> pssm = motif.pssm

14.6  模体比较
----------------------

当有多个模体时，我们就会想去比较它们。

在我们开始比较之前，应当要指出模体的边界通常比较模糊。这也就是说我们需要比较不同长度的模体，因此这些比较也要涉及到相关的比对。所以我们需要考虑两个东西：

-   模体比对
-   比较比对后模体的相关函数

为了比对模体，我们使用PSSMs的不含间隔的比对，并且用0来代替矩阵开始和结束位置缺失的列。这说明我们能够有效地利用背景分布来代替PSSM中缺失的列。距离函数然后可以返回模体间最小的距离，以及比对中相应的偏移量。

举个例子，先导入和测试模体 ``m`` 相似的模体：

.. code:: python

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

为了让模体能够进行相互比较，选择和模体 ``m`` 相同伪计数和背景值：

.. code:: python

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

我们将用皮尔逊相关（Pearson correlation）来比较这些模体。由于我们想要让它偏向于一个距离长度，我们实际上取1−\ *r* ，其中 *r* 是皮尔逊相关系数（Pearson correlation coefficient，PCC）：

.. code:: python

    >>> distance, offset = pssm.dist_pearson(pssm_reb1)
    >>> print "distance = %5.3g" % distance
    distance = 0.239
    >>> print offset
    -2

这意味着模体 ``m`` 和 ``m_reb1`` 间最佳PCC可以从下面的比对中获得：

.. code:: python

    m:      bbTACGCbb
    m_reb1: GTTACCCGG

其中 ``b`` 代表背景分布。PCC值大概为1−0.239=0.761。


14.7  查找 *De novo* 模体
-----------------------------

如今，Biopython对 *De novo* 模体查找的支持是有限的。也就是说，我们支持AlignAce和MEME的运行和读取。由于模体查找工具如雨后春笋般出现，所以很欢迎新的分析程序加入进来。

14.7.1  MEME
~~~~~~~~~~~~

假设用MEME以你喜欢的参数设置来跑序列，并把结果保存在文件 ``meme.out`` 中。你可以用以下的命令来得到MEME输出的模体：

.. code:: python

    >>> from Bio import motifs
    >>> motifsM = motifs.parse(open("meme.out"), "meme")

.. code:: python

    >>> motifsM
    [<Bio.motifs.meme.Motif object at 0xc356b0>]

除了最想要的一系列模体外，结果中还包含了很多有用的信息，可以通过那些一目了然的属性名获得：

-  ``.alphabet``
-  ``.datafile``
-  ``.sequence_names``
-  ``.version``
-  ``.command``

由MEME解析得到的模体可以像平常的模体对象（有实例）一样处理，它们也提供了一些额外的功能，可以为实例增加额外的信息。

.. code:: python

    >>> motifsM[0].consensus
    Seq('CTCAATCGTA', IUPACUnambiguousDNA())
    >>> motifsM[0].instances[0].sequence_name
    'SEQ10;'
    >>> motifsM[0].instances[0].start
    3
    >>> motifsM[0].instances[0].strand
    '+'

.. code:: python

    >>> motifsM[0].instances[0].pvalue
    8.71e-07

14.7.2  AlignAce
~~~~~~~~~~~~~~~~

我们可以用AlignACE程序实现类似的效果。假如，你把结果保存在 ``alignace.out`` 文件中。你可以用下面的代码读取结果：

.. code:: python

    >>> from Bio import motifs
    >>> motifsA = motifs.parse(open("alignace.out"),"alignace")

同样，你的模体也和正常的模体对象有相同的属性：

.. code:: python

    >>> motifsA[0].consensus
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

事实上，你甚至可以观察到，AlignAce找到了一个和MEME非常相似的模体。下面只是MEME模体互补链的一个较长版本：

.. code:: python

    >>> motifsM[0].reverse_complement().consensus
    Seq('TACGATTGAG', IUPACUnambiguousDNA())

如果你的机器上安装了AlignAce，你可以直接从Biopython中运行AlignAce。下面就是一个如何运行AlignAce的简单例子（其他参数可以用关键字参数来调用）：

.. code:: python

    >>> command="/opt/bin/AlignACE"
    >>> input_file="test.fa"
    >>> from Bio.motifs.applications import AlignAceCommandline
    >>> cmd = AlignAceCommandline(cmd=command,input=input_file,gcback=0.6,numcols=10)
    >>> stdout,stderr= cmd()

由于AlignAce把所有的结果输出到标准输出，所以你可以通过读取结果的第一部分来获得模体：

.. code:: python

    >>> motifs = motifs.parse(stdout,"alignace")

.. _sec-links:

14.8  相关链接
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

14.9  旧版Bio.Motif模块
-------------------------------

本章剩下部分将介绍Biopython 1.61版本前的 ``Bio.Motifs`` 模块，该模块取代了Biopython 1.50版本中基于两个早期Biopython模块—— ``Bio.AlignAce`` 和 ``Bio.MEME`` 的 ``Bio.Motif`` 模块。

为了平滑的过渡，早期的 ``Bio.Motif`` 模块将会和它的取代者 ``Bio.Motifs`` 一同维护到至少发行两个版本，并且持续至少一年。

14.9.1  模体对象
~~~~~~~~~~~~~~~~~~~~~

由于我们对模体分析感兴趣，不过让我们首先看看 ``Motif`` 对象。第一步要先导入模体库：

.. code:: python

    >>> from Bio import Motif

然后可以开始创建第一个模体对象。创建一个DNA模体：

.. code:: python

    >>> from Bio.Alphabet import IUPAC
    >>> m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)

现在这里面什么也没有，往新建的模体加入一些序列：

.. code:: python

    >>> from Bio.Seq import Seq
    >>> m.add_instance(Seq("TATAA",m.alphabet))
    >>> m.add_instance(Seq("TATTA",m.alphabet))
    >>> m.add_instance(Seq("TATAA",m.alphabet))
    >>> m.add_instance(Seq("TATAA",m.alphabet))

现在我们有了一个完整的 ``Motif`` 实例，我们可以试着从中获取一些基本信息。先看看长度和一致序列：

.. code:: python

    >>> len(m)
    5
    >>> m.consensus()
    Seq('TATAA', IUPACUnambiguousDNA())

对于DNA模体，我们还可以获得一个模体的反向互补序列：

.. code:: python

    >>> m.reverse_complement().consensus()
    Seq('TTATA', IUPACUnambiguousDNA())
    >>> for i in m.reverse_complement().instances:
    ...     print i
    TTATA
    TAATA
    TTATA
    TTATA

我们也可以简单的调取模体的信息容量：

.. code:: python

    >>> print "%0.2f" % m.ic()
    5.27

这给我们提供了模体中信息容量的比特数，这指出和背景有多少不同。

展示模体最常用的就是PWM（位置仅重矩阵）。它概括了在模体上任意位置出现一个符号（这里指核苷酸）的概率。这个可以用 ``.pwm()`` 方法来计算：

.. code:: python

    >>> m.pwm()
    [{'A': 0.05, 'C': 0.05, 'T': 0.85, 'G': 0.05}, 
     {'A': 0.85, 'C': 0.05, 'T': 0.05, 'G': 0.05}, 
     {'A': 0.05, 'C': 0.05, 'T': 0.85, 'G': 0.05}, 
     {'A': 0.65, 'C': 0.05, 'T': 0.25, 'G': 0.05}, 
     {'A': 0.85, 'C': 0.05, 'T': 0.05, 'G': 0.05}]

模体的PWM中的概率是基于实例中的计数，但我们发现，虽然模体中没有出现G和C，可是它们的概率仍然是非0的。这主要是因为有伪计数的存在，简单地说，就是一种常用的方式来承认我们认知的不完备以及为了避免使用0进行对数运算而出现的技术问题。

我可以调整伪计数添加到模体对象两个属性的方式。 ``.background`` 是我们假设代表背景分布的所有字符的概率分布，是非模体序列（通常基于各自基因组的GC含量）。在模体创建的时候，就默认的设置为一个统一分布：

.. code:: python

    >>> m.background  
    {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25}

另一个就是 ``.beta`` ，这个参数可以说明我们应该给伪计数设定为何值。默认设定为1.0。

.. code:: python

    >>> m.beta
    1.0

所以输入伪计数的总量等于一个实例的输入总量。

使用背景分布和附加了伪计数的pwm，可以很容易的计算log-odd比率，这告诉我们在背景下，一个来自模体特定碱基的log-odd值。我们可以使用 ``.log_odds()`` 方法：

.. code:: python

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

此处，我们可以看出如果模体中的碱基比背景中出现频率更高，其值为正值，反之则为负值。0.0说明在背景和模体中出现的概率是相同的（如第二个位置的“T”）。

14.9.1.1  模体读写
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

手动从实例创建一个模体确实没什么技术含量，所以很有必要有一些读写功能来读取和写出模体。对于如何存储模体还没有一个固定的标准，但是有一些格式比其他格式更流行。这些格式的主要区别在于模体的创建是基于实例还是一些PWM矩阵。其中一个最流行的模体数据库就是 `JASPAR <http://jaspar.genereg.net>`__ ，该数据库保存了上述两种类型的格式，所以让我们看看是如何从实例中导入JASPAR模体：

.. code:: python

    >>> from Bio import Motif
    >>> arnt = Motif.read(open("Arnt.sites"),"jaspar-sites")

从一个计数矩阵中导入：

.. code:: python

    >>> srf = Motif.read(open("SRF.pfm"),"jaspar-pfm")

``arnt`` 和 ``srf`` 模体可以为我们做相同的事情，但是它们使用不同的内部表现形式来展现模体。我们可以用 ``has_counts`` 和 ``has_instances`` 属性来区分它们：

.. code:: python

    >>> arnt.has_instances
    True
    >>> srf.has_instances
    False
    >>> srf.has_counts
    True

.. code:: python

    >>> srf.counts
    {'A': [2, 9, 0, 1, 32, 3, 46, 1, 43, 15, 2, 2],
     'C': [1, 33, 45, 45, 1, 1, 0, 0, 0, 1, 0, 1],
     'G': [39, 2, 1, 0, 0, 0, 0, 0, 0, 0, 44, 43],
     'T': [4, 2, 0, 0, 13, 42, 0, 45, 3, 30, 0, 0]}

对于模体的不同表现形式，可以用转换功能来实现相互转换：

.. code:: python

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

在这里需要注意的是 ``make_instances_from_counts()`` 方法创建的是假实例，因为按照相同的pwm能够得到许多不同的实例，所以不能反过来重建原来的矩阵。不过这对我们利用PWM来展现模体没有什么影响，但从基于计数的模体中导出实例时要小心。

说到导出，让我们看看导出函数。我们可以按fasta的格式导出：

.. code:: python

    >>> print m.format("fasta")
    >instance0
    TATAA
    >instance1
    TATTA
    >instance2
    TATAA
    >instance3
    TATAA

或者是按TRANSFAC样的矩阵格式导出（能被一些处理软件识别）

.. code:: python

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

最后，如果能够联网，我们可以创建一个 `weblogo <http://weblogo.berkeley.edu>`__ ：

.. code:: python

    >>> arnt.weblogo("Arnt.png")

我们可以把得到的标识图以png的格式保存到特定的文件中。

14.9.2  查找实例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

模体中最常用的就是在一些序列中查找实例。为解释这部分，我们将手动创建一个如下的序列：

.. code:: python

    test_seq=Seq("TATGATGTAGTATAATATAATTATAA",m.alphabet)

查找实例最简单的方法就是在模体中查找具体匹配的实例：

.. code:: python

    >>> for pos,seq in m.search_instances(test_seq):
    ...     print pos,seq.tostring()
    ... 
    10 TATAA
    15 TATAA
    21 TATAA

对于互补序列，也可以用相同的方法（为了找到互补链上的实例）：

.. code:: python

    >>> for pos,seq in m.reverse_complement().search_instances(test_seq):
    ...     print pos,seq.tostring()
    ... 
    12 TAATA
    20 TTATA

提高模体的log-odds值能让查为位置更加简单:


.. code:: python

    >>> for pos,score in m.search_pwm(test_seq,threshold=5.0):
    ...     print pos,score
    ... 
    10 8.44065060871
    -12 7.06213898545
    15 8.44065060871
    -20 8.44065060871
    21 8.44065060871

你可能注意到阀值参数，在这里随意地设为5.0。按 *log*\ :sub:`2` 来算，我们应当查找那些在模体中出现概率为背景中出现概率32倍的序列。默认的阀值是0.0,在些阀值下，会把所有比背景中出现概率大的模体实例都找出来。

如果不想那么随意的选择一个阀值，你可以研究一下 ``Motif.score_distribution`` 类，它为模体提供一个相应的得分分布。由于得分的空间分布随着模体长度而成倍增长，我们正用一个近似于给定精度值计算，从而使计算成本易于控制：

.. code:: python

    >>> sd = Motif.score_distribution(m,precision=10**4)

上面那个sd对象可以用来决定许多不同的阀值。

我们可以设定一个需要的假阳性率（找到一个由此序列在这个背景下产生的模体实例的概率）：

.. code:: python

    >>> sd.threshold_fpr(0.01)
    4.3535838726139886

或者假阴性率（找不到模体产生的实例的概率）：

.. code:: python

    >>> sd.threshold_fnr(0.1)
    0.26651713652234044

或者一个阀值（近似），满足假阳性率和假阴性率之间的关系（fnr/fpr≃ *t*)：

.. code:: python

    >>> sd.threshold_balanced(1000)
    8.4406506087056368

或者一个阀值能够大体满足假阳性率和信息含量的 −\ *log* 值之间的相等关系（像Hertz和Stormo的Patser软件所用的一样）：

在我们这个例子中，当以1000比率的平衡阀值查找实例时，你可以得到一个让你获得相同结果（对于这个序列来说）的阀值：

.. code:: python

    >>> for pos,score in m.search_pwm(test_seq,threshold=sd.threshold_balanced(1000)):
    ...     print pos,score
    ... 
    10 8.44065060871
    15 8.44065060871
    -20 8.44065060871
    21 8.44065060871

14.9.3  模体比较
~~~~~~~~~~~~~~~~~~~~~~~~

当有多个模体时，我们就会想去比较他们。对此， ``Bio.Motif`` 有三种不同的方法来进行模体比较。

在我们开始比较之前，应当指出模体的边界通常是相当模糊的。也就是说我们经常需要比较不同长度的模体，因此这些比较涉及到相关的比对。所以我们需要考虑两个要点：

-   模体比对
-   比较比对后模体的相关函数

在 ``Bio.Motif`` 中有三种比较方法，这些方法都是基于来源于模体比对的想法，而采用不同方式。简单来说，我们使用不含间隔的PSSMs比对，并且用0来代替矩阵同背景相比，在开始和结束位置出现缺失的列。这三种比较方法都可以解释成距离估量，但是只有一个（ ``dist——dpq`` ）满足三角不等式。这些方法都返回距离的最小值和模体相应的偏移量。

为了展示这些比较功能是如何实现的，导入和测试模体 ``m`` 相似的其他模体：

.. code:: python

    >>> ubx=Motif.read(open("Ubx.pfm"),"jaspar-pfm")
    <Bio.Motif.Motif.Motif object at 0xc29b90>
    >>> ubx.consensus()
    Seq('TAAT', IUPACUnambiguousDNA())

第一个展示的功能是基于皮尔逊相关（Pearson correlation）的。因为我们想让它类似于一个距离估量，所以我们实际上取 1−\ *r* ，其中的 *r* 是皮尔逊相关系数（Pearson correlation coefficient，PCC）：

.. code:: python

    >>> m.dist_pearson(ubx)
    (0.41740393308237722, 2)

这意味着模体 ``m`` 和 ``Ubx`` 间最佳的PCC可以从下面的比对中获得：

.. code:: python

    bbTAAT
    TATAAb

其中 ``b`` 代表背景分布。PCC值大概为 1-0.42=0.58.如果我们尝试计算Ubx模体的互补序列：

.. code:: python

    >>> m.dist_pearson(ubx.reverse_complement())
    (0.25784180151584823, 1)

我们可以发现更好的PCC值（大概为0.75），并且比对也是不同的：

.. code:: python

    bATTA
    TATAA

还有两个其他的功能函数： ``dist_dpq`` ,这是基于Kullback-Leibler散度的真正度量（满足三角不等式）。

.. code:: python

    >>> m.dist_dpq(ubx.reverse_complement())
    (0.49292358382899853, 1)

还有 ``dist_product`` 方法，它是基于概率的方法，这概率可以看成是两个模体独立产生两个相同实例的概率。

.. code:: python

    >>> m.dist_product(ubx.reverse_complement())
    (0.16224587301064275, 1)

14.9.4  *De novo* 模体查找
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

目前，Biopython对 *de novo* 模体查找只有一些有限的支持。也就是说，我们只支持AlignAce和MEME的运行和读取。由于现模体查找工具发展如雨后春笋般，我们很欢迎有新的贡献者加入。

14.9.4.1  MEME
^^^^^^^^^^^^^^

假如你以中意的参数用MEME来跑你自己的序列，并把得到的结果保存在 ``meme.out`` 文件中。你可以用以下代码读取MEME产生的文件获得那些模体：

.. code:: python

    >>> motifsM = list(Motif.parse(open("meme.out"),"MEME"))
    >>> motifsM
    [<Bio.Motif.MEMEMotif.MEMEMotif object at 0xc356b0>]

除了那一系列想要的模体外，结果对象中还有很多有用的信息，可以用那些一目了然的属性名来获取：

-  ``.alphabet``
-  ``.datafile``
-  ``.sequence_names``
-  ``.version``
-  ``.command``

MEME解析器得到的模体可以像通常模体（含有实例）一样进行处理，它们也可以通过对实例添加附加信息而提供一些额外的功能。

.. code:: python

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

对于AlignACE程序也可以做相同的事情。假如你把结果存储于文件 ``alignace.out`` 文件中。你可以用以下代码读取结果：

.. code:: python

    >>> motifsA=list(Motif.parse(open("alignace.out"),"AlignAce"))

同样，得到的模体也和平常的模体一样：

.. code:: python

    >>> motifsA[0].consensus()
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

事实上，你甚至可以发现AlignAce和MEME得到的模体十分相似，只不过AlignAce模体是MEME模体反向互补序列的加长版本而已：

.. code:: python

    >>> motifsM[0].reverse_complement().consensus()
    Seq('TACGATTGAG', IUPACUnambiguousDNA())

如果你的机器上安装了AlignAce，你也可以直接从Biopython中启动。下面就是一个如何启动的小例子（其他参数可以用关键字参数指定）：

.. code:: python

    >>> command="/opt/bin/AlignACE"
    >>> input_file="test.fa"
    >>> from Bio.Motif.Applications import AlignAceCommandline
    >>> cmd = AlignAceCommandline(cmd=command,input=input_file,gcback=0.6,numcols=10)
    >>> stdout,stderr= cmd()

由于AlignAce把结果打印到标准输出，因此你可以通过读取结果的第一部分来获得你想要的模体：

.. code:: python

    motifs=list(Motif.parse(stdout,"AlignAce"))


