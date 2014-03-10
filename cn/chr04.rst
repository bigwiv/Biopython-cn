.. _chapter-SeqRecord:

第4章  序列注释对象
======================================

第 :ref:`3 <chapter-Bio.Seq>` 章介绍了序列对象的基本情况。紧接上章的 ``Seq`` 类，这章主要讲Sequence record 或称之为 ``SeqRecord`` 类, 该类在 ``Bio.SeqRecord`` 模块中有定义。 它（见 ``SeqFeature`` 对象）可使序列与高级属性（如identifiers 和 features）关联。其应用贯穿序列输入/输出的交互界面 ``Bio.SeqIO`` 过程中 （详见第 :ref:`5 <chapter-Bio.SeqIO>` 章）。

如读者只需处理FASTA格式的序列文件等简单数据,可略过本章。如涉及带注释内容的数据（如 GenBank或EMBL格式文件）, 本章内容则非常重要。

尽管本章内容涵盖了 ``SeqRecord`` 和 ``SeqFeature`` 对象的大部分内容，但如需了解更多，读者可自行查阅 ``SeqRecord`` wiki (`http://biopython.org/wiki/SeqRecord <http://biopython.org/wiki/SeqRecord>`__ ),和内置帮助文档 (或在线文档 `SeqRecord <http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html>`__ 和 `SeqFeature <http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html>`__ )，获取更多信息:

.. code:: python

    >>> from Bio.SeqRecord import SeqRecord
    >>> help(SeqRecord)
    ...

4.1  SeqRecord对象
-------------------------

``SeqRecord`` (Sequence Record) 类包含在 ``Bio.SeqRecord`` 模块中。该类是 ``Bio.SeqIO`` 序列输入/输出交互界面 (详见第 :ref:`5 <chapter-Bio.SeqIO>` 章)的基本数据类型。可以把identifiers 和features等高级属性与序列关联起来 (参见第 :ref:`3 <chapter-Bio.Seq>` 章)。

``SeqRecord`` 类非常简单,包括下列属性:

**.seq**
    – 序列自身（即 ``Seq`` 对象）。
**.id**
    – 序列主ID（-字符串类型）。通常类同于accession number。
**.name**
    – 序列名/id （-字符串类型）。 可以是accession number, 也可是clone名（类似GenBank record中的LOCUS id）。
**.description**
    – 序列描述（-字符串类型）。
**.letter\_annotations**
    – 对照序列的每个字母逐字注释（per-letter-annotations），以信息名为键（keys），信息内容为值（value）所构成的字典。值与序列等长，用Python列表、元组或字符串表示。.letter\_annotations可用于质量分数(如第 :ref:`18.1.6 <sec-FASTQ-filtering-example>` 节) 或二级结构信息 (如 Stockholm/PFAM 比对文件)等数据的存储。
**.annotations**
    – 用于储存附加信息的字典。信息名为键（keys），信息内容为值（value）。用于保存序列的零散信息（如unstructured information）。
**.features**
    – ``SeqFeature`` 对象列表，储存序列的结构化信息（structured information），如：基因位置, 蛋白结构域。features 详见本章第三节（ 第 :ref:`4.3 <sec-seq_features>` 节）。
**.dbxrefs**
    – 储存数据库交叉引用信息（cross-references）的字符串列表。

4.2  创建 SeqRecord
-------------------------

使用 ``SeqRecord`` 对象非常简单，因为所有的信息都存储在该类的属性中；通常不必手动新建，用 ``Bio.SeqIO`` 从序列文件读取即可（见第 :ref:`5 <chapter-Bio.SeqIO>` 章）。 当然新建 ``SeqRecord`` 也不复杂。

4.2.1  从头新建SeqRecord
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SeqRecord`` 最少只需包含 ``Seq`` 对象:

.. code:: python

    >>> from Bio.Seq import Seq
    >>> simple_seq = Seq("GATC")
    >>> from Bio.SeqRecord import SeqRecord
    >>> simple_seq_r = SeqRecord(simple_seq)

还可以通过初始化函数给 id, name和description赋值；反之，它们被设为默认值“unknown”（可随后编辑）:

.. code:: python

    >>> simple_seq_r.id
    '<unknown id>'
    >>> simple_seq_r.id = "AC12345"
    >>> simple_seq_r.description = "Made up sequence I wish I could write a paper about"
    >>> print simple_seq_r.description
    Made up sequence I wish I could write a paper about
    >>> simple_seq_r.seq
    Seq('GATC', Alphabet())

标识符对输出 ``SeqRecord`` 内容到文件很重要，可随SeqRecord同时建立:

.. code:: python

    >>> from Bio.Seq import Seq
    >>> simple_seq = Seq("GATC")
    >>> from Bio.SeqRecord import SeqRecord
    >>> simple_seq_r = SeqRecord(simple_seq, id="AC12345")

上述章节已提到，``SeqRecord`` 含有一个 ``annotations`` 属性，用于储存各种杂乱注释的字典。添加annotations示例如下:

.. code:: python

    >>> simple_seq_r.annotations["evidence"] = "None. I just made it up."
    >>> print simple_seq_r.annotations
    {'evidence': 'None. I just made it up.'}
    >>> print simple_seq_r.annotations["evidence"]
    None. I just made it up.

``letter_annotations`` 也是字典，其值为与序列等长的内置Python字符串、列表或元组:

.. code:: python

    >>> simple_seq_r.letter_annotations["phred_quality"] = [40,40,38,30]
    >>> print simple_seq_r.letter_annotations
    {'phred_quality': [40, 40, 38, 30]}
    >>> print simple_seq_r.letter_annotations["phred_quality"]
    [40, 40, 38, 30]

``dbxrefs`` 和 ``features`` 分别是字符串和 ``SeqFeature`` 对象的Python列表，将在后续章节讨论。

4.2.2  根据FASTA文件创建SeqRecord对象
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

本节以鼠疫耶尔森菌株（*Yersinia pestis biovar Microtus* str. 91001 ）的pPCP1质粒全长序列为例,说明从FASTA文件创建SeqRecord的过程。该序列原始文件来自NCBI，可在Biopython单元测试GenBank文件夹下找到，也可点击 `NC_005816.fna <http://biopython.org/SRC/biopython/Tests/GenBank/NC_005816.fna>`__ 下载。

序列以大于号开头，该文件只包含一条序列:

.. code:: python

    >gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus ... pPCP1, complete sequence
    TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGGGGGTAATCTGCTCTCC
    ...

回顾第 :ref:`2 <chapter-quick-start>` 章的内容，我们已经遇到过 ``Bio.SeqIO.parse(...)`` 函数，用于遍历 ``SeqRecord`` 对象中的所有记录。 此处，我们介绍 ``Bio.SeqIO`` 模块中的另一个类似函数Bio.SeqIO.read()，用于读取单条序列的文件 （详见第 :ref:`5 <chapter-Bio.SeqIO>` 章）:

.. code:: python

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("NC_005816.fna", "fasta")
    >>> record
    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG',
    SingleLetterAlphabet()), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|',
    description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus ... sequence',
    dbxrefs=[])

现在让我们逐个介绍 ``SeqRecord`` 对象中的主要属性，从给予我们序列属性的 ``Seq`` 对象 开始:

.. code:: python

    >>> record.seq
    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG', SingleLetterAlphabet())

此处 ``Bio.SeqIO`` 默认为通用字母表（generic alphabet）, 而非判断是否DNA序列。如果FASTA文件中序列类型已知，也可通过 ``Bio.SeqIO`` 自行设定 (见第 :ref:`5 <chapter-Bio.SeqIO>` 章用法)。

接下来介绍 identifiers 和 description:

.. code:: python

    >>> record.id
    'gi|45478711|ref|NC_005816.1|'
    >>> record.name
    'gi|45478711|ref|NC_005816.1|'
    >>> record.description
    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus ... pPCP1, complete sequence'

FASTA文件中序列名所在行的第一个单词(去除大于号后) 被当作 ``id`` 和 ``name`` ；而将整行 (去除大于号后) 作为 description。这样设定是为了向后兼容，同时也为了便于处理如下序列:

.. code:: python

    >Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1
    TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGGGGGTAATCTGCTCTCC
    ...

Note: 读取FASTA 文件时其他注释属性为空:

.. code:: python

    >>> record.dbxrefs
    []
    >>> record.annotations
    {}
    >>> record.letter_annotations
    {}
    >>> record.features
    []

本例中FASTA文件源于NCBI，其规范的格式，意味着我们可以方便的解析这些信息并选择提取GI和accession number等信息。然后，对于从其他来源获得的FASTA文件，并不能确保能获得这些信息。

4.2.3  从 GenBank文件创建 SeqRecord
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

仍以疫耶尔森菌株pPCP1质粒全长序列（*Yersinia pestis biovar Microtus* str. 91001 plasmid pPCP1）为例，不同的是这次使用Genbank格式的文件，该文件同样包含在Biopython单元测试/GenBank文件夹下, 也可点击 `NC_005816.gb <http://biopython.org/SRC/biopython/Tests/GenBank/NC_005816.gb>`__
下载。

该文件只含一条记录 (只有一个 LOCUS 行):

.. code:: python

    LOCUS       NC_005816               9609 bp    DNA     circular BCT 21-JUL-2008
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete
                sequence.
    ACCESSION   NC_005816
    VERSION     NC_005816.1  GI:45478711
    PROJECT     GenomeProject:10638
    ...

同样使用 ``Bio.SeqIO`` 读取文件，代码跟处理FASTA 文件类似 (详见第 :ref:`5 <chapter-Bio.SeqIO>` 章):

.. code:: python

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("NC_005816.gb", "genbank")
    >>> record
    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG',
    IUPACAmbiguousDNA()), id='NC_005816.1', name='NC_005816',
    description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.',
    dbxrefs=['Project:10638'])

你可能已经发现了一些不同之处，逐个环顾各个属性，序列字符串和上述类似，但此处 ``Bio.SeqIO`` 可自动识别序列类型 （详见第 :ref:`5 <chapter-Bio.SeqIO>` 章）:

.. code:: python

    >>> record.seq
    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG', IUPACAmbiguousDNA())

``name`` 源于 LOCUS行, ``id`` 附加了版本后缀。description源于DEFINITION 行:

.. code:: python

    >>> record.id
    'NC_005816.1'
    >>> record.name
    'NC_005816'
    >>> record.description
    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.'

GenBank 文件中per-letter annotations为空:

.. code:: python

    >>> record.letter_annotations
    {}

多数注释信息储存在 ``annotations`` 字典中，例如:

.. code:: python

    >>> len(record.annotations)
    11
    >>> record.annotations["source"]
    'Yersinia pestis biovar Microtus str. 91001'

``dbxrefs`` 列表中的数据来自 PROJECT 或DBLINK行:

.. code:: python

    >>> record.dbxrefs
    ['Project:10638']

最后也许也可能是最有意思的， ``features`` 列表以 ``SeqFeature`` 对象的形式保存了features table中的所有entries（如genes和CDS等）。

.. code:: python

    >>> len(record.features)
    29

接下来，我们将在 第 :ref:`4.3 <sec-seq_features>` 节介绍 ``SeqFeature`` 对象。

.. _sec-seq_features:

4.3  Feature, location 和 position对象
-------------------------------------------

4.3.1  SeqFeature对象
~~~~~~~~~~~~~~~~~~~~~~~~~

序列特征是描述一条序列不可或缺的部分。抛开序列本身，你需要一种方式去组织和获取关于这条序列的 “抽象” 信息。 尽管设计一个通用的类囊括序列的所有特征看似是不可能的，但是Biopython的 ``SeqFeature``
类试图尽可能多的囊括序列的所有特征。Biopython主要依据 GenBank/EMBL 特征表来设计相应的对象，认识到这一点，将有助于你更快更好的理解Biopython ``SeqFeature`` 对象。

``SeqFeature`` 对象的关键目的在于描述其相对于父序列（parent sequence，通常为 ``SeqRecord`` 对象）所处的位置（location）, 通常是介于两个positions间的一个区域（region），后续第 :ref:`4.3.2 <sec-locations>` 节将详细说明。

``SeqFeature`` 对象含大量属性，首先一一例出，然后在后续章节举例说明其用法:

**.type**
    – 用文字描述的feature类型 (如 ‘CDS’ 或 ‘gene’)。
**.location**
    – ``SeqFeature`` 在序列中所处的位置。见第 :ref:`4.3.2 <sec-locations>` 节。 ``SeqFeature`` 设计了众多针对location对象的功能，包含一系列简写的属性。

    **.ref**
        – ``.location.ref`` 简写 --location对象相关的参考序列。通常为空（None）。
    **.ref\_db**
        – ``.location.ref_db`` 简写 -- 指定 ``.ref`` 相关数据库名称。通常为空（None）。
    **.strand**
        – ``.location.strand`` 简写 -- 表示feature所处序列的strand。在双链核酸序列中，1表示正链, -1表示负链, 0 表示strand信息很重要但未知, None表示strand信息未知且不重要。蛋白和单链核酸序列为None。 

**.qualifiers**
    – 存储feature附加信息（Python字典）。键（key）为值（value）所存信息的单字简要描述，值为实际信息。比如，键为 “evidence” ，而值为 “computational (non-experimental)”。 这只是为了提醒人们注意，该feature没有被实验所证实（湿实验）。Note：为与GenBank/EMBL文件中的feature tables对应，规定.qualifiers 中值为字符串数组（即使只有一个字符串）。
**.sub\_features**
    – 只有在描述复杂位置时才使用，如 GenBank/EMBL文件中的 ‘joins’ 位置。 已被 ``CompoundLocation`` 对象取代，因此略过不提。

.. _sec-locations:

4.3.2  Positions和locations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SeqFeature`` 对象主要用于描述相对于父序列中的位置（region）信息。Region用location对象表示，通常是两个position间的范围。为了区分location和position，我们定义如下:

**position**
    – 表示位于序列中的单一位置, 可以是精确的也可以是不确定的位置（如5, 20, ``<100`` 和 ``>200`` ）。
**location**
    – 介于两个positions间的区域。比如5..20 (5到20)。

之所以特意提及这两个概念是因为我经常混淆两者。

4.3.2.1  FeatureLocation 对象
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

多数 ``SeqFeature`` 特别简单（真核基因例外），只需起点、终点以及strand信息。最基本的 ``FeatureLocation`` 对象中通常包括上述三点信息。

但实际情况未必如此简单，因为我们还需处理包含几个区域的复合locations，而且position本身很可能是不精确的。

4.3.2.2  CompoundLocation 对象
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

为了更方便的处理EMBL/GenBank文件中的 ‘join’ locations，Biopython 1.62引入 ``CompoundLocation`` 对象。

4.3.2.3  模糊Positions
^^^^^^^^^^^^^^^^^^^^^^^^

目前，我们只处理过简单position，feature location复杂因素之一就是由position本身
不准确所致。生物学中许多问题都是不确定的，比如：你通过双核苷酸priming证明了mRNA
的转录起始位点是这两个位点中的一个。这是十分有价值的发现，但困难来自于怎样表述这
个位点信息。为了处理类似情况，我们用模糊位点（fuzzy position）表示。根据fuzzy 
position的不同，我们用5个类分别描述:

**ExactPosition**
    – 精确位点，用一个数字表示。从该对象的 ``position`` 属性可得知精确位点信息。
**BeforePosition**
    – 位于某个特定位点前。如 ```<13'`` , 在GenBank/EMBL中代表实际位点位于13之前。从该对象的 ``position`` 属性可得知上边界信息。 
**AfterPosition**
    – 与 ``BeforePosition`` 相反,如 ```>13'`` , 在GenBank/EMBL中代表实际位点位于13以后。从该对象的 ``position`` 属性可获知下边界信息。
**WithinPosition**
    – 介于两个特定位点之间，偶尔在GenBank/EMBL locations用到。如 ‘(1.5)’, GenBank/EMBL中代表实际位点位于1到5之间。该对象需要两个position属性表示，第一个 ``position`` 表示下边界（本例为1）， ``extension`` 表示上边界与下边界的差值（本例为4）。因此在WithinPosition中， ``object.position`` 表示下边界， ``object.position + object.extension`` 表示上边界。
**OneOfPosition**
    – 表示几个位点中的一个（GenBank/EMBL文件中偶尔能看到），比如在基因起始位点不明确或者有两个候选位点的时候可以使用，或者用于明确表示两个相关基因特征时使用。 
**UnknownPosition**
    – 代表未知位点。在GenBank/EMBL文件中没有使用，对应 UniProt中的 ‘?’ feature坐标。

举例说明创建一个fuzzy end points:

.. code:: python

    >>> from Bio import SeqFeature
    >>> start_pos = SeqFeature.AfterPosition(5)
    >>> end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
    >>> my_location = SeqFeature.FeatureLocation(start_pos, end_pos)

Note：Biopython 1.59以后，fuzzy-locations有修改, 特别是BetweenPosition和
WithinPosition，现在必须显示用整数表示。起点为较小值，终点则为较大值。

print输出 ``FeatureLocation`` 对象，可看到简洁的结果:

.. code:: python

    >>> print my_location
    [>5:(8^9)]

也可通过start和end属性得到fuzzy position的起始/终止位点:

.. code:: python

    >>> my_location.start
    AfterPosition(5)
    >>> print my_location.start
    >5
    >>> my_location.end
    BetweenPosition(9, left=8, right=9)
    >>> print my_location.end
    (8^9)

如果你只想获取数字，不理会模糊positions，则可将fuzzy position强制转换成一个整数:

.. code:: python

    >>> int(my_location.start)
    5
    >>> int(my_location.end)
    9

为了兼容旧版Biopython，保留了整数形式的 ``nofuzzy_start`` and ``nofuzzy_end`` :

.. code:: python

    >>> my_location.nofuzzy_start
    5
    >>> my_location.nofuzzy_end
    9

Notice：上述例子只是为了帮助你理解fuzzy locations。

相似的，如果要建立一个精确location，只需将整数传递给 ``FeaturePosition``
构造函数, 即可建立 ``ExactPosition`` 对象:

.. code:: python

    >>> exact_location = SeqFeature.FeatureLocation(5, 9)
    >>> print exact_location
    [5:9]
    >>> exact_location.start
    ExactPosition(5)
    >>> int(exact_location.start)
    5
    >>> exact_location.nofuzzy_start
    5

以上是Biopython处理fuzzy position的实现方法。希望读者能体会之所以这样设计，都是为了使用上的方便（至少不比精确位点复杂）

4.3.2.4  Location testing
^^^^^^^^^^^^^^^^^^^^^^^^^

可用Python关键词 ``in`` 检验某个碱基或氨基酸残基的父坐标是否位于
feature/location中。

假定你想知道某个SNP位于哪个feature里，并知道该SNP的索引位置是4350（Python 计数）。一个简单的实现方案是用循环遍历所有features:

.. code:: python

    >>> from Bio import SeqIO
    >>> my_snp = 4350
    >>> record = SeqIO.read("NC_005816.gb", "genbank")
    >>> for feature in record.features:
    ...     if my_snp in feature:
    ...         print feature.type, feature.qualifiers.get('db_xref')
    ...
    source ['taxon:229193']
    gene ['GeneID:2767712']
    CDS ['GI:45478716', 'GeneID:2767712']

Note： GenBank /EMBL 文件中的 gene 和CDS features（ ``join`` ）只包含外显子，不含内含子。

4.3.3  使用feature 或 location描述序列
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SeqFeature`` 或 location object对象并没有直接包含任何序列，只是可根据储存的location (见第 :ref:`4.3.2 <sec-locations>` 节)，从父序列中取得。例如：某一短基因位于负链5:18位置，由于GenBank/EMBL文件以1开始计数，Biopython中表示为 ``complement(6..18)`` :

.. code:: python

    >>> from Bio.Seq import Seq
    >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
    >>> example_parent = Seq("ACCGAGACGGCAAAGGCTAGCATAGGTATGAGACTTCCTTCCTGCCAGTGCTGAGGAACTGGGAGCCTAC")
    >>> example_feature = SeqFeature(FeatureLocation(5, 18), type="gene", strand=-1)

你可以用切片从父序列截取5:18,然后取反向互补序列。如果是Biopython 1.59或以后版本，可使用如下方法:

.. code:: python

    >>> feature_seq = example_parent[example_feature.location.start:example_feature.location.end].reverse_complement()
    >>> print feature_seq
    AGCCTTTGCCGTC

不过在处理复合 features (joins)时，此法相当繁琐。此时可以使用 ``SeqFeature`` 对象的 ``extract`` 方法处理:

.. code:: python

    >>> feature_seq = example_feature.extract(example_parent)
    >>> print feature_seq
    AGCCTTTGCCGTC

``SeqFeature`` 或 location对象的长度等同于所表示序列的长度。

.. code:: python

    >>> print example_feature.extract(example_parent)
    AGCCTTTGCCGTC
    >>> print len(example_feature.extract(example_parent))
    13
    >>> print len(example_feature)
    13
    >>> print len(example_feature.location)
    13

简单 ``FeatureLocation`` 对象的长度等于终止osition减去起始position的差值；而 ``CompoundLocation`` 的长度则为各片段长度之和。

4.4  References
---------------

对一条序列的注释还包括参考文献（reference），Biopython通过
``Bio.SeqFeature.Reference`` 对象来储存相关的文献信息。

References属性储存了 ``期刊名`` 、 ``题名`` 、 ``作者`` 等信息。此外还包括 ``medline_id`` 、 ``pubmed_id`` 以及 ``comment`` 。

通常reference 也有 ``location`` 对象，便于文献涉及研究对象在序列中的定位。该 ``location`` 有可能是一个fuzzy location（见第 :ref:`4.3.2 <sec-locations>` 节）。

文献对象都以列表储存在 ``SeqRecord`` 对象的 ``annotations`` 字典中。 字典的键为 “references”。reference对象也是为了方便处理文献而设计，希望能满足各种使用需求。

.. _sec-SeqRecord-format:

4.5  格式化方法
----------------------

``SeqRecord`` 类中的 ``format()`` 能将字符串转换成被 ``Bio.SeqIO`` 支持的格式，如FASTA:

.. code:: python

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_protein

    record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD" \
                          +"GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK" \
                          +"NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM" \
                          +"SSAC", generic_protein),
                       id="gi|14150838|gb|AAK54648.1|AF376133_1",
                       description="chalcone synthase [Cucumis sativus]")
                       
    print record.format("fasta")

输出为:

.. code:: python

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC

``format`` 方法接收单个必选参数，小写字母字符串是 ``Bio.SeqIO`` 模块支持的输出格式 (见第 :ref:`5 <chapter-Bio.SeqIO>` 章)。然而，此 ``format()`` 方法并不适用于包含多条序列的文件格式 (如多序列比对格式)（详见第 :ref:`5.5.4 <sec-Bio.SeqIO-and-StringIO>` 节）。

.. _sec-SeqRecord-slicing:

4.6  SeqRecord切片
------------------------

通过切片截取 ``SeqRecord`` 的部分序列可得到一条新的 ``SeqRecord`` 。此处需引起注意的是per-letter annotations也被取切片, 但新序列中的features保持不变 (locations相应调整)。

以前述Genbank文件为例:

.. code:: python

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("NC_005816.gb", "genbank")

.. code:: python

    >>> record
    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG',
    IUPACAmbiguousDNA()), id='NC_005816.1', name='NC_005816',
    description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.',
    dbxrefs=['Project:10638'])

.. code:: python

    >>> len(record)
    9609
    >>> len(record.features)
    41

本例中，我们关注 ``YP_pPCP05`` 质粒上的 ``pim`` 基因。从GenBank文件可直接看出 ``pim`` gene/CDS location是 ``4343..4780`` （相应的Python 位置是 ``4342:4780`` ）。Location信息位于GenBank文件第12和13 entries中, 由于python以0开始计数，因此python中，它们是 ``features`` 列表中的 entries 11和12:

.. code:: python

    >>> print record.features[20]
    type: gene
    location: [4342:4780](+)
    qualifiers: 
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    <BLANKLINE>

.. code:: python

    >>> print record.features[21]
    type: CDS
    location: [4342:4780](+)
    qualifiers: 
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity ...']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSH...']

从父记录中取切片（4300 到 4800），观测所得到的features数量:

.. code:: python

    >>> sub_record = record[4300:4800]

.. code:: python

    >>> sub_record
    SeqRecord(seq=Seq('ATAAATAGATTATTCCAAATAATTTATTTATGTAAGAACAGGATGGGAGGGGGA...TTA',
    IUPACAmbiguousDNA()), id='NC_005816.1', name='NC_005816',
    description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.',
    dbxrefs=[])

.. code:: python

    >>> len(sub_record)
    500
    >>> len(sub_record.features)
    2

子记录（sub_record）只包括两个features, 分别是 ``YP_pPCP05`` 质粒的gene和CDS:

.. code:: python

    >>> print sub_record.features[0]
    type: gene
    location: [42:480](+)
    qualifiers: 
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    <BLANKLINE>

.. code:: python

    >>> print sub_record.features[20]
    type: CDS
    location: [42:480](+)
    qualifiers: 
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity ...']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSH...']

注意：locations已被调整至对应生成的新父序列!

尽可能灵敏和直观地获取子记录的相关特征（和任意的per-letter annotation），但是对于其余注释，Biopython无法判断是否仍然适用于子记录。因此子记录忽略了 ``annotations`` 和 ``dbxrefs`` 以避免引起歧义。

.. code:: python

    >>> sub_record.annotations
    {}
    >>> sub_record.dbxrefs
    []

为了便于实际操作，子记录保留了 ``id`` , ``name`` 和 ``description`` :

.. code:: python

    >>> sub_record.id
    'NC_005816.1'
    >>> sub_record.name
    'NC_005816'
    >>> sub_record.description
    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.'

上述例子很好的展示了问题，由于子记录不包括完整的质粒序列，因此description是错的。我们可以将子记录看做是截短版的GenBank文件，可用第 :ref:`4.5 <sec-SeqRecord-format>` 节中所述 ``format`` 方法纠正：
:

.. code:: python

    >>> sub_record.description = "Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial."
    >>> print sub_record.format("genbank")
    ...

FASTQ例子参见第 :ref:`18.1.7 <sec-FASTQ-slicing-off-primer>` 节
和第 :ref:`18.1.8 <sec-FASTQ-slicing-off-adaptor>` 节（此例中per-letter annotations (read质量分数) 也被取切片）。

.. _sec-SeqRecord-addition:

4.7  SeqRecord对象相加
-----------------------------

``SeqRecord`` 对象可相加得到一个新的 ``SeqRecord`` 。注意：per-letter annotations也相加, features (locations 调整)；而其它annotation 保持不变(如id、name和description)。

以FASTQ 文件中的第一条记录为例说明per-letter annotation （第 :ref:`5 <chapter-Bio.SeqIO>` 章详细介绍 ``SeqIO`` 函数）:

.. code:: python

    >>> from Bio import SeqIO
    >>> record = SeqIO.parse("example.fastq", "fastq").next()
    >>> len(record)
    25
    >>> print record.seq
    CCCTTCTTGTCTTCAGCGTTTCTCC

.. code:: python

    >>> print record.letter_annotations["phred_quality"]
    [26, 26, 18, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 22, 26, 26, 26, 26,
    26, 26, 26, 23, 23]

假设上述序列数据来自Roche 454测序, 你根据其它信息得知 ``TTT`` 应该是 ``TT`` 。此时可分别用切片提取第三个 ``T`` 前后的序列（ ``SeqRecord`` ）:

.. code:: python

    >>> left = record[:20]
    >>> print left.seq
    CCCTTCTTGTCTTCAGCGTT
    >>> print left.letter_annotations["phred_quality"]
    [26, 26, 18, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 22, 26, 26, 26, 26]
    >>> right = record[21:]
    >>> print right.seq
    CTCC
    >>> print right.letter_annotations["phred_quality"]
    [26, 26, 23, 23]

两部分相加:

.. code:: python

    >>> edited = left + right
    >>> len(edited)
    24
    >>> print edited.seq
    CCCTTCTTGTCTTCAGCGTTCTCC

.. code:: python

    >>> print edited.letter_annotations["phred_quality"]
    [26, 26, 18, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 22, 26, 26, 26, 26,
    26, 26, 23, 23]

很容易和直观吧！上述两步可合并:

.. code:: python

    >>> edited = record[:20] + record[21:]

现在以GenBank文件（假定是环状基因组）为例说明features:

.. code:: python

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("NC_005816.gb", "genbank")

.. code:: python

    >>> record
    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG',
    IUPACAmbiguousDNA()), id='NC_005816.1', name='NC_005816',
    description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.',
    dbxrefs=['Project:10638'])

.. code:: python

    >>> len(record)
    9609
    >>> len(record.features)
    41
    >>> record.dbxrefs
    ['Project:58037']

.. code:: python

    >>> record.annotations.keys()
    ['comment', 'sequence_version', 'source', 'taxonomy', 'keywords', 'references',
    'accessions', 'data_file_division', 'date', 'organism', 'gi']

可改变起点:

.. code:: python

    >>> shifted = record[2000:] + record[:2000]

.. code:: python

    >>> shifted
    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA',
    IUPACAmbiguousDNA()), id='NC_005816.1', name='NC_005816',
    description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.',
    dbxrefs=[])

.. code:: python

    >>> len(shifted)
    9609

Note: 上述方法并不完美（丢失了数据库交叉引用dbxrefs 和源feature）:

.. code:: python

    >>> len(shifted.features)
    40
    >>> shifted.dbxrefs
    []
    >>> shifted.annotations.keys()
    []

这是因为 ``SeqRecord`` 切片对 annotation 保留非常谨慎 (错误保留 annotation 可能引起大问题)。如果你想保留数据库的交叉引用dbxrefs和其余annotations 字典必须明确说明，才能得以保留:

.. code:: python

    >>> shifted.dbxrefs = record.dbxrefs[:]
    >>> shifted.annotations = record.annotations.copy()
    >>> shifted.dbxrefs
    ['Project:10638']
    >>> shifted.annotations.keys()
    ['comment', 'sequence_version', 'source', 'taxonomy', 'keywords', 'references',
    'accessions', 'data_file_division', 'date', 'organism', 'gi']

Note: 此例中序列record的identifiers也应调整（因为NCBI的reference链接的是未经修改的 *原始* 序列）。

.. _sec-SeqRecord-reverse-complement:

4.8  反向互补SeqRecord对象
--------------------------------------------

为消除序列反向互补后annotation改变带来的困难，Biopython 1.57 ``SeqRecord`` 对象加入了 ``reverse_complement`` 方法。这也成为Biopython 1.57的新特性之一。

序列用Seq对象中的reverse_complement方法反向互补。Features随location而改变，strand也被重新计算。复制并反转per-letter-annotation（通常情况下这种做法比较合适，如对质量分数注释的反转）。然而多数annotation的转变却存有问题。

比如record ID是accession号，该accession不应被用于反向互补序列。默认identifier转换可导致后续分析中的轻度数据损坏。因此 ``SeqRecord`` 的id、
name、description、annotations和dbxrefs默认不变。

``SeqRecord`` 对象的 ``reverse_complement`` 法用多个可选参数以对应record的属性。将这些参数设为 ``True`` 表示复制旧值；而 ``False`` 意为用缺省值替换旧值。当然也可自定义新值。

举例:

.. code:: python

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("NC_005816.gb", "genbank")
    >>> print record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)
    NC_005816.1 9609 41 1 11

反向互补该record并给ID赋予新值 - 注意：多数annotation丢失，而features仍在:

.. code:: python

    >>> rc = record.reverse_complement(id="TESTING")
    >>> print rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)
    TESTING 9609 41 0 0

