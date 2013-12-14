第10章 Swiss-Prot和ExPASy
=================================

10.1  解析Swiss-Prot文件
------------------------------

Swiss-Prot
( `http://www.expasy.org/sprot <http://www.expasy.org/sprot>`__ )是一个
蛋白质序列数据库。 Biopython能够解析纯文本的Swiss-Prot文件,
这种格式也被Swiss-Prot、TrEMBL和PIRPSD的UniProt数据库使用。然而我们并
不支持UniProKB的XML格式文件。

10.1.1  Parsing Swiss-Prot records
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

在 \ `5.3.2 <#sec:SeqIO_ExPASy_and_SwissProt>`__ 章节中, 我们描述过怎样将一个
Swiss-Prot记录中的序列提出来作为一个 ``SeqRecord`` 对象。 此外，你可以将
Swiss-Prot记录存到  ``Bio.SwissProt.Record`` 对象, 这实际上存储了Swiss-Prot记录
中所包含的的全部信息。在这部分我们将介绍怎样从一个Swiss-Prot文件中提
取 ``Bio.SwissProt.Record`` 对象。

为了解析Swiss-Prot记录，我们首先需要得到一个Swiss-Prot记录文件。根据该Swiss-Prot
记录的储存位置和储存方式，获取该记录文件的方式也有所不同：

-  本地打开Swiss-Prot文件：

   .. code:: python
      
       >>> handle = open("myswissprotfile.dat")

-  打开使用gzip压缩的Swiss-Prot文件：

   .. code:: python

       >>> import gzip
       >>> handle = gzip.open("myswissprotfile.dat.gz")

-  在线打开Swiss-Prot文件：

   .. code:: python

       >>> import urllib
       >>> handle = urllib.urlopen("http://www.somelocation.org/data/someswissprotfile.dat")

-  从ExPASy数据库在线打开Swiss-Prot文件
   (见 `10.5.1 <#subsec:expasy_swissprot>`__ 章节):

   .. code:: python

       >>> from Bio import ExPASy
       >>> handle = ExPASy.get_sprot_raw(myaccessionnumber)

对于解析来说，关键点在于Swiss-Prot格式的数据，而不是获取它的方式。

我们可以用 \ `5.3.2 <#sec:SeqIO_ExPASy_and_SwissProt>`__ 章节中描述的方式，
通过 ``Bio.SeqIO`` 来获取格式未知的 ``SeqRecord`` 对象。此外，我们也可以
用 ``Bio.SwissProt`` 来获取更加匹配基本文件格式的 ``Bio.SwissProt.Record`` 对象。

我们使用 ``read()`` 函数来从文件中读取一个Swiss-Prot记录：

.. code:: python

    >>> from Bio import SwissProt
    >>> record = SwissProt.read(handle)

该函数只适用于仅存储了一个Swiss-Prot记录的文件，而当文件中没有或存在
多个记录时使用该函数，会出现 ``ValueError`` 提示。

现在我们可以输出一些与这些记录相关的信息：

.. code:: python

    >>> print record.description
    'RecName: Full=Chalcone synthase 3; EC=2.3.1.74; AltName: Full=Naringenin-chalcone synthase 3;'
    >>> for ref in record.references:
    ...     print "authors:", ref.authors
    ...     print "title:", ref.title
    ...
    authors: Liew C.F., Lim S.H., Loh C.S., Goh C.J.;
    title: "Molecular cloning and sequence analysis of chalcone synthase cDNAs of
    Bromheadia finlaysoniana.";
    >>> print record.organism_classification
    ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', ..., 'Bromheadia']

为了解析包含多个Swiss-Prot记录的文件，我们使用 ``parse`` 函数。这个函数能够让我们对
文件中的记录进行循环迭代操作。

比如，我们要解析整个Swiss-Prot数据库并且收集所有的描述。你可以从
`ExPAYs FTP site <ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz>`__ 
下载这些gzip压缩文件 ``uniprot_sprot.dat.gz`` (大约 300MB)。文件中含有
``uniprot_sprot.dat`` 一个文件(至少1.5GB)。

如同这一部分刚开始所描述的，你可以按照如下所示的方法使用python
的 ``gzip`` 模块打开并解压 ``.gz`` 文件:

.. code:: python

    >>> import gzip
    >>> handle = gzip.open("uniprot_sprot.dat.gz")

然而，解压一个大文件比较耗时，而且每次用这种方式打开一个
文件都是比较慢的。所以，如果你有空闲的硬盘空间并且在
最开始就在硬盘里通过解压到来得到 ``uniprot_sprot.dat`` ，这样能够在以后就可以像平常那样来打开文件：

.. code:: python

    >>> handle = open("uniprot_sprot.dat")

到2009年6月为止，从ExPASy下载下来的整个Swiss-Prot数据库一共
有468851个Swiss-Prot记录，一种建立关于这些记录的描述列表的
间接方式就是使用一种列表解析：

.. code:: python

    >>> from Bio import SwissProt
    >>> handle = open("uniprot_sprot.dat")
    >>> descriptions = [record.description for record in SwissProt.parse(handle)]
    >>> len(descriptions)
    468851
    >>> descriptions[:5]
    ['RecName: Full=Protein MGF 100-1R;',
     'RecName: Full=Protein MGF 100-1R;',
     'RecName: Full=Protein MGF 100-1R;',
     'RecName: Full=Protein MGF 100-1R;',
     'RecName: Full=Protein MGF 100-2L;']

或者对记录迭代器使用for循环：

.. code:: python

    >>> from Bio import SwissProt
    >>> descriptions = []
    >>> handle = open("uniprot_sprot.dat")
    >>> for record in SwissProt.parse(handle):
    ...     descriptions.append(record.description)
    ...
    >>> len(descriptions)
    468851

由于输入文件太大，这两种方法在我的新台式机上花费大约十一分钟（用解压好的
``uniprot_sprot.dat`` 作为输入文件）。

从Swiss-Prot记录中提取任何你想要的信息也同样简单。比如你想看看一个
Swiss-Prot记录中的成员，就输入：

.. code:: python

    >>> dir(record)
    ['__ doc__ ', '__ init__ ', '__ module__ ', 'accessions', 'annotation_update',
    'comments', 'created', 'cross_references', 'data_class', 'description',
    'entry_name', 'features', 'gene_name', 'host_organism', 'keywords',
    'molecule_type', 'organelle', 'organism', 'organism_classification',
    'references', 'seqinfo', 'sequence', 'sequence_length',
    'sequence_update', 'taxonomy_id']

10.1.2  解析Swiss-Prot关键词和分类列表
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Swiss-Prot也会提供一个 ``keywlist.txt`` 文件，该文件列出了Swiss-Prot中所用到
的关键词和分类。其中所包含的词条形式如下：

.. code:: python

    ID   2Fe-2S.
    AC   KW-0001
    DE   Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron
    DE   atoms complexed to 2 inorganic sulfides and 4 sulfur atoms of
    DE   cysteines from the protein.
    SY   Fe2S2; [2Fe-2S] cluster; [Fe2S2] cluster; Fe2/S2 (inorganic) cluster;
    SY   Di-mu-sulfido-diiron; 2 iron, 2 sulfur cluster binding.
    GO   GO:0051537; 2 iron, 2 sulfur cluster binding
    HI   Ligand: Iron; Iron-sulfur; 2Fe-2S.
    HI   Ligand: Metal-binding; 2Fe-2S.
    CA   Ligand.
    //
    ID   3D-structure.
    AC   KW-0002
    DE   Protein, or part of a protein, whose three-dimensional structure has
    DE   been resolved experimentally (for example by X-ray crystallography or
    DE   NMR spectroscopy) and whose coordinates are available in the PDB
    DE   database. Can also be used for theoretical models.
    HI   Technical term: 3D-structure.
    CA   Technical term.
    //
    ID   3Fe-4S.
    ...

文件中的词条可以通过使用 ``Bio.SwissProt.KeyWList`` 模块中的 ``parse`` 函数
来解析，并且每一个词条都会被存储在名为 ``Bio.SwissProt.KeyWList.Record`` 的
python字典里。


.. code:: python

    >>> from Bio.SwissProt import KeyWList
    >>> handle = open("keywlist.txt")
    >>> records = KeyWList.parse(handle)
    >>> for record in records:
    ...     print record['ID']
    ...     print record['DE']

这些命令行将会输出：

.. code:: python

    2Fe-2S.
    Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron atoms
    complexed to 2 inorganic sulfides and 4 sulfur atoms of cysteines from the
    protein.
    ...

10.2  解析Prosite记录
-----------------------------

Prosite是一个包含了蛋白质结构域、蛋白家族、功能位点以及识别它们的模式和图
谱，而且它是和Swiss-Prot同时开发出来的。
在Biopython中，Prosite记录是由 ``Bio.ExPASy.Prosite.Record`` 类来表示的，
其中的成员与该Prosite记录中的不同区域相对应。

一般来说，一个Prosite文件可以包含多个Prosite记录。比如，从 `ExPASy FTP
site <ftp://ftp.expasy.org/databases/prosite/prosite.dat>`__ 网站下载
下来的、容纳了整个Prosite记录的 ``prosite.dat`` 文件，含有2073条记录（2007年12月发布的第20.24版本）。
为了解析这样一个文件，我们再次使用一个迭代器：

.. code:: python

    >>> from Bio.ExPASy import Prosite
    >>> handle = open("myprositefile.dat")
    >>> records = Prosite.parse(handle)

现在我们可以逐个提取这些记录并输出其中一些信息。比如，使用包含整个Prosite数据库的
文件将会使我们找到如下等信息：

.. code:: python

    >>> from Bio.ExPASy import Prosite
    >>> handle = open("prosite.dat")
    >>> records = Prosite.parse(handle)
    >>> record = records.next()
    >>> record.accession
    'PS00001'
    >>> record.name
    'ASN_GLYCOSYLATION'
    >>> record.pdoc
    'PDOC00001'
    >>> record = records.next()
    >>> record.accession
    'PS00004'
    >>> record.name
    'CAMP_PHOSPHO_SITE'
    >>> record.pdoc
    'PDOC00004'
    >>> record = records.next()
    >>> record.accession
    'PS00005'
    >>> record.name
    'PKC_PHOSPHO_SITE'
    >>> record.pdoc
    'PDOC00005'

如果你想知道有多少条Prosite记录，你可以输入：

.. code:: python

    >>> from Bio.ExPASy import Prosite
    >>> handle = open("prosite.dat")
    >>> records = Prosite.parse(handle)
    >>> n = 0
    >>> for record in records: n+=1
    ...
    >>> print n
    2073

为了从这些数据中读取某一条特定的记录，可以使用 ``read`` 函数：

.. code:: python

    >>> from Bio.ExPASy import Prosite
    >>> handle = open("mysingleprositerecord.dat")
    >>> record = Prosite.read(handle)

如果并不存在或存在多个你想要找的Prosite记录时，这个函数将会输出一个“ValueError”提示。

10.3  解析Prosite文件记录
-------------------------------------------

在上述的Prosite示例中，像 ``'PDOC00001'`` 、 ``'PDOC00004'`` 、 ``'PDOC00005'`` 等这样的编号指的就
是Prosite文件。Prosite文件记录可以以单个文件（ ``prosite.doc`` ）的形式从ExPASy获取，并
且该文件包含了所有Prosite文档记录。

我们使用 ``Bio.ExPASy.Prodoc`` 中的解析器来解析这些Prosite文档记录。比如，为了生成一个包含所有
Prosite文档记录的编号列表，你可以使用：

.. code:: python

    >>> from Bio.ExPASy import Prodoc
    >>> handle = open("prosite.doc")
    >>> records = Prodoc.parse(handle)
    >>> accessions = [record.accession for record in records]

进一步可以使用 ``read()`` 函数来对这些数据中具体某一条文档记录来进行查询。

10.4  解析酶记录
----------------------------

ExPASy的酶数据库是一个关于酶的系统命名信息的数据库。如下所示是一个比较典型的酶的记录

.. code:: python

    ID   3.1.1.34
    DE   Lipoprotein lipase.
    AN   Clearing factor lipase.
    AN   Diacylglycerol lipase.
    AN   Diglyceride lipase.
    CA   Triacylglycerol + H(2)O = diacylglycerol + a carboxylate.
    CC   -!- Hydrolyzes triacylglycerols in chylomicrons and very low-density
    CC       lipoproteins (VLDL).
    CC   -!- Also hydrolyzes diacylglycerol.
    PR   PROSITE; PDOC00110;
    DR   P11151, LIPL_BOVIN ;  P11153, LIPL_CAVPO ;  P11602, LIPL_CHICK ;
    DR   P55031, LIPL_FELCA ;  P06858, LIPL_HUMAN ;  P11152, LIPL_MOUSE ;
    DR   O46647, LIPL_MUSVI ;  P49060, LIPL_PAPAN ;  P49923, LIPL_PIG   ;
    DR   Q06000, LIPL_RAT   ;  Q29524, LIPL_SHEEP ;
    //

在这个例子中，第一行显示了脂蛋白脂肪酶（第二行）的酶编号(EC, Enzyme Commission)。
脂蛋白脂肪酶其他的名称有“清除因子脂肪酶”和“甘油二脂脂肪酶”（第三行至第五行）。
开头为“CA”的那一行显示了该酶的催化活性。评论行开头为“CC”。“PR”行显示了对应Prosite
文档记录的参考，以及“DR”行显示了Swiss-Prot记录的参考。
然而并不是所有的词条都必需出现在酶记录当中。

在Biopython中，一个酶记录由 ``Bio.ExPASy.Enzyme.Record`` 类来代表。这个记录源于对应
于酶相关文件中所用到的双字母编码的python字典和哈希键。为了阅读含有一个酶记录的酶文件，
你可以使用 ``Bio.ExPASy.Enzyme`` 中的 ``read`` 函数：

.. code:: python

    >>> from Bio.ExPASy import Enzyme
    >>> handle = open("lipoprotein.txt")
    >>> record = Enzyme.read(handle)
    >>> record["ID"]
    '3.1.1.34'
    >>> record["DE"]
    'Lipoprotein lipase.'
    >>> record["AN"]
    ['Clearing factor lipase.', 'Diacylglycerol lipase.', 'Diglyceride lipase.']
    >>> record["CA"]
    'Triacylglycerol + H(2)O = diacylglycerol + a carboxylate.'
    >>> record["PR"]
    ['PDOC00110']

.. code:: python

    >>> record["CC"]
    ['Hydrolyzes triacylglycerols in chylomicrons and very low-density lipoproteins
    (VLDL).', 'Also hydrolyzes diacylglycerol.']
    >>> record["DR"]
    [['P11151', 'LIPL_BOVIN'], ['P11153', 'LIPL_CAVPO'], ['P11602', 'LIPL_CHICK'],
    ['P55031', 'LIPL_FELCA'], ['P06858', 'LIPL_HUMAN'], ['P11152', 'LIPL_MOUSE'],
    ['O46647', 'LIPL_MUSVI'], ['P49060', 'LIPL_PAPAN'], ['P49923', 'LIPL_PIG'],
    ['Q06000', 'LIPL_RAT'], ['Q29524', 'LIPL_SHEEP']]

如果没有找到或者找到多个酶记录时， ``read`` 函数会反馈一个ValueError提示。

所有酶记录都可以从 `ExPASy FTP site <ftp://ftp.expasy.org/databases/enzyme/enzyme.dat>`__ 网站下载
为单个文件（ ``enzyme.dat`` ），该文件包含了4877个记录（2009年3月发布的第三版）。为了打开含有多个
酶记录的文件，你可以使用 ``Bio.ExPASy.Enzyme`` 中的 ``parse`` 函数来获得一个迭代器：

.. code:: python

    >>> from Bio.ExPASy import Enzyme
    >>> handle = open("enzyme.dat")
    >>> records = Enzyme.parse(handle)

我们现在每次都可以对这些记录进行迭代。比如我们可以对那些已有的酶记录做一个EC编号列表：

.. code:: python

    >>> ecnumbers = [record["ID"] for record in records]

10.5  Accessing the ExPASy server
---------------------------------

Swiss-Prot、Prosite和Prosite文档记录可以从
`http://www.expasy.org <http://www.expasy.org>`__ 的ExPASy网络服务器下载到。在ExPASy服
务器上可以进行六种查询：

**get\_prodoc\_entry**
    下载一个HTML格式的Prosite文档记录
**get\_prosite\_entry**
    下载一个HTML格式的Prosite记录
**get\_prosite\_raw**
    下载一个原始格式的Prosite或Prosite文档记录
**get\_sprot\_raw**
    下载一个原始格式的Swiss-Prot记录
**sprot\_search\_ful**
    搜索一个Swiss-Prot记录
**sprot\_search\_de**
    搜索一个Swiss-Prot记录

为了从python脚本来访问该网络服务器，我们可以使用 ``Bio.ExPASy`` 模块。

10.5.1  获取一个Swiss-Prot记录
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

现在让我们来寻找一个关于兰花的查儿酮合成酶（对于寻找和兰花相关的有趣东西的理由
请看 \ `2.3 <#sec:orchids>`__ 章节）。查儿酮合成酶参与了植物中类黄酮的生物合成，
类黄酮能够合成包含色素和UV保护分子等物质。

如果你要对Swiss-Prot进行搜索，你可以找到三个关于查儿酮合成酶的兰花蛋白，id编号
为O23729, O23730, O23731。现在我们要写一个能够获取这些蛋白并能够找到一些有趣
的信息的脚本。

首先，我们使用 ``Bio.ExPASy`` 中的 ``get_sprot_raw()`` 函数来获取这些记录。这个函
数非常棒，因为你可以给它提供一个id然后得到一个原始文本记录（不会受到HTML的干扰）。
然后我们可以使用 ``Bio.SwissProt.read`` 来提取对应的Swiss-Prot记录，也可以使用 ``Bio.SeqIO.read`` 来
得到一个序列记录SeqRecord。下列代码能够实现我刚刚提到的任务：

.. code:: python

    >>> from Bio import ExPASy
    >>> from Bio import SwissProt

    >>> accessions = ["O23729", "O23730", "O23731"]
    >>> records = []

    >>> for accession in accessions:
    ...     handle = ExPASy.get_sprot_raw(accession)
    ...     record = SwissProt.read(handle)
    ...     records.append(record)

如果你提供给 ``ExPASy.get_sprot_raw`` 的编号并不存在，那么 ``SwissProt.read(handle)`` 会反
馈一个 ``ValueError`` 提示。你可以根据 ``ValueException`` 异常来找到无效的编号：

.. code:: python

    >>> for accession in accessions:
    ...     handle = ExPASy.get_sprot_raw(accession)
    ...     try:
    ...         record = SwissProt.read(handle)
    ...     except ValueException:
    ...         print "WARNING: Accession %s not found" % accession
    ...     records.append(record)

10.5.2  搜索Swiss-Prot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

现在，你可以察觉到我已经提前知道了这个记录的编号。的确， ``get_sprot_raw()`` 需要一个词条或者编号。
当你并没有编号或者词条的时候，你可使用 ``sprot_search_de()`` 或者 ``sprot_search_ful()`` 函数来解决问题。

``sprot_search_de()`` 在ID, DE, GN, OS和OG行进行搜索；
``sprot_search_ful()`` 则在所有行进行搜索。具体相关细节分别在
`http://www.expasy.org/cgi-bin/sprot-search-de <http://www.expasy.org/cgi-bin/sprot-search-de>`__ 
和
`http://www.expasy.org/cgi-bin/sprot-search-ful <http://www.expasy.org/cgi-bin/sprot-search-ful>`__ 上有说明。
注意它们的默认情况下并不搜索TrEMBL（参数为 ``trembl`` ）。还要注意它们返回的是html网页，然而编号却可以很容易从中得到：

.. code:: python

    >>> from Bio import ExPASy
    >>> import re

    >>> handle = ExPASy.sprot_search_de("Orchid Chalcone Synthase")
    >>> # or:
    >>> # handle = ExPASy.sprot_search_ful("Orchid and {Chalcone Synthase}")
    >>> html_results = handle.read()
    >>> if "Number of sequences found" in html_results:
    ...     ids = re.findall(r'HREF="/uniprot/(\w+)"', html_results)
    ... else:
    ...     ids = re.findall(r'href="/cgi-bin/niceprot\.pl\?(\w+)"', html_results)

10.5.3  获取Prosite和Prosite文档记录
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

我们可以得到HTML格式和原始格式的Prosite和Prosite文档记录。为了用biopython解析Prosite和Prosite文档记录，
你应该使用原始格式的记录。而对于其他的目的，你或许会对HTML格式感兴趣。

为了获取一个原始格式的Prosite或者Prosite文档的记录，请使用 ``get_prosite_raw()`` 。
例如，为了下载一个prosite记录并以原始格式输出，你可以使用：

.. code:: python

    >>> from Bio import ExPASy
    >>> handle = ExPASy.get_prosite_raw('PS00001')
    >>> text = handle.read()
    >>> print text

为了获取一个Prosite记录并将其解析成一个 ``Bio.Prosite.Record`` 对象，请使用：

.. code:: python

    >>> from Bio import ExPASy
    >>> from Bio import Prosite
    >>> handle = ExPASy.get_prosite_raw('PS00001')
    >>> record = Prosite.read(handle)

该函数也可以用于获取Prosite文档记录并解析到一个 ``Bio.ExPASy.Prodoc.Record`` 对象：

.. code:: python

    >>> from Bio import ExPASy
    >>> from Bio.ExPASy import Prodoc
    >>> handle = ExPASy.get_prosite_raw('PDOC00001')
    >>> record = Prodoc.read(handle)

对于不存在的编号， ``ExPASy.get_prosite_raw`` 返回一个空字符串。当遇到空字符
串， ``Prosite.read`` 和 ``Prodoc.read`` 会反馈一个ValueError错误。你可以
根据这些错误异常提示来找到无效的编号。

``get_prosite_entry()`` 和 ``get_prodoc_entry()`` 函数可用于下载HTML格式的Prosite和Prosite文档记录。
为了生成展示单个Prosite记录的网页，你可以使用：

.. code:: python

    >>> from Bio import ExPASy
    >>> handle = ExPASy.get_prosite_entry('PS00001')
    >>> html = handle.read()
    >>> output = open("myprositerecord.html", "w")
    >>> output.write(html)
    >>> output.close()

类似地，Prosite文档文本的网页展示如下：

.. code:: python

    >>> from Bio import ExPASy
    >>> handle = ExPASy.get_prodoc_entry('PDOC00001')
    >>> html = handle.read()
    >>> output = open("myprodocrecord.html", "w")
    >>> output.write(html)
    >>> output.close()

对于这些函数，无效的编号会返回一个HTML格式的错误信息。

10.6  浏览Prosite数据库
-----------------------------------

`ScanProsite <http://www.expasy.org/tools/scanprosite/>`__  允许你通过向Prosite数据库提供一个
Uniprot或者PDB序列编号或序列来在线浏览蛋白序列。关于ScanProsite更多的信息，请阅
读 `ScanProsite文档 <http://www.expasy.org/tools/scanprosite/scanprosite-doc.html>`__ 以及
`程序性访问ScanProsite说明文档 <http://www.expasy.org/tools/scanprosite/ScanPrositeREST.html>`__ 。

你也可以使用Biopython的 ``Bio.ExPASy.ScanProsite`` 模块来从python浏览Prosite数据库，这个模块既
能够帮你安全访问ScanProsite，也可以对ScanProsite返回的结果进行解析。为了查看下边序列中
的Prosite模式（pattern）：

.. code:: python

    MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFT
    CRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN

你可以使用下边的代码：

.. code:: python

    >>> sequence = "MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFT
    CRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN"
    >>> from Bio.ExPASy import ScanProsite
    >>> handle = ScanProsite.scan(seq=sequence)

你可以通过执行 ``handle.read()`` 获取原始XML格式的搜索结果。此外，我们可以使用 ``Bio.ExPASy.ScanProsite.read``
来将原始的XML数据解析到一个python对象：

.. code:: python

    >>> result = ScanProsite.read(handle)
    >>> type(result)
    <class 'Bio.ExPASy.ScanProsite.Record'>

 ``Bio.ExPASy.ScanProsite.Record`` 对象源自一个由ScanProsite返回的包含了ScanProsite hits的列表，这个对象也能够存储hits的数量以及所找到序列的数量。本次ScanProsite搜索找到了6个hits：

.. code:: python

    >>> result.n_seq
    1
    >>> result.n_match
    6
    >>> len(result)
    6
    >>> result[0]
    {'signature_ac': u'PS50948', 'level': u'0', 'stop': 98, 'sequence_ac': u'USERSEQ1', 'start': 16, 'score': u'8.873'}
    >>> result[1]
    {'start': 37, 'stop': 39, 'sequence_ac': u'USERSEQ1', 'signature_ac': u'PS00005'}
    >>> result[2]
    {'start': 45, 'stop': 48, 'sequence_ac': u'USERSEQ1', 'signature_ac': u'PS00006'}
    >>> result[3]
    {'start': 60, 'stop': 62, 'sequence_ac': u'USERSEQ1', 'signature_ac': u'PS00005'}
    >>> result[4]
    {'start': 80, 'stop': 83, 'sequence_ac': u'USERSEQ1', 'signature_ac': u'PS00004'}
    >>> result[5]
    {'start': 106, 'stop': 111, 'sequence_ac': u'USERSEQ1', 'signature_ac': u'PS00008'}

其他的ScanProsite参数可以以关键词参数的形式被传递，更多的信息详见 `程序性访问
ScanProsite说明文档 <http://www.expasy.org/tools/scanprosite/ScanPrositeREST.html>`__ 。
比如，传递 ``lowscore=1`` 可以帮我们找到一个新的低分值hit：

.. code:: python

    >>> handle = ScanProsite.scan(seq=sequence, lowscore=1)
    >>> result = ScanProsite.read(handle)
    >>> result.n_match
    7

