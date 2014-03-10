.. _chapter-entrez:

第9章  访问NCBI Entrez数据库
============================================

Entrez
(`http://www.ncbi.nlm.nih.gov/Entrez <http://www.ncbi.nlm.nih.gov/Entrez>`__)
是一个给客户提供NCBI各个数据库（如PubMed, GeneBank, GEO等等）访问的检索系统。
用户可以通过浏览器手动输入查询条目访问Entrez，也可以使用Biopython的 ``Bio.Entrez`` 模块以编程方式访问来访问Entrez。
如果使用第二种方法，用户用一个Python脚本就可以实现在PubMed里面搜索或者从GenBank下载数据。

``Bio.Entrez`` 模块利用了Entrez Programming Utilities（也称作EUtils），包含八个工具，详情请见NCBI的网站：
`http://www.ncbi.nlm.nih.gov/entrez/utils/ <http://www.ncbi.nlm.nih.gov/entrez/utils/>`__.
每个工具都能在Python的 ``Bio.Entrez`` 模块中找到对应函数，后面会详细讲到。这个模块可以保证用来查询的URL
的正确性，并且向NCBI要求的一样，每三秒钟查询的次数不超过一。

EUtils返回的输出结果通常是XML格式，我们有以下不同的方法来解析这种类型的输入文件：

#. 使用 ``Bio.Entrez``\ 解析器将XML输出的解析成Python对象;
#. 使用Python标准库中的DOM (Document Object Model)解析器;
#. 使用Python标准库中的SAX (Simple API for XML)解析器;
#. 把XML输出当做原始的文本文件，通过字符串查找和处理来进行解析；

对于DOM和SAX解析器，可以查看Python的文档. ``Bio.Entrez`` 中使用到的解析器将会在下面讨论.

NCBI使用DTD (Document Type Definition)文件来描述XML文件中所包含信息的结构. 大多数NCBI使用的DTD文件
格式都包含在了Biopython发行包里。当NCBI Entrez读入一个XML格式的文件的时候，``Bio.Entrez``
将会使用DTD文件。

有时候，你可能会发现与某种特殊的XML相关的DTD文件在Biopython发行包里面不存在。当NCBI升级它的
DTD文件的时候，这种情况可能发生。如果发生这种情况，``Entrez.read`` 将会显示丢失的DTD文件名字和URL的
警示信息。解析器会通过互联网获取缺失的DTD文件，让XML的分析继续正常进行。如果本地存在对应的DTD文件的
话，处理起来会更快。因此，为了更快的处理，我们可以通过警示信息里面的URL来下载对应的DTD文件，将文件放在DTD
文件默认存放的文件夹 ``...site-packages/Bio/Entrez/DTDs`` 。如果你没有权限进入这个文件夹，你也可以把
DTD文件放到 ``~/.biopython/Bio/Entrez/DTDs`` 这个目录，``~`` 表示的是你的Home目录。因为这个目录会先于
``...site-packages/Bio/Entrez/DTDs`` 被解析器读取，所以当 ``...site-packages/Bio/Entrez/DTDs`` 
下面的DTD文件过时的时候，你也可以将最新版本的DTD文件放到Home目录的那个文件夹下面。当然也有其他方案，如果你
是通过源码来安装的Biopython，你可以将DTD文件放到源码的 ``Bio/Entrez/DTDs`` 文件夹下，然后重新安装Biopython。
这样会将新的DTD文件和之前的一样地安装到正确的位置。

Entrez Programming Utilities也可以生成其他格式的输出文件，比如Fasta、序列数据库里面的GenBank文件格式
或者文献数据库里面的MedLine格式，更多内容将会在章节 :ref:`9.12 <sec-entrez-specialized-parsers>` 中讨论。

.. _sec-entrez-guidelines:

9.1  Entrez 简介
----------------------

在我们通过Biopython访问NCBI的线上资源（通过 ``Bio.Entrez`` 或者其他模块）的时候，请先阅读 `NCBI的Entrez
用户规范 <http://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen>`__.
如果NCBI发现你在滥用他们的系统，他们会禁止你的访问。

详细规范如下：

-  对任何连续超过100次的访问请求，请在周末时间或者避开美国的使用高峰时间。这个取决于你是否遵从。
-  使用这个网址 `http://eutils.ncbi.nlm.nih.gov <http://eutils.ncbi.nlm.nih.gov>`__ ，
   而不是通常的NCBI网址。Biopython使用的是这个网址。
-  每秒钟不要超过三次请求（比2009年年初的每三秒钟最多一次请求要宽松）。这个由Biopython自动强制实行。
-  使用email参数，这样如果遇到什么问题，NCBI可以通过邮件联系到你。你可以在每次请求Entrez的时候明确的设置
   这个参数（例如，在参数列表中包含 ``email="A.N.Other@example.com"`` ），或者你也可以设置一个全局的email
   地址：

   .. code:: python

       >>> from Bio import Entrez
       >>> Entrez.email = "A.N.Other@example.com"

   ``Bio.Entrez`` 将会在每次向Entrez请求的时候使用这个邮件地址。请千万不要胡乱的填写邮件地址，不填写都比
   这要好。邮件的参数从2010年6月1日将是强制的参数。在过度使用的情况下，NCBI会在封锁用户访问E-utilities之前尝试通过
   用户提供的邮件地址联系。

-  如果你是在一个大的软件包里面使用Biopython的，请通过tool这个参数明确说明。你既可以在每次请求访问Entrez
   的时候通过参数明确地指明使用的工具（例如，在参数列表中包含 ``tool="MyLocalScript"`` ），或者你也可以
   设置一个全局的tool名称：

   .. code:: python

       >>> from Bio import Entrez
       >>> Entrez.tool = "MyLocalScript"

   默认的tool名称是Biopython。

-  对于大规模的查询请求，NCBI也推荐使用他们的会话历史特性（ WebEnv会话cookie字符串，见
   章节 :ref:`9.15 <sec-entrez-webenv>` ）。 只是这个稍微有点复杂。
   

最后，根据你的使用情况选择不同的策略。如果你打算下载大量的数据，最好使用其他的方法。比如，你想得到所有人的
基因的数据，那么考虑通过FTP得到每个染色体的GenBank文件，然后将这些文件导入到你自己的BioSQL数据库里面去。
(请见章节 :ref:`18.5 <sec-BioSQL>` ).

.. _sec-entrez-einfo:

9.2  EInfo: 获取Entrez数据库的信息
------------------------------------------------------------

EInfo为每个NCBI的数据库提供了条目索引，最近更新的时间以及可用的链接。此外，你可以很容易的使用EInfo通过
Entrez获取所有数据库名字的列表：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.einfo()
    >>> result = handle.read()
    
变量 ``result`` 现在包含了XML格式的数据库列表：

.. code:: python

    >>> print result
    <?xml version="1.0"?>
    <!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD eInfoResult, 11 May 2002//EN"
     "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eInfo_020511.dtd">
    <eInfoResult>
    <DbList>
            <DbName>pubmed</DbName>
            <DbName>protein</DbName>
            <DbName>nucleotide</DbName>
            <DbName>nuccore</DbName>
            <DbName>nucgss</DbName>
            <DbName>nucest</DbName>
            <DbName>structure</DbName>
            <DbName>genome</DbName>
            <DbName>books</DbName>
            <DbName>cancerchromosomes</DbName>
            <DbName>cdd</DbName>
            <DbName>gap</DbName>
            <DbName>domains</DbName>
            <DbName>gene</DbName>
            <DbName>genomeprj</DbName>
            <DbName>gensat</DbName>
            <DbName>geo</DbName>
            <DbName>gds</DbName>
            <DbName>homologene</DbName>
            <DbName>journals</DbName>
            <DbName>mesh</DbName>
            <DbName>ncbisearch</DbName>
            <DbName>nlmcatalog</DbName>
            <DbName>omia</DbName>
            <DbName>omim</DbName>
            <DbName>pmc</DbName>
            <DbName>popset</DbName>
            <DbName>probe</DbName>
            <DbName>proteinclusters</DbName>
            <DbName>pcassay</DbName>
            <DbName>pccompound</DbName>
            <DbName>pcsubstance</DbName>
            <DbName>snp</DbName>
            <DbName>taxonomy</DbName>
            <DbName>toolkit</DbName>
            <DbName>unigene</DbName>
            <DbName>unists</DbName>
    </DbList>
    </eInfoResult>

因为这是一个相当简单的XML文件，我们可以简单的通过字符串查找提取里面所包含的信息。使用 ``Bio.Entrez`` 的解析器，
我们可以直接将这个XML读入到一个Python对象里面去：

.. code:: python

    >>> from Bio import Entrez
    >>> handle = Entrez.einfo()
    >>> record = Entrez.read(handle)

现在 ``record`` 是拥有一个确定键值的字典：

.. code:: python

    >>> record.keys()
    [u'DbList']

这个键对应的值存储了上面XML文件里面包含的数据库名字的列表：

.. code:: python

    >>> record["DbList"]
    ['pubmed', 'protein', 'nucleotide', 'nuccore', 'nucgss', 'nucest',
     'structure', 'genome', 'books', 'cancerchromosomes', 'cdd', 'gap',
     'domains', 'gene', 'genomeprj', 'gensat', 'geo', 'gds', 'homologene',
     'journals', 'mesh', 'ncbisearch', 'nlmcatalog', 'omia', 'omim', 'pmc',
     'popset', 'probe', 'proteinclusters', 'pcassay', 'pccompound',
     'pcsubstance', 'snp', 'taxonomy', 'toolkit', 'unigene', 'unists']

对于这些数据库，我们可以使用EInfo获得更多的信息：

.. code:: python

    >>> handle = Entrez.einfo(db="pubmed")
    >>> record = Entrez.read(handle)
    >>> record["DbInfo"]["Description"]
    'PubMed bibliographic record'
    >>> record["DbInfo"]["Count"]
    '17989604'
    >>> record["DbInfo"]["LastUpdate"]
    '2008/05/24 06:45'

通过 ``record["DbInfo"].keys()`` 可以获取存储在这个记录里面的其他信息。这里面最有用的信息之一是一个ESearch可用的
搜索值列表：

.. code:: python

    >>> for field in record["DbInfo"]["FieldList"]:
    ...     print "%(Name)s, %(FullName)s, %(Description)s" % field
    ALL, All Fields, All terms from all searchable fields
    UID, UID, Unique number assigned to publication
    FILT, Filter, Limits the records
    TITL, Title, Words in title of publication
    WORD, Text Word, Free text associated with publication
    MESH, MeSH Terms, Medical Subject Headings assigned to publication
    MAJR, MeSH Major Topic, MeSH terms of major importance to publication
    AUTH, Author, Author(s) of publication
    JOUR, Journal, Journal abbreviation of publication
    AFFL, Affiliation, Author's institutional affiliation and address
    ...

这是一个很长的列表，但是间接的告诉你在使用PubMed的时候，你可以通过 ``Jones[AUTH]`` 搜索作者，或者通过
``Sanger[AFFL]`` 将作者范围限制在Sanger Centre。这个会非常方便，特别是在你对某个数据库不太熟悉的时候。

9.3  ESearch: 搜索Entrez数据库
--------------------------------------------

我们可以使用 ``Bio.Entrez.esearch()`` 来搜索任意的数据库。例如，我们在PubMed中搜索跟Biopython相关的文献：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.esearch(db="pubmed", term="biopython")
    >>> record = Entrez.read(handle)
    >>> record["IdList"]
    ['19304878', '18606172', '16403221', '16377612', '14871861', '14630660', '12230038']

在输出的结果中，我们可以看到七个PubMed IDs（包括19304878，这个是Biopython应用笔记的PMID），你可以通过
EFetch来获取这些文献（请见章节 :ref:`9.6 <sec-efetch>` ）。

你也可以通过ESearch来搜索GenBank。我们将以在*Cypripedioideae* orchids中搜索*matK*基因为例，快速展示
一下（请见章节 :ref:`9.2 <sec-entrez-einfo>` 关于EInfo：一种查明你可以在哪个Entrez数据库中搜索的方法）。

.. code:: python

    >>> handle = Entrez.esearch(db="nucleotide",term="Cypripedioideae[Orgn] AND matK[Gene]")
    >>> record = Entrez.read(handle)
    >>> record["Count"]
    '25'
    >>> record["IdList"]
    ['126789333', '37222967', '37222966', '37222965', ..., '61585492']

每个IDs(126789333, 37222967, 37222966, …)是GenBank的一个标识。请见章节 :ref:`9.6 <sec-efetch>`
此章包含了怎样下载这些GenBank的记录的信息。

注意，不是像 ``Cypripedioideae[Orgn]`` 这样在搜索的时候加上特定的物种名字，而是需要在搜索的时候使用NCBI的
taxon ID，像 ``txid158330[Orgn]`` 这样。这个并没有记录在ESearch的帮助页面上，NCBI通过邮件回复解释了这个
问题。你可以通过经常和Entrez的网站接口互动，来推断搜索条目的格式。例如，在基因组搜索的时候加上 ``complete[prop]`` 
可以把结果限制在完成的基因组上。

作为最后一个例子，让我们获取一个computational journal名字的列表：

.. code:: python

    >>> handle = Entrez.esearch(db="journals", term="computational")
    >>> record = Entrez.read(handle)
    >>> record["Count"]
    '16'
    >>> record["IdList"]
    ['30367', '33843', '33823', '32989', '33190', '33009', '31986',
     '34502', '8799', '22857', '32675', '20258', '33859', '32534',
     '32357', '32249']

同样，我们可以通过EFetch来获得关于每个journal IDs更多的消息。

ESearch有很多有用的参数——参见 `ESearch 帮助页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html>`__
来获取更多信息.

9.4  EPost: 上传identifiers的列表
-------------------------------------------

EPost上传在后续搜索中将会用到的IDs的列表，参见 `EPost 帮助页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/epost_help.html>`__
来获取更多信息. 通过 ``Bio.Entrez.epost()`` 函数可以在Biopython中实现。

为了举一个关于此用法的例子，假设你有一个想通过EFetch下载的IDs的长长的列表（可能是序列，也有可能是引用的
其他内容）。当你通过EFetch发出下载请求的时候，你的IDs列表、数据库等，将会被转变成一个长的URL，然后被发送
到服务器。如果IDs列表很长，URL也会很长，长的URL可能会断掉（比如，一些代理不能复制全部的内容）。

另外，你也可以把以上分成两步来完成，首先用EPost来上传IDs的列表（这个使用了一个内部的 “HTML post” ，而不是
“HTML get” ， 避开了long URL可能产生的问题）。由于历史记录的支持，你可以使用EFetch来指向这个长的IDs列表，
并且下载相关的数据。

让我们通过下面一个简单的例子来看看EPost是如何工作的——上传了一些PubMed的IDs：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> id_list = ["19304878", "18606172", "16403221", "16377612", "14871861", "14630660"]
    >>> print Entrez.epost("pubmed", id=",".join(id_list)).read()
    <?xml version="1.0"?>
    <!DOCTYPE ePostResult PUBLIC "-//NLM//DTD ePostResult, 11 May 2002//EN"
     "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/ePost_020511.dtd">
    <ePostResult>
     <QueryKey>1</QueryKey>
     <WebEnv>NCID_01_206841095_130.14.22.101_9001_1242061629</WebEnv>
    </ePostResult>

返回的XML包含了两个重要的字符串， ``QueryKey`` 和 ``WebEnv`` ，两个字符串一起确定了之前的历史记录。你可以
使用其他的Entrez工具，例如EFetch，来提取这些值：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> id_list = ["19304878", "18606172", "16403221", "16377612", "14871861", "14630660"]
    >>> search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(id_list)))
    >>> webenv = search_results["WebEnv"]
    >>> query_key = search_results["QueryKey"] 

第 :ref:`9.15 <sec-entrez-webenv>` 章节讲述了如何使用历史的特性。

9.5  ESummary: 通过主要的IDs来获取摘要
----------------------------------------------------

ESummary可以通过一个primary IDs来获取文章的摘要（参见 `ESummary 帮助页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/esummary_help.html>`__
来获取更多信息）。在Biopython中，ESummary以 ``Bio.Entrez.esummary()`` 的形式出现。根据上面的搜索结果，
我们可以获得ID为30367杂志相关的更多信息：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.esummary(db="journals", id="30367")
    >>> record = Entrez.read(handle)
    >>> record[0]["Id"]
    '30367'
    >>> record[0]["Title"]
    'Computational biology and chemistry'
    >>> record[0]["Publisher"]
    'Pergamon,'

.. _sec-efetch:

9.6  EFetch: 从Entrez下载更多的记录
-------------------------------------------------

当你想要从Entrez中提取完整的记录的时候，你可以使用EFetch。 在 `EFetch的帮助页面 <http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html>`__
可以查到EFetch可以起作用的数据库。

NCBI大部分的数据库都支持多种不同的文件格式。当使用 ``Bio.Entrez.efetch()`` 从Entrez下载特定的某种格式的时候，
需要 ``rettype`` 和或者 ``retmode`` 这些可选的参数。对于不同数据库类型不同的搭配在下面的网页中有描述：
`NCBI efetch webpage <http://www.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html>`__
(例如：
`literature <http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetchlit_help.html>`__,
`sequences <http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetchseq_help.html>`__
and
`taxonomy <http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetchtax_help.html>`__).

一种常用的用法是下载FASTA或者GenBank/GenPept的文本格式 (接着可以使用 ``Bio.SeqIO`` 来解析, 参见 :ref:`5.3.1 <sec-SeqIO_GenBank_Online>`
和 :ref:`9.6 <sec-efetch>` ）。从上面 *Cypripedioideae* 的例子,我们可以通过 ``Bio.Entrez.efetch`` 
从GenBank下载记录186972394。

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
    >>> print handle.read()
    LOCUS       EU490707                1302 bp    DNA     linear   PLN 05-MAY-2008
    DEFINITION  Selenipedium aequinoctiale maturase K (matK) gene, partial cds;
                chloroplast.
    ACCESSION   EU490707
    VERSION     EU490707.1  GI:186972394
    KEYWORDS    .
    SOURCE      chloroplast Selenipedium aequinoctiale
      ORGANISM  Selenipedium aequinoctiale
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliophyta; Liliopsida; Asparagales; Orchidaceae;
                Cypripedioideae; Selenipedium.
    REFERENCE   1  (bases 1 to 1302)
      AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A.,
                Endara,C.L., Williams,N.H. and Moore,M.J.
      TITLE     Phylogenetic utility of ycf1 in orchids
      JOURNAL   Unpublished
    REFERENCE   2  (bases 1 to 1302)
      AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A.,
                Endara,C.L., Williams,N.H. and Moore,M.J.
      TITLE     Direct Submission
      JOURNAL   Submitted (14-FEB-2008) Department of Botany, University of
                Florida, 220 Bartram Hall, Gainesville, FL 32611-8526, USA
    FEATURES             Location/Qualifiers
         source          1..1302
                         /organism="Selenipedium aequinoctiale"
                         /organelle="plastid:chloroplast"
                         /mol_type="genomic DNA"
                         /specimen_voucher="FLAS:Blanco 2475"
                         /db_xref="taxon:256374"
         gene            <1..>1302
                         /gene="matK"
         CDS             <1..>1302
                         /gene="matK"
                         /codon_start=1
                         /transl_table=11
                         /product="maturase K"
                         /protein_id="ACC99456.1"
                         /db_xref="GI:186972395"
                         /translation="IFYEPVEIFGYDNKSSLVLVKRLITRMYQQNFLISSVNDSNQKG
                         FWGHKHFFSSHFSSQMVSEGFGVILEIPFSSQLVSSLEEKKIPKYQNLRSIHSIFPFL
                         EDKFLHLNYVSDLLIPHPIHLEILVQILQCRIKDVPSLHLLRLLFHEYHNLNSLITSK
                         KFIYAFSKRKKRFLWLLYNSYVYECEYLFQFLRKQSSYLRSTSSGVFLERTHLYVKIE
                         HLLVVCCNSFQRILCFLKDPFMHYVRYQGKAILASKGTLILMKKWKFHLVNFWQSYFH
                         FWSQPYRIHIKQLSNYSFSFLGYFSSVLENHLVVRNQMLENSFIINLLTKKFDTIAPV
                         ISLIGSLSKAQFCTVLGHPISKPIWTDFSDSDILDRFCRICRNLCRYHSGSSKKQVLY
                         RIKYILRLSCARTLARKHKSTVRTFMRRLGSGLLEEFFMEEE"
    ORIGIN      
            1 attttttacg aacctgtgga aatttttggt tatgacaata aatctagttt agtacttgtg
           61 aaacgtttaa ttactcgaat gtatcaacag aattttttga tttcttcggt taatgattct
          121 aaccaaaaag gattttgggg gcacaagcat tttttttctt ctcatttttc ttctcaaatg
          181 gtatcagaag gttttggagt cattctggaa attccattct cgtcgcaatt agtatcttct
          241 cttgaagaaa aaaaaatacc aaaatatcag aatttacgat ctattcattc aatatttccc
          301 tttttagaag acaaattttt acatttgaat tatgtgtcag atctactaat accccatccc
          361 atccatctgg aaatcttggt tcaaatcctt caatgccgga tcaaggatgt tccttctttg
          421 catttattgc gattgctttt ccacgaatat cataatttga atagtctcat tacttcaaag
          481 aaattcattt acgccttttc aaaaagaaag aaaagattcc tttggttact atataattct
          541 tatgtatatg aatgcgaata tctattccag tttcttcgta aacagtcttc ttatttacga
          601 tcaacatctt ctggagtctt tcttgagcga acacatttat atgtaaaaat agaacatctt
          661 ctagtagtgt gttgtaattc ttttcagagg atcctatgct ttctcaagga tcctttcatg
          721 cattatgttc gatatcaagg aaaagcaatt ctggcttcaa agggaactct tattctgatg
          781 aagaaatgga aatttcatct tgtgaatttt tggcaatctt attttcactt ttggtctcaa
          841 ccgtatagga ttcatataaa gcaattatcc aactattcct tctcttttct ggggtatttt
          901 tcaagtgtac tagaaaatca tttggtagta agaaatcaaa tgctagagaa ttcatttata
          961 ataaatcttc tgactaagaa attcgatacc atagccccag ttatttctct tattggatca
         1021 ttgtcgaaag ctcaattttg tactgtattg ggtcatccta ttagtaaacc gatctggacc
         1081 gatttctcgg attctgatat tcttgatcga ttttgccgga tatgtagaaa tctttgtcgt
         1141 tatcacagcg gatcctcaaa aaaacaggtt ttgtatcgta taaaatatat acttcgactt
         1201 tcgtgtgcta gaactttggc acggaaacat aaaagtacag tacgcacttt tatgcgaaga
         1261 ttaggttcgg gattattaga agaattcttt atggaagaag aa
    //

参数 ``rettype="gb"`` 和 ``retmode="text"`` 让我们下载的数据为GenBank格式。

需要注意的是直到2009年，Entrez EFetch API要求使用 “genbank” 作为返回类型，然而现在NCBI坚持使用官方的
“gb” 或 “gbwithparts” （或者针对蛋白的“gp”) 返回类型。同样需要注意的是，直到2012年2月，
Entrez EFetch API默认的返回格式为纯文本格式文件，现在默认的为XML格式。

作为另外的选择，你也可以使用 ``rettype="fasta"`` 来获取Fasta格式的文件；参见 `EFetch Sequences 帮助页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/efetchseq_help.html>`__ 。
记住，可选的数据格式决定于你要下载的数据库——请参见 `EFetch 帮助页面 <http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html>`__.

如果你要获取记录的格式是 ``Bio.SeqIO`` 所接受的一种格式(见第 :ref:`5 <chapter-Bio.SeqIO>` 章),
你可以直接将其解析为一个 ``SeqRecord`` ：

.. code:: python

    >>> from Bio import Entrez, SeqIO
    >>> handle = Entrez.efetch(db="nucleotide", id="186972394",rettype="gb", retmode="text")
    >>> record = SeqIO.read(handle, "genbank")
    >>> handle.close()
    >>> print record
    ID: EU490707.1
    Name: EU490707
    Description: Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast.
    Number of features: 3
    ...
    Seq('ATTTTTTACGAACCTGTGGAAATTTTTGGTTATGACAATAAATCTAGTTTAGTA...GAA', IUPACAmbiguousDNA())

需要注意的是，一种更加典型的用法是先把序列数据保存到一个本地文件，*然后* 使用 ``Bio.SeqIO`` 来解析。这样就避免了
在运行脚本的时候需要重复的下载同样的文件，并减轻NCBI服务器的负载。例如：

.. code:: python

    import os
    from Bio import SeqIO
    from Bio import Entrez
    Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    filename = "gi_186972394.gbk"
    if not os.path.isfile(filename):
        # Downloading...
        net_handle = Entrez.efetch(db="nucleotide",id="186972394",rettype="gb", retmode="text")
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print "Saved"

    print "Parsing..."
    record = SeqIO.read(filename, "genbank")
    print record

为了得到XML格式的输出，你可以使用 ``Bio.Entrez.read()`` 函数和参数 ``retmode="xml"`` 进行解析，：

.. code:: python

    >>> from Bio import Entrez
    >>> handle = Entrez.efetch(db="nucleotide", id="186972394", retmode="xml")
    >>> record = Entrez.read(handle)
    >>> handle.close()
    >>> record[0]["GBSeq_definition"] 
    'Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast'
    >>> record[0]["GBSeq_source"] 
    'chloroplast Selenipedium aequinoctiale'

就像这样处理数据。例如解析其他数据库特异的文件格式（例如，PubMed中用到的 ``MEDLINE`` 格式），请参见章节 :ref:`9.12 <sec-entrez-specialized-parsers>` .

如果你想使用 ``Bio.Entrez.esearch()`` 进行搜索，然后用 ``Bio.Entrez.efetch()`` 下载数据，那么你需要用到
WebEnv的历史特性，请参加见章节 :ref:`9.15 <sec-entrez-webenv>` .

.. _sec-elink:

9.7  ELink: 在NCBI Entrez中搜索相关的条目
------------------------------------------------------

ELink，在Biopython中是 ``Bio.Entrez.elink()`` ，可以用来在NCBI Entrez数据库中寻找相关的条目。例如，你
可以使用它在gene数据库中寻找核苷酸条目，或者其他很酷的事情。

让我们使用ELink来在2009年的 *Bioinformatics* 杂志中寻找与Biopython应用相关的文章。这篇文章的PubMed ID
是19304878：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"
    >>> pmid = "19304878"
    >>> record = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid))

变量 ``record`` 包含了一个Python列表，列出了已经搜索过的数据库。因为我们特指了一个PubMed ID来搜索，所以
``record`` 只包含了一个条目。这个条目是一个字典变量，包含了我们需要寻找的条目的信息，以及能搜索到的所有相关
的内容：

.. code:: python

    >>> record[0]["DbFrom"]
    'pubmed'
    >>> record[0]["IdList"]
    ['19304878']

键 ``"LinkSetDb"`` 包含了搜索结果，将每个目标数据库保存为一个列表。在我们这个搜索中，我们只在PubMed数据库
中找到了结果（尽管已经被分到了不同的分类）：

.. code:: python

    >>> len(record[0]["LinkSetDb"])
    5
    >>> for linksetdb in record[0]["LinkSetDb"]:
    ...     print linksetdb["DbTo"], linksetdb["LinkName"], len(linksetdb["Link"])
    ... 
    pubmed pubmed_pubmed 110
    pubmed pubmed_pubmed_combined 6
    pubmed pubmed_pubmed_five 6
    pubmed pubmed_pubmed_reviews 5
    pubmed pubmed_pubmed_reviews_five 5

实际的搜索结果被保存在键值为 ``"Link"`` 的字典下。在标准搜索下，总共找到了110个条目。让我们现在看看我们第一个
搜索结果：

.. code:: python

    >>> record[0]["LinkSetDb"][0]["Link"][0]
    {u'Id': '19304878'}

这个就是我们搜索的文章，从中并不能看到更多的结果，所以让我们来看看我们的第二个搜索结果：

.. code:: python

    >>> record[0]["LinkSetDb"][0]["Link"][1]
    {u'Id': '14630660'}

这个PubMed ID为14530660的文章是关于Biopython PDB解析器的。

我们通过一个循环来打印出所有的PubMed IDs：

.. code:: python

    >>> for link in record[0]["LinkSetDb"][0]["Link"] : print link["Id"]
    19304878
    14630660
    18689808
    17121776
    16377612
    12368254
    ......

现在漂亮极了，但是对我个人而言，我对某篇文章是否被引用过更感兴趣。好吧，ELink也可以完成这个——至少对PubMed
Central的杂志来说是这样的（请见章节 :ref:`9.15.3 <sec-elink-citations>` ）。

关于ELink的帮助，请见 `ELink 帮助页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/elink_help.html>`__ 。
这是一个关于 `link names <http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/entrezlinks.html>`__
的整个的子页面， 描述了不同的数据库可以怎样交叉的索引。

9.8  EGQuery: 全局搜索- 统计搜索的条目
----------------------------------------------------

EGQuery提供搜索字段在每个Entrez数据库中的数目。当我们只需要知道在每个数据库中能找到的条目的个数，
而不需要知道具体搜索结果的时候，这个非常的有用（请见例子 :ref:`9.14.2 <sec-entrez_example_genbank>` ）。

在这个例子中，我们使用 ``Bio.Entrez.egquery()`` 来获取跟 “Biopython” 相关的数目：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.egquery(term="biopython")
    >>> record = Entrez.read(handle)
    >>> for row in record["eGQueryResult"]: print row["DbName"], row["Count"]
    ...
    pubmed 6
    pmc 62
    journals 0
    ...

请见 `EGQuery 帮助页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/egquery_help.html>`__
获得更多信息.

9.9  ESpell: 获得拼写建议
-------------------------------------------

ESpell可以检索拼写建议。在这个例子中，我们使用 ``Bio.Entrez.espell()`` 来获得Biopython正确的拼写：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.espell(term="biopythooon")
    >>> record = Entrez.read(handle)
    >>> record["Query"]
    'biopythooon'
    >>> record["CorrectedQuery"]
    'biopython'

请见 `ESpell 帮助页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/espell_help.html>`__
获得更多信息. 这个的主要用法是在使用GUI工具的时候为搜索的条目自动的提供拼写建议。

9.10  解析大的Entrez XML文件
-----------------------------------

``Entrez.read`` 函数将Entrez返回的结果读取到一个Python对象里面去，这个对象被保存在内存中。对于解析太大的
XML文件而内存不够时，可以使用 ``Entrez.parse`` 这个函数。这是一个生成器函数，它将一个一个的读取XML文件里面的内容。只有XML
文件是一个列表对象的时候，这个函数才有用（换句话说，如果在一个内存无限的计算机上 ``Entrez.read`` 将返回一个
Python列表）。

例如，你可以通过NCBI的FTP站点从Entrez Gene 数据库中下载某个物种全部的条目作为一个文件。这个文件可能很大。
作为一个例子，在2009年9月4日，文件 ``Homo_sapiens.ags.gz`` 包含了Entrez Gene数据库中人的序列，文件大小
有116576kB。这个文件是 ``ASN`` 格式，可以通过NCBI的 ``gene2xml`` 程序转成XML格式（请到NCBI的FTP站点获取
更多的信息）：

.. code:: python

    gene2xml -b T -i Homo_sapiens.ags -o Homo_sapiens.xml

XML结果文件有6.1GB. 在大多数电脑上尝试 ``Entrez.read`` 都会导致 ``MemoryError`` 。

XML文件 ``Homo_sapiens.xml`` 包含了一个Entrez gene记录的列表，每个对应于人的一个Entrez基因信息。 ``Entrez.parse`` 
将一个一个的读取这些记录。这样你可以通过遍历每个记录的方式打印或者存储每个记录相关的信息。例如，下面这个脚本
遍历了Entrez基因里面的记录，打印了每个基因的数目和名字：

.. code:: python

    >>> from Bio import Entrez
    >>> handle = open("Homo_sapiens.xml")
    >>> records = Entrez.parse(handle)

    >>> for record in records:
    ...     status = record['Entrezgene_track-info']['Gene-track']['Gene-track_status']
    ...     if status.attributes['value']=='discontinued':
    ...         continue
    ...     geneid = record['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
    ...     genename = record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
    ...     print geneid, genename

将会打印以下内容:

.. code:: python

    1 A1BG
    2 A2M
    3 A2MP
    8 AA
    9 NAT1
    10 NAT2
    11 AACP
    12 SERPINA3
    13 AADAC
    14 AAMP
    15 AANAT
    16 AARS
    17 AAVS1
    ...

9.11  错误处理
---------------------

当解析XML文件的时候，可能出现一下三个错误：

-  这个文件可能不是以常规的 XML 文件格式开头；
-  这个文件可能不完整或者包含一些非 XML 格式的内容；
-  这个文件是正常的 XML 文件，但是包含和相关 DTD 文件无关的条目。

第一种情况会在，例如，你尝试把一个 Fasta 文件当做 XML 文件来处理时发生：

.. code:: python

    >>> from Bio import Entrez
    >>> handle = open("NC_005816.fna") # a Fasta file
    >>> record = Entrez.read(handle)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/site-packages/Bio/Entrez/__init__.py", line 257, in read
        record = handler.read(handle)
      File "/usr/local/lib/python2.7/site-packages/Bio/Entrez/Parser.py", line 164, in read
        raise NotXMLError(e)
    Bio.Entrez.Parser.NotXMLError: Failed to parse the XML data (syntax error: line 1, column 0). Please make sure that the input data are in XML format.

这时候，解析器找不到 ``<?xml ...`` 标签，而这是一个 XML 文件开始的标志，那么可以确定这个文件不是 XML 文件。

当你的文件是XML格式，但是是不完整的（例如，提前结束了），那么解析器会报CorruptedXMLError错误。下面
这个是一个XML文件提前结束的例子：

.. code:: python

    <?xml version="1.0"?>
    <!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD eInfoResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eInfo_020511.dtd">
    <eInfoResult>
    <DbList>
            <DbName>pubmed</DbName>
            <DbName>protein</DbName>
            <DbName>nucleotide</DbName>
            <DbName>nuccore</DbName>
            <DbName>nucgss</DbName>
            <DbName>nucest</DbName>
            <DbName>structure</DbName>
            <DbName>genome</DbName>
            <DbName>books</DbName>
            <DbName>cancerchromosomes</DbName>
            <DbName>cdd</DbName>

这个会生成以下的日志文件：

.. code:: python

    >>> Entrez.read(handle)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/site-packages/Bio/Entrez/__init__.py", line 257, in read
        record = handler.read(handle)
      File "/usr/local/lib/python2.7/site-packages/Bio/Entrez/Parser.py", line 160, in read
        raise CorruptedXMLError(e)
    Bio.Entrez.Parser.CorruptedXMLError: Failed to parse the XML data (no element found: line 16, column 0). Please make sure that the input data are not corrupted.

    >>>

注意，报错信息告诉你在XML文件的什么位置检测到了错误。

如果XML文件当中包含有对应DTD文件中没有描述的标签的时候，会发生第三类错误。以下是这样一个XML文件的例子：

.. code:: python

    <?xml version="1.0"?>
    <!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD eInfoResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eInfo_020511.dtd">
    <eInfoResult>
            <DbInfo>
            <DbName>pubmed</DbName>
            <MenuName>PubMed</MenuName>
            <Description>PubMed bibliographic record</Description>
            <Count>20161961</Count>
            <LastUpdate>2010/09/10 04:52</LastUpdate>
            <FieldList>
                    <Field>
    ...
                    </Field>
            </FieldList>
            <DocsumList>
                    <Docsum>
                            <DsName>PubDate</DsName>
                            <DsType>4</DsType>
                            <DsTypeName>string</DsTypeName>
                    </Docsum>
                    <Docsum>
                            <DsName>EPubDate</DsName>
    ...
            </DbInfo>
    </eInfoResult>

在这个文件里面，因为一些原因，``<DocsumList>`` （还有一些其他的）标签没有在DTD文件 ``eInfo_020511.dtd`` 
中列出来，XML文件对应DTD文件的第二行会特别的描述出来。默认情况下，如果没有找到DTD文件中的标签，解析器
会中止并报ValidationError错误。

.. code:: python

    >>> from Bio import Entrez
    >>> handle = open("einfo3.xml")
    >>> record = Entrez.read(handle)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/site-packages/Bio/Entrez/__init__.py", line 257, in read
        record = handler.read(handle)
      File "/usr/local/lib/python2.7/site-packages/Bio/Entrez/Parser.py", line 154, in read
        self.parser.ParseFile(handle)
      File "/usr/local/lib/python2.7/site-packages/Bio/Entrez/Parser.py", line 246, in startElementHandler
        raise ValidationError(name)
    Bio.Entrez.Parser.ValidationError: Failed to find tag 'DocsumList' in the DTD. To skip all tags that are not represented in the DTD, please call Bio.Entrez.read or Bio.Entrez.parse with validate=False.

可选地，你可以让解析器跳过这样的标签，而不是报ValidationError错误。通过调用 ``Entrez.read`` 或者
``Entrez.parse`` 并使参数 ``validate`` 等于False可以实现这个功能：

.. code:: python

    >>> from Bio import Entrez
    >>> handle = open("einfo3.xml")
    >>> record = Entrez.read(handle,validate=False)
    >>>

当然，XML文件中的tag没有出现在对应DTD文件中的信息，将不会在 ``Entrez.read`` 的返回记录中出现。

.. _sec-entrez-specialized-parsers:

9.12  专用的解析器
-------------------------

函数 ``Bio.Entrez.read()`` 可以处理大部分（如果不是所有的话）Entrez返回的XML文件。Entrez也可以
允许你通过其他格式来获取数据，有时候，这种方式在可读性上比XML文件格式更具优势（或者下载文件的大小）。

为了使用 ``Bio.Entrez.efetch()`` 函数从Entrez中提取一种特有的文件格式，需要指明 ``rettype`` 和或者或 ``retmode`` 
等可选参数。不同的组合在 `NCBI efetch的页面 <http://www.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html>`__ 。
有对不同数据库的描述。

一个显然的例子是，你可能更想以FASTA或者 GenBank/GenPept ( 这些可以通过 ``Bio.SeqIO`` 来处理, 请见 :ref:`5.3.1 <sec-SeqIO_GenBank_Online>`
和 :ref:`9.6 <sec-efetch>` ） 纯文本形式下载序列。对于文献数据库，Biopython包含了一个处理PubMed中
使用的 ``MEDLINE`` 格式的解析器。

.. _sec-entrez-and-medline:

9.12.1  解析Medline记录
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

你可以在 ``Bio.Medline`` 中找到Medline的解析器。假设你想处理包含一个Medline记录的 ``pubmed_result1.txt`` 
文件。你可以在Biopython的 ``Tests\Medline`` 目录下找到这个文件，这个文件内容如下所示：

.. code:: python

    PMID- 12230038
    OWN - NLM
    STAT- MEDLINE
    DA  - 20020916
    DCOM- 20030606
    LR  - 20041117
    PUBM- Print
    IS  - 1467-5463 (Print)
    VI  - 3
    IP  - 3
    DP  - 2002 Sep
    TI  - The Bio* toolkits--a brief overview.
    PG  - 296-302
    AB  - Bioinformatics research is often difficult to do with commercial software. The
          Open Source BioPerl, BioPython and Biojava projects provide toolkits with
    ...

我们首先打开文件，然后解析它：

.. code:: python

    >>> from Bio import Medline
    >>> input = open("pubmed_result1.txt")
    >>> record = Medline.read(input)


现在 ``record`` 将 Medline记录以Python字典的形式保存起来：

.. code:: python

    >>> record["PMID"]
    '12230038'

.. code:: python

    >>> record["AB"]
    'Bioinformatics research is often difficult to do with commercial software.
    The Open Source BioPerl, BioPython and Biojava projects provide toolkits with
    multiple functionality that make it easier to create customised pipelines or
    analysis. This review briefly compares the quirks of the underlying languages
    and the functionality, documentation, utility and relative advantages of the
    Bio counterparts, particularly from the point of view of the beginning
    biologist programmer.'

用于Medline记录的键值可以相当模糊，使用

.. code:: python

    >>> help(record)

可以做一个简单的总结。

为了解析包含多个Medline记录的文件，你可以使用 ``parse`` 函数来代替：

.. code:: python

    >>> from Bio import Medline
    >>> input = open("pubmed_result2.txt")
    >>> records = Medline.parse(input)
    >>> for record in records:
    ...     print record["TI"]
    A high level interface to SCOP and ASTRAL implemented in python.
    GenomeDiagram: a python package for the visualization of large-scale genomic data.
    Open source clustering software.
    PDB file parser and structure class implemented in Python.

你可以通过 ``Bio.Entrez.efetch`` 来下载Medline记录，而不是保存在某个文件里。例如，让我们来查看PubMed
里面跟Biopython相关的所有所有Medline记录：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.esearch(db="pubmed",term="biopython")
    >>> record = Entrez.read(handle)
    >>> record["IdList"]
    ['19304878', '18606172', '16403221', '16377612', '14871861', '14630660', '12230038']

现在我们使用 ``Bio.Entrez.efetch`` 来下载这些Medline记录:

.. code:: python

    >>> idlist = record["IdList"]
    >>> handle = Entrez.efetch(db="pubmed",id=idlist,rettype="medline",retmode="text")

这里，我们使 ``rettype="medline", retmode="text"`` 来以纯文本形式的Medline格式来得到这些记录。现在
我们使用 ``Bio.Medline`` 来解析这些记录：

.. code:: python

    >>> from Bio import Medline
    >>> records = Medline.parse(handle)
    >>> for record in records:
    ...     print record["AU"]
    ['Cock PJ', 'Antao T', 'Chang JT', 'Chapman BA', 'Cox CJ', 'Dalke A', ..., 'de Hoon MJ']
    ['Munteanu CR', 'Gonzalez-Diaz H', 'Magalhaes AL']
    ['Casbon JA', 'Crooks GE', 'Saqi MA']
    ['Pritchard L', 'White JA', 'Birch PR', 'Toth IK']
    ['de Hoon MJ', 'Imoto S', 'Nolan J', 'Miyano S']
    ['Hamelryck T', 'Manderick B']
    ['Mangalam H']

为了比对，我们展示了一个XML格式的例子：

.. code:: python

    >>> idlist = record["IdList"]
    >>> handle = Entrez.efetch(db="pubmed",id=idlist,rettype="medline",retmode="xml")
    >>> records = Entrez.read(handle)
    >>> for record in records:
    ...     print record["MedlineCitation"]["Article"]["ArticleTitle"]
    Biopython: freely available Python tools for computational molecular biology and
     bioinformatics.
    Enzymes/non-enzymes classification model complexity based on composition, sequence,
     3D and topological indices.
    A high level interface to SCOP and ASTRAL implemented in python.
    GenomeDiagram: a python package for the visualization of large-scale genomic data.
    Open source clustering software.
    PDB file parser and structure class implemented in Python.
    The Bio* toolkits--a brief overview.

需要注意的是，在上面这两个例子当中，为了简便我们混合使用了 ESearch 和 EFetch。在这种情形下，NCBI 希望你
使用他们的历史记录特性，在下面章节中会讲到Section :ref:`9.15 <sec-entrez-webenv>` .

9.12.2  解析GEO记录
~~~~~~~~~~~~~~~~~~~~~~~~~~~

GEO ( `Gene Expression Omnibus <http://www.ncbi.nlm.nih.gov/geo/>`__ ) 是高通量基因表达和杂交芯片
数据的数据库。 ``Bio.Geo`` 模块可以用来解析GEO格式的数据。

下面的代码展示了怎样将一个名称为 ``GSE16.txt`` 的GEO文件存进一个记录，并打印该记录：

.. code:: python

    >>> from Bio import Geo
    >>> handle = open("GSE16.txt")
    >>> records = Geo.parse(handle)
    >>> for record in records:
    ...     print record

你可以使用 ESearch 来搜索 “gds” 数据库 (GEO 数据集) :

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com" # Always tell NCBI who you are
    >>> handle = Entrez.esearch(db="gds",term="GSE16")
    >>> record = Entrez.read(handle)
    >>> record["Count"]
    2
    >>> record["IdList"]
    ['200000016', '100000028']

通过Entrez网站，UID “200000016” 是GDS16，其他的hit “100000028” 是相关的平台。不幸的是，在写
这份指南的时候，NCBI貌似还不支持通过Entrez下载GEO文件（不论XML文件，还是SOFT格式的文件）。

然而，可以相当直接的通过 FTP `ftp://ftp.ncbi.nih.gov/pub/geo/ <ftp://ftp.ncbi.nih.gov/pub/geo/>`__ 来下载 GEO 文件。
在这个例子当中，你需要的文件应该是 `ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/GSE16/GSE16_family.soft.gz <ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/GSE16/GSE16_family.soft.gz>`__
（一个压缩文件，参见Python的gzip 模块）。

9.12.3  解析UniGene记录
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UniGene是NCBI的转录组数据库，每个UniGene记录展示了该转录本在某个特定物种中相关的基因。一个典型的UniGene
记录如下所示：

.. code:: python

    ID          Hs.2
    TITLE       N-acetyltransferase 2 (arylamine N-acetyltransferase)
    GENE        NAT2
    CYTOBAND    8p22
    GENE_ID     10
    LOCUSLINK   10
    HOMOL       YES
    EXPRESS      bone| connective tissue| intestine| liver| liver tumor| normal| soft tissue/muscle tissue tumor| adult
    RESTR_EXPR   adult
    CHROMOSOME  8
    STS         ACC=PMC310725P3 UNISTS=272646
    STS         ACC=WIAF-2120 UNISTS=44576
    STS         ACC=G59899 UNISTS=137181
    ...
    STS         ACC=GDB:187676 UNISTS=155563
    PROTSIM     ORG=10090; PROTGI=6754794; PROTID=NP_035004.1; PCT=76.55; ALN=288
    PROTSIM     ORG=9796; PROTGI=149742490; PROTID=XP_001487907.1; PCT=79.66; ALN=288
    PROTSIM     ORG=9986; PROTGI=126722851; PROTID=NP_001075655.1; PCT=76.90; ALN=288
    ...
    PROTSIM     ORG=9598; PROTGI=114619004; PROTID=XP_519631.2; PCT=98.28; ALN=288

    SCOUNT      38
    SEQUENCE    ACC=BC067218.1; NID=g45501306; PID=g45501307; SEQTYPE=mRNA
    SEQUENCE    ACC=NM_000015.2; NID=g116295259; PID=g116295260; SEQTYPE=mRNA
    SEQUENCE    ACC=D90042.1; NID=g219415; PID=g219416; SEQTYPE=mRNA
    SEQUENCE    ACC=D90040.1; NID=g219411; PID=g219412; SEQTYPE=mRNA
    SEQUENCE    ACC=BC015878.1; NID=g16198419; PID=g16198420; SEQTYPE=mRNA
    SEQUENCE    ACC=CR407631.1; NID=g47115198; PID=g47115199; SEQTYPE=mRNA
    SEQUENCE    ACC=BG569293.1; NID=g13576946; CLONE=IMAGE:4722596; END=5'; LID=6989; SEQTYPE=EST; TRACE=44157214
    ...
    SEQUENCE    ACC=AU099534.1; NID=g13550663; CLONE=HSI08034; END=5'; LID=8800; SEQTYPE=EST
    //

这个记录展示了这个转录本（如 ``SEQUENCE`` 行展示）是来自人的NAT2基因，编码en N-acetyltransferase。
``PROTSIM`` 显示的是和NAT2显著相似的蛋白质， ``STS`` 展示的是基因组当中的STS位点。

我们使用 ``Bio.UniGene`` 模块来解析UniGene文件：

.. code:: python

    >>> from Bio import UniGene
    >>> input = open("myunigenefile.data")
    >>> record = UniGene.read(input)
    
``UniGene.read`` 返回的是一个包含一些和UniGene记录的字段相对应属性的Python对象。例如，

.. code:: python

    >>> record.ID
    "Hs.2"
    >>> record.title
    "N-acetyltransferase 2 (arylamine N-acetyltransferase)"


``EXPRESS`` 和 ``RESTR_EXPR`` 两行被存储为字符串的Python列表：

.. code:: python

    ['bone', 'connective tissue', 'intestine', 'liver', 'liver tumor', 'normal', 'soft tissue/muscle tissue tumor', 'adult']

跟 ``STS`` , ``PROTSIM`` , 和 ``SEQUENCE`` 相关的特有的对象被保存在如下键所对应的字典中：

.. code:: python

    >>> record.sts[0].acc
    'PMC310725P3'
    >>> record.sts[0].unists
    '272646'

和 ``PROTSIM`` 、 ``SEQUENCE`` 这两行相似。

我们使用 ``Bio.UniGene`` 中的 ``parse`` 函数来处理一个文件中包含多个UniGene记录的情况：

.. code:: python

    >>> from Bio import UniGene
    >>> input = open("unigenerecords.data")
    >>> records = UniGene.parse(input)
    >>> for record in records:
    ...     print record.ID

9.13  使用代理
-------------------

通常状况下，你不需要使用代理，但是如果你的网络有问题的时候，我们有以下应对方法。在内部， ``Bio.Entrez`` 使用
一个标准的 Python 库 ``urllib`` 来访问 NCBI的服务器。这个将检查叫做 ``http_proxy`` 的环境变量来自动配置简单
的代理服务。不幸的是，这个模块不支持需要认证的代理。

你可以选择设定环境变量 ``http_proxy`` 。同样，你可以在Python脚本开头的地方设置这个参数，例如：

.. code:: python

    import os
    os.environ["http_proxy"] = "http://proxyhost.example.com:8080"

参见 `urllib
文档 <http://www.python.org/doc/lib/module-urllib.html>`__ 获得更多信息。

9.14  实例
--------------

9.14.1  PubMed和Medline
~~~~~~~~~~~~~~~~~~~~~~~~~~

如果你是在医药领域或者对人类的问题感兴趣（或者尽管并你不感兴趣，大多数情况下也适用!），PubMed(`http://www.ncbi.nlm.nih.gov/PubMed/ <http://www.ncbi.nlm.nih.gov/PubMed/>`__)
是一个包含了各方面的非常优秀的资源。像其他的一样，我们希望能够通过 Python 脚本从中抓取一些信息。

在这个例子当中，我们要查询PubMed当中所有跟Orchids相关的文章(见 :ref:`2.3 <sec-orchids>` 我们的动机)。
我们首先看看有多少这样的文章：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.egquery(term="orchid")
    >>> record = Entrez.read(handle)
    >>> for row in record["eGQueryResult"]:
    ...     if row["DbName"]=="pubmed":
    ...         print row["Count"]
    463

现在我们使用 ``Bio.Entrez.efetch`` 这个函数来下载这463篇文章的PubMed IDs：

.. code:: python

    >>> handle = Entrez.esearch(db="pubmed", term="orchid", retmax=463)
    >>> record = Entrez.read(handle)
    >>> idlist = record["IdList"]
    >>> print idlist

返回值是一个Python列表，包含了所有和orchids相关文章的PubMed IDs：

.. code:: python

    ['18680603', '18665331', '18661158', '18627489', '18627452', '18612381',
    '18594007', '18591784', '18589523', '18579475', '18575811', '18575690',
    ...

这样我们就得到了这些信息，显然我们想要得到对应的Medline records和更多额外的信息。这里，我们将以纯文本的
形式下载和Medline records相关的信息，然后使用 ``Bio.Medline`` 模块来解析他们：

.. code:: python

    >>> from Bio import Medline
    >>> handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",
                               retmode="text")
    >>> records = Medline.parse(handle)

注意 - 我们完成了一次搜索和获取，NCBI更希望你在这种情况下使用他们的历史记录支持。请见章节 :ref:`9.15 <sec-entrez-webenv>` .

请记住 ``records`` 是一个迭代器，所以你只能访问这些records一次。如果你想保存这些records，你需要把他们转成列表：

.. code:: python

    >>> records = list(records)

现在让我们迭代这些records，然后分别打印每一个record的信息：

.. code:: python

    >>> for record in records:
    ...     print "title:", record.get("TI", "?")
    ...     print "authors:", record.get("AU", "?")
    ...     print "source:", record.get("SO", "?")
    ...     print

这个的输出结果是这样的:

.. code:: python

    title: Sex pheromone mimicry in the early spider orchid (ophrys sphegodes):
    patterns of hydrocarbons as the key mechanism for pollination by sexual
    deception [In Process Citation]
    authors: ['Schiestl FP', 'Ayasse M', 'Paulus HF', 'Lofstedt C', 'Hansson BS',
    'Ibarra F', 'Francke W']
    source: J Comp Physiol [A] 2000 Jun;186(6):567-74

特别有意思的是作者的列表，作者的列表会作为一个标准的Python列表返回。这使得用标准的Python工具操作和搜索
变得简单。例如，我们可以像下面的代码这样循环读取所有条目来查找某个特定的作者：

.. code:: python

    >>> search_author = "Waits T"

    >>> for record in records:
    ...     if not "AU" in record:
    ...         continue
    ...     if search_author in record["AU"]:
    ...         print "Author %s found: %s" % (search_author, record["SO"])

希望这个章节可以让你知道Entrez和Medline借口的能力和便利性和怎样同时使用他们。

.. _sec-entrez_example_genbank:

9.14.2  搜索，下载，和解析Entrez核酸记录
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

这里我们将展示一个关于远程Entrez查询的简单例子。在 :ref:`2.3 <sec-orchids>` 节，我们讲到了使用NCBI
的Entrez网站来搜索 NCBI 的核酸数据库来获得关于Cypripedioideae的信息。现在我们看看如何使用Python脚本
自动的处理。在这个例子当中，我们仅仅展示如何使用Entrez模块来连接，获取结果，解析他们。

首先，我们在下载这些结果之前，使用EGQuery来计算结果的数目。EGQuery 将会告诉我们在每个数据库中分别有多少
搜索结果，但在我们这个例子当中，我们只对核苷酸感兴趣：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.egquery(term="Cypripedioideae")
    >>> record = Entrez.read(handle)
    >>> for row in record["eGQueryResult"]:
    ...     if row["DbName"]=="nuccore":
    ...         print row["Count"]
    814

所以，我们预期能找到814个 Entrez 核酸记录（这是我在2008年得到的结果；在未来这个结果应该会增加）。如果你得
到了高的不可思议的结果数目时，你可能得重新考虑是否需要下载所有的这些结果，下载是我们的下一步：

.. code:: python

    >>> from Bio import Entrez
    >>> handle = Entrez.esearch(db="nucleotide", term="Cypripedioideae", retmax=814)
    >>> record = Entrez.read(handle)

在这里, ``record`` 是一个包含了搜索结果和一些辅助信息的Python字典。仅仅作为参考信息，让我们看看在这些字典当中
究竟存储了些什么内容： 

.. code:: python

    >>> print record.keys()
    [u'Count', u'RetMax', u'IdList', u'TranslationSet', u'RetStart', u'QueryTranslation']

首先, 让我们检查看看我们得到了多少个结果:

.. code:: python

    >>> print record["Count"]
    '814'

这个结果是我们所期望的。这814个结果被存在了 ``record['IdList']`` 中:

.. code:: python

    >>> print len(record["IdList"])
    814

让我们看看前五个结果:

.. code:: python

    >>> print record["IdList"][:5]
    ['187237168', '187372713', '187372690', '187372688', '187372686']

我们可以使用 ``efetch`` 来下载这些结果. 尽管你可以一个一个的下载这些记录，但为了减少 NCBI 服务器的负载，最好呢还是
一次性的下载所有的结果。然而在这个情况下，你应该完美的使用在后面章节 :ref:`9.15 <sec-entrez-webenv>` 中会要讲到的历史记录特性。

.. code:: python

    >>> idlist = ",".join(record["IdList"][:5])
    >>> print idlist
    187237168,187372713,187372690,187372688,187372686
    >>> handle = Entrez.efetch(db="nucleotide", id=idlist, retmode="xml")
    >>> records = Entrez.read(handle)
    >>> print len(records)
    5

每个这样的records对应一个GenBank record.

.. code:: python

    >>> print records[0].keys()
    [u'GBSeq_moltype', u'GBSeq_source', u'GBSeq_sequence',
     u'GBSeq_primary-accession', u'GBSeq_definition', u'GBSeq_accession-version',
     u'GBSeq_topology', u'GBSeq_length', u'GBSeq_feature-table',
     u'GBSeq_create-date', u'GBSeq_other-seqids', u'GBSeq_division',
     u'GBSeq_taxonomy', u'GBSeq_references', u'GBSeq_update-date',
     u'GBSeq_organism', u'GBSeq_locus', u'GBSeq_strandedness']

    >>> print records[0]["GBSeq_primary-accession"]
    DQ110336

    >>> print records[0]["GBSeq_other-seqids"]
    ['gb|DQ110336.1|', 'gi|187237168']

    >>> print records[0]["GBSeq_definition"]
    Cypripedium calceolus voucher Davis 03-03 A maturase (matR) gene, partial cds;
    mitochondrial

    >>> print records[0]["GBSeq_organism"]
    Cypripedium calceolus

你可以用这个来快速的开始搜索 —— 但是对于频繁的使用请见 :ref:`9.15 <sec-entrez-webenv>` .

.. _sec-entrez-search-fetch-genbank:

9.14.3  搜索、下载和解析GenBank record
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GenBank record 格式是保存序列信息、序列特征和其他相关信息非常普遍的一种方法。这种格式是从 NCBI 数据库 
`http://www.ncbi.nlm.nih.gov/ <http://www.ncbi.nlm.nih.gov/>`__ 获取信息非常好的一种方式 .

在这个例子当中，我们将展示怎样去查询 NCBI 数据库，根据query提取记录，然后使用 ``Bio.SeqIO`` 解析他们 ——
在 :ref:`5.3.1 <sec-SeqIO_GenBank_Online>` 中提到过这些。简单起见，这个例子*不会*使用 WebEnv 历史记录特性
—— 请到 :ref:`9.15 <sec-entrez-webenv>` 查看。

首先，我们想要查询找出要获取的记录的ID。这里我们快速的检索我们最喜欢的一个物种 *Opuntia* (多刺的梨型仙人掌)。我们
可以做一个快速的检索来获得所有满足要求的GIs（GenBank标志符）。首先我们看看有多少个记录：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.egquery(term="Opuntia AND rpl16")
    >>> record = Entrez.read(handle)
    >>> for row in record["eGQueryResult"]:
    ...     if row["DbName"]=="nuccore":
    ...         print row["Count"]
    ...
    9

现在我们下载GenBank identifiers的列表：

.. code:: python

    >>> handle = Entrez.esearch(db="nuccore", term="Opuntia AND rpl16")
    >>> record = Entrez.read(handle)
    >>> gi_list = record["IdList"]
    >>> gi_list
    ['57240072', '57240071', '6273287', '6273291', '6273290', '6273289', '6273286',
    '6273285', '6273284']

现在我们使用这些GIs来下载GenBank records —— 注意在老的Biopython版本中，你必须将GI号用逗号隔开传递给Entrez，例如
在 Biopython 1.59中，你可以传递一个列表，下面的内容会为你做转换：

.. code:: python

    >>> gi_str = ",".join(gi_list)
    >>> handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")

如果你想看原始的 GenBank 文件，你可以从这个句柄中读取并打印结果：

.. code:: python

    >>> text = handle.read()
    >>> print text
    LOCUS       AY851612                 892 bp    DNA     linear   PLN 10-APR-2007
    DEFINITION  Opuntia subulata rpl16 gene, intron; chloroplast.
    ACCESSION   AY851612
    VERSION     AY851612.1  GI:57240072
    KEYWORDS    .
    SOURCE      chloroplast Austrocylindropuntia subulata
      ORGANISM  Austrocylindropuntia subulata
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons;
                Caryophyllales; Cactaceae; Opuntioideae; Austrocylindropuntia.
    REFERENCE   1  (bases 1 to 892)
      AUTHORS   Butterworth,C.A. and Wallace,R.S.
    ...

在这个例子当中，我们只是得到了原始的记录。为了得到对Python友好的格式，我们可以使用 ``Bio.SeqIO`` 将GenBank
数据转化成 ``SeqRecord`` 对象，包括 ``SeqFeature`` 对象 (请见第 :ref:`5 <chapter-Bio.SeqIO>` 章):

.. code:: python

    >>> from Bio import SeqIO
    >>> handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
    >>> records = SeqIO.parse(handle, "gb")

我们现在可以逐个查看这些record来寻找我们感兴趣的信息：

.. code:: python

    >>> for record in records: 
    >>> ...    print "%s, length %i, with %i features" \
    >>> ...           % (record.name, len(record), len(record.features))
    AY851612, length 892, with 3 features
    AY851611, length 881, with 3 features
    AF191661, length 895, with 3 features
    AF191665, length 902, with 3 features
    AF191664, length 899, with 3 features
    AF191663, length 899, with 3 features
    AF191660, length 893, with 3 features
    AF191659, length 894, with 3 features
    AF191658, length 896, with 3 features

使用这些自动的查询提取功能相对于手动处理是一个很大的进步。尽管这些模块需要遵守NCBI每秒钟最多三次的规则，然而NCBI
有其他像避开高峰时刻的建议。请见章节 :ref:`9.1 <sec-entrez-guidelines>` 。尤其需要注意的是，这个例子没有
用到 WebEnv 历史记录特性。你应该使用这个来完成一些琐碎的搜索和下载的工作，请见章节 :ref:`9.15 <sec-entrez-webenv>` 。

最后，如果你计划重复你的分析，你应该下载这些record *一次* ，然后将他们保存在你的硬盘里，在本地进行分析；而不是
从 NCBI 下载之后就马上进行分析（像这个例子一样）。

9.14.4  查看物种的谱系关系
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

仍然以植物为例子，让我们找出Cyripedioideae兰花家族的谱系。首先让我们在Taxonomy数据库中查找跟Cypripedioideae
相关的记录，确实找到了一个确切的 NCBI taxonomy 标识号：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.esearch(db="Taxonomy", term="Cypripedioideae")
    >>> record = Entrez.read(handle)
    >>> record["IdList"]
    ['158330']
    >>> record["IdList"][0]
    '158330'

现在，我们使用 ``efetch`` 从 Taxonomy 数据库中下载这些条目，然后解析它：

.. code:: python

    >>> handle = Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    >>> records = Entrez.read(handle)

再次，这个record保存了许多的信息：

.. code:: python

    >>> records[0].keys()
    [u'Lineage', u'Division', u'ParentTaxId', u'PubDate', u'LineageEx',
     u'CreateDate', u'TaxId', u'Rank', u'GeneticCode', u'ScientificName',
     u'MitoGeneticCode', u'UpdateDate']

我们可以直接从这个record获得谱系信息：

.. code:: python

    >>> records[0]["Lineage"]
    'cellular organisms; Eukaryota; Viridiplantae; Streptophyta; Streptophytina;
     Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliophyta;
     Liliopsida; Asparagales; Orchidaceae'

这个record数据包含的信息远远超过在这里显示的 —— 例如查看 ``"LineageEx"`` 而不是 ``"Lineage"`` 相关的
信息，你也可以得到谱系里面的 NCBI taxon 标识号信息。

.. _sec-entrez-webenv:

9.15  使用历史记录和WebEnv
----------------------------------

通常，你想做一系列相关的查询。最典型的是，进行一个搜索，精炼搜索，然后提取详细的搜索结果。你 *可以* 通过一系列
独立的调用Entrez来完成这些工作。然而，NCBI更希望你利用历史记录支持的优势来完成这个 - 例如将ESearch
和EFetch结合起来。

另外一个关于历史记录支持，典型的使用是结合EPost和EFetch。你可以使用EPost来上传一个标识号的
列表，这样就开始一些新的history session。接下来你就可以用EFetch指向这个session来下载这些数据（而不是那些
标识号）。

9.15.1  利用 history 来搜索和下载序列
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

假设我们想搜索和下载所有的 *Opuntia* rpl16核酸序列，然后将它们保存到一个FASTA文件里。就像章节 :ref:`9.14.3 <sec-entrez-search-fetch-genbank>`
里一样, 我们可以简单的用 ``Bio.Entrez.esearch()`` 得到一个GI号的列表，然后调用 ``Bio.Entrez.efetch()`` 
来下载他们。

然而，被认同的方法是使用历史记录特性来进行搜索。然后，我们可以通过指向这些搜索结果的引用来获取他们 - 
NCBI 将会提前进行缓冲。

为此，调用 ``Bio.Entrez.esearch()`` 是正常的，但是需要额外的 ``usehistory="y"`` 参数，

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "history.user@example.com"
    >>> search_handle = Entrez.esearch(db="nucleotide",term="Opuntia[orgn] and rpl16",
                                       usehistory="y")
    >>> search_results = Entrez.read(search_handle)
    >>> search_handle.close()

当你得到XML输出的时候，它仍然包括了常见的搜索结果：

.. code:: python

    >>> gi_list = search_results["IdList"]
    >>> count = int(search_results["Count"])
    >>> assert count == len(gi_list)

然而，你将得到两个额外的信息， ``WebEnv`` 会话cookie 和 ``QueryKey`` :

.. code:: python

    >>> webenv = search_results["WebEnv"]
    >>> query_key = search_results["QueryKey"] 

将这些值保存到 ``session_cookie`` 和 ``query_key`` 后，我们可以使用它们作为 ``Bio.Entrez.efetch()`` 的
参数，而不用提供GI numbers的identifiers。

对于小数据量你一次下载所有的数据也没有关系，但是最好能够分批下载。你可以使用 ``restart`` 和 ``retmax`` 来说明
哪一部分搜索结果是你想得到的（条目以0开始计算，返回结果的最大数目）。例如：

.. code:: python

    batch_size = 3
    out_handle = open("orchid_rpl16.fasta", "w")
    for start in range(0,count,batch_size):
        end = min(count, start+batch_size)
        print "Going to download record %i to %i" % (start+1, end)
        fetch_handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text",
                                     retstart=start, retmax=batch_size,
                                     webenv=webenv, query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()

我们以此为例来说明，这个例子分三次来下载FASTA records。除非你是要下载基因组或者染色体数目，你最好选取一个
比较大的batch大小。

9.15.2  利用history来搜索和下载综述
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

这是另外一个history的例子，搜索过去几年当中发表的关于 *Opuntia* 的文章，然后下载到一个MedLine格式的文件里：

.. code:: python

    from Bio import Entrez
    Entrez.email = "history.user@example.com"
    search_results = Entrez.read(Entrez.esearch(db="pubmed",
                                                term="Opuntia[ORGN]",
                                                reldate=365, datetype="pdat",
                                                usehistory="y"))
    count = int(search_results["Count"])
    print "Found %i results" % count

    batch_size = 10
    out_handle = open("recent_orchid_papers.txt", "w")
    for start in range(0,count,batch_size):
        end = min(count, start+batch_size)
        print "Going to download record %i to %i" % (start+1, end)
        fetch_handle = Entrez.efetch(db="pubmed",
                                     rettype="medline", retmode="text",
                                     retstart=start, retmax=batch_size,
                                     webenv=search_results["WebEnv"],
                                     query_key=search_results["QueryKey"])
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()

在写这份文档的时候，这个搜索返回了28个匹配结果 - 但是因为这个是跟时间相关的搜索，因此返回结果会发生变化。
像在上面 :ref:`9.12.1 <sec-entrez-and-medline>` 讲到的一样, 你可以使用 ``Bio.Medline`` 来解析
保存下来的记录。

.. _sec-elink-citations:

9.15.3  搜索引用文章
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

回到 Section :ref:`9.7 <sec-elink>` 我们提到可以使用ELink来搜索制定文章的引用。不幸的是，这个只包含
PubMed Central（为PubMed中所有文献来做这个事情，意味这NIH将要付出更多的工作）包含的那些杂志。让我们以
Biopython PDB parser文章为例来试试看， PubMed ID 14630660：

.. code:: python

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"
    >>> pmid = "14630660"
    >>> results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc",
    ...                                    LinkName="pubmed_pmc_refs", from_uid=pmid))
    >>> pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
    >>> pmc_ids
    ['2744707', '2705363', '2682512', ..., '1190160']

好极了 - 11篇文章。但是为什么没有Biopython应用笔记（PubMed ID 19304878）呢？好吧，你可能已经从变量的
名称中猜到了，实际上他们不是PubMed IDs，而是PubMed Central IDs。我们的应用笔记是列表当中第三个引用的
文章， PMCID 2682512。

那么，如果（像我）你希望得到的是PubMed IDs的列表的话，该怎么做呢？好吧，你可以使用再次使用ELink来更改他们。
这将成为两步处理，所以你应该使用历史记录特性来完成这个工作（章节 :ref:`9.15 <sec-entrez-webenv>` ）。

但是首先，让我们使用更直接的方法来进行第二次调用ELink：

.. code:: python

    >>> results2 = Entrez.read(Entrez.elink(dbfrom="pmc", db="pubmed", LinkName="pmc_pubmed",
    ...                                     from_uid=",".join(pmc_ids)))
    >>> pubmed_ids = [link["Id"] for link in results2[0]["LinkSetDb"][0]["Link"]]
    >>> pubmed_ids
    ['19698094', '19450287', '19304878', ..., '15985178']

这次，你可以立即的看到Biopython应用笔记作为第三个hit（PubMed ID 19304878）。

现在，让我们重新使用历史记录再试一遍 … *TODO*.

最终，不要忘记在Entrez调用的时候，加上你 *自己* 的电子邮箱地址。

