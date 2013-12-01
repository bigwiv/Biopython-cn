第一章 介绍
=======================

1.1  什么是Biopython？
-----------------------

Biopython工程是一个使用Python来开发计算分子生物学工具的国际团体。(`http://www.python.org <http://www.python.org>`__)
Python是一种面向对象的、解释型的、灵活的语言，在计算机科学中日益流行。Python易学，语法明晰，并且能很容易的使用以C，C++或
者FORTRAN编写的模块实现扩展。

Biopython官网(`http://www.biopython.org <http://www.biopython.org>`__)为使用和研究生物信息学的开发者提供了一个在线的
资源库，包括模块、脚本以及一些基于Python的软件的网站链接。一般来讲，Biopython致力于通过创造高质量的和可重复利用的模块及
类，从而使得Python在生物信息学中的应用变得更加容易。Biopython的特点包括解析各种生物信息学格式的文件(BLAST， Clustalw， FASTA，
Genbank...)，访问在线的服务器(NCBI，Expasy...)，常见和不那么常见程序的接口(Clustalw， DSSP，MSMS...)，标准的序列类，各
种收集的模块，KD树数据结构等等，还有一些文档。

基本来说，我们喜欢使用Python来编程，并且希望通过创建高质量、可复用的模块和脚本来使得Python在生物信息学中的应用变得容易。

1.2  在Biopython包中我能发现什么？
---------------------------------------------

主要的Biopython发行版本有很多种功能，包括：

-  将生物信息学文件解析为Python可用的数据结构，包含以下支持的格式：

   -  Blast输出结果 – standalone和在线Blast
   -  Clustalw
   -  FASTA
   -  GenBank
   -  PubMed和Medline
   -  ExPASy文件, 如Enzyme和Prosite
   -  SCOP, 包括‘dom’和‘lin’文件
   -  UniGene
   -  SwissProt

-  被支持格式的文件可以通过记录来重复或者通过字典界面来索引。
-  处理常见的生物信息学在线数据库的代码：

   -  NCBI – Blast, Entrez和PubMed服务
   -  ExPASy – Swiss-Prot和Prosite条目, 包括Prosite搜索

-  常见生物信息学程序的接口，例如：

   -  NCBI的Standalone Blast
   -  Clustalw比对程序
   -  EMBOSS命令行工具

-  一个能处理序列、ID和序列特征的标准序列类。
-  对序列实现常规操作的工具，如翻译，转录和权重计算。
-  利用k最近邻接、Bayes或SVM对数据进行分类的代码。
-  处理比对的代码，包括创建和处理替换矩阵的标准方法。
-  分发并行任务到不同进程的代码。
-  实现序列的基本操作，翻译以及BLAST等功能的GUI程序。
-  使用这些模块的详细文档和帮助，包括此文件，在线的wiki文档，网站和邮件列表。
-  整合BioSQL，一个也被BioPerl和BioJava支持的数据库架构。

我们希望这些能给你足够的理由去下载并开始使用Biopython！

1.3  安装Biopython
-------------------------

Biopython的所有安装信息在此文档中分开，以便于更容易保持更新。

简短的版本去我们的下载页面(`http://biopython.org/wiki/Download <http://biopython.org/wiki/Download>`__),
下载并安装所列举的dependencies，然后下载并安装Biopython。Biopython能在多种平台上运行（Windows，Mac，各种版本的Linux和Unix）。
对于Windows我们提供预编译的一键式安装程序，而对Unix和其他操作系统，你必须按照附带的README文件从源开始安装。这通常
是很简单的，只要标准命令：

.. code:: verbatim

    python setup.py build
    python setup.py test
    sudo python setup.py install

（事实上你可以跳过build和test，直接install。但最好是确保所有的东西看起来都没问题。）

我们的安装说明的详细版本包括了Python的安装，Biopython dependencies的安装以及Biopython本身的安装。可从PDF
(`http://biopython.org/DIST/docs/install/Installation.pdf <http://biopython.org/DIST/docs/install/Installation.pdf>`__)
和HTML格式获得。
(`http://biopython.org/DIST/docs/install/Installation.html <http://biopython.org/DIST/docs/install/Installation.html>`__)。

1.4  常见问答（FAQ）
-------------------------------------

#. *在科学出版中我怎样引用Biopython？*
   请引用我们的应用笔记 [`1 <#cock2009>`__, Cock *et al.* ,  2009] 作为主要的Biopython参考文献。另外，如果可以，请
   引用以下任意出版物，特别是作为Biopython特定模块的参考文献的话。
   （更多信息可在我们网站上获得）：

   -  对于官方项目声明: [`13 <#chapman2000>`__,
      Chapman and Chang, 2000];
   -  对于 ``Bio.PDB``: [`18 <#hamelryck2003a>`__, Hamelryck and
      Manderick, 2003];
   -  对于 ``Bio.Cluster``: [`14 <#dehoon2004>`__, De Hoon *et al.*,
      2004];
   -  对于 ``Bio.Graphics.GenomeDiagram``: [`2 <#pritchard2006>`__,
      Pritchard *et al.*, 2006];
   -  对于 ``Bio.Phylo`` 和 ``Bio.Phylo.PAML``: [`9 <#talevich2012>`__,
      Talevich *et al.*, 2012];
   -  对于在Biopython，BioPerl，BioRuby，BioJava和EMBOSS支持的FASTQ格
      式文件：[`7 <#cock2010>`__, Cock *et al.*, 2010].

#. *我该怎样以大写字母写“Biopython”？写成“BioPython”可以吗？*
    正确的大写是“Biopython”而不是“BioPython”（虽然对于BioPerl，BioRuby
    和BioJava是这样）。

#. *我怎样查看自己安装的Biopython的版本？*
    使用以下代码：

   .. code:: verbatim

         >>> import Bio
         >>> print Bio.__version__
         ...
         

    如果 “\ ``import Bio``\ ” 这行报错，说明Biopython未被安装。如果第二行报错，
    你的版本已经很过时了。如果版本号以“+”号结束，说明你用的并不是官方版本，而
    是开发代码的快照。

#. *此文档的最新版本在哪里？*
    如果你下载的是一个Biopython源代码包，那么它将包含此文档HTML和PDF两种格式
    的相应版本。此文档最新出版的版本可通过在线获得（每个版本的更新）：

   -  `http://biopython.org/DIST/docs/tutorial/Tutorial.html <http://biopython.org/DIST/docs/tutorial/Tutorial.html>`__
   -  `http://biopython.org/DIST/docs/tutorial/Tutorial.pdf <http://biopython.org/DIST/docs/tutorial/Tutorial.pdf>`__

    如果你使用的是从我们库中获得的尚未发布的最新代码，你可以在这里找到还在开发中
    的教程的拷贝：

   -  `http://biopython.org/DIST/docs/tutorial/Tutorial-dev.html <http://biopython.org/DIST/docs/tutorial/Tutorial-dev.html>`__
   -  `http://biopython.org/DIST/docs/tutorial/Tutorial-dev.pdf <http://biopython.org/DIST/docs/tutorial/Tutorial-dev.pdf>`__

#. *我需要哪一个“Numerical Python”？*
    对于Biopython 1.48或更早的版本，你需要老的Numeric模块。对于Biopython 1.49
    及更高的版本，你需要更新的NumPy来代替。Numeric和NumPy都可以在同一台机器上安
    装。也可以访问： `http://numpy.scipy.org/ <http://numpy.scipy.org/>`__
#. *为什么* ``Seq`` *对象缺少了这篇教程里的（反向）transcription和translation方法？*
    你需要Biopython 1.49或更新的版本。或者，使用以下 \ `3.14 <#sec:seq-module-functions>`__ 部分中的 ``Bio.Seq`` 模块
    功能。
#. *为什么* ``Seq`` *对象缺少了这篇教程中的upper和lower方法？*
    你需要Biopython 1.53或更新版本。或者，使用 ``str(my_seq).upper()`` 来获得
    大写字符串。如果你需要一个Seq对象，试试 ``Seq(str(my_seq).upper())`` ，但是
    要小心重用相同的字母。
#. *为什么* ``Seq`` *对象的translation方法不支持本教程中描述的* ``cds`` *选项？*
    你需要Biopython 1.51或更新版本。
#. *为什么* ``Bio.SeqIO`` *不能正常工作？它导入正常但是没有解析函数等。*
    你需要Biopython 1.43或更新版本。较老的版本确实包含了一些相关的代码在 ``Bio.SeqIO`` 下面但是后来就被移除了——这就是为什么import是正常的。
#. *为什么* ``Bio.SeqIO.read()`` *不能正常工作？该模块导入正常但是并没有read函数！*
    你需要Biopython 1.45或更新的版本。或者，使用 ``Bio.SeqIO.parse(...).next()`` 来代替。
#. *为什么没有* ``Bio.AlignIO`` *？模块导入失败！*
    你需要Biopython 1.46或更新的版本。 
#. ``Bio.SeqIO`` *和* ``Bio.AlignIO`` *读写什么样的文件格式？*
    请检查内建文档（``from Bio import SeqIO``，然后 ``help(SeqIO)`` ），或见wiki上的最
    新条目：
    `http://biopython.org/wiki/SeqIO <http://biopython.org/wiki/SeqIO>`__
    以及
    `http://biopython.org/wiki/AlignIO <http://biopython.org/wiki/AlignIO>`__
#. *为什么* ``Bio.SeqIO`` *和* ``Bio.AlignIO`` *的input函数不让我提供一个序列字母？*
    你需要Biopython 1.49或更新版本。
#. *为什么* ``Bio.SeqIO`` *和* ``Bio.AlignIO`` *函数* ``parse`` *，* ``read`` *和* ``write`` *不能使用文件名？它们坚持句柄！*
    你需要Biopython 1.54或更新的版本。或者明确使用句柄。
    (见 Section \ `22.1 <#sec:appendix-handles>`__). 一定要记得当你写完数据后关闭输
    出句柄。
#. *为什么* ``Bio.SeqIO.write()`` *和* ``Bio.AlignIO.write()`` *函数不接受单个记录
   或比对？它们坚持需要一个列表或迭代器！*
    你需要Biopython 1.54或更新版本，或将该条目以 ``[...]`` 包起来形成一个单元素的列表。
#. *为什么* ``str(...)`` *不给我一个* ``Seq`` *对象的全序列？*
    你需要Biopython 1.45或更新的版本。或者，与其使用 ``str(my_seq)``，不如试试 ``my_seq.tostring()`` 这也能在最近的Biopython版本上工作）。
#. *为什么* ``Bio.Blast`` *不能处理最新的NCBI blast输出文本文件结果？*
    NCBI在不断的调整BLAST工具的纯文本输出，导致我们的解析器需要不断更新。
    如果你没使用最新版本的Biopython，你可以试试升级。但是，我们（还有NCBI）推荐你使用
    HTML格式输出来代替，因为HTML是设计给电脑程序读取的。
#. *为什么* ``Bio.Entrez.read()`` *不能正常工作？模块导入正常但是没有read函数！*
    你需要Biopython 1.46或更新的版本。
#. *为什么* ``Bio.Entrez.parse()`` *不能正常工作？模块导入正常但是没有parse函数！*
    你需要Biopython 1.52或更新的版本。
#. *为什么我的脚本使用了* ``Bio.Entrez.efetch()`` *便停止工作了？*
    这可能是由于NCBI在2012年2月引进EFetch 2.0后发生了改变。首先，他们改变了默认的返回方式——
    你可能想添加 ``retmode="text"`` 到你的call。其次，他们对于怎么提供一个ID列表变得更加严格——
    Biopython 1.59及之后版本或自动将一个列表转换成逗号分隔的字符串。
#. *为什么* ``Bio.Blast.NCBIWWW.qblast()`` *没有给出与NCBI BLAST网站上相同的结果？*
    你需要指定相同的选项——NCBI经常调整网站上的默认设置，并且他们不再匹配QBLAST的默认设置了。
    请检查gap罚分和期望值阈值。
#. *为什么* ``Bio.Blast.NCBIXML.read()`` *不正常工作？模块导入了但是没有read函数！*
    你需要Biopython 1.50或更新的版本。或者，使用 ``Bio.Blast.NCBIXML.parse(...).next()`` 代替。
#. *为什么我的* ``SeqRecord`` *对象没有一个* ``letter_annotations`` *的属性？*
    Per-letter-annotation已经被加入到Biopython 1.50中。
#. *为什么我无法切片我的* ``SeqRecord`` *来获取一个子记录？*
    你需要Biopython 1.50或更新版本。
#. *为什么我无法一起添加* ``SeqRecord`` *对象？*
    你需要Biopython 1.53或更新版本。
#. *为什么* ``Bio.SeqIO.convert()`` *或* ``Bio.AlignIO.convert()`` *不能正常工作？模块导入
    正常但是没有convert函数！*
    你需要Biopython 1.52或更新版本。或者，按以下教程中描述的结合 ``parse`` 和 ``write`` 函数。
    （见 Sections \ `5.5.2 <#sec:SeqIO-conversion>`__ 和 \ `6.2.1 <#sec:converting-alignments>`__）。
#. *为什么* ``Bio.SeqIO.index()`` *不能正常工作？模块导入正常但是没有index函数！*
    你需要Biopython 1.52或更新版本。
#. *为什么* ``Bio.SeqIO.index_db()`` *不能正常工作？模块导入正常但是没有*\ * ``index_db`` *\ *函数！*
    你需要Biopython 1.57或更新版本。（有SQLite3的Python支持）
#. ``MultipleSeqAlignment`` *对象在哪里？* ``Bio.Align`` *模块导入正常但是这个类不在那里！*
    你需要Biopython 1.54或更新版本。或者，较早的 ``Bio.Align.Generic.Alignment`` 类支持它的一些功能，
    但是现在不推荐使用这个。
#. *为什么我不能直接从应用程序包装器上运行命令行工具？*
    你需要Biopython 1.55或更新版本。或者，直接使用Python的 ``subprocess`` 模块。
#. *我看到过一个代码的目录，但是我找不到那个能干嘛的代码了。它藏在哪儿了？*
    我们知道，我们的代码存放在 ``__init__.py`` 文件里。如果你此前没有在这个文件里寻找代码那么这可能会
    让人困惑。我们这样做的原因是为了让用户更容易导入。比如，不一定要像 ``from Bio.GenBank import GenBank``
    来导入一个“repetitive”，你仅需使用 ``from Bio import GenBank`` 就行。
#. *为什么CVS的代码貌似过期了？*
    2009年9月下旬，在Biopython 1.52发布之后，我们从使用CVS转变为使用git，git是一个分散式的版本控制系统。
    旧的CVS服务仍可作为静态和只读备份，但是如果你想获取最新的代码，你需要使用git。详见我们的网站获取更多
    信息：

对于更一般的问题，Python FAQ页面 `http://www.python.org/doc/faq/ <http://www.python.org/doc/faq/>`__
可能会有帮助。
