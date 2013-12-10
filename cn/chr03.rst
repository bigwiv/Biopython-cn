第3章	生物序列对象
===========================

生物学序列以绝对的优势成为生物信息学研究的重点对象。这一章我们将简要介绍
Biopython处理这些序列的机制-- ``Seq`` 对象。第 `4 <#chapter:SeqRecord>`__ 
章将要引入与此相关的 ``SeqRecord`` 对象，此对象可以将序列信息和注释结合起来，
用于第 \ `5 <#chapter:Bio.SeqIO>`__ 章序列的输入/输出。

序列实质上就是由字母构成的字符串，比如 ``AGTACACTGGT`` ，看起来很自然，因为
这就是序列在生物学文件中的常用代表格式。

 ``Seq`` 对象和标准的Python字符串有两个明显的不同。首先，它们使用不同的方法。
尽管``Seq``对象支持常规字符串的很多方法，但是它的 ``translate()`` 方法在做
生物学翻译时是不同的。相似的还有其他的生物学相关的方法，比如 ``reverse_complement()`` 。
其次， ``Seq`` 对象具有一个重要的属性-- ``alphabet`` ，这一对象用于描述由单个
字母构成的序列字符串的 “mean” （意义），以及如何解释这一字符串。例如， ``AGTACACTGGT`` 
序列是个DNA序列还是一段富含Alanines, Glycines, Cysteines and Threonines
的蛋白质序列？

3.1  序列和字母表
----------------------------

字母对象可能是使得``Seq``对象不仅仅是字符串的重要因素。目前，Biopython的字母表
定义在``Bio.Alphabet``模块。我们将会使用IUPAC字母表
（ `http://www.chem.qmw.ac.uk/iupac/ <http://www.chem.qmw.ac.uk/iupac/>`__ ）
来处理我们比较青睐的对象：DNA、RNA和蛋白质对象。

 ``Bio.Alphabet.IUPAC`` 提供了蛋白质、DNA和RNA的基本定义，并且提供了扩展
和定制基本定义的功能。例如蛋白质，有一个基本的IUPACProtein类，另外还有
一个ExtendedIUPACProtein类。这个类包含除20种常见氨基酸外的其他氨基酸元素，
比如 “U” （或 “Sec”代表硒代半胱氨酸），“O”（或“Pyl” 代表吡咯赖氨酸），
还有歧意字母“B” （或“Asx”代表天冬酰胺或者天冬氨酸），“Z”（或“Glx”代表谷
氨酰胺或者谷氨酸），“J” （或“Xle”代表亮氨酸或异亮氨酸），“X” （或“Xxx”代
表未知氨基酸）。同理，对于DNA有IUPACUnambiguousDNA、IUPACAmbiguousDNA和
ExtendedIUPACDNA类，分别提供基本字母，每种可能下的歧意字母和修饰后的碱基。
同样地，RNA可以使用IUPACAmbiguousRNA和IUPACUnambiguousRNA类。

使用字母表类有两方面的优势。首先，明确了Seq对象包含的信息的类型；其次通过
类型检查，它提供了约束信息的工具。

我们已经知道了将要处理的对象，现在让我们看看怎么使用这些类做一些有意思的事情。
你可以创建一条有通用字母组成的模糊序列，如下：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> my_seq = Seq("AGTACACTGGT")
    >>> my_seq
    Seq('AGTACACTGGT', Alphabet())
    >>> my_seq.alphabet
    Alphabet()

然而，如果可能，你要在创建序列对象的时候就尽量明确指定字母的类型，如下创建一条
明确的DNA字母表对象。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
    >>> my_seq
    Seq('AGTACACTGGT', IUPACUnambiguousDNA())
    >>> my_seq.alphabet
    IUPACUnambiguousDNA()

当然,除非这真的是一个氨基酸序列：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_prot = Seq("AGTACACTGGT", IUPAC.protein)
    >>> my_prot
    Seq('AGTACACTGGT', IUPACProtein())
    >>> my_prot.alphabet
    IUPACProtein()

3.2  序列表现的像字符串一样
-------------------------------

在许多时候，我们可以讲Seq对象处理成正常的Python字符串，比如取序列长度，迭代
元素：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = Seq("GATCG", IUPAC.unambiguous_dna)
    >>> for index, letter in enumerate(my_seq):
    ...     print index, letter
    0 G
    1 A
    2 T
    3 C
    4 G
    >>> print len(my_seq)
    5

你可以像字符串那样获取序列的元素（但是请记住，Python计数从0开始）：

.. code:: verbatim

    >>> print my_seq[0] #first letter
    G
    >>> print my_seq[2] #third letter
    T
    >>> print my_seq[-1] #last letter
    G

``Seq`` 对象有一个 ``.count()`` 方法，类似于字符串。记住这意味就像Python的
字符串一样进行着非重叠的计数。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> "AAAA".count("AA")
    2
    >>> Seq("AAAA").count("AA")
    2

但是在某些生物学上，你可能需要使用重叠计数（就像上面的例子中如果重复计
数结果将为3）。当计算耽搁字母出现的次数时，重叠计数和非重叠计数没有差别。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPAC.unambiguous_dna)
    >>> len(my_seq)
    32
    >>> my_seq.count("G")
    9
    >>> 100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq)
    46.875

你当然可以使用上面的代码段计算GC含量，但是记住 ``Bio.SeqUtils`` 模块已经
建立了好几个GC函数，类如：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> from Bio.SeqUtils import GC
    >>> my_seq = Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPAC.unambiguous_dna)
    >>> GC(my_seq)
    46.875

注意在使用 ``Bio.SeqUtils.GC()`` 函数时会自动处理序列和可代表G或者C的歧意核苷酸
字母S混合的情况。

然后还要注意，就像正常的Python字符串， ``Seq`` 对象在某些方式下是只读的。如果需要
编辑序列，比如模拟点突变，请看后续的 \ `3.12 <#sec:mutable-seq>`__章节中讲述的
 ``MutableSeq`` 对象。

3.3  切取序列
-----------------------

一个较为复杂的例子，让我们切取序列。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
    >>> my_seq[4:12]
    Seq('GATGGGCC', IUPACUnambiguousDNA())

要注意两个有意思的地方。首先，序列第一个元素从0开始，这是符合Python字符串的规则的。
这在计算机科学上是普遍现象，但在生物学上不是这样。当你做切片的时候，第一项包含了
（比如例子中的4），而最后一项去除了（例子中的12）。这是Python的规则，但当然这不是
世界上所有人都希望的。主要是为了和Python保持一致。

第二个需要注意的地方是，切片是在序列数据字符串上执行的，但是产生的新对象保留了原始
 ``Seq`` 对象的字母表信息。

同样和Python字符串一样，你可以通过设置起始位置、终止位置和 *步幅* （间隔数，默认为1）
进行切片。例如，我们可以分别获取下面DNA序列密码子第一、第二、第三位的碱基。

.. code:: verbatim

    >>> my_seq[0::3]
    Seq('GCTGTAGTAAG', IUPACUnambiguousDNA())
    >>> my_seq[1::3]
    Seq('AGGCATGCATC', IUPACUnambiguousDNA())
    >>> my_seq[2::3]
    Seq('TAGCTAAGAC', IUPACUnambiguousDNA())

你可能已经注意到Python字符串中的另一个奇特步幅设定：使用-1返回倒序字符串切片。
当然以也可以使用 ``Seq`` 对象来完成。

.. code:: verbatim

    >>> my_seq[::-1]
    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG', IUPACUnambiguousDNA())

3.4  将序列对象转换成字符串
-------------------------------------

如果你仅仅需要一个单纯的字符串，就像写入文件或者插入数据库，这事很容易就
可以实现的：

.. code:: verbatim

    >>> str(my_seq)
    'GATCGATGGGCCTATATAGGATCGAAAATCGC'

尽管对 ``Seq`` 对象调用 ``str()`` 方法将以字符串的形式返回全长序列，但是你经常不需要
特地做这个转换。当使用print打印声明是，Python会自动转换。

.. code:: verbatim

    >>> print my_seq
    GATCGATGGGCCTATATAGGATCGAAAATCGC

当你进行Python字符串格式化或者插入操作符（ ``%`` ）时，
可以直接把 ``Seq`` 对象和 ``%s`` 占位符一起使用：

.. code:: verbatim

    >>> fasta_format_string = ">Name\n%s\n" % my_seq
    >>> print fasta_format_string
    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    <BLANKLINE>

这一行代码展示的是一个简单的FASTA格式的记录（不用关心自动换行）。
\ `4.5 <#sec:SeqRecord-format>`__ 部分将介绍一个简洁的方式从 ``SeqRecord`` 
对象中获取FASTA格式的字符串，更详细的读写FASTA格式的序列文件将在第
\ `5 <#chapter:Bio.SeqIO>`__ 章介绍。

 *注意：*  如果你使用Biopython 1.44或者更旧的版本，使用 ``str(my_seq)`` 只会
返回一个截短了的序列。这时候可以使用 ``my_seq.tostring()`` ，为了保持向后兼
容性，这一方法在当前的Biopython版本中还有保留。

.. code:: verbatim

    >>> my_seq.tostring()
    'GATCGATGGGCCTATATAGGATCGAAAATCGC'

3.5  连接或添加序列
--------------------------------------

当然，原则上你可以将任何两个Seq对象加在一起，就像Python字符串一样去连接它们。
但是你不能将两个不相容的字母表加在一起，比如蛋白质序列和核苷酸序列就不能简单
叠加。

.. code:: verbatim

    >>> from Bio.Alphabet import IUPAC
    >>> from Bio.Seq import Seq
    >>> protein_seq = Seq("EVRNAK", IUPAC.protein)
    >>> dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
    >>> protein_seq + dna_seq
    Traceback (most recent call last):
    ...
    TypeError: Incompatible alphabets IUPACProtein() and IUPACUnambiguousDNA()

如果你 *真的* 想这么做，你必须首先将两个序列转换成通用字母表。

.. code:: verbatim

    >>> from Bio.Alphabet import generic_alphabet
    >>> protein_seq.alphabet = generic_alphabet
    >>> dna_seq.alphabet = generic_alphabet
    >>> protein_seq + dna_seq
    Seq('EVRNAKACGT', Alphabet())

这里有个例子是将通用核苷酸序列加到明确的IUPAC DNA序列上，最后生成一段
模糊的核苷酸序列。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_nucleotide
    >>> from Bio.Alphabet import IUPAC
    >>> nuc_seq = Seq("GATCGATGC", generic_nucleotide)
    >>> dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
    >>> nuc_seq
    Seq('GATCGATGC', NucleotideAlphabet())
    >>> dna_seq
    Seq('ACGT', IUPACUnambiguousDNA())
    >>> nuc_seq + dna_seq
    Seq('GATCGATGCACGT', NucleotideAlphabet())

3.6  改变大小写
------------------

Python字符串具有很有用的转换大小写的 ``upper`` 和 ``lower`` 方法。从
Biopython 1.53起， ``Seq`` 对象也获取了类似的方法应用于字母表。例如:

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna
    >>> dna_seq = Seq("acgtACGT", generic_dna)
    >>> dna_seq
    Seq('acgtACGT', DNAAlphabet())
    >>> dna_seq.upper()
    Seq('ACGTACGT', DNAAlphabet())
    >>> dna_seq.lower()
    Seq('acgtacgt', DNAAlphabet())

这在不区分大小写进行匹配的时候很有用。

.. code:: verbatim

    >>> "GTAC" in dna_seq
    False
    >>> "GTAC" in dna_seq.upper()
    True

注意，严格地说IUPAC字母表仅仅是对于大写字母构成的序列的，因此：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
    >>> dna_seq
    Seq('ACGT', IUPACUnambiguousDNA())
    >>> dna_seq.lower()
    Seq('acgt', DNAAlphabet())

3.7  核苷酸序列和（反向）互补序列
---------------------------------------------------

对于核苷酸序列，你可以使用 ``Seq`` 对象内置的方法很容易地获得 ``Seq`` 
的互补或反向互补序列。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
    >>> my_seq
    Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPACUnambiguousDNA())
    >>> my_seq.complement()
    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG', IUPACUnambiguousDNA())
    >>> my_seq.reverse_complement()
    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC', IUPACUnambiguousDNA())

在前面的方法中，使用切片的-1的步长可以很容易的获取一个 ``Seq`` 对象的反向序列。

.. code:: verbatim

    >>> my_seq[::-1]
    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG', IUPACUnambiguousDNA())

在所有这些操作中，字母的属性一直保留着。这是非常有用的，以防你不小心做一些
奇怪的事情，比如获取蛋白质序列的（反向）互补序列。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> protein_seq = Seq("EVRNAK", IUPAC.protein)
    >>> protein_seq.complement()
    Traceback (most recent call last):
    ...
    ValueError: Proteins do not have complements!

  \ `5.5.3 <#sec:SeqIO-reverse-complement>`__ 部分的例子将 ``Seq`` 对象的反向互补
 方法和 ``Bio.SeqIO`` 对于序列的输入/输出方法结合起来。

3.8  转录
------------------

在谈论转录之前，我想先说明一下链的问题。考虑以下（编造的）编码短肽的双链DNA的延伸：

+------+------------------------------------------------------+------+
|      |
+------+------------------------------------------------------+------+
|      | DNA coding strand (aka Crick strand, strand +1)      |      |
+------+------------------------------------------------------+------+
| 5’   | ``ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG``          | 3’   |
+------+------------------------------------------------------+------+
|      | ``|||||||||||||||||||||||||||||||||||||||``          |      |
+------+------------------------------------------------------+------+
| 3’   | ``TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC``          | 5’   |
+------+------------------------------------------------------+------+
|      | DNA template strand (aka Watson strand, strand −1)   |      |
+------+------------------------------------------------------+------+
|      |
+------+------------------------------------------------------+------+
|      | \|                                                   |      |
+------+------------------------------------------------------+------+
|      | Transcription                                        |      |
+------+------------------------------------------------------+------+
|      | ↓                                                    |      |
+------+------------------------------------------------------+------+
|      |
+------+------------------------------------------------------+------+
| 5’   | ``AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG``          | 3’   |
+------+------------------------------------------------------+------+
|      | Single stranded messenger RNA                        |      |
+------+------------------------------------------------------+------+
|      |
+------+------------------------------------------------------+------+

实际的生物学上的转录过程是将模板链反向互补（TCAG → CUGA）生成mRNA。但是，
在Biopython和生物信息学领域，我们通常会直接利用编码链，因为我们可以通过
T → U的转换获得mRNA。

现在让我们着手真实地使用Biopython做一个转录。首先，让我们分别创建DNA序列的
编码链和模板链的 ``Seq`` 对象：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
    >>> coding_dna
    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())
    >>> template_dna = coding_dna.reverse_complement()
    >>> template_dna
    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT', IUPACUnambiguousDNA())

这是和上面的图表相一致的，记住按照惯例核苷酸序列通常是从5’到3’端方向的，
而图中所示的模板链是反向的。

现在让我们使用 ``Seq`` 对象内置的 ``transcribe`` 方法将编码链转录成对应的mRNA：

.. code:: verbatim

    >>> coding_dna
    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())
    >>> messenger_rna = coding_dna.transcribe()
    >>> messenger_rna
    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

就如你看到的，这里做的全部工作是将T → U转换，并调整字母表。

如果你确实想从模板链去做一个真正的生物学上的转录，需要两步：

.. code:: verbatim

    >>> template_dna.reverse_complement().transcribe()
    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

``Seq`` 对象还包含了从mRNA逆向转录为DNA编码链的方法。同样，这仅仅是从U
→ T的替代并伴随着字母表的变化：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", IUPAC.unambiguous_rna)
    >>> messenger_rna
    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())
    >>> messenger_rna.back_transcribe()
    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())

*注意：* ``Seq`` 对象的 ``transcribe`` 和 ``back_transcribe`` 方法直到
Biopython 1.49版本才出现，在较早的版本中你需要使用 ``Bio.Seq`` 模块的函
数替代，详见 \ `3.14 <#sec:seq-module-functions>`__ 部分。


3.9  翻译
----------------

继续使用在转录那个小节中的例子，现在让我们将这个mRNA翻译成相对应的
蛋白质序列，利用的是 ``Seq`` 对象众多生物学方法中的一个：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", IUPAC.unambiguous_rna)
    >>> messenger_rna
    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())
    >>> messenger_rna.translate()
    Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

你也可以直接从编码DNA链进行翻译：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
    >>> coding_dna
    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())
    >>> coding_dna.translate()
    Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

你应该注意到在上面的蛋白质序列中，除了末尾的终止符外，在序列中间还有一个终止符。
其实这是一个精心选择的例子，因为由它我们可以引申讲一下可选参数，包括不同的翻译
表（遗传密码）。

Biopython上可用的翻译表是基于 `NCBI <http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`__ 
（参考这个教程的下一个部分）。默认情况下，翻译使用的是 *标准* 遗传密码（NCBI上table id 1)。
假设我们需要翻译一个线粒体序列，我们就需要告诉翻译函数使用相关的遗传密码：

.. code:: verbatim

    >>> coding_dna.translate(table="Vertebrate Mitochondrial")
    Seq('MAIVMGRWKGAR*', HasStopCodon(IUPACProtein(), '*'))

你也可以利用NCBI上表格的标号来指定所使用的遗传密码，这样更简洁一些，
在GenBank文件的特征注释中经常包含表格的标号：

.. code:: verbatim

    >>> coding_dna.translate(table=2)
    Seq('MAIVMGRWKGAR*', HasStopCodon(IUPACProtein(), '*'))

现在你可能想将上面的核苷酸序列仅翻译到阅读框的第一个终止密码子，然后停止
（这更符合自然现象）。

.. code:: verbatim

    >>> coding_dna.translate()
    Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))
    >>> coding_dna.translate(to_stop=True)
    Seq('MAIVMGR', IUPACProtein())
    >>> coding_dna.translate(table=2)
    Seq('MAIVMGRWKGAR*', HasStopCodon(IUPACProtein(), '*'))
    >>> coding_dna.translate(table=2, to_stop=True)
    Seq('MAIVMGRWKGAR', IUPACProtein())

注意到当你使用 ``to_stop`` 参数时，终止密码子本身是不翻译的，终止的符号也是
不显现在蛋白质序列中的。

如果你不喜欢默认的星号作为终止符号，你也可以自己指定终止符。

.. code:: verbatim

    >>> coding_dna.translate(table=2, stop_symbol="@")
    Seq('MAIVMGRWKGAR@', HasStopCodon(IUPACProtein(), '@'))

现在假设你有一条完整的编码序列CDS，这是一种核苷酸序列（例如mRNA剪切以后），
序列全长都是密码子（也就是长度是3的倍数），开始于起始密码子，终止于终止密
码子，阅读框内没有内部的终止密码子。通常情况下，给你一条完整的CDS，默认的
翻译方法即可以翻译出你想要的（有时使用 ``to_stop`` 选项）。但是，如果序列使
用的是非标准的起始密码子呢？这种情况在细菌中很常见，比如 ``E. coli`` 
K12中的基因yaaX：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna
    >>> gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" + \
    ...            "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT" + \
    ...            "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT" + \
    ...            "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT" + \
    ...            "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA",
    ...            generic_dna)
    >>> gene.translate(table="Bacterial")
    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR*',
    HasStopCodon(ExtendedIUPACProtein(), '*')
    >>> gene.translate(table="Bacterial", to_stop=True)
    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR',
    ExtendedIUPACProtein())

在细菌遗传密码中 ``GTG`` 是个有效的起始密码子。 *正常情况下* 编码缬氨酸，
如果作为起始密码子，则翻译成甲硫氨酸。当你告诉Biopython你的序列是完整CDS时，
这事将会发生。

.. code:: verbatim

    >>> gene.translate(table="Bacterial", cds=True)
    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR',
    ExtendedIUPACProtein())

除了告诉Biopython翻译时使用另一种起始密码子编码甲硫氨酸外，使用这一选项同样能
确保你的序列是个真实有效的CDS（如果不是将会抛出异常）。

第 \ `18.1.3 <#sec:SeqIO-translate>`__ 章的例子将把 ``Seq`` 对象的翻译方法和
 ``Bio.SeqIO`` 对象的对于序列的输入/输出方法结合起来。 

3.10  翻译表
------------------------

在前面的章节中我们讨论了 ``Seq`` 对象的转录方法（并且提到了 ``Bio.Seq`` 模块
中的等效函数--参见第 \ `3.14 <#sec:seq-module-functions>`__ 章节）。实质上
使用的这些密码子表对象来自与NCBI的 ```ftp://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt`` 
<ftp://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt>`__ ，还有
 `http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi <http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`__ 
以一种更易读的形式呈现。

和前面一样，让我们仅仅关注两个选择：标准的翻译表和脊椎动物线粒体DNA的翻译表。

.. code:: verbatim

    >>> from Bio.Data import CodonTable
    >>> standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    >>> mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

另一种方式，这些表也可以分别以标号1和2来标识：

.. code:: verbatim

    >>> from Bio.Data import CodonTable
    >>> standard_table = CodonTable.unambiguous_dna_by_id[1]
    >>> mito_table = CodonTable.unambiguous_dna_by_id[2]

你可以在打印后直观地比较这些实际的翻译表：

.. code:: verbatim

    >>> print standard_table
    Table 1 Standard, SGC0

      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--

和

.. code:: verbatim

    >>> print mito_table
    Table 2 Vertebrate Mitochondrial, SGC1

      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--

你会发现下面的特性很有用，比如当你查找新基因时：

.. code:: verbatim

    >>> mito_table.stop_codons
    ['TAA', 'TAG', 'AGA', 'AGG']
    >>> mito_table.start_codons
    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
    >>> mito_table.forward_table["ACG"]
    'T'

3.11  比较Seq对象
---------------------------

序列之间的比较实际上是一个比较复杂的话题，没有简单的方法来判断两个序列是等同的。
核心的问题是字母的意义是依赖于上下文的。字母 “A” 既可以是DNA、RNA也可以使蛋白质序
列的一部分。 Biopython在 ``Seq`` 对象中包含了字母表对象，以此尝试获得这些信息。所
以比较两个 ``Seq`` 对象意味着既要考虑两个序列的字符串 *又要* 考虑字母表。

举个例子，你可能会觉得 ``Seq("ACGT", IUPAC.unambiguous_dna)`` 和
 ``Seq("ACGT", IUPAC.ambiguous_dna)`` 这两个DNA ``Seq`` 对象是一样的，尽管它们确实具
有不同的字母表。根据上下文来判断是很重要的。

下面这种情况更遭：假设你认为 ``Seq("ACGT", IUPAC.unambiguous_dna)`` 和
 ``Seq("ACGT")`` （也就是默认的通用字母表）是等同的。那么依照逻辑，
  ``Seq("ACGT", IUPAC.protein)`` 和 ``Seq("ACGT")`` 也是等同的。现在从理
 论上讲，如果 *A*\ =\ *B* ， *B*\ =\ *C* ，那么通过递延性，我们会期望
 *A*\ =\ *C* 。因此遵从逻辑上的一致性我们需要将 ``Seq("ACGT", IUPAC.unambiguous_dna)`` 
和 ``Seq("ACGT", IUPAC.protein)`` 等同起来，虽然大部分人会同意这一递延，
但是这是错误的。这一递延性的问题也会影响使用 ``Seq`` 对象作为Python字典
的键值。

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> seq1 = Seq("ACGT", IUPAC.unambiguous_dna)
    >>> seq2 = Seq("ACGT", IUPAC.unambiguous_dna)

那么接下来Biopython会怎么做？等同性测试是Python对象默认要做的测试。经过
检验查看内存中的对象是不是同一个。这是一个非常严格的测试：

.. code:: verbatim

    >>> seq1 == seq2
    False
    >>> seq1 == seq1
    True

如果你真想这么做，你可以更明确地使用Python中的 ``id`` 函数，

.. code:: verbatim

    >>> id(seq1) == id(seq2)
    False
    >>> id(seq1) == id(seq1)
    True

在日常使用中，你的所有序列可能都是同一个字母表，或者至少都是同一类型的序列
（都是DNA、RNA或者都是蛋白质）。你可能想要的只是以字符串的形式比较这些序列，
那么直接这么做：

.. code:: verbatim

    >>> str(seq1) == str(seq2)
    True
    >>> str(seq1) == str(seq1)
    True

作为一个扩展，你可以建立一个Python字典，以 ``Seq`` 对象作为键值。一般情况下，
将序列作为字符串赋予键值更有用。详见 \ `3.4 <#sec:seq-to-string>`__ 部分。

3.12  MutableSeq对象
------------------------

就像正常的Python字符串， ``Seq`` 对象是 “只读的” ，在Python术语上就是不可变的。
除了想要 ``Seq`` 对象表现得向一个字符串之外，这是一个很有用的默认，因为在生
物学应用上你往往需要确保你没有改动你的序列数据：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)

当你尝试编辑序列是你看看会发生什么：

.. code:: verbatim

    >>> my_seq[5] = "G"
    Traceback (most recent call last):
    ...
    TypeError: 'Seq' object does not support item assignment

但是你可以使用 ``MutableSeq`` 对象将它转换成可变的序列，然后做任何你想要做的。

.. code:: verbatim

    >>> mutable_seq = my_seq.tomutable()
    >>> mutable_seq
    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA', IUPACUnambiguousDNA())

或者你可以直接从字符串建立一个 ``MutableSeq`` 对象：

.. code:: verbatim

    >>> from Bio.Seq import MutableSeq
    >>> from Bio.Alphabet import IUPAC
    >>> mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)

这两种方式都可以将序列对象转换成可变的：

.. code:: verbatim

    >>> mutable_seq
    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA', IUPACUnambiguousDNA())
    >>> mutable_seq[5] = "C"
    >>> mutable_seq
    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA', IUPACUnambiguousDNA())
    >>> mutable_seq.remove("T")
    >>> mutable_seq
    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA', IUPACUnambiguousDNA())
    >>> mutable_seq.reverse()
    >>> mutable_seq
    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG', IUPACUnambiguousDNA())

请注意与 ``Seq`` 对象不同的是， ``MutableSeq`` 对象的各种方法都是实时呈现的，比如
 ``reverse_complement()`` 和 ``reverse()`` 方法！

Python中可变对象和不可变对象的一个重要的技术差别就是 ``MutableSeq`` 对象不可以作为
字典的键值 ，但是Python字符串或者 ``Seq`` 对象就可以。

一旦你的 ``MutableSeq`` 对象编辑完成，很容易将它变回到只读的 ``Seq`` 对象，你只需：

.. code:: verbatim

    >>> new_seq = mutable_seq.toseq()
    >>> new_seq
    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG', IUPACUnambiguousDNA())

就像你从 ``Seq`` 对象中获取字符串一样，你也可以从 ``MutableSeq`` 获得（参见
 `3.4 <#sec:seq-to-string>`__ 章节）。

3.13  UnknownSeq对象
------------------------

 ``UnknownSeq`` 对象是基本的 ``Seq`` 对象中的一个子类，其目的是一个已知长度的
序列，但序列并不是由实际的字母组成的。在这种情况下，你当然可以将其作为一个
正常的 ``Seq`` 对象，但是存储由一百万个 “N” 字母组成的字符串会浪费相当大量的内
存，这时你可以只存储一个 “N” 和序列所需的长度（整数）。

.. code:: verbatim

    >>> from Bio.Seq import UnknownSeq
    >>> unk = UnknownSeq(20)
    >>> unk
    UnknownSeq(20, alphabet = Alphabet(), character = '?')
    >>> print unk
    ????????????????????
    >>> len(unk)
    20

当然你也可以指定一个字母，而不仅仅是 “?” 。一般核苷酸序列默认为 “N” ，蛋白质
序列默认为 “X” 。

.. code:: verbatim

    >>> from Bio.Seq import UnknownSeq
    >>> from Bio.Alphabet import IUPAC
    >>> unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna)
    >>> unk_dna
    UnknownSeq(20, alphabet = IUPACAmbiguousDNA(), character = 'N')
    >>> print unk_dna
    NNNNNNNNNNNNNNNNNNNN

你可以使用所有常规的 ``Seq`` 对象，记住这些可以节省内存的 ``UnknownSeq`` 对象，
如你所希望的那样在恰当的地方使用。

.. code:: verbatim

    >>> unk_dna
    UnknownSeq(20, alphabet = IUPACAmbiguousDNA(), character = 'N')
    >>> unk_dna.complement()
    UnknownSeq(20, alphabet = IUPACAmbiguousDNA(), character = 'N')
    >>> unk_dna.reverse_complement()
    UnknownSeq(20, alphabet = IUPACAmbiguousDNA(), character = 'N')
    >>> unk_dna.transcribe()
    UnknownSeq(20, alphabet = IUPACAmbiguousRNA(), character = 'N')
    >>> unk_protein = unk_dna.translate()
    >>> unk_protein
    UnknownSeq(6, alphabet = ProteinAlphabet(), character = 'X')
    >>> print unk_protein
    XXXXXX
    >>> len(unk_protein)
    6

你也许能够在自己的代码中找到 ``UnknownSeq`` 对象的应用，但你更可能首先在由
 ``Bio.SeqIO`` 创建的 ``SeqRecord`` 对象中遇到 ``UnknownSeq`` 对象（参见第
 \ `5 <#chapter:Bio.SeqIO>`__ 章）。一些序列格式的文件不总是由实际的序列组成，
像GenBank和EMBL文件就可能包含各种特征的列表，而序列部分仅展示contig信息。
又或者在测序工作中的QUAL文件仅包含质量分数，而 *从未* 包含序列，取而代之的
和QUAL文件同时生成的FASTA格式文件 *确实* 是由序列构成。

3.14  直接使用字符串
-----------------------------------

在这一章的结尾，对于那些 *真的* 不想使用序列对象的人（或者那些更喜欢面向
对象的函数式编程风格的人）， ``Bio.Seq`` 的模块级别的函数可以接受普通的
Python字符串，比如 ``Seq`` 对象（包括 ``UnknownSeq`` 对象）或者 ``MutableSeq`` 对象：

.. code:: verbatim

    >>> from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
    >>> my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
    >>> reverse_complement(my_string)
    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'
    >>> transcribe(my_string)
    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'
    >>> back_transcribe(my_string)
    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'
    >>> translate(my_string)
    'AVMGRWKGGRAAG*'

尽管如此，我们鼓励你使用默认的 ``Seq`` 对象。
