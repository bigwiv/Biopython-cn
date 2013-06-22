第6章 多序列比对
==============================================

多序列比对（Multiple Sequence Alignment, MSA）是指对多个序列进行对位排列。 这通常需要保证序列间的等同位点处在同一列上，并通过引进小横线（-）以保证最终的序列具有相同的长度。这种序列排列可以视作是由字符组成的矩阵。在Biopython中，多序列排列中每一个序列是以 ``SeqRecord`` 对象来表示的。

这里我们介绍一种新的对象 -- ``MultipleSeqAlignment`` 来表示这样一类数据，我们还将介绍 ``Bio.AlignIO`` 模块来读写不同格式的多序列比对数据（ ``Bio.AlignIO`` 在设计上与之前介绍的 ``Bio.SeqIO`` 模块是类似的）。Biopython中， ``Bio.SeqIO`` 和 ``Bio.AlignIO`` 都能读写各种格式的多序列排列数据。在实际处理中，使用哪一个模块取决于用户需要对数据进行何种操作。

本章的第一部分是关于各种常用多序列排列程序（ClustalW和MUSCLE）的Biopython命令行封装。

6.1 读取多序列排列数据
-------------------------------------------

在Biopython中，有两种方法读取多序列排列数据， ``Bio.AlignIO.read()`` 和 ``Bio.AlignIO.parse()`` 。这两种方法跟 ``Bio.SeqIO`` 处理一个和多个数据的设计方式是一样的。 ``Bio.AlignIO.read()`` 只能读取一个多序列排列而 ``Bio.AlignIO.parse()`` 可以依次读取多个序列排列数据。 

使用 ``Bio.AlignIO.parse()`` 将会返回一个 ``MultipleSeqAlignment`` 的 *迭代器（iterator）* 。迭代器往往在循环中使用。在实际数据分析过程中会时常处理包含有多个序列排列的文件。例如PHYLIP中的 ``seqboot`` ，EMBOSS工具箱中的 ``water`` 和 ``needle``, 以及Bill Pearson的FASTA程序。

然而在大多数情况下，你所遇到的文件仅仅包括一个多序列排列。这时，你应该使用 ``Bio.AlignIO.read()`` ，这将返回一个 ``MultipleSeqAlignment`` 对象。

这两个函数都接受两个必须参数。

#. 第一个参数为包含有多序列排列数据的 *句柄（handle）* 。在实际操作中，这往往是一个具有可读权限的句柄对象（详细信息请见 `22.1 <#sec:appendix-handles>`__ ）或者一个储存数据的文件名。

#. 第二个参数为文件格式（小写）。与 ``Bio.SeqIO`` 模块一样，Biopython不会对将读取的文件格式进行猜测。所有 ``Bio.AlignIO`` 模块支持的多序列排列数据格式可以在 ```http://biopython.org/wiki/AlignIO`` <http://biopython.org/wiki/AlignIO>`__ 中找到。

``Bio.AlignIO`` 模块还接受一个可选参数 ``seq_count`` 。这一参数将在 `6.1.3 <#sec:AlignIO-count-argument>`__ 中具体讨论。它可以处理不确定的多序列排列格式，或者包含有多个序列的排列。

另一个可选参数 ``alphabet`` 允许用户指定序列排列文件的字符（alphabet），它可以用来说明序列排列的类型（DNA，RNA或蛋白质）。因为大多数序列排列格式并不区别序列的类型，因此指定这一参数可能会对后期的分析产生帮助。 ``Bio.AlignIO`` 默认将使用一般字符（generic alphabet），这将不区分各种序列排列类型。

6.1.1 单一的序列排列
~~~~~~~~~~~~~~~~~~~~~~~~

例如，请见以下PFAM（或者Stockholm）格式的蛋白序列排列文件。

.. code:: verbatim

    # STOCKHOLM 1.0
    #=GS COATB_BPIKE/30-81  AC P03620.1
    #=GS COATB_BPIKE/30-81  DR PDB; 1ifl ; 1-52;
    #=GS Q9T0Q8_BPIKE/1-52  AC Q9T0Q8.1
    #=GS COATB_BPI22/32-83  AC P15416.1
    #=GS COATB_BPM13/24-72  AC P69541.1
    #=GS COATB_BPM13/24-72  DR PDB; 2cpb ; 1-49;
    #=GS COATB_BPM13/24-72  DR PDB; 2cps ; 1-49;
    #=GS COATB_BPZJ2/1-49   AC P03618.1
    #=GS Q9T0Q9_BPFD/1-49   AC Q9T0Q9.1
    #=GS Q9T0Q9_BPFD/1-49   DR PDB; 1nh4 A; 1-49;
    #=GS COATB_BPIF1/22-73  AC P03619.2
    #=GS COATB_BPIF1/22-73  DR PDB; 1ifk ; 1-50;
    COATB_BPIKE/30-81             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
    #=GR COATB_BPIKE/30-81  SS    -HHHHHHHHHHHHHH--HHHHHHHH--HHHHHHHHHHHHHHHHHHHHH----
    Q9T0Q8_BPIKE/1-52             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
    COATB_BPI22/32-83             DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
    COATB_BPM13/24-72             AEGDDP...AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
    #=GR COATB_BPM13/24-72  SS    ---S-T...CHCHHHHCCCCTCCCTTCHHHHHHHHHHHHHHHHHHHHCTT--
    COATB_BPZJ2/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
    Q9T0Q9_BPFD/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
    #=GR Q9T0Q9_BPFD/1-49   SS    ------...-HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH--
    COATB_BPIF1/22-73             FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA
    #=GR COATB_BPIF1/22-73  SS    XX-HHHH--HHHHHH--HHHHHHH--HHHHHHHHHHHHHHHHHHHHHHH---
    #=GC SS_cons                  XHHHHHHHHHHHHHHHCHHHHHHHHCHHHHHHHHHHHHHHHHHHHHHHHC--
    #=GC seq_cons                 AEssss...AptAhDSLpspAT-hIu.sWshVsslVsAsluIKLFKKFsSKA
    //

这是PFAM数据库中Phage\_Coat\_Gp8的种子排列（PF05371）。该排列下载于一个已经过期的PFAM数据库版本（ ```http://pfam.sanger.ac.uk/`` <http://pfam.sanger.ac.uk/>`__ ），但这并不影响我们的例子。假设你已经将以上内容下载到一个名为''PF05371\_seed.sth''的文件中，并在Python的当前工作目录下。

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")

这段代码将在屏幕上打印出该序列排列的概要信息：

.. code:: verbatim

    >>> print alignment
    SingleLetterAlphabet() alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73

你会注意到，以上输出截短了中间一部分序列的内容。你也可以很容易地通过控制多序列排列中每一个序列（为 ``SeqRecord`` 对象）来输出你所喜欢的格式。例如：

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    >>> print "Alignment length %i" % alignment.get_alignment_length()
    Alignment length 52
    >>> for record in alignment:
    ...     print "%s - %s" % (record.seq, record.id)
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73

你也可以使用上面alignment对象的 ``format`` 方法来以指定的格式显示它。具体信息可以参见 `6.2.2 <#sec:alignment-format-method>`__ 。

你是否已经注意到以上原始数据文件中包含有引用蛋白数据库（PDB）以及相关二级结构的信息？你可以尝试一下代码：

.. code:: verbatim

    >>> for record in alignment:
    ...     if record.dbxrefs:
    ...         print record.id, record.dbxrefs
    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']

如果你希望显示所有的序列注释信息，请使用以下例子：

.. code:: verbatim

    >>> for record in alignment:
    ...     print record

Sanger网站
```http://pfam.sanger.ac.uk/family?acc=PF05371`` <http://pfam.sanger.ac.uk/family?acc=PF05371>`__
可以让你下载各种不同的序列排列的格式。以下例子为FASTA格式：

.. code:: verbatim

    >COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
    >Q9T0Q8_BPIKE/1-52
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
    >COATB_BPI22/32-83
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
    >COATB_BPM13/24-72
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
    >COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
    >Q9T0Q9_BPFD/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
    >COATB_BPIF1/22-73
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA

注意Sanger网站有一个选项可以将序列排列中的间隔（gap）用小圆点或者是小横线表示。在以上例子中，序列间隔由小横线表示。假设你已经下载该文件，并保存为 “PF05371\_seed.faa”。你可以使用以下代码来读入该序列排列。

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.faa", "fasta")
    print alignment

你可能已经发现，以上代码中唯一的变化只是指定格式的参数。所返回的alignment对象将会包含同样的序列和序列名字。但是仔细的读者会发现，每一个alignment的SeqRecord中并不包含数据的引用注释。这是因为FASTA格式本身并没有包含这一类信息。

此外，除了使用Sanger网站，你也可以利用 ``Bio.AlignIO`` 来将原始的Stockholm格式转化成FASTA文件格式（见以下代码）。

对于任何一种Biopython支持的格式，你都可以用一样的方式读取它（通过指定文件的格式）。例如，你可以使用“phylip”来表示PHYLIP格式文件，用"nexus"来指定NEXUS格式文件或者用“emboss”来指定EMBOSS工具箱的输出文件。读者可以在以下链接中找到所有支持的格式。```http://biopython.org/wiki/AlignIO`` <http://biopython.org/wiki/AlignIO>`__ 和 `online <http://biopython.org/DIST/docs/api/Bio.AlignIO-module.html>`__:

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> help(AlignIO)
    ...

6.1.2  多个序列排列
~~~~~~~~~~~~~~~~~~~~~~~~~~

在前一章中，我们旨在读取的文件仅包含有一个序列排列。然而，在很多情况下，文件可能包含有多个序列排列。这时，你可以使用 ``Bio.AlignIO.parse()`` 来读取它们。

假设我们有一个PHYLIP格式的很小的序列排列：

.. code:: verbatim

        5    6
    Alpha     AACAAC
    Beta      AACCCC
    Gamma     ACCAAC
    Delta     CCACCA
    Epsilon   CCAAAC

如果你想用PHYLIP工具包来bootstrap一个系统发生树，其中的一个步骤是用 ``bootseq`` 程序来产生许多序列排列。这将给出类似于以下格式的序列排列：

.. code:: verbatim

        5     6
    Alpha     AAACCA
    Beta      AAACCC
    Gamma     ACCCCA
    Delta     CCCAAC
    Epsilon   CCCAAA
        5     6
    Alpha     AAACAA
    Beta      AAACCC
    Gamma     ACCCAA
    Delta     CCCACC
    Epsilon   CCCAAA
        5     6
    Alpha     AAAAAC
    Beta      AAACCC
    Gamma     AACAAC
    Delta     CCCCCA
    Epsilon   CCCAAC
    ...
        5     6
    Alpha     AAAACC
    Beta      ACCCCC
    Gamma     AAAACC
    Delta     CCCCAA
    Epsilon   CAAACC

如果你想用 ``Bio.AlignIO`` 来读取这个文件，你可以使用：

.. code:: verbatim

    from Bio import AlignIO
    alignments = AlignIO.parse("resampled.phy", "phylip")
    for alignment in alignments:
        print alignment
        print

这将给出以下的输出（这时只显示缩略的一部分）：

.. code:: verbatim

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAACCA Alpha
    AAACCC Beta
    ACCCCA Gamma
    CCCAAC Delta
    CCCAAA Epsilon

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAACAA Alpha
    AAACCC Beta
    ACCCAA Gamma
    CCCACC Delta
    CCCAAA Epsilon

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAAAAC Alpha
    AAACCC Beta
    AACAAC Gamma
    CCCCCA Delta
    CCCAAC Epsilon

    ...

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAAACC Alpha
    ACCCCC Beta
    AAAACC Gamma
    CCCCAA Delta
    CAAACC Epsilon

与 ``Bio.SeqIO.parse`` 一样， ``Bio.SeqIO.parse()`` 将返回一个迭代器（iterator）。如果你希望把所有的序列排列都读取到内存中，以下代码将把它们储存在一个列表对象里。

.. code:: verbatim

    from Bio import AlignIO
    alignments = list(AlignIO.parse("resampled.phy", "phylip"))
    last_align = alignments[-1]
    first_align = alignments[0]

6.1.3  含糊的序列排列
~~~~~~~~~~~~~~~~~~~~~~~~~~~

许多序列排列的文件格式可以非常明确地储存多个序列排列。然而，例如FASTA一类的普通序列文件格式并没有很直接的分隔符来分开多个序列排列。读者可以见以下例子：

.. code:: verbatim

    >Alpha
    ACTACGACTAGCTCAG--G
    >Beta
    ACTACCGCTAGCTCAGAAG
    >Gamma
    ACTACGGCTAGCACAGAAG
    >Alpha
    ACTACGACTAGCTCAGG--
    >Beta
    ACTACCGCTAGCTCAGAAG
    >Gamma
    ACTACGGCTAGCACAGAAG

以上FASTA格式文件可以认为是一个包含有6条序列的序列排列（有重复序列名）。或者从文件名来看，这很可能是两个序列排列，每一个包含有三个序列，只是这两个序列排列恰好具有相同的长度。

以下是另一个例子：

.. code:: verbatim

    >Alpha
    ACTACGACTAGCTCAG--G
    >Beta
    ACTACCGCTAGCTCAGAAG
    >Alpha
    ACTACGACTAGCTCAGG--
    >Gamma
    ACTACGGCTAGCACAGAAG
    >Alpha
    ACTACGACTAGCTCAGG--
    >Delta
    ACTACGGCTAGCACAGAAG

同样，这也可能是一个包含有六个序列的序列排列。然而，根据序列名判断，这很可能是三个两两间的序列比较，而且恰好有同样的长度。

最后一个例子也类似：

.. code:: verbatim

    >Alpha
    ACTACGACTAGCTCAG--G
    >XXX
    ACTACCGCTAGCTCAGAAG
    >Alpha
    ACTACGACTAGCTCAGG
    >YYY
    ACTACGGCAAGCACAGG
    >Alpha
    --ACTACGAC--TAGCTCAGG
    >ZZZ
    GGACTACGACAATAGCTCAGG

在这一个例子中，由于序列有不同的长度，这不能被当作是一个包含六个序列的单独的序列排列。很显然，这可以被看成是三个两两间的序列排列。

很明显，将多个序列排列以FASTA格式储存并不方便。然而，在某些情况下，如果你一定要这么做， ``Bio.AlignIO`` 依然能够处理上述情形（但是所有的序列排列必须都含有相同的序列）。一个很常见的例子是，我们经常会使用EMBOSS工具箱中的 ``needle`` 和 ``water`` 来产生许多两两间的序列排列（尽管在这种情况下，你可以指定数据格式为“emboss”给 ``Bio.AlignIO`` ）。

为了处理这样的FASTA格式的数据，我们可以指定 ``Bio.AlignIO.parse()`` 的第三个可选参数 ``seq_count`` ，这一参数将告诉Biopython你所期望的每个序列排列中序列的个数。例如：

.. code:: verbatim

    for alignment in AlignIO.parse(handle, "fasta", seq_count=2):
        print "Alignment length %i" % alignment.get_alignment_length()
        for record in alignment:
            print "%s - %s" % (record.seq, record.id)
        print

这将给出：

.. code:: verbatim

    Alignment length 19
    ACTACGACTAGCTCAG--G - Alpha
    ACTACCGCTAGCTCAGAAG - XXX

    Alignment length 17
    ACTACGACTAGCTCAGG - Alpha
    ACTACGGCAAGCACAGG - YYY

    Alignment length 21
    --ACTACGAC--TAGCTCAGG - Alpha
    GGACTACGACAATAGCTCAGG - ZZZ

如果你使用 ``Bio.AlignIO.read()`` 或者 ``Bio.AlignIO.parse()`` 而不指定 ``seq_count`` ，这将返回一个包含有六条序列的序列排列。对于上面的第三个例子，由于序列长度不同，Biopython将会报告一个错误。

如果数据格式本身包含有分割符， ``Bio.AlignIO`` 可以很聪明地自动确定文件中每一个序列排列而无需指定 ``seq_count`` 选项。如果你仍然指定 ``seq_count`` 但是却与数据本身的分隔符相冲突，Biopython也将报告一个错误。

注意指定这一可选的 ``seq_count`` 参数将假设文件中所有的序列排列都包含相同数目的序列。假如你真的遇到每一个序列排列都有不同数目的序列， ``Bio.AlignIO`` 将无法读取。这时，我们建议你使用 ``Bio.SeqIO`` 来读取数据，然后将序列转化为序列排列。

6.2  序列排列数据的写出
-----------------------

我们已经讨论了 ``Bio.AlignIO.read()`` 和 ``Bio.AlignIO.parse()`` 来读取各种格式的序列排列，现在让我们来使用 ``Bio.AlignIO.write()`` 写出序列排列文件。

这一函数接受三个参数：一个 ``MultipleSeqAlignment`` 对象（或者是一个 ``Alignment`` 对象），一个可写的文件句柄（handle）或者期望写出的文件名，以及写出文件的格式。

这里有一个手动构造一个 ``MultipleSeqAlignment`` 对象的例子（注意 ``MultipleSeqAlignment`` 是由若干个 ``SeqRecord`` 组成的）：

.. code:: verbatim

    from Bio.Alphabet import generic_dna
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    align1 = MultipleSeqAlignment([
                 SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
                 SeqRecord(Seq("ACT-CTAGCTAG", generic_dna), id="Beta"),
                 SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),
             ])

    align2 = MultipleSeqAlignment([
                 SeqRecord(Seq("GTCAGC-AG", generic_dna), id="Delta"),
                 SeqRecord(Seq("GACAGCTAG", generic_dna), id="Epsilon"),
                 SeqRecord(Seq("GTCAGCTAG", generic_dna), id="Zeta"),
             ])

    align3 = MultipleSeqAlignment([
                 SeqRecord(Seq("ACTAGTACAGCTG", generic_dna), id="Eta"),
                 SeqRecord(Seq("ACTAGTACAGCT-", generic_dna), id="Theta"),
                 SeqRecord(Seq("-CTACTACAGGTG", generic_dna), id="Iota"),
             ])

    my_alignments = [align1, align2, align3]

现在我们有一个包含三个 ``MultipleSeqAlignment`` 对象的列表（ ``my_alignments`` ），现在我们将它写出为PHYLIP格式：

.. code:: verbatim

    from Bio import AlignIO
    AlignIO.write(my_alignments, "my_example.phy", "phylip")

如果你用你喜欢的文本编辑器在你当前的工作目录下找到 ``my_example.phy`` 文件，你会看到以下内容：

.. code:: verbatim

     3 12
    Alpha      ACTGCTAGCT AG
    Beta       ACT-CTAGCT AG
    Gamma      ACTGCTAGDT AG
     3 9
    Delta      GTCAGC-AG
    Epislon    GACAGCTAG
    Zeta       GTCAGCTAG
     3 13
    Eta        ACTAGTACAG CTG
    Theta      ACTAGTACAG CT-
    Iota       -CTACTACAG GTG

在更多情况下，你希望读取一个已经含有序列排列的文件，经过某些操作（例如去掉一些行和列）然后将它重新储存起来。

假如你希望知道有多少序列排列被 ``Bio.AlignIO.write()`` 函数写入句柄中。如果你的序列排列都被放在一个列表中（如同以上的例子），你可以很容易地使用 ``len(my_alignments)`` 来获得这一信息。然而，如果你的序列排列在一个迭代器对象中，你无法轻松地完成这件事情。为此， ``Bio.AlignIO.write()`` 将会返回它所写出的序列排列个数。

*注意* - 如果你所指定给 ``Bio.AlignIO.write()`` 的文件已经存在在当前目录下，这一文件将被直接覆盖掉而不会有任何警告。

6.2.1  序列排列的格式间转换
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Bio.AlignIO`` 模块中的序列排列格式转化功能与 ``Bio.SeqIO`` （见 `5.5.2 <#sec:SeqIO-conversion>`__ ）模块的格式转化是一样的。在通常情况下，我们建议使用 ``Bio.AlignIO.parse()`` 来读取序列排列数据，然后使用 ``Bio.AlignIO.write()`` 函数来写出。或者你也可以直接使用 ``Bio.AlignIO.convert()`` 函数来实现格式的转换。

在本例中，我们将读取PFAM/Stockholm格式的序列排列，然后将其保存为Clustal格式。

.. code:: verbatim

    from Bio import AlignIO
    count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
    print "Converted %i alignments" % count

或者，使用 ``Bio.AlignIO.parse()`` 和 ``Bio.AlignIO.write()`` ：

.. code:: verbatim

    from Bio import AlignIO
    alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
    count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
    print "Converted %i alignments" % count

``Bio.AlignIO.write()`` 函数默认处理的情形是一个包括有多个序列排列的对象。在以上例子中，我们给予 ``Bio.AlignIO.write()`` 的参数是一个由 ``Bio.AlignIO.parse()`` 函数返回的一个迭代器。

在以下例子中，我们知道序列排列文件中仅包含有一个序列排列，因此我们使用 ``Bio.AlignIO.read()`` 函数来读取数据，然后使用 ``Bio.AlignIO.write()`` 来将保存数据保存为另一种格式。

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    AlignIO.write([alignment], "PF05371_seed.aln", "clustal")

使用以上两个例子，你都可以将PFAM/Stockholm格式的序列排列数据转化为Clustal格式。

.. code:: verbatim

    CLUSTAL X (1.81) multiple sequence alignment


    COATB_BPIKE/30-81                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS
    Q9T0Q8_BPIKE/1-52                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS
    COATB_BPI22/32-83                   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS
    COATB_BPM13/24-72                   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
    COATB_BPZJ2/1-49                    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFAS
    Q9T0Q9_BPFD/1-49                    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
    COATB_BPIF1/22-73                   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVS

    COATB_BPIKE/30-81                   KA
    Q9T0Q8_BPIKE/1-52                   RA
    COATB_BPI22/32-83                   KA
    COATB_BPM13/24-72                   KA
    COATB_BPZJ2/1-49                    KA
    Q9T0Q9_BPFD/1-49                    KA
    COATB_BPIF1/22-73                   RA

另外，你也可以使用以下代码将它保存为PHYLIP格式。

.. code:: verbatim

    from Bio import AlignIO
    AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")

你可以获得以下PHYLIP格式的文件输出：

.. code:: verbatim

     7 52
    COATB_BPIK AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
    Q9T0Q8_BPI AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
    COATB_BPI2 DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
    COATB_BPM1 AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
    COATB_BPZJ AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
    Q9T0Q9_BPF AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
    COATB_BPIF FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

               KA
               RA
               KA
               KA
               KA
               KA
               RA

PHYLIP格式最大的一个缺陷就是它严格地要求每一条序列的ID是都为10个字符（ID中多出的字符将被截短）。在这一个例子中，截短的序列ID依然是唯一的（只是缺少了可读性）。在某些情况下，我们并没有一个好的方式去压缩序列的ID。以下例子提供了另一种解决方案 —— 利用自定义的序列ID来代替原本的序列ID。

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    name_mapping = {}
    for i, record in enumerate(alignment):
        name_mapping[i] = record.id
        record.id = "seq%i" % i
    print name_mapping

    AlignIO.write([alignment], "PF05371_seed.phy", "phylip")

以上代码将会建立一个字典对象实现自定义的ID和原始ID的映射。

.. code:: verbatim

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', ...}

以下为PHYLIP的格式输出：

.. code:: verbatim

     7 52
    seq0       AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
    seq1       AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
    seq2       DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
    seq3       AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
    seq4       AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
    seq5       AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
    seq6       FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

               KA
               RA
               KA
               KA
               KA
               KA
               RA

由于序列ID的限制性，PHYLIP格式不是储存序列排列的理想格式。我们建议你将数据储存成PFAM/Stockholm或者其它能对序列排列进行注释的格式来保存你的数据。

6.2.2  Getting your alignment objects as formatted strings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Bio.AlignIO`` interface is based on handles, which means if you
want to get your alignment(s) into a string in a particular file format
you need to do a little bit more work (see below). However, you will
probably prefer to take advantage of the alignment object’s ``format()``
method. This takes a single mandatory argument, a lower case string
which is supported by ``Bio.AlignIO`` as an output format. For example:

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    print alignment.format("clustal")

As described in Section \ `4.5 <#sec:SeqRecord-format>`__), the
``SeqRecord`` object has a similar method using output formats supported
by ``Bio.SeqIO``.

Internally the ``format()`` method is using the ``StringIO`` string
based handle and calling ``Bio.AlignIO.write()``. You can do this in
your own code if for example you are using an older version of
Biopython:

.. code:: verbatim

    from Bio import AlignIO
    from StringIO import StringIO

    alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")

    out_handle = StringIO()
    AlignIO.write(alignments, out_handle, "clustal")
    clustal_data = out_handle.getvalue()

    print clustal_data

6.3  Manipulating Alignments
----------------------------

Now that we’ve covered loading and saving alignments, we’ll look at what
else you can do with them.

6.3.1  Slicing alignments
~~~~~~~~~~~~~~~~~~~~~~~~~

First of all, in some senses the alignment objects act like a Python
``list`` of ``SeqRecord`` objects (the rows). With this model in mind
hopefully the actions of ``len()`` (the number of rows) and iteration
(each row as a ``SeqRecord``) make sense:

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    >>> print "Number of rows: %i" % len(alignment)
    Number of rows: 7
    >>> for record in alignment:
    ...     print "%s - %s" % (record.seq, record.id)
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73

You can also use the list-like ``append`` and ``extend`` methods to add
more rows to the alignment (as ``SeqRecord`` objects). Keeping the list
metaphor in mind, simple slicing of the alignment should also make sense
- it selects some of the rows giving back another alignment object:

.. code:: verbatim

    >>> print alignment
    SingleLetterAlphabet() alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73
    >>> print alignment[3:7]
    SingleLetterAlphabet() alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73

What if you wanted to select by column? Those of you who have used the
NumPy matrix or array objects won’t be surprised at this - you use a
double index.

.. code:: verbatim

    >>> print alignment[2,6]
    T

Using two integer indices pulls out a single letter, short hand for
this:

.. code:: verbatim

    >>> print alignment[2].seq[6]
    T

You can pull out a single column as a string like this:

.. code:: verbatim

    >>> print alignment[:,6]
    TTT---T

You can also select a range of columns. For example, to pick out those
same three rows we extracted earlier, but take just their first six
columns:

.. code:: verbatim

    >>> print alignment[3:6,:6]
    SingleLetterAlphabet() alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49

Leaving the first index as ``:`` means take all the rows:

.. code:: verbatim

    >>> print alignment[:,:6]
    SingleLetterAlphabet() alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73

This brings us to a neat way to remove a section. Notice columns 7, 8
and 9 which are gaps in three of the seven sequences:

.. code:: verbatim

    >>> print alignment[:,6:9]
    SingleLetterAlphabet() alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73

Again, you can slice to get everything after the ninth column:

.. code:: verbatim

    >>> print alignment[:,9:]
    SingleLetterAlphabet() alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73

Now, the interesting thing is that addition of alignment objects works
by column. This lets you do this as a way to remove a block of columns:

.. code:: verbatim

    >>> edited = alignment[:,:6] + alignment[:,9:]
    >>> print edited
    SingleLetterAlphabet() alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73

Another common use of alignment addition would be to combine alignments
for several different genes into a meta-alignment. Watch out though -
the identifiers need to match up (see
Section \ `4.7 <#sec:SeqRecord-addition>`__ for how adding ``SeqRecord``
objects works). You may find it helpful to first sort the alignment rows
alphabetically by id:

.. code:: verbatim

    >>> edited.sort()
    >>> print edited
    SingleLetterAlphabet() alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49

Note that you can only add two alignments together if they have the same
number of rows.

6.3.2  Alignments as arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on what you are doing, it can be more useful to turn the
alignment object into an array of letters – and you can do this with
NumPy:

.. code:: verbatim

    >>> import numpy as np
    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    >>> align_array = np.array([list(rec) for rec in alignment], np.character)
    >>> align_array.shape
    (7, 52)

If you will be working heavily with the columns, you can tell NumPy to
store the array by column (as in Fortran) rather then its default of by
row (as in C):

.. code:: verbatim

    >>> align_array = np.array([list(rec) for rec in alignment], np.character, order="F")

Note that this leaves the original Biopython alignment object and the
NumPy array in memory as separate objects - editing one will not update
the other!

6.4  Alignment Tools
--------------------

There are *lots* of algorithms out there for aligning sequences, both
pairwise alignments and multiple sequence alignments. These calculations
are relatively slow, and you generally wouldn’t want to write such an
algorithm in Python. Instead, you can use Biopython to invoke a command
line tool on your behalf. Normally you would:

#. Prepare an input file of your unaligned sequences, typically this
   will be a FASTA file which you might create using ``Bio.SeqIO`` (see
   Chapter \ `5 <#chapter:Bio.SeqIO>`__).
#. Call the command line tool to process this input file, typically via
   one of Biopython’s command line wrappers (which we’ll discuss here).
#. Read the output from the tool, i.e. your aligned sequences, typically
   using ``Bio.AlignIO`` (see earlier in this chapter).

All the command line wrappers we’re going to talk about in this chapter
follow the same style. You create a command line object specifying the
options (e.g. the input filename and the output filename), then invoke
this command line via a Python operating system call (e.g. using the
``subprocess`` module).

Most of these wrappers are defined in the ``Bio.Align.Applications``
module:

.. code:: verbatim

    >>> import Bio.Align.Applications
    >>> dir(Bio.Align.Applications)
    ...
    ['ClustalwCommandline', 'DialignCommandline', 'MafftCommandline', 'MuscleCommandline',
    'PrankCommandline', 'ProbconsCommandline', 'TCoffeeCommandline' ...]

(Ignore the entries starting with an underscore – these have special
meaning in Python.) The module ``Bio.Emboss.Applications`` has wrappers
for some of the `EMBOSS suite <http://emboss.sourceforge.net/>`__,
including ``needle`` and ``water``, which are described below in
Section \ `6.4.5 <#seq:emboss-needle-water>`__, and wrappers for the
EMBOSS packaged versions of the PHYLIP tools (which EMBOSS refer to as
one of their EMBASSY packages - third party tools with an EMBOSS style
interface). We won’t explore all these alignment tools here in the
section, just a sample, but the same principles apply.

6.4.1  ClustalW
~~~~~~~~~~~~~~~

ClustalW is a popular command line tool for multiple sequence alignment
(there is also a graphical interface called ClustalX). Biopython’s
``Bio.Align.Applications`` module has a wrapper for this alignment tool
(and several others).

Before trying to use ClustalW from within Python, you should first try
running the ClustalW tool yourself by hand at the command line, to
familiarise yourself the other options. You’ll find the Biopython
wrapper is very faithful to the actual command line API:

.. code:: verbatim

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> help(ClustalwCommandline)
    ...

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <http://biopython.org/DIST/docs/tutorial/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). This is a small FASTA file containing seven prickly-pear
DNA sequences (from the cactus family *Opuntia*).

By default ClustalW will generate an alignment and guide tree file with
names based on the input FASTA file, in this case ``opuntia.aln`` and
``opuntia.dnd``, but you can override this or make it explicit:

.. code:: verbatim

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
    >>> print cline
    clustalw2 -infile=opuntia.fasta

Notice here we have given the executable name as ``clustalw2``,
indicating we have version two installed, which has a different filename
to version one (``clustalw``, the default). Fortunately both versions
support the same set of arguments at the command line (and indeed,
should be functionally identical).

You may find that even though you have ClustalW installed, the above
command doesn’t work – you may get a message about “command not found”
(especially on Windows). This indicated that the ClustalW executable is
not on your PATH (an environment variable, a list of directories to be
searched). You can either update your PATH setting to include the
location of your copy of ClustalW tools (how you do this will depend on
your OS), or simply type in the full path of the tool. For example:

.. code:: verbatim

    >>> import os
    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> clustalw_exe = r"C:\Program Files\new clustal\clustalw2.exe"
    >>> clustalw_cline = ClustalwCommandline(clustalw_exe, infile="opuntia.fasta")

.. code:: verbatim

    >>> assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
    >>> stdout, stderr = clustalw_cline()

Remember, in Python strings ``\n`` and ``\t`` are by default interpreted
as a new line and a tab – which is why we’re put a letter “r” at the
start for a raw string that isn’t translated in this way. This is
generally good practice when specifying a Windows style file name.

Internally this uses the ``subprocess`` module which is now the
recommended way to run another program in Python. This replaces older
options like the ``os.system()`` and the ``os.popen*`` functions.

Now, at this point it helps to know about how command line tools “work”.
When you run a tool at the command line, it will often print text output
directly to screen. This text can be captured or redirected, via two
“pipes”, called standard output (the normal results) and standard error
(for error messages and debug messages). There is also standard input,
which is any text fed into the tool. These names get shortened to stdin,
stdout and stderr. When the tool finishes, it has a return code (an
integer), which by convention is zero for success.

When you run the command line tool like this via the Biopython wrapper,
it will wait for it to finish, and check the return code. If this is non
zero (indicating an error), an exception is raised. The wrapper then
returns two strings, stdout and stderr.

In the case of ClustalW, when run at the command line all the important
output is written directly to the output files. Everything normally
printed to screen while you wait (via stdout or stderr) is boring and
can be ignored (assuming it worked).

What we care about are the two output files, the alignment and the guide
tree. We didn’t tell ClustalW what filenames to use, but it defaults to
picking names based on the input file. In this case the output should be
in the file ``opuntia.aln``. You should be able to work out how to read
in the alignment using ``Bio.AlignIO`` by now:

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("opuntia.aln", "clustal")
    >>> print align
    SingleLetterAlphabet() alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191

In case you are interested (and this is an aside from the main thrust of
this chapter), the ``opuntia.dnd`` file ClustalW creates is just a
standard Newick tree file, and ``Bio.Phylo`` can parse these:

.. code:: verbatim

    >>> from Bio import Phylo
    >>> tree = Phylo.read("opuntia.dnd", "newick")
    >>> Phylo.draw_ascii(tree)
                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658

Chapter `13 <#sec:Phylo>`__ covers Biopython’s support for phylogenetic
trees in more depth.

6.4.2  MUSCLE
~~~~~~~~~~~~~

MUSCLE is a more recent multiple sequence alignment tool than ClustalW,
and Biopython also has a wrapper for it under the
``Bio.Align.Applications`` module. As before, we recommend you try using
MUSCLE from the command line before trying it from within Python, as the
Biopython wrapper is very faithful to the actual command line API:

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> help(MuscleCommandline)
    ...

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <http://biopython.org/DIST/docs/tutorial/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). You can then tell MUSCLE to read in this FASTA file, and
write the alignment to an output file:

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.txt")
    >>> print cline
    muscle -in opuntia.fasta -out opuntia.txt

Note that MUSCLE uses “-in” and “-out” but in Biopython we have to use
“input” and “out” as the keyword arguments or property names. This is
because “in” is a reserved word in Python.

By default MUSCLE will output the alignment as a FASTA file (using
gapped sequences). The ``Bio.AlignIO`` module should be able to read
this alignment using ``format="fasta"``. You can also ask for
ClustalW-like output:

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clw=True)
    >>> print cline
    muscle -in opuntia.fasta -out opuntia.aln -clw

Or, strict ClustalW output where the original ClustalW header line is
used for maximum compatibility:

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clwstrict=True)
    >>> print cline
    muscle -in opuntia.fasta -out opuntia.aln -clwstrict

The ``Bio.AlignIO`` module should be able to read these alignments using
``format="clustal"``.

MUSCLE can also output in GCG MSF format (using the ``msf`` argument),
but Biopython can’t currently parse that, or using HTML which would give
a human readable web page (not suitable for parsing).

You can also set the other optional parameters, for example the maximum
number of iterations. See the built in help for details.

You would then run MUSCLE command line string as described above for
ClustalW, and parse the output using ``Bio.AlignIO`` to get an alignment
object.

6.4.3  MUSCLE using stdout
~~~~~~~~~~~~~~~~~~~~~~~~~~

Using a MUSCLE command line as in the examples above will write the
alignment to a file. This means there will be no important information
written to the standard out (stdout) or standard error (stderr) handles.
However, by default MUSCLE will write the alignment to standard output
(stdout). We can take advantage of this to avoid having a temporary
output file! For example:

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
    >>> print muscle_cline
    muscle -in opuntia.fasta

If we run this via the wrapper, we get back the output as a string. In
order to parse this we can use ``StringIO`` to turn it into a handle.
Remember that MUSCLE defaults to using FASTA as the output format:

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
    >>> stdout, stderr = muscle_cline()
    >>> from StringIO import StringIO
    >>> from Bio import AlignIO
    >>> align = AlignIO.read(StringIO(stdout), "fasta")
    >>> print align
    SingleLetterAlphabet() alignment with 7 rows and 906 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191663
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191665
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191664
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191661
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191660
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191659
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191658

The above approach is fairly simple, but if you are dealing with very
large output text the fact that all of stdout and stderr is loaded into
memory as a string can be a potential drawback. Using the ``subprocess``
module we can work directly with handles instead:

.. code:: verbatim

    >>> import subprocess
    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
    >>> child = subprocess.Popen(str(muscle_cline),
    ...                          stdout=subprocess.PIPE,
    ...                          stderr=subprocess.PIPE,
    ...                          shell=(sys.platform!="win32"))
    >>> from Bio import AlignIO
    >>> align = AlignIO.read(child.stdout, "fasta")
    >>> print align
    SingleLetterAlphabet() alignment with 7 rows and 906 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191663
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191665
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191664
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191661
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191660
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191659
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191658

6.4.4  MUSCLE using stdin and stdout
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We don’t actually *need* to have our FASTA input sequences prepared in a
file, because by default MUSCLE will read in the input sequence from
standard input! Note this is a bit more advanced and fiddly, so don’t
bother with this technique unless you need to.

First, we’ll need some unaligned sequences in memory as ``SeqRecord``
objects. For this demonstration I’m going to use a filtered version of
the original FASTA file (using a generator expression), taking just six
of the seven sequences:

.. code:: verbatim

    >>> from Bio import SeqIO
    >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)

Then we create the MUSCLE command line, leaving the input and output to
their defaults (stdin and stdout). I’m also going to ask for strict
ClustalW format as for the output.

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(clwstrict=True)
    >>> print muscle_cline
    muscle -clwstrict

Now for the fiddly bits using the ``subprocess`` module, stdin and
stdout:

.. code:: verbatim

    >>> import subprocess
    >>> import sys
    >>> child = subprocess.Popen(str(cline),
    ...                          stdin=subprocess.PIPE,
    ...                          stdout=subprocess.PIPE,
    ...                          stderr=subprocess.PIPE,
    ...                          shell=(sys.platform!="win32"))                     

That should start MUSCLE, but it will be sitting waiting for its FASTA
input sequences, which we must supply via its stdin handle:

.. code:: verbatim

    >>> SeqIO.write(records, child.stdin, "fasta")
    6
    >>> child.stdin.close()

After writing the six sequences to the handle, MUSCLE will still be
waiting to see if that is all the FASTA sequences or not – so we must
signal that this is all the input data by closing the handle. At that
point MUSCLE should start to run, and we can ask for the output:

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> align = AlignIO.read(child.stdout, "clustal")
    >>> print align
    SingleLetterAlphabet() alignment with 6 rows and 900 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF19166
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF19166
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF19165
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF19165

Wow! There we are with a new alignment of just the six records, without
having created a temporary FASTA input file, or a temporary alignment
output file. However, a word of caution: Dealing with errors with this
style of calling external programs is much more complicated. It also
becomes far harder to diagnose problems, because you can’t try running
MUSCLE manually outside of Biopython (because you don’t have the input
file to supply). There can also be subtle cross platform issues (e.g.
Windows versus Linux), and how you run your script can have an impact
(e.g. at the command line, from IDLE or an IDE, or as a GUI script).
These are all generic Python issues though, and not specific to
Biopython.

If you find working directly with ``subprocess`` like this scary, there
is an alternative. If you execute the tool with ``muscle_cline()`` you
can supply any standard input as a big string,
``muscle_cline(stdin=...)``. So, provided your data isn’t very big, you
can prepare the FASTA input in memory as a string using ``StringIO``
(see Section \ `22.1 <#sec:appendix-handles>`__):

.. code:: verbatim

    >>> from Bio import SeqIO
    >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)
    >>> from StringIO import StringIO
    >>> handle = StringIO()
    >>> SeqIO.write(records, handle, "fasta")
    6
    >>> data = handle.getvalue()

You can then run the tool and parse the alignment as follows:

.. code:: verbatim

    >>> stdout, stderr = muscle_cline(stdin=data)
    >>> from Bio import AlignIO
    >>> align = AlignIO.read(StringIO(stdout), "clustal")
    >>> print align
    SingleLetterAlphabet() alignment with 6 rows and 900 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF19166
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF19166
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF19165
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF19165

You might find this easier, but it does require more memory (RAM) for
the strings used for the input FASTA and output Clustal formatted data.

6.4.5  EMBOSS needle and water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `EMBOSS <http://emboss.sourceforge.net/>`__ suite includes the
``water`` and ``needle`` tools for Smith-Waterman algorithm local
alignment, and Needleman-Wunsch global alignment. The tools share the
same style interface, so switching between the two is trivial – we’ll
just use ``needle`` here.

Suppose you want to do a global pairwise alignment between two
sequences, prepared in FASTA format as follows:

.. code:: verbatim

    >HBA_HUMAN
    MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
    KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
    AVHASLDKFLASVSTVLTSKYR

in a file ``alpha.fasta``, and secondly in a file ``beta.fasta``:

.. code:: verbatim

    >HBB_HUMAN
    MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
    VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
    KEFTPPVQAAYQKVVAGVANALAHKYH

Let’s start by creating a complete ``needle`` command line object in one
go:

.. code:: verbatim

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline(asequence="alpha.faa", bsequence="beta.faa",
    ...                                  gapopen=10, gapextend=0.5, outfile="needle.txt")
    >>> print needle_cline
    needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5

Why not try running this by hand at the command prompt? You should see
it does a pairwise comparison and records the output in the file
``needle.txt`` (in the default EMBOSS alignment file format).

Even if you have EMBOSS installed, running this command may not work –
you might get a message about “command not found” (especially on
Windows). This probably means that the EMBOSS tools are not on your PATH
environment variable. You can either update your PATH setting, or simply
tell Biopython the full path to the tool, for example:

.. code:: verbatim

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline(r"C:\EMBOSS\needle.exe",
    ...                                  asequence="alpha.faa", bsequence="beta.faa",
    ...                                  gapopen=10, gapextend=0.5, outfile="needle.txt")

Remember in Python that for a default string ``\n`` or ``\t`` means a
new line or a tab – which is why we’re put a letter “r” at the start for
a raw string.

At this point it might help to try running the EMBOSS tools yourself by
hand at the command line, to familiarise yourself the other options and
compare them to the Biopython help text:

.. code:: verbatim

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> help(NeedleCommandline)
    ...

Note that you can also specify (or change or look at) the settings like
this:

.. code:: verbatim

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline()
    >>> needle_cline.asequence="alpha.faa"
    >>> needle_cline.bsequence="beta.faa"
    >>> needle_cline.gapopen=10
    >>> needle_cline.gapextend=0.5
    >>> needle_cline.outfile="needle.txt"
    >>> print needle_cline
    needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5
    >>> print needle_cline.outfile
    needle.txt

Next we want to use Python to run this command for us. As explained
above, for full control, we recommend you use the built in Python
``subprocess`` module, but for simple usage the wrapper object usually
suffices:

.. code:: verbatim

    >>> stdout, stderr = needle_cline()
    >>> print stdout + stderr
    Needleman-Wunsch global alignment of two sequences

Next we can load the output file with ``Bio.AlignIO`` as discussed
earlier in this chapter, as the ``emboss`` format:

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("needle.txt", "emboss")
    >>> print align
    SingleLetterAlphabet() alignment with 2 rows and 149 columns
    MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR HBA_HUMAN
    MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRF...KYH HBB_HUMAN

In this example, we told EMBOSS to write the output to a file, but you
*can* tell it to write the output to stdout instead (useful if you don’t
want a temporary output file to get rid of – use ``stdout=True`` rather
than the ``outfile`` argument), and also to read *one* of the one of the
inputs from stdin (e.g. ``asequence="stdin"``, much like in the MUSCLE
example in the section above).

This has only scratched the surface of what you can do with ``needle``
and ``water``. One useful trick is that the second file can contain
multiple sequences (say five), and then EMBOSS will do five pairwise
alignments.

Note - Biopython includes its own pairwise alignment code in the
``Bio.pairwise2`` module (written in C for speed, but with a pure Python
fallback available too). This doesn’t work with alignment objects, so we
have not covered it within this chapter. See the module’s docstring
(built in help) for details.
