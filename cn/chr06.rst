第6章 多序列比对
==============================================

多序列比对（Multiple Sequence Alignment, MSA）是指对多个序列进行对位排列。 这通常需要保证序列间的等同位点处在同一列上，并通过引进小横线（-）以保证最终的序列具有相同的长度。这种序列比对可以视作是由字符组成的矩阵。在Biopython中，多序列比对中每一个序列是以 ``SeqRecord`` 对象来表示的。

这里我们介绍一种新的对象 -- ``MultipleSeqAlignment`` 来表示这样一类数据，我们还将介绍 ``Bio.AlignIO`` 模块来读写不同格式的多序列比对数据（ ``Bio.AlignIO`` 在设计上与之前介绍的 ``Bio.SeqIO`` 模块是类似的）。Biopython中， ``Bio.SeqIO`` 和 ``Bio.AlignIO`` 都能读写各种格式的多序列比对数据。在实际处理中，使用哪一个模块取决于用户需要对数据进行何种操作。

本章的最后一部分是关于各种常用多序列比对程序（ClustalW和MUSCLE）的Biopython命令行封装。

6.1 读取多序列比对数据
-------------------------------------------

在Biopython中，有两种方法读取多序列比对数据， ``Bio.AlignIO.read()`` 和 ``Bio.AlignIO.parse()`` 。这两种方法跟 ``Bio.SeqIO`` 处理一个和多个数据的设计方式是一样的。 ``Bio.AlignIO.read()`` 只能读取一个多序列比对而 ``Bio.AlignIO.parse()`` 可以依次读取多个序列比对数据。 

使用 ``Bio.AlignIO.parse()`` 将会返回一个 ``MultipleSeqAlignment`` 的 *迭代器（iterator）* 。迭代器往往在循环中使用。在实际数据分析过程中会时常处理包含有多个多序列比对的文件。例如PHYLIP中的 ``seqboot`` ，EMBOSS工具箱中的 ``water`` 和 ``needle``, 以及Bill Pearson的FASTA程序。

然而在大多数情况下，你所遇到的文件仅仅包括一个多序列比对。这时，你应该使用 ``Bio.AlignIO.read()`` ，这将返回一个 ``MultipleSeqAlignment`` 对象。

这两个函数都接受两个必须参数：

#. 第一个参数为包含有多序列比对数据的 *句柄（handle）* 。在实际操作中，这往往是一个打开的文件（详细信息请见 `22.1 <#sec:appendix-handles>`__ ）或者文件名。

#. 第二个参数为多序列比对文件格式（小写）。与 ``Bio.SeqIO`` 模块一样，Biopython不会对将读取的文件格式进行猜测。所有 ``Bio.AlignIO`` 模块支持的多序列比对数据格式可以在 `http://biopython.org/wiki/AlignIO <http://biopython.org/wiki/AlignIO>`__ 中找到。

``Bio.AlignIO`` 模块还接受一个可选参数 ``seq_count`` 。这一参数将在 `6.1.3 <#sec:AlignIO-count-argument>`__ 中具体讨论。它可以处理不确定的多序列比对格式，或者包含有多个序列的排列。

另一个可选参数 ``alphabet`` 允许用户指定序列比对文件的字符（alphabet），它可以用来说明序列比对的类型（DNA，RNA或蛋白质）。因为大多数序列比对格式并不区别序列的类型，因此指定这一参数可能会对后期的分析产生帮助。 ``Bio.AlignIO`` 默认将使用一般字符（generic alphabet），这将不区分各种序列比对类型。

6.1.1 单一的序列比对
~~~~~~~~~~~~~~~~~~~~~~~~

例如，请见以下PFAM（或者Stockholm）格式的蛋白序列比对文件。

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

这是PFAM数据库中Phage\_Coat\_Gp8的种子排列（PF05371）。该排列下载于一个已经过期的PFAM数据库版本（ `http://pfam.sanger.ac.uk/ <http://pfam.sanger.ac.uk/>`__ ），但这并不影响我们的例子（假设你已经将以上内容下载到一个名为''PF05371\_seed.sth''的文件中，并在Python的当前工作目录下）：

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")

这段代码将在屏幕上打印出该序列比对的概要信息：

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

你会注意到，以上输出截短了中间一部分序列的内容。你也可以很容易地通过控制多序列比对中每一条序列（作为 ``SeqRecord`` 对象）来输出你所喜欢的格式。例如：

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

你是否已经注意到以上原始数据文件中包含有蛋白数据库（PDB）交叉引用以及相关二级结构的信息？你可以尝试以下代码：

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
`http://pfam.sanger.ac.uk/family?acc=PF05371 <http://pfam.sanger.ac.uk/family?acc=PF05371>`__
可以让你下载各种不同的序列比对的格式。以下例子为FASTA格式：

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

注意Sanger网站有一个选项可以将序列比对中的间隔（gap）用小圆点或者是小横线表示。在以上例子中，序列间隔由小横线表示。假设你已经下载该文件，并保存为 “PF05371\_seed.faa”。你可以使用以下代码来读入该序列比对。

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.faa", "fasta")
    print alignment

你可能已经发现，以上代码中唯一的变化只是指定格式的参数。所返回的alignment对象将会包含同样的序列和序列名字。但是仔细的读者会发现，每一个alignment的SeqRecord中并不包含数据的引用注释。这是因为FASTA格式本身并没有包含这一类信息。

此外，除了使用Sanger网站，你也可以利用 ``Bio.AlignIO`` 来将原始的Stockholm格式转换成FASTA文件格式（见下文）。

对于任何一种Biopython支持的格式，你都可以用同样的方式读取它（通过指定文件的格式）。例如，你可以使用“phylip”来表示PHYLIP格式文件，用"nexus"来指定NEXUS格式文件或者用“emboss”来指定EMBOSS工具箱的输出文件。读者可以在以下链接中找到所有支持的格式（ `http://biopython.org/wiki/AlignIO <http://biopython.org/wiki/AlignIO>`__ ），或者内置的帮助中（以及在线文档 `online <http://biopython.org/DIST/docs/api/Bio.AlignIO-module.html>`__ ）：

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> help(AlignIO)
    ...

6.1.2  多个序列比对
~~~~~~~~~~~~~~~~~~~~~~~~~~

在前一章中，我们旨在读取仅包含有一个序列比对的文件。然而，在很多情况下，文件可能包含有多个序列比对。这时，你可以使用 ``Bio.AlignIO.parse()`` 来读取它们。

假设我们有一个PHYLIP格式的很小的序列比对：

.. code:: verbatim

        5    6
    Alpha     AACAAC
    Beta      AACCCC
    Gamma     ACCAAC
    Delta     CCACCA
    Epsilon   CCAAAC

如果你想用PHYLIP工具包来bootstrap一个系统发生树，其中的一个步骤是用 ``bootseq`` 程序来产生许多序列比对。这将给出类似于以下格式的序列比对：

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

与 ``Bio.SeqIO.parse`` 一样， ``Bio.SeqIO.parse()`` 将返回一个迭代器（iterator）。如果你希望把所有的序列比对都读取到内存中，以下代码将把它们储存在一个列表对象里。

.. code:: verbatim

    from Bio import AlignIO
    alignments = list(AlignIO.parse("resampled.phy", "phylip"))
    last_align = alignments[-1]
    first_align = alignments[0]

6.1.3  含糊的序列比对
~~~~~~~~~~~~~~~~~~~~~~~~~~~

许多序列比对的文件格式可以非常明确地储存多个序列比对。然而，例如FASTA一类的普通序列文件格式并没有很直接的分隔符来分开多个序列比对。读者可以见以下例子：

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

以上FASTA格式文件可以认为是一个包含有6条序列的序列比对（有重复序列名）。或者从文件名来看，这很可能是两个序列比对，每一个包含有三个序列，只是这两个序列比对恰好具有相同的长度。

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

同样，这也可能是一个包含有六个序列的序列比对。然而，根据序列名判断，这很可能是三个两两间的序列比较，而且恰好有同样的长度。

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

在这一个例子中，由于序列有不同的长度，这不能被当作是一个包含六个序列的单独的序列比对。很显然，这可以被看成是三个两两间的序列比对。

很明显，将多个序列比对以FASTA格式储存并不方便。然而，在某些情况下，如果你一定要这么做， ``Bio.AlignIO`` 依然能够处理上述情形（但是所有的序列比对必须都含有相同的序列）。一个很常见的例子是，我们经常会使用EMBOSS工具箱中的 ``needle`` 和 ``water`` 来产生许多两两间的序列比对 —— 然而在这种情况下，你可以指定数据格式为“emboss”，``Bio.AlignIO`` 仍然能够识别这些原始输出。

为了处理这样的FASTA格式的数据，我们可以指定 ``Bio.AlignIO.parse()`` 的第三个可选参数 ``seq_count`` ，这一参数将告诉Biopython你所期望的每个序列比对中序列的个数。例如：

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

如果你使用 ``Bio.AlignIO.read()`` 或者 ``Bio.AlignIO.parse()`` 而不指定 ``seq_count`` ，这将返回一个包含有六条序列的序列比对。对于上面的第三个例子，由于序列长度不同，导致它们不能被解析为一个序列比对，Biopython将会抛出一个异常。

如果数据格式本身包含有分割符， ``Bio.AlignIO`` 可以很聪明地自动确定文件中每一个序列比对，而无需指定 ``seq_count`` 选项。如果你仍然指定 ``seq_count`` 但是却与数据本身的分隔符相冲突，Biopython将产生一个错误。

注意指定这一可选的 ``seq_count`` 参数将假设文件中所有的序列比对都包含相同数目的序列。假如你真的遇到每一个序列比对都有不同数目的序列， ``Bio.AlignIO`` 将无法读取。这时，我们建议你使用 ``Bio.SeqIO`` 来读取数据，然后将序列转换为序列比对。

6.2  序列比对的写出
-----------------------

我们已经讨论了 ``Bio.AlignIO.read()`` 和 ``Bio.AlignIO.parse()`` 来读取各种格式的序列比对，现在让我们来使用 ``Bio.AlignIO.write()`` 写出序列比对文件。

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

如果你用你喜欢的文本编辑器在你当前的工作目录下打开 ``my_example.phy`` 文件，你会看到以下内容：

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

在更多情况下，你希望读取一个已经含有序列比对的文件，经过某些操作（例如去掉一些行和列）然后将它重新储存起来。

假如你希望知道有多少序列比对被 ``Bio.AlignIO.write()`` 函数写入句柄中。如果你的序列比对都被放在一个列表中（如同以上的例子），你可以很容易地使用 ``len(my_alignments)`` 来获得这一信息。然而，如果你的序列比对在一个生成器/迭代器对象中，你无法轻松地完成这件事情。为此， ``Bio.AlignIO.write()`` 将会返回它所写出的序列比对个数。

*注意* - 如果你所指定给 ``Bio.AlignIO.write()`` 的文件已经存在在当前目录下，这一文件将被直接覆盖掉而不会有任何警告。

6.2.1  序列比对的格式间转换
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Bio.AlignIO`` 模块中的序列比对格式转换功能与 ``Bio.SeqIO`` （见 `5.5.2 <#sec:SeqIO-conversion>`__ ）模块的格式转换是一样的。在通常情况下，我们建议使用 ``Bio.AlignIO.parse()`` 来读取序列比对数据，然后使用 ``Bio.AlignIO.write()`` 函数来写出。或者你也可以直接使用 ``Bio.AlignIO.convert()`` 函数来实现格式的转换。

在本例中，我们将读取PFAM/Stockholm格式的序列比对，然后将其保存为Clustal格式：

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

``Bio.AlignIO.write()`` 函数默认处理的情形是一个包括有多个序列比对的对象。在以上例子中，我们给予 ``Bio.AlignIO.write()`` 的参数是一个由 ``Bio.AlignIO.parse()`` 函数返回的一个迭代器。

在以下例子中，我们知道序列比对文件中仅包含有一个序列比对，因此我们使用 ``Bio.AlignIO.read()`` 函数来读取数据，然后使用 ``Bio.AlignIO.write()`` 来将数据保存为另一种格式：

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    AlignIO.write([alignment], "PF05371_seed.aln", "clustal")

使用以上两个例子，你都可以将PFAM/Stockholm格式的序列比对数据转换为Clustal格式：

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

另外，你也可以使用以下代码将它保存为PHYLIP格式：

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

PHYLIP格式最大的一个缺陷就是它严格地要求每一条序列的ID是都为10个字符（ID中多出的字符将被截短）。在这一个例子中，截短的序列ID依然是唯一的（只是缺少了可读性）。在某些情况下，我们并没有一个好的方式去压缩序列的ID。以下例子提供了另一种解决方案 —— 利用自定义的序列ID来代替原本的序列ID：

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    name_mapping = {}
    for i, record in enumerate(alignment):
        name_mapping[i] = record.id
        record.id = "seq%i" % i
    print name_mapping

    AlignIO.write([alignment], "PF05371_seed.phy", "phylip")

以上代码将会建立一个字典对象实现自定义的ID和原始ID的映射：

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

由于序列ID的限制性，PHYLIP格式不是储存序列比对的理想格式。我们建议你将数据储存成PFAM/Stockholm或者其它能对序列比对进行注释的格式来保存你的数据。

6.2.2  将序列比对对象转换为格式化字符串（formatted strings）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

因为 ``Bio.AlignIO`` 模块是基于文件句柄的，因此你如果想将序列比对读入为一个字符串对象，你需要做一些额外的工作。然而，我们提供一个 ``format()`` 方法来帮助你实现这项任务。 ``format()`` 方法需要用户提供一个小写的格式参数（这可以是任何 ``AlignIO`` 支持的序列比对格式）。例如：

.. code:: verbatim

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    print alignment.format("clustal")

我们在 `4.5 <#sec:SeqRecord-format>`__ 中讲到， ``Bio.SeqIO`` 也有一个对 ``SeqRecord`` 输出的方法。

``format()`` 方法是利用 ``StringIO`` 以及 ``Bio.AlignIO.write()`` 来实现以上输出的。如果你使用的是较老版本的Biopython，你可以使用以下代码来完成相同的工作：

.. code:: verbatim

    from Bio import AlignIO
    from StringIO import StringIO

    alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")

    out_handle = StringIO()
    AlignIO.write(alignments, out_handle, "clustal")
    clustal_data = out_handle.getvalue()

    print clustal_data

6.3  序列比对的操纵
-------------------

现在我们已经了解了如何读入和写出序列比对。让我们继续看看如何对读入的序列比对进行操作。

6.3.1  序列比对的切片（slice）操作
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

首先，用户可以认为读入的序列比对是一个由 ``SeqRecord`` 对象构成的Python列表（list）。有了这样一个印象以后，你可以使用 ``len()`` 方法来得到行数（序列比对的个数），你也可以对序列比对进行迭代。

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

你可以使用列表所拥有的 ``append`` 和 ``extend`` 方法来给序列比对增加序列。请读者一定要正确理解序列比对与其包含的序列的关系，这样你就可以使用切片操作来获得其中某些序列比对。

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

假如你需要获得特定的列该怎么办呢？如果你接触过Numpy矩阵那么一定对下面的语法非常熟悉，使用双切片：

.. code:: verbatim

    >>> print alignment[2,6]
    T

使用两个整数来获得序列比对中的一个字符，这其实是以下操作的简化方式：

.. code:: verbatim

    >>> print alignment[2].seq[6]
    T

你可以用下面的代码来获取整列：

.. code:: verbatim

    >>> print alignment[:,6]
    TTT---T

你也可以同时选择特定的行和列。例如，以下代码将打印出第3到6行的前6列：

.. code:: verbatim

    >>> print alignment[3:6,:6]
    SingleLetterAlphabet() alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49

使用 ``:`` 将打印出所有行：

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

切片给我们提供了一个简单的方式来去除一部分序列比对。在以下例子中，有三条序列的7，8，9三列为间隔（-）。

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

你也可以通过切片来获得第9列以后的所有序列：

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

现在，你可以通过列来操纵序列比对。这也是你能够去除序列比对中的许多列。例如：

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

另一个经常使用的序列比对操作是将多个基因的序列比对拼接成一个大的序列比对（meta-alignment）。
在进行这种操作时一定要注意序列的ID需要匹配（具体请见 `4.7 <#sec:SeqRecord-addition>`__ 关于 ``SeqRecord``
的说明)。为了达到这种目的，用 ``sort()`` 方法将序列ID按照字母顺序进行排列可能会有所帮助：

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

注意：只有当两个序列比对拥有相同的行的时候才能进行序列比对的拼接。

6.3.2  序列比对作为数组
~~~~~~~~~~~~~~~~~~~~~~~~~~~

根据你的需要，有时将序列比对转换为字符数组是非常方便的。你可以用 ``Numpy`` 来实现这一目的：

.. code:: verbatim

    >>> import numpy as np
    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    >>> align_array = np.array([list(rec) for rec in alignment], np.character)
    >>> align_array.shape
    (7, 52)

如果你需要频繁地使用列操作，你可以让 ``Numpy`` 将序列比对以列的形式进行储存（与Fortran一样），而不是 ``Numpy`` 默认形式（与C一样以行储存）：

.. code:: verbatim

    >>> align_array = np.array([list(rec) for rec in alignment], np.character, order="F")

注意， ``Numpy`` 的数组和Biopython默认的序列比对对象是分别储存在内存中的，编辑其中的一个不会更新另一个的值。

6.4  构建序列比对的工具
--------------------

目前有非常多的算法来帮助你构建一个序列比对，包括两两间的比对和多序列比对。这些算法在计算上往往是非常慢的，你一定不会希望用Python来实现他们。然而，你可以使用Biopython来运行命令行程序。通常你需要：

#. 准备一个包含未比对序列的输入文件，一般为FASTA格式的序列。你可以使用 ``Bio.SeqIO`` 来创建一个 (具体见第5章 <#chapter:Bio.SeqIO>`__).
#. 在Biopython中运行一个命令行程序来构建序列比对（我们将在这里详细介绍）。这需要通过Biopython的打包程序（wrapper）来实现。
#. 读取以上程序的输出，也就是排列好的序列比对。这往往可以通过 ``Bio.AlignIO`` 来实现（请看本章前部分内容）。

本章所介绍的所有的命令行打包程序都将以同样的方式使用。你创造一个命令行对象来指定各种参数（例如：输入文件名，输出文件名等），然后通过Python的系统命令模块来运行这一程序（例如：使用 ``subprocess`` 进程）。

大多数的打包程序都在 ``Bio.Align.Applications`` 中定义：

.. code:: verbatim

    >>> import Bio.Align.Applications
    >>> dir(Bio.Align.Applications)
    ...
    ['ClustalwCommandline', 'DialignCommandline', 'MafftCommandline', 'MuscleCommandline',
    'PrankCommandline', 'ProbconsCommandline', 'TCoffeeCommandline' ...]

（以下划线开头的记录不是Biopython打包程序，这些变量在Python中有特殊的含义。） ``Bio.Emboss.Applications`` 中包含对 `EMBOSS  <http://emboss.sourceforge.net/>`__ 的打包程序（包括 ``needle`` 和 ``water`` ）。EMBOSS和PHYLIP的打包程序将在 `6.4.5 <#seq:emboss-needle-water>`__ 节中详细介绍。在本章中，我们并不打算将所有的序列比对程序都予以介绍，但是Biopython中各种序列比对程序都具有相同的使用方式。

6.4.1  ClustalW
~~~~~~~~~~~~~~~

ClustalW是一个非常流行的进行多序列比对的命令行程序（其还有一个图形化的版本称之为ClustalX）。Biopython的 ``Bio.Align.Applications`` 模块包含这一多序列比对程序的打包程序。

我们建议你在Python中使用ClustalW之前在命令行界面下手动使用ClustalW，这样能使你更清楚这一程序的参数。你会发现Biopython打包程序非常严格地遵循实际的命令行API：

.. code:: verbatim

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> help(ClustalwCommandline)
    ...

作为最简单的一个例子，你仅仅需要一个FASTA格式的序列文件作为输入，例如： `opuntia.fasta <http://biopython.org/DIST/docs/tutorial/examples/opuntia.fasta>`__ （你可以在线或者在Biopython/Doc/examples文件夹中找到该序列）。 `opuntia.fasta` 包含着7个prickly-pear的DNA序列（来自仙人掌科）。

ClustalW在默认情况下会产生一个包括所有输入序列的序列比对以及一个由输入序列名字构成的指导树（guide tree）。例如，用上述文件作为输入，ClustalW将会输出 ``opuntia.aln`` 和 ``opuntia.dnd`` 两个文件：

.. code:: verbatim

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
    >>> print cline
    clustalw2 -infile=opuntia.fasta

注意这里我们给出的执行文件名是 ``clustalw2`` ，这是ClustalW的第二个版本（第一个版本的文件名为 ``clustalw`` ）。ClustalW的这两个版本具有相同的参数，并且在功能上也是一致的。

你可能会发现，尽管你安装了ClustalW，以上的命令行却无法正确运行。你可能会得到“command not found”的错误信息（尤其是在Windows上）。这往往是由于ClustalW的运行程序并不在系统的工作目录PATH下（一个包含着运行程序路径的环境变量）。你既可以修改PATH，使其包括ClustalW的运行程序（不同系统需要以不同的方式修改），或者你也可以直接指定程序的绝对路径。例如：

.. code:: verbatim

    >>> import os
    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> clustalw_exe = r"C:\Program Files\new clustal\clustalw2.exe"
    >>> clustalw_cline = ClustalwCommandline(clustalw_exe, infile="opuntia.fasta")

.. code:: verbatim

    >>> assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
    >>> stdout, stderr = clustalw_cline()

注意，Python中 ``\n`` 和 ``\t`` 会被解析为一个新行和制表空白（tab）。然而，如果你将一个小写的“r”放在字符串的前面，这一字符串就将保留原始状态，而不被解析。这种方式对于指定Windows风格的文件名来说是一种良好的习惯。

Biopython在内部使用较新的 ``subprocess`` 模块来实现打包程序，而不是 ``os.system()`` 和 ``os.popen*`` 。

现在，我们有必要去了解命令行工具是如何工作的。当你使用一个命令行时，它往往会在屏幕上输出一些内容。这一输出可以被保存或重定向。在系统输出中，有两种管道（pipe）来区分不同的输出信息--标准输出（standard output）包含正常的输出内容，标准错误（standard error）显示错误和调试信息。同时，系统也接受标准输入（standard input）。这也是命令行工具如何读取数据文件的。当程序运行结束以后，它往往会返回一个整数。一般返回值为0意味着程序正常结束。

当你使用Biopython打包程序来调用命令行工具的时候，它将会等待程序结束，并检查程序的返回值。如果返回值不为0，Biopython将会提示一个错误信息。Biopython打包程序将会输出两个字符串，标准输出和标准错误。

在ClustalW的例子中，当你使用程序时，所有重要的输出都被保存到输出文件中。所有打印在屏幕上的内容（通过 stdout or stderr）可以被忽略掉（假设它已经成功运行）。

当运行ClustalW的时候，我们所关心的往往是输出的序列比对文件和指导树文件。ClustalW会自动根据输入数据的文件名来命名输出文件。在本例中，输出文件将是 ``opuntia.aln`` 。当你成功运行完ClustalW以后，你可以使用 ``Bio.AlignIO`` 来读取输出结果：

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

另一个输出文件 ``opuntia.dnd`` 中包含有一个newick格式的指导树，你可以使用Biopython中的 ``Bio.Phylo`` 来读取它：

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

`13 <#sec:Phylo>`__  章中详细介绍了如何使用Biopython对进化树数据进行处理。

6.4.2  MUSCLE
~~~~~~~~~~~~~

MUSCLE是另一个较新的序列比对工具，Biopython的 ``Bio.Align.Applications`` 中也有针对Muscle的打包程序。与ClustalW一样，我们也建议你先在命令行界面下使用MUSCLE以后再使用Biopython打包程序。你会发现，Biopython的打包程序非常严格地包括了所有命令行输入参数：

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> help(MuscleCommandline)
    ...

作为最简单的例子，你只需要一个Fasta格式的数据文件作为输入。例如： `opuntia.fasta <http://biopython.org/DIST/docs/tutorial/examples/opuntia.fasta>`__ 然后你可以告诉MUSCLE来读取该FASTA文件，并将序列比对写出：

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.txt")
    >>> print cline
    muscle -in opuntia.fasta -out opuntia.txt

注意，MUSCLE使用“-in”和“-out”来指定输入和输出文件，而在Biopython中，我们使用“input”和“out”作为关键字来指定输入输出。这是由于“in”是Python的一个关键词而被保留。

默认情况下，MUSCLE的输出文件将是包含间隔（gap）的FASTA格式文件。 当你指定 ``format=fasta`` 时， ``Bio.AlignIO`` 能够读取该FASTA文件。你也可以告诉MUSCLE来输出ClustalW-like的文件结果：

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clw=True)
    >>> print cline
    muscle -in opuntia.fasta -out opuntia.aln -clw

或者，严格的ClustalW的输出文件（这将输出原始的ClustalW的文件标签）。例如：

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clwstrict=True)
    >>> print cline
    muscle -in opuntia.fasta -out opuntia.aln -clwstrict

你可以使用 ``Bio.AlignIO`` 的 ``format="clustal"`` 参数来读取这些序列比对输出。

MUSCLE也可以处理GCG和MSF（使用 ``msf`` 参数）甚至HTML格式，但是目前Biopython并不能读取它们。

你也可以设置MUSCLE其它的可选参数，例如最大数目的迭代数。具体信息请查阅Biopython的内部帮助文档。

6.4.3  MUSCLE标准输出
~~~~~~~~~~~~~~~~~~~~~~~~~~

使用以上的MUSCLE命令行将会把序列比对结果写出到一个文件中。然而MUSCLE也允许你将序列比对结果作为系统的标准输出。Biopython打包程序可以利用这一特性来避免创建一个临时文件。例如：

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
    >>> print muscle_cline
    muscle -in opuntia.fasta

如果你使用打包程序运行上述命令，程序将返回一个字符串对象。为了读取它，我们可以使用 ``StringIO`` 模块。记住MUSCLE将默认以FASTA格式输出序列比对：

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

以上是一个非常简单的例子，如果你希望处理较大的输出数据，我们并不建议你将它们全部读入内存中。对于这种情况， ``subprocess`` 模块可以非常方便地处理。例如：

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

6.4.4  以标准输入和标准输出使用MUSCLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

事实上，我们并不需要将序列放在一个文件里来使用MUSCLE。MUSCLE可以读取系统标准输入的内容。注意，这有一点高级和繁琐，若非必须，你可以不用关心这个技术。

为了让MUSCLE读取标准输入的内容，我们首先需要将未排列的序列以 ``SeqRecord`` 对象的形式读入到内存。在这里，我们将以一个规则来选择特定的序列（序列长度小于900bp的），使用生成器表达式。

.. code:: verbatim

    >>> from Bio import SeqIO
    >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)

随后，我们需要建立一个MUSCLE命令行，但是不指定输入和输出（MUSCLE默认为标准输入和标准输出）。这里，我们将指定输出格式为严格的Clustal格式：

.. code:: verbatim

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(clwstrict=True)
    >>> print muscle_cline
    muscle -clwstrict

我们使用Python的内置模块 ``subprocess`` 来实现这一目的：

.. code:: verbatim

    >>> import subprocess
    >>> import sys
    >>> child = subprocess.Popen(str(cline),
    ...                          stdin=subprocess.PIPE,
    ...                          stdout=subprocess.PIPE,
    ...                          stderr=subprocess.PIPE,
    ...                          shell=(sys.platform!="win32"))                     

这一命令将启动MUSCLE，但是它将会等待FASTA格式的输入数据。我们可以通过标准输入句柄来提供给它：

.. code:: verbatim

    >>> SeqIO.write(records, child.stdin, "fasta")
    6
    >>> child.stdin.close()

在将6条序列写入句柄后，MUSCLE仍将会等待，判断是否所有的FASTA序列全部输入完毕了。我们可以关闭句柄来提示给MUSCLE。这时，MUSCLE将开始运行。最后，我们可以在标准输出中获得结果：

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

现在我们在没有创造一个FASTA文件的情况下获得了一个序列比对。然而，由于你没有在Biopython外运行MUSCLE，这会使调试程序的难度增大，而且存在程序跨平台使用的问题（Windows和Linux）。

如果你觉得 ``subprocess`` 不方便使用，Biopython提供了另一种方式。如果你用 ``muscle_cline()`` 来运行外部程序（如MUSCLE），你可以用一个字符串对象作为输入。例如，你可以以这种方式使用： ``muscle_cline(stdin=...)`` 。假如你的序列文件不大，你可以将其储存为 ``StringIO`` 对象（具体见 `22.1 <#sec:appendix-handles>`__)：

.. code:: verbatim

    >>> from Bio import SeqIO
    >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)
    >>> from StringIO import StringIO
    >>> handle = StringIO()
    >>> SeqIO.write(records, handle, "fasta")
    6
    >>> data = handle.getvalue()

你可以以下方式运行外部程序和读取结果：

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

你可能觉得这种方式更便捷，但它需要更多的内存（这是由于我们是以字符串对象来储存输入的FASTA文件和输出的Clustal排列）。

6.4.5  EMBOSS包的序列比对工具——needle和water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


`EMBOSS <http://emboss.sourceforge.net/>`__ 包有两个序列比对程序—— ``water`` 和 ``needle`` 来实现Smith-Waterman做局部序列比对（local alignment）和Needleman-Wunsch算法来做全局排列（global alignment）。这两个程序具有相同的使用方式，因此我们仅以 ``needle`` 为例。

假设你希望做全局的序列两两排列，你可以将FASTA格式序列以如下方式储存：

.. code:: verbatim

    >HBA_HUMAN
    MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
    KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
    AVHASLDKFLASVSTVLTSKYR

以上内容在 ``alpha.fasta`` 文件中，另一个在 ``beta.fasta`` 中如下：

.. code:: verbatim

    >HBB_HUMAN
    MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
    VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
    KEFTPPVQAAYQKVVAGVANALAHKYH

让我们开始使用一个完整的 ``needle`` 命令行对象：

.. code:: verbatim

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline(asequence="alpha.faa", bsequence="beta.faa",
    ...                                  gapopen=10, gapextend=0.5, outfile="needle.txt")
    >>> print needle_cline
    needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5

你可能会有疑问，为什么不直接在终端里运行这一程序呢？你会发现，它将进行一个序列两两间的排列，并把结果记录在 ``needle.txt`` 中（以EMBOSS默认的序列比对格式）。

即使你安装了EMBOSS，使用以上命令仍可能会出错，你可能获得一个错误消息“command not found”，尤其是在Windows环境中。这很可能是由于EMBOSS工具的安装目录并不在系统的PATH中。遇到这种情况，你既可以更新系统的环境变量，也可以在Biopython中指定EMBOSS的安装路径。例如：

.. code:: verbatim

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline(r"C:\EMBOSS\needle.exe",
    ...                                  asequence="alpha.faa", bsequence="beta.faa",
    ...                                  gapopen=10, gapextend=0.5, outfile="needle.txt")

在Python中， ``\n`` 和 ``\t`` 分别意味着换行符和制表符。而在字符串前有一个“r”代表着raw字符串（ ``\n`` 和 ``\t`` 将保持它们本来的状态）。

现在你可以自己尝试着手动运行EMBOSS工具箱中的程序，比较一下各个参数以及其对应的Biopython打包程序帮助文档：

.. code:: verbatim

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> help(NeedleCommandline)
    ...

提示：你也可以指定特定的参数设置。例如：

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

现在我们获得了一个 ``needle`` 命令行，并希望在Python中运行它。我们在之前解释过，如果你希望完全地控制这一过程， ``subprocess`` 是最好的选择，但是如果你只是想尝试使用打包程序，以下命令足以达到目的：

.. code:: verbatim

    >>> stdout, stderr = needle_cline()
    >>> print stdout + stderr
    Needleman-Wunsch global alignment of two sequences

随后，我们需要载入 ``Bio.AlignIO`` 模块来读取needle输出（ ``emboss`` 格式）：

.. code:: verbatim

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("needle.txt", "emboss")
    >>> print align
    SingleLetterAlphabet() alignment with 2 rows and 149 columns
    MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR HBA_HUMAN
    MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRF...KYH HBB_HUMAN

在这个例子中，我们让EMBOSS将结果保存到一个输出文件中，但是你也可以让其写入标准输出中（这往往是在不需要临时文件的情况下的选择，你可以使用 ``stdout=True`` 参数而不是 ``outfile`` 参数）。与MUSCLE的例子一样，你也可以从标准输入里读取序列（ ``asequence="stdin"`` 参数）。

以上例子仅仅介绍了 ``needle`` 和 ``water`` 最简单的使用。一个有用的小技巧是，第二个序列文件可以包含有多个序列，EMBOSS工具将将每一个序列与第一个文件进行两两序列比对。

注意，Biopython有它自己的两两比对模块 ``Bio.pairwise2`` （用C语言编写）。但是它无法与序列比对对象一起工作，因此我们不在本章讨论它。具体信息请查阅模块的docstring（内部帮助文档）。
