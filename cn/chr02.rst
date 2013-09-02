第二章 快速开始 —— 你能用Biopython做什么？
========================================================

此部分旨在能让你快速开始Biopython，并给你一个大概的了解什么能用以及如何使用它。
此部分的所有例子都会假设你有个大致的Python使用知识，并且前提是你已经在你系统上
安装了Biopython。如果你认为你需要认真复习Python，主流的Python网站提供了相当多
的免费文档，你可以从以下网站开始(```http://www.python.org/doc/`` <http://www.python.org/doc/>`__)。

由于计算机上大量的生物学工作涉及到网上的数据库，某些例子也会需要联网操作才能
完成。

OK了，这就是所有的了，让我们看看我们能用Biopython做什么吧。

2.1  Biopython提供概览
------------------------------------------------

正如介绍中提到的，Biopython是一个库的集合，这个库能给在计算机上解决生物学家
感兴趣的事情。基本上这意味着你需要至少一点编程经验（当然是Python！）或至少
有兴趣学习编程。Biopython的任务就是通过提供可重复利用的库，让你作为程序员，
从而使得你的工作变得更加容易，以便于你能集中回答你所感兴趣的问题，而不必去
关注解析特殊文件格式的内部问题（当然，如果你想帮我们写一个原本不存在的解析
器并把它贡献给Biopython，请继续！）。所以Biopython的职责是让你更加轻松！

值得一提的是，Biopython通常能给出多种方式来解决“相同的事情”。在最近的版本中，
情况有所改善，但这仍可让人沮丧，因为在理想的Python中应该只有一种正确的方式
去解决问题。但是，这也可以成为一个真正的好处，因为它给了你很多灵活性和对库的
控制。本教程给你展示普通的或简单的方式去处理问题以便于你能自己处理事情。想要
学习更多替代的可能性，请查看Cookbook（Chapter `18 <#chapter:cookbook>`__,
这里有一些很酷的技巧和提示），进阶部分(Chapter `20 <#chapter:advanced>`__)，
内建“文档”（通过Python help命令），或者`API 文档 <http://biopython.org/DIST/docs/api/>`__)
或者最根本上代码本身。

2.2  处理序列
---------------------------

有讨论余地的（理所当然！），生物信息学的中心对象是序列。因此，我们先快速开始介绍
一下Biopython处理序列的机制，主要是``Seq``对象，这个我们也将会在第三章中详细讨论 
\ `3 <#chapter:Bio.Seq>`__。

大多数时候当我们想到一条序列时，在我们脑海中都会有一串类似于‘\ ``AGTACACTGGT``\ ’的
字母串。你可以按以下步骤创建一个``Seq``对象——“\ ``>>>``\表示Python提示符后紧跟你要
输入的内容：

.. code:: verbatim

    >>> from Bio.Seq import Seq
    >>> my_seq = Seq("AGTACACTGGT")
    >>> my_seq
    Seq('AGTACACTGGT', Alphabet())
    >>> print my_seq
    AGTACACTGGT
    >>> my_seq.alphabet
    Alphabet()

我们这里是一个带有*通用*字母的序列对象——说明我们还*没有*指定其为一条DNA还是蛋白质
序列（OK，一个蛋白质序列应该有很多的Ala，Gly，Cys和Thr！）。在第三章我们将讨论更
多关于字母的\ `3 <#chapter:Bio.Seq>`__。

除了有一个字母，``Seq``对象支持不同于Python的字符串方法。你不能对一个纯字符串做以
下处理：

.. code:: verbatim

    >>> my_seq
    Seq('AGTACACTGGT', Alphabet())
    >>> my_seq.complement()
    Seq('TCATGTGACCA', Alphabet())
    >>> my_seq.reverse_complement()
    Seq('ACCAGTGTACT', Alphabet())

另一个最重要的类是``SeqRecord``或Sequence Record。它保留了一条序列（作为``Seq``对象）
额外的注释信息，包括ID，name和description。用于读写序列文件格式的``Bio.SeqIO``模块
能与``SeqRecord``对象一起工作，稍后我们将会介绍，详细内容在第五章也有\ `5 <#chapter:Bio.SeqIO>`__。

这涵盖了基本的功能和Biopython序列类的使用。现在你应该有一些想法像怎么和Biopython库互动，
是时候去钻研它的乐趣，探讨处理生物学文件格式的有趣世界了！

2.3  用法示例
--------------------

在我们跳到语法解析器和Biopython处理的其他所有事之前，让我们先建立一个例子来激励我们
所做的每一件事情并让生活变得更加有趣。毕竟，如果没有任何生物学在这篇教程里，那你为
什么还要读它呢？

因为我喜欢植物，我想我们就来一个基于植物的例子吧（对其他生物的爱好者对不起了！）。
刚刚结束了去我们当地的一个温室的旅程，我们对Lady Slipper Orchids突然有一个难以置信
的观念（如果你想知道为什么，请看看一些Lady Slipper Orchids的照片<http://www.flickr.com/search/?q=lady+slipper+orchid&s=int&z=t>`__，
并试试Google图片搜索<http://images.google.com/images?q=lady%20slipper%20orchid>`__）。

当然，兰花不仅仅只有外观好看，它们也极力吸引着人们研究其进化和系统分类学。
假设我们正在考虑写一份关于Lady Slipper Orchids进化研究的基金方案，我们就会想
看看别人已经做了什么样的研究，然后我们能够增加一些什么内容。

经过一些阅读之后，我们发现Lady Slipper Orchids属于兰科拖鞋兰亚科并且由5个属组成：
*Cypripedium*，*Paphiopedilum*，*Phragmipedium*，*Selenipedium*和*Mexipedium*。

这已经给了我们足够多的信息来探究更多的东西。现在，让我们看看Biopython工具能起到怎样的作用。
我们从一条2.4部分解析出来的序列开始\ `2.4 <#sec:sequence-parsing>`__， 但是我们稍后还是回到
兰花上来——比如我们将在PubMed上搜索有关兰花的文章然后在第九章中提取GenBank上的序列\ `9 <#chapter:entrez>`__，
第十章中从Swiss-Prot上从特定的兰花蛋白中提取数据\ `10 <#chapter:swiss_prot>`__，最后在6.4.1
部分我们用ClustalW对兰花蛋白进行多序列比对\ `6.4.1 <#sec:align_clustal>`__。 

2.4  解析序列文件格式
----------------------------------

很多生物信息学工作的一大部分都会涉及到处理各种包含有生物学数据的文件格式类型。这些文件加载了
有趣的生物学数据，因而一个特殊的挑战是需要将这些文件解析成你能使用某种编程语言操作的格式。然而
这些解析工作有时会让人感到失望，因为这些格式有可能经常改变，而一个细微的改变也有可能打破甚至是
设计得最好的解析器。

我们现在开始简明的介绍``Bio.SeqIO``模块——你可以在第五章中查看更多\ `5 <#chapter:Bio.SeqIO>`__。
我们从在线搜索我们的朋友——Lady Slipper Orchids——开始。为保持这篇介绍简单一点，我们仅仅手动使用
NCBI网站。我们先看看NCBI上的nucleotide库，使用在线的Entrez搜索
(```http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?db=Nucleotide`` <http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?db=Nucleotide>`__)
包含Cypripedioideae所有东西（这是Lady Slipper Orchids的亚科）。

当本教程最初编写时，这个搜索仅给我们找到了94条匹配的信息，我们将结果保存为FASTA格式文本文件和
GenBank格式文本文件（文件```ls_orchid.fasta`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.fasta>`__
和```ls_orchid.gbk`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.gbk>`__，
也包含在Biopython源代码包下``docs/tutorial/examples/``）。

如果你今天搜索，你将会获得几百个的匹配结果，跟随者此教程，如果你想要看看包含基因的相同的列表，
请下载这两个文件或者从Biopython源代码中拷贝 ``docs/examples/`` 。在2.5部分\ `2.5 <#sec:connecting-with-biological-databases>`__
我们将会看到怎样使用Python做类似的搜索。

2.4.1  简单的FASTA解析示例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

如果你用你喜好的文本编辑器打开了lady slipper orchids的FASTA文件```ls_orchid.fasta`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.fasta>`__，
你会看到文件开头类似这样：

.. code:: verbatim

    >gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTG
    AATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGTGACCCTGATTTGTTGTTGGG
    ...

它包含有94条查询记录，每一行都以“\ ``>``\ ”开头，（大于符号）紧随其后的是一行或多行序列。现在
试试以下Python代码：

.. code:: verbatim

    from Bio import SeqIO
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print seq_record.id
        print repr(seq_record.seq)
        print len(seq_record)

你应该会得到类似这样的一些东西出现在屏幕上：

.. code:: verbatim

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC', SingleLetterAlphabet())
    740
    ...
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC', SingleLetterAlphabet())
    592

注意FASTA文件并没有指定字母表，因此``Bio.SeqIO``被默认为一般的``SingleLetterAlphabet()``而不是特定的DNA序列。

2.4.2  简单的GenBank解析示例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

现在我们来加载一个GenBank文件```ls_orchid.gbk`` <http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.gbk>`__
来代替——注意此次的代码与上面处理FASTA文件的一小段几乎完全相同——仅有的不同之处是我们改变了文件名和格式的字符串：

.. code:: verbatim

    from Bio import SeqIO
    for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
        print seq_record.id
        print repr(seq_record.seq)
        print len(seq_record)

这段代码应该会给出：

.. code:: verbatim

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC', IUPACAmbiguousDNA())
    740
    ...
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC', IUPACAmbiguousDNA())
    592

这一次``Bio.SeqIO``能够明智的选择一个字母表，IUPAC Ambiguous DNA。你也应该注意到了这个例子中有一个较短的字符串
被用作``seq_record.id``。

2.4.3  我爱解析——请别停止讨论它！
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Biopython有很多的解析器，基于它们所解析的文件格式，每一个都有自己独特的地位。第五章\ `5 <#chapter:Bio.SeqIO>`__ covers ``Bio.SeqIO``
将详细探讨，而第六章\ `6 <#chapter:Bio.AlignIO>`__将介绍用于序列比对的``Bio.AlignIO``。

由于最主流的文件格式都有解析器整合在``Bio.SeqIO``和/或 ``Bio.AlignIO``中，对于一些比较罕见的或者不被人们喜爱的
文件格式，要么根本就没有解析器，要么就是一些没有链接的老的解析器。请也到wiki页面```http://biopython.org/wiki/SeqIO`` <http://biopython.org/wiki/SeqIO>`__
以及```http://biopython.org/wiki/AlignIO`` <http://biopython.org/wiki/AlignIO>`__查看最新信息，或者咨询邮件列
表。wiki页面上应该包含了支持文件类型的最新列表，还有一些附加的例子。

另一个查找特定解析器信息和如何处理它们的地就是Cookbook（本教程的第十八章 `18 <#chapter:cookbook>`__）。如果你
没有找到你要的信息，请考虑及时帮帮你那可怜的过劳的文档机，并提交一份关于它的cookbook的entry号！（一旦你知道怎
么做了，那就是了！）

2.5  与生物学数据库结合
-----------------------------------------

在生物信息学中你需要做的很普通的事情之一是从生物学数据库中提取信息。手动访问这些数据库可能会非常枯燥乏味，尤其
是当你有很多重复的工作要做的时候。Biopython试图通过用Python脚本做一些可用的在线数据库来节省你的时间和精力。当前，
Biopython拥有从以下数据库中获取信息的代码：

-  NCBI的`Entrez <http://www.ncbi.nlm.nih.gov/Entrez/>`__ （和 `PubMed <http://www.ncbi.nlm.nih.gov/PubMed/>`__）
   ——见第九章\ `9 <#chapter:entrez>`__。
-  `ExPASy <http://www.expasy.org/>`__ ——见第十章\ `10 <#chapter:swiss_prot>`__。
-  `SCOP <http://scop.mrc-lmb.cam.ac.uk/scop/>`__——见``Bio.SCOP.search()``函数。

这些模块里的代码基本上使很容易地写出与这些页面中CGI脚本交互的Python代码，因此你能很方便的获得需要处理的格式。在某些
情况下，结果可能会紧密结合到Biopython解析器中从而使得提取信息更加简单。

2.6  下一步该做什么
--------------------

现在你已经做到这一步，你应该对基本的Biopython有一个很好的了解了，并准备好开始用它完成一些有用的工作。现在最好先完成
阅读本教程，然后如果你想开始看看源代码，并看看自动生成文档。

一旦你有一份你要做的工作的蓝图，以及Biopython能完成它的库，你应该看看Cookbook（第十八章`18 <#chapter:cookbook>`__），
在这里可能会有一些类似你工作的示例代码。

如果你知道你想要做什么，但是还没弄明白怎么去做，请随时将你的问题贴出到主要的Biopython列表中（见
```http://biopython.org/wiki/Mailing_lists`` <http://biopython.org/wiki/Mailing_lists>`__）。这不仅方便我们回答你的
问题，也有助于我们改进文档以便于它能帮到下一个做你工作的人。

请享受代码吧！