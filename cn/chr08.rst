第8章  BLAST和其他序列搜索工具(*实验性质的代码*)
======================================================================

*WARNING*: 这章教程介绍了Biopython中一个 *实验的* 模块。它正在被加入到
Biopython中，并且以一个预尾期的状态整理到教程当中，这样在我们发布稳定版
的之前可以收到一系列的反馈和并作改进。那时有些细节可能会改变，并且用到
当前 ``Bio.SearchIO`` 模块的脚步也需要更新。切记！为了与NCBI BLAST有关的
代码可以稳定工作，请继续使用第 \ `7 <#chapter:blast>`__ 章介绍的 Bio.Blast。

生物序列的鉴定是生物信息工作的主要部分。有几个工具像BLAST（可能是最流行
的），FASTA ,HMMER还有许多其它的都有这个功能，每个工具都有独特的算法和
途径。一般来说，这些工具都是用你的序列在数据库中搜索可能的匹配。随着序列
数量的增加（匹配数也会随之增加），将会有成百上千的可能匹配，解析这些结果
无疑变得越来越困难。理所当然，人工解析搜索结果变得不可能。而且你可能会同
时用几种不同的搜索工具，每种工具都有独特的统计方法、规则和输出格式。可以
想象，同时用多种工具搜索多条序列是多么恐怖的事。

我们对此非常了解，所以我们在Biopython创造了 ``Bio.SearchIO`` 亚模块。
``Bio.SearchIO`` 模块使从搜索结果中提取信息变得简单，并且可以处理不同工具
的不同标准和规则。``SearchIO`` 和BioPerl中模块名字一致。

在本章中，我们将学习 ``Bio.SearchIO`` 的主要功能，知道它可以为你做什么。我
们将使用两个主要的搜索工具：BLAST和FASTA。它们只是用来阐明思路，你可以轻
易地把工作流程应用到 ``Bio.SearchIO`` 支持的其他工具中。欢迎你使用我们将要
用到的搜索结果文件。BLAST搜索结果文件可以在
`here <http://biopython.org/SRC/Tests/Tutorial/my_blast.xml>`__ 下载。
BLAT输出结果文件可以在
`here <http://biopython.org/SRC/Tests/Tutorial/my_blat.psl>`__ 下载。两个结果
文件都是用下面这条序列搜索产生的：

.. code:: verbatim

    >mystery_seq
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

BLAST的XML结果是用 ``blastn`` 搜索NCBI的 ``refseq_rna`` 数据库得到的。对于
BLAT，数据库是2009年2月份的 ``hg19`` 人类基因组草图，输出格式是PSL。

我们从 ``Bio.SearchIO`` 的对象模型的介绍开始。这个模型代表你的搜索结果，因
此它是 ``Bio.SearchIO`` 的核心。然后，我们会介绍 ``Bio.SearchIO`` 常用的主要
功能。

现在一切就绪，让我们开始第一步：介绍核心对象模型。

8.1  SearchIO对象模型
------------------------------

尽管多数搜索工具的输出风格极为不同，但是它们包含的概念很相似：


-  输出文件可能包含一条或更多的搜索查询的结果。
-  在每次查询中，你会在给定的数据库中得到一个或更多的hits。
-  在每个hit中，你会得到一个或更多包含查询序列和数据库序列实际比对的区域。
-  一些工具如BLAT和Exonerate可能会把这些区域分成几个比对片段（或在BLAT中
   称为区块，在Exonerate中称为可能外显子）。这并不是很常见，像BLAST和
   HMMER就不这么做。

知道这些共性之后，我们决定它们作为创造 ``Bio.SearchIO`` 对象模型的基础。对
象模型包括一个Python对象的嵌套分层，每个都代表一个上面列出的概念。这些对
象是：

-  ``QueryResult``，代表单个查询序列。
-  ``Hit``，代表单个的数据库hit。``Hit`` 对象包含在 ``QueryResult`` 中，
   每个 ``QueryResult`` 中有0个或更多 ``Hit`` 对象。
-  ``HSP`` (high-scoring pair（高分片段）的缩写)，代表查询和匹配序列中有
   意义比对的区域。``HSP`` 对象包含在 ``Hit`` 对象中，而且每个 ``Hit`` 有一个
   或更多的 ``HSP`` 对象。   
-  ``HSPFragment``，代表查询和匹配序列中单个的邻近比对。 ``HSPFragment``
   对象包含在 ``HSP`` 对象中。多数的搜索工具如BLAST和HMMER把 ``HSP`` 和
   ``HSPFragment`` 合并，因为一个 ``HSP`` 只含有一个 ``HSPFragment``。但是
   像BLAT和Exonerate会产生含有多个 ``HSPFragment`` 的 ``HSP`` 。似乎有些困
   惑？不要紧，稍后我们将详细介绍这两个对象。

这四个对象是当你用 ``Bio.SearchIO`` 会碰到的。 ``Bio.SearchIO`` 四
个主要方法： ``read`` ，``parse``，``index``，or ``index_db`` 中任意一个都可
以产生这四个对象。这些方法的会在后面的部分详细介绍。这部分只会用到 ``read`` 和
``parse`` ，这两个方法和 ``Bio.SeqIO`` 以及 ``Bio.AlignIO`` 中的 ``read`` 和 ``parse`` 方法功
能相似：

-  ``read`` 用于搜索有单个查询的输出文件并且返回一个 ``QueryResult`` 对象。
-  ``parse`` 用于搜索有多个查询的输出文件并且返回一个可以yield
   ``QueryResult`` 对象的generator。

完成这些之后，让我们开始学习每个 ``Bio.SearchIO`` 对象，从 ``QueryResult``
开始。

8.1.1  QueryResult
~~~~~~~~~~~~~~~~~~

``QueryResult``，代表单个查询序列，每个 ``QueryResult`` 中有0个或更多 ``Hit``
对象。我们来看看BLAST文件时什么样的：

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> blast_qresult = SearchIO.read('my_blast.xml', 'blast-xml')
    >>> print blast_qresult
    Program: blastn (2.2.27+)
      Query: 42291 (61)
             mystery_seq
     Target: refseq_rna
       Hits: ----  -----  ----------------------------------------------------------
                #  # HSP  ID + description                                          
             ----  -----  ----------------------------------------------------------
                0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
                1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
                2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...
                3      2  gi|301171322|ref|NR_035857.1|  Pan troglodytes microRNA...
                4      1  gi|301171267|ref|NR_035851.1|  Pan troglodytes microRNA...
                5      2  gi|262205330|ref|NR_030198.1|  Homo sapiens microRNA 52...
                6      1  gi|262205302|ref|NR_030191.1|  Homo sapiens microRNA 51...
                7      1  gi|301171259|ref|NR_035850.1|  Pan troglodytes microRNA...
                8      1  gi|262205451|ref|NR_030222.1|  Homo sapiens microRNA 51...
                9      2  gi|301171447|ref|NR_035871.1|  Pan troglodytes microRNA...
               10      1  gi|301171276|ref|NR_035852.1|  Pan troglodytes microRNA...
               11      1  gi|262205290|ref|NR_030188.1|  Homo sapiens microRNA 51...
               12      1  gi|301171354|ref|NR_035860.1|  Pan troglodytes microRNA...
               13      1  gi|262205281|ref|NR_030186.1|  Homo sapiens microRNA 52...
               14      2  gi|262205298|ref|NR_030190.1|  Homo sapiens microRNA 52...
               15      1  gi|301171394|ref|NR_035865.1|  Pan troglodytes microRNA...
               16      1  gi|262205429|ref|NR_030218.1|  Homo sapiens microRNA 51...
               17      1  gi|262205423|ref|NR_030217.1|  Homo sapiens microRNA 52...
               18      1  gi|301171401|ref|NR_035866.1|  Pan troglodytes microRNA...
               19      1  gi|270133247|ref|NR_032574.1|  Macaca mulatta microRNA ...
               20      1  gi|262205309|ref|NR_030193.1|  Homo sapiens microRNA 52...
               21      2  gi|270132717|ref|NR_032716.1|  Macaca mulatta microRNA ...
               22      2  gi|301171437|ref|NR_035870.1|  Pan troglodytes microRNA...
               23      2  gi|270133306|ref|NR_032587.1|  Macaca mulatta microRNA ...
               24      2  gi|301171428|ref|NR_035869.1|  Pan troglodytes microRNA...
               25      1  gi|301171211|ref|NR_035845.1|  Pan troglodytes microRNA...
               26      2  gi|301171153|ref|NR_035838.1|  Pan troglodytes microRNA...
               27      2  gi|301171146|ref|NR_035837.1|  Pan troglodytes microRNA...
               28      2  gi|270133254|ref|NR_032575.1|  Macaca mulatta microRNA ...
               29      2  gi|262205445|ref|NR_030221.1|  Homo sapiens microRNA 51...
               ~~~
               97      1  gi|356517317|ref|XM_003527287.1|  PREDICTED: Glycine ma...
               98      1  gi|297814701|ref|XM_002875188.1|  Arabidopsis lyrata su...
               99      1  gi|397513516|ref|XM_003827011.1|  PREDICTED: Pan panisc...

虽然我们才接触对象模型的皮毛，但是你已经可以看到一些 有用的信息了。通过调用
``QueryResult`` 对象的 ``print`` 方法，你可以看到：

-  程序的名称和版本 (blastn version 2.2.27+)
-  查询的ID，描述和序列的长度(ID是42291，描述是 ‘mystery\_seq’，长度是61)
-  搜索的目标数据库 (refseq\_rna)
-  hits结果的快速预览。对于我们的查询序列，有100个可能的hits（表格中标记
   0-99）对于每个hit，我们可以看到它包含的高分比对片段（HSP)，ID和一个片
   段的描述。注意， ``Bio.SearchIO`` 截断了表格，只显示0-29，然后是97-99。
 
现在让我们用同样的步骤来检查BLAT的结果：

.. code:: verbatim

    >>> blat_qresult = SearchIO.read('my_blat.psl', 'blat-psl')
    >>> print blat_qresult
    Program: blat (<unknown version>)
      Query: mystery_seq (61)
             <unknown description>
     Target: <unknown target>
       Hits: ----  -----  ----------------------------------------------------------
                #  # HSP  ID + description                                          
             ----  -----  ----------------------------------------------------------
                0     17  chr19  <unknown description>                              

马上可以看到有些不同点。有些是由于BLAT使用PSL格式储存它的信息，稍后会看
到。其余是由于BLAST和BLAT搜索的程序和数据库之间明显的差异造成的：

-  程序名称和版本。 ``Bio.SearchIO`` 知道程序是BLAST，但是在输出文件中没
   有信息显示程序版本，所以默认是 ‘<unknown version>’。
-  查询的ID，描述和序列的长度。注意，这些细节和BLAST的细节只有细小的差别，
   ID是 ‘mystery\_seq’ 而不是42991，这是未知描述，但是序列长度仍是61。这
   实际上是文件格式本身导致的差异。BLAST有时创建自己的查询ID并且用你的原
   始ID作为序列描述。
-  目标数据库是未知的，因为BLAT输出文件没提到相关信息。
-  最后，hits列表完全不同，这里，我们的查询序列只hit到 ‘chr19’ 数据库条
   目，但是我们可以看到它含有17个HSP区域。这真是让人诧异，但是考虑到我们
   使用的是不同的程序，并且这些程序都有自己的数据库。

所有通过调用 ``print``方法看到的信息都可以单独地用Python的对象属性入
口标记获得（又叫点标记法）。同样还可以用相同的方法获得其他格式特有的属性。

.. code:: verbatim

    >>> print "%s %s" % (blast_qresult.program, blast_qresult.version)
    blastn 2.2.27+
    >>> print "%s %s" % (blat_qresult.program, blat_qresult.version)
    blat <unknown version>
    >>> blast_qresult.param_evalue_threshold    # blast-xml specific
    10.0

想获得一个可访问属性的完整列表，可以查询每个格式特有的文档。这些是 `for
BLAST <http://biopython.org/DIST/docs/api/Bio.SearchIO.BlastIO-module.html>`__
and for
`BLAT <http://biopython.org/DIST/docs/api/Bio.SearchIO.BlatIO-module.html>`__.

已经看到了在 ``QueryResult`` 对象上调用 ``print`` 方法，让我们研究的更深
一些。 ``QueryResult``到底是什么？就Python对象来说， ``QueryResult`` 混合
了列表和字典的特性。换句话说，也就是一个包含了列表和字典方便功能的容器对象。

和列表以及字典一样， ``QueryResult`` 对象是可迭代的。每次迭代返回一个hit
对象：

.. code:: verbatim

    >>> for hit in blast_qresult:
    ...     hit
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)
    Hit(id='gi|301171311|ref|NR_035856.1|', query_id='42291', 1 hsps)
    Hit(id='gi|270133242|ref|NR_032573.1|', query_id='42291', 1 hsps)
    Hit(id='gi|301171322|ref|NR_035857.1|', query_id='42291', 2 hsps)
    Hit(id='gi|301171267|ref|NR_035851.1|', query_id='42291', 1 hsps)
    ...

要得到 ``QueryResult`` 对象有多少条目(hits)，可以简单调用Python的 ``len`` 
方法：
.. code:: verbatim

    >>> len(blast_qresult)
    100
    >>> len(blat_qresult)
    1

同列表类似，你可以用切片来获得 ``QueryResult``对象的条目(hits)：

.. code:: verbatim

    >>> blast_qresult[0]        # retrieves the top hit
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)
    >>> blast_qresult[-1]       # retrieves the last hit
    Hit(id='gi|397513516|ref|XM_003827011.1|', query_id='42291', 1 hsps)

要得到多个条目，你同样可以对 ``QueryResult`` 对象作切片。这种情况下，切片
一个包含被切hits的新 ``QueryResult`` 对象：

.. code:: verbatim

    >>> blast_slice = blast_qresult[:3]     # slices the first three hits
    >>> print blast_slice
    Program: blastn (2.2.27+)
      Query: 42291 (61)
             mystery_seq
     Target: refseq_rna
       Hits: ----  -----  ----------------------------------------------------------
                #  # HSP  ID + description                                          
             ----  -----  ----------------------------------------------------------
                0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
                1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
                2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...

同字典类似，可以通过hit的ID获取hits。如果你知道一个特定的hit ID存在于一个
搜索结果中时，这将特别有用：

.. code:: verbatim

    >>> blast_qresult['gi|262205317|ref|NR_030195.1|']
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)

你可以用 ``hits`` 方法获得完整的 ``Hit`` 对象，也可以用 ``hit_keys``方法
获得完整的``Hit`` IDs：

.. code:: verbatim

    >>> blast_qresult.hits
    [...]       # list of all hits
    >>> blast_qresult.hit_keys
    [...]       # list of all hit IDs

如果你想确定一个特殊的hit是否存在于查询对象中该怎么做呢？可以用 ``in`` 
关键字作一个简单的成员检验：

.. code:: verbatim

    >>> 'gi|262205317|ref|NR_030195.1|' in blast_qresult
    True
    >>> 'gi|262205317|ref|NR_030194.1|' in blast_qresult
    False

有时候，只知道一个hit是否存在是不够的；你可能也会想知道hit的排名。 ``index`` 
方法可以帮助你：

.. code:: verbatim

    >>> blast_qresult.index('gi|301171437|ref|NR_035870.1|')
    22

记住，我们用的是Python风格的索引，是从0开始。这代表hit的排名是23而不是22。

同样，注意你看的hit排名是基于原始搜索输出文件的本来顺序。不同的搜索工具可
能会基于不同的标准排列这些hits。

如果原本的hit排序不合你意，可以用 ``QueryResult`` 对象的 ``sort`` 方法。
它和Python的 ``list.sort`` 方法很相似，只是有个是否创建一个新的排序后的
``QueryResult`` 对象的选项。

这里有个用 ``QueryResult.sort`` 方法排序hits的例子，这个方法基于每个hit
的完整序列长度。对于这个特殊的排序，我们设置 ``in_place`` 参数等于 ``False`` ，
这样排序方法会返回一个新的 ``QueryResult`` 对象，而原来的对象是未排序的。
我们同样可以设置 ``reverse`` 参数等于True以递减排序。

.. code:: verbatim

    >>> for hit in blast_qresult[:5]:   # id and sequence length of the first five hits
    ...     print hit.id, hit.seq_len
    ...
    gi|262205317|ref|NR_030195.1| 61
    gi|301171311|ref|NR_035856.1| 60
    gi|270133242|ref|NR_032573.1| 85
    gi|301171322|ref|NR_035857.1| 86
    gi|301171267|ref|NR_035851.1| 80

    >>> sort_key = lambda hit: hit.seq_len
    >>> sorted_qresult = blast_qresult.sort(key=sort_key, reverse=True, in_place=False)
    >>> for hit in sorted_qresult[:5]:
    ...     print hit.id, hit.seq_len
    ...
    gi|397513516|ref|XM_003827011.1| 6002
    gi|390332045|ref|XM_776818.2| 4082
    gi|390332043|ref|XM_003723358.1| 4079
    gi|356517317|ref|XM_003527287.1| 3251
    gi|356543101|ref|XM_003539954.1| 2936

有 ``in_place`` 参数的好处是可以保留原本的顺序，后面可能会用到。注意这不
是 ``QueryResult.sort`` 的默认行为，需要我们明确地设置 ``in_place`` 为True。

现在，你已经知道使用 ``QueryResult`` 对象。但是，在我们学习 ``Bio.SearchIO`` 
模块下个对象前，先了解下可以使 ``QueryResult`` 对象更易使用的两个方法：
``filter`` 和 ``map`` 方法。

如果你对Python的列表推导式、generator表达式或内建的 ``filter`` 和 ``map`` 
很熟悉，就知道（不知道就是看看吧!)它们在处理list-like的对象时有多有用。
你可以用这些内建的方法来操作 ``QueryResult`` 对象，这将止于常规的list，
并且你会丧失作更多有趣操作的能力。

这就是为什么 ``QueryResult`` 对象提供自己特有的 ``filter`` 和 ``map`` 
方法。对于 ``filter`` 有相似的 ``hit_filter`` 和 ``hsp_filter`` 方法，
从名称就可以看出，这些方法过滤 ``QueryResult`` 对象的 ``Hit`` 对象或者
``HSP`` 对象。同样的，对于 ``map`` ， ``QueryResult`` 对象同样提供相似
的  ``hit_map`` 和 ``hsp_map`` 方法。这些方法分别应用于 ``QueryResult`` 
对象的所有hits或者HSPs。 

让我们来看看这些方法的功能，从 ``hit_filter`` 开始。这个方法接受一个回调
函数，这个函数检验给定的 ``Hit`` 是否符合你设定的条件。换句话说，这个方法
必须接受一个单独 ``Hit`` 对象作为参数并且返回True或False。 

这里有个用 ``hit_filter`` 筛选出只有一个HSP的 ``Hit`` 对象的例子：

.. code:: verbatim

    >>> filter_func = lambda hit: len(hit.hsps) > 1     # the callback function
    >>> len(blast_qresult)      # no. of hits before filtering
    100
    >>> filtered_qresult = blast_qresult.hit_filter(filter_func)
    >>> len(filtered_qresult)   # no. of hits after filtering
    37
    >>> for hit in filtered_qresult[:5]:    # quick check for the hit lengths
    ...     print hit.id, len(hit.hsps)
    gi|301171322|ref|NR_035857.1| 2
    gi|262205330|ref|NR_030198.1| 2
    gi|301171447|ref|NR_035871.1| 2
    gi|262205298|ref|NR_030190.1| 2
    gi|270132717|ref|NR_032716.1| 2

``hsp_filter`` 和 ``hit_filter``功能相同，只是它过滤每个hit中的 ``HSP`` 对象，
而不是 ``Hit`` 。

对于 ``map`` 方法，同样接受一个回调函数作为参数。但是回调函数返回修改过的
 ``Hit`` 或 ``HSP``对象（取决于你是否使用 ``hit_map`` 或 ``hsp_map``方法），
 而不是返回 ``True`` 或 ``False``。

来看一个用 ``hit_map`` 方法来重命名hit ID的例子：

.. code:: verbatim

    >>> def map_func(hit):
    ...     hit.id = hit.id.split('|')[3]   # renames 'gi|301171322|ref|NR_035857.1|' to 'NR_035857.1'
    ...     return hit
    ...
    >>> mapped_qresult = blast_qresult.hit_map(map_func)
    >>> for hit in mapped_qresult[:5]:
    ...     print hit.id
    NR_030195.1
    NR_035856.1
    NR_032573.1
    NR_035857.1
    NR_035851.1

同样的， ``hsp_map`` 和 ``hit_map``作用相似, 但是作用于 ``HSP`` 对象而不
是 ``Hit`` 对象。

8.1.2  Hit
~~~~~~~~~~

``Hit`` 对象代表从单个数据库获得所有查询结果。在 ``Bio.SearchIO``对象等级
中是二级容器。它们被包含在 ``QueryResult``对象中，同时它们又包含 ``HSP`` 
对象。

看看它们是什么样的，从我们的BLAST搜索开始：

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> blast_qresult = SearchIO.read('my_blast.xml', 'blast-xml')
    >>> blast_hit = blast_qresult[3]    # fourth hit from the query result

.. code:: verbatim

    >>> print blast_hit
    Query: 42291
           mystery_seq
      Hit: gi|301171322|ref|NR_035857.1| (86)
           Pan troglodytes microRNA mir-520c (MIR520C), microRNA
     HSPs: ----  --------  ---------  ------  ---------------  ---------------------
              #   E-value  Bit score    Span      Query range              Hit range
           ----  --------  ---------  ------  ---------------  ---------------------
              0   8.9e-20     100.47      60           [1:61]                [13:73]
              1   3.3e-06      55.39      60           [0:60]                [13:73]

可以看到我们获得了必要的信息：

-  查询ID和描述信息。一个hit总是和一个查询绑定，所有我们同样希望记录原始
   查询。这些值可以通过 ``query_id`` 和  ``query_description`` 属性从hit
   中获取。
-  我们同样得到了hit ID、描述和序列全长。它们可以分别通过 ``id``，
   ``description``，和 ``seq_len`` 获取。
-  最后，有一个含有这个hit的HSPs的简短信息的表。在每行中，HSP重要信息被
   列出来：HSP索引，e值，得分，长度（包括gap），查询序列坐标和hit坐标。

现在，和BLAT结果作对比。记住，在BLAT搜索结果中，我们发现有一个含有17HSP的
hit。

.. code:: verbatim

    >>> blat_qresult = SearchIO.read('my_blat.psl', 'blat-psl')
    >>> blat_hit = blat_qresult[0]      # the only hit
    >>> print blat_hit
    Query: mystery_seq
           <unknown description>
      Hit: chr19 (59128983)
           <unknown description>
     HSPs: ----  --------  ---------  ------  ---------------  ---------------------
              #   E-value  Bit score    Span      Query range              Hit range
           ----  --------  ---------  ------  ---------------  ---------------------
              0         ?          ?       ?           [0:61]    [54204480:54204541]
              1         ?          ?       ?           [0:61]    [54233104:54264463]
              2         ?          ?       ?           [0:61]    [54254477:54260071]
              3         ?          ?       ?           [1:61]    [54210720:54210780]
              4         ?          ?       ?           [0:60]    [54198476:54198536]
              5         ?          ?       ?           [0:61]    [54265610:54265671]
              6         ?          ?       ?           [0:61]    [54238143:54240175]
              7         ?          ?       ?           [0:60]    [54189735:54189795]
              8         ?          ?       ?           [0:61]    [54185425:54185486]
              9         ?          ?       ?           [0:60]    [54197657:54197717]
             10         ?          ?       ?           [0:61]    [54255662:54255723]
             11         ?          ?       ?           [0:61]    [54201651:54201712]
             12         ?          ?       ?           [8:60]    [54206009:54206061]
             13         ?          ?       ?          [10:61]    [54178987:54179038]
             14         ?          ?       ?           [8:61]    [54212018:54212071]
             15         ?          ?       ?           [8:51]    [54234278:54234321]
             16         ?          ?       ?           [8:61]    [54238143:54238196]

我们得到了和前面看到的BLAST hit详细程度相似的结果。但是有些不同需要解释：

-  e-value和bit score列的值。因为BLAT HSP没有e-values和bit scores，默
   认显示‘?’.
-  span列是怎么回事呢？span值本来是显示完整的比对长度，包含所有的残基和
   gap。但是PSL格式目前还不支持这些信息并且 ``Bio.SearchIO`` 也不打算去
   猜它到底是多少，所有我们得到了和e-value以及bit score列相同的 ‘?’。 

就Python对象来说， ``Hit`` 和列表行为最相似，但是额外含有 ``HSP`` 。如果
你对列表熟悉，在使用 ``Hit``对象是不会遇到困难。

和列表一样， ``Hit`` 对象是可迭代的，并且每次迭代返回一个 ``HSP`` 对象：

.. code:: verbatim

    >>> for hsp in blast_hit:
    ...     hsp
    HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='42291', 1 fragments)
    HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='42291', 1 fragments)

你可以对 ``Hit`` 对象调用 ``len`` 方法查看它含有多少个 ``HSP`` 对象：

.. code:: verbatim

    >>> len(blast_hit)
    2
    >>> len(blat_hit)
    17

你可以对 ``Hit``对象使用切片取得单个或多个 ``HSP`` 对象，和 ``QueryResult``
一样，如果切取多个 ``HSP``  ，会返回包含被切 ``HSP``  的一个新 ``Hit``对象。

.. code:: verbatim

    >>> blat_hit[0]                 # retrieve single items
    HSP(hit_id='chr19', query_id='mystery_seq', 1 fragments)
    >>> sliced_hit = blat_hit[4:9]  # retrieve multiple items
    >>> len(sliced_hit)
    5
    >>> print sliced_hit
    Query: mystery_seq
           <unknown description>
      Hit: chr19 (59128983)
           <unknown description>
     HSPs: ----  --------  ---------  ------  ---------------  ---------------------
              #   E-value  Bit score    Span      Query range              Hit range
           ----  --------  ---------  ------  ---------------  ---------------------
              0         ?          ?       ?           [0:60]    [54198476:54198536]
              1         ?          ?       ?           [0:61]    [54265610:54265671]
              2         ?          ?       ?           [0:61]    [54238143:54240175]
              3         ?          ?       ?           [0:60]    [54189735:54189795]
              4         ?          ?       ?           [0:61]    [54185425:54185486]

你同样可以对一个 ``Hit`` 里的 ``HSP``  排序，和你在 ``QueryResult`` 对象
中看到的方法一样。

最后，同样可以对 ``Hit`` 对象使用 ``filter`` 和 ``map``方法。和 ``QueryResult`` 
不同， ``Hit`` 对象只有一种 ``filter`` (``Hit.filter``) 和一种 ``map`` (``Hit.map``)。

8.1.3  HSP
~~~~~~~~~~

``HSP`` (高分片段)代表hit序列中的一个区域，该区域包含对于查询序列有意义的
比对。它包含了你的查询序列和一个数据库条目之间精确的匹配。由于匹配取决于
序列搜索工具的算法， ``HSP``  含有大部分统计信息，这些统计是由搜索工具计
算得到的。这使得不同搜索工具的 ``HSP``  对象之间的差异和你在 ``QueryResult`` 
以及 ``Hit`` 对象看到的差异更加明显，

我们来看看BLAST和BLAT搜索的例子。先看BLAST HSP：

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> blast_qresult = SearchIO.read('my_blast.xml', 'blast-xml')
    >>> blast_hsp = blast_qresult[0][0]    # first hit, first hsp

.. code:: verbatim

    >>> print blast_hsp
          Query: 42291 mystery_seq
            Hit: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520...
    Query range: [0:61] (1)
      Hit range: [0:61] (1)
    Quick stats: evalue 4.9e-23; bitscore 111.29
      Fragments: 1 (61 columns)
         Query - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
           Hit - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

和 ``QueryResult`` 以及 ``Hit``类似，调用 ``HSP``  的 ``print`` 方法,
显示细节：

-  有query和hit ID以及描述。我们需要这些来辨识我买的 ``HSP``  。
-  我们同样得到了query和hit序列的匹配范围。这里用的的切片标志着范围的表示
   是使用Python的索引风格（从0开始，半开区间）。圆括号里的数字表示正负链。
   这里，两条序列都是正链。
-  还有一些简短统计：e-value和bitscore。
-  还有一些HSP片段的信息。现在可以忽略，稍后会解释。
-  最后，还有query和hit的比对本身。

这些信息可以用点标记从它们本身获得，和 ``Hit`` 以及 ``QueryResult``相同： 

.. code:: verbatim

    >>> blast_hsp.query_range
    (0, 61)

.. code:: verbatim

    >>> blast_hsp.evalue
    4.91307e-23

They’re not the only attributes available, though. ``HSP`` objects come
with a default set of properties that makes it easy to probe their
various details. Here are some examples:

.. code:: verbatim

    >>> blast_hsp.hit_start         # start coordinate of the hit sequence
    0
    >>> blast_hsp.query_span        # how many residues in the query sequence
    61
    >>> blast_hsp.aln_span          # how long the alignment is
    61

Check out the ``HSP``
`documentation <http://biopython.org/DIST/docs/api/Bio.SearchIO._model.hsp-module.html>`__
for a full list of these predefined properties.

Furthermore, each sequence search tool usually computes its own
statistics / details for its ``HSP`` objects. For example, an XML BLAST
search also outputs the number of gaps and identical residues. These
attributes can be accessed like so:

.. code:: verbatim

    >>> blast_hsp.gap_num       # number of gaps
    0
    >>> blast_hsp.ident_num     # number of identical residues
    61

These details are format-specific; they may not be present in other
formats. To see which details are available for a given sequence search
tool, you should check the format’s documentation in ``Bio.SearchIO``.
Alternatively, you may also use ``.__dict__.keys()`` for a quick list of
what’s available:

.. code:: verbatim

    >>> blast_hsp.__dict__.keys()
    ['bitscore', 'evalue', 'ident_num', 'gap_num', 'bitscore_raw', 'pos_num', '_items']

Finally, you may have noticed that the ``query`` and ``hit`` attributes
of our HSP are not just regular strings:

.. code:: verbatim

    >>> blast_hsp.query
    SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet()), id='42291', name='aligned query sequence', description='mystery_seq', dbxrefs=[])
    >>> blast_hsp.hit
    SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet()), id='gi|262205317|ref|NR_030195.1|', name='aligned hit sequence', description='Homo sapiens microRNA 520b (MIR520B), microRNA', dbxrefs=[])

They are ``SeqRecord`` objects you saw earlier in
Section \ `4 <#chapter:SeqRecord>`__! This means that you can do all
sorts of interesting things you can do with ``SeqRecord`` objects on
``HSP.query`` and/or ``HSP.hit``.

It should not surprise you now that the ``HSP`` object has an
``alignment`` property which is a ``MultipleSeqAlignment`` object:

.. code:: verbatim

    >>> print blast_hsp.aln
    DNAAlphabet() alignment with 2 rows and 61 columns
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG 42291
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG gi|262205317|ref|NR_030195.1|

Having probed the BLAST HSP, let’s now take a look at HSPs from our BLAT
results for a different kind of HSP. As usual, we’ll begin by invoking
``print`` on it:

.. code:: verbatim

    >>> blat_qresult = SearchIO.read('my_blat.psl', 'blat-psl')
    >>> blat_hsp = blat_qresult[0][0]       # first hit, first hsp
    >>> print blat_hsp
          Query: mystery_seq <unknown description>
            Hit: chr19 <unknown description>
    Query range: [0:61] (1)
      Hit range: [54204480:54204541] (1)
    Quick stats: evalue ?; bitscore ?
      Fragments: 1 (? columns)

Some of the outputs you may have already guessed. We have the query and
hit IDs and descriptions and the sequence coordinates. Values for evalue
and bitscore is ‘?’ as BLAT HSPs do not have these attributes. But The
biggest difference here is that you don’t see any sequence alignments
displayed. If you look closer, PSL formats themselves do not have any
hit or query sequences, so ``Bio.SearchIO`` won’t create any sequence or
alignment objects. What happens if you try to access ``HSP.query``,
``HSP.hit``, or ``HSP.aln``? You’ll get the default values for these
attributes, which is ``None``:

.. code:: verbatim

    >>> blat_hsp.hit is None
    True
    >>> blat_hsp.query is None
    True
    >>> blat_hsp.aln is None
    True

This does not affect other attributes, though. For example, you can
still access the length of the query or hit alignment. Despite not
displaying any attributes, the PSL format still have this information so
``Bio.SearchIO`` can extract them:

.. code:: verbatim

    >>> blat_hsp.query_span     # length of query match
    61
    >>> blat_hsp.hit_span       # length of hit match
    61

Other format-specific attributes are still present as well:

.. code:: verbatim

    >>> blat_hsp.score          # PSL score
    61
    >>> blat_hsp.mismatch_num   # the mismatch column
    0

So far so good? Things get more interesting when you look at another
‘variant’ of HSP present in our BLAT results. You might recall that in
BLAT searches, sometimes we get our results separated into ‘blocks’.
These blocks are essentially alignment fragments that may have some
intervening sequence between them.

Let’s take a look at a BLAT HSP that contains multiple blocks to see how
``Bio.SearchIO`` deals with this:

.. code:: verbatim

    >>> blat_hsp2 = blat_qresult[0][1]      # first hit, second hsp
    >>> print blat_hsp2
          Query: mystery_seq <unknown description>
            Hit: chr19 <unknown description>
    Query range: [0:61] (1)
      Hit range: [54233104:54264463] (1)
    Quick stats: evalue ?; bitscore ?
      Fragments: ---  --------------  ----------------------  ----------------------
                   #            Span             Query range               Hit range
                 ---  --------------  ----------------------  ----------------------
                   0               ?                  [0:18]     [54233104:54233122]
                   1               ?                 [18:61]     [54264420:54264463]

What’s happening here? We still some essential details covered: the IDs
and descriptions, the coordinates, and the quick statistics are similar
to what you’ve seen before. But the fragments detail is all different.
Instead of showing ‘Fragments: 1’, we now have a table with two data
rows.

This is how ``Bio.SearchIO`` deals with HSPs having multiple fragments.
As mentioned before, an HSP alignment may be separated by intervening
sequences into fragments. The intervening sequences are not part of the
query-hit match, so they should not be considered part of query nor hit
sequence. However, they do affect how we deal with sequence coordinates,
so we can’t ignore them.

Take a look at the hit coordinate of the HSP above. In the
``Hit range:`` field, we see that the coordinate is
``[54233104:54264463]``. But looking at the table rows, we see that not
the entire region spanned by this coordinate matches our query.
Specifically, the intervening region spans from ``54233122`` to
``54264420``.

Why then, is the query coordinates seem to be contiguous, you ask? This
is perfectly fine. In this case it means that the query match is
contiguous (no intervening regions), while the hit match is not.

All these attributes are accessible from the HSP directly, by the way:

.. code:: verbatim

    >>> blat_hsp2.hit_range         # hit start and end coordinates of the entire HSP
    (54233104, 54264463)
    >>> blat_hsp2.hit_range_all     # hit start and end coordinates of each fragment
    [(54233104, 54233122), (54264420, 54264463)]
    >>> blat_hsp2.hit_span          # hit span of the entire HSP
    31359
    >>> blat_hsp2.hit_span_all      # hit span of each fragment
    [18, 43]
    >>> blat_hsp2.hit_inter_ranges  # start and end coordinates of intervening regions in the hit sequence
    [(54233122, 54264420)]
    >>> blat_hsp2.hit_inter_spans   # span of intervening regions in the hit sequence
    [31298]

Most of these attributes are not readily available from the PSL file we
have, but ``Bio.SearchIO`` calculates them for you on the fly when you
parse the PSL file. All it needs are the start and end coordinates of
each fragment.

What about the ``query``, ``hit``, and ``aln`` attributes? If the HSP
has multiple fragments, you won’t be able to use these attributes as
they only fetch single ``SeqRecord`` or ``MultipleSeqAlignment``
objects. However, you can use their ``*_all`` counterparts:
``query_all``, ``hit_all``, and ``aln_all``. These properties will
return a list containing ``SeqRecord`` or ``MultipleSeqAlignment``
objects from each of the HSP fragment. There are other attributes that
behave similarly, i.e. they only work for HSPs with one fragment. Check
out the ``HSP``
`documentation <http://biopython.org/DIST/docs/api/Bio.SearchIO._model.hsp-module.html>`__
for a full list.

Finally, to check whether you have multiple fragments or not, you can
use the ``is_fragmented`` property like so:

.. code:: verbatim

    >>> blat_hsp2.is_fragmented     # BLAT HSP with 2 fragments
    True
    >>> blat_hsp.is_fragmented      # BLAT HSP from earlier, with one fragment
    False

Before we move on, you should also know that we can use the slice
notation on ``HSP`` objects, just like ``QueryResult`` or ``Hit``
objects. When you use this notation, you’ll get an ``HSPFragment``
object in return, the last component of the object model.

8.1.4  HSPFragment
~~~~~~~~~~~~~~~~~~

``HSPFragment`` represents a single, contiguous match between the query
and hit sequences. You could consider it the core of the object model
and search result, since it is the presence of these fragments that
determine whether your search have results or not.

In most cases, you don’t have to deal with ``HSPFragment`` objects
directly since not that many sequence search tools fragment their HSPs.
When you do have to deal with them, what you should remember is that
``HSPFragment`` objects were written with to be as compact as possible.
In most cases, they only contain attributes directly related to
sequences: strands, reading frames, alphabets, coordinates, the
sequences themselves, and their IDs and descriptions.

These attributes are readily shown when you invoke ``print`` on an
``HSPFragment``. Here’s an example, taken from our BLAST search:

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> blast_qresult = SearchIO.read('my_blast.xml', 'blast-xml')
    >>> blast_frag = blast_qresult[0][0][0]    # first hit, first hsp, first fragment
    >>> print blast_frag
          Query: 42291 mystery_seq
            Hit: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520...
    Query range: [0:61] (1)
      Hit range: [0:61] (1)
      Fragments: 1 (61 columns)
         Query - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
           Hit - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

At this level, the BLAT fragment looks quite similar to the BLAST
fragment, save for the query and hit sequences which are not present:

.. code:: verbatim

    >>> blat_qresult = SearchIO.read('my_blat.psl', 'blat-psl')
    >>> blat_frag = blat_qresult[0][0][0]    # first hit, first hsp, first fragment
    >>> print blat_frag
          Query: mystery_seq <unknown description>
            Hit: chr19 <unknown description>
    Query range: [0:61] (1)
      Hit range: [54204480:54204541] (1)
      Fragments: 1 (? columns)

In all cases, these attributes are accessible using our favorite dot
notation. Some examples:

.. code:: verbatim

    >>> blast_frag.query_start      # query start coordinate
    0
    >>> blast_frag.hit_strand       # hit sequence strand
    1
    >>> blast_frag.hit              # hit sequence, as a SeqRecord object
    SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet()), id='gi|262205317|ref|NR_030195.1|', name='aligned hit sequence', description='Homo sapiens microRNA 520b (MIR520B), microRNA', dbxrefs=[])

8.2  A note about standards and conventions
-------------------------------------------

Before we move on to the main functions, there is something you ought to
know about the standards ``Bio.SearchIO`` uses. If you’ve worked with
multiple sequence search tools, you might have had to deal with the many
different ways each program deals with things like sequence coordinates.
It might not have been a pleasant experience as these search tools
usually have their own standards. For example, one tools might use
one-based coordinates, while the other uses zero-based coordinates. Or,
one program might reverse the start and end coordinates if the strand is
minus, while others don’t. In short, these often creates unnecessary
mess must be dealt with.

We realize this problem ourselves and we intend to address it in
``Bio.SearchIO``. After all, one of the goals of ``Bio.SearchIO`` is to
create a common, easy to use interface to deal with various search
output files. This means creating standards that extend beyond the
object model you just saw.

Now, you might complain, "Not another standard!". Well, eventually we
have to choose one convention or the other, so this is necessary. Plus,
we’re not creating something entirely new here; just adopting a standard
we think is best for a Python programmer (it is Biopython, after all).

There are three implicit standards that you can expect when working with
``Bio.SearchIO``:

-  The first one pertains to sequence coordinates. In ``Bio.SearchIO``,
   all sequence coordinates follows Python’s coordinate style:
   zero-based and half open. For example, if in a BLAST XML output file
   the start and end coordinates of an HSP are 10 and 28, they would
   become 9 and 28 in ``Bio.SearchIO``. The start coordinate becomes 9
   because Python indices start from zero, while the end coordinate
   remains 28 as Python slices omit the last item in an interval.
-  The second is on sequence coordinate orders. In ``Bio.SearchIO``,
   start coordinates are always less than or equal to end coordinates.
   This isn’t always the case with all sequence search tools, as some of
   them have larger start coordinates when the sequence strand is minus.
-  The last one is on strand and reading frame values. For strands,
   there are only four valid choices: ``1`` (plus strand), ``-1`` (minus
   strand), ``0`` (protein sequences), and ``None`` (no strand). For
   reading frames, the valid choices are integers from ``-3`` to ``3``
   and ``None``.

Note that these standards only exist in ``Bio.SearchIO`` objects. If you
write ``Bio.SearchIO`` objects into an output format, ``Bio.SearchIO``
will use the format’s standard for the output. It does not force its
standard over to your output file.

8.3  Reading search output files
--------------------------------

There are two functions you can use for reading search output files into
``Bio.SearchIO`` objects: ``read`` and ``parse``. They’re essentially
similar to ``read`` and ``parse`` functions in other submodules like
``Bio.SeqIO`` or ``Bio.AlignIO``. In both cases, you need to supply the
search output file name and the file format name, both as Python
strings. You can check the documentation for a list of format names
``Bio.SearchIO`` recognizes.

``Bio.SearchIO.read`` is used for reading search output files with only
one query and returns a ``QueryResult`` object. You’ve seen ``read``
used in our previous examples. What you haven’t seen is that ``read``
may also accept additional keyword arguments, depending on the file
format.

Here are some examples. In the first one, we use ``read`` just like
previously to read a BLAST tabular output file. In the second one, we
use a keyword argument to modify so it parses the BLAST tabular variant
with comments in it:

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('tab_2226_tblastn_003.txt', 'blast-tab')
    >>> qresult
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)
    >>> qresult2 = SearchIO.read('tab_2226_tblastn_007.txt', 'blast-tab', comments=True)
    >>> qresult2
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

These keyword arguments differs among file formats. Check the format
documentation to see if it has keyword arguments that modifies its
parser’s behavior.

As for the ``Bio.SearchIO.parse``, it is used for reading search output
files with any number of queries. The function returns a generator
object that yields a ``QueryResult`` object in each iteration. Like
``Bio.SearchIO.read``, it also accepts format-specific keyword
arguments:

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('tab_2226_tblastn_001.txt', 'blast-tab')
    >>> for qresult in qresults:
    ...     print qresult.id
    gi|16080617|ref|NP_391444.1|
    gi|11464971:4-101
    >>> qresults2 = SearchIO.parse('tab_2226_tblastn_005.txt', 'blast-tab', comments=True)
    >>> for qresult in qresults2:
    ...     print qresult.id
    random_s00
    gi|16080617|ref|NP_391444.1|
    gi|11464971:4-101

8.4  Dealing with large search output files with indexing
---------------------------------------------------------

Sometimes, you’re handed a search output file containing hundreds or
thousands of queries that you need to parse. You can of course use
``Bio.SearchIO.parse`` for this file, but that would be grossly
inefficient if you need to access only a few of the queries. This is
because ``parse`` will parse all queries it sees before it fetches your
query of interest.

In this case, the ideal choice would be to index the file using
``Bio.SearchIO.index`` or ``Bio.SearchIO.index_db``. If the names sound
familiar, it’s because you’ve seen them before in
Section \ `5.4.2 <#sec:SeqIO-index>`__. These functions also behave
similarly to their ``Bio.SeqIO`` counterparts, with the addition of
format-specific keyword arguments.

Here are some examples. You can use ``index`` with just the filename and
format name:

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> idx = SearchIO.index('tab_2226_tblastn_001.txt', 'blast-tab')
    >>> sorted(idx.keys())
    ['gi|11464971:4-101', 'gi|16080617|ref|NP_391444.1|']
    >>> idx['gi|16080617|ref|NP_391444.1|']
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

Or also with the format-specific keyword argument:

.. code:: verbatim

    >>> idx = SearchIO.index('tab_2226_tblastn_005.txt', 'blast-tab', comments=True)
    >>> sorted(idx.keys())
    ['gi|11464971:4-101', 'gi|16080617|ref|NP_391444.1|', 'random_s00']
    >>> idx['gi|16080617|ref|NP_391444.1|']
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

Or with the ``key_function`` argument, as in ``Bio.SeqIO``:

.. code:: verbatim

    >>> key_function = lambda id: id.upper()    # capitalizes the keys
    >>> idx = SearchIO.index('tab_2226_tblastn_001.txt', 'blast-tab', key_function=key_function)
    >>> sorted(idx.keys())
    ['GI|11464971:4-101', 'GI|16080617|REF|NP_391444.1|']
    >>> idx['GI|16080617|REF|NP_391444.1|']
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

``Bio.SearchIO.index_db`` works like as ``index``, only it writes the
query offsets into an SQLite database file.

8.5  Writing and converting search output files
-----------------------------------------------

It is occasionally useful to be able to manipulate search results from
an output file and write it again to a new file. ``Bio.SearchIO``
provides a ``write`` function that lets you do exactly this. It takes as
its arguments an iterable returning ``QueryResult`` objects, the output
filename to write to, the format name to write to, and optionally some
format-specific keyword arguments. It returns a four-item tuple, which
denotes the number or ``QueryResult``, ``Hit``, ``HSP``, and
``HSPFragment`` objects that were written.

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('mirna.xml', 'blast-xml')     # read XML file
    >>> SearchIO.write(qresults, 'results.tab', 'blast-tab')    # write to tabular file
    (3, 239, 277, 277)

You should note different file formats require different attributes of
the ``QueryResult``, ``Hit``, ``HSP`` and ``HSPFragment`` objects. If
these attributes are not present, writing won’t work. In other words,
you can’t always write to the output format that you want. For example,
if you read a BLAST XML file, you wouldn’t be able to write the results
to a PSL file as PSL files require attributes not calculated by BLAST
(e.g. the number of repeat matches). You can always set these attributes
manually, if you really want to write to PSL, though.

Like ``read``, ``parse``, ``index``, and ``index_db``, ``write`` also
accepts format-specific keyword arguments. Check out the documentation
for a complete list of formats ``Bio.SearchIO`` can write to and their
arguments.

Finally, ``Bio.SearchIO`` also provides a ``convert`` function, which is
simply a shortcut for ``Bio.SearchIO.parse`` and ``Bio.SearchIO.write``.
Using the convert function, our example above would be:

.. code:: verbatim

    >>> from Bio import SearchIO
    >>> SearchIO.convert('mirna.xml', 'blast-xml', 'results.tab', 'blast-tab')
    (3, 239, 277, 277)

As ``convert`` uses ``write``, it is only limited to format conversions
that have all the required attributes. Here, the BLAST XML file provides
all the default values a BLAST tabular file requires, so it works just
fine. However, other format conversions are less likely to work since
you need to manually assign the required attributes first.
