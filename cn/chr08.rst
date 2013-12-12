第8章  BLAST和其他序列搜索工具(*实验性质的代码*)
======================================================================

*WARNING*: 这章教程介绍了Biopython中一个 *实验的* 模块。它正在被加入到
Biopython中，并且以一个预尾期的状态整理到教程当中，这样在我们发布稳定版
的之前可以收到一系列的反馈和并作改进。那时有些细节可能会改变，并且用到
当前 ``Bio.SearchIO`` 模块的脚本也需要更新。切记！为了与NCBI BLAST有关的
代码可以稳定工作，请继续使用第 \ `7 <#chapter:blast>`__ 章介绍的Bio.Blast。

生物序列的鉴定是生物信息工作的主要部分。有几个工具像BLAST（可能是最流行
的），FASTA ，HMMER还有许多其它的都有这个功能，每个工具都有独特的算法和
途径。一般来说，这些工具都是用你的序列在数据库中搜索可能的匹配。随着序列
数量的增加（匹配数也会随之增加），将会有成百上千的可能匹配，解析这些结果
无疑变得越来越困难。理所当然，人工解析搜索结果变得不可能。而且你可能会同
时用几种不同的搜索工具，每种工具都有独特的统计方法、规则和输出格式。可以
想象，同时用多种工具搜索多条序列是多么恐怖的事。

我们对此非常了解，所以我们在Biopython构建了 ``Bio.SearchIO`` 亚模块。``Bio.SearchIO`` 模块使从搜索结果中提取信息变得简单，并且可以处理不同工具
的不同标准和规则。``SearchIO`` 和BioPerl中模块名称一致。

在本章中，我们将学习 ``Bio.SearchIO`` 的主要功能，了解它可以做什么。我
们将使用两个主要的搜索工具：BLAST和FASTA。它们只是用来阐明思路，你可以轻
易地把工作流程应用到 ``Bio.SearchIO`` 支持的其他工具中。欢迎你使用我们将要
用到的搜索结果文件。BLAST搜索结果文件可以在
`这里 <http://biopython.org/SRC/Tests/Tutorial/my_blast.xml>`__ 下载。
BLAT输出结果文件可以在
`这里 <http://biopython.org/SRC/Tests/Tutorial/my_blat.psl>`__ 下载。两个结果
文件都是用下面这条序列搜索产生的：

.. code:: python

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

尽管多数搜索工具的输出风格极为不同，但是它们蕴含的理念很相似：


-  输出文件可能包含一条或更多的搜索查询的结果。
-  在每次查询中，你会在给定的数据库中得到一个或更多的hit（命中）。
-  在每个hit中，你会得到一个或更多包含query（查询)序列和数据库序列实际比对的区域。
-  一些工具如BLAT和Exonerate可能会把这些区域分成几个比对片段（或在BLAT中
   称为区块，在Exonerate中称为可能外显子）。这并不是很常见，像BLAST和
   HMMER就不这么做。

知道这些共性之后，我们决定把它们作为创造 ``Bio.SearchIO`` 对象模型的基础。对
象模型是Python对象组成的嵌套分级系统，每个对象都代表一个上面列出的概念。这些对
象是：

-  ``QueryResult``，代表单个搜索查询。
-  ``Hit``，代表单个的数据库hit。``Hit`` 对象包含在 ``QueryResult`` 中，
   每个 ``QueryResult`` 中有0个或多个 ``Hit`` 对象。
-  ``HSP`` (high-scoring pair（高分片段）)，代表query和hit序列中有
   意义比对的区域。``HSP`` 对象包含在 ``Hit`` 对象中，而且每个 ``Hit`` 有一个
   或更多的 ``HSP`` 对象。   
-  ``HSPFragment``，代表query和hit序列中单个的邻近比对。 ``HSPFragment``
   对象包含在 ``HSP`` 对象中。多数的搜索工具如BLAST和HMMER把 ``HSP`` 和
   ``HSPFragment`` 合并，因为一个 ``HSP`` 只含有一个 ``HSPFragment``。但是
   像BLAT和Exonerate会产生含有多个 ``HSPFragment`` 的 ``HSP`` 。似乎有些困
   惑？不要紧，稍后我们将详细介绍这两个对象。

这四个对象是当你用 ``Bio.SearchIO`` 会碰到的。 ``Bio.SearchIO`` 四
个主要方法： ``read`` ，``parse``，``index`` 或 ``index_db`` 中任意一个都可
以产生这四个对象。这些方法的会在后面的部分详细介绍。这部分只会用到 ``read`` 和
``parse`` ，这两个方法和 ``Bio.SeqIO`` 以及 ``Bio.AlignIO`` 中的 ``read`` 和 ``parse`` 方法功
能相似：

-  ``read`` 用于单query对输出文件进行搜索并且返回一个 ``QueryResult`` 对象。
-  ``parse`` 用于多query对输出文件进行搜索并且返回一个可以yield ``QueryResult`` 对象的generator。

完成这些之后，让我们开始学习每个 ``Bio.SearchIO`` 对象，从 ``QueryResult`` 开始。

8.1.1  QueryResult
~~~~~~~~~~~~~~~~~~

``QueryResult``，代表单query搜索，每个 ``QueryResult`` 中有0个或多个 ``Hit`` 对象。我们来看看BLAST文件是什么样的：

.. code:: python

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

虽然我们才接触对象模型的皮毛，但是你已经可以看到一些有用的信息了。通过调用``QueryResult`` 对象的 ``print`` 方法，你可以看到：

-  程序的名称和版本 (blastn version 2.2.27+)
-  query的ID，描述和序列长度(ID是42291，描述是 ‘mystery\_seq’，长度是61)
-  搜索的目标数据库 (refseq\_rna)
-  hit结果的快速预览。对于我们的查询序列，有100个可能的hit（表格中表示的是
   0-99）对于每个hit，我们可以看到它包含的高分比对片段（HSP)，ID和一个片
   段描述。注意， ``Bio.SearchIO`` 截断了表格，只显示0-29，然后是97-99。
 
现在让我们用同样的步骤来检查BLAT的结果：

.. code:: python

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

马上可以看到一些不同点。有些是由于BLAT使用PSL格式储存它的信息，稍后会看
到。其余是由于BLAST和BLAT搜索的程序和数据库之间明显的差异造成的：

-  程序名称和版本。 ``Bio.SearchIO`` 知道程序是BLAST，但是在输出文件中没
   有信息显示程序版本，所以默认是 ‘<unknown version>’。
-  query的ID，描述和序列的长度。注意，这些细节和BLAST的细节只有细小的差别，
   ID是 ‘mystery\_seq’ 而不是42991，缺少描述，但是序列长度仍是61。这
   实际上是文件格式本身导致的差异。BLAST有时创建自己的query ID并且用你的原
   始ID作为序列描述。
-  目标数据库是未知的，因为BLAT输出文件没提到相关信息。
-  最后，hit列表完全不同。这里，我们的查询序列只命中到 ‘chr19’ 数据库条
   目，但是我们可以看到它含有17个HSP区域。这并不让人诧异，考虑到我们
   使用的是不同的程序，并且这些程序都有自己的数据库。

所有通过调用 ``print`` 方法看到的信息都可以单独地用Python的对象属性获得（又叫点标记法）。同样还可以用相同的方法获得其他格式特有的属性。

.. code:: python

    >>> print "%s %s" % (blast_qresult.program, blast_qresult.version)
    blastn 2.2.27+
    >>> print "%s %s" % (blat_qresult.program, blat_qresult.version)
    blat <unknown version>
    >>> blast_qresult.param_evalue_threshold    # blast-xml specific
    10.0

想获得一个可访问属性的完整列表，可以查询每个格式特有的文档。这些是 
`BLAST <http://biopython.org/DIST/docs/api/Bio.SearchIO.BlastIO-module.html>`__
`BLAT <http://biopython.org/DIST/docs/api/Bio.SearchIO.BlatIO-module.html>`__.

已经知道了在 ``QueryResult`` 对象上可以调用 ``print`` 方法，让我们研究的更深
一些。 ``QueryResult`` 到底是什么？就Python对象来说， ``QueryResult`` 混合
了列表和字典的特性。换句话说，也就是一个包含了列表和字典功能的容器对象。

和列表以及字典一样， ``QueryResult`` 对象是可迭代的。每次迭代返回一个hit
对象：

.. code:: python

    >>> for hit in blast_qresult:
    ...     hit
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)
    Hit(id='gi|301171311|ref|NR_035856.1|', query_id='42291', 1 hsps)
    Hit(id='gi|270133242|ref|NR_032573.1|', query_id='42291', 1 hsps)
    Hit(id='gi|301171322|ref|NR_035857.1|', query_id='42291', 2 hsps)
    Hit(id='gi|301171267|ref|NR_035851.1|', query_id='42291', 1 hsps)
    ...

要得到 ``QueryResult`` 对象有多少hit，可以简单调用Python的 ``len`` 方法：

.. code:: python

    >>> len(blast_qresult)
    100
    >>> len(blat_qresult)
    1

同列表类似，你可以用切片来获得 ``QueryResult`` 对象的hit：

.. code:: python

    >>> blast_qresult[0]        # retrieves the top hit
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)
    >>> blast_qresult[-1]       # retrieves the last hit
    Hit(id='gi|397513516|ref|XM_003827011.1|', query_id='42291', 1 hsps)

要得到多个hit，你同样可以对 ``QueryResult`` 对象作切片。这种情况下，返回一个包含被切hit的新 ``QueryResult`` 对象：

.. code:: python

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

同字典类似，可以通过ID获取hit。如果你知道一个特定的hit ID存在于一个搜索结果中时，特别有用：

.. code:: python

    >>> blast_qresult['gi|262205317|ref|NR_030195.1|']
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)

你可以用 ``hits`` 方法获得完整的 ``Hit`` 对象，也可以用 ``hit_keys`` 方法获得完整的 ``Hit`` ID：

.. code:: python

    >>> blast_qresult.hits
    [...]       # list of all hits
    >>> blast_qresult.hit_keys
    [...]       # list of all hit IDs

如果你想确定一个特定的hit是否存在于查询结果中该怎么做呢？可以用 ``in`` 关键字作一个简单的成员检验：

.. code:: python

    >>> 'gi|262205317|ref|NR_030195.1|' in blast_qresult
    True
    >>> 'gi|262205317|ref|NR_030194.1|' in blast_qresult
    False

有时候，只知道一个hit是否存在是不够的，你可能也会想知道hit的排名。 ``index`` 方法可以帮助你：

.. code:: python

    >>> blast_qresult.index('gi|301171437|ref|NR_035870.1|')
    22

记住，我们用的是Python风格的索引，是从0开始。这代表hit的排名是23而不是22。

同样，注意你看的hit排名是基于原始搜索输出文件的本来顺序。不同的搜索工具可
能会基于不同的标准排列hit。

如果原本的hit排序不合你意，可以用 ``QueryResult`` 对象的 ``sort`` 方法。
它和Python的 ``list.sort`` 方法很相似，只是有个是否创建一个新的排序后的
``QueryResult`` 对象的选项。

这里有个用 ``QueryResult.sort`` 方法对hit排序的例子，这个方法基于每个hit
的完整序列长度。对于这个特殊的排序，我们设置 ``in_place`` 参数等于 ``False`` ，
这样排序方法会返回一个新的 ``QueryResult`` 对象，而原来的对象是未排序的。
我们同样可以设置 ``reverse`` 参数等于 `` True `` 以递减排序。

.. code:: python

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

有 ``in_place`` 参数的好处是可以保留原本的顺序，后面会用到。注意这不是 ``QueryResult.sort`` 的默认行为，需要我们明确地设置 ``in_place`` 为 ``True`` 。

现在，你已经知道使用 ``QueryResult`` 对象。但是，在我们学习 ``Bio.SearchIO`` 
模块下个对象前，先了解下可以让 ``QueryResult`` 对象更易使用的两个方法：
``filter`` 和 ``map`` 方法。

如果你对Python的列表推导式、generator表达式或内建的 ``filter`` 和 ``map`` 
很熟悉，就知道（不知道就是看看吧!)它们在处理list-like的对象时有多有用。
你可以用这些内建的方法来操作 ``QueryResult`` 对象，但是这只对正常list有效，并且可操作性也会受到限制。

这就是为什么 ``QueryResult`` 对象提供自己特有的 ``filter`` 和 ``map`` 
方法。对于 ``filter`` 有相似的 ``hit_filter`` 和 ``hsp_filter`` 方法，
从名称就可以看出，这些方法过滤 ``QueryResult`` 对象的 ``Hit`` 对象或者
``HSP`` 对象。同样的，对于 ``map`` ， ``QueryResult`` 对象同样提供相似
的  ``hit_map`` 和 ``hsp_map`` 方法。这些方法分别应用于 ``QueryResult`` 对象 ``hit`` 或者 ``HSP`` 对象。 

让我们来看看这些方法的功能，从 ``hit_filter`` 开始。这个方法接受一个回调
函数，这个函数检验给定的 ``Hit`` 是否符合你设定的条件。换句话说，这个方法
必须接受一个单独 ``Hit`` 对象作为参数并且返回  ``True`` 或 ``False`` 。 

这里有个用 ``hit_filter`` 筛选出只有一个HSP的 ``Hit`` 对象的例子：

.. code:: python

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

``hsp_filter`` 和 ``hit_filter`` 功能相同，只是它过滤每个hit中的 ``HSP`` 对象，
而不是 ``Hit`` 。

对于 ``map`` 方法，同样接受一个回调函数作为参数。但是回调函数返回修改过的 ``Hit`` 或 ``HSP`` 对象（取决于你是否使用 ``hit_map`` 或 ``hsp_map`` 方法），
而不是返回 ``True`` 或 ``False``。

来看一个用 ``hit_map`` 方法来重命名hit ID的例子：

.. code:: python

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

同样的， ``hsp_map`` 和 ``hit_map`` 作用相似, 但是作用于 ``HSP`` 对象而不是 ``Hit`` 对象。

8.1.2  Hit
~~~~~~~~~~

``Hit`` 对象代表从单个数据库获得所有查询结果。在 ``Bio.SearchIO`` 对象等级中是二级容器。它们被包含在 ``QueryResult`` 对象中，同时它们又包含 ``HSP`` 对象。

看看它们是什么样的，从我们的BLAST搜索开始：

.. code:: python

    >>> from Bio import SearchIO
    >>> blast_qresult = SearchIO.read('my_blast.xml', 'blast-xml')
    >>> blast_hit = blast_qresult[3]    # fourth hit from the query result

.. code:: python

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

-  query的ID和描述信息。一个hit总是和一个query绑定，所以我们同样希望记录原始
   query。这些值可以通过 ``query_id`` 和  ``query_description`` 属性从hit
   中获取。
-  我们同样得到了hit的ID、描述和序列全长。它们可以分别通过 ``id``，
   ``description``，和 ``seq_len`` 获取。
-  最后，有一个hit含有的HSP的简短信息表。在每行中，HSP重要信息被
   列出来：HSP索引，e值，得分，长度（包括gap），query坐标和hit坐标。

现在，和BLAT结果作对比。记住，在BLAT搜索结果中，我们发现有一个含有17HSP的hit。

.. code:: python

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

我们得到了和前面看到的BLAST hit详细程度相似的结果。但是有些不同点需要解释：

-  e-value和bit score列的值。因为BLAT HSP没有e-values和bit scores，默
   认显示‘?’.
-  span列是怎么回事呢？span值本来是显示完整的比对长度，包含所有的残基和
   gap。但是PSL格式目前还不支持这些信息并且 ``Bio.SearchIO`` 也不打算去猜它到底是多少，所有我们得到了和e-value以及bit score列相同的 ‘?’。 

就Python对象来说， ``Hit`` 和列表行为最相似，但是额外含有 ``HSP`` 。如果
你对列表熟悉，在使用 ``Hit`` 对象是不会遇到困难。

和列表一样， ``Hit`` 对象是可迭代的，并且每次迭代返回一个 ``HSP`` 对象：

.. code:: python

    >>> for hsp in blast_hit:
    ...     hsp
    HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='42291', 1 fragments)
    HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='42291', 1 fragments)

你可以对 ``Hit`` 对象调用 ``len`` 方法查看它含有多少个 ``HSP`` 对象：

.. code:: python

    >>> len(blast_hit)
    2
    >>> len(blat_hit)
    17

你可以对 ``Hit`` 对象作切片取得单个或多个 ``HSP`` 对象，和 ``QueryResult``
一样，如果切取多个 ``HSP``  ，会返回包含被切片 ``HSP``  的一个新 ``Hit`` 对象。

.. code:: python

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

最后，同样可以对 ``Hit`` 对象使用 ``filter`` 和 ``map`` 方法。和 ``QueryResult`` 
不同， ``Hit`` 对象只有一种 ``filter`` (``Hit.filter``) 和一种 ``map`` (``Hit.map``)。

8.1.3  HSP
~~~~~~~~~~

``HSP`` (高分片段)代表hit序列中的一个区域，该区域包含对于查询序列有意义的
比对。它包含了你的查询序列和一个数据库条目之间精确的匹配。由于匹配取决于
序列搜索工具的算法， ``HSP``  含有大部分统计信息，这些统计是由搜索工具计
算得到的。这使得不同搜索工具的 ``HSP``  对象之间的差异和你在 ``QueryResult`` 
以及 ``Hit`` 对象看到的差异更加明显。

我们来看看BLAST和BLAT搜索的例子。先看BLAST HSP：

.. code:: python

    >>> from Bio import SearchIO
    >>> blast_qresult = SearchIO.read('my_blast.xml', 'blast-xml')
    >>> blast_hsp = blast_qresult[0][0]    # first hit, first hsp

.. code:: python

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

和 ``QueryResult`` 以及 ``Hit`` 类似，调用 ``HSP``  的 ``print`` 方法,
显示细节：

-  有query和hit的ID以及描述。我们需要这些来识别我们的 ``HSP``  。
-  我们同样得到了query和hit序列的匹配范围。这里用的的切片标志着范围的表示
   是使用Python的索引风格（从0开始，半开区间）。圆括号里的数字表示正负链。
   这里，两条序列都是正链。
-  还有一些简短统计：e-value和bitscore。
-  还有一些HSP片段的信息。现在可以忽略，稍后会解释。
-  最后，还有query和hit的比对。

这些信息可以用点标记从它们本身获得，和 ``Hit`` 以及 ``QueryResult`` 相同： 

.. code:: python

    >>> blast_hsp.query_range
    (0, 61)

.. code:: python

    >>> blast_hsp.evalue
    4.91307e-23

它们并不是仅有的属性， ``HSP``  对象有一系列的属性，使得获得它们的具体信
息更加容易。下面是一些例子：

.. code:: python

    >>> blast_hsp.hit_start         # start coordinate of the hit sequence
    0
    >>> blast_hsp.query_span        # how many residues in the query sequence
    61
    >>> blast_hsp.aln_span          # how long the alignment is
    61

查看 ``HSP``
`文档 <http://biopython.org/DIST/docs/api/Bio.SearchIO._model.hsp-module.html>`__
获取完整的属性列表。

不仅如此，每个搜索工具通常会对它的 ``HSP``  对象作统计学或其他细节计算。例如，一个
XML BLAST搜索同样输出gap以及相同的残基数量。这些属性可以像这样被获取：

.. code:: python

    >>> blast_hsp.gap_num       # number of gaps
    0
    >>> blast_hsp.ident_num     # number of identical residues
    61

这些细节是格式特异的；它们可能不会出现在其他的格式中。要知道哪些细节在给
定的序列搜索工具中是存在的，你应该查看那种格式的在 ``Bio.SearchIO`` 中的
文档。或者可以用 ``.__dict__.keys()`` 获得快速列表：

.. code:: python

    >>> blast_hsp.__dict__.keys()
    ['bitscore', 'evalue', 'ident_num', 'gap_num', 'bitscore_raw', 'pos_num', '_items']

最后，你可能已经注意到了，我们HSP的 ``query`` 和 ``hit`` 属性不只是规律字符串： 


.. code:: python

    >>> blast_hsp.query
    SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet()), id='42291', name='aligned query sequence', description='mystery_seq', dbxrefs=[])
    >>> blast_hsp.hit
    SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet()), id='gi|262205317|ref|NR_030195.1|', name='aligned hit sequence', description='Homo sapiens microRNA 520b (MIR520B), microRNA', dbxrefs=[])

它们是你已经在第 \ `4 <#chapter:SeqRecord>`__ 章看到过的 ``SeqRecord`` 对象！
意味着你可以对 ``SeqRecord`` 对象做的各种有趣的事同样适用于 ``HSP.query`` 和 ``HSP.hit`` 对象。

现在 ``HSP``  对象有个 ``alignment`` 属性（一个 ``MultipleSeqAlignment`` 
对象）应该不会让你感到惊讶：

.. code:: python

    >>> print blast_hsp.aln
    DNAAlphabet() alignment with 2 rows and 61 columns
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG 42291
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG gi|262205317|ref|NR_030195.1|


探索完BLAST HSP对象，让我们看看来自BLAT结果的不一样的HSP。我们将对它调用 ``print`` 方法： 

.. code:: python

    >>> blat_qresult = SearchIO.read('my_blat.psl', 'blat-psl')
    >>> blat_hsp = blat_qresult[0][0]       # first hit, first hsp
    >>> print blat_hsp
          Query: mystery_seq <unknown description>
            Hit: chr19 <unknown description>
    Query range: [0:61] (1)
      Hit range: [54204480:54204541] (1)
    Quick stats: evalue ?; bitscore ?
      Fragments: 1 (? columns)

一些输出你应该已经猜到了。我们得到了查询序列、hit ID、描述以及序列坐标。
evalue和bitscore的值是 ‘?’ ，因为BLAT HSP并没有这些属性。但是最大的不同
是你看不到任何的序列比对展示。如果你看的更仔细，PSL格式本身并没有任何的
hit和query序列，所以 ``Bio.SearchIO`` 不会创建任何序列或者比对对象。如果
你尝试获取 ``HSP.query`` ，``HSP.hit`` ， 或者 ``HSP.aln`` 属性会怎么样
呢？你会得到这些属性的默认值 ``None`` ：

.. code:: python

    >>> blat_hsp.hit is None
    True
    >>> blat_hsp.query is None
    True
    >>> blat_hsp.aln is None
    True

这并不影响其他的属性。例如，你仍然可以获取query和hit比对的长度。尽管不显
示任何的属性，但是PSL格式还是有这些信息的，所以 ``Bio.SearchIO`` 可以抽
提出这些信息。

.. code:: python

    >>> blat_hsp.query_span     # length of query match
    61
    >>> blat_hsp.hit_span       # length of hit match
    61

其他格式特异的属性同样被展示出来：

.. code:: python

    >>> blat_hsp.score          # PSL score
    61
    >>> blat_hsp.mismatch_num   # the mismatch column
    0

到目前为止，一切还不错？当你看到BLAT结果中不同的HSP时，事情变得更有趣了。
你可能会回想起在BLAT搜索中，有时我们把结果分成 ‘blocks’ 。这些区块是必需比对片段，可能会有些内含子在它们之间。

让我们看看 ``Bio.SearchIO`` 怎么处理包含多个区块的BLAT HSP：

.. code:: python

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

怎么回事？我们仍然得到了一些必要的信息：ID，描述信息，坐标和快速统计，和
你前面看到的一样。但是片段信息完全不同。我们得到了有两行数据的表格，而不是显示 ‘Fragment: 1’。

这就是 ``Bio.SearchIO`` 处理含有多片段HSP的方式。和前面提到的一样，一个
HSP比对可能会被内含子分成多个片段。内含子不是query-hit匹配的一部分，所以
它们不能被当成query或hit序列的一部分。但是，它们确实影响我们处理序列坐标，
所以我们不能忽视。

看看上面的HSP的hit坐标。在 ``Hit range`` 区域，我们看到坐标是
``[54233104:54264463]``。但是看看表格中的行，我们发现不是坐标跨度的所有区域
都能匹配我们的query。特殊的是，间断区域从 ``54233122`` 到 ``54264420`` 。

你可能会问，为什么query坐标好像是邻近的?这是很好的。在这个例子中，query是连续的（无间断区域），但是hit却不是。

所有的这些属性都是可以直接从HSP获取的，通过这样的方式：

.. code:: python

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

这些属性中大多数都不能简单地从PSL文件获得，但是当你分析PSL文件时，
``Bio.SearchIO`` 会动态地帮你计算。它需要的只是每个片段的开始和结束坐标。

``query``， ``hit``， 和 ``aln`` 属性又是什么情况？如果HSP含有多个片段，
你就不能使用这些属性，因为它们只取回单个 ``SeqRecord`` 或
``MultipleSeqAlignment`` 对象。但是，你可以用相应的 ``*_all`` 方法：
``query_all``， ``hit_all``， 和 ``aln_all``。 这些属性会返回包含每个HSP
片段的 ``SeqRecord`` 或 ``MultipleSeqAlignment`` 对象的列表。还有其他相同
功能的属性，也就是只对只有一个片段的HSP有效。查看 ``HSP``
`文档 <http://biopython.org/DIST/docs/api/Bio.SearchIO._model.hsp-module.html>`__
获得完整的列表。

最后，想要检查是否是多片段HSP，你可以用 ``is_fragmented`` 属性：

.. code:: python

    >>> blat_hsp2.is_fragmented     # BLAT HSP with 2 fragments
    True
    >>> blat_hsp.is_fragmented      # BLAT HSP from earlier, with one fragment
    False

在进入下部分之前，你只需要了解我们可以对 ``HSP`` 对象使用切片，和
``QueryResult`` 或 ``Hit`` 对象一样。当你使用切片的时候，会返回一个
``HSPFragment`` 对象。

8.1.4  HSP片段
~~~~~~~~~~~~~~~~~~

``HSPFragment`` 代表query和hit之间单个连续匹配。应该把它当作对象模型
和搜索结果的核心，因为它决定你的搜索是否有结果。

在多数情况下，你不必直接处理 ``HSPFragment`` 对象，因为没有那么多搜索工具
断裂它们的HSP。当你确实需要处理它们时，需要记住的是 ``HSPFragment`` 对象
要被写地尽量压缩。在多数情况下，它们仅仅包含直接与序列有关的属性：正负链，
阅读框，字母表，位置坐标，序列本身以及它们的ID和描述。

当你对 ``HSPFragment`` 对象调用 ``print`` 方法时，这些属性可以非常简单地显示
出来。这里有个从我们BLAST搜索得到的例子：

.. code:: python

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

在这个水平上，BLAT和BLAST片段看起来非常相似，除了没有出现的query和hit序列：

.. code:: python

    >>> blat_qresult = SearchIO.read('my_blat.psl', 'blat-psl')
    >>> blat_frag = blat_qresult[0][0][0]    # first hit, first hsp, first fragment
    >>> print blat_frag
          Query: mystery_seq <unknown description>
            Hit: chr19 <unknown description>
    Query range: [0:61] (1)
      Hit range: [54204480:54204541] (1)
      Fragments: 1 (? columns)

在所有情况下，这些属性都可以通过我们最爱的点标记访问。一些例子：

.. code:: python

    >>> blast_frag.query_start      # query start coordinate
    0
    >>> blast_frag.hit_strand       # hit sequence strand
    1
    >>> blast_frag.hit              # hit sequence, as a SeqRecord object
    SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet()), id='gi|262205317|ref|NR_030195.1|', name='aligned hit sequence', description='Homo sapiens microRNA 520b (MIR520B), microRNA', dbxrefs=[])

8.2  一个关于标准和惯例的注意事项
-------------------------------------------

在我们进入到主要功能前，你需要知道 ``Bio.SearchIO`` 使用的一些标准。如果
你已经接触过多序列搜索工具，你可能必须面对每个程序处理事情方式不同的问题，
如序列位置坐标。这可能不是一个令人高兴的经历，因为这些搜索工具通常有它们
自己的标准。例如，一种工具可能使用“从1开始”(one-based)的坐标，而其他工具
使用“从0开始”(zero-based)的坐标。或者，一种程序在处理负链时，可能会反转
开始和结束坐标，而其他程序确不会。简而言之，会产生一些必须要处理的混乱。

我们意识到这种问题，并且打算在 ``Bio.SearchIO`` 中解决。毕竟， ``Bio.SearchIO`` 的目标之一就是创建一个通用简单的接口来处理多种不同的搜索
输出文件。意味着要制定一个超越你所见的对象模型的标准。

现在，你可能抱怨，”不要又来一个标准“。好吧，最后我们必须选择一个标准，这
是必须的。并且，我们并不是创造一个全新的事物；只是采用一个我们觉得对Python
使用者最好的标准（这是Biopython，毕竟）。

在使用 ``Bio.SearchIO`` 时你可以认为有个三个隐含的标准：

-  第一个适用于序列坐标。在 ``Bio.SearchIO`` 模块中，所有序列坐标遵循Python
   的坐标风格：
   从0开始，半开区间。例如，在一个BLAST XML输出文件中，HSP的起始和结束坐标
   是10和28，它们在 ``Bio.SearchIO`` 中将变成9和28。起始坐标变成9因为Python
   中索引是从0开始，而结束坐标仍然是28因为Python索引删除了区间中最后一个
   项目。
-  第二个是关于序列坐标顺序。在 ``Bio.SearchIO`` 中，开始坐标总是小于或
   等于结束坐标。但是这不是在所有的序列搜索工具中都始终适用。因为当序列
   为负链时，起始坐标会更大一些。
-  最后一个标准是关于链和阅读框的值。对于链值，只有四个可选值： ``1`` (正链)， ``-1`` (负链)， ``0`` (蛋白序列)， 和 ``None`` (无链)。对于阅读框，
   可选值是从 ``-3`` 至 ``3`` 的整型以及 ``None`` 。
   
注意，这些标准只是存在于 ``Bio.SearchIO`` 对象中。如果你把 ``Bio.SearchIO`` 
对象写入一种输出格式， ``Bio.SearchIO`` 会使用该格式的标准来输出。它并不
强加它的标准到你的输出文件。

8.3  读取搜索输出文件
--------------------------------

有两个方法，你可以用来读取搜索输出文件到 ``Bio.SearchIO`` 对象： ``read`` 和 ``parse``。
它们和其他亚模块如 ``Bio.SeqIO`` 或 ``Bio.AlignIO`` 中的 ``read`` 和 ``parse`` 方法在
本质上是相似的。你都需要提供搜索输出文件名和文件格式名，都是Python字符串类型。你可以
查阅文档来获得 ``Bio.SearchIO`` 可以识别的格式清单。

``Bio.SearchIO.read`` 用于读取单query的搜索输出文件并且返回一个 ``QueryResult`` 对象。你在前面的例子中已经看到过 ``read`` 的使用了。
你没看到的是， ``read`` 同样接受额外的关键字参数，取决于文件的格式。

这里有一些例子。在第一个例子中，我们和前面一样用 ``read`` 读BLAST表格输出
文件。在第二个例子中，我们用一个关键字来修饰，所以它分析带有注释的BLAST
表格变量。

.. code:: python

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('tab_2226_tblastn_003.txt', 'blast-tab')
    >>> qresult
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)
    >>> qresult2 = SearchIO.read('tab_2226_tblastn_007.txt', 'blast-tab', comments=True)
    >>> qresult2
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

这些关键字在不同的文件格式中是不一样的。查看格式文档，看看它是否有关键字参数来控制它的分析器行为。

对于 ``Bio.SearchIO.parse``，是用来读取含有任意数量query的搜索输出文件。
这个方法返回一个generator对象，在每次迭代中yield一个 ``QueryResult`` 对象。
和 ``Bio.SearchIO.read`` 一样，它同样接受格式特异的关键字参数：

.. code:: python

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

8.4  用索引处理含有大量搜索输出的文件
---------------------------------------------------------

有时，你得到了一个包含成百上千个query的搜索输出文件要分析，你当然可以使用
``Bio.SearchIO.parse`` 来处理，但是如果你仅仅需要访问少数query的话，效率
是及其低下的。这是因为 ``parse`` 会分析所有的query，直到找到你感兴趣。

在这种情况下，理想的选择是用 ``Bio.SearchIO.index`` 或 ``Bio.SearchIO.index_db`` 
来索引文件。如果名字听起来很熟悉，是因为你之前已经见过了，在
Section \ `5.4.2 <#sec:SeqIO-index>`__。这些方法和 ``Bio.SeqIO`` 
中相应的方法行为很相似，只是多了些格式特异的关键字参数。

这里有一些例子。你可以只用文件名和格式名来 ``index`` 

.. code:: python

    >>> from Bio import SearchIO
    >>> idx = SearchIO.index('tab_2226_tblastn_001.txt', 'blast-tab')
    >>> sorted(idx.keys())
    ['gi|11464971:4-101', 'gi|16080617|ref|NP_391444.1|']
    >>> idx['gi|16080617|ref|NP_391444.1|']
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

或者依旧使用格式特异的关键字参数：

.. code:: python

    >>> idx = SearchIO.index('tab_2226_tblastn_005.txt', 'blast-tab', comments=True)
    >>> sorted(idx.keys())
    ['gi|11464971:4-101', 'gi|16080617|ref|NP_391444.1|', 'random_s00']
    >>> idx['gi|16080617|ref|NP_391444.1|']
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

或者使用 ``key_function`` 参数，和 ``Bio.SeqIO`` 中一样：

.. code:: python

    >>> key_function = lambda id: id.upper()    # capitalizes the keys
    >>> idx = SearchIO.index('tab_2226_tblastn_001.txt', 'blast-tab', key_function=key_function)
    >>> sorted(idx.keys())
    ['GI|11464971:4-101', 'GI|16080617|REF|NP_391444.1|']
    >>> idx['GI|16080617|REF|NP_391444.1|']
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

``Bio.SearchIO.index_db`` 和 ``index`` 作用差不多，不同的只是它把query
偏移量写入一个SQLite数据库文件中。

8.5  写入和转换搜索输出文件
-----------------------------------------------

有时候，读取一个搜索输出文件，作些调整并写到一个新的文件是很有用的。
``Bio.SearchIO`` 提供了一个 ``write`` 方法，让你可以准确地完成这种工作。
它需要的参数是：一个可迭代返回 ``QueryResult`` 的对象，输出文件名，输出文件
格式和一些可选的格式特异的关键字参数。它返回一个4项目的元组，分别代表
被写入的 ``QueryResult``， ``Hit``， ``HSP``， 和 ``HSPFragment`` 对象的数量。 

.. code:: python

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('mirna.xml', 'blast-xml')     # read XML file
    >>> SearchIO.write(qresults, 'results.tab', 'blast-tab')    # write to tabular file
    (3, 239, 277, 277)

你应该注意，不同的文件格式需要 ``QueryResult``， ``Hit``， ``HSP`` 和
``HSPFragment`` 对象的不同属性。如果这些属性不存在，那么将不能写入。
也就是，你想写入的格式可能有时也会失效。举个例子，如果你读取一个BLASTXML文件，
你就不能将结果写入PSL文件，因为PSL文件需要一些属性，而这些属性BLAST却不能
提供（如重复匹配的数量）。如果你确实想写到PSL，可以手工设置这些属性。

和 ``read``， ``parse``， ``index`` 和 ``index_db`` 相似， ``write`` 同
样接受格式特异的关键字参数。查阅文档获得 ``Bio.SearchIO`` 可写格式和这些
格式的参数的完整清单。

最后， ``Bio.SearchIO`` 同样提供一个 ``convert`` 方法，可以理解为 ``Bio.SearchIO.parse`` 和 ``Bio.SearchIO.write`` 的简单替代方法。使用 ``convert`` 方法的例子如下：

.. code:: python

    >>> from Bio import SearchIO
    >>> SearchIO.convert('mirna.xml', 'blast-xml', 'results.tab', 'blast-tab')
    (3, 239, 277, 277)

因为 ``convert`` 使用 ``write`` 方法，所以只有所有需要的属性都存在时，格式
转换才能正常工作。这里由于BLAST XML文件提供BLAST 表格文件所需的所有默认值，
格式转换才能正常完成。但是，其他格式转换就可能不会正常工作，因为你需要先手工指定所需的属性。

