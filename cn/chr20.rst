第20章 高级
====================

20.1  解析器的设计
-------------------

过去很多Biopython解析器都是根据面向事件设计出来的，包括Scanner和Consumer。

Scanners是将输入的数据源进行逐行分析，只要识别出数据中的信息就会发送一个事件。
例如，如果数据中包含生物名的信息，scanner只要读到某行包含名称信息时就会产生一个“生物名”的事件。

Consumers是用来接收Scanners所发出事件的对象。
接着上面的例子，当consumer收到了“生物名”事件，在当前应用程序中无论以何种方式都会运行。

这是一个非常灵活的构架，如果你想要将一个文件解析成多种其他格式的话，这会是很有优势的。
例如，“Bio.GenBank”模块可以运用这种方式构建“SeqRecord”或者其他独特的文件格式记录对象。

最近，很多添加了“Bio.SeqIO”和“Bio.AlignIO”的解析器在使用一种更为简单的方法，
但是只能产生单一形式的文件格式（分别是“SeqRecord”和“MultipleSeqAlignment”）。
在某些情况，“Bio.SeqIO”解析器实际上包含了另一种Biopython解析器。
例如，“Bio.SwissProt”解析器产生了特定的SwissProt格式对象，又转换成了“SeqRecord”格式对象。


20.2  替换矩阵
---------------------------

20.2.1  SubsMat
~~~~~~~~~~~~~~~

这个模块提供了一个类和一些固定的方法来产生替换矩阵，类似于BLOSUM或者PAM矩阵，但是是基于用户提供的数据。
此外，你还可以从已建立的替换矩阵集合MatrixInfo.py中选择一个矩阵。
“SeqMat”类来自于一个dictionary库:

.. code:: verbatim

    class SeqMat(dict)

这个dictionary的格式是
``{(i1,j1):n1, (i1,j2):n2,...,(ik,jk):nk}`` i和j是字母，而n是一个值。

#. 属性

   #. ``self.alphabet``: Bio.Alphabet中定义的一个类
   #. ``self.ab_list``: 排列好的字母列表。主要是内部原因而需要。

#. 方法

   #. .. code:: verbatim

          __init__(self,data=None,alphabet=None, mat_name='', build_later=0):

      #. ``data``: 可以是一个dictionary，也可以是另一个SeqMat实例。
      #. ``alphabet``: 一个Bio.Alphabet的实例。如果没有提供，可以从数据构建一个alphabet。
      #. ``mat_name``: 矩阵名，例如 "BLOSUM62" 或者 "PAM250"
      #. ``build_later``: 默认值为假。如果是真，用户应该只提供alphabet和空的dictionary。
	     如果想要之后再构建矩阵，这样会跳过alphabet大小和矩阵大小的检查。
	     
   #. .. code:: verbatim

          entropy(self,obs_freq_mat)

      #. ``obs_freq_mat``: 一个观测频率矩阵。基于“obs_freq_mat”的频率返回矩阵的熵值。矩阵实例应该是LO或者SUBS。

   #. .. code:: verbatim

          sum(self)

      计算矩阵的alphabet中每个字母值的总和，返回值是以dictionary的形式
      ``{i1: s1, i2: s2,...,in:sn}``, 其中:

      -  i: 一个字母;
      -  s: 半矩阵中某个字母值的总和;
      -  n: alphabet中字母的个数。

   #. .. code:: verbatim

          print_mat(self,f,format="%4d",bottomformat="%4s",alphabet=None)

      将矩阵打印到文件句柄f。“format”是矩阵值的格式；“bottomformat”是底部行的格式，包括矩阵字母。
      下面是一个3字母矩阵的例子：

      .. code:: verbatim

          A 23
          B 12 34
          C 7  22  27
            A   B   C

      ``alphabet``可选的自变量是alphabet中所有字符的一个字符串。
	  如果用户提供了数据，则轴上字母的顺序是根据字符串中的顺序，而不是字母表的顺序。

#. 用法
   
   安排下面这个部分是因为大多数的读者希望能够知道如何产生一个对数机率矩阵。
   当然，也有人研究产生的中间过渡矩阵。
   但是大部分的人只是想要一个对数机率矩阵，仅此而已。

   #. 产生一个可接受的代替矩阵

      首先，你应该从数据中产生出一个可接受代替矩阵（ARM）。
	  ARM中的数值是根据你的数据中替换的个数决定的。数据可以是一对或者多对的序列比对结果。
	  例如，丙氨酸被半胱氨酸替换了10次，而半胱氨酸被丙氨酸替换了12次，其相对应的ARM为：

      .. code:: verbatim

          ('A','C'): 10, ('C','A'): 12

      由于顺序并不重要，用户也可以只用一个输入:

      .. code:: verbatim

          ('A','C'): 22

      一个SeqMat实例的初始化可以用全矩阵（第一种计数方法：10,12），也可以用半矩阵（后一种方法，22）。
	  一个蛋白字母全矩阵的大小应该是20x20 = 400。而一个这样的半矩阵大小是20x20/2 + 20/2 = 210。
	  这是因为相同字母的输入并没有改变（矩阵的对角线）。如果一个大小为N的alphabet：

      #. 全矩阵大小: N\*N
      #. 半矩阵大小: N(N+1)/2

      如果已经通过了全矩阵，SeqMat的构造器会自动产生半矩阵。
	  如果通过了半矩阵，则键的字母将按照字母表顺序排列('A','C')，而不是('C','A')。
	  
      讲到这里，如果你想知道的仅仅只是怎样产生一个对数机率矩阵的话，请直接看用法示例那个章节。
	  对于想要更加深入地知道核苷酸/氨基酸频率数据的读者，接下来要讲的是内部功能的细节。


   #. 产生观察频率矩阵(OFM)

      用法:

      .. code:: verbatim

          OFM = SubsMat._build_obs_freq_mat(ARM)

      观察频率矩阵是由可接受代替矩阵产生的，只是将替换的个数换成了替换频率。

   #. 产生一个期望频率矩阵(EFM)

      用法:

      .. code:: verbatim

          EFM = SubsMat._build_exp_freq_mat(OFM,exp_freq_table)

      #. ``exp_freq_table``: 应该是一个FreqTable的实例。
	     更加详细的信息可以查看章节\ `20.2.2 <#sec:freq_table>`__ 
	     简单地说，期望频率表表示字母表中每个元素显示的频率。
		 这个表相当于一个dictionary，字母是键，字母对应的频率是值。总值是1。
		
      期望频率表可以（理论上说应该）从观察频率矩阵得到。
	  所以大多数情况你可以用下面的代码产生``exp_freq_table``:

      .. code:: verbatim

          >>> exp_freq_table = SubsMat._exp_freq_table_from_obs_freq(OFM)
          >>> EFM = SubsMat._build_exp_freq_mat(OFM,exp_freq_table)

      如果你想的话，也可以使用自己提供的``exp_freq_table``

   #. 产生一个替换频率矩阵(SFM)

      用法:

      .. code:: verbatim

          SFM = SubsMat._build_subs_mat(OFM,EFM)

      使用观察频率矩阵(OFM)和期望频率矩阵(EFM)。
	  得到相应值的除法结果。

   #. 产生一个对数机率矩阵(LOM)

      用法:

      .. code:: verbatim

          LOM=SubsMat._build_log_odds_mat(SFM[,logbase=10,factor=10.0,round_digit=1])

      #. 使用一个替换频率矩阵(SFM)。
      #. ``logbase``: 用来产生对数机率值的对数的底。
      #. ``factor``: 因子是对数机率值的乘数。
	     每个数通过log(LOM[key])\*factor产生，如果需要还可以通过``round_digit``四舍五入到小数点后。
	   
#. 用法示例

   正如大部分的人都想用最简单的方法产生对数机率矩阵，SubsMat提供了一个函数可以完成所有的事情：
 
   .. code:: verbatim

       make_log_odds_matrix(acc_rep_mat,exp_freq_table=None,logbase=10,
                             factor=10.0,round_digit=0):

   #. ``acc_rep_mat``: 用户提供可接受代替矩阵
   #. ``exp_freq_table``: 期望频率表。如果提供了就使用期望频率表，否则就从``acc_rep_mat``产生。
   #. ``logbase``: 对数机率矩阵的对数的底。默认底为10。
   #. ``round_digit``: 四舍五入的小数点位数。默认为0。

20.2.2  FreqTable
~~~~~~~~~~~~~~~~~

.. code:: verbatim

    FreqTable.FreqTable(UserDict.UserDict)

#. 属性:

   #. ``alphabet``: 一个Bio.Alphabet的实例。
   #. ``data``: 频率dictionary
   #. ``count``: 计数dictionary(如果有计数的话)。

#. 功能:

   #. ``read_count(f)``: 从字符串f读入一个计数文件。然后将其转换成频率。
   #. ``read_freq(f)``: 从字符串f读入一个频率数据文件。
      当然，我们不用计数，我们感兴趣的是字母频率。

#. 用法示例: 文件中有残基的个数，用空格分格，形式如下（以3个字母为例）：

   .. code:: verbatim

       A   35
       B   65
       C   100

   用``FreqTable.read_count(file_handle)``方法读入。
   
   一个等价的频率文件:

   .. code:: verbatim

       A  0.175
       B  0.325
       C  0.5

   相反地，残基频率或者计数也可以作为dictionary输入。
   一个计数dictionary的例子（3个字母）：

   .. code:: verbatim

       {'A': 35, 'B': 65, 'C': 100}

   这也意味着'C'的频率是0.5，'B'的频率是0.325，'A'的频率是0.175。A、B、C的总和为200。

   一个相同数据的频率dictionary如下：

   .. code:: verbatim

       {'A': 0.175, 'B': 0.325, 'C': 0.5}

   总和为1。

   当传入一个dictionary数据作为参数，应该指出这是一个计数还是频率的dictionary。
   因此FreqTable类的构建需要两个参数：
   dictionary本身和FreqTable.COUNT或者FreqTable.FREQ，分别代表计数或者频率。

   读入期望的计数，readCount会产生频率。下面的任意一个都可以用来产生频率表（ftab）：

   .. code:: verbatim

       >>> from SubsMat import *
       >>> ftab = FreqTable.FreqTable(my_frequency_dictionary,FreqTable.FREQ)
       >>> ftab = FreqTable.FreqTable(my_count_dictionary,FreqTable.COUNT)
       >>> ftab = FreqTable.read_count(open('myCountFile'))
       >>> ftab = FreqTable.read_frequency(open('myFrequencyFile'))
