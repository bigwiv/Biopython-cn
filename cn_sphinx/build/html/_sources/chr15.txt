第15章 聚类分析
============================

聚类分析是根据元素相似度，进行分组的过程。在生物信息学中，聚类分析广泛
用于基因表达数据分析，用来对具有相似表达谱的基因归类；从而鉴定功能相关的基
因，或预测未知基因的功能。

Biopython中的 ``Bio.Cluster`` 模块提供了常用的聚类算法。虽然Bio.Cluster被设计用于
基因表达数据，它也可用于其他类型数据的聚类。 ``Bio.Cluster`` 
和其使用的C聚类库的说明见De Hoon *et al.* [:ref:`14 <dehoon2004>`].

``Bio.Cluster`` 包含了以下四种聚类算法：

-  系统聚类（成对重心法，最短距离，最大距离和平均连锁法);
-  *k*-means, *k*-medians, 和 *k*-medoids 聚类;
-  自组织映射（Self-Organizing Maps）;
-  主成分分析

数据表示法 
------------------------

用于聚类的输入为一个 *n* x *m* 的Python 数值矩阵 ``data``。在基因表达数据聚类中，
每一行表示不同的基因，每一列表示不同的实验条件。 ``Bio.Cluster`` 既可以
针对每行（基因），也可以针对每列（实验条件）进行聚类。

缺失值
------------------------

在芯片实验中，经常会有些缺失值，可以用一个额外的 *n* × *m* Numerical Python
整型矩阵 ``mask`` 表示。例如 ``mask[i,j]==0`` ，表示 ``data[i,j]`` 是个缺失值，
并且在分析中被忽略。

随机数生成器
------------------------

*k*-means/medians/medoids 聚类和 Self-Organizing 
Maps (SOMs) 需要调用随机数生成器。在 ``Bio.Cluster`` 中，正态分布随机数
生成器的算法是基于L’Ecuyer [:ref:`25 <lecuyer1988>`] ，二项分布的随机数
生成器算法是基于Kachitvichyanukul and Schmeiser [:ref:`23 <kachitvichyanukul1988>`] 
开发的BTPE算法。随机数生成器在调用时会首先进行初始化。由于随机数生成器使用了
两个乘同余发生器（multiplicative linear congruential generators），所以初始化时需要两个整型的
种子。这两个种子可以调用系统提供的 ``rand`` （C标准库）函数生成。在 ``Bio.Cluster`` 中，
我们首先调用 ``srand`` 使用以秒为单位的时间戳的值初始值，再用 ``rand`` 随机产生两
个随机数作为种子来产生正态分布的随机数。

.. _sec-distancefunctions:

15.1 距离函数
------------------------

为了对元素根据相似度进行聚类，第一步需要定义相似度。``Bio.Cluster`` 提供了八种不同
的距离函数来衡量相似度或者距离，分别用不同的字母代表：

-  ``'e'``: Euclidean 距离;
-  ``'b'``: City-block 距离.
-  ``'c'``: Pearson 相关系数;
-  ``'a'``: Pearson相关系数的绝对值;
-  ``'u'``: Uncentered Pearson correlation （相当于两个数据向量的夹角余弦值）
-  ``'x'``: uncentered Pearson correlation的绝对值;
-  ``'s'``: Spearman’s 秩相关系数;
-  ``'k'``: Kendall’s τ.

前两个距离函数满足三角形的两边和大于第三边的特点：

.. math::

  d\left(\underline{u},\underline{v}\right) \leq d\left(\underline{u},\underline{w}\right) + d\left(\underline{w},\underline{v}\right) \textrm{ for all } \underline{u}, \underline{v}, \underline{w},

所以称之为 *metrics*。 在任何语言中，这个意味着两点之间直线最短。

剩余的六种距离函数跟相关系数有关，距离 *d* 是由相关性 *r* 确定： *d*\ =1−\ *r*。
请注意这类距离函数是 *semi-metrics* ，因此不满足三角形的两边之和大于第三边的
性质。例如

.. math::
  
  \underline{u}=\left(1,0,-1\right);

  \underline{v}=\left(1,1,0\right);

  \underline{w}=\left(0,1,1\right);

通过计算Pearson距离，可以得到 *d*\ (*u*,\ *w*) = 1.8660, 而
*d*\ (*u*,\ *v*)+\ *d*\ (*v*,\ *w*) = 1.6340.

Euclidean 距离
~~~~~~~~~~~~~~~~~~

在 ``Bio.Cluster`` 中, Euclidean 距离被定义为

.. math::
  
  d = {\frac{1} {n}} \sum_{i=1}^{n} \left(x_i-y_i\right)^{2}.

求和时只考虑*x*\ :sub:`*i*` 和 *y*\ :sub:`*i*` 都存在的值, 分母 *n* 
也相应的做出调整。当分析表达谱数据时，由于 *x*\ :sub:`*i*` 和 *y*\ :sub:`*i*` 
会直接相减, 因此在使用Euclidean距离前，请对表达谱数据标准化处理。

City-block distance
~~~~~~~~~~~~~~~~~~~

city-block distance也称之为Manhattan 距离，跟Euclidean距离有一定的相似性。Euclidean距离
表示的是两点间最短的距离，而city-block距离则是两点在所有维度中距离的和。由于基因表达的数据
经常会有缺失数据，在 ``Bio.Cluster`` 中，city-block距离定义为总距离除以总维度：

.. math::
  
  d = {\frac{1} {n}} \sum_{i=1}^n \left|x_i-y_i\right|.

City-block distance类似于当你在从城市里一个位置到另一个位置时，所经过街道的距离。
而在Euclidean 距离中，表达谱的数据会直接相减，因此必须先对数据进行标准化。

Pearson 相关系数
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pearson相关系数定义为：

.. math::

  r = \frac{1}{n} \sum_{i=1}^n \left( \frac{x_i -\bar{x}}{\sigma_x} \right) \left(\frac{y_i -\bar{y}}{\sigma_y} \right),

其中 x, ȳ 分别是 *x* 和 *y* 的样品均值, σ\ :sub:`*x*`, σ\ :sub:`*y*` 
是 *x* 和 *y* 的样品标准差. Pearson相关系数是用于测量 *x* and *y* 散点图对直线的
拟合程度。如果所有的点都在直线上，那么Pearson相关系数为 +1 or -1, 取决于直线的斜率
是正还是负。如果Pearson 相关系数等于0，表明 *x* 和 *y* 之间没有相关性。

*Pearson distance* 定义为

.. math::
  
  d_{\textrm{P}} \equiv 1 - r.

由于Pearson 相关性的值介于 -1 和 1之间, Pearson 距离的范围为 0 和 2 之间.

Absolute Pearson correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

通过对Pearson相关系数取绝对值，可以得到一个0和1之间的数。如果绝对值是1，
所有的点都位于一条斜率为正或负直线上。当绝对值为0时，表明 *x* and *y* 没有相关性。

对应的距离定义为：

.. math::

  d_{\textrm A} \equiv 1 - \left|r\right|,

其中 *r* 是 Pearson 相关系数. 由于Pearson的相关系数的绝对值介于 0 和 1之间, 对应的
距离也位于0和1之间。

在基因表达数据分析中，应当注意，当相关性的绝对值等于1时，表明两组基因的表达情况完全一样或者完全
相反。

Uncentered correlation (夹角余弦)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

在某些情况下，使用 *uncentered correlation* 比常规的Pearson相关系数更合适。
uncentered correlation 定义为：

.. math::

  r_{\textrm U} = \frac{1}{n} \sum_{i=1}^{n} \left(\frac{x_i}{\sigma_x^{(0)}} \right) \left(\frac{y_i}{\sigma_y^{(0)}} \right),

其中

.. math::     

  \begin{eqnarray}
  \sigma_x^{(0)} & = & \sqrt{{\frac{1}{n}} \sum_{i=1}^{n}x_i^2}; \nonumber \\
  \sigma_y^{(0)} & = & \sqrt{{\frac{1}{n}} \sum_{i=1}^{n}y_i^2}. \nonumber 
  \end{eqnarray}

这个公式同Pearson相关系数的公式形式一样，只是把样本均值 x, ȳ 设为0 。
uncentered correlation 适用于表达量基准为0的情况。例如，在对基因表达分析中，使用
比值对数时，当log-ratio 等于0 表明红绿信号强度相等，也意味着实验处理
不影响基因的表达量。

uncentered correlation 系数对应的距离计算方法为：

.. math::
  
  d_{\mbox{U}} \equiv 1 - r_{\mbox{U}},

其中 *r*\ :sub:`U` 是uncentered 相关性系数。 由于uncentered系数位于-1 和 1
之间，对应的距离范围为 0 与 2之间。

由于 uncentered 相关系数值等同于两个数据向量在 *n* 维空间里的夹角余弦，因此也常称为夹角余弦。

Absolute uncentered correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

与 Pearson 相关性类似, 也可以用uncentered correlation的绝对值来定义距离:

.. math::

  d_{\mbox{AU}} \equiv 1 - \left|r_{\mbox{U}}\right|,

其中 *r*\ :sub:`U` 是 uncentered相关系数。由于uncentered 相关系数的
绝对值位于 0 和 1 之间，对应的距离也为位于 0 和 1之间。

从几何学上来讲，uncentered相关系数的绝对值等于两个数据所在向量的支持线（supporting lines）
的角度余弦值（即不考虑向量的方向性）。

Spearman rank correlation
~~~~~~~~~~~~~~~~~~~~~~~~~

Spearman秩相关系数是一种非参的相关性测量方法，对于数据中的离群点，比Pearson相关系数
有更好的稳健性。

为了计算Spearman秩相关系数，首先对每个数据集里的数据按值排序，得到每个数据的对应的
秩。然后，计算对两个数据的秩集合计算Pearson相关系数，得到Spearman的相关系数。

同Pearson相关性类似，Spearman秩相关系数对应的距离定义为：

.. math::

  d_{\mbox{S}} \equiv 1 - r_{\mbox{S}},

其中 *r*\ :sub:`S` 是Spearman秩相关系数。

Kendall’s τ
~~~~~~~~~~~

Kendall’s τ 是另一个非参的计算相关性的方法。它同Spearman秩相关系数类似，但它不对数据进行排序，
而是使用相对秩来计算  τ (see Snedecor & Cochran [:ref:`29 <snedecor1989>`] ) 。

Kendall’s τ 对应的距离计算为：

.. math::

  d_{\mbox{K}} \equiv 1 - \tau.

因为 Kendall’s τ 位于 -1 和 1之间, 对应的距离位于 0 和 2之间。

Weighting
~~~~~~~~~

对于 ``Bio.Cluster`` 中大部分距离函数，都可以使用加权向量。加权向量包含着
数据集中每个元素的权重。如果元素 *i* 的权重为 *w*\ :sub:`*i*`，那么将会认为该元素
出现了 *w*\ :sub:`*i*` 次 。权重值可以不为整数。对于 Spearman 秩相关系数
和Kendall’s τ, 权重没有太大的意义，因此不适用于这两个函数。

.. _subsec-distancematrix:

计算距离矩阵
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

距离矩阵是 ``data`` 中，所有元素的两两间的距离的平方矩阵，可以用 ``Bio.Cluster`` 模块中 ``distancematrix`` 函数计算：
 
.. code:: python

    >>> from Bio.Cluster import distancematrix
    >>> matrix = distancematrix(data)

其中，包含以下参数：

-  ``data`` (必选)
    包含所有元素的矩阵
-  ``mask`` (默认: ``None``)
    缺失数据矩阵。若 ``mask[i,j]==0``, 则 ``data[i,j]`` 缺失。若 ``mask==None``, 表明没有缺失数据。
-  ``weight`` (默认: ``None``)
    权重矩阵。若 ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择使用 ``data`` 中的行 (``transpose==0``), 或者列 (``transpose==1``)来计算距离.
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择距离函数 (具体见 :ref:`15.1 <sec-distancefunctions>` ).

为了节省内存，函数返回的距离矩阵是一个一维数组的列表。每行的列数等于行号。
因此，第一行有0个元素。例如：

.. code:: python

    [array([]),
     array([1.]),
     array([7., 3.]),
     array([4., 2., 6.])]

对应的距离矩阵为：

.. math::

  \left(
  \begin{array}{cccc}
  0 & 1 & 7 & 4  \\
  1 & 0 & 3 & 2  \\
  7 & 3 & 0 & 6  \\
  4 & 2 & 6 & 0
  \end{array}
  \right).

15.2  计算类的相关性质
------------------------------------

.. _subsec-clustercentroids:

计算类中心
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

类中心可以定义为该类中在每个维度上所有元素的平均值或者中值，可以用 ``Bio.Cluster`` 中的 ``clustercentroids`` 
函数计算：
 
.. code:: python

    >>> from Bio.Cluster import clustercentroids
    >>> cdata, cmask = clustercentroids(data)

包含了以下参数:

-  ``data`` (必选)
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    缺失数据矩阵。若 ``mask[i,j]==0``, 则 ``data[i,j]`` 缺失。若 ``mask==None``, 则明没有缺失数据。
-  ``clusterid`` (默认: ``None``)
    一个表示每个元素的所属类的整型向量。如果 ``clusterid`` 是 ``None``, 表明所有的元素属于相同的类。
-  ``method`` (默认: ``'a'``)
    指定使用算术平方根 (``method=='a'``) 或者中值(``method=='m'``) 来计算类中心。
-  ``transpose`` (默认: ``0``)
    选择使用 ``data`` 中的行 (``transpose==0``), 或者列 (``transpose==1``) 来计算类中心.

这个函数返回值为元组 ``(cdata, cmask)``。 类中心的数据存储在一个二维的Numerical Python 
数组 ``cdata`` 中, 缺失值的结果存储在二维的Numerical Python整型数组 ``cmask`` 中。 当 ``transpose`` = ``0`` 时，
这两个数组的维度是（类数，列数），当 ``transpose`` = ``1`` 时，数组的长度为 （行数，类数）。
其中每一行（当 ``transpose`` = ``0``) 或者 每一列（当 ``transpose`` = ``1`` ）
包含着对应每类对应的数据的平均值。

计算类间距离
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

根据每个 *items* 的距离函数，我们可以计算出两个 *clusters* 的距离。两个类别的
算术平均值之间的距离通常用于重心法聚类和 *k*-means 聚类，而 *k*-medoids
聚类中，通常利用两类的中值进行计算。最短距离法利用的是两类间最近的元素之间的距离，
而最大距离法利用最长的元素之间的距离。在两两平均连锁聚类法中，
类间的距离定义为类内所有对应元素两两间距离的平均值。

为了计算两类之间的距离，可以利用:

.. code:: python

    >>> from Bio.Cluster import clusterdistance
    >>> distance = clusterdistance(data)

其中，包含的参数有：

-  ``data`` (必选)
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    缺失数据矩阵。若 ``mask[i,j]==0``, 则 ``data[i,j]`` 缺失。若 ``mask==None``, 则明没有缺失数据。
-  ``weight`` (默认: ``None``)
    权重矩阵。若 ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``index1`` (默认: ``0``)
    第一个类所包含的元素索引的列表。如果一个类别只包含一个元素 *i* ，则数据类型
    可以为一个列表 ``[i]``, 或者整数 ``i``.
-  ``index2`` (默认: ``0``)
    第二个类所包含的元素的列表。如果一个类别只包含一个元素 *i* ，则数据类型
    可以为一个列表 ``[i]``, 或者整数 ``i``.
-  ``method`` (默认: ``'a'``)
    选择计算类别间距离的方法:

   -  ``'a'``: 使用两个类中心的距离 (算术平均值);
   -  ``'m'``: 使用两个类中心的距离 (中值);
   -  ``'s'``: 使用两类中最短的两个元素之间的距离;
   -  ``'x'``: 使用两类中最长的两个元素之间的距离;
   -  ``'v'``: 使用两类中对应元素间的距离的平均值作为距离。

-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择距离函数 (具体见 :ref:`15.1 <sec-distancefunctions>` ).
-  ``transpose`` (默认: ``0``)
    选择使用 ``data`` 中的行 (``transpose==0``), 或者列 (``transpose==1``)来计算距离.

15.3  划分算法
-----------------------------

划分算法依据所有元素到各自聚类中心距离之和最小化原则，
将元素分为 *k* 类。类的个数 *k* 由用户定义。 ``Bio.Cluster`` 提供了三种不同
的算法:

-  *k*-means 聚类
-  *k*-medians 聚类
-  *k*-medoids 聚类

这些算法的区别在于如何定义聚类中心。在 *k*-means 中, 聚类中心定义为该类中所有
元素的平均值。 在 *k*-medians 聚类中， 利用每个维度的中间值来计算。
最后， *k*-medoids 聚类中，聚类中心定义为该类中，距离其他所有元素距离之和最小的元素所在的位置。
这个方法适用于已知距离矩阵，但是原始数据矩阵未知的情况，例如根据结构相似度对蛋白进行聚类。

expectation-maximization (EM) 算法通常用于将数据分成 *k* 组。在 EM算法的起始阶段,
随机的把元素分配到不同的组。为了保证所有的类都包含元素，可以利用二项分布的方法随机
为每类挑选元素。然后，随机的对分组进行排列，保证每个元素有相同的概率被分到任何一个类别。
最终，保证每类中至少含有一个元素。

之后进行迭代:

-  利用均值，中值或者medoid计算每类的中心;
-  计算每类的元素离各自中心的距离;
-  对于每个元素，判别其离哪个聚类中心最近;
-  将元素重新分配到最近的聚类，当不能进行调整时，迭代终止。

为了避免迭代中产生空的类别，在 *k*-means 和 *k*-medians 聚类中，算法始终记录着每类中元素的
个数，并且阻止最后一个元素被分到其他的类别中。对于 *k*-medoids 聚类, 这种检查就是没有必要的，
因为当只剩最后一个元素时，它离中心的距离为0，所以不会被分配到其他的类别中。

由于起始阶段的每类中的元素分配是随机的，而通常当EM算法执行时，可能产生不同的聚类结果。为了找到最优的聚类结果，
可以对进行 *k*-means 算法重复多次，每次都以不同的随机分配作为起始。每次运行后，都会保存所有元素距离
其中心距离之和，并且选择总距离最小的运行结果最为最终的结果。

EM算法运行的次数取决于需要聚类元素的多少。一般而言，我们可以根据最优解被发现的次数来选择。
这个次数会作为划分算法的返回值。如果最优解被多次返回，那么不太可能存在比这个
更优的解。然后，如果最优解只被发现一次，那么可能存在着距离更小的解。但是，如果需要聚类的
元素过多的话（多余几百），那么很难找到一个全局最优解。

EM算法会在不能进行任何分配的时候停止。我们注意到，在某些随机的起始分配中，由于
相同的解会在迭代中周期性的重复，从而导致EM算法的失败。因此，我们在迭代中也会
检查是否有周期性出现的解存在。首先，在给定数目的迭代后，当前的聚类结果会保存作为一个参考。之后
继续迭代一定次数，比较该结果同之前保存的结果，可以确定之前的结果是否重复出现。
如果有重复出现，迭代会终止。如果没有出现，那么再次迭代后的结果会保存作为新的参考。
通常，会首先重复10次迭代，再保存结果为新的参考。之后，迭代的次数会翻倍，保证在长的周期中也可以
检测到该解。

*k*-means and *k*-medians
~~~~~~~~~~~~~~~~~~~~~~~~~

*k*-means 和 *k*-medians 算法可以利用 ``Bio.Cluster``中的 ``kcluster`` 实现:

.. code:: python

    >>> from Bio.Cluster import kcluster
    >>> clusterid, error, nfound = kcluster(data)

其中，包含的参数有：

-  ``data`` (必选)
    包含所有元素的矩阵。
-  ``nclusters`` (默认: ``2``)
    期望的类的数目 *k*.
-  ``mask`` (默认: ``None``)
    缺失数据矩阵。若 ``mask[i,j]==0``, 则 ``data[i,j]`` 缺失。若 ``mask==None``, 则明没有缺失数据。
-  ``weight`` (默认: ``None``)
    权重矩阵。若 ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择使用 ``data`` 中的行 (``transpose==0``), 或者列 (``transpose==1``)来计算距离.
    -  ``npass`` (默认: ``1``)
    *k*-means/-medians 聚类算法运行的次数，每次运行使用不同的随机的起始值。
    如果指定了 ``initialid`` , 程序会忽略``npass`` 的值，并且聚类算法只会运行一次。
-  ``method`` (默认: ``a``)
    指定聚类中心计算方法:

   -  ``method=='a'``: 算数平均值 (*k*-means clustering);
   -  ``method=='m'``: 中值 (*k*-medians clustering).

   当指定 ``method`` 使用其他值时，算法会采用算数平均值。
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择距离函数 (具体见 :ref:`15.1 <sec-distancefunctions>` ).
    尽管八种距离都可以用于 ``kcluster`` 计算,
    但从经验上来讲，Euclidean 距离适合 *k*-means 算法, city-block 距离适合 *k*-medians.
-  ``initialid`` (默认: ``None``)
    指定EM算法运行初始的聚类类别。如果 ``initialid==None``, 那么每运行一次EM算法时，
    都会采取不同的随机初始聚类，总共运行的次数由 ``npass`` 决定。如果 ``initialid`` 不是 ``None``,
    那么它应该为一个长度为类别数的1维数组，每类中至少含有一个元素。通常当初始分类确定后，EM算法的结果也就确定了。

这个函数的返回值为一个包含 ``(clusterid, error, nfound)`` 的元组，其中 ``clusterid`` 是
一个整型矩阵，为每行或列所在的类。 ``error`` 是最优聚类解中，每类内距离的总和，
``nfound`` 指的是最优解出现的次数。

*k*-medoids 聚类
~~~~~~~~~~~~~~~~~~~~~~

``kmedoids`` 函数根据提供的距离矩阵和聚类数，来运行 *k*-medoids 聚类：

.. code:: python

    >>> from Bio.Cluster import kmedoids
    >>> clusterid, error, nfound = kmedoids(distance)

其中，包含的参数有: , nclusters=2, npass=1,
initialid=None)\|

-  ``distance`` (必选)
    两两元素间的距离矩阵，可以通过三种不同的方法提供：

   -  提供一个2D的 Numerical Python 数组 (函数只会使用矩阵里左下角数据):

      .. code:: python

          distance = array([[0.0, 1.1, 2.3],
                            [1.1, 0.0, 4.5],
                            [2.3, 4.5, 0.0]])

   -  输入一个一维的 Numerical Python 数组，包含了距离矩阵左下角的数据：

      .. code:: python

          distance = array([1.1, 2.3, 4.5])

   -  输入一个列表，包含距离矩阵左下角的数据：

      .. code:: python

          distance = [array([]|,
                      array([1.1]),
                      array([2.3, 4.5])
                     ]

   三种方法对应着同样的距离矩阵。
-  ``nclusters`` (默认: ``2``)
    期望的类的数目 *k*.
-  ``npass`` (默认: ``1``)
    *k*-medoids 聚类算法运行的次数，每次运行使用不同的随机的起始值。
    如果指定了 ``initialid`` , ``npass`` 的值会忽略，并且聚类算法只会运行一次。
-  ``initialid`` (默认: ``None``)
    指定EM算法运行初始的聚类类别。如果 ``initialid==None``, 那么每运行一次EM算法时，
    都会采取不同的随机初始聚类，总共运行的次数由 ``npass`` 决定。如果 ``initialid`` 不是 ``None``,
    那么它应该为一个长度为类别数的1维数组，每类中至少含有一个元素。通常当初始分类确定后，EM算法的结果也就确定了。

函数返回值为一个 包含 ``(clusterid, error, nfound)`` 的元组, 其中
``clusterid`` 一个整型矩阵，为每行或列类所在的类。``error`` 是在最优解中，类内距离的总和，
``nfound`` 指的是最优解出现的次数。需要注意的是， ``clusterid`` 中的类号是指的是代表聚类中心的元素号。

15.4  系统聚类
-----------------------------

系统聚类同 *k*-means 聚类有本质的不同。在系统聚类中，基因间或者实验条件间的相似度是通过
树的形式展现出来的。由于可以利用Treeview或者Java Treeview来查看这些树的结构，因此系统聚类在基因表达谱数据中得到普遍应用。

系统聚类的第一步是计算所有元素间的距离矩阵。之后，融合两个最近的元素成为一个节点。然后，不断的
通过融合相近的元素或者节点来形成新的节点，直到所有的元素都属于同一个节点。在追溯元素和节点融合
的过程的同时形成了树的结构。不同于 *k*-means 使用的EM算法，系统聚类的过程是固定的。

系统聚类也存在着几个不同的方法，他们区别在于如何计算子节点间的距离。在
``Bio.Cluster`` 中，提供了最短距离法（ pairwise single）,最长距离法（maximum）, 类平均法（average）,
和重心法（centroid linkage）。

-  在最短距离法中，节点间的距离被定义两个节点最近样品间距离。
-  在最短距离法中，节点间的距离被定义两个节点最远样品间距离。
-  在类平均法中，节点间的距离被定义为所有样品对之间的平均距离。
-  在重心法中，节点间的距离被定义为两个节点重心间的距离。重心的计算是通过对
   每类中所有元素进行计算的。由于每次都要计算新的节点与 其他元素和已存在节点的距离，
   因此重心法的运行时间比其他系统聚类的方法更长。该方法另外一个特性是，当聚类树的
   长大的时候，距离并不会增加，有时候反而减少。这是由于使用Pearson相关系数作为距离时，
   对重心的计算和距离的计算不一致产生:因为Pearson相关系数在计算距离时会对数据进行有效归一化，，
   但是重心的计算不会存在该种归一化。

对于最短距离法，最长距离法和类平均法时，两个节点之间的距离是直接对类别里的元素计算得到的。
因此，聚类的算法在得到距离矩阵后，不一定需要提供最开始的基因表达数据。而对于重心法而言，
新生成的节点的中心必须依靠原始的数据，而不是仅仅依靠距离矩阵。

最短距离法的实现是根据 SLINK algorithm (R. Sibson, 1973), 这个算法具有快速和高效的特点。
并且这个方法聚类的结果同传统的方法结果一致。并且该算法，也可以有效的运用于大量的数据，而传统的
算法则需要大量的内存需求和运行时间。

展示系统聚类的结果
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

系统聚类的结果是用树的结构展示所有节点，每个节点包含两个元素或者子节点。通常，我们既关心那个元素
或者哪个子节点互相融合，也关心二者之间的距离（或者相似度）。我们可以调用 ``Bio.Cluster`` 中的
``Node`` 类，来存储聚类树的一个节点。 ``Node`` 的实例包含以下三个属性：

-  ``left``
-  ``right``
-  ``distance``

其中, ``left`` 和 ``right`` 是合并到该节点两个元素或子节点的编号。
``distance`` 指的是二者间的距离。其中元素的编号是从0到（元素数目-1），
而聚类的组别是从-1到-（元素数目-1）。请注意，节点的数目比元素的数目少一。


为了创建一个新的 ``Node`` 对象,我们需要指定 ``left`` 和 ``right``; 
``distance`` 是可选的。

.. code:: python

    >>> from Bio.Cluster import Node
    >>> Node(2,3)
    (2, 3): 0
    >>> Node(2,3,0.91)
    (2, 3): 0.91

对于已存在 ``Node`` 对象的 ``left``, ``right``, 和 ``distance`` 都是可以直接修改的：

.. code:: python

    >>> node = Node(4,5)
    >>> node.left = 6
    >>> node.right = 2
    >>> node.distance = 0.73
    >>> node
    (6, 2): 0.73

当 ``left`` 和 ``right`` 不是整数的时候，或者 ``distance`` 不能被转化成浮点值，会抛出错误。

 Python的类 ``Tree`` 包含着整个系统聚类的结果。 ``Tree`` 的对象可以通过
 一个 ``Node`` 的列表创建:

.. code:: python

    >>> from Bio.Cluster import Node, Tree
    >>> nodes = [Node(1,2,0.2), Node(0,3,0.5), Node(-2,4,0.6), Node(-1,-3,0.9)]
    >>> tree = Tree(nodes)
    >>> print tree
    (1, 2): 0.2
    (0, 3): 0.5
    (-2, 4): 0.6
    (-1, -3): 0.9

 ``Tree`` 的初始器会检查包含节点的列表是否是一个正确的系统聚类树的结果:

.. code:: python

    >>> nodes = [Node(1,2,0.2), Node(0,2,0.5)]
    >>> Tree(nodes)
    Traceback (most recent call last):
      File "<stdin>", line 1, in ?
    ValueError: Inconsistent tree

也可以使用中括号来对 ``Tree`` 对象进行检索：

.. code:: python

    >>> nodes = [Node(1,2,0.2), Node(0,-1,0.5)]
    >>> tree = Tree(nodes)
    >>> tree[0]
    (1, 2): 0.2
    >>> tree[1]
    (0, -1): 0.5
    >>> tree[-1]
    (0, -1): 0.5

因为 ``Tree`` 对象是只读的，我们不能对 ``Tree`` 对象中任何一个节点进行改变。然而，我们可以将其
转换成一个节点的列表，对列表进行操作，最后创建新的树。

.. code:: python

    >>> tree = Tree([Node(1,2,0.1), Node(0,-1,0.5), Node(-2,3,0.9)])
    >>> print tree
    (1, 2): 0.1
    (0, -1): 0.5
    (-2, 3): 0.9
    >>> nodes = tree[:]
    >>> nodes[0] = Node(0,1,0.2)
    >>> nodes[1].left = 2
    >>> tree = Tree(nodes)
    >>> print tree
    (0, 1): 0.2
    (2, -1): 0.5
    (-2, 3): 0.9

这个性质保证了``Tree`` 结果的正确性。

为了利用可视化工具，例如Java Treeview，来查看系统聚类树，最好对所有节点的距离进行标准化，
使其位于0和1之间。可以通过对 ``Tree`` 对象调用 ``scale`` 方法来实现这个功能：

.. code:: python

    >>> tree.scale()

这个方法不需要任何参数，返回值是 ``None``.

经过系统聚类后，可以对 ``Tree`` 对象进行剪接，将所有的元素分为 *k* 类：

.. code:: python

    >>> clusterid = tree.cut(nclusters=1)

其中 ``nclusters`` (默认是 ``1``) 是期望的类别数 *k*。这个方法会忽略树结构里面的
最高的 *k*\ −1 节点，最终形成 *k* 个独立的类别。对于 *k* 必须为正数，并且小于或者等于
元素的数目。这个方法会返回一个数组 ``clusterid`` ,包含着每个元素对应的类。

运行系统聚类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了进行系统聚类，可以用 ``Bio.Cluster`` 中的 ``treecluster`` 函数。

.. code:: python

    >>> from Bio.Cluster import treecluster
    >>> tree = treecluster(data)

包括以下参数:

-  ``data``
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    缺失数据矩阵。若 ``mask[i,j]==0``, 则 ``data[i,j]`` 缺失。若 ``mask==None``, 则明没有缺失数据。
-  ``weight`` (默认: ``None``)
    权重矩阵。若 ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择使用 ``data`` 中的行 (``transpose==0``), 或者列 (``transpose==1``)来计算距离.
-  ``method`` (默认: ``'m'``)
    选择节点间距离计算方法:

   -  ``method=='s'``: 最小距离法
   -  ``method=='m'``: 最大距离法
   -  ``method=='c'``: 重心法
   -  ``method=='a'``: 类平均法

-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择距离函数 (具体见 :ref:`15.1 <sec-distancefunctions>` ).

为了对距离矩阵进行系统聚类，可以在调用 ``treecluster`` 时，
用 ``distancematrix`` 参数来代替 ``data`` 参数：

.. code:: python

    >>> from Bio.Cluster import treecluster
    >>> tree = treecluster(distancematrix=distance)

这种情况下，需要定义下列参数：

-  ``distancematrix``
    元素两两间的距离矩阵，可以通过三种不同的方法提供：

   -  提供一个2D的 Numerical Python 数组 (函数只会使用矩阵里左下角数据):

      .. code:: python

          distance = array([[0.0, 1.1, 2.3], 
                            [1.1, 0.0, 4.5],
                            [2.3, 4.5, 0.0]])

   -  输入一个一维的 Numerical Python 数组，包含了距离矩阵左下角的数据：

      .. code:: python

          distance = array([1.1, 2.3, 4.5])

   -  输入一个列表，包含距离矩阵左下角的数据：

      .. code:: python

          distance = [array([]),
                      array([1.1]),
                      array([2.3, 4.5])

      三种方法对应着同样的距离矩阵。由于 ``treecluster`` 会对距离矩阵中的值进行随机洗牌，
      如果后面需要调用这个距离矩阵，请在调用 ``treecluster`` 之情，事先存到一个新的变量

-  ``method``
    选择节点间距离计算方法:

   -  ``method=='s'``: 最小距离法
   -  ``method=='m'``: 最大距离法
   -  ``method=='a'``: 类平均法

   其中，最小距离法、最大距离法和类平均法可以只通过距离矩阵计算，而重心法却不行。

当调用 ``treecluster``时,  ``data`` 或者 ``distancematrix`` 总有一个必须为 ``None``。

函数返回一个 ``Tree`` 对象，该对象包含着 (元素数目-1）个节点，当选择行作为聚类时，元素的
数目同行数一致；当使用列作为聚类时，元素的数目同列数一致。每个节点都意味着一对相邻连锁的
事件，其中节点的性质 ``left`` 和 ``right`` 包含着每个合并的元素或者子节点的编号， ``distance`` 
是两个合并元素或者子节点的距离。元素编号是从 0 到 (元素数目 − 1) , 而类别是从 -1 到 −(元素
数目 -1 ）

15.5  Self-Organizing Maps
--------------------------

Self-Organizing Maps (SOMs) 是由 Kohonen 在描述神经网络的时候发明的 (see for instance Kohonen, 1997 [:ref:`24 <kohonen1997>`] ).
Tamayo (1999) 第一次讲 Self-Organizing Maps 应用到基因表达数据上。
[:ref:`30 <tamayo1999>`].

SOMs 根据某种拓扑结果将元素进行分类。通常选用的是矩形的拓扑结构。在SOMs生成的类别中，相邻的
两个类的拓扑结构相似度高于他们对其他的相似度。

计算SOM的第一步是随机分配数据向量到每个类别中，如果使用行进行聚类，那么每个数据向量中的元素
个数等于列数。

一个SOM 会一次读入一行，并且找到该向量最近的拓扑聚类结构。之后利用找到的数据向量对
这个类别的数据向量和相邻的类别的数据向量进行调整。调整如下：

.. math::

  \Delta \underline{x}_{\textrm{cell}} = \tau \cdot \left(\underline{x}_{\textrm{row}} - \underline{x}_{\textrm{cell}} \right).

参数 τ 会随着迭代次数增加而减少。可以用一个简单的线性函数来定义其与迭代次数的关系：

.. math::

  \tau = \tau_{\textrm{init}} \cdot \left(1 - {\frac{1}{n}}\right),

τ\ :sub:`init` 是指定的起始的 τ 值， *i* 是当前迭代的次数， *n* 是总的需要迭代的次数。
在迭代开始时，τ变化很快，然而在迭代末尾，变化越来越小。

所有在半径 *R* 内的类别都会在每次迭代中进行调整。半径也会随着迭代的增加而减小：

.. math::

  R = R_{\textrm{max}} \cdot \left(1 - {\frac{1}{n}}\right),

其中最大的半径定义为：

.. math::

  R_{\textrm{max}} = \sqrt{N_x^2 + N_y^2},

其中 (*N*\ :sub:`*x*`, *N*\ :sub:`*y*`) 是定义拓扑结构的矩形维度。

函数 ``somcluster`` 可以用来在一个矩形的网格里计算 Self-Organizing Map。
首先，初始化一个随机数产生器。利用随机化产生器来对节点数据进行初始化。在SOM中，
基因或者芯片的调整顺序同样是随机的。用户可以定义总的SOM迭代的次数。

运行 ``somcluster``, 例如：

.. code:: python

    >>> from Bio.Cluster import somcluster
    >>> clusterid, celldata = somcluster(data)

其中，可以定义一下参数:

-  ``data`` (required)
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    缺失数据矩阵。若 ``mask[i,j]==0``, 则 ``data[i,j]`` 缺失。若 ``mask==None``, 则明没有缺失数据。
-  ``weight`` (默认: ``None``)
    权重矩阵。若 ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择使用 ``data`` 中的行 (``transpose==0``), 或者列 (``transpose==1``)来聚类.
-  ``nxgrid, nygrid`` (默认: ``2, 1``)
    当Self-Organizing Map计算的时候，矩形的网格所包含的横向和纵向的格子。
-  ``inittau`` (默认: ``0.02``)
    SOM算法中，参数 τ 的初始值，默认是 0.02。 这个初始值同Michael Eisen’s Cluster/TreeView 一致。
-  ``niter`` (默认: ``1``)
    迭代运行的次数。
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择距离函数 (具体见 :ref:`15.1 <sec-distancefunctions>` ).

这个函数返回的是一个元组 ``(clusterid, celldata)``:

-  ``clusterid``:
    一个两列的数组，行的数目等于待聚类元素的个数。每行包含着在矩形SOM网格中，将每个元素分配到的
    格子的 *x* 和 *y* 的坐标。
-  ``celldata``:
    当以行进行聚类时，生成的矩阵维度为 (``nxgrid``, ``nygrid``, number of columns)；
    当以列进行聚类时，生成的矩阵维度为 (``nxgrid``, ``nygrid``, number of  rows)。
    在这个矩阵里， ``[ix][iy]`` 表示着一个一维向量，其中用于计算该类中心的这基因的表达谱数据.

15.6  主成分分析
----------------------------------

主成分分析 (PCA) 被广泛的用于分析多维数据，一个将主成分分析应用于表达谱数据的请见
Yeung and Ruzzo (2001) [:ref:`33 <yeung2001>`].

简而言之，PCA是一种坐标转换的方法，转换后的基础向量成为主成分，变换前的每行可以用主成分的
线性关系显示。主成分的选择是基于是残差尽可能的小的原则。例如，一个 *n* × 3 的数据矩阵可以表示为三维
空间内的一个椭圆球形的点的云。第一主成分是这个椭圆球形的最长轴，第二主成分是次长轴，第三主成分
是最短的轴。矩阵中，每一行都可以用主成分的线性关系展示。一般而言，为了对数据进行降维，只保留最
重要的几个主成分。剩余的残差认为是不可解释的方差。

可以通过计算数据的协方差矩阵的特征向量来得到主成分。每个主成分对应的特征值决定了
其在数据中代表的方差的大小。

在进行主成分分析前，矩阵的数据每一列都要减去其平均值。在上面椭圆球形云的例子中，数据在3D
空间中，围绕着其中心分布，而主成分则显示着每个点对其中心的变化。

函数 ``pca`` 首先使用奇异值分解（singular value decomposition）来计算矩阵的特征值和
特征向量。奇异值分解使用的是Algol写的C语言的 ``svd`` [:ref:`16 <golub1971>`] , 利用的是
Householder bidiagonalization 和 QR 算法的变异。主成分，每个数据在主成分上的坐标和主成分
对应的特征值都会被计算出来，并按照特征值的降序排列。如果需要数据中心，则需要在调用 ``pca`` 
前，对每列数据减去其平均值。

将主成分分析应用于二维矩阵 ``data``,可以：

.. code:: python

    >>> from Bio.Cluster import pca
    >>> columnmean, coordinates, components, eigenvalues = pca(data)

函数会返回一个元组 ``columnmean, coordinates, components, eigenvalues`` :

-  ``columnmean``
    包含 ``data`` 每列均值的数组 .
-  ``coordinates``
    ``data`` 中每行数据在主成分上对应的坐标。
-  ``components``
    主成分
-  ``eigenvalues``
    每个主成分对应的特征值

原始的数据 ``data`` 可以通过计算 ``columnmean +  dot(coordinates, components)`` 得到。

15.7  处理 Cluster/TreeView-type 文件
------------------------------------------

Cluster/TreeView 是一个对基因表达数据可视化的工具。他们最初由 `Michael
Eisen <http://rana.lbl.gov>`__ 在 Stanford University 完成。``Bio.Cluster`` 
包含着读写 Cluster/TreeView 对应的文件格式的函数。因此，将结果保存为该格式后，
可以用Treeview对结果进行直接的查看。我们推荐使用 Alok Saldanha 的
`http://jtreeview.sourceforge.net/ <http://jtreeview.sourceforge.net/>`__\ Java
TreeView 程序。这个软件可以显示系统聚类和 *k*-means 聚类的结果。

类 ``Record`` 的一个对象包含着一个 Cluster/TreeView-type数据文件需要的所有信息。
为了将结果保存到一个 ``Record`` 对象中，首先需要打开一个文件，并读取：

.. code:: python

    >>> from Bio import Cluster
    >>> handle = open("mydatafile.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()

两步操作使得你可以较灵活地操作不同来源的数据，例如：

.. code:: python

    >>> import gzip # Python standard library
    >>> handle = gzip.open("mydatafile.txt.gz")

来打开一个gzipped文件，或者利用

.. code:: python

    >>> import urllib # Python standard library
    >>> handle = urllib.urlopen("http://somewhere.org/mydatafile.txt")

来打开一个网络文件，然后调用 ``read``.

``read`` 命令会读取一个由制表符分割的文本文件 ``mydatafile.txt``，文件包含着
符合Michael Eisen’s Cluster/TreeView格式的基因表达数据。具体的格式说明，可以参见
Cluster/TreeView手册，链接见 `Michael Eisen’s lab
website <http://rana.lbl.gov/manuals/ClusterTreeView.pdf>`__ 或者 `our
website <http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/cluster3.pdf>`__.

一个 ``Record`` 对象有以下的性质:

-  ``data``
    包含基因表达数据的矩阵，每行为基因，每列为芯片。
-  ``mask``
    缺失值的整型数组。如果 ``mask[i,j]==0``, 则 ``data[i,j]`` 是缺失的. 如果 ``mask==None``,
    那么没有数据缺失。
-  ``geneid``
    包含每个基因的独特说明的列表 (例如 ORF 数目).
-  ``genename``
    包含每个基因说明的列表（例如基因名）。如果文件中不包含该数据，
    那么 ``genename`` 被设为 ``None``.
-  ``gweight``
    计算表达谱数据中，基因间的距离使用的权重。如果文件中不含该信息，则
    ``gweight`` 为 ``None``.
-  ``gorder``
    期望输出文件中基因的排列的顺序。如果文件中不含该信息，则
    ``gorder`` 为``None``.
-  ``expid``
    包含每个芯片说明的列表，例如实验条件。
-  ``eweight``
    计算表达谱数据中，不同芯片间的距离使用的权重。如果文件中不含该信息，则
    ``eweight`` 为 ``None``.
-  ``eorder``
    期望输出文件中基因的排列的顺序。如果文件中不含该信息，则 ``eorder`` 为  ``None``.
-  ``uniqid``
    用于代替文件中 UNIQID 的字符串.

在载入 ``Record`` 对象后，上述的每个性质可以直接读取和修改。例如，可以对
``record.data`` 直接取对数来对数据进行log转换。

计算距离矩阵
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了计算record中存储元素的距离矩阵，可以用：

.. code:: python

    >>> matrix = record.distancematrix()

其中，包含以下参数：

-  ``transpose`` (默认: ``0``)
    选择对 ``data`` 的行 (``transpose==0``), 或者列 (``transpose==1``)计算距离。
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择合适的元素距离算法 (见 :ref:`15.1 <sec-distancefunctions>` ).

函数会返回一个距离矩阵，每行的列数等于行数。(见 :ref:`15.1 <subsec-distancematrix>` ).

计算聚类中心
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了计算存储在record中的元素的聚类中心，利用：

.. code:: python

    >>> cdata, cmask = record.clustercentroids()

-  ``clusterid`` (默认: ``None``)
    展示每个元素所属类的整型向量。如果缺少 ``clusterid``,默认所有的元素属于同一类。
-  ``method`` (默认: ``'a'``)
    选择使用算术平均值 (``method=='a'``) 或者中值 (``method=='m'``)来计算聚类中心。
-  ``transpose`` (默认: ``0``)
    选择计算``data`` 的行 (``transpose==0``), 或者列 (``transpose==1``)计算中心。

函数返回元组 ``cdata, cmask`` ; 见 :ref:`15.2 <subsec-clustercentroids>` .

计算两类间的距离
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了计算存储在record中的两类的距离，利用：

.. code:: python

    >>> distance = record.clusterdistance()

其中，包含以下参数：

-  ``index1`` (默认: ``0``)
    第一个类别所包含的元素的列表。如果一个类别只包含一个元素 *i* 
    可以为一个列表 ``[i]``, 或者整数 ``i``.
-  ``index2`` (默认: ``0``)
    第二个类别所包含的元素的列表。如果一个类别只包含一个元素 *i* 
    可以为一个列表 ``[i]``, 或者整数 ``i``.
-  ``method`` (默认: ``'a'``)
    选择计算类别间距离的方法:

   -  ``'a'``: 使用两个聚类中心的距离 (算术平均值);
   -  ``'m'``: 使用两个聚类中心的距离 (中值);
   -  ``'s'``: 使用两类中最短的两个元素之间的距离;
   -  ``'x'``: 使用两类中最长的两个元素之间的距离;
   -  ``'v'``: 使用两类中两两元素距离的平均值作为距离。

-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择使用的距离函数 (见 :ref:`15.1 <sec-distancefunctions>` ).
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行 ( ``transpose==0`` ), 或者列 ( ``transpose==1`` )计算距离。

进行系统聚类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了对存储在record中的数据进行系统聚类，利用：

.. code:: python

    >>> tree = record.treecluster()

包含以下参数:

-  ``transpose`` (默认: ``0``)
    选择使用行 ( ``transpose==0`` ) 或者列 ( ``transpose==1`` ) 用于聚类
-  ``method`` (默认: ``'m'``)
    选择合适的节点距离计算方法:

   -  ``method=='s'``: 最小距离法
   -  ``method=='m'``: 最大距离法
   -  ``method=='c'``: 重心法
   -  ``method=='a'``: 类平均法

-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择使用的距离函数(见 :ref:`15.1 <sec-distancefunctions>` ).
-  ``transpose``
    选择使用基因或者芯片进行聚类，如果是 ``transpose==0`` , 则使用基因 (行) 进行聚类，如果使用
    ``transpose==1``, 芯片 (列) 用于聚类.

函数返回 ``Tree`` 对象。对象包含 (元素数目 − 1） 节点, 如果使用行进行聚类时，元素数目为总行数；
当使用列进行聚类时，元素数目为总列数。每个节点描述着一对节点连接，然而节点的性质 ``left`` 和
``right`` 包含着相邻节点所有的元素和子节点数， ``distance`` 显示着左右节点的距离。
元素从 0 到 (元素数目 − 1) 进行索引, 而类别从 -1 to −(元素数目−1)进行索引。

进行 *k*-means or *k*-medians 聚类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了对存储在record中的元素进行 *k*-means 或者 *k*-medians 聚类，可以使用：

.. code:: python

    >>> clusterid, error, nfound = record.kcluster()

包含以下参数:

-  ``nclusters`` (默认: ``2``)
    类的数目 *k*.
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行 ( ``transpose==0`` ), 或者列 ( ``transpose==1`` )计算距离。
-  ``npass`` (默认: ``1``)
    *k*-means/-medians 聚类算法运行的次数，每次运行使用不同的随机的起始值。
    如果指定了 ``initialid`` , ``npass`` 的值会忽略，并且聚类算法只会运行一次。
-  ``method`` (默认: ``a``)
    指定确定聚类中心的方法:

   -  ``method=='a'``: 算数平均值 (*k*-means clustering);
   -  ``method=='m'``: 中间值 (*k*-medians clustering).

   当指定 ``method`` 使用其他值时，算法会采用算数平均值。
-  ``dist`` (默认: ``'e'`` , Euclidean distance)
    选择使用的距离函数 (见 :ref:`15.1 <sec-distancefunctions>` ).

这个函数返回的是一个元组 ``(clusterid, error, nfound)`` , 其中 ``clusterid`` 是一个每行或则列对应的类的编号。
``error`` 是最优解的类内的距离和， ``nfound`` 是最优解被发现的次数。

计算Self-Organizing Map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

可以利用以下命令，计算对存储在record中元素计算 Self-Organizing Map ：

.. code:: python

    >>> clusterid, celldata = record.somcluster()

包含以下参数:

-  ``transpose`` (默认: ``0`` )
    选择 使用 ``data`` 的行 ( ``transpose==0`` ), 或者列 ( ``transpose==1`` )计算距离.
-  ``nxgrid, nygrid`` (默认: ``2, 1``)
    当Self-Organizing Map计算时，在矩形网格里的横向和纵向格子数目
-  ``inittau`` (默认: ``0.02``)
    用于SOM算法的参数 τ 的初始值。默认的 ``inittau`` 是0.02，同Michael Eisen’s Cluster/TreeView 程序中
    使用的参数一致。
-  ``niter`` (默认: ``1`` )
    迭代运行的次数。
-  ``dist`` (默认: ``'e'`` , Euclidean distance)
    选择使用的距离函数(见 :ref:`15.1 <sec-distancefunctions>` ).

函数返回一个元组 ``(clusterid, celldata)`` :

-  ``clusterid``:
    一个二维数组，行数同待聚类的元素数目相同。每行的内容对应着该元素在矩形SOM方格内 *x* 和 *y* 的坐标。
-  ``celldata``:
    格式为一个矩阵，如果是对行聚类，内容为 ( ``nxgrid`` , ``nygrid`` , 列数)，如果是对列聚类，
    那么内容为 ( ``nxgrid`` , ``nygrid`` , 行数) 。矩阵中，坐标 ``[ix][iy]`` 对应的是该坐标的网格里的
    基因表达数据的聚类中心的一维向量。

保存聚类结果
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了保存聚类结果，可以利用：

.. code:: python

    >>> record.save(jobname, geneclusters, expclusters)

包含以下参数:

-  ``jobname``
    字符串 ``jobname`` 作为保存的文件名。
-  ``geneclusters``
    这个参数指的是基因（以行聚类）的结果。在 *k*-means 聚类中，这个参数是一个一维的数组，包含着
    每个基因对应的类别，可以通过 ``kcluster`` 得到。在系统聚类中， ``geneclusters`` 是一个 ``Tree`` 对象。
-  ``expclusters``
    这个参数指的是实验条件（以列聚类）的结果。在 *k*-means 聚类中，这个参数是一个一维的数组，包含着
    每个实验条件对应的类别，可以通过 ``kcluster`` 得到。在系统聚类中， ``geneclusters`` 是一个``Tree`` 对象。

这个方法会生成文本文件 ``jobname.cdt``, ``jobname.gtr``, ``jobname.atr``, ``jobname*.kgg``, 
和/或 ``jobname*.kag`` 。 这些文件可以用于后续分析。如果 ``geneclusters`` 和 ``expclusters`` 
都是 ``None`` , 那这个方法只会生成 ``jobname.cdt`` ; 这个文件可以被读取，生成一个新的 ``Record`` 对象.

15.8  示例
-------------------------

以下是一个系统聚类的例子，其中使用最短距离法对基因进行聚类，用最大距离法对实验条件进行聚类。
由于使用 Euclidean 距离对基因进行聚类，因此需要将节点距离 ``genetree`` 进行调整，使其处于0和1之间。
这种调整对于Java TreeView正确显示树结构也是很必须的。同时使用 uncentered correlation 对实验条件进行聚类。
在这种情况下，不需要任何的调整，因为 ``exptree`` 中的结果已经位于0和2之间。 示例中使用的
文件 ``cyano.txt`` 可以从 ``data`` 文件夹中找到。

.. code:: python

    >>> from Bio import Cluster
    >>> handle = open("cyano.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()
    >>> genetree = record.treecluster(method='s')
    >>> genetree.scale()
    >>> exptree = record.treecluster(dist='u', transpose=1)
    >>> record.save("cyano_result", genetree, exptree)

这个命令会生成 ``cyano_result.cdt`` , ``cyano_result.gtr`` , 和 ``cyano_result.atr`` 等文件。

同样的，也可以保存一个 *k*-means 聚类的结果:

.. code:: python

    >>> from Bio import Cluster
    >>> handle = open("cyano.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()
    >>> (geneclusters, error, ifound) = record.kcluster(nclusters=5, npass=1000)
    >>> (expclusters, error, ifound) = record.kcluster(nclusters=2, npass=100, transpose=1)
    >>> record.save("cyano_result", geneclusters, expclusters)

上述代码将生成文件 ``cyano_result_K_G2_A2.cdt`` , ``cyano_result_K_G2.kgg`` , 和 ``cyano_result_K_A2.kag`` 。

15.9  附加函数
-------------------------

``median(data)`` 返回一维数组 ``data`` 的中值

``mean(data)`` 返回一维数组 ``data`` 的均值。

``version()`` 返回使用的C聚类库的版本号。

