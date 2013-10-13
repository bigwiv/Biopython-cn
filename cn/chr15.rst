第十五章 聚类分析
============================

聚类分析是对一堆元素按照相似度来进行分组的过程。在生物信息中，聚类分析经常
用于分析基因表达的数据，用来发现有相同的表达图谱的基因，以鉴定功能相关的基
因，或者发现未知的功能相关基因。

Biopython中的 ``Bio.Cluster`` 模块提供了常用的聚类算法，并且设计的时候考虑了
基因表达数据的应用。然而，这个模块也可以用于其他类型数据。 ``Bio.Cluster`` 
和其使用的C聚类库的说明见De Hoon *et al.* [`14 <#dehoon2004>`__\ ].

``Bio.Cluster`` 包含了以下四种聚类算法：

-  分层聚类（两两聚类，single-，complete-, and average-linkage);
-  *k*-means, *k*-medians, and *k*-medoids clustering;
-  Self-Organizing Maps;
-  主成分分析

数据结构 
~~~~~~~~~~~~~~~~~~~

用于聚类的数据为一个 *n* x *m* 的Numerical Python 矩阵 ``data``。
其中，每一行表示不同的基因，每一列表示不同的实验条件。 ``Bio.Cluster`` 既可以
针对每行（基因），也可以针对每列（实验条件）进行聚类。

缺失值
~~~~~~~~~~~~~~


在芯片实验中，经常会有些缺失值，通常用一个额外的 *n* × *m* Numerical Python
整型矩阵 ``mask`` 表示。例如 ``mask[i,j]==0`` ，表示 ``data[i,j]`` 是个缺失值，
并且在分析中忽略。

随机数生成器
~~~~~~~~~~~~~~~~~~~~~~~
*k*-means/medians/medoids clustering algorithms and Self-Organizing
Maps (SOMs) 需要调用随机数生成器。 ``Bio.Cluster`` 中使用的正态分布随机数
生成器使用的算法是基于L’Ecuyer [`25 <#lecuyer1988>`__\ ],二项分布的随机数
生成是根据Kachitvichyanukul and Schmeiser [`23 <#kachitvichyanukul1988>`__\ ]
开发的BTPE算法。随机数生成器在调用时首先进行初始化。由于随机数生成器使用了
两个multiplicative linear congruential generators，所以初始化时需要两个整型的
种子，并调用系统提供的 ``rand`` （C标注库）生成随机数。在 ``Bio.Cluster`` 中，
也可以调用 ``srand`` 使用当前时间的秒作为初始值，用 ``rand`` 随机产生的头两
个随机数作为种子来产生正态分布的随机数。


15.1 距离函数
------------------------
为了根据相似度进行聚类，首先需要定义相似度。``Bio.Cluster``提供了八种不同
的距离函数来计算相似度或者距离，分别用不同的字母代表：

-  ``'e'``: Euclidean 距离;
-  ``'b'``: City-block 距离.
-  ``'c'``: Pearson 相关系数;
-  ``'a'``: Pearson相关系数的绝对值;
-  ``'u'``: Uncentered Pearson correlation (相当于两个数据向量的形成角度的cos值
-  ``'x'``: 绝对的uncentered Pearson correlation;
-  ``'s'``: Spearman’s rank correlation;
-  ``'k'``: Kendall’s τ.

前两个距离函数满足三角形的两边和大于第三边的特点：


*d*

| ⎛
|  ⎜
|  ⎝

+-------+
| *u*   |
+-------+
+-------+

,

+-------+
| *v*   |
+-------+
+-------+

| ⎞
|  ⎟
|  ⎠

≤ \ *d*

| ⎛
|  ⎜
|  ⎝

+-------+
| *u*   |
+-------+
+-------+

,

+-------+
| *w*   |
+-------+
+-------+

| ⎞
|  ⎟
|  ⎠

+ \ *d*

| ⎛
|  ⎜
|  ⎝

+-------+
| *w*   |
+-------+
+-------+

,

+-------+
| *v*   |
+-------+
+-------+

| ⎞
|  ⎟
|  ⎠

for all  

+-------+
| *u*   |
+-------+
+-------+

, 

+-------+
| *v*   |
+-------+
+-------+

, 

+-------+
| *w*   |
+-------+
+-------+

,

所以称之为 *metrics*. 在任何语言中，这个意味着两点之间直线最短。

剩余的六中距离函数同相关系数有关，距离 *d* 是有相关性*r* 确定： *d*\ =1−\ *r*。
请注意这类距离函数是 *semi-metrics* ，因此不满足三角形的两边之和大于第三边的
性质。例如


+-------+
| *u*   |
+-------+
+-------+

=

| ⎛
|  ⎝

1,0,−1

| ⎞
|  ⎠

;

+-------+
| *v*   |
+-------+
+-------+

=

| ⎛
|  ⎝

1,1,0

| ⎞
|  ⎠

;

+-------+
| *w*   |
+-------+
+-------+

=

| ⎛
|  ⎝

0,1,1

| ⎞
|  ⎠

;

计算Pearson距离 *d*\ (*u*,\ *w*) = 1.8660, 而
*d*\ (*u*,\ *v*)+\ *d*\ (*v*,\ *w*) = 1.6340.

Euclidean 距离
~~~~~~~~~~~~~~~~~~

在 ``Bio.Cluster`` 中, 定义 Euclidean 距离为

*d* = 

+-------+
| 1     |
+-------+
+-------+
| *n*   |
+-------+

 

+-----------+
| *n*       |
+-----------+
| ∑         |
+-----------+
| *i*\ =1   |
+-----------+

 

| ⎛
|  ⎝

*x*\ :sub:`*i*`\ −\ *y*\ :sub:`*i*`

| ⎞
|  ⎠

:sup:`2`.

计算时，求和时只考虑*x*\ :sub:`*i*` 和 *y*\ :sub:`*i*` 都存在的值, 分母 *n* 
也相应的做出调整。当分析表达谱数据时，由于 *x*\ :sub:`*i*` 和 *y*\ :sub:`*i*` 
会直接相减, 在使用Euclidean距离前，请对表达谱数据归一化处理.

City-block distance
~~~~~~~~~~~~~~~~~~~

city-block distance也称之为Manhattan 距离，跟Euclidean距离以相关性。Euclidean距离
表示的是两点间最短的距离，而city-block距离是所有维度中距离的和。由于基因表达的数据
经常会有缺失数据，在 ``Bio.Cluster`` 中，city-block距离定义为总距离除以
总维度：

*d* = 

+-------+
| 1     |
+-------+
+-------+
| *n*   |
+-------+

 

+-----------+
| *n*       |
+-----------+
| ∑         |
+-----------+
| *i*\ =1   |
+-----------+

 

| ⎪
|  ⎪

*x*\ :sub:`*i*`\ −\ *y*\ :sub:`*i*`

| ⎪
|  ⎪


这个相当于当你在从城市里一个位置到另一个位置时，所经过街道的距离。
跟Euclidean 距离类似，表达谱的数据会直接相减，因此必须先对数据进行归一化才能
使用。

Pearson 相关系数
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pearson相关系数定义为：

*r* = 

+-------+
| 1     |
+-------+
+-------+
| *n*   |
+-------+

 

+-----------+
| *n*       |
+-----------+
| ∑         |
+-----------+
| *i*\ =1   |
+-----------+

 

| ⎛
|  ⎜
|  ⎜
|  ⎝

+----------------------+
| *x*\ :sub:`*i*` −x   |
+----------------------+
+----------------------+
| σ\ :sub:`*x*`        |
+----------------------+

 

| ⎞
|  ⎟
|  ⎟
|  ⎠

| ⎛
|  ⎜
|  ⎜
|  ⎝

+----------------------+
| *y*\ :sub:`*i*` −ȳ   |
+----------------------+
+----------------------+
| σ\ :sub:`*y*`        |
+----------------------+

 

| ⎞
|  ⎟
|  ⎟
|  ⎠


其中 x, ȳ 分别是 *x* 和 *y* 的样品均值, σ\ :sub:`*x*`, σ\ :sub:`*y*` 
是 *x* 和 *y* 的样品标准差. Pearson相关系数是用于测量 *x* and *y* 散点图对直线的
拟合程度。如果所有的点都在直线上，那么Pearson相关系数为 +1 or -1, 取决于直线的斜率
是正还是负。如果Pearson 相关系数等于0，表明 *x* 和 *y* 之间没有相关性。

*Pearson distance* 定义为

+----------------------------+
| *d*\ :sub:`P` ≡ 1 − *r*.   |
+----------------------------+

由于the Pearson 相关性介于 -1 和 1之间, Pearson 距离的范围为 0 和 2 之间.

Absolute Pearson correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

通过对Pearson相关系数取绝对值，可以得到一个0和1之间的数。如果绝对值是1，所有的点
都位于一条直线上，无论斜率为正还是负。当绝对值为0时，表明 *x* and *y* 没有相关性。

对应的距离定义为：

+------------------------+------+-------+------+-----+
| *d*\ :sub:`A` ≡ 1 −    | ⎪    | *r*   | ⎪    | ,   |
|                        |  ⎪   |       |  ⎪   |     |
+------------------------+------+-------+------+-----+

其中 *r* 是 Pearson 相关系数. 由于Pearson的相关系数介于 0 和 1之间, 对应的
距离也位于0和1之间。

在基因表达数据中，绝对相关性等于1，表明两组基因的表达情况完全一样或者完全
相反，在使用时，应该注意这一点。

Uncentered correlation (cosine of the angle)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

在某些情况下，使用 *uncentered correlation* 比常规的Pearson相关系数更合适。
uncentered correlation 定义为：

*r*\ :sub:`U` = 

+-------+
| 1     |
+-------+
+-------+
| *n*   |
+-------+

 

+-----------+
| *n*       |
+-----------+
| ∑         |
+-----------+
| *i*\ =1   |
+-----------+

 

| ⎛
|  ⎜
|  ⎜
|  ⎝

+-----------------------------+
| *x*\ :sub:`*i*`             |
+-----------------------------+
+-----------------------------+
| σ\ :sub:`*x*`\ :sup:`(0)`   |
+-----------------------------+

 

| ⎞
|  ⎟
|  ⎟
|  ⎠

| ⎛
|  ⎜
|  ⎜
|  ⎝

+-----------------------------+
| *y*\ :sub:`*i*`             |
+-----------------------------+
+-----------------------------+
| σ\ :sub:`*y*`\ :sup:`(0)`   |
+-----------------------------+

 

| ⎞
|  ⎟
|  ⎟
|  ⎠

,

其中

     

σ\ :sub:`*x*`\ :sup:`(0)`

 =

 

√

+-------+
| 1     |
+-------+
+-------+
| *n*   |
+-------+

 

+-----------+
| *n*       |
+-----------+
| ∑         |
+-----------+
| *i*\ =1   |
+-----------+

*x*\ :sub:`*i*`\ :sup:`2`

;  

 

σ\ :sub:`*y*`\ :sup:`(0)`

 =

 

√

+-------+
| 1     |
+-------+
+-------+
| *n*   |
+-------+

 

+-----------+
| *n*       |
+-----------+
| ∑         |
+-----------+
| *i*\ =1   |
+-----------+

*y*\ :sub:`*i*`\ :sup:`2`

.  

 
这个公式同Pearson相关系数的公式一样，只是把样本均值 x, ȳ 设为0 。
uncentered correlation 适用于表达量基准为0的情况。例如，在对基因表达情况计算
比值后取对数，当log-ratio 等于0 表明红色或绿色信号强度相等，也意味着实验处理
不影响基因的表达量。

uncentered correlation 系数对应的距离计算方法为：

+--------------------------------------+
| *d*\ :sub:`U` ≡ 1 − *r*\ :sub:`U`,   |
+--------------------------------------+

其中 *r*\ :sub:`U` 是uncentered 系数。 由于uncentered系数位于-1 和 1
之间，对应的距离范围为 0 与 2之间。

由于 uncentered 系数同 *n* 维空间里的两个数据向量所成角度的cosine值相同，因此
得名。

Absolute uncentered correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

与 Pearson 相关性类似, 也可以用uncentered correlation的绝对值来定义距离:

+-------------------------+------+-----------------+------+-----+
| *d*\ :sub:`AU` ≡ 1 −    | ⎪    | *r*\ :sub:`U`   | ⎪    | ,   |
|                         |  ⎪   |                 |  ⎪   |     |
+-------------------------+------+-----------------+------+-----+

其中 *r*\ :sub:`U` 是 uncentered相关系数。由于uncentered 相关系数的
绝对值位于 0 和 1 之间，对应的距离也为位于 0 和 1之间。

从几何学上来讲，uncentered相关系数的绝对值等于两个数据组成的向量的supporting lines
的角度的cosine值（即不考虑向量的方向性）。

Spearman rank correlation
~~~~~~~~~~~~~~~~~~~~~~~~~

Spearman秩相关系数是一种非参的相关性测量方法，同Pearson相关系数相比，对于离群点
有更好的稳健性。

为了计算Spearman秩相关系数，首先对每个数据集里的数据按值排序，得到每个数据的
秩。然后，计算两个数据集的秩的Pearson相关系数，得到Spearson的相关系数。

同Pearson相关性类似，Spearman秩相关系数对应的距离定义为：

+--------------------------------------+
| *d*\ :sub:`S` ≡ 1 − *r*\ :sub:`S`,   |
+--------------------------------------+

其中 *r*\ :sub:`S` 是Spearman秩相关系数。

Kendall’s τ
~~~~~~~~~~~

Kendall’s τ 是另一个非参的计算相关性的方法。它同Spearman秩相关系数类似，但它使用秩来计算
 τ (see Snedecor & Cochran [`29 <#snedecor1989>`__\ ]) 。

Kendall’s τ 对应的距离计算为：

+--------------------------+
| *d*\ :sub:`K` ≡ 1 − τ.   |
+--------------------------+

因为 Kendall’s τ 位于 -1 和 1之间, 对应的距离位于 0 和 2之间。

Weighting
~~~~~~~~~

对于 ``Bio.Cluster`` 中大部分距离函数，都可以使用权重矩阵。权重矩阵包含着
数据集中每个元素的权重。如果元素 *i* 的权重为 *w*\ :sub:`*i*`，那么这个元素
计算为元素的值乘以 *w*\ :sub:`*i*` 。权重值不需要为整数。对于 Spearman 秩相关系数
和Kendall’s τ, 权重没有很好的定义，因此不能用于这两个函数。

计算距离矩阵
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
距离矩阵是 ``data`` 每个元素的两两间的距离的平方矩阵，可以用 ``Bio.Cluster`` 模块中 ``distancematrix`` 函数计算：
 
.. code:: verbatim

    >>> from Bio.Cluster import distancematrix
    >>> matrix = distancematrix(data)

其中，包含以下参数：

-  ``data`` (必选)
    包含所有元素的矩阵
-  ``mask`` (默认: ``None``)
    显示是否为缺失数据的矩阵。若
   ``mask[i,j]==0``, 那么 ``data[i,j]`` 缺失。若 ``mask==None``,
   那么表明没有缺失数据。
-  ``weight`` (默认: ``None``)
    计算距离时使用的权重矩阵。若
   ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    定义使用的距离函数 (具体见
   `15.1 <#sec:distancefunctions>`__).

为了节省内存，函数运行返回的距离矩阵是一个1D 数组的列表。每一行的列数等于
行号。因此，第一行有0个元素。例如一个返回值为：

.. code:: verbatim

    [array([]),
     array([1.]),
     array([7., 3.]),
     array([4., 2., 6.])]

对应的距离矩阵为：

| ⎛
|  ⎜
|  ⎜
|  ⎜
|  ⎝

+-----+-----+-----+-------+
| 0   | 1   | 7   | 4     |
+-----+-----+-----+-------+
| 1   | 0   | 3   | 2     |
+-----+-----+-----+-------+
| 7   | 3   | 0   | 6     |
+-----+-----+-----+-------+
| 4   | 2   | 6   | 0     |
+-----+-----+-----+-------+

| ⎞
|  ⎟
|  ⎟
|  ⎟
|  ⎠

.

15.2  计算聚类的相关性质
------------------------------------

计算聚类中心
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

聚类中心可以是所有聚类元素的在每个维度上的平均值或者中值，可以用 ``Bio.Cluster`` 中的 ``clustercentroids`` 
函数计算：
 
.. code:: verbatim

    >>> from Bio.Cluster import clustercentroids
    >>> cdata, cmask = clustercentroids(data)

包含了一下参数:

-  ``data`` (必选)
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    用来表示数据是否缺失的整型数组。如果
   ``mask[i,j]==0``, 那么 ``data[i,j]`` 是缺失的. 如果 ``mask==None``,
   那么没有数据缺失.
-  ``clusterid`` (默认: ``None``)
    一个整型向量，用来表示每个元素属于那个类别。如果
   ``clusterid`` 是 ``None``, 表明所有的元素属于相同的类别。
-  ``method`` (默认: ``'a'``)
    指定使用算术平方根 (``method=='a'``) 或者中值
   median (``method=='m'``) 来计算聚类中心。
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).

这个函数返回值为元组 ``(cdata, cmask)``。 聚类中心的数据存储在一个二维的Numerical Python 
数组 ``cdata`` 中, 缺失值的结果存储在二维的Numerical Python整型数组 ``cmask`` 中。 当 ``transpose`` 
为 ``0`` 时，这些数组的长度为（聚类数，列数），当 ``transpose`` 是 ``1`` 时，数组的
长度为 （行数，聚类数）。每一行（当 ``transpose`` = ``0``) 或者 每一列（当 ``transpose`` = ``1`` ）
包含着对应每一聚类中心对应的数据。

计算每类之间的距离
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

根据每个 *items* 的距离函数，我们可以计算出两个 *clusters* 的距离。两个类别的
数学平均值之间的距离通常用于两两间的centroid-linkage聚类和 *k*-means聚类，而 *k*-medoids
聚类中，通常利用两类的中值进行计算。两类中，最短的元素之间的距离用于pairwise single-linkage的聚类，
而最长的元素之间的距离用于计算pairwise maximum-linkage 聚类。在pairwise average-linkage聚类中，
两类之间的距离定义为两两元素间距离的平均值。

为了计算两类之间的距离，可以利用:

.. code:: verbatim

    >>> from Bio.Cluster import clusterdistance
    >>> distance = clusterdistance(data)

其中，包含的参数有：

-  ``data`` (必选)
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    用来表示数据是否缺失的整型数组。如果
   ``mask[i,j]==0``, 那么 ``data[i,j]`` 是缺失的. 如果 ``mask==None``,
   那么没有数据缺失。
-  ``weight`` (默认: ``None``)
    计算距离时使用的权重矩阵。若
   ``weight==None``, 则假设所有的数据使用相同的权重。
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
    选择使用的距离函数 (see
   `15.1 <#sec:distancefunctions>`__).
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).

15.3  Partitioning algorithms
-----------------------------

Partitioning algorithms 依据所有元素到各自聚类中心距离之和最小化原则，
将元素分为 *k* 类。类别的个数 *k* 由用户定义。 ``Bio.Cluster`` 提供了三种不同
的算法:

-  *k*-means 聚类
-  *k*-medians 聚类
-  *k*-medoids 聚类

这些算法的区别在于如何定义聚类中心。在 *k*-means 中, 聚类中心定义为该类中所有
元素的mean data vector。 在 *k*-medians 聚类中， 利用每个维度的中间值来计算。
最后， *k*-medoids 聚类中，聚类中心定义为该类中，离其他所有元素距离之和最小的元素的位置。
这个方法适用于已知距离矩阵，但是原始数据矩阵未知的情况，例如根据结构相似度对蛋白进行聚类
的情况。

expectation-maximization (EM) 算法通常用于将数据分成 *k* 组。在 EM算法的起始阶段,
随机的把元素分配到不同的组。为了保证不存在空元素的类别，可以利用二项分布的方法随机
为每类挑选元素。然后，随机的对分组进行permute，保证每个元素有相同的概率去任何一个类别。
最终，保证每类中至少含有一个元素。

之后进行迭代:

-  利用均值，中值或者medoid计算每类的中心;
-  计算每类的元素离各自中心的距离;
-  对于每个元素，判别其离那个聚类中心最近;
-  对元素重新进行聚类，当不能进行调整时，迭代终止。

为了避免迭代中产生空的类别，在 *k*-means 和 *k*-medians 聚类中，算法始终记录着每类中元素的
个数，并且阻止最后一个元素被分到其他的类别中。对于 *k*-medoids 聚类, 这种检查就是没有必要的，
因为当只剩最后一个元素时，它里中心的距离就为0，所以不会被分配到其他的类别中。

由于起始阶段的元素是随机的，通常当EM算法执行时，会产生不同的聚类结果。为了找到最优的聚类结果，
 *k*-means 算法会重复很多次，每次都以不同的随机分配作为起始。每次运行后，都会保存所有元素距离
 其中心距离之和，并且选择距离最小的那次运行结果最为最终的结果。

EM算法运行的次数取决于需要聚类元素的多少。一般而言，我们可以考虑最优解被发现的次数，
这个次数会作为partitioning算法的返回值。如果最优解被多次返回，那么不太可能存在比这个
更优的解。然后，如果最优解只被发现一次，那么可能存在着距离更小的解。但是，如果需要聚类的
元素过多的话（多余几百），那么很难找到一个全局最优解。

EM算法会在不能进行任何分配的时候停止。我们注意到，对于某些起始的分配，由于
相同的解会在迭代中周期性的重复，从而导致EM算法的失败。因此，我们在迭代中也会
检查这种周期性出现的解。在给定数目的迭代后，当前的聚类结果会保存作为参考。之后
继续迭代一定次数，比较该结果同之前保存的结果，可以确定之前的结果是否重复出现。
如果有重复出现，迭代会终止。如果没有出现，那么再次迭代后的结果会保存作为新的参考。
首先，会重复10次迭代，再保存新的参考。之后，迭代的次数会翻倍，保证在长的周期中也可以
检测到该解。

*k*-means and *k*-medians
~~~~~~~~~~~~~~~~~~~~~~~~~

*k*-means and *k*-medians 算法可以利用 ``Bio.Cluster``中的 ``kcluster`` 实现:

.. code:: verbatim

    >>> from Bio.Cluster import kcluster
    >>> clusterid, error, nfound = kcluster(data)

其中，包含的参数有：

-  ``data`` (必选)
    包含所有元素的矩阵。
-  ``nclusters`` (默认: ``2``)
    聚类的数目 *k*.
-  ``mask`` (默认: ``None``)
    用来表示数据是否缺失的整型数组。如果
   ``mask[i,j]==0``, 那么 ``data[i,j]`` 是缺失的. 如果 ``mask==None``,
   那么没有数据缺失。
-  ``weight`` (默认: ``None``)
    计算距离时使用的权重矩阵。若
   ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).
-  ``npass`` (默认: ``1``)
    *k*-means/-medians 聚类算法运行的次数，每次运行使用不同的随机的起始值。
    如果指定了 ``initialid`` , ``npass`` 的值会忽略，并且聚类算法只会运行一次。
-  ``method`` (默认: ``a``)
    指定确定聚类中心的方法:

   -  ``method=='a'``: 算数平均值 (*k*-means clustering);
   -  ``method=='m'``: 中间值 (*k*-medians clustering).

   当指定 ``method`` 使用其他值时，算法会采用算数平均值。
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择使用的距离函数 (see
   `15.1 <#sec:distancefunctions>`__). 尽管八种距离都可以用于 ``kcluster`` 计算,
   但从经验上来讲，Euclidean 距离适合 *k*-means 算法, city-block 距离适合 *k*-medians.
-  ``initialid`` (默认: ``None``)
    指定EM算法运行初始的聚类类别。如果
   ``initialid==None``, 那么每运行一次EM算法时，都会采取不同的随机初始聚类，总共
   运行的次数由 ``npass`` 决定。如果
   ``initialid`` 不是 ``None``, 那么它应该为一个长度为聚类数的1维数组，每类中至少含有
   一个元素。当初始类别给定后，EM算法的结果也就确定了。

这个函数的返回值为一个包含 ``(clusterid, error, nfound)`` 的元组，其中 ``clusterid`` 是
一个整型矩阵，包含着每类中所包含的行或列的数目。 ``error`` 是最优聚类解中，每类内距离的总和，
``nfound`` 指的是最优解出现的次数。

*k*-medoids 聚类
~~~~~~~~~~~~~~~~~~~~~~

``kmedoids`` 函数根据提供的距离矩阵和聚类数，来运行 *k*-medoids 聚类：

.. code:: verbatim

    >>> from Bio.Cluster import kmedoids
    >>> clusterid, error, nfound = kmedoids(distance)

其中，包含的参数有: , nclusters=2, npass=1,
initialid=None)\|

-  ``distance`` (必选)
    元素两两间的距离矩阵，可以通过三种不同的方法提供：

   -  提供一个2D的 Numerical Python 数组 (函数只会使用矩阵里左下角数据):

      .. code:: verbatim

          distance = array([[0.0, 1.1, 2.3],
                            [1.1, 0.0, 4.5],
                            [2.3, 4.5, 0.0]])

   -  输入一个1D的 Numerical Python 数组，包含了距离矩阵左下角的数据：

      .. code:: verbatim

          distance = array([1.1, 2.3, 4.5])

   -  输入一个列表，包含距离矩阵左下角的数据：

      .. code:: verbatim

          distance = [array([]|,
                      array([1.1]),
                      array([2.3, 4.5])
                     ]

   三种方法对应着同样的距离矩阵。
-  ``nclusters`` (默认: ``2``)
    聚类的数目 *k*.
-  ``npass`` (默认: ``1``)
    *k*-means/-medians 聚类算法运行的次数，每次运行使用不同的随机的起始值。
    如果指定了 ``initialid`` , ``npass`` 的值会忽略，并且聚类算法只会运行一次。
-  ``initialid`` (默认: ``None``)
    指定EM算法运行初始的聚类类别。如果
   ``initialid==None``, 那么每运行一次EM算法时，都会采取不同的随机初始聚类，总共
   运行的次数由 ``npass`` 决定。如果
   ``initialid`` 不是 ``None``, 那么它应该为一个长度为聚类数的1维数组，每类中至少含有
   一个元素。当初始类别给定后，EM算法的结果也就确定了。

函数返回值为一个 包含 ``(clusterid, error, nfound)`` 的元组, 其中
``clusterid`` 一个整型矩阵，包含着每类中所包含的行或列的数目。``error`` 是最优聚类解中，每类内距离的总和，
``nfound`` 指的是最优解出现的次数。需要注意的是， ``clusterid`` 中的聚类数指的是代表聚类中的元素的个数。

15.4  系统聚类
-----------------------------

系统聚类同 *k*-means 聚类有本质的不同。在系统聚类中，基因间或者实验条件间的相似度是通过
树的形式展现出来的。可以利用Treeview或者Java Treeview来查看这些树的结构，因此很便于基因表达
数据中系统聚类的运用。

系统聚类的第一步是计算所有元素间的距离矩阵。之后，融合两个最近的元素成为一个节点。然后，不断的
通过融合相近的元素或者节点来形成新的节点，直到所有的元素都属于同一个节点。在追溯元素和节点融合
的过程的同时形成了树的结构。不同于 *k*-means 使用的EM算法，系统聚类的过程是固定的。

系统聚类也存在着几个不同的方法，他们区别在于如何计算子节点内的距离。在
``Bio.Cluster`` 中，提供了最短距离法（ pairwise single）,最长距离法（maximum）, 类平均法（average）,
和重心法（centroid linkage）。

-  在最短距离法中，节点间的距离被定义两个节点最近样品间距离。
-  在最短距离法中，节点间的距离被定义两个节点最远样品间距离。
-  在类平均法中，节点间的距离被定义为所有样品对之间的平均距离。
-  在重心法中，节点间的距离被定义为两个节点重心间的距离。重心的计算是通过对
每类中所有元素进行计算的。由于每次都要计算新的节点对各个元素和已存在节点的距离，
因此重心法的运行时间比其他系统聚类的方法更长。该方法另外一个特性是，当聚类树的
长大的时候，距离并不会增加，有时候反而减少。由于使用Pearson相关对重心的计算和距离的计算不一致，导致了
这种特性：因为Pearson相关对于计算距离很有效，但是对于重心的计算不是很好normalization。

对于最短距离法，最长距离法和类平均法时，两个节点之间的距离是直接对类别里的元素计算得到的。
因此，聚类的算法在得到距离矩阵后，不一定需要提供最开始的基因表达数据。而对于重心法而言，
新生成的节点的中心必须依靠原始的数据，而不是仅仅依靠距离矩阵。

最短距离法的实现是根据 SLINK algorithm (R. Sibson, 1973), 这个算法具有快速和高效的特点。
并且这个方法聚类的结果同传统的方法结果一致。并且这个算法，可以运用于大量的数据，而传统的
算法则需要大量的内存需求和运行时间。

Representing a hierarchical clustering solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

分层聚类的结果是用树的结构展示所有节点，每个节点包含两个元素或者子节点。通常，我们既关心那个元素
或者哪个子节点互相融合，也关心二者之间的距离（或者相似度）。我们可以调用 ``Bio.Cluster``中的
``Node``函数，来存储聚类树的一个节点。 ``Node``的实例包含以下三个性质：

-  ``left``
-  ``right``
-  ``distance``

其中, ``left`` 和 ``right`` 是两个需要合并的节点所包含的元素或者子节点的个数。
``distance`` 指的是两个节点的距离。需要聚类的元素的编号是从0到（元素数目-1），
而聚类的组别是从-1到-（元素数目-1）。请注意，节点的数目比元素的数目少一。


为了创建一个新的 ``Node`` 对象,我们需要指定 ``left`` 和 ``right``; 
``distance`` 是可选的。

.. code:: verbatim

    >>> from Bio.Cluster import Node
    >>> Node(2,3)
    (2, 3): 0
    >>> Node(2,3,0.91)
    (2, 3): 0.91

一个已存在 ``Node`` 对象的 ``left``, ``right``, 和 ``distance`` 都是可以直接修改的：

.. code:: verbatim

    >>> node = Node(4,5)
    >>> node.left = 6
    >>> node.right = 2
    >>> node.distance = 0.73
    >>> node
    (6, 2): 0.73

当 ``left`` 和 ``right`` 不是整数的时候，或者 ``distance`` 不能被转化成浮点值，会抛出错误。

 Python的类 ``Tree`` 包含着整个系统聚类的结果。 ``Tree`` 的对象可以通过
 一个 ``Node`` 的列表创建:

.. code:: verbatim

    >>> from Bio.Cluster import Node, Tree
    >>> nodes = [Node(1,2,0.2), Node(0,3,0.5), Node(-2,4,0.6), Node(-1,-3,0.9)]
    >>> tree = Tree(nodes)
    >>> print tree
    (1, 2): 0.2
    (0, 3): 0.5
    (-2, 4): 0.6
    (-1, -3): 0.9

 ``Tree`` 的初始器会检查包含节点的list是否是一个有效的系统聚类树的结果:

.. code:: verbatim

    >>> nodes = [Node(1,2,0.2), Node(0,2,0.5)]
    >>> Tree(nodes)
    Traceback (most recent call last):
      File "<stdin>", line 1, in ?
    ValueError: Inconsistent tree

也可以用中括号来对 ``Tree`` 对象进行检索：

.. code:: verbatim

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

.. code:: verbatim

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

这保证了``Tree`` 都具有良好的结果。

为了利用可视化工具，例如Java Treeview，来查看系统聚类树，最好对所有节点的距离进行归一化，
使其位于0和1之间。可以通过对 ``Tree`` 对象调用 ``scale`` 方法来实现这个功能：

.. code:: verbatim

    >>> tree.scale()

这个方法不需要任何参数，返回值是 ``None``.

经过系统聚类后，可以对 ``Tree`` 对象进行剪接，将所有的元素分为 *k* 类：

.. code:: verbatim

    >>> clusterid = tree.cut(nclusters=1)

其中 ``nclusters`` (默认是 ``1``) 是期望的类别数 *k*。这个方法会忽略树结构里面的
最高的 *k*\ −1 节点，最终形成 *k* 个独立的类别。对于 *k* 必须为正数，并且小于或者等于
元素的数目。这个方法会返回一个数组 ``clusterid`` ,包含着每类中所包含元素的个数。

运行系统聚类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了进行系统聚类，可以用 ``Bio.Cluster`` 中的 ``treecluster`` 函数。

.. code:: verbatim

    >>> from Bio.Cluster import treecluster
    >>> tree = treecluster(data)

包含着以下参数:

-  ``data``
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    用来表示数据是否缺失的整型数组。如果
   ``mask[i,j]==0``, 那么 ``data[i,j]`` 是缺失的. 如果 ``mask==None``,
   那么没有数据缺失。
-  ``weight`` (默认: ``None``)
    计算距离时使用的权重矩阵。若
   ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).
-  ``method`` (默认: ``'m'``)
    选择合适的节点距离计算方法:

   -  ``method=='s'``: 最小距离法
   -  ``method=='m'``: 最大距离法
   -  ``method=='c'``: 重心法
   -  ``method=='a'``: 类平均法

-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择合适的元素距离算法 (see
   `15.1 <#sec:distancefunctions>`__).

为了对一个事先计算好的距离矩阵进行系统聚类，可以在调用 ``treecluster`` 时候，
用 ``distancematrix`` 参数来代替 ``data`` 参数：

.. code:: verbatim

    >>> from Bio.Cluster import treecluster
    >>> tree = treecluster(distancematrix=distance)

这种情况下，需要定义下面的参数：

-  ``distancematrix``
    元素两两间的距离矩阵，可以通过三种不同的方法提供：

   -  提供一个2D的 Numerical Python 数组 (函数只会使用矩阵里左下角数据):

      .. code:: verbatim

          distance = array([[0.0, 1.1, 2.3], 
                            [1.1, 0.0, 4.5],
                            [2.3, 4.5, 0.0]])

   -  输入一个1D的 Numerical Python 数组，包含了距离矩阵左下角的数据：

      .. code:: verbatim

          distance = array([1.1, 2.3, 4.5])

   -  输入一个列表，包含距离矩阵左下角的数据：

      .. code:: verbatim

          distance = [array([]),
                      array([1.1]),
                      array([2.3, 4.5])

      三种方法对应着同样的距离矩阵。由于 ``treecluster`` 会对距离矩阵中的值进行随机洗牌，
      如果后面需要调用这个距离矩阵，请在调用 ``treecluster`` 之情，事先存到一个新的变量

-  ``method``
    选择合适的节点距离计算方法:

   -  ``method=='s'``: 最小距离法
   -  ``method=='m'``: 最大距离法
   -  ``method=='a'``: 类平均法

   其中，最小距离法、最大距离法和类平均法可以只通过距离矩阵计算，而重心法却不行。

当调用 ``treecluster``时,  ``data`` 或者 ``distancematrix`` 总有一个必须为 ``None``。

函数返回一个 ``Tree`` 对象，该对象包含着 (元素数目-1）个节点，当选择行作为聚类时，元素的
数目同行数一致；当使用列作为聚类时，元素的数目同列数一致。每个节点都意味着一对相邻连锁的
事件，其中节点的性质 ``left`` 和 ``right`` 包含着每个合并的元素或者子节点的元素数， ``distance`` 
是两个合并元素或者子节点的距离。元素是从 0 to (元素数目 − 1) 进行标记, 而类别是从 -1 到 −(元素
数目 -1 ）

15.5  Self-Organizing Maps
--------------------------

Self-Organizing Maps (SOMs) 是由 Kohonen 在描述神经网络的时候发明的 (see for instance Kohonen, 1997 [`24 <#kohonen1997>`__\ ]).
Tamayo (1999) 第一次讲 Self-Organizing Maps 应用到基因表达数据上。
[`30 <#tamayo1999>`__\ ].

SOMs 根据某种拓扑结果将元素进行分类。通常选用的是矩形的拓扑结构。在SOMs生成的类别中，相邻的
两个类的拓扑结构相似度高于他们对其他的相似度。

计算SOM的第一步是随机分配数据向量到每个类别中，如果使用行进行聚类，那么每个数据向量中的元素
个数等于列数。

一个SOM 会一次读入整行，并且找到该向量最近的拓扑聚类结构。通过对整行里的数据向量对
这个类别的数据向量，同相邻的类别的数据向量进行调整。调整如下：

Δ 

+-------+
| *x*   |
+-------+
+-------+

:sub:`cell` = τ · 

| ⎛
|  ⎜
|  ⎝

+-------+
| *x*   |
+-------+
+-------+

:sub:`row` − 

+-------+
| *x*   |
+-------+
+-------+

:sub:`cell` 

| ⎞
|  ⎟
|  ⎠

.

参数 τ 会随着迭代次数增加而减少。可以用一个简单的线性函数来定义其与迭代次数的关系：

τ = τ\ :sub:`init` · 

| ⎛
|  ⎜
|  ⎜
|  ⎝

1 − 

+--------+
| *i*    |
+--------+
+--------+
| *n*    |
+--------+

| ⎞
|  ⎟
|  ⎟
|  ⎠

,

τ\ :sub:`init` 是可以指定的起始的 τ 值， *i* 是当前迭代的次数， *n* 是总的需要迭代的次数。
在迭代开始时，τ变化很快，然而在迭代末尾，变化越来越小。

所有在半径 *R* 内的类别都会在每次迭代中进行调整。半径也会随着迭代的增加而减小：

*R* = *R*\ :sub:`max` · 

| ⎛
|  ⎜
|  ⎜
|  ⎝

1 − 

+--------+
| *i*    |
+--------+
+--------+
| *n*    |
+--------+

| ⎞
|  ⎟
|  ⎟
|  ⎠

,

其中最大的半径定义为：

*R*\ :sub:`max` = 

√

+---------------------------------------------------------+
+---------------------------------------------------------+
| *N*\ :sub:`*x*`\ :sup:`2` + *N*\ :sub:`*y*`\ :sup:`2`   |
+---------------------------------------------------------+

,

其中 (*N*\ :sub:`*x*`, *N*\ :sub:`*y*`) 是定义拓扑结构的矩形维度。

函数 ``somcluster`` 可以用来在一个矩形的网格里计算 Self-Organizing Map。
首先，初始化一个随机数产生器。利用随机化产生器来对节点数据进行初始化。
基因或者芯片的使用顺序同样是随机的。用户可以定义总的SOM迭代的次数。

运行 ``somcluster``, 例如：

.. code:: verbatim

    >>> from Bio.Cluster import somcluster
    >>> clusterid, celldata = somcluster(data)

其中，可以定义一下参数:

-  ``data`` (required)
    包含所有元素的矩阵。
-  ``mask`` (默认: ``None``)
    用来表示数据是否缺失的整型数组。如果
   ``mask[i,j]==0``, 那么 ``data[i,j]`` 是缺失的. 如果 ``mask==None``,
   那么没有数据缺失。
-  ``weight`` (默认: ``None``)
    计算距离时使用的权重矩阵。若
   ``weight==None``, 则假设所有的数据使用相同的权重。
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).
-  ``nxgrid, nygrid`` (默认: ``2, 1``)
    当Self-Organizing Map计算的时候，矩形的网格所包含的横向和纵向的格子。
-  ``inittau`` (默认: ``0.02``)
    SOM算法中，参数 τ 的初始值，默认是 0.02。 这个初始值同Michael Eisen’s Cluster/TreeView 一致。
-  ``niter`` (默认: ``1``)
    迭代运行的次数。
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择合适的元素距离算法 (见
   `15.1 <#sec:distancefunctions>`__).

这个函数返回的是一个元组 ``(clusterid, celldata)``:

-  ``clusterid``:
    一个两列的数组，行的数目等于待聚类元素的个数。每行包含着在矩形SOM网格中，将每个元素分配到的
    格子的 *x* 和 *y* 的坐标。
-  ``celldata``:
    当以行进行聚类时，生成的矩阵维度为 (``nxgrid``, ``nygrid``, number of columns)；
   当以列进行聚类时，生成的矩阵维度为 (``nxgrid``, ``nygrid``, number of  rows)。
   在这个矩阵里， ``[ix][iy]`` 表示着一个1D向量，其中用于计算该类中心的这基因的表达谱数据.

15.6  主成分分析
----------------------------------

主成分分析 (PCA) 被广泛的用于分析多维数据，一个将主成分分析应用于表达谱数据的请见
Yeung and Ruzzo (2001) [`33 <#yeung2001>`__\ ].

简而言之，PCA是一种坐标转换的方法，转换后的基础向量成为主成分，变换前的每行可以用主成分利用
线性关系显示。主成分的选择是基于是残差尽可能的小。例如，一个 *n* × 3 的数据矩阵可以表示为三维
空间内的一个椭圆球形的点的云。第一主成分是这个椭圆球形的最长轴，第二主成分是次长轴，第三主成分
是最短的轴。矩阵中，每一行都可以用主成分的线性关系展示。一般而言，为了对数据进行降维，只保留最
重要的几个主成分。剩余的残差认为是不可解释的方差。

可以通过计算数据的协方差矩阵的特征向量来得到主成分。每个主成分对应的特征值决定了
其在数据中代表的方差的大小。

在进行主成分分析前，矩阵的数据每一列都要减去其平均值。在上面的例子中，椭圆球形云在3D
空间中，围绕着其中心分布，而主成分则显示着每个点对其中心的变化。

函数 ``pca`` 首先使用奇异值分解（singular value decomposition）来计算矩阵的特征值和
特征向量。奇异值分解使用的是Algol写的C语言的 ``svd`` [`16 <#golub1971>`__\ ], 利用的是
Householder bidiagonalization 和 QR 算法的变异。主成分，每个数据在主成分上的坐标和主成分
对应的特征值都会被计算出来，并按照特征值的降序排列。如果需要数据中心，则需要在调用 ``pca`` 
前，对每列数据减去其平均值。

将主成分分析应用于二维矩阵 ``data``,可以：

.. code:: verbatim

    >>> from Bio.Cluster import pca
    >>> columnmean, coordinates, components, eigenvalues = pca(data)

函数会返回一个元组：
``columnmean, coordinates, components, eigenvalues``:

-  ``columnmean``
    包含 ``data`` 每列均值的数组 .
-  ``coordinates``
     数据 ``data`` 中每行在主成分上对应的坐标。
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
```http://jtreeview.sourceforge.net/`` <http://jtreeview.sourceforge.net/>`__\ Java
TreeView 程序。这个软件可以显示系统聚类和 *k*-means 聚类的结果。

类 ``Record`` 的一个对象包含着一个 Cluster/TreeView-type数据文件需要的所有信息。
为了将结果保存到一个 ``Record`` 对象中，首先需要打开一个文件，并读取：

.. code:: verbatim

    >>> from Bio import Cluster
    >>> handle = open("mydatafile.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()

分成两步操作可以对数据的来源有更好的可变性，例如可以利用：

.. code:: verbatim

    >>> import gzip # Python standard library
    >>> handle = gzip.open("mydatafile.txt.gz")

来打开一个gzipped文件，或者利用

.. code:: verbatim

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
    用来表示数据是否缺失的整型数组。如果
   ``mask[i,j]==0``, 那么 ``data[i,j]`` 是缺失的. 如果 ``mask==None``,
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

.. code:: verbatim

    >>> matrix = record.distancematrix()

其中，包含以下参数：

-  ``transpose`` (默认: ``0``)
       选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择合适的元素距离算法 (见
   `15.1 <#sec:distancefunctions>`__).

函数会返回一个距离矩阵，每行的列数等于行数。(see section
`15.1 <#subsec:distancematrix>`__).

计算聚类中心
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了计算存储在record中的元素的聚类中心，利用：

.. code:: verbatim

    >>> cdata, cmask = record.clustercentroids()

-  ``clusterid`` (默认: ``None``)
    用于显示每个元素属于哪类的整型的向量。如果缺少 ``clusterid``,
    那么所有的元素属于同一类。
-  ``method`` (默认: ``'a'``)
    选择是否使用算术平均值 (``method=='a'``) 或者中值 (``method=='m'``) 
    来计算聚类中心。
-  ``transpose`` (默认: ``0``)
    选择计算``data`` 的行的中心 (``transpose==0``), 或者计算列的中心 (``transpose==1``).

函数返回元组 ``cdata, cmask``; 见 section
`15.2 <#subsec:clustercentroids>`__ for a description.

计算两类间的距离
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了计算存储在record中，两类的距离，利用：

.. code:: verbatim

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
    选择使用的距离函数 (见
   `15.1 <#sec:distancefunctions>`__).
-  ``transpose`` (默认: ``0``)
        选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``)..

进行系统聚类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了对存储在record中的数据进行系统聚类，利用：

.. code:: verbatim

    >>> tree = record.treecluster()

包含以下参数:

-  ``transpose`` (默认: ``0``)
    选择使用行 (``transpose==0``) 或者列 (``transpose==1``) 用于聚类
-  ``method`` (默认: ``'m'``)
    选择合适的节点距离计算方法:

   -  ``method=='s'``: 最小距离法
   -  ``method=='m'``: 最大距离法
   -  ``method=='c'``: 重心法
   -  ``method=='a'``: 类平均法

-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择使用的距离函数(见
   `15.1 <#sec:distancefunctions>`__).
-  ``transpose``
    选择使用基因或者芯片进行聚类，如果是 ``transpose==0``,则使用基因 (行) 进行聚类，如果使用
   ``transpose==1``, 芯片 (列) 用于聚类.

函数返回 ``Tree`` 对象。对象包含 (元素数目 − 1） 节点, 如果使用行进行聚类时，元素数目为总行数；
当使用列进行聚类时，元素数目为总列数。每个节点描述着一对节点连接，然而节点的性质 ``left`` 和
``right`` 包含着相邻节点所有的元素和子节点数， ``distance`` 显示着左右节点的距离。
元素从 0 到 (元素数目 − 1) 进行索引, 而类别从 -1 to −(元素数目−1)进行索引。

进行 *k*-means or *k*-medians 聚类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了对存储在record中的元素进行 *k*-means 或者 *k*-medians 聚类，可以使用：

.. code:: verbatim

    >>> clusterid, error, nfound = record.kcluster()

包含以下参数:

-  ``nclusters`` (默认: ``2``)
    聚类的数目 *k*.
-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).
-  ``npass`` (默认: ``1``)
    *k*-means/-medians 聚类算法运行的次数，每次运行使用不同的随机的起始值。
    如果指定了 ``initialid`` , ``npass`` 的值会忽略，并且聚类算法只会运行一次。
-  ``method`` (默认: ``a``)
    指定确定聚类中心的方法:

   -  ``method=='a'``: 算数平均值 (*k*-means clustering);
   -  ``method=='m'``: 中间值 (*k*-medians clustering).

   当指定 ``method`` 使用其他值时，算法会采用算数平均值。
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择使用的距离函数 (see
   `15.1 <#sec:distancefunctions>`__).

这个函数返回的是一个元组 ``(clusterid, error, nfound)``, 其中 ``clusterid`` 是一个包含
该类中所有的行或者类别的数目（还是id）， ``error`` 是最优解的类内的距离和， ``nfound`` 
是最优解被发现的次数。

计算Self-Organizing Map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

可以利用以下命令，计算对存储在record中元素计算 Self-Organizing Map ：

.. code:: verbatim

    >>> clusterid, celldata = record.somcluster()

包含以下参数:

-  ``transpose`` (默认: ``0``)
    选择 使用 ``data`` 的行行之间计算距离 (``transpose==0``), 或者列与列计算距离 (``transpose==1``).
-  ``nxgrid, nygrid`` (默认: ``2, 1``)
    当Self-Organizing Map计算时，在矩形网格里的横向和纵向格子数目
-  ``inittau`` (默认: ``0.02``)
    用于SOM算法的参数 τ 的初始值。默认的 ``inittau`` 是0.02，同Michael Eisen’s Cluster/TreeView 程序中
    使用的参数一致。
-  ``niter`` (默认: ``1``)
    迭代运行的次数。
-  ``dist`` (默认: ``'e'``, Euclidean distance)
    选择使用的距离函数(see
   `15.1 <#sec:distancefunctions>`__).

函数返回一个元组 ``(clusterid, celldata)``:

-  ``clusterid``:
    一个二维数组，行数同待聚类的元素数目相同。每行的内容对应着该元素在矩形SOM方格内 *x* 和
   *y* 的坐标。
-  ``celldata``:
    格式为一个矩阵，如果是对行聚类，内容为 (``nxgrid``, ``nygrid``, 列数)，如果是对列聚类，
    那么内容为 (``nxgrid``, ``nygrid``, 行数) 。矩阵中，坐标 ``[ix][iy]`` 对应的是该坐标的网格里的
    基因表达数据的聚类中心的1D向量。

保存聚类结果
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

为了保存聚类结果，可以利用：

.. code:: verbatim

    >>> record.save(jobname, geneclusters, expclusters)

包含以下参数:

-  ``jobname``
    字符串 ``jobname`` 作为保存的文件名。
-  ``geneclusters``
    这个参数指的是基因（以行聚类）的结果。在 *k*-means 聚类中，这个参数是一个1D的数组，包含着
    每个基因对应的类别，可以通过 ``kcluster`` 得到。在系统聚类中， ``geneclusters`` 是一个``Tree`` 对象。
-  ``expclusters``
    这个参数指的是实验条件（以列聚类）的结果。在 *k*-means 聚类中，这个参数是一个1D的数组，包含着
    每个实验条件对应的类别，可以通过 ``kcluster`` 得到。在系统聚类中， ``geneclusters`` 是一个``Tree`` 对象。

这个方法会生成文本文件 ``jobname.cdt``, ``jobname.gtr``, ``jobname.atr``, ``jobname*.kgg``, 
和/或 ``jobname*.kag``。 这些文件可以用于后续分析。如果 ``geneclusters`` 和 ``expclusters`` 
都是 ``None``, 那这个方法只会生成 ``jobname.cdt``; 这个文件可以被读取，生成一个新的 ``Record`` 对象.

15.8  示例
-------------------------

以下是一个系统聚类的例子，其中使用最短距离法对基因进行聚类，用最大距离法对实验条件进行聚类。
由于使用 Euclidean 距离对基因进行聚类，因此需要将节点距离 ``genetree`` 进行调整，使其处于0和1之间。
这种调整对于Java TreeView正确显示树结构也是很必须的。为了对实验条件进行聚类，使用了 非中心的关联分析（uncentered 
correlation）。在这种情况下，不需要任何的调整，因为 ``exptree`` 中的结果已经位于0和2之间。 示例中使用的
文件 ``cyano.txt`` 可以再 ``data`` 文件夹中找到。

.. code:: verbatim

    >>> from Bio import Cluster
    >>> handle = open("cyano.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()
    >>> genetree = record.treecluster(method='s')
    >>> genetree.scale()
    >>> exptree = record.treecluster(dist='u', transpose=1)
    >>> record.save("cyano_result", genetree, exptree)

这个命令会生成 ``cyano_result.cdt``, ``cyano_result.gtr``, 和 ``cyano_result.atr``等文件。

同样的，也可以保存一个 *k*-means 聚类的结果:

.. code:: verbatim

    >>> from Bio import Cluster
    >>> handle = open("cyano.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()
    >>> (geneclusters, error, ifound) = record.kcluster(nclusters=5, npass=1000)
    >>> (expclusters, error, ifound) = record.kcluster(nclusters=2, npass=100, transpose=1)
    >>> record.save("cyano_result", geneclusters, expclusters)

这个会生成文件 ``cyano_result_K_G2_A2.cdt``,``cyano_result_K_G2.kgg``, 和 ``cyano_result_K_A2.kag``.

15.9  其他函数
-------------------------

``median(data)`` 返回1D数组 ``data``的中值

``mean(data)`` 返回1D数组``data``的均值。

``version()`` 返回使用的C聚类库的版本号。

