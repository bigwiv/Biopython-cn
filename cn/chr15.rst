第十五章 聚类分析
============================

聚类分析是对一堆元素按照相似度来进行分组的过程。在生物信息中，聚类分析经常
用于分析基因表达的数据，用来发现有相同的表达图谱的基因，以鉴定功能相关的基
因，或者发现未知的功能相关基因。

Biopython中的 ``Bio.Cluster`` 模块提供了常用的聚类算法，并且设计的时候考虑了
基因表达数据的应用。然而，这个模块也可以用于其他类型数据。 ``Bio.Cluster`` 
和其使用的C聚类库的说明见De Hoon *et al.* [`14 <#dehoon2004>`__\ ].

The following four clustering approaches are implemented in
``Bio.Cluster``:
``Bio.Cluster``包含了以下四种聚类算法：

-  Hierarchical clustering (pairwise centroid-, single-, complete-, and
   average-linkage);
-  *k*-means, *k*-medians, and *k*-medoids clustering;
-  Self-Organizing Maps;
-  Principal Component Analysis.

-  分层聚类（两两聚类，single-，complete-, and average-linkage);
-  *k*-means, *k*-medians, and *k*-medoids clustering;
-  Self-Organizing Maps;
-  主成分分析

Data representation
数据结构
~~~~~~~~~~~~~~~~~~~

The data to be clustered are represented by a *n* × *m* Numerical Python
array ``data``. Within the context of gene expression data clustering,
typically the rows correspond to different genes whereas the columns
correspond to different experimental conditions. The clustering
algorithms in ``Bio.Cluster`` can be applied both to rows (genes) and to
columns (experiments).
用于聚类的数据为一个*n* x *m*的Numerical Python 矩阵 ``data``。
其中，每一行表示不同的基因，每一列表示不同的实验条件。``Bio.Cluster``既可以
针对每行（基因），也可以针对每列（实验条件）进行聚类。

Missing values
缺失值
~~~~~~~~~~~~~~

Often in microarray experiments, some of the data values are missing,
which is indicated by an additional *n* × *m* Numerical Python integer
array ``mask``. If ``mask[i,j]==0``, then ``data[i,j]`` is missing and
is ignored in the analysis.
在芯片实验中，经常会有些缺失值，通常用一个额外的 *n* × *m* Numerical Python
integer矩阵``mask``表示。例如``mask[i,j]==0``，表示``data[i,j]``是个缺失值，
并且在分析中忽略

Random number generator
随机数生成器
~~~~~~~~~~~~~~~~~~~~~~~
*k*-means/medians/medoids clustering algorithms and Self-Organizing
Maps (SOMs) 需要调用随机数生成器。``Bio.Cluster``中使用的正态分布随机数
生成器使用的算法是基于L’Ecuyer [`25 <#lecuyer1988>`__\ ],二项分布的随机数
生成是根据Kachitvichyanukul and Schmeiser [`23 <#kachitvichyanukul1988>`__\ ]
开发的BTPE算法。随机数生成器在调用时首先进行初始化。由于随机数生成器使用了
两个multiplicative linear congruential generators，所以初始化时需要两个整型的
种子，并调用系统提供的``rand``（C标注库）生成随机数。在``Bio.Cluster``中，
也可以调用``srand``使用当前时间的秒作为初始值，用``rand``随机产生的头两
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

and are therefore refered to as *metrics*. In everyday language, this
means that the shortest distance between two points is a straight line.

The remaining six distance measures are related to the correlation
coefficient, where the distance *d* is defined in terms of the
correlation *r* by *d*\ =1−\ *r*. Note that these distance functions are
*semi-metrics* that do not satisfy the triangle inequality. For example,
for

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

we find a Pearson distance *d*\ (*u*,\ *w*) = 1.8660, while
*d*\ (*u*,\ *v*)+\ *d*\ (*v*,\ *w*) = 1.6340.

Euclidean distance
~~~~~~~~~~~~~~~~~~

In ``Bio.Cluster``, we define the Euclidean distance as

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

Only those terms are included in the summation for which both
*x*\ :sub:`*i*` and *y*\ :sub:`*i*` are present, and the denominator *n*
is chosen accordingly. As the expression data *x*\ :sub:`*i*` and
*y*\ :sub:`*i*` are subtracted directly from each other, we should make
sure that the expression data are properly normalized when using the
Euclidean distance.

City-block distance
~~~~~~~~~~~~~~~~~~~

The city-block distance, alternatively known as the Manhattan distance,
is related to the Euclidean distance. Whereas the Euclidean distance
corresponds to the length of the shortest path between two points, the
city-block distance is the sum of distances along each dimension. As
gene expression data tend to have missing values, in ``Bio.Cluster`` we
define the city-block distance as the sum of distances divided by the
number of dimensions:

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

.

This is equal to the distance you would have to walk between two points
in a city, where you have to walk along city blocks. As for the
Euclidean distance, the expression data are subtracted directly from
each other, and we should therefore make sure that they are properly
normalized.

The Pearson correlation coefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Pearson correlation coefficient is defined as

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

,

in which x, ȳ are the sample mean of *x* and *y* respectively, and
σ\ :sub:`*x*`, σ\ :sub:`*y*` are the sample standard deviation of *x*
and *y*. The Pearson correlation coefficient is a measure for how well a
straight line can be fitted to a scatterplot of *x* and *y*. If all the
points in the scatterplot lie on a straight line, the Pearson
correlation coefficient is either +1 or -1, depending on whether the
slope of line is positive or negative. If the Pearson correlation
coefficient is equal to zero, there is no correlation between *x* and
*y*.

The *Pearson distance* is then defined as

+----------------------------+
| *d*\ :sub:`P` ≡ 1 − *r*.   |
+----------------------------+

As the Pearson correlation coefficient lies between -1 and 1, the
Pearson distance lies between 0 and 2.

Absolute Pearson correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By taking the absolute value of the Pearson correlation, we find a
number between 0 and 1. If the absolute value is 1, all the points in
the scatter plot lie on a straight line with either a positive or a
negative slope. If the absolute value is equal to zero, there is no
correlation between *x* and *y*.

The corresponding distance is defined as

+------------------------+------+-------+------+-----+
| *d*\ :sub:`A` ≡ 1 −    | ⎪    | *r*   | ⎪    | ,   |
|                        |  ⎪   |       |  ⎪   |     |
+------------------------+------+-------+------+-----+

where *r* is the Pearson correlation coefficient. As the absolute value
of the Pearson correlation coefficient lies between 0 and 1, the
corresponding distance lies between 0 and 1 as well.

In the context of gene expression experiments, the absolute correlation
is equal to 1 if the gene expression profiles of two genes are either
exactly the same or exactly opposite. The absolute correlation
coefficient should therefore be used with care.

Uncentered correlation (cosine of the angle)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, it may be preferable to use the *uncentered correlation*
instead of the regular Pearson correlation coefficient. The uncentered
correlation is defined as

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

where

     

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

 

This is the same expression as for the regular Pearson correlation
coefficient, except that the sample means x, ȳ are set equal to zero.
The uncentered correlation may be appropriate if there is a zero
reference state. For instance, in the case of gene expression data given
in terms of log-ratios, a log-ratio equal to zero corresponds to the
green and red signal being equal, which means that the experimental
manipulation did not affect the gene expression.

The distance corresponding to the uncentered correlation coefficient is
defined as

+--------------------------------------+
| *d*\ :sub:`U` ≡ 1 − *r*\ :sub:`U`,   |
+--------------------------------------+

where *r*\ :sub:`U` is the uncentered correlation. As the uncentered
correlation coefficient lies between -1 and 1, the corresponding
distance lies between 0 and 2.

The uncentered correlation is equal to the cosine of the angle of the
two data vectors in *n*-dimensional space, and is often referred to as
such.

Absolute uncentered correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As for the regular Pearson correlation, we can define a distance measure
using the absolute value of the uncentered correlation:

+-------------------------+------+-----------------+------+-----+
| *d*\ :sub:`AU` ≡ 1 −    | ⎪    | *r*\ :sub:`U`   | ⎪    | ,   |
|                         |  ⎪   |                 |  ⎪   |     |
+-------------------------+------+-----------------+------+-----+

where *r*\ :sub:`U` is the uncentered correlation coefficient. As the
absolute value of the uncentered correlation coefficient lies between 0
and 1, the corresponding distance lies between 0 and 1 as well.

Geometrically, the absolute value of the uncentered correlation is equal
to the cosine between the supporting lines of the two data vectors
(i.e., the angle without taking the direction of the vectors into
consideration).

Spearman rank correlation
~~~~~~~~~~~~~~~~~~~~~~~~~

The Spearman rank correlation is an example of a non-parametric
similarity measure, and tends to be more robust against outliers than
the Pearson correlation.

To calculate the Spearman rank correlation, we replace each data value
by their rank if we would order the data in each vector by their value.
We then calculate the Pearson correlation between the two rank vectors
instead of the data vectors.

As in the case of the Pearson correlation, we can define a distance
measure corresponding to the Spearman rank correlation as

+--------------------------------------+
| *d*\ :sub:`S` ≡ 1 − *r*\ :sub:`S`,   |
+--------------------------------------+

where *r*\ :sub:`S` is the Spearman rank correlation.

Kendall’s τ
~~~~~~~~~~~

Kendall’s τ is another example of a non-parametric similarity measure.
It is similar to the Spearman rank correlation, but instead of the ranks
themselves only the relative ranks are used to calculate τ (see Snedecor
& Cochran [`29 <#snedecor1989>`__\ ]).

We can define a distance measure corresponding to Kendall’s τ as

+--------------------------+
| *d*\ :sub:`K` ≡ 1 − τ.   |
+--------------------------+

As Kendall’s τ is always between -1 and 1, the corresponding distance
will be between 0 and 2.

Weighting
~~~~~~~~~

For most of the distance functions available in ``Bio.Cluster``, a
weight vector can be applied. The weight vector contains weights for the
items in the data vector. If the weight for item *i* is *w*\ :sub:`*i*`,
then that item is treated as if it occurred *w*\ :sub:`*i*` times in the
data. The weight do not have to be integers. For the Spearman rank
correlation and Kendall’s τ, weights do not have a well-defined meaning
and are therefore not implemented.

Calculating the distance matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The distance matrix is a square matrix with all pairwise distances
between the items in ``data``, and can be calculated by the function
``distancematrix`` in the ``Bio.Cluster`` module:

.. code:: verbatim

    >>> from Bio.Cluster import distancematrix
    >>> matrix = distancematrix(data)

where the following arguments are defined:

-  ``data`` (required)
    Array containing the data for the items.
-  ``mask`` (default: ``None``)
    Array of integers showing which data are missing. If
   ``mask[i,j]==0``, then ``data[i,j]`` is missing. If ``mask==None``,
   then all data are present.
-  ``weight`` (default: ``None``)
    The weights to be used when calculating distances. If
   ``weight==None``, then equal weights are assumed.
-  ``transpose`` (default: ``0``)
    Determines if the distances between the rows of ``data`` are to be
   calculated (``transpose==0``), or between the columns of ``data``
   (``transpose==1``).
-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).

To save memory, the distance matrix is returned as a list of 1D arrays.
The number of columns in each row is equal to the row number. Hence, the
first row has zero elements. An example of the return value is

.. code:: verbatim

    [array([]),
     array([1.]),
     array([7., 3.]),
     array([4., 2., 6.])]

This corresponds to the distance matrix

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

15.2  Calculating cluster properties
------------------------------------

Calculating the cluster centroids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The centroid of a cluster can be defined either as the mean or as the
median of each dimension over all cluster items. The function
``clustercentroids`` in ``Bio.Cluster`` can be used to calculate either:

.. code:: verbatim

    >>> from Bio.Cluster import clustercentroids
    >>> cdata, cmask = clustercentroids(data)

where the following arguments are defined:

-  ``data`` (required)
    Array containing the data for the items.
-  ``mask`` (default: ``None``)
    Array of integers showing which data are missing. If
   ``mask[i,j]==0``, then ``data[i,j]`` is missing. If ``mask==None``,
   then all data are present.
-  ``clusterid`` (default: ``None``)
    Vector of integers showing to which cluster each item belongs. If
   ``clusterid`` is ``None``, then all items are assumed to belong to
   the same cluster.
-  ``method`` (default: ``'a'``)
    Specifies whether the arithmetic mean (``method=='a'``) or the
   median (``method=='m'``) is used to calculate the cluster center.
-  ``transpose`` (default: ``0``)
    Determines if the centroids of the rows of ``data`` are to be
   calculated (``transpose==0``), or the centroids of the columns of
   ``data`` (``transpose==1``).

This function returns the tuple ``(cdata, cmask)``. The centroid data
are stored in the 2D Numerical Python array ``cdata``, with missing data
indicated by the 2D Numerical Python integer array ``cmask``. The
dimensions of these arrays are (number of clusters, number of columns)
if ``transpose`` is ``0``, or (number of rows, number of clusters) if
``transpose`` is ``1``. Each row (if ``transpose`` is ``0``) or column
(if ``transpose`` is ``1``) contains the averaged data corresponding to
the centroid of each cluster.

Calculating the distance between clusters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a distance function between *items*, we can define the distance
between two *clusters* in several ways. The distance between the
arithmetic means of the two clusters is used in pairwise
centroid-linkage clustering and in *k*-means clustering. In *k*-medoids
clustering, the distance between the medians of the two clusters is used
instead. The shortest pairwise distance between items of the two
clusters is used in pairwise single-linkage clustering, while the
longest pairwise distance is used in pairwise maximum-linkage
clustering. In pairwise average-linkage clustering, the distance between
two clusters is defined as the average over the pairwise distances.

To calculate the distance between two clusters, use

.. code:: verbatim

    >>> from Bio.Cluster import clusterdistance
    >>> distance = clusterdistance(data)

where the following arguments are defined:

-  ``data`` (required)
    Array containing the data for the items.
-  ``mask`` (default: ``None``)
    Array of integers showing which data are missing. If
   ``mask[i,j]==0``, then ``data[i,j]`` is missing. If ``mask==None``,
   then all data are present.
-  ``weight`` (default: ``None``)
    The weights to be used when calculating distances. If
   ``weight==None``, then equal weights are assumed.
-  ``index1`` (default: ``0``)
    A list containing the indices of the items belonging to the first
   cluster. A cluster containing only one item *i* can be represented
   either as a list ``[i]``, or as an integer ``i``.
-  ``index2`` (default: ``0``)
    A list containing the indices of the items belonging to the second
   cluster. A cluster containing only one items *i* can be represented
   either as a list ``[i]``, or as an integer ``i``.
-  ``method`` (default: ``'a'``)
    Specifies how the distance between clusters is defined:

   -  ``'a'``: Distance between the two cluster centroids (arithmetic
      mean);
   -  ``'m'``: Distance between the two cluster centroids (median);
   -  ``'s'``: Shortest pairwise distance between items in the two
      clusters;
   -  ``'x'``: Longest pairwise distance between items in the two
      clusters;
   -  ``'v'``: Average over the pairwise distances between items in the
      two clusters.

-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).
-  ``transpose`` (default: ``0``)
    If ``transpose==0``, calculate the distance between the rows of
   ``data``. If ``transpose==1``, calculate the distance between the
   columns of ``data``.

15.3  Partitioning algorithms
-----------------------------

Partitioning algorithms divide items into *k* clusters such that the sum
of distances over the items to their cluster centers is minimal. The
number of clusters *k* is specified by the user. Three partitioning
algorithms are available in ``Bio.Cluster``:

-  *k*-means clustering
-  *k*-medians clustering
-  *k*-medoids clustering

These algorithms differ in how the cluster center is defined. In
*k*-means clustering, the cluster center is defined as the mean data
vector averaged over all items in the cluster. Instead of the mean, in
*k*-medians clustering the median is calculated for each dimension in
the data vector. Finally, in *k*-medoids clustering the cluster center
is defined as the item which has the smallest sum of distances to the
other items in the cluster. This clustering algorithm is suitable for
cases in which the distance matrix is known but the original data matrix
is not available, for example when clustering proteins based on their
structural similarity.

The expectation-maximization (EM) algorithm is used to find this
partitioning into *k* groups. In the initialization of the EM algorithm,
we randomly assign items to clusters. To ensure that no empty clusters
are produced, we use the binomial distribution to randomly choose the
number of items in each cluster to be one or more. We then randomly
permute the cluster assignments to items such that each item has an
equal probability to be in any cluster. Each cluster is thus guaranteed
to contain at least one item.

We then iterate:

-  Calculate the centroid of each cluster, defined as either the mean,
   the median, or the medoid of the cluster;
-  Calculate the distances of each item to the cluster centers;
-  For each item, determine which cluster centroid is closest;
-  Reassign each item to its closest cluster, or stop the iteration if
   no further item reassignments take place.

To avoid clusters becoming empty during the iteration, in *k*-means and
*k*-medians clustering the algorithm keeps track of the number of items
in each cluster, and prohibits the last remaining item in a cluster from
being reassigned to a different cluster. For *k*-medoids clustering,
such a check is not needed, as the item that functions as the cluster
centroid has a zero distance to itself, and will therefore never be
closer to a different cluster.

As the initial assignment of items to clusters is done randomly, usually
a different clustering solution is found each time the EM algorithm is
executed. To find the optimal clustering solution, the *k*-means
algorithm is repeated many times, each time starting from a different
initial random clustering. The sum of distances of the items to their
cluster center is saved for each run, and the solution with the smallest
value of this sum will be returned as the overall clustering solution.

How often the EM algorithm should be run depends on the number of items
being clustered. As a rule of thumb, we can consider how often the
optimal solution was found; this number is returned by the partitioning
algorithms as implemented in this library. If the optimal solution was
found many times, it is unlikely that better solutions exist than the
one that was found. However, if the optimal solution was found only
once, there may well be other solutions with a smaller within-cluster
sum of distances. If the number of items is large (more than several
hundreds), it may be difficult to find the globally optimal solution.

The EM algorithm terminates when no further reassignments take place. We
noticed that for some sets of initial cluster assignments, the EM
algorithm fails to converge due to the same clustering solution
reappearing periodically after a small number of iteration steps. We
therefore check for the occurrence of such periodic solutions during the
iteration. After a given number of iteration steps, the current
clustering result is saved as a reference. By comparing the clustering
result after each subsequent iteration step to the reference state, we
can determine if a previously encountered clustering result is found. In
such a case, the iteration is halted. If after a given number of
iterations the reference state has not yet been encountered, the current
clustering solution is saved to be used as the new reference state.
Initially, ten iteration steps are executed before resaving the
reference state. This number of iteration steps is doubled each time, to
ensure that periodic behavior with longer periods can also be detected.

*k*-means and *k*-medians
~~~~~~~~~~~~~~~~~~~~~~~~~

The *k*-means and *k*-medians algorithms are implemented as the function
``kcluster`` in ``Bio.Cluster``:

.. code:: verbatim

    >>> from Bio.Cluster import kcluster
    >>> clusterid, error, nfound = kcluster(data)

where the following arguments are defined:

-  ``data`` (required)
    Array containing the data for the items.
-  ``nclusters`` (default: ``2``)
    The number of clusters *k*.
-  ``mask`` (default: ``None``)
    Array of integers showing which data are missing. If
   ``mask[i,j]==0``, then ``data[i,j]`` is missing. If ``mask==None``,
   then all data are present.
-  ``weight`` (default: ``None``)
    The weights to be used when calculating distances. If
   ``weight==None``, then equal weights are assumed.
-  ``transpose`` (default: ``0``)
    Determines if rows (``transpose`` is ``0``) or columns
   (``transpose`` is ``1``) are to be clustered.
-  ``npass`` (default: ``1``)
    The number of times the *k*-means/-medians clustering algorithm is
   performed, each time with a different (random) initial condition. If
   ``initialid`` is given, the value of ``npass`` is ignored and the
   clustering algorithm is run only once, as it behaves
   deterministically in that case.
-  ``method`` (default: ``a``)
    describes how the center of a cluster is found:

   -  ``method=='a'``: arithmetic mean (*k*-means clustering);
   -  ``method=='m'``: median (*k*-medians clustering).

   For other values of ``method``, the arithmetic mean is used.
-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__). Whereas all eight distance
   measures are accepted by ``kcluster``, from a theoretical viewpoint
   it is best to use the Euclidean distance for the *k*-means algorithm,
   and the city-block distance for *k*-medians.
-  ``initialid`` (default: ``None``)
    Specifies the initial clustering to be used for the EM algorithm. If
   ``initialid==None``, then a different random initial clustering is
   used for each of the ``npass`` runs of the EM algorithm. If
   ``initialid`` is not ``None``, then it should be equal to a 1D array
   containing the cluster number (between ``0`` and ``nclusters-1``) for
   each item. Each cluster should contain at least one item. With the
   initial clustering specified, the EM algorithm is deterministic.

This function returns a tuple ``(clusterid, error, nfound)``, where
``clusterid`` is an integer array containing the number of the cluster
to which each row or cluster was assigned, ``error`` is the
within-cluster sum of distances for the optimal clustering solution, and
``nfound`` is the number of times this optimal solution was found.

*k*-medoids clustering
~~~~~~~~~~~~~~~~~~~~~~

The ``kmedoids`` routine performs *k*-medoids clustering on a given set
of items, using the distance matrix and the number of clusters passed by
the user:

.. code:: verbatim

    >>> from Bio.Cluster import kmedoids
    >>> clusterid, error, nfound = kmedoids(distance)

where the following arguments are defined: , nclusters=2, npass=1,
initialid=None)\|

-  ``distance`` (required)
    The matrix containing the distances between the items; this matrix
   can be specified in three ways:

   -  as a 2D Numerical Python array (in which only the left-lower part
      of the array will be accessed):

      .. code:: verbatim

          distance = array([[0.0, 1.1, 2.3],
                            [1.1, 0.0, 4.5],
                            [2.3, 4.5, 0.0]])

   -  as a 1D Numerical Python array containing consecutively the
      distances in the left-lower part of the distance matrix:

      .. code:: verbatim

          distance = array([1.1, 2.3, 4.5])

   -  as a list containing the rows of the left-lower part of the
      distance matrix:

      .. code:: verbatim

          distance = [array([]|,
                      array([1.1]),
                      array([2.3, 4.5])
                     ]

   These three expressions correspond to the same distance matrix.
-  ``nclusters`` (default: ``2``)
    The number of clusters *k*.
-  ``npass`` (default: ``1``)
    The number of times the *k*-medoids clustering algorithm is
   performed, each time with a different (random) initial condition. If
   ``initialid`` is given, the value of ``npass`` is ignored, as the
   clustering algorithm behaves deterministically in that case.
-  ``initialid`` (default: ``None``)
    Specifies the initial clustering to be used for the EM algorithm. If
   ``initialid==None``, then a different random initial clustering is
   used for each of the ``npass`` runs of the EM algorithm. If
   ``initialid`` is not ``None``, then it should be equal to a 1D array
   containing the cluster number (between ``0`` and ``nclusters-1``) for
   each item. Each cluster should contain at least one item. With the
   initial clustering specified, the EM algorithm is deterministic.

This function returns a tuple ``(clusterid, error, nfound)``, where
``clusterid`` is an array containing the number of the cluster to which
each item was assigned, ``error`` is the within-cluster sum of distances
for the optimal *k*-medoids clustering solution, and ``nfound`` is the
number of times the optimal solution was found. Note that the cluster
number in ``clusterid`` is defined as the item number of the item
representing the cluster centroid.

15.4  Hierarchical clustering
-----------------------------

Hierarchical clustering methods are inherently different from the
*k*-means clustering method. In hierarchical clustering, the similarity
in the expression profile between genes or experimental conditions are
represented in the form of a tree structure. This tree structure can be
shown graphically by programs such as Treeview and Java Treeview, which
has contributed to the popularity of hierarchical clustering in the
analysis of gene expression data.

The first step in hierarchical clustering is to calculate the distance
matrix, specifying all the distances between the items to be clustered.
Next, we create a node by joining the two closest items. Subsequent
nodes are created by pairwise joining of items or nodes based on the
distance between them, until all items belong to the same node. A tree
structure can then be created by retracing which items and nodes were
merged. Unlike the EM algorithm, which is used in *k*-means clustering,
the complete process of hierarchical clustering is deterministic.

Several flavors of hierarchical clustering exist, which differ in how
the distance between subnodes is defined in terms of their members. In
``Bio.Cluster``, pairwise single, maximum, average, and centroid linkage
are available.

-  In pairwise single-linkage clustering, the distance between two nodes
   is defined as the shortest distance among the pairwise distances
   between the members of the two nodes.
-  In pairwise maximum-linkage clustering, alternatively known as
   pairwise complete-linkage clustering, the distance between two nodes
   is defined as the longest distance among the pairwise distances
   between the members of the two nodes.
-  In pairwise average-linkage clustering, the distance between two
   nodes is defined as the average over all pairwise distances between
   the items of the two nodes.
-  In pairwise centroid-linkage clustering, the distance between two
   nodes is defined as the distance between their centroids. The
   centroids are calculated by taking the mean over all the items in a
   cluster. As the distance from each newly formed node to existing
   nodes and items need to be calculated at each step, the computing
   time of pairwise centroid-linkage clustering may be significantly
   longer than for the other hierarchical clustering methods. Another
   peculiarity is that (for a distance measure based on the Pearson
   correlation), the distances do not necessarily increase when going up
   in the clustering tree, and may even decrease. This is caused by an
   inconsistency between the centroid calculation and the distance
   calculation when using the Pearson correlation: Whereas the Pearson
   correlation effectively normalizes the data for the distance
   calculation, no such normalization occurs for the centroid
   calculation.

For pairwise single-, complete-, and average-linkage clustering, the
distance between two nodes can be found directly from the distances
between the individual items. Therefore, the clustering algorithm does
not need access to the original gene expression data, once the distance
matrix is known. For pairwise centroid-linkage clustering, however, the
centroids of newly formed subnodes can only be calculated from the
original data and not from the distance matrix.

The implementation of pairwise single-linkage hierarchical clustering is
based on the SLINK algorithm (R. Sibson, 1973), which is much faster and
more memory-efficient than a straightforward implementation of pairwise
single-linkage clustering. The clustering result produced by this
algorithm is identical to the clustering solution found by the
conventional single-linkage algorithm. The single-linkage hierarchical
clustering algorithm implemented in this library can be used to cluster
large gene expression data sets, for which conventional hierarchical
clustering algorithms fail due to excessive memory requirements and
running time.

Representing a hierarchical clustering solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The result of hierarchical clustering consists of a tree of nodes, in
which each node joins two items or subnodes. Usually, we are not only
interested in which items or subnodes are joined at each node, but also
in their similarity (or distance) as they are joined. To store one node
in the hierarchical clustering tree, we make use of the class ``Node``,
which defined in ``Bio.Cluster``. An instance of ``Node`` has three
attributes:

-  ``left``
-  ``right``
-  ``distance``

Here, ``left`` and ``right`` are integers referring to the two items or
subnodes that are joined at this node, and ``distance`` is the distance
between them. The items being clustered are numbered from 0 to (number
of items − 1), while clusters are numbered from -1 to −(number of
items−1). Note that the number of nodes is one less than the number of
items.

To create a new ``Node`` object, we need to specify ``left`` and
``right``; ``distance`` is optional.

.. code:: verbatim

    >>> from Bio.Cluster import Node
    >>> Node(2,3)
    (2, 3): 0
    >>> Node(2,3,0.91)
    (2, 3): 0.91

The attributes ``left``, ``right``, and ``distance`` of an existing
``Node`` object can be modified directly:

.. code:: verbatim

    >>> node = Node(4,5)
    >>> node.left = 6
    >>> node.right = 2
    >>> node.distance = 0.73
    >>> node
    (6, 2): 0.73

An error is raised if ``left`` and ``right`` are not integers, or if
``distance`` cannot be converted to a floating-point value.

The Python class ``Tree`` represents a full hierarchical clustering
solution. A ``Tree`` object can be created from a list of ``Node``
objects:

.. code:: verbatim

    >>> from Bio.Cluster import Node, Tree
    >>> nodes = [Node(1,2,0.2), Node(0,3,0.5), Node(-2,4,0.6), Node(-1,-3,0.9)]
    >>> tree = Tree(nodes)
    >>> print tree
    (1, 2): 0.2
    (0, 3): 0.5
    (-2, 4): 0.6
    (-1, -3): 0.9

The ``Tree`` initializer checks if the list of nodes is a valid
hierarchical clustering result:

.. code:: verbatim

    >>> nodes = [Node(1,2,0.2), Node(0,2,0.5)]
    >>> Tree(nodes)
    Traceback (most recent call last):
      File "<stdin>", line 1, in ?
    ValueError: Inconsistent tree

Individual nodes in a ``Tree`` object can be accessed using square
brackets:

.. code:: verbatim

    >>> nodes = [Node(1,2,0.2), Node(0,-1,0.5)]
    >>> tree = Tree(nodes)
    >>> tree[0]
    (1, 2): 0.2
    >>> tree[1]
    (0, -1): 0.5
    >>> tree[-1]
    (0, -1): 0.5

As a ``Tree`` object is read-only, we cannot change individual nodes in
a ``Tree`` object. However, we can convert the tree to a list of nodes,
modify this list, and create a new tree from this list:

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

This guarantees that any ``Tree`` object is always well-formed.

To display a hierarchical clustering solution with visualization
programs such as Java Treeview, it is better to scale all node distances
such that they are between zero and one. This can be accomplished by
calling the ``scale`` method on an existing ``Tree`` object:

.. code:: verbatim

    >>> tree.scale()

This method takes no arguments, and returns ``None``.

After hierarchical clustering, the items can be grouped into *k*
clusters based on the tree structure stored in the ``Tree`` object by
cutting the tree:

.. code:: verbatim

    >>> clusterid = tree.cut(nclusters=1)

where ``nclusters`` (defaulting to ``1``) is the desired number of
clusters *k*. This method ignores the top *k*\ −1 linking events in the
tree structure, resulting in *k* separated clusters of items. The number
of clusters *k* should be positive, and less than or equal to the number
of items. This method returns an array ``clusterid`` containing the
number of the cluster to which each item is assigned.

Performing hierarchical clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform hierarchical clustering, use the ``treecluster`` function in
``Bio.Cluster``.

.. code:: verbatim

    >>> from Bio.Cluster import treecluster
    >>> tree = treecluster(data)

where the following arguments are defined:

-  ``data``
    Array containing the data for the items.
-  ``mask`` (default: ``None``)
    Array of integers showing which data are missing. If
   ``mask[i,j]==0``, then ``data[i,j]`` is missing. If ``mask==None``,
   then all data are present.
-  ``weight`` (default: ``None``)
    The weights to be used when calculating distances. If
   ``weight==None``, then equal weights are assumed.
-  ``transpose`` (default: ``0``)
    Determines if rows (``transpose==0``) or columns (``transpose==1``)
   are to be clustered.
-  ``method`` (default: ``'m'``)
    defines the linkage method to be used:

   -  ``method=='s'``: pairwise single-linkage clustering
   -  ``method=='m'``: pairwise maximum- (or complete-) linkage
      clustering
   -  ``method=='c'``: pairwise centroid-linkage clustering
   -  ``method=='a'``: pairwise average-linkage clustering

-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).

To apply hierarchical clustering on a precalculated distance matrix,
specify the ``distancematrix`` argument when calling ``treecluster``
function instead of the ``data`` argument:

.. code:: verbatim

    >>> from Bio.Cluster import treecluster
    >>> tree = treecluster(distancematrix=distance)

In this case, the following arguments are defined:

-  ``distancematrix``
    The distance matrix, which can be specified in three ways:

   -  as a 2D Numerical Python array (in which only the left-lower part
      of the array will be accessed):

      .. code:: verbatim

          distance = array([[0.0, 1.1, 2.3], 
                            [1.1, 0.0, 4.5],
                            [2.3, 4.5, 0.0]])

   -  as a 1D Numerical Python array containing consecutively the
      distances in the left-lower part of the distance matrix:

      .. code:: verbatim

          distance = array([1.1, 2.3, 4.5])

   -  as a list containing the rows of the left-lower part of the
      distance matrix:

      .. code:: verbatim

          distance = [array([]),
                      array([1.1]),
                      array([2.3, 4.5])

   These three expressions correspond to the same distance matrix. As
   ``treecluster`` may shuffle the values in the distance matrix as part
   of the clustering algorithm, be sure to save this array in a
   different variable before calling ``treecluster`` if you need it
   later.
-  ``method``
    The linkage method to be used:

   -  ``method=='s'``: pairwise single-linkage clustering
   -  ``method=='m'``: pairwise maximum- (or complete-) linkage
      clustering
   -  ``method=='a'``: pairwise average-linkage clustering

   While pairwise single-, maximum-, and average-linkage clustering can
   be calculated from the distance matrix alone, pairwise
   centroid-linkage cannot.

When calling ``treecluster``, either ``data`` or ``distancematrix``
should be ``None``.

This function returns a ``Tree`` object. This object contains (number of
items − 1) nodes, where the number of items is the number of rows if
rows were clustered, or the number of columns if columns were clustered.
Each node describes a pairwise linking event, where the node attributes
``left`` and ``right`` each contain the number of one item or subnode,
and ``distance`` the distance between them. Items are numbered from 0 to
(number of items − 1), while clusters are numbered -1 to −(number of
items−1).

15.5  Self-Organizing Maps
--------------------------

Self-Organizing Maps (SOMs) were invented by Kohonen to describe neural
networks (see for instance Kohonen, 1997 [`24 <#kohonen1997>`__\ ]).
Tamayo (1999) first applied Self-Organizing Maps to gene expression data
[`30 <#tamayo1999>`__\ ].

SOMs organize items into clusters that are situated in some topology.
Usually a rectangular topology is chosen. The clusters generated by SOMs
are such that neighboring clusters in the topology are more similar to
each other than clusters far from each other in the topology.

The first step to calculate a SOM is to randomly assign a data vector to
each cluster in the topology. If rows are being clustered, then the
number of elements in each data vector is equal to the number of
columns.

An SOM is then generated by taking rows one at a time, and finding which
cluster in the topology has the closest data vector. The data vector of
that cluster, as well as those of the neighboring clusters, are adjusted
using the data vector of the row under consideration. The adjustment is
given by

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

The parameter τ is a parameter that decreases at each iteration step. We
have used a simple linear function of the iteration step:

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

τ\ :sub:`init` is the initial value of τ as specified by the user, *i*
is the number of the current iteration step, and *n* is the total number
of iteration steps to be performed. While changes are made rapidly in
the beginning of the iteration, at the end of iteration only small
changes are made.

All clusters within a radius *R* are adjusted to the gene under
consideration. This radius decreases as the calculation progresses as

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

in which the maximum radius is defined as

*R*\ :sub:`max` = 

√

+---------------------------------------------------------+
+---------------------------------------------------------+
| *N*\ :sub:`*x*`\ :sup:`2` + *N*\ :sub:`*y*`\ :sup:`2`   |
+---------------------------------------------------------+

,

where (*N*\ :sub:`*x*`, *N*\ :sub:`*y*`) are the dimensions of the
rectangle defining the topology.

The function ``somcluster`` implements the complete algorithm to
calculate a Self-Organizing Map on a rectangular grid. First it
initializes the random number generator. The node data are then
initialized using the random number generator. The order in which genes
or microarrays are used to modify the SOM is also randomized. The total
number of iterations in the SOM algorithm is specified by the user.

To run ``somcluster``, use

.. code:: verbatim

    >>> from Bio.Cluster import somcluster
    >>> clusterid, celldata = somcluster(data)

where the following arguments are defined:

-  ``data`` (required)
    Array containing the data for the items.
-  ``mask`` (default: ``None``)
    Array of integers showing which data are missing. If
   ``mask[i,j]==0``, then ``data[i,j]`` is missing. If ``mask==None``,
   then all data are present.
-  ``weight`` (default: ``None``)
    contains the weights to be used when calculating distances. If
   ``weight==None``, then equal weights are assumed.
-  ``transpose`` (default: ``0``)
    Determines if rows (``transpose`` is ``0``) or columns
   (``transpose`` is ``1``) are to be clustered.
-  ``nxgrid, nygrid`` (default: ``2, 1``)
    The number of cells horizontally and vertically in the rectangular
   grid on which the Self-Organizing Map is calculated.
-  ``inittau`` (default: ``0.02``)
    The initial value for the parameter τ that is used in the SOM
   algorithm. The default value for ``inittau`` is 0.02, which was used
   in Michael Eisen’s Cluster/TreeView program.
-  ``niter`` (default: ``1``)
    The number of iterations to be performed.
-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).

This function returns the tuple ``(clusterid, celldata)``:

-  ``clusterid``:
    An array with two columns, where the number of rows is equal to the
   number of items that were clustered. Each row contains the *x* and
   *y* coordinates of the cell in the rectangular SOM grid to which the
   item was assigned.
-  ``celldata``:
    An array with dimensions (``nxgrid``, ``nygrid``, number of columns)
   if rows are being clustered, or (``nxgrid``, ``nygrid``, number of
   rows) if columns are being clustered. Each element ``[ix][iy]`` of
   this array is a 1D vector containing the gene expression data for the
   centroid of the cluster in the grid cell with coordinates
   ``[ix][iy]``.

15.6  Principal Component Analysis
----------------------------------

Principal Component Analysis (PCA) is a widely used technique for
analyzing multivariate data. A practical example of applying Principal
Component Analysis to gene expression data is presented by Yeung and
Ruzzo (2001) [`33 <#yeung2001>`__\ ].

In essence, PCA is a coordinate transformation in which each row in the
data matrix is written as a linear sum over basis vectors called
principal components, which are ordered and chosen such that each
maximally explains the remaining variance in the data vectors. For
example, an *n* × 3 data matrix can be represented as an ellipsoidal
cloud of *n* points in three dimensional space. The first principal
component is the longest axis of the ellipsoid, the second principal
component the second longest axis of the ellipsoid, and the third
principal component is the shortest axis. Each row in the data matrix
can be reconstructed as a suitable linear combination of the principal
components. However, in order to reduce the dimensionality of the data,
usually only the most important principal components are retained. The
remaining variance present in the data is then regarded as unexplained
variance.

The principal components can be found by calculating the eigenvectors of
the covariance matrix of the data. The corresponding eigenvalues
determine how much of the variance present in the data is explained by
each principal component.

Before applying principal component analysis, typically the mean is
subtracted from each column in the data matrix. In the example above,
this effectively centers the ellipsoidal cloud around its centroid in 3D
space, with the principal components describing the variation of points
in the ellipsoidal cloud with respect to their centroid.

The function ``pca`` below first uses the singular value decomposition
to calculate the eigenvalues and eigenvectors of the data matrix. The
singular value decomposition is implemented as a translation in C of the
Algol procedure ``svd`` [`16 <#golub1971>`__\ ], which uses Householder
bidiagonalization and a variant of the QR algorithm. The principal
components, the coordinates of each data vector along the principal
components, and the eigenvalues corresponding to the principal
components are then evaluated and returned in decreasing order of the
magnitude of the eigenvalue. If data centering is desired, the mean
should be subtracted from each column in the data matrix before calling
the ``pca`` routine.

To apply Principal Component Analysis to a rectangular matrix ``data``,
use

.. code:: verbatim

    >>> from Bio.Cluster import pca
    >>> columnmean, coordinates, components, eigenvalues = pca(data)

This function returns a tuple
``columnmean, coordinates, components, eigenvalues``:

-  ``columnmean``
    Array containing the mean over each column in ``data``.
-  ``coordinates``
    The coordinates of each row in ``data`` with respect to the
   principal components.
-  ``components``
    The principal components.
-  ``eigenvalues``
    The eigenvalues corresponding to each of the principal components.

The original matrix ``data`` can be recreated by calculating
``columnmean +  dot(coordinates, components)``.

15.7  Handling Cluster/TreeView-type files
------------------------------------------

Cluster/TreeView are GUI-based codes for clustering gene expression
data. They were originally written by `Michael
Eisen <http://rana.lbl.gov>`__ while at Stanford University.
``Bio.Cluster`` contains functions for reading and writing data files
that correspond to the format specified for Cluster/TreeView. In
particular, by saving a clustering result in that format, TreeView can
be used to visualize the clustering results. We recommend using Alok
Saldanha’s
```http://jtreeview.sourceforge.net/`` <http://jtreeview.sourceforge.net/>`__\ Java
TreeView program, which can display hierarchical as well as *k*-means
clustering results.

An object of the class ``Record`` contains all information stored in a
Cluster/TreeView-type data file. To store the information contained in
the data file in a ``Record`` object, we first open the file and then
read it:

.. code:: verbatim

    >>> from Bio import Cluster
    >>> handle = open("mydatafile.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()

This two-step process gives you some flexibility in the source of the
data. For example, you can use

.. code:: verbatim

    >>> import gzip # Python standard library
    >>> handle = gzip.open("mydatafile.txt.gz")

to open a gzipped file, or

.. code:: verbatim

    >>> import urllib # Python standard library
    >>> handle = urllib.urlopen("http://somewhere.org/mydatafile.txt")

to open a file stored on the Internet before calling ``read``.

The ``read`` command reads the tab-delimited text file
``mydatafile.txt`` containing gene expression data in the format
specified for Michael Eisen’s Cluster/TreeView program. For a
description of this file format, see the manual to Cluster/TreeView. It
is available at `Michael Eisen’s lab
website <http://rana.lbl.gov/manuals/ClusterTreeView.pdf>`__ and at `our
website <http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/cluster3.pdf>`__.

A ``Record`` object has the following attributes:

-  ``data``
    The data array containing the gene expression data. Genes are stored
   row-wise, while microarrays are stored column-wise.
-  ``mask``
    This array shows which elements in the ``data`` array, if any, are
   missing. If ``mask[i,j]==0``, then ``data[i,j]`` is missing. If no
   data were found to be missing, ``mask`` is set to ``None``.
-  ``geneid``
    This is a list containing a unique description for each gene (i.e.,
   ORF numbers).
-  ``genename``
    This is a list containing a description for each gene (i.e., gene
   name). If not present in the data file, ``genename`` is set to
   ``None``.
-  ``gweight``
    The weights that are to be used to calculate the distance in
   expression profile between genes. If not present in the data file,
   ``gweight`` is set to ``None``.
-  ``gorder``
    The preferred order in which genes should be stored in an output
   file. If not present in the data file, ``gorder`` is set to ``None``.
-  ``expid``
    This is a list containing a description of each microarray, e.g.
   experimental condition.
-  ``eweight``
    The weights that are to be used to calculate the distance in
   expression profile between microarrays. If not present in the data
   file, ``eweight`` is set to ``None``.
-  ``eorder``
    The preferred order in which microarrays should be stored in an
   output file. If not present in the data file, ``eorder`` is set to
   ``None``.
-  ``uniqid``
    The string that was used instead of UNIQID in the data file.

After loading a ``Record`` object, each of these attributes can be
accessed and modified directly. For example, the data can be
log-transformed by taking the logarithm of ``record.data``.

Calculating the distance matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To calculate the distance matrix between the items stored in the record,
use

.. code:: verbatim

    >>> matrix = record.distancematrix()

where the following arguments are defined:

-  ``transpose`` (default: ``0``)
    Determines if the distances between the rows of ``data`` are to be
   calculated (``transpose==0``), or between the columns of ``data``
   (``transpose==1``).
-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).

This function returns the distance matrix as a list of rows, where the
number of columns of each row is equal to the row number (see section
`15.1 <#subsec:distancematrix>`__).

Calculating the cluster centroids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To calculate the centroids of clusters of items stored in the record,
use

.. code:: verbatim

    >>> cdata, cmask = record.clustercentroids()

-  ``clusterid`` (default: ``None``)
    Vector of integers showing to which cluster each item belongs. If
   ``clusterid`` is not given, then all items are assumed to belong to
   the same cluster.
-  ``method`` (default: ``'a'``)
    Specifies whether the arithmetic mean (``method=='a'``) or the
   median (``method=='m'``) is used to calculate the cluster center.
-  ``transpose`` (default: ``0``)
    Determines if the centroids of the rows of ``data`` are to be
   calculated (``transpose==0``), or the centroids of the columns of
   ``data`` (``transpose==1``).

This function returns the tuple ``cdata, cmask``; see section
`15.2 <#subsec:clustercentroids>`__ for a description.

Calculating the distance between clusters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To calculate the distance between clusters of items stored in the
record, use

.. code:: verbatim

    >>> distance = record.clusterdistance()

where the following arguments are defined:

-  ``index1`` (default: ``0``)
    A list containing the indices of the items belonging to the first
   cluster. A cluster containing only one item *i* can be represented
   either as a list ``[i]``, or as an integer ``i``.
-  ``index2`` (default: ``0``)
    A list containing the indices of the items belonging to the second
   cluster. A cluster containing only one item *i* can be represented
   either as a list ``[i]``, or as an integer ``i``.
-  ``method`` (default: ``'a'``)
    Specifies how the distance between clusters is defined:

   -  ``'a'``: Distance between the two cluster centroids (arithmetic
      mean);
   -  ``'m'``: Distance between the two cluster centroids (median);
   -  ``'s'``: Shortest pairwise distance between items in the two
      clusters;
   -  ``'x'``: Longest pairwise distance between items in the two
      clusters;
   -  ``'v'``: Average over the pairwise distances between items in the
      two clusters.

-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).
-  ``transpose`` (default: ``0``)
    If ``transpose==0``, calculate the distance between the rows of
   ``data``. If ``transpose==1``, calculate the distance between the
   columns of ``data``.

Performing hierarchical clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform hierarchical clustering on the items stored in the record,
use

.. code:: verbatim

    >>> tree = record.treecluster()

where the following arguments are defined:

-  ``transpose`` (default: ``0``)
    Determines if rows (``transpose==0``) or columns (``transpose==1``)
   are to be clustered.
-  ``method`` (default: ``'m'``)
    defines the linkage method to be used:

   -  ``method=='s'``: pairwise single-linkage clustering
   -  ``method=='m'``: pairwise maximum- (or complete-) linkage
      clustering
   -  ``method=='c'``: pairwise centroid-linkage clustering
   -  ``method=='a'``: pairwise average-linkage clustering

-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).
-  ``transpose``
    Determines if genes or microarrays are being clustered. If
   ``transpose==0``, genes (rows) are being clustered. If
   ``transpose==1``, microarrays (columns) are clustered.

This function returns a ``Tree`` object. This object contains (number of
items − 1) nodes, where the number of items is the number of rows if
rows were clustered, or the number of columns if columns were clustered.
Each node describes a pairwise linking event, where the node attributes
``left`` and ``right`` each contain the number of one item or subnode,
and ``distance`` the distance between them. Items are numbered from 0 to
(number of items − 1), while clusters are numbered -1 to −(number of
items−1).

Performing *k*-means or *k*-medians clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform *k*-means or *k*-medians clustering on the items stored in
the record, use

.. code:: verbatim

    >>> clusterid, error, nfound = record.kcluster()

where the following arguments are defined:

-  ``nclusters`` (default: ``2``)
    The number of clusters *k*.
-  ``transpose`` (default: ``0``)
    Determines if rows (``transpose`` is ``0``) or columns
   (``transpose`` is ``1``) are to be clustered.
-  ``npass`` (default: ``1``)
    The number of times the *k*-means/-medians clustering algorithm is
   performed, each time with a different (random) initial condition. If
   ``initialid`` is given, the value of ``npass`` is ignored and the
   clustering algorithm is run only once, as it behaves
   deterministically in that case.
-  ``method`` (default: ``a``)
    describes how the center of a cluster is found:

   -  ``method=='a'``: arithmetic mean (*k*-means clustering);
   -  ``method=='m'``: median (*k*-medians clustering).

   For other values of ``method``, the arithmetic mean is used.
-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).

This function returns a tuple ``(clusterid, error, nfound)``, where
``clusterid`` is an integer array containing the number of the cluster
to which each row or cluster was assigned, ``error`` is the
within-cluster sum of distances for the optimal clustering solution, and
``nfound`` is the number of times this optimal solution was found.

Calculating a Self-Organizing Map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To calculate a Self-Organizing Map of the items stored in the record,
use

.. code:: verbatim

    >>> clusterid, celldata = record.somcluster()

where the following arguments are defined:

-  ``transpose`` (default: ``0``)
    Determines if rows (``transpose`` is ``0``) or columns
   (``transpose`` is ``1``) are to be clustered.
-  ``nxgrid, nygrid`` (default: ``2, 1``)
    The number of cells horizontally and vertically in the rectangular
   grid on which the Self-Organizing Map is calculated.
-  ``inittau`` (default: ``0.02``)
    The initial value for the parameter τ that is used in the SOM
   algorithm. The default value for ``inittau`` is 0.02, which was used
   in Michael Eisen’s Cluster/TreeView program.
-  ``niter`` (default: ``1``)
    The number of iterations to be performed.
-  ``dist`` (default: ``'e'``, Euclidean distance)
    Defines the distance function to be used (see
   `15.1 <#sec:distancefunctions>`__).

This function returns the tuple ``(clusterid, celldata)``:

-  ``clusterid``:
    An array with two columns, where the number of rows is equal to the
   number of items that were clustered. Each row contains the *x* and
   *y* coordinates of the cell in the rectangular SOM grid to which the
   item was assigned.
-  ``celldata``:
    An array with dimensions (``nxgrid``, ``nygrid``, number of columns)
   if rows are being clustered, or (``nxgrid``, ``nygrid``, number of
   rows) if columns are being clustered. Each element ``[ix][iy]`` of
   this array is a 1D vector containing the gene expression data for the
   centroid of the cluster in the grid cell with coordinates
   ``[ix][iy]``.

Saving the clustering result
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To save the clustering result, use

.. code:: verbatim

    >>> record.save(jobname, geneclusters, expclusters)

where the following arguments are defined:

-  ``jobname``
    The string ``jobname`` is used as the base name for names of the
   files that are to be saved.
-  ``geneclusters``
    This argument describes the gene (row-wise) clustering result. In
   case of *k*-means clustering, this is a 1D array containing the
   number of the cluster each gene belongs to. It can be calculated
   using ``kcluster``. In case of hierarchical clustering,
   ``geneclusters`` is a ``Tree`` object.
-  ``expclusters``
    This argument describes the (column-wise) clustering result for the
   experimental conditions. In case of *k*-means clustering, this is a
   1D array containing the number of the cluster each experimental
   condition belongs to. It can be calculated using ``kcluster``. In
   case of hierarchical clustering, ``expclusters`` is a ``Tree``
   object.

This method writes the text file ``jobname.cdt``, ``jobname.gtr``,
``jobname.atr``, ``jobname*.kgg``, and/or ``jobname*.kag`` for
subsequent reading by the Java TreeView program. If ``geneclusters`` and
``expclusters`` are both ``None``, this method only writes the text file
``jobname.cdt``; this file can subsequently be read into a new
``Record`` object.

15.8  Example calculation
-------------------------

This is an example of a hierarchical clustering calculation, using
single linkage clustering for genes and maximum linkage clustering for
experimental conditions. As the Euclidean distance is being used for
gene clustering, it is necessary to scale the node distances
``genetree`` such that they are all between zero and one. This is needed
for the Java TreeView code to display the tree diagram correctly. To
cluster the experimental conditions, the uncentered correlation is being
used. No scaling is needed in this case, as the distances in ``exptree``
are already between zero and two. The example data ``cyano.txt`` can be
found in the ``data`` subdirectory.

.. code:: verbatim

    >>> from Bio import Cluster
    >>> handle = open("cyano.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()
    >>> genetree = record.treecluster(method='s')
    >>> genetree.scale()
    >>> exptree = record.treecluster(dist='u', transpose=1)
    >>> record.save("cyano_result", genetree, exptree)

This will create the files ``cyano_result.cdt``, ``cyano_result.gtr``,
and ``cyano_result.atr``.

Similarly, we can save a *k*-means clustering solution:

.. code:: verbatim

    >>> from Bio import Cluster
    >>> handle = open("cyano.txt")
    >>> record = Cluster.read(handle)
    >>> handle.close()
    >>> (geneclusters, error, ifound) = record.kcluster(nclusters=5, npass=1000)
    >>> (expclusters, error, ifound) = record.kcluster(nclusters=2, npass=100, transpose=1)
    >>> record.save("cyano_result", geneclusters, expclusters)

This will create the files ``cyano_result_K_G2_A2.cdt``,
``cyano_result_K_G2.kgg``, and ``cyano_result_K_A2.kag``.

15.9  Auxiliary functions
-------------------------

``median(data)`` returns the median of the 1D array ``data``.

``mean(data)`` returns the mean of the 1D array ``data``.

``version()`` returns the version number of the underlying C Clustering
Library as a string.

