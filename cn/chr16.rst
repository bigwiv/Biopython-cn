第十六章 监督学习方法
=======================================

注意本章介绍的所有监督学习方法都需要先安装Numerical Python （numpy）。

16.1 Logistic 回归模型
-----------------------------------

16.1.1 背景和目的
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Logistic回归是一种监督学习方法，通过若干预测变量 *x*\ :sub:`*i*` 的加权和尝试将样本划分为 *K* 个不同类别。Logistic回归模型可用来计算预测变量的权重β\ :sub:`*i*`。在Biopython中，logistic回归模型目前只实现了二类别（*K*=2）分类，而预测变量的数量没有限制。

作为一个例子，我们试着预测细菌中的操纵子结构。一个操纵子是在一条DNA链上许多相邻基因组成的一个集合，可以被共同转录为一条mRNA分子，这条mRNA分子经翻译后产生多个不同的蛋白质。我们将以枯草芽孢杆菌的操纵子数据进行说明，它的一个操纵子平均包含2.4个基因。

首先为了理解细菌的基因调节，我们需要知道其操纵子的结构。枯草芽孢杆菌大约10%的基因操纵子结构已经通过实验获知。剩下的90%的基因操纵子结构可以通过一种监督学习方法来预测。

在这种方法中，我们需要选择某些与操纵子结构有关的容易计算的预测变量 *x*\ :sub:`*i*` 。例如可以选择基因间碱基对距离来来作为其中一个预测变量。同一个操纵子中的相邻基因往往距离相对较近，而位于不同操纵子的相邻基因间通常具有更大的空间来容纳启动子和终止子序列。另一个预测变量可以基于基因表达量度。根据操纵子的定义，属于同一个操纵子的基因有相同的基因表达轮廓，而不同操纵子的两个基因的表达轮廓也不相同。在实际操作中，由于存在测量误差，对相同操纵子的基因表达轮廓的测量不会完全一致。为了测量基因表达轮廓的相似性，我们假设测量误差服从正态分布，然后计算对应的对数似然分。

现在我们有了两个预测变量，可以据此预测在同一条DNA链上两个相邻基因是否属于相同的操纵子：
-  *x*\ :sub:`1`：两基因间的碱基对数
-  *x*\ :sub:`2`：两基因表达轮廓的相似度

在logistic回归模型中，我们使用这两个预测变量的加权和来计算一个联合得分 *S*：

+-----------------------------------------------------------------------------------+
| *S* = β:sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`.     (16.1)   |
+-----------------------------------------------------------------------------------+

根据下面两组示例基因，logistic回归模型对参数 β\ :sub:`0` ， β\ :sub:`1`, β\ :sub:`2` 给出合适的值：
-  OP: 相邻基因，相同DNA链，属于相同操纵子
-  NOP: 相邻基因，相同DNA链，属于不同操纵子

在logistic回归模型中，属于某个类别的概率依赖于通过logistic函数得出的分数。对于这两类OP和NOP，相应概率可如下表述：

     

Pr(\ *OP*\ \|\ *x*\ :sub:`1`, \ *x*\ :sub:`2`)

 =

 

+--------------------------------------------------------------------------+
| exp(β\ :sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`)     |
+--------------------------------------------------------------------------+
+--------------------------------------------------------------------------+
| 1+exp(β\ :sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`)   |
+--------------------------------------------------------------------------+

   

    (16.2)

Pr(\ *NOP*\ \|\ *x*\ :sub:`1`, \ *x*\ :sub:`2`)

 =

 

+--------------------------------------------------------------------------+
| 1                                                                        |
+--------------------------------------------------------------------------+
+--------------------------------------------------------------------------+
| 1+exp(β\ :sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`)   |
+--------------------------------------------------------------------------+

   

    (16.3)

使用一组已知是否属于相同操纵子（OP类别）或不同操纵子（NOP类别）的基因对，通过最大化相应概率函数的对数似然值，我们可以计算权重 β\ :sub:`0`, β\ :sub:`1`, β\ :sub:`2`。
(`16.2 <#eq:OP>`__) and (`16.3 <#eq:NOP>`__).

16.1.2 训练logistic回归模型
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    --------------

    +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | 表16.1： 已知类别(OP or NOP)的相邻基因对.如果两个基因相重叠，其基因间距离为负值   |
    +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | 基因对           | 基因间距离 (*x*\ :sub:`1`)   | 基因表达得分 (*x*\ :sub:`2`)   | 类别   |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *cotJA* — *cotJB*   | -53                                  | -200.78                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yesK* — *yesL*     | 117                                  | -267.14                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplA* — *lplB*     | 57                                   | -163.47                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplB* — *lplC*     | 16                                   | -190.30                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplC* — *lplD*     | 11                                   | -220.94                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplD* — *yetF*     | 85                                   | -193.94                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfmT* — *yfmS*     | 16                                   | -182.71                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfmF* — *yfmE*     | 15                                   | -180.41                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *citS* — *citT*     | -26                                  | -181.73                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *citM* — *yflN*     | 58                                   | -259.87                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfiI* — *yfiJ*     | 126                                  | -414.53                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lipB* — *yfiQ*     | 191                                  | -249.57                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfiU* — *yfiV*     | 113                                  | -265.28                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfhH* — *yfhI*     | 145                                  | -312.99                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *cotY* — *cotX*     | 154                                  | -213.83                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yjoB* — *rapA*     | 147                                  | -380.85                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *ptsI* — *splA*     | 93                                   | -291.13                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+

    --------------

表`16.1 <#table:training>`__ 列出了枯草芽孢杆菌的一些基因对，这些基因的操纵子结构已知。让我们根据表中的这些数据来计算其logistic回归模型：

.. code:: verbatim

    >>> from Bio import LogisticRegression
    >>> xs = [[-53, -200.78],
              [117, -267.14],
              [57, -163.47],
              [16, -190.30],
              [11, -220.94],
              [85, -193.94],
              [16, -182.71],
              [15, -180.41],
              [-26, -181.73],
              [58, -259.87],
              [126, -414.53],
              [191, -249.57],
              [113, -265.28],
              [145, -312.99],
              [154, -213.83],
              [147, -380.85],
              [93, -291.13]]
    >>> ys = [1,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0]
    >>> model = LogisticRegression.train(xs, ys)

这里， ``xs`` 和 ``ys`` 是训练数据： ``xs`` 包含每个基因对的预测变量， ``ys`` 指定是否这个基因对属于相同操纵子（ ``1`` ，类别OP）或不同操纵子（``0``，类别NOP）。Logistic回归模型结果存储在 ``model`` 中，包含权重 β\ :sub:`0`, β\ :sub:`1`, and β\ :sub:`2`:

.. code:: verbatim

    >>> model.beta
    [8.9830290157144681, -0.035968960444850887, 0.02181395662983519]

注意 β\ :sub:`1` 是负的，这是因为具有更短基因间距离的基因对有更高的概率属于相同操纵子（类别OP）。另一方面， β\ :sub:`2` 为正，因为属于相同操纵子的基因对通常有更高的基因表达轮廓相似性得分。参数 β\ :sub:`0` 是正值是因为在这个训练数据中操纵子基因对占据大多数。

函数 ``train`` 有两个可选参数： ``update_fn`` 和 ``typecode`` 。 ``update_fn`` 可用来指定一个回调函数，以迭代数和对数似然值做参数。在这个例子中，我们可以使用这个回调函数追踪模型计算（使用Newton-Raphson迭代来最大化logistic回归模型的对数似然函数）进度：

.. code:: verbatim

    >>> def show_progress(iteration, loglikelihood):
            print "Iteration:", iteration, "Log-likelihood function:", loglikelihood
    >>>
    >>> model = LogisticRegression.train(xs, ys, update_fn=show_progress)
    Iteration: 0 Log-likelihood function: -11.7835020695
    Iteration: 1 Log-likelihood function: -7.15886767672
    Iteration: 2 Log-likelihood function: -5.76877209868
    Iteration: 3 Log-likelihood function: -5.11362294338
    Iteration: 4 Log-likelihood function: -4.74870642433
    Iteration: 5 Log-likelihood function: -4.50026077146
    Iteration: 6 Log-likelihood function: -4.31127773737
    Iteration: 7 Log-likelihood function: -4.16015043396
    Iteration: 8 Log-likelihood function: -4.03561719785
    Iteration: 9 Log-likelihood function: -3.93073282192
    Iteration: 10 Log-likelihood function: -3.84087660929
    Iteration: 11 Log-likelihood function: -3.76282560605
    Iteration: 12 Log-likelihood function: -3.69425027154
    Iteration: 13 Log-likelihood function: -3.6334178602
    Iteration: 14 Log-likelihood function: -3.57900855837
    Iteration: 15 Log-likelihood function: -3.52999671386
    Iteration: 16 Log-likelihood function: -3.48557145163
    Iteration: 17 Log-likelihood function: -3.44508206139
    Iteration: 18 Log-likelihood function: -3.40799948447
    Iteration: 19 Log-likelihood function: -3.3738885624
    Iteration: 20 Log-likelihood function: -3.3423876581
    Iteration: 21 Log-likelihood function: -3.31319343769
    Iteration: 22 Log-likelihood function: -3.2860493346
    Iteration: 23 Log-likelihood function: -3.2607366863
    Iteration: 24 Log-likelihood function: -3.23706784091
    Iteration: 25 Log-likelihood function: -3.21488073614
    Iteration: 26 Log-likelihood function: -3.19403459259
    Iteration: 27 Log-likelihood function: -3.17440646052
    Iteration: 28 Log-likelihood function: -3.15588842703
    Iteration: 29 Log-likelihood function: -3.13838533947
    Iteration: 30 Log-likelihood function: -3.12181293595
    Iteration: 31 Log-likelihood function: -3.10609629966
    Iteration: 32 Log-likelihood function: -3.09116857282
    Iteration: 33 Log-likelihood function: -3.07696988017
    Iteration: 34 Log-likelihood function: -3.06344642288
    Iteration: 35 Log-likelihood function: -3.05054971191
    Iteration: 36 Log-likelihood function: -3.03823591619
    Iteration: 37 Log-likelihood function: -3.02646530573
    Iteration: 38 Log-likelihood function: -3.01520177394
    Iteration: 39 Log-likelihood function: -3.00441242601
    Iteration: 40 Log-likelihood function: -2.99406722296
    Iteration: 41 Log-likelihood function: -2.98413867259

一旦对数似然函数得分增加值小于0.01，迭代将终止。如果在500次迭代后还没有到达收敛， ``train`` 函数返回并抛出一个 ``AssertionError`` 。

可选的关键字 ``typecode`` 几乎可以一直忽略。这个关键字允许用户选择要使用的数值矩阵类型。当为了避免大数据计算的内存问题时，可能有必要使用单精度浮点数（Float8，Float16等等）而不是默认的double型。

16.1.3 使用logistic回归模型进行分类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

调用 ``classify`` 函数可以进行分类。给定一个logistic回归模型和 *x*\ :sub:`1` 和 *x*\ :sub:`2` 的值（例如，操纵子结构未知的基因对）， ``classify`` 函数返回 ``1`` 或 ``0`` ，分别对应类别OP和NOP。例如，考虑基因对 *yxcE* ， *yxcD* 和 *yxiB* ， *yxiA* ：

    --------------

    +-------------------------------------------------------------+
    | 表16.2：操纵子状态未知的相邻基因对   |
    +-------------------------------------------------------------+

    +-------------------+------------------------------------+---------------------------------------+
    | 基因对         | 基因间距离 *x*\ :sub:`1`   | 基因表达得分 *x*\ :sub:`2`   |
    +-------------------+------------------------------------+---------------------------------------+
    | *yxcE* — *yxcD*   | 6                                  | -173.143442352                        |
    +-------------------+------------------------------------+---------------------------------------+
    | *yxiB* — *yxiA*   | 309                                | -271.005880394                        |
    +-------------------+------------------------------------+---------------------------------------+

    --------------

Logistic回归模型预测 *yxcE* ， *yxcD* 属于相同操纵子（类别OP），而 *yxiB* ， *yxiA* 属于相同操纵子:

.. code:: verbatim

    >>> print "yxcE, yxcD:", LogisticRegression.classify(model, [6,-173.143442352])
    yxcE, yxcD: 1
    >>> print "yxiB, yxiA:", LogisticRegression.classify(model, [309, -271.005880394])
    yxiB, yxiA: 0

（这个结果和生物学文献报道的一致）。

为了确定这个预测的可信度，我们可以调用 ``calculate`` 函数来获得类别OP和NOP的概率(equations
(`16.2 <#eq:OP>`__) and `16.3 <#eq:NOP>`__)。对于*yxcE*, *yxcD*我们发现

.. code:: verbatim

    >>> q, p = LogisticRegression.calculate(model, [6,-173.143442352])
    >>> print "class OP: probability =", p, "class NOP: probability =", q
    class OP: probability = 0.993242163503 class NOP: probability = 0.00675783649744

对于 *yxiB* ， *yxiA*

.. code:: verbatim

    >>> q, p = LogisticRegression.calculate(model, [309, -271.005880394])
    >>> print "class OP: probability =", p, "class NOP: probability =", q
    class OP: probability = 0.000321211251817 class NOP: probability = 0.999678788748

为了确定回归模型的预测精确性，我们可以把模型应用到训练数据上：

.. code:: verbatim

    >>> for i in range(len(ys)):
            print "True:", ys[i], "Predicted:", LogisticRegression.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0

表示除一个基因对外所有预测都是正确的。Leave-one-out分析可以对预测精确性给出一个更可信的估计，这是通过从训练数据中移除要预测的基因，再重新计算模型实现：

.. code:: verbatim

    >>> for i in range(len(ys)):
            model = LogisticRegression.train(xs[:i]+xs[i+1:], ys[:i]+ys[i+1:])
            print "True:", ys[i], "Predicted:", LogisticRegression.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0

Leave-one-out分析显示这个logistic回归模型的预测只对两个基因对不正确，对应预测精确度为88%。

16.1.4 Logistic回归，线性判别分析和支持向量机
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Logistic回归模型类似于线性判别分析。在线性判别分析中，类别概率同样可由方程(`16.2 <#eq:OP>`__) and (`16.3 <#eq:NOP>`__)给出。但是，不是直接估计系数β，我们首先对预测变量 *x* 拟合一个正态分布。然后通过这个正态分布的平均值和方差计算系数β。如果 *x* 的分布确实是正态的，线性判别分析将比logistic回归模型有更好的性能。另一方面，logistic回归模型对于偏态到正态的广泛分布更加强健。

另一个相似的方法是应用线性核函数的支持向量机。这样的SVM也使用一个预测变量的线性组合，但是是从靠近类别之间的边界区域的预测变量 *x* 来估计系数β。如果logistic回归模型(equations (`16.2 <#eq:OP>`__) and (`16.3 <#eq:NOP>`__))很好的描述了远离边界区域的 *x* ，我们可以期望logistic回归模型优于线性核函数SVM，因为它应用了更多数据。如果不是，SVM可能更好。

16.2 *k*-最近邻居（ *KNN* ）
---------------------------

16.2.1 背景和目的
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

点是基于训练数据集的 *k* 个最近邻居类别进行分类的。

在Biopython中， *KNN* 方法可在 ``Bio.KNN`` 中获得。我们使用同样的操纵子数据集来说明Biopython中 *KNN* 方法的用法。

16.2.2 初始化一个 *KNN* 模型
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

使用表`16.1 <#table:training>`__中的数据，我们创建和初始化一个*KNN*模型：

.. code:: verbatim

    >>> from Bio import kNN
    >>> k = 3
    >>> model = kNN.train(xs, ys, k)

这里 ``xs`` 和 ``ys`` 和`16.1.2 <#subsec:LogisticRegressionTraining>`__中的相同。 ``k`` 是分类中的邻居数 *k* 。对于二分类，为 *k* 选择一个奇数可以避免tied votes。函数名 ``train`` 在这里有点不合适，因为就没有训练模型：这个函数仅仅是用来存储模型变量 ``xs`` ， ``ys`` 和 ``k`` 。

16.2.3 使用*KNN*模型来分类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

应用*KNN*模型对新数据进行分类，我们使用 ``classify`` 函数。这个函数以一个数据点(*x*\ :sub:`1`,\ *x*\ :sub:`2`)为参数并在训练数据集 ``xs`` 中寻找 *k* -最近邻居。然后基于在这*k*个邻居中出现最多的类别（ ``ys`` ）来对数据点(*x*\ :sub:`1`,\ *x*\ :sub:`2`)进行分类。

对于基因对 *yxcE* 、 *yxcD* 和 *yxiB* 、 *yxiA* 的例子，我们发现：

.. code:: verbatim

    >>> x = [6, -173.143442352]
    >>> print "yxcE, yxcD:", kNN.classify(model, x)
    yxcE, yxcD: 1
    >>> x = [309, -271.005880394]
    >>> print "yxiB, yxiA:", kNN.classify(model, x)
    yxiB, yxiA: 0

和logistic回归模型一致，*yxcE*,*yxcD*被归为一类（类别OP），*yxiB*,*yxiA*属于不同操纵子。

函数 ``classify`` 可以指定距离函数和权重函数作为可选参数。距离函数影响作为最近邻居的 *k* 个类别的选择，因为这些到查询点（ *x* ， *y* ）有最小距离的类别被定义为邻居。默认使用欧几里德距离。另外，我们也可以如示例中的使用曼哈顿距离：

.. code:: verbatim

    >>> def cityblock(x1, x2):
    ...    assert len(x1)==2
    ...    assert len(x2)==2
    ...    distance = abs(x1[0]-x2[0]) + abs(x1[1]-x2[1])
    ...    return distance
    ...
    >>> x = [6, -173.143442352]
    >>> print "yxcE, yxcD:", kNN.classify(model, x, distance_fn = cityblock)
    yxcE, yxcD: 1

权重函数可以用于权重投票。例如，相比于相邻较远的邻居，我们可能想给更近的邻居一个更高的权重：

.. code:: verbatim

    >>> def weight(x1, x2):
    ...    assert len(x1)==2
    ...    assert len(x2)==2
    ...    return exp(-abs(x1[0]-x2[0]) - abs(x1[1]-x2[1]))
    ...
    >>> x = [6, -173.143442352]
    >>> print "yxcE, yxcD:", kNN.classify(model, x, weight_fn = weight)
    yxcE, yxcD: 1

默认所有邻居有相同权重。

为了确定这些预测的置信度，我们可以调用函数 ``calculate`` 来计算分配到类别OP和NOP的总权重。对于默认的加权方案，这样减少了每个分类的邻居数量。对于 *yxcE* ， *yxcD* ， 我们发现

.. code:: verbatim

    >>> x = [6, -173.143442352]
    >>> weight = kNN.calculate(model, x)
    >>> print "class OP: weight =", weight[0], "class NOP: weight =", weight[1]
    class OP: weight = 0.0 class NOP: weight = 3.0

这意味着 ``x1`` ， ``x2`` 的所有三个邻居都属于NOP类别。对另一个例子 *yesK* ， *yesL* 我们发现

.. code:: verbatim

    >>> x = [117, -267.14]
    >>> weight = kNN.calculate(model, x)
    >>> print "class OP: weight =", weight[0], "class NOP: weight =", weight[1]
    class OP: weight = 2.0 class NOP: weight = 1.0

这意思是两个邻居是操纵子对，另一个是非操纵子对

对于*KNN*方法的预测精确性，我们对训练数据应用模型：

.. code:: verbatim

    >>> for i in range(len(ys)):
            print "True:", ys[i], "Predicted:", kNN.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0

显示除了两个基因对所有预测都是正确的。Leave-one-out分析可以对预测精确性给出一个更可信的估计，这是通过从训练数据中移除要预测的基因，再重新计算模型实现：

.. code:: verbatim

    >>> for i in range(len(ys)):
            model = kNN.train(xs[:i]+xs[i+1:], ys[:i]+ys[i+1:])
            print "True:", ys[i], "Predicted:", kNN.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 1

Leave-one-out分析显示这个 *KNN* 模型的预测正确17个基因对中的13个，对应预测精确度为76%。

16.3 Naive贝叶斯
-----------------

这部分将描述模块 ``Bio.NaiveBayes`` .

16.4 最大熵
---------------------

这部分将描述模块 ``Bio.MaximumEntropy``.

16.5  马尔科夫模型
-------------------

这部分将描述模块 ``Bio.MarkovModel`` 和/或
``Bio.HMM.MarkovModel`` .

