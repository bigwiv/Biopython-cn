Chapter 13  Bio.Phylo系统发育分析
========================================
Biopython1.54开始引入了Bio.Phylo模块，与SeqIO和AlignIO类似，它的目的是提供
一个通用的独立于源数据格式的方法来使用系统进化树，同时提供一致的API来进行
I/O操作。

Bio.Phylo在一篇开放获取的期刊文章中有介绍
[`9 <#talevich2012>`__, Talevich *et al.*, 2012], 这可能对您也有所帮助。


13.1  示例: 树中有什么？ 
-----------------------------

为了熟悉这个模块，让我们首先从一个已经创建好的树开始，从几个不同的角度来审视
它。接着我们将给树的分支上颜色，并使用一个特殊的phyloXML特性，最终保存树。

在终端中使用你喜欢的编辑器创建一个简单的Newick文件：

.. code:: verbatim

    % cat > simple.dnd <<EOF
    > (((A,B),(C,D)),(E,F,G));
    > EOF

这棵树没有分支长度，只有一个拓扑结构和标记的端点。（如果你有一个真正的树文件，
你也可以使用它来替代进行示例操作。）

选择启动你的Python解释器：

.. code:: verbatim

    % ipython -pylab

对于交互式操作，使用参数``-pylab``启动IPython解释器能启用**matplotlib**整合
功能，这样图像就能自动弹出来。我们将在这个示例中使用该功能。

现在，在Python终端中，读取树文件，给定文件名和格式名。

.. code:: verbatim

    >>> from Bio import Phylo
    >>> tree = Phylo.read("simple.dnd", "newick")

以字符串打印该树对象我们将得到整个对象的层次结构概况。

.. code:: verbatim

    >>> print tree

    Tree(weight=1.0, rooted=False, name="")
        Clade(branch_length=1.0)
            Clade(branch_length=1.0)
                Clade(branch_length=1.0)
                    Clade(branch_length=1.0, name="A")
                    Clade(branch_length=1.0, name="B")
                Clade(branch_length=1.0)
                    Clade(branch_length=1.0, name="C")
                    Clade(branch_length=1.0, name="D")
            Clade(branch_length=1.0)
                Clade(branch_length=1.0, name="E")
                Clade(branch_length=1.0, name="F")
                Clade(branch_length=1.0, name="G")

``Tree`` 对象包含树的全局信息，如树是有根树还是无根树。它包含一个根进化枝，
和以此往下以列表嵌套的所有进化枝，直至叶子分支。

函数 ``draw_ascii`` 创建一个简单的ASCII-art(纯文本)系统发生图。在没有更好
图形工具的情况下，这对于交互研究来说是一个方便的可视化展示方式。

.. code:: verbatim

    >>> Phylo.draw_ascii(tree)
                                                        ________________________ A
                               ________________________|
                              |                        |________________________ B
      ________________________|
     |                        |                         ________________________ C
     |                        |________________________|
    _|                                                 |________________________ D
     |
     |                         ________________________ E
     |                        |
     |________________________|________________________ F
                              |
                              |________________________ G

如果你安装有 **matplotlib** 或者 **pylab**, 你可以使用 ``draw`` 函数一个图像(见 Fig.
`13.1 <#fig:phylo-simple-draw>`__):

.. code:: verbatim

    >>> tree.rooted = True
    >>> Phylo.draw(tree)

|image5|

13.1.1  给树的分支上颜色
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
函数 ``draw`` 和 ``draw_graphviz`` 支持在树中显示不同的颜色和分支宽度。
从Biopython 1.59开始，Clade对象就开始支持 ``color`` 和 ``width`` 属性，
且使用他们不需要额外支持。这两个属性都表示导向给定的进化枝前面的分支的
属性，并依次往下作用，所以所有的后代分支在显示时也都继承相同的宽度和颜
色。

在早期的Biopython版本中，PhyloXML树有些特殊的特性，使用这些属性需要首先
将这个树转换为一个基本树对象的子类Phylogeny，该类在Bio.Phylo.PhyloXML模
块中。

在Biopython 1.55和之后的版本中，这是一个很方便的树方法：

.. code:: verbatim

    >>> tree = tree.as_phyloxml()

在Biopython 1.54中, 你能通过导入一个额外的模块实现相同的事情：

.. code:: verbatim

    >>> from Bio.Phylo.PhyloXML import Phylogeny
    >>> tree = Phylogeny.from_tree(tree)

注意Newick和Nexus文件类型并不支持分支颜色和宽度，如果你在Bio.Phylo中使用
这些属性，你只能保存这些值到PhyloXML格式中。（你也可以保存成Newick或Nexus
格式，但是颜色和宽度信息在输出的文件时会被忽略掉。）

现在我们开始赋值颜色。首先，我们将设置根进化枝为灰色。我们能通过赋值24位
的颜色值来实现，用三位数的RGB值、HTML格式的十六进制字符串、或者预先设置好的
颜色名称。

.. code:: verbatim

    >>> tree.root.color = (128, 128, 128)

Or:

.. code:: verbatim

    >>> tree.root.color = "#808080"

Or:

.. code:: verbatim

    >>> tree.root.color = "gray"

一个进化枝的颜色会被当作从上而下整个进化枝的颜色，所以我们这里设置根的
的颜色会将整个树的颜色变为灰色。我们能通过在树中下面分支赋值不同的颜色
来重新定义某个分支的颜色。

让我们先定位“E”和“F”最近祖先（MRCA）节点。方法 ``common_ancestor`` 返回
原始树中这个进化枝的引用，所以当我们设置该进化枝为“salmon”颜色时，这个颜
色则会在原始的树中显示出来。

.. code:: verbatim

    >>> mrca = tree.common_ancestor({"name": "E"}, {"name": "F"})
    >>> mrca.color = "salmon"

当我没碰巧明确地知道某个进化枝在树中的位置，以嵌套列表的形式，我们就能
通过索引的方式直接跳到那个位置。这里，索引 ``[0,1]`` 表示根节点的第一个
子代节点的第二个子代。

.. code:: verbatim

    >>> tree.clade[0,1].color = "blue"

最后，展示一下我们的工作结果 (see Fig. `13.1.1 <#fig:phylo-color-draw>`__):

.. code:: verbatim

    >>> Phylo.draw(tree)

|image6|

注意进化枝的颜色包括导向它的分支和它的子代的分支。E和F的共同祖先结果刚好
在根分支下面，而通过这样上色，我们能清楚的看出这个树的根在哪里。

我们已经完成了很多！现在让我们休息一下，保存一下我们的工作。使用一个文件
名或句柄（这里我们使用标准输出来查看将会输出什么）和 ``phyloxml`` 格式来
调用 ``write`` 函数。PhyloXML格式保存了我们设置的颜色，所以你能通过其他树
查看工具，如Archaeopteryx，打开这个phyloXML文件，这些颜色也会显示出来。

.. code:: verbatim

    >>> import sys
    >>> Phylo.write(tree, sys.stdout, "phyloxml")

    <phy:phyloxml xmlns:phy="http://www.phyloxml.org">
      <phy:phylogeny rooted="true">
        <phy:clade>
          <phy:branch_length>1.0</phy:branch_length>
          <phy:color>
            <phy:red>128</phy:red>
            <phy:green>128</phy:green>
            <phy:blue>128</phy:blue>
          </phy:color>
          <phy:clade>
            <phy:branch_length>1.0</phy:branch_length>
            <phy:clade>
              <phy:branch_length>1.0</phy:branch_length>
              <phy:clade>
                <phy:name>A</phy:name>
                ...

本章的其余部分将更加细致的介绍Bio.Phylo核心功能。关于Bio.Phylo的更多例
子，请参见Biopython.org上的Cookbook手册页面。

```http://biopython.org/wiki/Phylo_cookbook`` <http://biopython.org/wiki/Phylo_cookbook>`__

13.2  I/O 函数
-------------------

和SeqIO、AlignIO类似, Phylo使用四个函数处理文件的输入输出： ``parse`` 、
``read`` 、 ``write`` 和 ``convert`` ，所有的函数都支持Newick、NEXUS、
phyloXML和NeXML等树文件格式。

``read`` 函数解析并返回给定文件中的单个树。注意，如果文件中包含多个或不包含任何树，它将抛出一个错误。

.. code:: verbatim

    >>> from Bio import Phylo
    >>> tree = Phylo.read("Tests/Nexus/int_node_labels.nwk", "newick")
    >>> print tree

（Biopython发布包的 ``Tests/Nexus/`` 和 ``Tests/PhyloXML/`` 文件夹中有相应的例子）

处理多个（或者未知个数）的树文件，需要使用 ``parse`` 函数迭代给定文件中的每一个树。

.. code:: verbatim

    >>> trees = Phylo.parse("Tests/PhyloXML/phyloxml_examples.xml", "phyloxml")
    >>> for tree in trees:
    ...     print tree

使用 ``write‵‵ 函数输出一个或多个可迭代的树。

.. code:: verbatim

    >>> trees = list(Phylo.parse("phyloxml_examples.xml", "phyloxml"))
    >>> tree1 = trees[0]
    >>> others = trees[1:]
    >>> Phylo.write(tree1, "tree1.xml", "phyloxml")
    1
    >>> Phylo.write(others, "other_trees.xml", "phyloxml")
    12

使用 ``convert`` 函数转换任何支持的树格式。

.. code:: verbatim

    >>> Phylo.convert("tree1.dnd", "newick", "tree1.xml", "nexml")
    1
    >>> Phylo.convert("other_trees.xml", "phyloxml", "other_trees.nex", 'nexus")
    12

和SeqIO和AlignIO类似，当使用字符串而不是文件作为输入输出时，需要使用 ‵‵StringIO`` 函数。

.. code:: verbatim

    >>> from Bio import Phylo
    >>> from StringIO import StringIO
    >>> handle = StringIO("(((A,B),(C,D)),(E,F,G));")
    >>> tree = Phylo.read(handle, "newick")

13.3  查看和导出树
---------------------------

了解一个 ``Tree`` 对象概况的最简单的方法是用 ``print`` 函数将它打印出来：

.. code:: verbatim

    >>> tree = Phylo.read("Tests/PhyloXML/example.xml", "phyloxml")
    >>> print tree
    Phylogeny(rooted='True', description='phyloXML allows to use either a "branch_length"
    attribute...', name='example from Prof. Joe Felsenstein's book "Inferring Phyl...')
        Clade()
            Clade(branch_length='0.06')
                Clade(branch_length='0.102', name='A')
                Clade(branch_length='0.23', name='B')
            Clade(branch_length='0.4', name='C')

上面实际上是Biopython的树对象层次结构的一个概况。然而更可能的情况是，你希望见到
画出树的形状，这里有三个函数来做这件事情。

如我们在demo中看到的一样， ``draw_ascii`` 打印一个树的ascii-art图像（有根进化树）
到标准输出，或者一个打开的文件句柄，若有提供。不是所有关于树的信息被显示出来，但是它提供了一个
不依靠于任何外部依赖的快速查看树的方法。

.. code:: verbatim

    >>> tree = Phylo.read("example.xml", "phyloxml")
    >>> Phylo.draw_ascii(tree)
                 __________________ A
      __________|
    _|          |___________________________________________ B
     |
     |___________________________________________________________________________ C

``draw`` 函数则使用matplotlib类库画出一个更加好看的图像。查看API文档以获得关于它所接受的
用来定制输出的参数。

.. code:: verbatim

    >>> tree = Phylo.read("example.xml", "phyloxml")
    >>> Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

|image7|

``draw_graphviz`` 则画出一个无根的进化分枝图（cladogram），但是它要求你安装有Graphviz、
PyDot或PyGraphviz、Network和matplotlib（或pylab）。使用上面相同的例子，和Graphviz中的
 ``dot`` 程序，让我们来画一个有根树（见图. `13.3 <#fig:phylo-dot>`__ ）：

.. code:: verbatim

    >>> tree = Phylo.read("example.xml", "phyloxml")
    >>> Phylo.draw_graphviz(tree, prog='dot')
    >>> import pylab
    >>> pylab.show()                    # Displays the tree in an interactive viewer
    >>> pylab.savefig('phylo-dot.png')  # Creates a PNG file of the same graphic

|image8|

（提示：如果你使用 ``-pylab`` 选项执行IPython，调用 ``draw_graphviz`` 将导致matplotlib
查看器自动运行，而不需要手动的调用 ``show()`` 方法。）

这将输出树对象到一个NetworkX图中，使用Graphviz来布局节点的位置，并使用matplotlib来显示
它。这里有几个关键词参数来修改结果图像，包括大多数被NetworkX函数 ``networkx.draw`` 和
``networkx.draw_graphviz`` 所接受的参数。

最终的显示也受所提供的树对象的 ``rooted`` 属性的影响。有根树在每个分支（branch）上显示
一个“head”来表明它的方向（见图. `13.3 <#fig:phylo-rooted>`__ ）：

.. code:: verbatim

    >>> tree = Phylo.read("simple.dnd", "newick")
    >>> tree.rooted = True
    >>> Phylo.draw_graphiz(tree)

|image9|

“prog”参数指定Graphviz的用来布局的引擎。默认的引擎 ``twopi`` 对任何大小的树都表现很好，
很可靠的避免交叉的分支出现。``neato``程序可能画出更加好看的中等大小的树，但是有时候会
有交叉分支出现（见图. `13.3 <#fig:phylo-color>`__ ）。 ``dot`` 程序或许对小型的树有用，
但是对于大一点的树的布局易产生奇怪的事情。

.. code:: verbatim

    >>> Phylo.draw_graphviz(tree, prog="neato")

|image10|

这个查看方式非常方便研究大型的树，因为matplotlib查看器可以放大选择的区域，使得杂乱的图像
变得稀疏。

.. code:: verbatim

    >>> tree = Phylo.read("apaf.xml", "phyloxml")
    >>> Phylo.draw_graphviz(tree, prog="neato", node_size=0)

|image11| |image12|

注意，分支长度并没有被正确地显示，因为Graphviz在布局时忽略了他们。然而，分支长度可以在输出
树为NetworkX图对象（ ``to_networkx`` ）时重新获得。

查看Biopython维基的Phylo页面
(```http://biopython.org/wiki/Phylo`` <http://biopython.org/wiki/Phylo>`__)
以获得关于 ``draw_ascii`` 、 ``draw_graphviz`` 和 ``to_networkx`` 的更加高级的功能的描述
和例子。

13.4  使用Tree和Clade对象
----------------------------------

``parse`` 和 ``read`` 方法产生的 ``Tree`` 对象是一些包含递归的子树的容器，连接到 ``Tree``
对象的 ``root`` 属性（不管进化树实际上被认为是否有根）。一个 ``Tree`` 包含进化树的全局信息，
如有根性（rootedness）和指向一个单独的 ``Clade`` 的引用; 一个 ``Clade`` 包含节点和进化枝
特异性信息，如分支长度（branch length）和一个它自身后代 ``Clade`` 实例的列表，附着在 ``clades``
属性上。

所以，这里 ``tree`` 和 ``tree.root`` 间是有区别的. 然而，实际操作中，你几乎不需要担心它。为了
缓和这个不同，``Tree`` 和 ``Clade`` 两者都继承自 ``TreeMixin``，它包含常用的用来查找、审视和
修改树和任何它的进化枝的方法的实现。这意味着，所有 ``tree`` 所支持的方法在 ``tree.root`` 和
任何它下面的clade中都能用。（ ``Clade`` 也有一个 ``root`` 属性，它返回clade对象本身。）

13.4.1  查找和遍历类方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For convenience, we provide a couple of simplified methods that return
all external or internal nodes directly as a list:

 **``get_terminals``**
    makes a list of all of this tree’s terminal (leaf) nodes.
**``get_nonterminals``**
    makes a list of all of this tree’s nonterminal (internal) nodes.

These both wrap a method with full control over tree traversal,
``find_clades``. Two more traversal methods, ``find_elements`` and
``find_any``, rely on the same core functionality and accept the same
arguments, which we’ll call a “target specification” for lack of a
better description. These specify which objects in the tree will be
matched and returned during iteration. The first argument can be any of
the following types:

-  A **TreeElement instance**, which tree elements will match by
   identity — so searching with a Clade instance as the target will find
   that clade in the tree;
-  A **string**, which matches tree elements’ string representation — in
   particular, a clade’s ``name`` *(added in Biopython 1.56)*;
-  A **class** or **type**, where every tree element of the same type
   (or sub-type) will be matched;
-  A **dictionary** where keys are tree element attributes and values
   are matched to the corresponding attribute of each tree element. This
   one gets even more elaborate:

   -  If an ``int`` is given, it matches numerically equal attributes,
      e.g. 1 will match 1 or 1.0
   -  If a boolean is given (True or False), the corresponding attribute
      value is evaluated as a boolean and checked for the same
   -  ``None`` matches ``None``
   -  If a string is given, the value is treated as a regular expression
      (which must match the whole string in the corresponding element
      attribute, not just a prefix). A given string without special
      regex characters will match string attributes exactly, so if you
      don’t use regexes, don’t worry about it. For example, in a tree
      with clade names Foo1, Foo2 and Foo3,
      ``tree.find_clades({"name": "Foo1"})`` matches Foo1,
      ``{"name": "Foo.*"}`` matches all three clades, and
      ``{"name": "Foo"}`` doesn’t match anything.

   Since floating-point arithmetic can produce some strange behavior, we
   don’t support matching ``float``\ s directly. Instead, use the
   boolean ``True`` to match every element with a nonzero value in the
   specified attribute, then filter on that attribute manually with an
   inequality (or exact number, if you like living dangerously).

   If the dictionary contains multiple entries, a matching element must
   match each of the given attribute values — think “and”, not “or”.

-  A **function** taking a single argument (it will be applied to each
   element in the tree), returning True or False. For convenience,
   LookupError, AttributeError and ValueError are silenced, so this
   provides another safe way to search for floating-point values in the
   tree, or some more complex characteristic.

After the target, there are two optional keyword arguments:

 **terminal**
    — A boolean value to select for or against terminal clades (a.k.a.
    leaf nodes): True searches for only terminal clades, False for
    non-terminal (internal) clades, and the default, None, searches both
    terminal and non-terminal clades, as well as any tree elements
    lacking the ``is_terminal`` method.
**order**
    — Tree traversal order: ``"preorder"`` (default) is depth-first
    search, ``"postorder"`` is DFS with child nodes preceding parents,
    and ``"level"`` is breadth-first search.

Finally, the methods accept arbitrary keyword arguments which are
treated the same way as a dictionary target specification: keys indicate
the name of the element attribute to search for, and the argument value
(string, integer, None or boolean) is compared to the value of each
attribute found. If no keyword arguments are given, then any TreeElement
types are matched. The code for this is generally shorter than passing a
dictionary as the target specification:
``tree.find_clades({"name": "Foo1"})`` can be shortened to
``tree.find_clades(name="Foo1")``.

(In Biopython 1.56 or later, this can be even shorter:
``tree.find_clades("Foo1")``)

Now that we’ve mastered target specifications, here are the methods used
to traverse a tree:

 **``find_clades``**
    Find each clade containing a matching element. That is, find each
    element as with ``find_elements``, but return the corresponding
    clade object. (This is usually what you want.)

    The result is an iterable through all matching objects, searching
    depth-first by default. This is not necessarily the same order as
    the elements appear in the Newick, Nexus or XML source file!

**``find_elements``**
    Find all tree elements matching the given attributes, and return the
    matching elements themselves. Simple Newick trees don’t have complex
    sub-elements, so this behaves the same as ``find_clades`` on them.
    PhyloXML trees often do have complex objects attached to clades, so
    this method is useful for extracting those.
**``find_any``**
    Return the first element found by ``find_elements()``, or None. This
    is also useful for checking whether any matching element exists in
    the tree, and can be used in a conditional.

Two more methods help navigating between nodes in the tree:

 **``get_path``**
    List the clades directly between the tree root (or current clade)
    and the given target. Returns a list of all clade objects along this
    path, ending with the given target, but excluding the root clade.
**``trace``**
    List of all clade object between two targets in this tree. Excluding
    start, including finish.

13.4.2  信息类方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~

These methods provide information about the whole tree (or any clade).

 **``common_ancestor``**
    Find the most recent common ancestor of all the given targets. (This
    will be a Clade object). If no target is given, returns the root of
    the current clade (the one this method is called from); if 1 target
    is given, this returns the target itself. However, if any of the
    specified targets are not found in the current tree (or clade), an
    exception is raised.
**``count_terminals``**
    Counts the number of terminal (leaf) nodes within the tree.
**``depths``**
    Create a mapping of tree clades to depths. The result is a
    dictionary where the keys are all of the Clade instances in the
    tree, and the values are the distance from the root to each clade
    (including terminals). By default the distance is the cumulative
    branch length leading to the clade, but with the
    ``unit_branch_lengths=True`` option, only the number of branches
    (levels in the tree) is counted.
**``distance``**
    Calculate the sum of the branch lengths between two targets. If only
    one target is specified, the other is the root of this tree.
**``total_branch_length``**
    Calculate the sum of all the branch lengths in this tree. This is
    usually just called the “length” of the tree in phylogenetics, but
    we use a more explicit name to avoid confusion with Python
    terminology.

The rest of these methods are boolean checks:

 **``is_bifurcating``**
    True if the tree is strictly bifurcating; i.e. all nodes have either
    2 or 0 children (internal or external, respectively). The root may
    have 3 descendents and still be considered part of a bifurcating
    tree.
**``is_monophyletic``**
    Test if all of the given targets comprise a complete subclade —
    i.e., there exists a clade such that its terminals are the same set
    as the given targets. The targets should be terminals of the tree.
    For convenience, this method returns the common ancestor (MCRA) of
    the targets if they are monophyletic (instead of the value
    ``True``), and ``False`` otherwise.
**``is_parent_of``**
    True if target is a descendent of this tree — not required to be a
    direct descendent. To check direct descendents of a clade, simply
    use list membership testing: ``if subclade in clade: ...``
**``is_preterminal``**
    True if all direct descendents are terminal; False if any direct
    descendent is not terminal.

13.4.3  修改类方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

这些方法都在原地对树进行修改，所以如果你想保持原来的树不变，你首先要使用Python的
``copy`` 模块对树进行完整的拷贝：

.. code:: verbatim

    tree = Phylo.read('example.xml', 'phyloxml')
    import copy
    newtree = copy.deepcopy(tree)

 **``collapse``**
    从树中删除目标，重新连接它的子代（children）到它的父亲节点（parent）。
**``collapse_all``**
    删除这个树的所有后代（descendents），只保留末端节点（terminals）。
    分支长度被保留，即到每个末端节点的距离保持不变。如指定一个目标（见上），
    只坍塌（collapses）和指定匹配的内部节点。
**``ladderize``**
    根据末端节点的个数，在原地对进化枝（clades）进行排序。越深的进化枝默认被放到最后，
    使用 ``reverse=True`` 将其放到最前。
**``prune``**
    从树中修剪末端进化枝（terminal clade）。如果分类名（taxon）来自一个二叉枝（bifurcation），
    连接的节点将被坍塌，它的分支长度将被加到剩下的末端节点上。这可能不再是一个有意义的值。
**``root_with_outgroup``**
    使用包含给定目标的外群进化枝（outgroup clade）重新确定树的根节点，即外群的共同祖先。该方法
    只在Tree对象中能用，不能用于Clade对象。

    如果外群和self.root一致，将不发生改变。如果外群进化枝是末端（即一个末端节点被作为外群），一个
    新的二叉根进化枝将被创建，且到给定外群的分支长度为0。否则，外群根部的内部节点变为整个树的一个
    三叉根。如果原先的根是一个二叉，它将被从树中遗弃。

    在所有的情况下，树的分支长度总和保持不变。

**``root_at_midpoint``**
    重新选择树中两个最远的节点的中点作为树的根。（这实际上是使用 ``root_with_outgroup`` 函数。）
**``split``**
    产生 *n* （默认为2）个 新的后代。在一个物种树中，这是一个物种形成事件。新的进化枝拥有给定的
    ``branch_length`` 以及和这个进化枝的根相同的名字，名字后面包含一个整数后缀（从0开始计数）——
    例如，分割名为“A”的进化枝将生成子进化枝“A0”和“A1”。

查看Biopython维基的Phylo页面
(```http://biopython.org/wiki/Phylo`` <http://biopython.org/wiki/Phylo>`__)
以获得更多已有方法的使用示例。

13.4.4  PhyloXML树的特性
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phyloXML文件格式包含用来注释树的，采用额外数据格式和图像提示的字段。

参加Biopython维基上的PhyloXML页面
(```http://biopython.org/wiki/PhyloXML`` <http://biopython.org/wiki/PhyloXML>`__)
以查看关于使用PhyloXML提供的额外注释特性的描述和例子。

13.5  运行外部程序
-----------------------------------

尽管Bio.Phylo本身不从序列比对推断进化树，但这里有一些第三方的程序可以使用。
他们通过 ``Bio.Phylo.Applications`` 模块获得支持，使用和 ``Bio.Emboss.Applications`` 、
 ``Bio.Align.Applications`` 以及其他模块相同的通用框架。

Biopython 1.58引入了一个PhyML的打包程序（wrapper）
(```http://www.atgc-montpellier.fr/phyml/`` <http://www.atgc-montpellier.fr/phyml/>`__)。
该程序接受一个 ``phylip-relaxed`` 格式（它是Phylip格式，然而没有对分类名称的10个字符的限制）
的比对输入和多种参数。一个快速的例子是：

.. code:: verbatim

    >>> from Bio import Phylo
    >>> from Bio.Phylo.Applications import PhymlCommandline
    >>> cmd = PhymlCommandline(input='Tests/Phylip/random.phy')
    >>> out_log, err_log = cmd()

这生成一个树文件盒一个统计文件，名称为：
[*input filename*\ ]\ ``_phyml_tree.txt`` 和
[*input filename*\ ]\ ``_phyml_stats.txt``. 树文件的格式是Newick格式：

.. code:: verbatim

    >>> tree = Phylo.read('Tests/Phylip/random.phy_phyml_tree.txt', 'newick')
    >>> Phylo.draw_ascii(tree)

一个类似的RAxML打包程序
(```http://sco.h-its.org/exelixis/software.html`` <http://sco.h-its.org/exelixis/software.html>`__)
也已经被添加到Biopython 1.60中。

注意，如果你系统中已经安装了EMBOSS的Phylip扩展，一些常用的Phylip程序，包括 ``dnaml`` 和 ``protml`` 
已经通过 ``Bio.Emboss.Applications`` 中的EMBOSS打包程序被支持。参见章节 \ `6.4 <#sec:alignment-tools>`__
以查看使用这些程序的例子和提示。

13.6  PAML整合
----------------------

Biopython 1.58引入了对PAML的支持
(```http://abacus.gene.ucl.ac.uk/software/paml.html`` <http://abacus.gene.ucl.ac.uk/software/paml.html>`__),
它是一个采用最大似然法（maximum likelihood）进行系统进化分析的程序包。目前，对程序codeml、baseml和yn00的支持
已经实现。由于PAML使用控制文件而不是命令行参数来控制运行时选项，这个打包程序（wrapper）的使用格式和Biopython
的其他应用打包程序有些差异。

一个典型的流程是：初始化一个PAML对象，指定一个比对文件，一个树文件，一个输出文件和工作路径。下一步，运行时
选项通过 ``set_options()`` 方法或者读入一个已有的控制文件来设定。最后，程序通过 ``run()`` 方法来运行，输出文件
将自动被解析到一个结果目录。


下面是一个codeml典型用法的例子：

.. code:: verbatim

    >>> from Bio.Phylo.PAML import codeml
    >>> cml = codeml.Codeml()
    >>> cml.alignment = "Tests/PAML/alignment.phylip"
    >>> cml.tree = "Tests/PAML/species.tree"
    >>> cml.out_file = "results.out"
    >>> cml.working_dir = "./scratch"
    >>> cml.set_options(seqtype=1,
    ...         verbose=0,
    ...         noisy=0,
    ...         RateAncestor=0,
    ...         model=0,
    ...         NSsites=[0, 1, 2],
    ...         CodonFreq=2,
    ...         cleandata=1,
    ...         fix_alpha=1,
    ...         kappa=4.54006)
    >>> results = cml.run()
    >>> ns_sites = results.get("NSsites")
    >>> m0 = ns_sites.get(0)
    >>> m0_params = m0.get("parameters")
    >>> print m0_params.get("omega")

已有的输出文件也可以通过模块的 ``read()`` 方法来解析：

.. code:: verbatim

    >>> results = codeml.read("Tests/PAML/Results/codeml/codeml_NSsites_all.out")
    >>> print results.get("lnL max")

这个新模块的详细介绍目前在Biopython维基上可以看到：
```http://biopython.org/wiki/PAML`` <http://biopython.org/wiki/PAML>`__

13.7  未来计划
------------------

Bio.Phylo 目前还在开发中，下面是我们可能会在将来的发布版本中添加的特性：

 **新方法**
    通常用来操作Tree和Clade对象的有用方法会首先出现在Biopython维基上，这样常规用户
    就能在我们添加到Bio.Phylo之前测试这些方法，看看它们是否有用：
    ```http://biopython.org/wiki/Phylo_cookbook`` <http://biopython.org/wiki/Phylo_cookbook>`__

**Bio.Nexus port**
    这个模块的大部分是在2009年NESCent主办的谷歌编程夏令营中写的，作为实现Python对phyloXML数据格式（见
    `13.4.4 <#sec:PhyloXML>`__ ）支持的一个项目。对Newick和Nexus格式的支持，已经通过导入Bio.Nexus模块
    的一部分被添加到Bio.Phylo使用的新类中。

    目前，Bio.Nexus包含一些还没有导入到Bio.Phylo类中的有用的特性——特别是，计算一致树（consensus tree）。
    如果你发现某些功能Bio.Phylo中没有，试试在Bio.Nexus中能不能找到。

我们乐意接受任何增强该模块功能和使用性的建议；如果有，只需要通过邮件列表或我们的bug数据库让我们知道。
