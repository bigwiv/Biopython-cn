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
到标准输出 to standard output, or an open file handle if
given. Not all of the available information about the tree is shown, but
it provides a way to quickly view the tree without relying on any
external dependencies.

.. code:: verbatim

    >>> tree = Phylo.read("example.xml", "phyloxml")
    >>> Phylo.draw_ascii(tree)
                 __________________ A
      __________|
    _|          |___________________________________________ B
     |
     |___________________________________________________________________________ C

The ``draw`` function draws a more attractive image using the matplotlib
library. See the API documentation for details on the arguments it
accepts to customize the output.

.. code:: verbatim

    >>> tree = Phylo.read("example.xml", "phyloxml")
    >>> Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

|image7|

``draw_graphviz`` draws an unrooted cladogram, but requires that you
have Graphviz, PyDot or PyGraphviz, NetworkX, and matplotlib (or pylab)
installed. Using the same example as above, and the ``dot`` program
included with Graphviz, let’s draw a rooted tree (see
Fig. `13.3 <#fig:phylo-dot>`__):

.. code:: verbatim

    >>> tree = Phylo.read("example.xml", "phyloxml")
    >>> Phylo.draw_graphviz(tree, prog='dot')
    >>> import pylab
    >>> pylab.show()                    # Displays the tree in an interactive viewer
    >>> pylab.savefig('phylo-dot.png')  # Creates a PNG file of the same graphic

|image8|

(Tip: If you execute IPython with the ``-pylab`` option, calling
``draw_graphviz`` causes the matplotlib viewer to launch automatically
without manually calling ``show()``.)

This exports the tree object to a NetworkX graph, uses Graphviz to lay
out the nodes, and displays it using matplotlib. There are a number of
keyword arguments that can modify the resulting diagram, including most
of those accepted by the NetworkX functions ``networkx.draw`` and
``networkx.draw_graphviz``.

The display is also affected by the ``rooted`` attribute of the given
tree object. Rooted trees are shown with a “head” on each branch
indicating direction (see Fig. `13.3 <#fig:phylo-rooted>`__):

.. code:: verbatim

    >>> tree = Phylo.read("simple.dnd", "newick")
    >>> tree.rooted = True
    >>> Phylo.draw_graphiz(tree)

|image9|

The “prog” argument specifies the Graphviz engine used for layout. The
default, ``twopi``, behaves well for any size tree, reliably avoiding
crossed branches. The ``neato`` program may draw more attractive
moderately-sized trees, but sometimes will cross branches (see
Fig. `13.3 <#fig:phylo-color>`__). The ``dot`` program may be useful
with small trees, but tends to do surprising things with the layout of
larger trees.

.. code:: verbatim

    >>> Phylo.draw_graphviz(tree, prog="neato")

|image10|

This viewing mode is particularly handy for exploring larger trees,
because the matplotlib viewer can zoom in on a selected region, thinning
out a cluttered graphic.

.. code:: verbatim

    >>> tree = Phylo.read("apaf.xml", "phyloxml")
    >>> Phylo.draw_graphviz(tree, prog="neato", node_size=0)

|image11| |image12|

Note that branch lengths are not displayed accurately, because Graphviz
ignores them when creating the node layouts. The branch lengths are
retained when exporting a tree as a NetworkX graph object
(``to_networkx``), however.

See the Phylo page on the Biopython wiki
(```http://biopython.org/wiki/Phylo`` <http://biopython.org/wiki/Phylo>`__)
for descriptions and examples of the more advanced functionality in
``draw_ascii``, ``draw_graphviz`` and ``to_networkx``.

13.4  Using Tree and Clade objects
----------------------------------

The ``Tree`` objects produced by ``parse`` and ``read`` are containers
for recursive sub-trees, attached to the ``Tree`` object at the ``root``
attribute (whether or not the phylogenic tree is actually considered
rooted). A ``Tree`` has globally applied information for the phylogeny,
such as rootedness, and a reference to a single ``Clade``; a ``Clade``
has node- and clade-specific information, such as branch length, and a
list of its own descendent ``Clade`` instances, attached at the
``clades`` attribute.

So there is a distinction between ``tree`` and ``tree.root``. In
practice, though, you rarely need to worry about it. To smooth over the
difference, both ``Tree`` and ``Clade`` inherit from ``TreeMixin``,
which contains the implementations for methods that would be commonly
used to search, inspect or modify a tree or any of its clades. This
means that almost all of the methods supported by ``tree`` are also
available on ``tree.root`` and any clade below it. (``Clade`` also has a
``root`` property, which returns the clade object itself.)

13.4.1  Search and traversal methods
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

13.4.2  Information methods
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

13.4.3  Modification methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These methods modify the tree in-place. If you want to keep the original
tree intact, make a complete copy of the tree first, using Python’s
``copy`` module:

.. code:: verbatim

    tree = Phylo.read('example.xml', 'phyloxml')
    import copy
    newtree = copy.deepcopy(tree)

 **``collapse``**
    Deletes the target from the tree, relinking its children to its
    parent.
**``collapse_all``**
    Collapse all the descendents of this tree, leaving only terminals.
    Branch lengths are preserved, i.e. the distance to each terminal
    stays the same. With a target specification (see above), collapses
    only the internal nodes matching the specification.
**``ladderize``**
    Sort clades in-place according to the number of terminal nodes.
    Deepest clades are placed last by default. Use ``reverse=True`` to
    sort clades deepest-to-shallowest.
**``prune``**
    Prunes a terminal clade from the tree. If taxon is from a
    bifurcation, the connecting node will be collapsed and its branch
    length added to remaining terminal node. This might no longer be a
    meaningful value.
**``root_with_outgroup``**
    Reroot this tree with the outgroup clade containing the given
    targets, i.e. the common ancestor of the outgroup. This method is
    only available on Tree objects, not Clades.

    If the outgroup is identical to self.root, no change occurs. If the
    outgroup clade is terminal (e.g. a single terminal node is given as
    the outgroup), a new bifurcating root clade is created with a
    0-length branch to the given outgroup. Otherwise, the internal node
    at the base of the outgroup becomes a trifurcating root for the
    whole tree. If the original root was bifurcating, it is dropped from
    the tree.

    In all cases, the total branch length of the tree stays the same.

**``root_at_midpoint``**
    Reroot this tree at the calculated midpoint between the two most
    distant tips of the tree. (This uses ``root_with_outgroup`` under
    the hood.)
**``split``**
    Generate *n* (default 2) new descendants. In a species tree, this is
    a speciation event. New clades have the given ``branch_length`` and
    the same name as this clade’s root plus an integer suffix (counting
    from 0) — for example, splitting a clade named “A” produces the
    sub-clades “A0” and “A1”.

See the Phylo page on the Biopython wiki
(```http://biopython.org/wiki/Phylo`` <http://biopython.org/wiki/Phylo>`__)
for more examples of using the available methods.

13.4.4  Features of PhyloXML trees
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The phyloXML file format includes fields for annotating trees with
additional data types and visual cues.

See the PhyloXML page on the Biopython wiki
(```http://biopython.org/wiki/PhyloXML`` <http://biopython.org/wiki/PhyloXML>`__)
for descriptions and examples of using the additional annotation
features provided by PhyloXML.

13.5  Running external applications
-----------------------------------

While Bio.Phylo doesn’t infer trees from alignments itself, there are
third-party programs available that do. These are supported through the
module ``Bio.Phylo.Applications``, using the same general framework as
``Bio.Emboss.Applications``, ``Bio.Align.Applications`` and others.

Biopython 1.58 introduced a wrapper for PhyML
(```http://www.atgc-montpellier.fr/phyml/`` <http://www.atgc-montpellier.fr/phyml/>`__).
The program accepts an input alignment in ``phylip-relaxed`` format
(that’s Phylip format, but without the 10-character limit on taxon
names) and a variety of options. A quick example:

.. code:: verbatim

    >>> from Bio import Phylo
    >>> from Bio.Phylo.Applications import PhymlCommandline
    >>> cmd = PhymlCommandline(input='Tests/Phylip/random.phy')
    >>> out_log, err_log = cmd()

This generates a tree file and a stats file with the names
[*input filename*\ ]\ ``_phyml_tree.txt`` and
[*input filename*\ ]\ ``_phyml_stats.txt``. The tree file is in Newick
format:

.. code:: verbatim

    >>> tree = Phylo.read('Tests/Phylip/random.phy_phyml_tree.txt', 'newick')
    >>> Phylo.draw_ascii(tree)

A similar wrapper for RAxML
(```http://sco.h-its.org/exelixis/software.html`` <http://sco.h-its.org/exelixis/software.html>`__)
was added in Biopython 1.60.

Note that some popular Phylip programs, including ``dnaml`` and
``protml``, are already available through the EMBOSS wrappers in
``Bio.Emboss.Applications`` if you have the Phylip extensions to EMBOSS
installed on your system. See Section \ `6.4 <#sec:alignment-tools>`__
for some examples and clues on how to use programs like these.

13.6  PAML integration
----------------------

Biopython 1.58 brought support for PAML
(```http://abacus.gene.ucl.ac.uk/software/paml.html`` <http://abacus.gene.ucl.ac.uk/software/paml.html>`__),
a suite of programs for phylogenetic analysis by maximum likelihood.
Currently the programs codeml, baseml and yn00 are implemented. Due to
PAML’s usage of control files rather than command line arguments to
control runtime options, usage of this wrapper strays from the format of
other application wrappers in Biopython.

A typical workflow would be to initialize a PAML object, specifying an
alignment file, a tree file, an output file and a working directory.
Next, runtime options are set via the ``set_options()`` method or by
reading an existing control file. Finally, the program is run via the
``run()`` method and the output file is automatically parsed to a
results dictionary.

Here is an example of typical usage of codeml:

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

Existing output files may be parsed as well using a module’s ``read()``
function:

.. code:: verbatim

    >>> results = codeml.read("Tests/PAML/Results/codeml/codeml_NSsites_all.out")
    >>> print results.get("lnL max")

Detailed documentation for this new module currently lives on the
Biopython wiki:
```http://biopython.org/wiki/PAML`` <http://biopython.org/wiki/PAML>`__

13.7  Future plans
------------------

Bio.Phylo is under active development. Here are some features we might
add in future releases:

 **New methods**
    Generally useful functions for operating on Tree or Clade objects
    appear on the Biopython wiki first, so that casual users can test
    them and decide if they’re useful before we add them to Bio.Phylo:

    ```http://biopython.org/wiki/Phylo_cookbook`` <http://biopython.org/wiki/Phylo_cookbook>`__

**Bio.Nexus port**
    Much of this module was written during Google Summer of Code 2009,
    under the auspices of NESCent, as a project to implement Python
    support for the phyloXML data format (see
    `13.4.4 <#sec:PhyloXML>`__). Support for Newick and Nexus formats
    was added by porting part of the existing Bio.Nexus module to the
    new classes used by Bio.Phylo.

    Currently, Bio.Nexus contains some useful features that have not yet
    been ported to Bio.Phylo classes — notably, calculating a consensus
    tree. If you find some functionality lacking in Bio.Phylo, try
    poking throught Bio.Nexus to see if it’s there instead.

We’re open to any suggestions for improving the functionality and
usability of this module; just let us know on the mailing list or our
bug database.
