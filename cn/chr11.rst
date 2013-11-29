第11章  走向3D：PDB模块
=========================================================

Bio.PDB是Biopython中处理生物大分子晶体结构的模块。除了别的类之外，Bio.PDB包含PDBParser类，此类能够产生一个Structure对象，以一种较方便的方式获取文件中的原子数据。只是在处理PDB文件头所包含的信息时，该类有一定的局限性。


11.1  晶体结构文件的读与写 
----------------------

11.1.1  读取PDB文件 
~~~~~~~~~~~~~~~~~~~~~

首先，我们创建一个 ``PDBParser`` 对象：


    >>> from Bio.PDB.PDBParser import PDBParser
    >>> p = PDBParser(PERMISSIVE=1)

``PERMISSIV`` 标签表示一些与PDB文件相关的问题（见 `11.7.1 <#problem%20structures>`__ ）会被忽略（注意某些原子和/或残基会丢失）。如果没有这个标签，则会在解析器运行期间有问题被检测到的时候生成一个 ``PDBConstructionException`` 标签。


接着通过 ``PDBParser`` 解析PDB文件，就产生了Structure对象（在此例子中，PDB文件为'pdb1fat.ent'，'1fat'是用户定义的结构名称）:


    >>> structure_id = "1fat"
    >>> filename = "pdb1fat.ent"
    >>> s = p.get_structure(structure_id, filename)

你可以从PDBParser对象中用 ``get_header`` 和 ``get_trailer`` 方法来提取PDB文件中的文件头和文件尾（简单的字符串列表）。然而许多PDB文件头包含不完整或错误的信息。许多错误在等价的mmCIF格式文件中得到修正。* 因此，如果你对文件头信息感兴趣，可以用下面即将讲到的 ``MMCIF2Dict`` 来提取信息，而不用处理PDB文件文件头。* 


现在澄清了，让我们回到解析PDB文件头这件事上。结构对象有个属性叫 ``header`` ，这是一个将头记录映射到其相应值的Python字典。


例子：

    >>> resolution = structure.header['resolution']
    >>> keywords = structure.header['keywords']

在这个字典中可用的关键字有 ``name`` 、 ``head`` 、 ``deposition_date`` 、 ``release_date`` 、 ``structure_method`` 、 ``resolution`` 、 ``structure_reference`` （映射到一个参考文献列表）、 ``journal_reference`` 、 ``author`` 、和 ``compound`` （映射到一个字典，其中包含结晶化合物的各种信息）。

没有创建 ``Structure`` 对象的时候，也可以创建这个字典，比如直接从PDB文件创建:

    >>> file = open(filename,'r')
    >>> header_dict = parse_pdb_header(file)
    >>> file.close()

11.1.2  读取mmCIF文件 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

与PDB文件的情形类似，先创建一个 ``MMCIFParser`` 对象：

    >>> from Bio.PDB.MMCIFParser import MMCIFParser
    >>> parser = MMCIFParser()

然后用这个解析器从mmCIF文件创建一个结构对象：

    >>> structure = parser.get_structure('1fat', '1fat.cif')

为了尽量少访问mmCIF文件，可以用 ``MMCIF2Dict`` 类创建一个Python字典来将所有mmCIF文件中各种标签映射到其对应的值上。若有多个值（像 ``_atom_site.Cartn_y`` 标签，储存的是所有原子的*y*坐标值），则这个标签映射到一个值列表。从mmCIF文件创建字典如下：

    >>> from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    >>> mmcif_dict = MMCIF2Dict('1FAT.cif')

例：从mmCIF文件获取溶剂含量:

    >>> sc = mmcif_dict['_exptl_crystal.density_percent_sol']

例：获取包含所有原子*y*坐标的列表:

    >>> y_list = mmcif_dict['_atom_site.Cartn_y']

11.1.3  读取PDB XML格式的文件
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

这个功能暂时还不支持，不过我们确实计划在未来支持这个功能（这项任务并不大）。如果你需要的话联系Biopython开发人员（ `biopython-dev@biopython.org <mailto:biopython-dev@biopython.org>`__ ）。

11.1.4  写PDB文件
~~~~~~~~~~~~~~~~~~~~~~~~~

可以用PDBIO类实现。当然也可很方便地输出一个结构的特定部分。

例子：保存一个结构


    >>> io = PDBIO()
    >>> io.set_structure(s)
    >>> io.save('out.pdb')

如果你想写出结构的一部分，可以用 `Select` 类（也在 ``PDBIO`` 中）来实现。 `Select` 有如下四种方法：

-  ``accept_model(model)``
-  ``accept_chain(chain)``
-  ``accept_residue(residue)``
-  ``accept_atom(atom)``

在默认情况下，每种方法的返回值都为1（表示model/chain/residue/atom被包含在输出结果中）。通过子类化 ``Select`` 和返回值0，你可以从输出中排除model、chain等。也许麻烦，但很强大。接下来的代码将只输出甘氨酸残基：



    >>> class GlySelect(Select):
    ...     def accept_residue(self, residue):
    ...         if residue.get_name()=='GLY':
    ...             return True
    ...         else:
    ...             return False
    ...
    >>> io = PDBIO()
    >>> io.set_structure(s)
    >>> io.save('gly_only.pdb', GlySelect())

如果这部分对你来说太复杂，那么 ``Dice`` 模块有一个很方便的 ``extract`` 函数，它可以输出一条链中起始和终止氨基酸残基之间的所有氨基酸残基。

11.2  结构的表示 
-------------------------------------------

一个 ``Structure`` 对象的整体布局遵循称为SMCRA（Structure/Model/Chain/Residue/Atom，结构/模型/链/残基/原子）的体系架构：
 -  结构由模型组成
 -  模型由多条链组成
 -  链由残基组成
 -  多个原子构成残基

这是很多结构生物学家/生物信息学家看待结构的方法，也是处理结构的一种简单而有效的方法。在需要的时候加上额外的材料。一个 ``Structure`` 对象的UML图（暂时忘掉 ``Disordered``吧）如下图所示 `11.1 <#fig:smcra>`__ 。这样的数据结构不一定最适用于表示一个结构的生物大分子内容，但要很好地解释一个描述结构的文件中所呈现的数据（最典型的如PDB或MMCIF文件），这样的数据结构就是必要的了。如果这种层次结构不能表示一个结构文件的内容，那么可以相当确定是这个文件有错误或至少描述结构不够明确。一旦不能生成SMCRA数据结构，就有理由怀疑出了故障。因此，解析PDB文件可用于检测可能的故障。我们将在 `11.7.1 <#problem%20structures>`__ 小节给出关于这一点的一些例子。

    --------------

    |image3|
    +------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | 图11.1：用来表示大分子结构的 ``Structure`` 类的SMCRA体系的UML图。带方块的实线表示集合，带箭头的实线表示引用，带三角形的实线表示继承，带三角形的虚线表示接口实现。    |
    +------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

    --------------

结构，模型，链，残基都是实体基类的子类。原子类仅仅（部分）实现了实体接口（因为原子类没有子类）。

对于每个实体子类，你可以用该子类的一个唯一标识符作为键来提取子类（比如，可以用原子名称作为键从残基对象中提取一个原子对象；用链的标识符作为键从域对象中提取链）。

无序的原子和残基用DisorderedAtom和DisorderedResidue类来描述，二者都是DisorderedEntityWrapper基本类的子类。他们隐藏了和无序的复杂性，表现得与原子和残基对象无二。

一般地，一个实体子类（如原子，残基，链，域）能通过标识符键从父类中提取。

    >>> child_entity = parent_entity[child_id]

你可以从一个父Entity对象获得一个子Entity的列表。需要注意的是这个列表以一定的方式排列（例如根据在Model对象中链对象的链标识符）。



    >>> child_list = parent_entity.get_list()

你也可以从子类获得父类：

    >>> parent_entity = child_entity.get_parent()

在SMCRA的所有层次结构，你可以提取 *完整的ID* 。完整的id是一个含所有从顶层对象（结构分子）到当前对象的id字符串的元组。Residue对象的完整id可以这么获取：

    >>> full_id = residue.get_full_id()
    >>> print full_id
    ("1abc", 0, "A", ("", 10, "A"))

与这些相对应的：
 - id为iabc的结构
 - id为0的Model
 - id为A的链
 - id为(" ", 10, "A")的残基


Residue的id表示这个残基不是hetero残基（也不是水分子），因为hetero值为空；它的序列标识符为10，插入码为“A”。


欲获取entity的id，用 ``get_id`` 方法：

    >>> entity.get_id()

可以用 ``has_id`` 方法检查这个entity是否有给定id的子类：

    >>> entity.has_id(entity_id)

entity的长度等于子类的个数：


    >>> nr_children = len(entity)

它能删除，重命名，添加，比如从父entity得到的子entity并没有进行完整性检查（有可能添加两个相同id的残基到同一条链上）。这需要Decorator类来完成包括完整性检查，但是如果你想使用最初的接口的话可以查看源代码（Entity.py）。


11.2.1  结构
~~~~~~~~~~~~~~~~~

结构对象是成此结构中最高的一层。它的id是用户指定的字符串。结构包含一系列子结构域。大部分晶体结构（但不是全部）含有一个结构域，但是NMR结构经常含有多个结构域。晶体结构中大部分子的无序状态也能导致多个结构域。


11.2.2  结构域
~~~~~~~~~~~~~~~

结构域对象的id是一个整数，得自该结构域在被解析的文件中的位置（自动的从0开始标记）。晶体结构通常只有一个结构域（id为0），而NMR文件则含有多个结构域。然而许多PDB解析器都假定只有一个结构域， ``Bio.PDB`` 中的 ``Structure`` 就是为解决这个而设计的，它能轻易的处理含有多个结构域的PDB文件。


举例如下，从结构对象中获取第一个结构域：


    >>> first_model = structure[0]

结构域对象储存着一个链对象的列表。


11.2.3  链
~~~~~~~~~~~~~~~~

链对象的id取自PDB/mmCIF文件中链的标识符，是个单字符（通常是一个字母）。结构域中的链都有一个唯一的id。例如，从一个结构域对象取出标识符为“A”的链：


    >>> chain_A = model["A"]

链对象储存着一个残基对象的列表。


11.2.4  残基
~~~~~~~~~~~~~~~~~~~

一个残基标识符是三个元素组成的元组：
 - **hetero-field** (hetfield)，即：
	- ``'W'`` 代表水分子
	-  氨基酸残基名字后面紧跟的 ``'H_'`` 代表其他hetero残基（例如 ``'H_GLC'`` 表示一个葡萄糖分子） 
	- 空值表示标准的氨基酸和核酸

   图表在 `11.4.1 <#hetero%20problems>`__ 部分描述。
 -   **序列标识符** （resseq），一个描述改残基在链上位置的整数（如100）；
 -  **插入码** （icode），一个字符串，如“A”。插入码有时候用来保持合适residue numbering scheme。一个Ser80的点突变（在Thr80和Asn81间插入）会有如下的序列标识符和插入码：Thr 80 A, Ser 80 B, Asn 81。通过这种方式，residue numbering scheme保持与野生型结构一致。 

上述葡萄糖的id可以是 ``(’H_GLC’, 100, ’A’)`` 。如果hetero标签和插入码为空，序列标识符为：

    # Full id
    >>> residue=chain[(' ', 100, ' ')]
    # Shortcut id
    >>> residue=chain[100]

hetero标签的起因是许多许多PDB文件用相同的序列标识符表示氨基酸和hetero残基或水分子，这会产生一个很明显的问题如果不是用hetero标签的话。

不令人吃惊的是，一个Residue对象存储这一系列的子Atom，也包含一个表示残基名字的字符串（如 “ASN”）和残基的片段标识符（这对X-PLOR的用户来说很熟悉，但是在SMCRA数据结构的构建红却不使用）。


让我们来看一些例子。插入码为空的Asn 10的残基id可为 ``(’ ’, 10, ’ ’)`` ；Water 10的残基id为 ``(’W’, 10, ’ ’)``；一个序列标识符是10的葡萄糖分子（一个残基名字是GLC的hetero残基）的残基id为 ``(’H_GLC’, 10, ’ ’)`` 。在这种情况下，三个相同插入码和序列标识符的残基可以位于同一条链上，因为他们的残基id是不同的。


在大多数情况下，hetflag和插入码都会为空，如 ``(’ ’, 10, ’ ’)`` 。在那些情况下，序列标识符可以用来作为完整的id：


    # use full id
    >>> res10 = chain[(' ', 10, ' ')]
    # use shortcut
    >>> res10 = chain[10]

在Chain上的每个Residue对象都应该有一个唯一的id。然而，无序残基以一种特别的方式对待，详见 `11.3.3 <#point%20mutations>`__ 。


一个Residue对象还有如下额外的方法：



    >>> residue.get_resname()       # returns the residue name, e.g. "ASN"
    >>> residue.is_disordered()     # returns 1 if the residue has disordered atoms
    >>> residue.get_segid()         # returns the SEGID, e.g. "CHN1"
    >>> residue.has_id(name)        # test if a residue has a certain atom

你可以用 ``is_aa(residue)`` 来测试一个Residue对象是否是氨基酸。


11.2.5  原子Atom
~~~~~~~~~~~~

Atom对象储存这所有和原子有关的数据，没有子类。原子的id是它的原子名字（如，“OG”代表Ser残基侧链的氧原子）。在Residue中Atom的id需要是唯一的。而且，对于无序原子会产生异常，见 `11.3.2 <#disordered%20atoms>`__ 描述。


原子id是简单的原子名字（如 ``’CA’`` ）。在实际中，原子名字通过去除PDB文件中原子名字中的空格创建的。


可是，在PDB文件中空格可以是原子名字的一部分。通常，钙原子称为 ``’CA..’`` 是为了和Cα原子（叫做 ``’.CA.’`` ）区分开。在这种情况，如果空格去除则会产生问题（如统一个残基中的两个原子都叫做 ``’CA’`` ），所以空格保留。


在PDB文件中，一个原子名字有4个字符组成，通常头尾皆为空格。为了方便使用空格经常被移除（在PDB文件中氨基酸的Cα原子标记为“.CA.”，点表示空格）。为了生成原子名字（然后是原子id），空格会被移除，除非回造成名字冲突（如两个Atom对象有相同的名字和id）。对于后者，会尝试原子名字包含空格。这种情况可能会发生当残基含名字为“.CA.”和“CA..”的原子，尽管这不怎么可能。

储存的原子数据包括原子名字，原子坐标（如果有的话还包括标准差），B因子（包括anisotropicB因子和可能存在的标准差），altloc标识符和完整的包括空格的原子名字。少用有时没有在PDB文件中存储的原子序号和原子电荷。

为了操纵原子坐标，可以用 ``’CA’`` 对象的 ``transform`` 方法。用 ``set_coord`` 方法可以直接指定原子坐标。

一个Atom对象还有如下额外的方法：


    >>> a.get_name()       # atom name (spaces stripped, e.g. "CA")
    >>> a.get_id()         # id (equals atom name)
    >>> a.get_coord()      # atomic coordinates
    >>> a.get_vector()     # atomic coordinates as Vector object
    >>> a.get_bfactor()    # isotropic B factor
    >>> a.get_occupancy()  # occupancy
    >>> a.get_altloc()     # alternative location specifier
    >>> a.get_sigatm()     # standard deviation of atomic parameters
    >>> a.get_siguij()     # standard deviation of anisotropic B factor
    >>> a.get_anisou()     # anisotropic B factor
    >>> a.get_fullname()   # atom name (with spaces, e.g. ".CA.")

siguij，anisotropic B因子和sigatm Numpy可以用来表示原子坐标。

方法 ``get_vector`` 回返回一个代表 ``Atom``  对象坐标的 ``Vector`` 对象，可以对原子坐标做向量操作。 ``Vector`` 实现了完整的三维向量操作、矩阵乘法（包括左和右）和一些高级的旋转相关的操作。


作为Bio.PDB的 ``Vector`` 模块的性能展示的一个例子，假设你查找Gly氨基酸残基的Cβ原子的位置,如果存在的话。将Gly残基的N原子沿Cα-C旋转-120度，能大致将其放在一个虚拟的Cβ原子的位置上。下面是使用 ``Vector`` 模块中的``rotaxis`` 方法（能用来完成绕一个轴的旋转）如何实现的：


    # get atom coordinates as vectors
    >>> n = residue['N'].get_vector() 
    >>> c = residue['C'].get_vector() 
    >>> ca = residue['CA'].get_vector()
    # center at origin
    >>> n = n - ca 
    >>> c = c - ca 
    # find rotation matrix that rotates n 
    # -120 degrees along the ca-c vector
    >>> rot = rotaxis(-pi * 120.0/180.0, c)
    # apply rotation to ca-n vector
    >>> cb_at_origin = n.left_multiply(rot)
    # put on top of ca atom
    >>> cb = cb_at_origin+ca

这个例子展示了在原子数据上能进行非常不一样的向量操作的可能性，这些操作是非常有用的。另外，所有有用的向量操作（交叉（用 ``*``\ ``*`` ），点积（用 ``*`` ），angle, norm等）和上述提到的 ``rotaxis`` 函数，``Vector`` 模块也有方法实现将一个向量叠加（ ``rotmat`` ）或反射（ ``refmat`` ）到两外一个向量上。


11.2.6  从结构分子中提取特定的 ``Atom/Residue/Chain/Model`` 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

例子如下：



    >>> model = structure[0]
    >>> chain = model['A']
    >>> residue = chain[100]
    >>> atom = residue['CA']

还有比较简单的方式：

    >>> atom = structure[0]['A'][100]['CA']

11.3  Disorder
--------------

Bio.PDB能够处理无序原子和点突变（比如Gly和Ala残基在相同位置上）。


11.3.1  一般途径 
~~~~~~~~~~~~~~~~~~~~~~~~

杂乱性可以从两个角度来解决：原子和残基的角度。一般来说，我们尝试压缩所有的会增加杂乱性。如果你仅仅想遍历所有Cα原子，你不必关注一些含杂乱侧链的残基。另一方面，应该考虑在数据结构中描述杂乱性。因此，杂乱原子或残基存储在特定的对象中表现得毫无杂乱性。这可以通过描述杂乱原子或残基的子集来完成。至于哪个子集（例如用到Ser残基侧链上的两个杂乱的OG原子）被挑选出来由用户来决定。


11.3.2  杂乱的原子
~~~~~~~~~~~~~~~~~~~~~~~~

杂乱原子可以用普通的 ``Atom`` 对象来描述，但是所有描述相同物理物理原子的 ``Atom`` 对象储存在一个 ``DisorderedAtom`` 对象中（见图Fig. `11.1 <#fig:smcra>`__ ）。 ``DisorderedAtom`` 对象中每个 ``Atom`` 对象都可以被它的altloc标识符唯一的索引。 ``DisorderedAtom`` 对象转寄所有uncaught method calls到选定的Atom对象，用过它的altloc标识符。以这种方式，原子杂乱性哪呢个正确的呗描述而没有额外的复杂性。换言之，如果你对杂乱原子不敢兴趣，你不会被它困扰。


每个杂乱原子有个特有的altloc标识符。你可以指定一个 ``DisorderedAtom`` 对象表现得像 ``Atom`` 对象和特定的altloc标识符关联：



    >>> atom.disordered_select('A') # select altloc A atom
    >>> print atom.get_altloc()
    "A"
    >>> atom.disordered_select('B') # select altloc B atom
    >>> print atom.get_altloc()
    "B"

11.3.3  杂乱的残基
~~~~~~~~~~~~~~~~~~~~~~~~~~~

普通例子
^^^^^^^^^^^

最常见的例子是一个残基包含一个或多个杂乱原子。这显然可以通过用DisorderedAtom对象描述这些杂乱原子来解决，并将DisorderedAtom对象保存在Residue对象中想正常的Atom对象一样。DisorderedAtom对象通过转寄所有uncaught method calls到其中一个Atom对象（被选中的Atom对象）一个表现的非常像正常的原子（事实上这个原子有最高的使用率）。


点突变
^^^^^^^^^^^^^^^

当点突变导致杂乱会产生一个特殊的例子，例如，当多肽中有两个或更多的点突变展示在晶体结构中。这个例子可以在PDB结构1EN2中找到。


既然这些残基属于不同的残基类型（举例说Ser 60 和Cys 60），他们不应该存储在一个 ``Residue`` 对象中像普通情况一样。在这是，每个残基被描述成 ``Residue`` 对象，所有 ``Residue`` 对象保存在一个 ``DisorderedResidue`` 对象中（见 Fig. `11.1 <#fig:smcra>`__ ）。


``DisorderedResidue`` 对象转寄所有uncaught methods到选定的 ``Residue`` 对象（默认下最后一个 ``Residue`` 对象被添加），然后表现得的一个正常的残基一样。在 ``DisorderedResidue`` 中的每个 ``Residue`` 对象通过残基名字被势必人出来。在上述例子中，残基Ser60在 ``DisorderedResidue`` 对象中的id为“SER”，而残基Cys 60则是“CYS”。他们能在 ``DisorderedResidue`` 中通过这个id选择现行的 ``Residue`` 对象。


例子：假设一个子类在10位有一个由Ser和Cys构成的点突变。让这个链的10位残基表现为Cys残基。



    >>> residue = chain[10]
    >>> residue.disordered_select('CYS')

另外，你能获得所有 ``Atom`` 对象的列表（如所有 ``DisorderedAtom`` 对象从他们各自的 ``Atom`` 对象中’unpacked’），通过使用 ``(Disordered)Residue`` 对象的 ``get_unpacked_list`` 方法。


11.4  Hetero残基
---------------------

11.4.1  相关问题
~~~~~~~~~~~~~~~~~~~~~~~~~~~

关于hetero残基的一个很普遍的问题是若干hetero和非hetero残基在同一条链中有同样的序列标识符和插入码。因此，为了生成每个hetero残基唯一的id，水分子和其他hetero残基应该以不同的方式来对待。


记住Residue残基有一个元组（hetfield, resseq, icode）作为id。hetfield值为空(“ ”)表示为氨基酸和核酸；为一个字符串表示水分子和其他hetero残基。hetfield的内容将在下面解释。

11.4.2  水分子
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

water residue的hetfield字符串由字母“W”组成。所以一个典型的水分子的残基id为(“W”, 1, “ ”)。

11.4.3  其他hetero残基
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

其他hetero残基的hetfield字符以“H\_”起始，后接残基名字。一个葡萄糖分子的残基名称为“GLC”，则hetfield字符为“H\_GLC”；它的残基id可以是(“H\_GLC”, 1,
“ ”)。

11.5  操纵Structure对象
-------------------------------------------

解析PDB文件，提取若干Model、Chain、Residue和Atom对象 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



    >>> from Bio.PDB.PDBParser import PDBParser
    >>> parser = PDBParser()
    >>> structure = parser.get_structure("test", "1fat.pdb")
    >>> model = structure[0]
    >>> chain = model["A"]
    >>> residue = chain[1]
    >>> atom = residue["CA"]

遍历一个结构中的所有原子
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



    >>> p = PDBParser()
    >>> structure = p.get_structure('X', 'pdb1fat.ent')
    >>> for model in structure:
    ...     for chain in model:
    ...         for residue in chain:
    ...             for atom in residue:
    ...                 print atom
    ...

如果你想遍历一个结构中所有原子，这儿有个捷径可以走：


    >>> atoms = structure.get_atoms()
    >>> for atom in atoms:
    ...     print atom
    ...

类似地，遍历一条链中的原子，可以这么做：


    >>> atoms = chain.get_atoms()
    >>> for atom in atoms:
    ...     print atom
    ...

遍历结构域中的所有残基
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

或者，如果你想遍历在一条链上的所有残基：


    >>> residues = model.get_residues()
    >>> for residue in residues:
    ...     print residue
    ...

你也可以用 ``Selection.unfold_entities`` 函数来从一个结构中获取所有残基：


    >>> res_list = Selection.unfold_entities(structure, 'R')

或者从链上获得所有原子：


    >>> atom_list = Selection.unfold_entities(chain, 'A')

明显的是， ``A=atom, R=residue, C=chain, M=model, S=structure`` 。你可以用此返回上一层，如从一个 ``Atoms`` 列表回溯到一个唯一的 ``Residue`` 或 ``Chain`` 的列表：


    >>> residue_list = Selection.unfold_entities(atom_list, 'R')
    >>> chain_list = Selection.unfold_entities(atom_list, 'C')

更多信息详见API文档。

从链中提取hetero残基（如葡萄糖（GLC）resseq 10的那部分）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    >>> residue_id = ("H_GLC", 10, " ")
    >>> residue = chain[residue_id]

输出链上所有的hetero残基
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    >>> for residue in chain.get_list():
    ...    residue_id = residue.get_id()
    ...    hetfield = residue_id[0]
    ...    if hetfield[0]=="H":
    ...        print residue_id
    ...

输出一个结构分子中所有B因子大于50的CA原子的坐标
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    >>> for model in structure.get_list():
    ...     for chain in model.get_list():
    ...         for residue in chain.get_list():
    ...             if residue.has_id("CA"):
    ...                 ca = residue["CA"]
    ...                 if ca.get_bfactor() > 50.0:
    ...                     print ca.get_coord()
    ...

输出所有含无序原子的残基
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    >>> for model in structure.get_list():
    ...     for chain in model.get_list():
    ...         for residue in chain.get_list():
    ...             if residue.is_disordered():
    ...                 resseq = residue.get_id()[1]
    ...                 resname = residue.get_resname()
    ...                 model_id = model.get_id()
    ...                 chain_id = chain.get_id()
    ...                 print model_id, chain_id, resname, resseq
    ...

遍历所有无序原子，并选取所有altloc A的原子（如果有的话）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

这需要确定的是SMCRA数据结构表现得如同altloc A原子存在。

    >>> for model in structure.get_list():
    ...     for chain in model.get_list():
    ...         for residue in chain.get_list():
    ...             if residue.is_disordered():
    ...                 for atom in residue.get_list():
    ...                     if atom.is_disordered():
    ...                         if atom.disordered_has_id("A"):
    ...                             atom.disordered_select("A")
    ...

从 ``Structure`` 对象中提取多肽
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

为了从一个结构中提取多肽，需要用 ``PolypeptideBuilder`` 从 ``Structure`` 构建一个 ``Polypeptide`` 对象的列表，如下所示：

    >>> model_nr = 1
    >>> polypeptide_list = build_peptides(structure, model_nr)
    >>> for polypeptide in polypeptide_list:
    ...     print polypeptide
    ...

Polypeptide对象是Residue对象的一个简单UserList，总是从单结构域中创建（在此例中为结构域1）。你可以从刚生成的 ``Polypeptide`` 对象提取序列作为 ``Seq`` 对象，或获得Cα原子的列表。多肽可以通过C-N 或a Cα-Cα距离标准来建立。


例子：


    # Using C-N 
    >>> ppb=PPBuilder()
    >>> for pp in ppb.build_peptides(structure): 
    ...     print pp.get_sequence()
    ...
    # Using CA-CA
    >>> ppb=CaPPBuilder()
    >>> for pp in ppb.build_peptides(structure): 
    ...     print pp.get_sequence()
    ...

需要注意的是，在上述例子中这个结构中，只有结构域0被 ``PolypeptideBuilder`` 考虑。尽管如此，还是可以用 ``PolypeptideBuilder`` 从 ``Model`` 和 ``Chain`` 对象创建 ``Polypeptide`` 对象。


从结构中获取序列
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

要做的第一件事就是从结构中提取所有多肽（如上述）。每条多肽的序列可以轻易的通过 ``Polypeptide`` 对象获得。这个序列当作Biopython中 ``Seq`` 对象来描述，它的字母表由 ``ProteinAlphabet`` 对象来定义。

例子：



    >>> seq = polypeptide.get_sequence()
    >>> print seq
    Seq('SNVVE...', <class Bio.Alphabet.ProteinAlphabet>)

11.6  结构分析
--------------------------

11.6.1  测定距离
~~~~~~~~~~~~~~~~~~~~~~~~~~~

两个已叠加的原子通过减法运算可以返回两个原子之间的距离。


    # Get some atoms
    >>> ca1 = residue1['CA']
    >>> ca2 = residue2['CA']
    # Simply subtract the atoms to get their distance
    >>> distance = ca1-ca2

11.6.2  测定角度
~~~~~~~~~~~~~~~~~~~~~~~~

用向量表示原子坐标，然后用 ``Vector`` 模块中的 ``calc_angle`` 函数可以计算角度。


    >>> vector1 = atom1.get_vector()
    >>> vector2 = atom2.get_vector()
    >>> vector3 = atom3.get_vector()
    >>> angle = calc_angle(vector1, vector2, vector3)

11.6.3  测定扭转角Measuring torsion angles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

用向量表示原子坐标，然后用 ``Vector`` 模块中的 ``calc_dihedral`` 函数可以计算角度。



    >>> vector1 = atom1.get_vector()
    >>> vector2 = atom2.get_vector()
    >>> vector3 = atom3.get_vector()
    >>> vector4 = atom4.get_vector()
    >>> angle = calc_dihedral(vector1, vector2, vector3, vector4)

11.6.4  确定原子-原子接触Determining atom-atom contacts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

用 ``NeighborSearch`` 实现临近查询。临近查询能够使用用C语言写的KD树模块（见 ``Bio.KDTree`` ）来比较快速的实现。它也包含了一个比较快速的方法找出一定距离内所有成对的点。

11.6.5  叠加两个结构
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

可以用 ``Superimposer`` 对象将两组坐标数据重叠。这个对象计算叠加和转化矩阵，这个矩阵通过叠加两个原子的列表以他们的RMSD最小的方式叠加在一起得到的矩阵。当然这两个列表含有相同数目的原子。 ``Superimposer`` 对象将叠加/转化应用在原子上。叠加和转化作为元组储存在 ``Superimposer`` 对象的 ``rotran`` 属性（注意的是叠加是right multiplying），RMSD储存在属性 ``rmsd`` 中。


``Superimposer`` 使用的算法来自[`17 <#golub1989>`__,
Golub & Van Loan]，用有意义的分解值（这在 ``Bio.SVDSuperimposer`` 模块中实现）。

例子：



    >>> sup = Superimposer()
    # Specify the atom lists
    # 'fixed' and 'moving' are lists of Atom objects
    # The moving atoms will be put on the fixed atoms
    >>> sup.set_atoms(fixed, moving)
    # Print rotation/translation/rmsd
    >>> print sup.rotran
    >>> print sup.rms 
    # Apply rotation/translation to the moving atoms
    >>> sup.apply(moving)

为了将结构根据他们的现行的位点叠加在一起，用active位点的原子计算叠加和转化矩阵（类似上述），将那些应用到整个分子。


11.6.6  相互映射这两个结构的残基到对方
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

首先，创建一个FASTA格式的比对文件，然后用``StructureAlignment`` 类。这个类也可以用来映射多余两个结构和比对。

11.6.7  计算Half Sphere Exposure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Half Sphere Exposure (HSE)是一种新的，2D测量solvent exposure [`20 <#hamelryck2005>`__\ ]。基本上，它计算围绕一个残基的Cα原子的个数，在它侧链的方向上，及反方向（在13 Å范围内）。除了简单，它还比其他solvent exposure测量工具好。


HSE有两种HSEα和HSEβ。前者仅用Cα原子的位置，而后者都用Cα和Cβ原子位置。HSE测定的通过 ``HSExposure`` 类计算的，也能计算原子接触数目。后者能返回一个将``Residue`` 对象映射到相应的HSEα,HSEβ和接触数目值的字典。


例子：



    >>> model = structure[0]
    >>> hse = HSExposure()
    # Calculate HSEalpha
    >>> exp_ca = hse.calc_hs_exposure(model, option='CA3')
    # Calculate HSEbeta
    >>> exp_cb=hse.calc_hs_exposure(model, option='CB')
    # Calculate classical coordination number
    >>> exp_fs = hse.calc_fs_exposure(model)
    # Print HSEalpha for a residue
    >>> print exp_ca[some_residue]

11.6.8  确定二级结构
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

与这个功能，你需要安装DSSP（获得一个对学术免费的证书，参见 ```http://www.cmbi.kun.nl/gv/dssp/`` <http://www.cmbi.kun.nl/gv/dssp/>`__ ）。然后用 ``DSSP`` 类，可以 ``Residue`` 对象到他们的二级结构上（和accessible surface area）。DSSP代码如下表所列 Table `11.1 <#cap:DSSP-codes>`__ 。注意DSSP（程序及其相应的类）不能处理多个结构域！

    --------------

    +--------+-----------------------------+
    | Code   | Secondary structure         |
    +--------+-----------------------------+
    | H      | α-helix                     |
    +--------+-----------------------------+
    | B      | Isolated β-bridge residue   |
    +--------+-----------------------------+
    | E      | Strand                      |
    +--------+-----------------------------+
    | G      | 3-10 helix                  |
    +--------+-----------------------------+
    | I      | Π-helix                     |
    +--------+-----------------------------+
    | T      | Turn                        |
    +--------+-----------------------------+
    | S      | Bend                        |
    +--------+-----------------------------+
    | -      | Other                       |
    +--------+-----------------------------+

    +--------------------------------------+
    | Table 11.1: Bio.PDB中的DSSP代码。   |
    +--------------------------------------+

    --------------

``DSSP`` 类也可以用来计算残基的易接近的表面。但是也参考 `11.6.9 <#subsec:residue_depth>`__ 。


11.6.9  计算残基深度
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

残基深度残基的原子到solvent accessible surface的平均距离。它是个费城新颖和强大的描述solvent accessibility的参数。为了这个功能，你需要安装Michel Sanner的 MSMS程序（ ```http://www.scripps.edu/pub/olson-web/people/sanner/html/msms_home.html`` <http://www.scripps.edu/pub/olson-web/people/sanner/html/msms_home.html>`__ ）。然后用 ``ResidueDepth`` 类。这个类想字典一样将 ``Residue`` 对象映射到相应的元组（残基深度，Cα深度）。Cα深度是残基的Cα原子到solvent accessible surface的距离。


例子：



    >>> model = structure[0]
    >>> rd = ResidueDepth(model, pdb_file)
    >>> residue_depth, ca_depth=rd[some_residue]

你也可以获得分子表面本身（通过 ``get_surface`` 函数），以Python数组的形式的形式和表面点。You can also get access to the molecular surface itself (via the
``get_surface`` function), in the form of a Numeric Python array with
the surface points.

11.7  PDB文件中的常见问题
----------------------------------

都知道一些PDB文件有语义错误（不是结构本身的错误，而是在PDB文件中描述时出错）。Bio.PDB可以有两种途径来处理这个问题。PDBParser对象能以两种方式处理：严格（restrictive）方式和宽松（permissive）方式（默认方式）：


例子:



    # Permissive parser
    >>> parser = PDBParser(PERMISSIVE=1)
    >>> parser = PDBParser() # The same (default)
    # Strict parser
    >>> strict_parser = PDBParser(PERMISSIVE=0)

在宽松状态（默认），还有错误的PDB文件会被认为是“正确的”（比如说一些残基或原子丢失）。这些错误包括：
 - 多个残基使用同一个标识符
 - 多个原子使用统一个标识符（考虑altloc识别符）


这些错误暗示了PDB文件中确实在错误（详情见 [`18 <#hamelryck2003a>`__, Hamelryck and Manderick, 2003] ）。在严格模式，PDB文件中的这些错误会报错，这能帮助发现PDB文件的存在的错误。

有些错误能自动修正。正常情况下，每个无序原子应该会有一个非空altloc标识符，可是有些结构没有遵循这个惯例，在用同一个原子上会同时存在两个无序位置的空的和非空的标识符。这个错误能够以正确的方式被解析。


有时候一个结构会有这样的情况：一部分残基属于A链，接下来一部分残基属于B链，然后又有一部分残基属于A链，这种链称为“断链”，这也能被正确的解析。


11.7.1  例子
~~~~~~~~~~~~~~~~

PDBParser/Structure类经将近800个结构测试（每个都属于不同的SCOP超家族），总共花费20分钟左右，也就是说平均每个结构只需1.5秒。在一台1000 MHz的PC上解析巨大的包含近64000个原子的核糖体亚单位（1FKK）只需10秒。

当明确的数据结构不能建立的时候会有三个异常发生。在这三个异常中，可能的起因是PDB文件中本应修正的错误没有被修正。产生异常总要比冒险不正确描述一个在那样的数据结构中的结构分子好的多。


11.7.1.1  重复残基
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

一个结构分子在相同的序列标识符（resseq3）和icode的链上存在两个相同的氨基酸残基。仔细观察可以发现这条链包含残基：Thr A3, …, Gly A202, Leu A3；很明显第二个Leu A3应该是Leu A203。类似的情况也存在于1FFK结构上（含残基Gly B64, Met B65, Glu B65, Thr B67；第二个Glu B65应该是Glu B66）。


11.7.1.2  重复原子
^^^^^^^^^^^^^^^^^^^^^^^^^

结构1EJG含有在A链22位的Ser/Pro点突变。依次，Ser 22含一些杂乱原子。和期望的一样，所有属于 Ser 22的原子都有一个非空的altloc标识符（B或C）。所有Pro 22的原子都为altloc A，除了含空altloc的N原子。这会申城一个异常，因为所有属于一个点突变的两个残基的原子都应该有非空的altloc。这导致这个原子被Ser 和 Pro 22共用，Ser22丢失了这个N原子。再者，这意味着在文件中的问题：这个N原子应该在Ser和Pro残基中都有描述，都与合适的altloc标识符关联。


11.7.2  自动修正
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

一些错误比较普遍，能够在没有太大解释错误的风险下被修改。这些错误如下所述。

11.7.2.1  无序原子的空altloc 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

正常情况下，每个无序原子应该会有一个非空altloc标识符，可是有些结构没有遵循这个惯例，在用同一个原子上会同时存在两个无序位置的空的和非空的标识符。这个错误能够以正确的方式被解析。

11.7.2.2  断链
^^^^^^^^^^^^^^^^^^^^^^^

有时候一个结构会有这样的情况：一部分残基属于A链，接下来一部分残基属于B链，然后又有一部分残基属于A链，这种链称为“断链”，这也能被正确的解析。

11.7.3  致命错误Fatal errors
~~~~~~~~~~~~~~~~~~~~

有时候一个PDB文件不能被明确的解释，这会产生异常而不是猜测和冒出错的风险，等待用户去修正这个PDB文件。这些异常如下所述。

11.7.3.1  残基重复
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

所有在一条链上的残基应该有一个唯一的id，这id基于下述生成：

 - 序列标识符（resseq）
 - 插入码（icode）
 - hetfield字符（“W”代表水分子，残基名字后面的“H\_”代表其他异性残基）
 - 发生点突变的残基名字（在DisorderedResidue对象中保存Residue对象）


如果这样还不能生成一个唯一的id，则会生成异常。


11.7.3.2  原子重复
^^^^^^^^^^^^^^^^^^^^^^^^^

一个残基上所有原子应该有一个唯一的id，这个id基于下述产生：
 - 原子名称（没有空格，否则会报错）
 - altloc标识符

如果这不能生成一个唯一的id，则会生成异常。

11.8  访问Protein Data Bank
-------------------------------------

11.8.1  从Protein Data Bank下载结构
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

结构可以从PDB数据库（Protein Data Bank）通过 ``PDBList`` 对象的 ``retrieve_pdb_file`` 方法下载。这种方法的要点是结构在PDB中的唯一标识符。Structures can be downloaded from the PDB (Protein Data Bank) by using
the ``retrieve_pdb_file`` method on a ``PDBList`` object. The argument
for this method is the PDB identifier of the structure.

    >>> pdbl = PDBList()
    >>> pdbl.retrieve_pdb_file('1FAT')

``PDBList`` 类也能作为命令行工具来使用：

    python PDBList.py 1fat

下载的文件将以 ``pdb1fat.ent`` 为名保存在当前工作目录。 ``retrieve_pdb_file`` 方法有个可选参数 ``pdir`` 来指定路径保存所下载的PDB文件。

``retrieve_pdb_file`` 方法还有其他选项可以指定下载的压缩格式（默认的 ``.Z`` 格式和 ``gunzip`` 格式）。另外，在创建 ``PDBList`` 对象时还可以指定PDB ftp站点。一般使用Worldwide Protein Data Bank（ ```ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/`` <ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/>`__ ）。详细内容参见API文档。再次感谢Kristian Rother对此模块的所做的贡献。

11.8.2  下载完整PDB数据 Downloading the entire PDB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

下面的命令将会保存所有PDB文件至 ``/data/pdb`` 路径：


    python PDBList.py all /data/pdb

    python PDBList.py all /data/pdb -d

在API中这个方法叫做 ``download_entire_pdb`` 。添加 ``-d`` 能保存所有文件在相同的路径。否则将分别保存至和它们PDB ID相对应的子目录中。根据网速，完整的下载全部PDB文件大概需要2-4天。


11.8.3  保持本地拷贝与PDN数据库更新
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

这也能通过 ``PDBList`` 对象来完成。可以简单的创建一个 ``PDBList`` 对象（指定本地拷贝保存的路径），然后调用 ``update_pdb`` 方法：
This can also be done using the ``PDBList`` object. One simply creates a
``PDBList`` object (specifying the directory where the local copy of the
PDB is present) and calls the ``update_pdb`` method:


    >>> pl = PDBList(pdb='/data/pdb')
    >>> pl.update_pdb()

当然还可以用 ``cronjob`` 实现本地拷贝每周的自动更新。可以指定PDB ftp站点（详见API文档）。

``PDBList`` 有其他许多另外的方法可供调用。 ``get_all_obsolete`` 方法可以获取已经废弃不用的PDB entries；  ``changed_this_week``  方法可以获得当前一周内新增加、修改或废弃的PDB entries。更多 ``PDBList`` 的用法参见API文档。


11.9  一般疑问
-------------------------------

11.9.1  对Bio.PDB测试如何？ 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

事实上，相当好。Bio.PDB已经在从PDB获得的近5500个结构上广泛的测试过，所有文件都能正常的被解析。更多细节可以参考在Bioinformatics上发表的关于Bio.PDB的文章。Bio.PDB已经并且正在被作为很使用的工具用于许多研究项目中。我几乎每天都在用它，处于研究目的、提升其性能和增加新属性。


11.9.2  它有多快？
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``PDBParser`` 的性能经将近800个结构测试（每个都属于不同的SCOP超家族），总共花费20分钟左右，也就是说平均每个结构只需1.5秒。在一台1000 MHz的PC上解析巨大的包含近64000个原子的核糖体亚单位（1FKK）只需10秒。总而言之，它比其他的一些应用程序更快。


11.9.3  是否支持分子图形展示？
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

现在已经有很多应用直接或间接用Pyhton解决问题的应用，将来也可能回用到Bio.PDB。我的选择是Pymol，我在Pymol中使用Bio.PDB非常成功，将来会有特定的Bio.PDB中会有特定的PyMol模块。基于Python的分子图形展示的解决方案包括：

-  PyMol:
   ```http://pymol.sourceforge.net/`` <http://pymol.sourceforge.net/>`__
-  Chimera:
   ```http://www.cgl.ucsf.edu/chimera/`` <http://www.cgl.ucsf.edu/chimera/>`__
-  PMV:
   ```http://www.scripps.edu/~sanner/python/`` <http://www.scripps.edu/~sanner/python/>`__
-  Coot:
   ```http://www.ysbl.york.ac.uk/~emsley/coot/`` <http://www.ysbl.york.ac.uk/~emsley/coot/>`__
-  CCP4mg:
   ```http://www.ysbl.york.ac.uk/~lizp/molgraphics.html`` <http://www.ysbl.york.ac.uk/~lizp/molgraphics.html>`__
-  mmLib: ```http://pymmlib.sourceforge.net/`` <http://pymmlib.sourceforge.net/>`__ 
-  VMD:
   ```http://www.ks.uiuc.edu/Research/vmd/`` <http://www.ks.uiuc.edu/Research/vmd/>`__
-  MMTK:
   ```http://starship.python.net/crew/hinsen/MMTK/`` <http://starship.python.net/crew/hinsen/MMTK/>`__

11.9.4  谁在用Bio.PDB？ 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bio.PDB曾用于构建DISEMBL，一个能预测蛋白结构中的非规则区域( ```http://dis.embl.de/`` <http://dis.embl.de/>`__ )；COLUMBA，一个提供经注释的蛋白结构的站点( ```http://www.columba-db.de/`` <http://www.columba-db.de/>`__ )。Bio.PDB也用于PDB数据库中蛋白质间大规模活性位点相似性的搜索[`19 <#hamelryck2003b>`__, Hamelryck, 2003]，开发新的算法鉴别线性二级结构元件[`26 <#majumdar2005>`__, Majumdar *et al.*, 2005]。

基于新属性和信息的需求反馈，也可以得知Bio.PDB也在许多大型制药公司中使用。



