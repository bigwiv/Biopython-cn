第19章 Biopython测试框架
===========================================

Biopython具有一个基于Python标准单元测试框架 `unittest<http://docs.python.org/library/unittest.html>` 
的回归测试框架（文件 ``run_tests.py``）。而为模块提供全面测试，
是确保Biopython代码在使用期内尽可能无bug的一个最重要的方面。
也经常是最被轻视的方面之一。本章旨在使运行Biopython测试和编
写测试代码尽可能容易。理想情况下，进入Biopython的每个模块都
应该有一个测试（还应该有对应文档！）。强烈鼓励我们所有开发
者，以及任何从源码安装Biopython的人运行单元测试。

19.1  运行测试
-----------------------

在你下载Biopython源码或者从我们的源码仓库签出时，你会发现一
个子目录调用 ``Tests``。 这包括关键脚本 ``run_tests.py``、
名为 ``test_XXX.py`` 的很多独立脚本、一个叫 ``output`` 的子目录和
 很多其他包含测试套件输入文件的子目录。

作为构建和安装Biopython的一部分，你通常会在命令行上从Biopython
源码顶层目录运行整个测试套件如下：

.. code:: verbatim

    python setup.py test

这事实上等价于转到 ``Tests`` 子目录，并运行：

.. code:: verbatim

    python run_tests.py

你通常会想要只运行测试的一部分，这可以如下来操作：

.. code:: verbatim

    python run_tests.py test_SeqIO.py test_AlignIO.py

当给出测试列表时， ``.py`` 扩展名是可选的，所以你可以只需打字：

.. code:: verbatim

    python run_tests.py test_SeqIO test_AlignIO

要运行 docstring 测试（见 `19.3 <#section:doctest>`__ 节）的话，
你可以用

.. code:: verbatim

    python run_tests.py doctest

缺省情况下， ``run_tests.py`` 运行所有测试，包括docstring测试。

如果一单个测试失败了，你还可以尝试直接运行它，它会给出更多信息。

重要的是，要注意单个单元测试有两类作用：

-  简单打印和比较脚本。 这些单元测试本质上是简短的 Python 示例
   程序，它们会打印出各种输出文本。对于一个名为 ``test_XXX.py`` 
   的测试文件，在 ``output`` 子目录（包含期望的输出）下会有一个
   叫做 ``test_XXX`` 的匹配文本文件。测试框架所做的全部就是运行
   脚本并检查输出的一致性。
-  基于 ``unittest`` 的标准测试。 这些会 ``import unittest`` ，然
   后定义 ``unittest.TestCase`` 类，这些类的每一个都带有一或多个
   像以 ``test_`` 开头、检查代码的某些方面的方法那样的子测试。这
   些测试不应该直接打印任何输出。

目前，大约一半的 Biopython 测试是 ``unittest``风格的测试，另一半
是 print-and-compare 测试。

直接运行一个简单的 print-and-compare 测试通常会在屏幕上给出大量
输出，但是并不检查输出是否跟期望的一样。如果测试以一个意外的错误
而失败，那么应该很容易精确定位脚本失败的位置。例如，对于一个print-and-compare 测试，试一下：

.. code:: verbatim

    python test_SeqIO.py

基于 ``unittest`` 的测试反倒是精确地显示你测试的哪一个小块失败了。
例如，

.. code:: verbatim

    python test_Cluster.py

19.2  编写测试
-------------------

假如说你想为一个叫做 ``Biospam`` 的模块写一些测试。这可以是你写的
一个模块，或者是一个还没有任何测试的现存模块。在下面的例子中，我们
假设 ``Biospam`` 是一个做简单数学的模块。

每个 Biopython 测试都可以有三个重要的文件和相关目录：

#. ``test_Biospam.py`` – 关于你的模块的真正测试代码。
#. ``Biospam`` [optional]– 一个包含任何必要输入文件的目录。任何会
   生成的输出文件也应该写在这里（并且最好在测试结束后打扫干净)以防
   堵塞主 Tests 目录。
#. ``output/Biospam`` – [只针对 print-and-compare 测试] 这个文件
   包括运行 ``test_Biospam.py`` 的期望输出。这个文件对于 ``unittest`` 
   风格的测试不是必须的，因为测试脚本 ``test_Biospam.py`` 会自己做验证。
你要自己决定你是想编写一个 print-and-compare 测试脚本还是一个 ``unittest`` 
风格的测试脚本。重要的是你不能把这两种风格混合在一个
测试脚本中。尤其是，不要在一个 print-and-compare 测试中使用``unittest`` 
特征。

 ``Tests`` 目录中任何具有 ``test_`` 前缀的脚本都会被 ``run_tests.py`` 
找到并运行。下面，我们展示一个示例测试脚本 ``test_Biospam.py`` ，针对
一个 print-and-compare 测试和一个基于 ``unittest`` 的测试。如果你把这个
脚本放进 Biopython的 ``Tests`` 目录，那么 ``run_tests.py`` 就会找到它并
执行其中包含的测试：

.. code:: verbatim

    $ python run_tests.py     
    test_Ace ... ok
    test_AlignIO ... ok
    test_BioSQL ... ok
    test_BioSQL_SeqIO ... ok
    test_Biospam ... ok
    test_CAPS ... ok
    test_Clustalw ... ok

…

.. code:: verbatim

    ----------------------------------------------------------------------
    Ran 107 tests in 86.127 seconds

19.2.1  编写一个 print-and-compare 测试
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
一个 print-and-compare 风格的测试对于初学者和新手来说是很容易写的- 本质上它只是一个使用你的模块的示例脚本。

为了写一个关于 ``Biospam`` 的 print-and-compare 测试，这是你应该
做的：

#. 编写一个叫 ``test_Biospam.py`` 的脚本

   -  这个脚本应该位于 Tests 目录
   -  脚本应该测试模块的所有重要功能（当然，你测试的越多、你的测试就
      越好！）。
   -  尽量避免任何平台特异的东西，例如打印浮点数而不用显式格式字符串
      来避免有太多小数位（不同的平台会给出稍微不同的值）。

#. 如果脚本需要文件来进行测试，这些应转到目录 Tests/Biospam 中进行
   （如果你只需一些通用的东西，像一个 FASTA 序列文件，或者一条
    GenBank 记录，试着用一个现存的样品输入文件来代替）。
#. 写出测试输出并验证输出是正确的。

   有两种方法可以做到这一点：

   #. 长期方法：

      -  运行脚本并将输出写到一个文件中。在 UNIX （包括 Linux 和 Mac OS X 
         ）机器上，你可以这样做： ``python test_Biospam.py > test_Biospam`` 
         这会把输出写到文件 ``test_Biospam`` 中。
      -  手动查看文件 ``test_Biospam`` 来确保输出正确。当你确定都没问
         题、没有bug后，你需要快速编辑 ``test_Biospam`` 文件使其首行为：
          ‘\ ``test_Biospam``\ ’  （不含引号）。
      -  复制文件 ``test_Biospam`` 到目录 Tests/output 中。

   #. 快速方法:

      -  运行 ``python run_tests.py -g test_Biospam.py`` 。回归测试框架
         很聪明的会以他喜欢的方式把输出放在恰当的地方。
      -  转到输出（应该在 ``Tests/output/test_Biospam``）并复查输出以确
         保其完全正确。

#. 现在改换到 Tests 目录并运行 ``python run_tests.py`` 进行回归测试。
   这会运行所有测试，而你会看到你的测试也在运行（并通过）。
#. 好了！这样你就得到了可用于签入或提交到Biopython的、关于你的模块的
   一个友好的测试。恭喜你！

例如，测试 ``Biospam`` 模块中的 ``addition`` 和 ``multiplication`` 功
能的测试脚本 ``test_Biospam.py`` 也许看起来是下面这个样子：

.. code:: verbatim

    from Bio import Biospam

    print "2 + 3 =", Biospam.addition(2, 3)
    print "9 - 1 =", Biospam.addition(9, -1)
    print "2 * 3 =", Biospam.multiplication(2, 3)
    print "9 * (- 1) =", Biospam.multiplication(9, -1)

我们用 ``python run_tests.py -g test_Biospam.py`` 来生成对应的输出，
并检查输出文件 ``output/test_Biospam`` ：

.. code:: verbatim

    test_Biospam
    2 + 3 = 5
    9 - 1 = 8
    2 * 3 = 6
    9 * (- 1) = -9

通常，更大的 print-and-compare 测试的困难在于追踪输出行与测试脚本
命令之间的对应关系。为此，打印出一些标记是很重要的，这些标记帮助你
把输入脚本按行和产生的输出匹配起来。

19.2.2  编写一个基于 unittest 的测试
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

我们想要Biopython中的所有模块都具有单元测试，并且一个简单的 
print-and-compare 测试比一点儿测试都没有要好。不过，尽管有一个陡峭的
学习曲线，使用 ``unittest`` 框架能给出一个更结构化的结果，并且如果有
一个测试失败，这能够清晰准确地指出测试的哪部分出了问题。子测试也可以
单独运行，这对于测试和排错很有帮助。

从2.1版开始 ``unittest`` 框架就包含在Python中了，并且存档在 
Python Library Reference （就是所推荐的你的枕边书）。也有 `关于unittest
的在线文档 <http://docs.python.org/library/unittest.html>`__。如果你
熟悉 ``unittest`` 系统（或类似于某些噪音测试框架的东西），你应该不会有
什么麻烦。你也许发现，寻找Biopython中的现成例子很有帮助。

这是关于 ``Biospam`` 的一个 ``unittest`` 风格的极小测试脚本，你可以
复制粘贴过去运行它：

.. code:: verbatim

    import unittest
    from Bio import Biospam

    class BiospamTestAddition(unittest.TestCase):

        def test_addition1(self):
            result = Biospam.addition(2, 3)
            self.assertEqual(result, 5)

        def test_addition2(self):
            result = Biospam.addition(9, -1)
            self.assertEqual(result, 8)

    class BiospamTestDivision(unittest.TestCase):

        def test_division1(self):
            result = Biospam.division(3.0, 2.0)
            self.assertAlmostEqual(result, 1.5)

        def test_division2(self):
            result = Biospam.division(10.0, -2.0)
            self.assertAlmostEqual(result, -5.0)


    if __name__ == "__main__":
        runner = unittest.TextTestRunner(verbosity = 2)
        unittest.main(testRunner=runner)

在分割测试中，我们使用 ``assertAlmostEqual`` 而不是 ``assertEqual`` 
以免因舍入误差造成的测试失败；详情以及 ``unittest`` 中的其他可用功能
参见Python文档中的 ``unittest`` 章节（`在线参考 <http://docs.python.org/library/unittest.html>`__）。

这里是基于 ``unittest`` 的测试的一些关键点：

-  测试实例存储在 ``unittest.TestCase`` 的子类中并涵盖了你的代码
    的一个基本方面。
-  对于任何在每个测试方法前后都要运行的重复代码，你可以使用方法 
   ``setUp`` 和 ``tearDown`` 。例如 ``setUp`` 方法可用于创建你正在
   测试的对象的实例，或打开一个文件句柄。 ``tearDown`` 可做任何整理，
   例如关闭文件句柄。
-  测试以 ``test_`` 为前缀并且每项测试应覆盖你所想要测试的内容的一个
   具体部分。一个类中你想包含多少个测试都行。
-  在测试脚本的末尾，你可以用

   .. code:: verbatim

       if __name__ == "__main__":
           runner = unittest.TextTestRunner(verbosity = 2)
           unittest.main(testRunner=runner)

   来执行测试脚本，当脚本是从	自己运行（而不是从 ``run_tests.py`` 导入）时。
   如果你运行该脚本，那么你会见到类似下面的东西:

   .. code:: verbatim

       $ python test_BiospamMyModule.py
       test_addition1 (__main__.TestAddition) ... ok
       test_addition2 (__main__.TestAddition) ... ok
       test_division1 (__main__.TestDivision) ... ok
       test_division2 (__main__.TestDivision) ... ok

       ----------------------------------------------------------------------
       Ran 4 tests in 0.059s

       OK

-  为了更清晰地表明每个测试都干了什么，你可以给每个测试加上 docstrings 。
   它们会在运行测试的时候显示出来，如果一个测试失败这会是有用的信息。

   .. code:: verbatim

       import unittest
       from Bio import Biospam

       class BiospamTestAddition(unittest.TestCase):

           def test_addition1(self):
               """An addition test"""
               result = Biospam.addition(2, 3)
               self.assertEqual(result, 5)

           def test_addition2(self):
               """A second addition test"""
               result = Biospam.addition(9, -1)
               self.assertEqual(result, 8)

       class BiospamTestDivision(unittest.TestCase):

           def test_division1(self):
               """Now let's check division"""
               result = Biospam.division(3.0, 2.0)
               self.assertAlmostEqual(result, 1.5)

           def test_division2(self):
               """A second division test"""
               result = Biospam.division(10.0, -2.0)
               self.assertAlmostEqual(result, -5.0)


       if __name__ == "__main__":
           runner = unittest.TextTestRunner(verbosity = 2)
           unittest.main(testRunner=runner)

   运行脚本你就会看到：

   .. code:: verbatim

       $ python test_BiospamMyModule.py
       An addition test ... ok
       A second addition test ... ok
       Now let's check division ... ok
       A second division test ... ok

       ----------------------------------------------------------------------
       Ran 4 tests in 0.001s

       OK

如果你的模块包含 docstring 测试（见`19.3 <#section:doctest>`__小节），
你也许想在要运行的测试中包含这些。你可以修改 ``if __name__ == "__main__":`` 
下面的代码如下面这样：

.. code:: verbatim

    if __name__ == "__main__":
        unittest_suite = unittest.TestLoader().loadTestsFromName("test_Biospam")
        doctest_suite = doctest.DocTestSuite(Biospam)
        suite = unittest.TestSuite((unittest_suite, doctest_suite))
        runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
        runner.run(suite)

这只与你执行 ``python test_Biospam.py`` 时是否想要运行 docstring 测试
有关；运行 ``python run_tests.py`` ，docstring 测试会自动运行（假设他们
被包含在 ``run_tests.py`` 中的 docstring 测试列表中，见下面的小节）。

19.3  编写 doctests
----------------------

Python 模块、类和函数支持使用 docstrings 创建文档。 `doctest 框架
 <http://docs.python.org/library/doctest.html>`__ （包含在Python中）
 允许开发者将工作例子嵌入在 docstrings 中，并自动测试这些例子。

目前只有一小部分 Biopython 包含 doctests 。 ``run_tests.py`` 脚本
看护着 doctests 的运行。为此， ``run_tests.py`` 脚本开头是要测试
的模块的一个手动编译列表，该列表允许我们跳过那些可能没有安装可选
外部依赖库的模块（例如 Reportlab 和 NumPy 库）。所以如果你在 Biopython 
模块中加一些针对 dostrings 的 doctests ，为了把它们包含在 Biopython 
套件中，你必须更新 ``run_tests.py`` 以包含你的模块。现在，
 ``run_tests.py`` 的相关部分看起来像下面这样：

.. code:: verbatim

    # This is the list of modules containing docstring tests.
    # If you develop docstring tests for other modules, please add
    # those modules here.
    DOCTEST_MODULES = ["Bio.Seq",
                       "Bio.SeqRecord",
                       "Bio.SeqIO",
                       "...",
                      ]
    #Silently ignore any doctests for modules requiring numpy!
    try:
        import numpy
        DOCTEST_MODULES.extend(["Bio.Statistics.lowess"])
    except ImportError:
        pass

注意我们首先把 doctests 看做文档，所以你应该坚持典型用法。通常处理错
误条件等诸如此类的复杂例子最好留给一个专门的单元测试。

注意，如果你想编写涉及文件解析的 doctests ，定义文件位置复杂性是很要
紧的。理想情况下，假设代码会从 ``Tests`` 目录运行，使用相对路径即可，
关于这一点的一个例子参见 ``Bio.SeqIO`` doctests 。

要想只运行 docstring 测试，使用

.. code:: verbatim

    $ python run_tests.py doctest
