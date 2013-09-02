第 22 章 附录：Python之外
===============================================

如果你对于Python编程还不是非常熟练，那么你在使用Biopython过程中遇到的问题都跟Python有关。本章节主要向读者提供一些使用Biopython文库的过程中有用的建议和一些常用代码。关于这一章的内容，你要是有什么好的建议的话，请一定告诉我们。

22.1  到底什么是句柄（handle）？
--------------------------------

这个文档中，句柄常被提及，而且也比较难理解（至少对我来说）。一般说来，你可以把句柄想象成一个对文本信息的“封装”。

对于普通文本信息的处理，使用句柄至少有两个好处：
Handles provide (at least) two benefits over plain text information:

#. 对于以不同方式存储的信息，句柄提供了一个标准的处理方法。这些文本信息可能来自文件、内存中的一个字符串、命令行指令的输出或者来自于远端网站信息，但是句柄提供了一种通用的方式处理这些不同格式和来源的文本信息。
#. 句柄可以依次读取文本信息，而不是一次读取所有信息。这点在处理数据比较大的文件时尤为有用，因为一次载入一个很大的文件可能会占去你所有的内存。

不论是从文件读取文本信息还是将文本信息写入文件，句柄都能胜任。在读取文件时，常用的函数有read()和readline(), 前者可以通过句柄读取所有文本信息，而后者则每次读取一行；对于文本信息的写入，则通常使用write()函数。

句柄最常见的使用就是从文件读取信息，这可以通过Python内置函数open来完成。 下面示例中，我们打开一个指向文件m_cold.fasta（可通过网址http://biopython.org/DIST/docs/tutorial/examples/m_cold.fasta获取）的句柄：

.. code:: verbatim

    >>> handle = open("m_cold.fasta", "r")
    >>> handle.readline()
    ">gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum ...\n"

Biopython中句柄常用来向解析器（parsers）传递信息。比如说，自从Biopython1.54以来，Bio.SeqIO和Bio.AlignIO模块中的主要函数都可以使用文件名来代替句柄使用：

.. code:: verbatim

    from Bio import SeqIO
    for record in SeqIO.parse("m_cold.fasta", "fasta"):
        print record.id, len(record)

在比较早的BioPython版本中，必须使用句柄。

.. code:: verbatim

    from Bio import SeqIO
    handle = open("m_cold.fasta", "r")
    for record in SeqIO.parse(handle, "fasta"):
        print record.id, len(record)
    handle.close()

这种操作方式仍有其用武之地，比如在解析一个gzip压缩的FASTA文件中：

.. code:: verbatim

    import gzip
    from Bio import SeqIO
    handle = gzip.open("m_cold.fasta.gz")
    for record in SeqIO.parse(handle, "fasta"):
        print record.id, len(record)
    handle.close()

在5.2章节中SeqIO_compressed部分有更多此类示例可供参考，其中包括,bzip2压缩文件的读取。

22.1.1  从字符串创建句柄Creating a handle from a string
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

一个比较有用的工具是将字符串中包含的文本信息传递给一个句柄。以下示例中通过Python标准文库cStringIO来完成：

.. code:: verbatim

    >>> my_info = 'A string\n with multiple lines.'
    >>> print my_info
    A string
     with multiple lines.
    >>> from StringIO import StringIO
    >>> my_info_handle = StringIO(my_info)
    >>> first_line = my_info_handle.readline()
    >>> print first_line
    A string
    <BLANKLINE>
    >>> second_line = my_info_handle.readline()
    >>> print second_line
     with multiple lines.
