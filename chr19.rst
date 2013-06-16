Chapter?19??The Biopython testing framework
===========================================

Biopython has a regression testing framework (the file ``run_tests.py``)
based on `unittest <http://docs.python.org/library/unittest.html>`__,
the standard unit testing framework for Python. Providing comprehensive
tests for modules is one of the most important aspects of making sure
that the Biopython code is as bug-free as possible before going out. It
also tends to be one of the most undervalued aspects of contributing.
This chapter is designed to make running the Biopython tests and writing
good test code as easy as possible. Ideally, every module that goes into
Biopython should have a test (and should also have documentation!). All
our developers, and anyone installing Biopython from source, are
strongly encouraged to run the unit tests.

19.1??Running the tests
-----------------------

When you download the Biopython source code, or check it out from our
source code repository, you should find a subdirectory call ``Tests``.
This contains the key script ``run_tests.py``, lots of individual
scripts named ``test_XXX.py``, a subdirectory called ``output`` and lots
of other subdirectories which contain input files for the test suite.

As part of building and installing Biopython you will typically run the
full test suite at the command line from the Biopython source top level
directory using the following:

.. code:: verbatim

    python setup.py test

This is actually equivalent to going to the ``Tests`` subdirectory and
running:

.. code:: verbatim

    python run_tests.py

You¡¯ll often want to run just some of the tests, and this is done like
this:

.. code:: verbatim

    python run_tests.py test_SeqIO.py test_AlignIO.py

When giving the list of tests, the ``.py`` extension is optional, so you
can also just type:

.. code:: verbatim

    python run_tests.py test_SeqIO test_AlignIO

To run the docstring tests (see section `19.3 <#section:doctest>`__),
you can use

.. code:: verbatim

    python run_tests.py doctest

By default, ``run_tests.py`` runs all tests, including the docstring
tests.

If an individual test is failing, you can also try running it directly,
which may give you more information.

Importantly, note that the individual unit tests come in two types:

-  Simple print-and-compare scripts. These unit tests are essentially
   short example Python programs, which print out various output text.
   For a test file named ``test_XXX.py`` there will be a matching text
   file called ``test_XXX`` under the ``output`` subdirectory which
   contains the expected output. All that the test framework does to is
   run the script, and check the output agrees.
-  Standard ``unittest``- based tests. These will ``import unittest``
   and then define ``unittest.TestCase`` classes, each with one or more
   sub-tests as methods starting with ``test_`` which check some
   specific aspect of the code. These tests should not print any output
   directly.

Currently, about half of the Biopython tests are ``unittest``-style
tests, and half are print-and-compare tests.

Running a simple print-and-compare test directly will usually give lots
of output on screen, but does not check the output matches the expected
output. If the test is failing with an exception error, it should be
very easy to locate where exactly the script is failing. For an example
of a print-and-compare test, try:

.. code:: verbatim

    python test_SeqIO.py

The ``unittest``-based tests instead show you exactly which
sub-section(s) of the test are failing. For example,

.. code:: verbatim

    python test_Cluster.py

19.2??Writing tests
-------------------

Let¡¯s say you want to write some tests for a module called ``Biospam``.
This can be a module you wrote, or an existing module that doesn¡¯t have
any tests yet. In the examples below, we assume that ``Biospam`` is a
module that does simple math.

Each Biopython test can have three important files and directories
involved with it:

#. ``test_Biospam.py`` ¨C The actual test code for your module.
#. ``Biospam`` [optional]¨C A directory where any necessary input files
   will be located. Any output files that will be generated should also
   be written here (and preferably cleaned up after the tests are done)
   to prevent clogging up the main Tests directory.
#. ``output/Biospam`` ¨C [for print-and-compare tests only] This file
   contains the expected output from running ``test_Biospam.py``. This
   file is not needed for ``unittest``-style tests, since there the
   validation is done in the test script ``test_Biospam.py`` itself.

It¡¯s up to you to decide whether you want to write a print-and-compare
test script or a ``unittest``-style test script. The important thing is
that you cannot mix these two styles in a single test script.
Particularly, don¡¯t use ``unittest`` features in a print-and-compare
test.

Any script with a ``test_`` prefix in the ``Tests`` directory will be
found and run by ``run_tests.py``. Below, we show an example test script
``test_Biospam.py`` both for a print-and-compare test and for a
``unittest``-based test. If you put this script in the Biopython
``Tests`` directory, then ``run_tests.py`` will find it and execute the
tests contained in it:

.. code:: verbatim

    $ python run_tests.py     
    test_Ace ... ok
    test_AlignIO ... ok
    test_BioSQL ... ok
    test_BioSQL_SeqIO ... ok
    test_Biospam ... ok
    test_CAPS ... ok
    test_Clustalw ... ok

¡­

.. code:: verbatim

    ----------------------------------------------------------------------
    Ran 107 tests in 86.127 seconds

19.2.1??Writing a print-and-compare test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A print-and-compare style test should be much simpler for beginners or
novices to write - essentially it is just an example script using your
new module.

Here is what you should do to make a print-and-compare test for the
``Biospam`` module.

#. Write a script called ``test_Biospam.py``

   -  This script should live in the Tests directory
   -  The script should test all of the important functionality of the
      module (the more you test the better your test is, of course!).
   -  Try to avoid anything which might be platform specific, such as
      printing floating point numbers without using an explicit
      formatting string to avoid having too many decimal places
      (different platforms can give very slightly different values).

#. If the script requires files to do the testing, these should go in
   the directory Tests/Biospam (if you just need something generic, like
   a FASTA sequence file, or a GenBank record, try and use an existing
   sample input file instead).
#. Write out the test output and verify the output to be correct.

   There are two ways to do this:

   #. The long way:

      -  Run the script and write its output to a file. On UNIX
         (including Linux and Mac OS X) machines, you would do something
         like: ``python test_Biospam.py > test_Biospam`` which would
         write the output to the file ``test_Biospam``.
      -  Manually look at the file ``test_Biospam`` to make sure the
         output is correct. When you are sure it is all right and there
         are no bugs, you need to quickly edit the ``test_Biospam`` file
         so that the first line is: ¡®\ ``test_Biospam``\ ¡¯ (no quotes).
      -  copy the ``test_Biospam`` file to the directory Tests/output

   #. The quick way:

      -  Run ``python run_tests.py -g test_Biospam.py``. The regression
         testing framework is nifty enough that it¡¯ll put the output in
         the right place in just the way it likes it.
      -  Go to the output (which should be in
         ``Tests/output/test_Biospam``) and double check the output to
         make sure it is all correct.

#. Now change to the Tests directory and run the regression tests with
   ``python run_tests.py``. This will run all of the tests, and you
   should see your test run (and pass!).
#. That¡¯s it! Now you¡¯ve got a nice test for your module ready to check
   in, or submit to Biopython. Congratulations!

As an example, the ``test_Biospam.py`` test script to test the
``addition`` and ``multiplication`` functions in the ``Biospam`` module
could look as follows:

.. code:: verbatim

    from Bio import Biospam

    print "2 + 3 =", Biospam.addition(2, 3)
    print "9 - 1 =", Biospam.addition(9, -1)
    print "2 * 3 =", Biospam.multiplication(2, 3)
    print "9 * (- 1) =", Biospam.multiplication(9, -1)

We generate the corresponding output with
``python run_tests.py -g test_Biospam.py``, and check the output file
``output/test_Biospam``:

.. code:: verbatim

    test_Biospam
    2 + 3 = 5
    9 - 1 = 8
    2 * 3 = 6
    9 * (- 1) = -9

Often, the difficulty with larger print-and-compare tests is to keep
track which line in the output corresponds to which command in the test
script. For this purpose, it is important to print out some markers to
help you match lines in the input script with the generated output.

19.2.2??Writing a unittest-based test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We want all the modules in Biopython to have unit tests, and a simple
print-and-compare test is better than no test at all. However, although
there is a steeper learning curve, using the ``unittest`` framework
gives a more structured result, and if there is a test failure this can
clearly pinpoint which part of the test is going wrong. The sub-tests
can also be run individually which is helpful for testing or debugging.

The ``unittest``-framework has been included with Python since version
2.1, and is documented in the Python Library Reference (which I know you
are keeping under your pillow, as recommended). There is also `online
documentaion for
unittest <http://docs.python.org/library/unittest.html>`__. If you are
familiar with the ``unittest`` system (or something similar like the
nose test framework), you shouldn¡¯t have any trouble. You may find
looking at the existing example within Biopython helpful too.

Here¡¯s a minimal ``unittest``-style test script for ``Biospam``, which
you can copy and paste to get started:

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

In the division tests, we use ``assertAlmostEqual`` instead of
``assertEqual`` to avoid tests failing due to roundoff errors; see the
``unittest`` chapter in the Python documentation for details and for
other functionality available in ``unittest`` (`online
reference <http://docs.python.org/library/unittest.html>`__).

These are the key points of ``unittest``-based tests:

-  Test cases are stored in classes that derive from
   ``unittest.TestCase`` and cover one basic aspect of your code
-  You can use methods ``setUp`` and ``tearDown`` for any repeated code
   which should be run before and after each test method. For example,
   the ``setUp`` method might be used to create an instance of the
   object you are testing, or open a file handle. The ``tearDown``
   should do any ¡°tidying up¡±, for example closing the file handle.
-  The tests are prefixed with ``test_`` and each test should cover one
   specific part of what you are trying to test. You can have as many
   tests as you want in a class.
-  At the end of the test script, you can use

   .. code:: verbatim

       if __name__ == "__main__":
           runner = unittest.TextTestRunner(verbosity = 2)
           unittest.main(testRunner=runner)

   to execute the tests when the script is run by itself (rather than
   imported from ``run_tests.py``). If you run this script, then you¡¯ll
   see something like the following:

   .. code:: verbatim

       $ python test_BiospamMyModule.py
       test_addition1 (__main__.TestAddition) ... ok
       test_addition2 (__main__.TestAddition) ... ok
       test_division1 (__main__.TestDivision) ... ok
       test_division2 (__main__.TestDivision) ... ok

       ----------------------------------------------------------------------
       Ran 4 tests in 0.059s

       OK

-  To indicate more clearly what each test is doing, you can add
   docstrings to each test. These are shown when running the tests,
   which can be useful information if a test is failing.

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

   Running the script will now show you:

   .. code:: verbatim

       $ python test_BiospamMyModule.py
       An addition test ... ok
       A second addition test ... ok
       Now let's check division ... ok
       A second division test ... ok

       ----------------------------------------------------------------------
       Ran 4 tests in 0.001s

       OK

If your module contains docstring tests (see section
`19.3 <#section:doctest>`__), you may want to include those in the tests
to be run. You can do so as follows by modifying the code under
``if __name__ == "__main__":`` to look like this:

.. code:: verbatim

    if __name__ == "__main__":
        unittest_suite = unittest.TestLoader().loadTestsFromName("test_Biospam")
        doctest_suite = doctest.DocTestSuite(Biospam)
        suite = unittest.TestSuite((unittest_suite, doctest_suite))
        runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
        runner.run(suite)

This is only relevant if you want to run the docstring tests when you
execute ``python test_Biospam.py``; with ``python run_tests.py``, the
docstring tests are run automatically (assuming they are included in the
list of docstring tests in ``run_tests.py``, see the section below).

19.3??Writing doctests
----------------------

Python modules, classes and functions support built in documentation
using docstrings. The `doctest
framework <http://docs.python.org/library/doctest.html>`__ (included
with Python) allows the developer to embed working examples in the
docstrings, and have these examples automatically tested.

Currently only a small part of Biopython includes doctests. The
``run_tests.py`` script takes care of running the doctests. For this
purpose, at the top of the ``run_tests.py`` script is a manually
compiled list of modules to test, which allows us to skip modules with
optional external dependencies which may not be installed (e.g. the
Reportlab and NumPy libraries). So, if you¡¯ve added some doctests to the
docstrings in a Biopython module, in order to have them included in the
Biopython test suite, you must update ``run_tests.py`` to include your
module. Currently, the relevant part of ``run_tests.py`` looks as
follows:

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

Note that we regard doctests primarily as documentation, so you should
stick to typical usage. Generally complicated examples dealing with
error conditions and the like would be best left to a dedicated unit
test.

Note that if you want to write doctests involving file parsing, defining
the file location complicates matters. Ideally use relative paths
assuming the code will be run from the ``Tests`` directory, see the
``Bio.SeqIO`` doctests for an example of this.

To run the docstring tests only, use

.. code:: verbatim

    $ python run_tests.py doctest
