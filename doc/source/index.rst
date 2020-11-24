.. OpenSees-Tools documentation master file, created by
   sphinx-quickstart on Thu Jul 11 14:12:43 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OpenSees-Tools's documentation!
==========================================

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    basic
    SectionAnalysis


Setting which OpenSeesPy to use
-------------------------------

Since there's more than one way to install and use OpenSeesPy, it's necessary to
make sure openseestools is using the same one that the rest of your code uses.
By default, openseestools will try to import a locally-built version of
OpenSeesPy, falling back to pip-installed OpenSeesPy. You can set this by
overriding the ``opensees`` variable.

For example, if your locally-built ``opensees.so`` isn't on your Python path
when you import openseestools, you can tell openseestools to use your version
instead:

.. code:: python

   # Import openseestools, which falls back to pip-installed openseespy
   import openseestools

   # Import your locally-built module and tell openseestools to use that
   import opensees as ops
   openseestools.opensees = ops
