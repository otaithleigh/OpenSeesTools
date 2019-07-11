Section Analysis
++++++++++++++++

The |SectionAnalysis| class provides methods for inspecting OpenSees fiber
sections.


.. autoclass:: OpenSeesTools.SectionAnalysis
    :members:


Example
=======

.. code-block::

    >>> from OpenSeesTools import SectionAnalysis, fourFiberSectionGJ
    >>> def createSection():
    ...     ops.uniaxialMaterial('Elastic', 1, 29000.0)
    ...     fourFiberSectionGJ(1, 1, area=10.0, Iy=144.0, Iz=94.0, GJ=11000.0)
    >>> SA = SectionAnalysis(createSection)
    >>> SA.printMaterialInfo()
       Material |   # Fibers |   Area |      Iz |   Iy
    ------------+------------+--------+---------+------
              1 |          4 |     10 | 93.9999 |  144
    ------------+------------+--------+---------+------
          Total |          4 |     10 | 93.9999 |  144


.. |SectionAnalysis| replace:: :class:`OpenSeesTools.SectionAnalysis`
