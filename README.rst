A Whole Cell Model of a cocci Bacterium in Python
--------------------------------------------------------------------------------------------------------

.. |PyPI version| image:: https://img.shields.io/pypi/v/chemw.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/chemw/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/freiburgermsu/chemw/workflows/Test%20ChemW/badge.svg
   :target: https://github.com/freiburgermsu/chemw/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/chemw
   :target: https://pepy.tech/project/chemw
   :alt: Downloads


Installation
+++++++++++++

The following command installs ``WCMpy`` in a command prompt/terminal environment::
 
 pip install wcmpy

_________________

ChemMW
++++++++++++++++++

+++++++++++
__init__
+++++++++++

The data environment, in a `Python IDE <https://www.simplilearn.com/tutorials/python-tutorial/python-ide>`_, is defined: 

.. code-block:: python

  import chemw
  chem_mw = chemw.ChemMW(verbose = False, printing = True)

- *verbose* & *printing* ``bool``: specifies whether troubleshooting information or MW results will be printed, respectively.

++++++++++++++++
mass()
++++++++++++++++

The parameterized data is fitted to the Hill equation, with the following arguments and their default values:

.. code-block:: python

 chem_mw.mass(formula)

- ** ````: .


++++++++++++++++++++++++++
Accessible content
++++++++++++++++++++++++++
The ``WCMpy`` object retains numerous components that are accessible to the user: 

- ** ````: .
