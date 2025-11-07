Installation instructions
=========================

Dependencies
------------

The ``fermihalos`` package requires the use of `NumPy <https://numpy.org/>` and `SciPy <https://scipy.org/>` for its correct functionality.

Installation
------------

There are many ways to install ``fermihalos``.

Using ``pip``
*************

To install ``fermihalos`` using the ``pip`` package manager simply run:

.. code_block:: bash

    pip install fermihalos

which will install the package and its dependencies.

Using ``conda``
***************

To install ``fermihalos`` using the ``conda`` package manager simply run:

.. code_block:: bash

    conda install fermihalos

which will install the package and its dependencies.

Using the setup.py
******************

In case a manual installation is required, the first step is to clone the GitHub repository to the installation folder:

.. code_block:: bash

    git clone https://github.com/Santiq22/FermiHalos.git

Then, since it is needed to have ``setuptools`` installed, it is recomended to have it updated. To do this, run:

.. code_block:: bash

    python -m pip install --upgrade setuptools

Then, to install the package in your computer, run on the root directory of the cloned repository:

.. code_block:: bash

    python setup.py install

This command will install the package and its required dependencies into the corresponding Python environment. The package will then be available for import and use. For more details on how to clone a GitHub repository or how to install a python package see `Cloning a repository <https://docs.github.com/es/repositories/creating-and-managing-repositories/cloning-a-repository>` and `Setuptools Documentation <https://setuptools.pypa.io/en/latest/index.html>`.