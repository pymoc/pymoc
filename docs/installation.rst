.. _installation:

Installation
############

Requirements
============

The only requirements are

- Python 2.7. or 3
- numpy_ (1.13 or later)
- scipy_ (1.3 or later)


Instructions
============

There are several ways to install PyMOC and its dependencies
on your system. As PyMOC is intended primarily as an educational
instrument, we highly encourage you to download the source from
GitHub_ and either install the model locally or utilize our
pre-build Docker container. Regardless, instructions are provided
for installing a packaged version of PyMOC from pip for those
who are not interested in what's going on under the hood.

The Easy Way: Managed Installing with pip or conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PyMOC is available on the Python Package Index (PyPi), and can be
installed via pip in a terminal window. First, make sure you have
numpy_ and scipy_ installed:

.. code-block:: bash

  $ pip install numpy>=1.13 
  $ pip install scipy>=1.3 

Once these dependencies are installed, PyMOC can be installed via

.. code-block:: bash

  $ pip install py-moc

If you do not have Python3 or pip configured on your system,
instructions can be found on `Real Python`_ and in the PyPa_
documentation.

Alternatively, you may want to manage your Python environment using
the Anaconda_ distribution. Once you have installed Anaconda, you 
can follow analagous steps to those required to install via pip

.. code-block:: bash

  $ conda install -c anaconda numpy=1.13 
  $ conda install -c anaconda scipy=1.3 
  $ conda install -c anaconda py-moc


The Hard(er) Way: Building from Source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install PyMOC from source, you will need to ensure that you have
Python3, and the proper versions of numpy_ and scipy_ installed on
your system, as for the above managed installation instructions.

Rather than installing PyMoc via pip or conda, however, you will
instead directly clone the code repository using git into whatever
local directory you would like to work with the code in

.. code-block:: bash

  $ git clone https://github.com/pymoc/pymoc.git

You can then navigate into the pymoc directory, and install the model

.. code-block:: bash

  $ cd pymoc
  $ python setup.py install

If you plan on modifying the model code for your own purposes, you should
install in development mode instead

.. code-block:: bash

  $ python setup.py develop

Should you receive an error from either of the above installation commands
regarding system privileges, you will need to either install as root (via
``sudo``), or if you do not have root privileges install for only your user
by passing the ``--user`` flag to the setup script.

The Zero Configuration Way: Running in Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are just getting started with Python, and have not (or don't
want to) configure your system to run PyMOC by hand, a containerized
Docker_ environment is available to uses. A Docker container is a pre-
packaged software execution environment, with all software dependencies
installed and ready to go, meaning that end users are not required to
go through tedious system configuration to get up and running. For those
familiar with Git & GitHub, Docker serves a similar purpose, but for full
software environments rather than code. Our automated testing suite runs
in the same Docker container provided here,so for beginners containerization
is a great alternative to avoid issues surrounding package version and
configuration.

To get started with Docker_, sign up for a free account on `Docker Hub`_,
download and install the desktop software for your operating system
(available for Windows, MacOS, and Linux). Once you have installed the
desktop software, open your preferred terminal and download the latest
version of the PyMOC image:

.. code-block:: bash

    $ docker pull pymoc/pymoc:latest

Next, navigate to the directory where you want to work on the PyMOC code,
and clone the repository from GitHub:

.. code-block:: bash

    $ git clone https://github.com/pymoc/PyMOC.git

If you don't have git installed, you can also download a zipfile of the latest
source code from https://github.com/pymoc/PyMOC/archive/master.zip.

Once the code has been cloned or downloaded and unzipped, you'll need to create
and start a local container from the PyMOC image:

.. code-block:: bash

  $ docker create -it --name pymoc -v <Path to Your Code>:/pymoc/ pymoc:latest
  $ docker start pymoc

From there, you can start a bash session inside of the container via:

.. code-block:: bash

  $ docker exec -it pymoc bash

And proceed with the instructions for installing and running the model outlined in
`The Hard(er) Way`, following the dependency and system configuration steps.

.. _numpy:  http://www.numpy.org/
.. _scipy:  http://www.scipy.org/
.. _Docker: https://www.docker.com/products/docker-desktop
.. _`Docker Hub`: https://hub.docker.com/signup
.. _PyPa: https://pip.pypa.io/en/stable/installing/
.. _`Real Python`: https://realpython.com/installing-python/
.. _Anaconda: https://www.anaconda.com/distribution/
.. _GitHub: https://www.github.com/pymoc/PyMOC/
