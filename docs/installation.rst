.. _installation:

Installation
############

Requirements
============

The only requirements are

- Python 2.7. or 3
- numpy_ (1.13 or later)
- scipy_

.. _numpy:  http://www.numpy.org/
.. _scipy:  http://www.scipy.org/

Instructions
============

There are several ways to install PyMOC and its dependencies
on your system. As PyMOC is intended primarily as an educational
instrument, we highly encourage you to download the source from
GitHub_ and either install the model locally or utilize our
pre-build Docker container. Regardless, instructions are provided
for installing a packaged version of PyMOC from pip for those
who are not interested in what's going on under the hood.

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

The Codeless Way: Installing with pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Hard(er) Way: Building from Source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _Docker: https://www.docker.com/products/docker-desktop
.. _`Docker Hub`: https://hub.docker.com/signup
