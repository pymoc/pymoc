.. pymoc documentation master file, created by
   sphinx-quickstart on Fri Aug 30 10:54:58 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyMOC: Python Meridional Overturning Circulation Model
======================================================
.. image:: https://img.shields.io/badge/code-GitHub-informational
  :target: https://www.github.com/pymoc/PyMOC
.. image:: https://circleci.com/gh/pymoc/pymoc/tree/master.svg?style=shield
  :target: https://circleci.com/gh/pymoc/pymoc/tree/master
.. image:: https://api.codeclimate.com/v1/badges/b03ff00b5c86d7afc364/test_coverage
  :target: https://codeclimate.com/github/pymoc/PyMOC/test_coverage
.. image:: https://api.codeclimate.com/v1/badges/b03ff00b5c86d7afc364/maintainability
  :target: https://codeclimate.com/github/pymoc/PyMOC/maintainability
.. image:: https://badge.fury.io/py/py-moc.svg
   :target: https://badge.fury.io/py/py-moc
.. image:: https://img.shields.io/badge/license-MIT-informational
  :target: https://github.com/pymoc/PyMOC/blob/master/LICENSE

PyMOC is a suite of python modules to build simple "toy" models for ocean's Meridional Overturning Circulation (MOC).

The model suite consists of several independent modules representing various ocean regions and dynamics. Specifically, there are modules for calculating the 1-D advective-diffusive buoyancy tendencies averaged over ocean basins, given the net isopycnal transports in and out of the column. The isopycnal transports are computed as diagnostic relationships, with modules to calculate the wind- and eddy-driven residual circulation in a southern-ocean-like re-entrant channel, as well as the geostrophic exchange between different basins or basin regions. These modules may be coupled to study the circulation in a wide range of different configurations.

The intended audiences for this model are researchers, educators and students in the geophysical sciences. The goal is to provide an accessible model appropriate for both experts and newcomers to geophysical modeling, but with physics that reflect the current state of our theoretical understanding of the deep ocean overturning circulation.

Configuration and execution of the PyMOC suite requires relatively little technical knowledge or computational resources. All modules are written in pure Python, and the only core dependencies are the NumPy and SciPy libraries. If configuration of your base system environment is undesirable, a preconfigured Docker container has been made available with all required software libraries pre-installed.

Anyone is more than welcome to contibute to the development of PyMOC, but is asked to adhere to the goal of keeping PyMOC well tested, stable, maintainable, and documented. Further details on installation, configuration, contribution, and issue reporting are provided in this documentation.

.. include:: contents.rst

