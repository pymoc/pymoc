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
.. image:: https://img.shields.io/badge/license-MIT-informational
  :target: https://github.com/pymoc/PyMOC/blob/master/LICENSE

PyMOC is a simple, modular suite of python ocean column models for
use in studying the Meridional Overturning Circulation (MOC). The 
MOC plays a critical role in the uptake and redistribution of heat
and carbon by the ocean, and as such both mediates and is governed
by shifts in the climate regime. As such, understanding of the MOC
is crucial to understanding of climate change.

The model suite consists of several independent modules representing
various ocean regions and dynamics. Specifically, there are modules
for calculating the advective-diffusive buoyancy balance in zonally
constrained ocean basins (such as the Atlantic Ocean) in both
transient and equilibrium states, in re-entrant surface channel flow
(such as in the shallow Southern Ocean), for calculating a thermal
wind balance between basins, and for calculating residual wind and
eddy driven circulation in a deep channel. These modules may be 
coupled to study an arbitrary circulatory structure (such as the AMOC).

The intended audiences for this model are educators and students
of the geophysical sciences. While the goal is to provide an accessible
model appropriate for newcomers to geophysical modeling, the physics
captured in PyMOC are robust enough to support basic research as well.

Configuration and execution of the PyMOC suite requires little
administrative knowledge on the technical end. All modules are written
in pure Python, and the only core dependencies are the NumPy and SciPy
libraries. If configuration of your base system environment is undesirable,
a preconfigured Docker container has been made available with all required
software libraries pre-installed. Furthermore, a goal of the development
team is to keep PyMOC well tested, stable, and maintainable to reduce
pain to the end user. Further details on installation, configuration,
contribution, and issue reporting is available in the documentation.

.. include:: contents.rst

