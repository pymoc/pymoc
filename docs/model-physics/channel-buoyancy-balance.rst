Southern Ocean Surface Buoyancy Balance
=======================================

Computation of the :ref:`Southern Ocean overturning transport<Southern Ocean overturning transport>`
requires knowledges of the surface buoyancy structure in the Southern Ocean. This surface buoyancy
can either be prescribed or computed via a horizontal residual advection diffusion equation for
the surface layer of the Southern Ocean:

.. math::

  \partial_t b_{SO} = -v^\dagger_{SO}\partial_y b_{SO} + \partial_y\left(K_s\partial_y b_{SO}\right) + \mathcal{B}_{SO}

Here :math:`\mathcal{B}_{SO}` is the surface forcing, which can be prescribed either via a buoyancy
restoring to a specified profile or via a fixed flux condition, or via a combination of the two (as in `Jansen & Nadeau (2019)`_). :math:`K_s` is a horizontal eddy diffusivity in the mixed layer, and :math:`v^\dagger_{SO}` is the residual meridional advection defined as

.. math::

  v^\dagger_{SO} = \frac{\Psi_{SO}}{L_x h}

where :math:`h` is the Southern Ocean mixed layer depth.

.. _`Jansen & Nadeau (2019)`: https://doi.org/10.1175/JPO-D-18-0187.1
