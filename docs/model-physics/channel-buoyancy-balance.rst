Southern Ocean Surface Buoyancy Balance
=======================================

Computation of the :ref:`Southern Ocean overturning transport<Southern Ocean overturning transport>`
requires knowledges of the surface buoyancy structure in the Southern Ocean. This surface buoyancy
can either be prescribed or computed via a horizontal residual advection diffusion equation for
the surface layer of the Southern Ocean:

.. math::

  \partial_t b_{SO} = -v^\dagger_{SO}\partial_y b_{SO} + \partial_y\left(K_s\partial_y b_{SO}\right) + \mathcal{B}_{SO}

Here :math:`\mathcal{B}_{SO}` is the surface forcing, which can be prescribed either via a buoyancy
restoring to a specified profile or via a fixed flux condition, or via a combination of the two,
with a fixed flux applied over part of the latitude range and restoring applied elsewhere 
(as in `Jansen & Nadeau (2019)`_). :math:`K_s` is a horizontal eddy diffusivity in the mixed layer, and
:math:`v^\dagger_{SO}` is the residual meridional advection defined as

.. math::

  v^\dagger_{SO} = \frac{\Psi_{SO}}{L_x h}

where :math:`h` is the Southern Ocean mixed layer depth.

.. _`Jansen & Nadeau (2019)`: https://doi.org/10.1175/JPO-D-18-0187.1

.. I commented the text below out because this is only one specific possible way to configure
.. the module
.. In agreement with `Jansen et al. (2018)`_ a fixed buoyancy flux is applied in the
.. region of Antarctic Bottom Water formation, constrained to those latitudes south of
.. :math:`y_{fixed}` where ice formation takes place. North of this latitude, a surface
.. restoring boundary condition is applies, resulting in a surface buoyancy forcing that
.. can be described as:

.. .. math::

..  \begin{aligned}
..  \mathcal{B}_{SO}&=\begin{cases}
..    -\frac{\mathcal{F}_{SO}}{h} & y\le y_{fixed} \\
..    \frac{v_{pist}}{h}\left(b_{SO,0}-b_{SO}\right) & y > y_{fixed}
..  \end{cases}
..  \end{aligned}

.. where :math:`h` is the Southern Ocean mixed layer depth, :math:`\mathcal{F}_{SO}`
.. is the specified buoyancy flux in the Antarctic Bottom Water formation region south
.. of :math:`y_{fixed}`, :math:`v_{pist}` is a piston velocity, :math:`b_{SO}` is the
.. meridional surface buoyancy profile, and :math:`b_{SO,0}` is a prescribed surface
.. buoyancy profile.
