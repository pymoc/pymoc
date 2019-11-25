Southern Ocean Surface Buoyancy Balance
=======================================

Understanding of the :ref:`Southern Ocean overturning transport<Southern Ocean overturning transport>`
requires the caclulation of the surface buoyancy structure overlying the channel.

In agreement with `Jansen et al. (2018)`_ a fixed buoyancy flux is applied in the
region of Antarctic Bottom Water formation, constrained to those latitudes south of
:math:`y_{fixed}` where ice formation takes place. North of this latitude, a surface
restoring boundary condition is applies, resulting in a surface buoyancy forcing that
can be described as:

.. math::

  \begin{aligned}
  \mathcal{B}_{SO}&=\begin{cases}
    -\frac{\mathcal{F}_{SO}}{h} & y\le y_{fixed} \\
    \frac{v_{pist}}{h}\left(b_{SO,0}-b_{SO}\right) & y > y_{fixed}
  \end{cases}
  \end{aligned}

where :math:`h` is the Southern Ocean mixed layer depth, :math:`\mathcal{F}_{SO}`
is the specified buoyancy flux in the Antarctic Bottom Water formation region south
of :math:`y_{fixed}`, :math:`v_{pist}` is a piston velocity, :math:`b_{SO}` is the
meridional surface buoyancy profile, and :math:`b_{SO,0}` is a prescribed surface
buoyancy profile.

:math:`\mathcal{B}_{SO}` is then used as the forcing term in the advective-difussive
balance that maintains the meridional buoyancy structure at the surface of the Southern
Ocean

.. math::

  \partial_tb_{SO} = -v^\dagger_{SO}\partial_yb_{SO} + \partial_y\left(K_s\partial_yb_{SO}\right) + \mathcal{B}_{SO}

where :math:`K_s` is a horizontal eddy diffusivity in the mixed layer, and 
:math:`v^\dagger_{SO}` is the residual meridional advection defined as

.. math::

  v^\dagger_{SO} = \frac{\Psi_{SO}}{L_xh}

.. _`Jansen et al. (2018)`: https://doi.org/10.1175/JCLI-D-17-0797.1