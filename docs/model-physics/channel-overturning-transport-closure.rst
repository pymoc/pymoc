Southern Ocean Overturning Transport Closure
============================================

For flow in reentrant zonal channels, such as the Atlantic Circumpolar Current
in the Southern Ocean, the meridional overturning circulation is maintained by
a balanance between the wind and bottom friction driven Ekman transport, and 
the eddy transport in the channel:

.. math::
  \Psi_{SO}=-\frac{\tau L_x}{\rho_of}+KL_xs

where :math:`\tau` is the wind stress magnitude, :math:`L_x` is the zonal extent
of the channel, :math:`\rho_o` is the reference density, and :math:`f` is the
Coriolis frequency. :math:`K` is a parameterized eddy diffusivity, in this case
following the `Gent & McWilliams (1990)`_ (GM) parameterization.

:math:`s=s(b)` is the isopycnal slope in the channel for buoyancy class :math:`b`

.. math::

  \begin{aligned}
  s\left(b\right)&=\begin{cases}
    \frac{z_B\left(b\right)}{L_y-y_{SO}\left(b\right)} & b >= b_{min} \\
    \frac{\tau}{K\rho_o f} & b < b_{min}
  \end{cases}
  \end{aligned}

where :math:`L_y` is the meridional extent of the channel, :math:`z_B\left(b\right)` is the 
depth of the isopycnal with density :math:`b` in the adjoining basin, and :math:`y_{SO}(b)`
is the latitude at which the isopycnal of density :math:`b` outcrops in the channel.
:math:`b_{min}` is the minimum surface buoyancy in the channel, where isopycnals with
density less than :math:`b_{min}` do not outcrop, accounting for transient cases or
equilibria in which bottom water is not formed in the channel.

.. _`Gent & McWilliams (1990)`: https://doi.org/10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2