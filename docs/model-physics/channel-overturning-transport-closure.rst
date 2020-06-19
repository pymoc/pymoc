.. _Southern Ocean overturning transport: 

Southern Ocean Overturning Transport Closure
============================================

For the Southern Ocean, the meridional overturning circulation
is assumed to be  maintained by a balance between the wind-driven
Ekman transport and the eddy transport (following `Marshall & Radko, 2003`_
and `Marshall & Zanna, 2014`_):

.. math::
  \Psi_{SO}=-\frac{\tau L_x}{\rho_of}+L_x K s

where :math:`\tau` is the wind stress magnitude, :math:`L_x` is the zonal extent
of the channel, :math:`\rho_o` is the reference density, and :math:`f` is the
Coriolis frequency. :math:`K` is an eddy "diffusivity",
following the `Gent & McWilliams (1990)`_ (GM) parameterization.

:math:`s=s(b)` is the isopycnal slope in the Southern Ocean:

.. math::

  \begin{aligned}
  s\left(b\right)&=\begin{cases}
    \frac{z_B\left(b\right)}{L_y-y_{SO}\left(b\right)} & b >= b_{min} \\
    \frac{\tau}{K\rho_o f} & b < b_{min}
  \end{cases}
  \end{aligned}

where :math:`L_y` is the meridional extent of the Southern Ocean channel, :math:`z_B\left(b\right)`
is the  depth of the isopycnal with density :math:`b` in the adjoining basin, and :math:`y_{SO}(b)`
is the latitude at which the isopycnal of density :math:`b` outcrops in the channel.
:math:`b_{min}` is the minimum surface buoyancy in the channel. Isopycnals with
density less than :math:`b_{min}` do not outcrop at the surface in the Southern Ocean.
Since we approximate the interior of the Southern Ocean to be adiabatic, the residual circulation
has to vanish on isopycnals that do not outcrop, and their slope is set to be such that the eddy-driven
transport exactly cancels the wind-driven Ekman transport.

.. _`Marshall & Radko, 2003`: https://doi.org/10.1175/1520-0485(2003)033<2341:RSFTAC>2.0.CO;2
.. _`Marshall & Zanna, 2014`: https://doi.org/10.1175/JCLI-D-13-00344.1
.. _`Gent & McWilliams (1990)`: https://doi.org/10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2
