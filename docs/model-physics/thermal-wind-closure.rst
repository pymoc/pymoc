Thermal Wind Closure
====================

Circulation between basins, such as between the Atlantic interior and the
North Atlantic deep water formation region, is represented by a Thermal
wind closure

.. math::
  \partial_{zz}\Psi^z_N\left(z\right)=\frac{1}{f}\left[b_B\left(z\right)-b_N\left(z\right)\right]

Where :math:`b_B(z)` is the basin interior buoyancy in depth space, :math:`b_N(z)` is the northern
basin buoyancy in depth space, :math:`f` is the coriolis parameter, and :math:`\Psi^z_N(z)` is 
the **zonal** overturning streamfunction in the northern basin.

`Nikurashin & Vallis (2012)`_ argued that in the North Atlantic, eastward currents are subducted
at the eastern boundary and propagate westward, where they meet the direct westward currents.
Currents reaching the western boundary of the basin are then balanced by meridional western
boundary currents. This structure is consistent the large scale potential vorticity balance in
the ocean gyres, where basin interior flow returns in a frictional western boundary current 
(e.g. the Gulf Stream).

Thus, we are able to utilize the thermal wind balanced zonal flow at the boundary between the
basin interior and the northern deep water formation region as an approximate measure of the
meridional overturning circulation between those basins.

The thermal wind balance is subject to the boundary conditions:

.. math::
  \begin{aligned}
  \Psi^z_N\left(z=0\right)&=0 \\
  \Psi^z_N\left(z=H\right)&=0
  \end{aligned}

Given a depth :math:`z_o` where :math:`b_N\left(z_o\right)=b_B\left(z_o\right)`, we impose
the further conditions that:

.. math::
  \begin{aligned}
  b_N\left(z\right)&=b_N & z >= z_o \\
  \Psi_N\left(z\right)&=0 & z < z_o
  \end{aligned}

where :math:`b_N` is the uniform buoyancy over the convective depth range.

Finally, the overturning streamfunction is rewritten isopycnal space, in order to reconcile
flow between the model domains, as

.. math::
  \Psi^b\left(b\right) = \int_{-H}^0 \partial_z\Psi\left(z\right)\mathcal{H}\left[b - b_{up}\left(z\right)\right]

by computing upwind density classes

.. math::
  \begin{aligned}
  b_{up}\left(z\right) = 
  \begin{cases} 
    b_N\left(z\right), & \partial_z\Psi\left(z\right)  > 0 \\
    b_B\left(z\right), & \partial_z\Psi\left(z\right)  < 0
  \end{cases}
  \end{aligned}

where :math:`\mathcal{H}` is the Heaviside step function.

.. _`Nikurashin & Vallis (2012)`: https://doi.org/10.1175/JPO-D-11-0189.1