Thermal Wind Closure
====================

The circulation between basins or basin regions, such as between the Atlantic interior and the
North Atlantic deep water formation region, is represented by a Thermal
wind closure

.. math::
  \partial_{zz}\Psi^z_N\left(z\right)=\frac{1}{f}\left[b_B\left(z\right)-b_N\left(z\right)\right]

Where :math:`b_B(z)` is the basin interior buoyancy, :math:`b_N(z)` is the northern
basin buoyancy, :math:`f` is the coriolis parameter, and :math:`\Psi^z_N(z)` is 
the **zonal** overturning streamfunction in the northern basin.

`Nikurashin & Vallis (2012)`_ argued that in the North Atlantic, eastward currents are subducted
at the eastern boundary and propagate westward towards the western boundary.
By mass balance the zonal overturning in the north of the basin is then matched by a
meridional overturning in the western boundary. This argument, which is adopted in PyMOC
justifies the use of the thermal wind relation for the meridonal overturning in the North
Atlantic---that is for the overturning between the columns representing the basin interior
and the northern deep water formation region. 

Vanishing net mass flux between the columns yields the boundary conditions:

.. math::
  \begin{aligned}
  \Psi^z_N\left(z=0\right)&=0 \\
  \Psi^z_N\left(z=H\right)&=0
  \end{aligned}

.. Given a depth :math:`z_o` where :math:`b_N\left(z_o\right)=b_B\left(z_o\right)`, we impose
.. the further conditions that:

.. .. math::
..  \begin{aligned}
..  b_N\left(z\right)&=b_N & z >= z_o \\
..  \Psi_N\left(z\right)&=0 & z < z_o
..  \end{aligned}
..
.. where :math:`b_N` is the uniform buoyancy over the convective depth range.
.. Notice that the above only applies for the equi_colum module, which solves for both the
.. overturning streamfunction and buoyancy profles at once - albeit under special conditions. 

Finally, the isopycnal overturning streamfunction is obtained by mapping the z-coordinate
streamfunction obtained from the thermal wind relation onto the buoyancy in the respective 
up-stream column:

.. math::
  \Psi^b\left(b\right) = \int_{-H}^0 \partial_z\Psi\left(z\right)\mathcal{H}\left[b - b_{up}\left(z\right)\right]

where :math:`\mathcal{H}` is the Heaviside step function and

.. math::
  \begin{aligned}
  b_{up}\left(z\right) = 
  \begin{cases} 
    b_N\left(z\right), & \partial_z\Psi\left(z\right)  > 0 \\
    b_B\left(z\right), & \partial_z\Psi\left(z\right)  < 0 \,.
  \end{cases}
  \end{aligned}


.. _`Nikurashin & Vallis (2012)`: https://doi.org/10.1175/JPO-D-11-0189.1
