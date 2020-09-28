Thermal Wind Closure
====================

The circulation between basins or basin regions is represented by a Thermal
wind closure (c.f. `Nikurashin & Vallis (2012)`_, `Jansen & Nadeau (2019)`_, `Nadeau & Jansen (2020)`_ )

.. math::
  \partial_{zz}\Psi^z\left(z\right)=\frac{1}{f}\left[b_2\left(z\right)-b_1\left(z\right)\right]

Where :math:`b_1(z)` and :math:`b_2(z)` are the buoyancy profiles in the adjacent regions, :math:`f` is the coriolis parameter,
and :math:`\Psi^z(z)` is the overturning streamfunction in depth space.

Vanishing net mass flux between the columns yields the boundary conditions:

.. math::
  \begin{aligned}
  \Psi^z\left(z=0\right)&=0 \\
  \Psi^z\left(z=-H\right)&=0
  \end{aligned}

The isopycnal overturning streamfunction is then obtained by mapping the z-coordinate
streamfunction obtained from the thermal wind relation onto the buoyancy in the respective 
up-stream column:

.. math::
  \Psi^b\left(b\right) = \int_{-H}^0 \partial_z\Psi^z\left(z\right)\mathcal{H}\left[b - b_{up}\left(z\right)\right]

where :math:`\mathcal{H}` is the Heaviside step function and

.. math::
  \begin{aligned}
  b_{up}\left(z\right) = 
  \begin{cases} 
    b_1\left(z\right), & \partial_z\Psi^z\left(z\right)  > 0 \\
    b_2\left(z\right), & \partial_z\Psi^z\left(z\right)  < 0 \,.
  \end{cases}
  \end{aligned}


.. _`Nikurashin & Vallis (2012)`: https://doi.org/10.1175/JPO-D-11-0189.1
.. _`Jansen & Nadeau (2019)`: https://doi.org/10.1175/JPO-D-18-0187.1
.. _`Nadeau & Jansen (2020)`: https://doi.org/10.1175/JPO-D-20-0034.1
