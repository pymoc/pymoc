Vertical Buoyancy Balance
=========================

Residual-Mean Advection-Diffusion Equation
##########################################

The verical buoyancy profiles within ocean basins (e.g. the Atlantic 
Ocean interior, the North Atlantic Deepwater formation region) is 
maintained by a balance between vertical advection & diffusion, as first
described by `Munk (1966)`_. This balance is calculated in a horizontally
averaged sense, an approximation that is allowable given the flatness of
the buoyancy slopwe in the ocean interior, where isopycnal outcropping is
not a concern. As a note, this approximation makes this balance unsuitable
for modeling of the Southern Ocean, where distinct
:doc:`overturning transport <overturning-transport-closure>` &
:doc:`surface buoyancy <surface-buoyancy-balance>` formulations are utilized.

We derive our advection-diffusoing equation by starting with the Boussinesq 
continuity equation in isopycnal space

.. math::
  \partial_t\sigma + \nabla_h\cdot\left(\sigma\vec{u}\right)=-\partial_b\left(\sigma\mathcal{B}\right)

Here we define :math:`\mathcal{B}` as the diabatic buoyancy tendency, including
small scale processes (those taking place horizontal over distances less than ~1km),
and :math:`\sigma\equiv\partial_b z` as a measure of the depth range covered by a given
buyancy (:math:`b`) class that we call the  isopycnal "thickness". We can now transform
the continuity equation into depth space, by integrating from the minimum modeled buoyancy
(:math:`b_{min}`) to buoyancy :math:`b`

.. math::

  \partial_t z=-\nabla_h\cdot\int_{b_{min}}^b\sigma\vec{u}db^\prime - \sigma\mathcal{B}

We next express the diabatic forcing term via a vertical diffusive buoyancy flux (:math:`F^b`),
defined relative to the eddy diffusivity (:math:`\kappa`) and surface buoyancy flux (:math:`F_s^b`)

.. math::

  \begin{aligned}
  \sigma\mathcal{B}&=-\partial_b F^b \\
  \partial_t z&=-\nabla_h\cdot\int_{b_{min}}^b\sigma\vec{u}db^\prime + \partial_b F^b \\
  F^b&=
    \begin{cases}
      F_s^b, & b\ge b_s(x,y) \\
      -\kappa\partial_z b, & b_{bot}(x,y)\lt b\lt b_s(x,y) \\
      0, & b\le b_{bot}(x,y)
    \end{cases}
  \end{aligned}
  

.. math::

  \begin{aligned}
  \partial_t b \approx -w^\dagger\partial_zb+\partial_z(\kappa_{e\!f\!f}\partial_zb)+\mathcal{B}_s
  \end{aligned}


.. _`Munk (1966)`: https://doi.org/10.1016/0011-7471(66)90602-4