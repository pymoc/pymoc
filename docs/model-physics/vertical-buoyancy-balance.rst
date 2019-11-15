Vertical Buoyancy Balance
=========================

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
  
We can now integrate :math:`\partial_t z` zonally and meridionally across the basin, assuming zero zonal
flux (consistent with a basin bounded by continental landmasses):

.. math::
  
  \begin{aligned}
  \partial_t\int_{x_1}^{x_2}\!\int_{y_1}^{y_2}zdxdy &=
    \int_{x_1}^{x_2}\!\int_{b_{min}}^b \sigma vdb^\prime \left.dx\right|_{y_2}^{y_1} +
    \partial_b\int_{y_1}^{y_2}\!\int_{x_1}^{x_2}F^bdxdy \\
    \Psi\left(y, b\right) &\equiv -\int_{x_1}^{x_2}\!\int_{b_{min}}^b \sigma vdb^\prime dx \\
    \left<X\right> &\equiv \frac{1}{A_o}\int_{x_1}^{x_2}\!\int_{y_1}^{y_2}Xdxdy \\
    A_o &\equiv \frac{1}{A_o}\int_{x_1}^{x_2}\!\int_{y_1}^{y_2}1dxdy \\
    \partial_t\left<z\right> &= \frac{1}{A_o}\left[\Psi\left(y_2, b\right) - \Psi\left(y_1, b\right)\right] + \partial_b\left<F^b\right>
  \end{aligned}

where :math:`\Psi\left(y, b\right)` is the meridional overturning streamfunction in bouyancy space, and
:math:`\left<X\right>` is the horizontal spatial mean of an arbitrary quantity :math:`X`. Taking the spatial
average of the isopycnal thickness :math:`\sigma`, we can derive the identity
:math:`\left.\left<\sigma\right>=\partial_b\left<z\right>\right|_{b=b\left(\left<z\right>\right)}` where
:math:`b\left(\left<z\right>\right)` is the buoyancy of the isopycnal with average depth :math:`z`.

This allows us to rearrange the terms of our equation, so that:

.. math::

  \begin{aligned}
  \left<\sigma\right>^{-1}\partial_t\left<z\right> &= \left.-\partial_b\left<b\right>\right|_{\left<z\right>} \\
  \left.-\partial_t\left<b\right>\right|_{\left<z\right>} &= -\frac{1}{A_o}\left[\Psi\left(b\left(\left<z\right>\right), y_2\right) - \Psi\left(b, y_1\right)\right]\partial_{\left<z\right>}b-\partial_{\left<z\right>}\left<F^b\right> \\
  \left.-\partial_t\left<b\right>\right|_{\left<z\right>} &= -w^\dagger\partial_{\left<z\right>}b-\partial_{\left<z\right>}\left<F^b\right>
  \end{aligned}

where we've defined :math:`w^\dagger` as the residual upwelling, due to convergence of volume flux along isopycnals
below :math:`b\left(\left<z\right>\right)`.

If we re-express the buyoyancy forcing in terms of a diabatic surface forcing :math:`\mathcal{B}_s` and an effective
vertical diffusivity :math:`k_{eff}\equiv\left(\frac{A_1}{A_o}\right)^2\kappa` where :math:`A_1` is the area of non-outcropping 
isopycnals, and drop the spatial averaging notation, our equations can be reduced to:

.. math::

  \begin{aligned}
  \partial_t b \approx -w^\dagger\partial_zb+\partial_z(\kappa_{e\!f\!f}\partial_zb)+\mathcal{B}_s
  \end{aligned}

Where the vertical buoyancy profile is maintained by an advective-diffusive balance relative to average isopycnal depths,
and accounting for surface forcing.  

.. _`Munk (1966)`: https://doi.org/10.1016/0011-7471(66)90602-4