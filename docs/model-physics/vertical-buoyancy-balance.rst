Vertical Buoyancy Balance
=========================

The verical buoyancy profiles within ocean basins (e.g. the Atlantic 
Ocean interior, the North Atlantic Deepwater formation region) are 
maintained by a balance between advection & diffusion. If we  assume that isopycnals
(i.e. surface of constant buoyancy) are relatively flat, horizonatal advection
has a small contribution to the buoyancy tendencies, and the dominant balance 
is between vertical advection and diffusion, as first described by `Munk (1966)`_..
Following `Jansen and Nadeau (2019)`_, we here relax the the strict requirement for
flat isopycnals by implementing a residual-advection diffusion equation, which is
horizontally averaged along isopycnals, such that horizontal buoayncy advection
vanishes by construction. 

The residual advection-diffusion equation can be derived by starting with the Boussinesq 
continuity equation in isopycnal space

.. math::
  \partial_t\sigma + \nabla_h\cdot\left(\sigma\vec{u}\right)=-\partial_b\left(\sigma\mathcal{B}\right)

Here we define :math:`\mathcal{B}` as the diabatic buoyancy tendency, including
small scale diapycnal mixing processes that cannot be resolved by the model.
:math:`\sigma\equiv\partial_b z` is the isopycnal "thickness"---a measure of the depth range covered by a given
buyancy (:math:`b`) class. We can now transform
the continuity equation into depth space, by integrating from the minimum modeled buoyancy
(:math:`b_{min}`) to buoyancy :math:`b`

.. math::

  \partial_t z=-\nabla_h\cdot\int_{b_{min}}^b\sigma\vec{u}db^\prime - \sigma\mathcal{B}

We next express the diabatic forcing term via a vertical diffusive buoyancy flux (:math:`F^b`),
defined relative to the eddy diffusivity (:math:`\kappa`) and surface buoyancy flux (:math:`F_s^b`)

.. math::

  \begin{aligned}
  \sigma\mathcal{B}&=-\partial_b F^b \\
  F^b&=
    \begin{cases}
      F_s^b, & b\ge b_s(x,y) \\
      -\kappa\partial_z b, & b_{bot}(x,y)\lt b\lt b_s(x,y) \\
      0, & b\le b_{bot}(x,y)
    \end{cases}
  \end{aligned}
  
We can now integrate :math:`\partial_t z` zonally and meridionally across the basin, assuming zero
flux across teh zonal boundaries (consistent with a basin bounded by continental landmasses):

.. math::
  
  \partial_t\int_{x_1}^{x_2}\!\int_{y_1}^{y_2}zdxdy =
    \int_{x_1}^{x_2}\!\int_{b_{min}}^b \sigma vdb^\prime \left.dx\right|_{y_2}^{y_1} +
    \partial_b\int_{y_1}^{y_2}\!\int_{x_1}^{x_2}F^bdxdy
  
Defining

.. math::

  \begin{aligned} 
    \Psi\left(y, b\right) &\equiv -\int_{x_1}^{x_2}\!\int_{b_{min}}^b \sigma vdb^\prime dx \\
    \langle X \rangle &\equiv \frac{1}{A_o}\int_{x_1}^{x_2}\!\int_{y_1}^{y_2}Xdxdy \\
    A_o &\equiv \frac{1}{A_o}\int_{x_1}^{x_2}\!\int_{y_1}^{y_2}1dxdy \\
  \end{aligned}

where :math:`\Psi\left(y, b\right)` is the meridional overturning streamfunction in bouyancy space, and
:math:`\langle X \rangle` is the horizontal spatial mean of an arbitrary quantity :math:`X`, we get 

.. math::
  \begin{aligned} 
    \partial_t\langle z \rangle &= \frac{1}{A_o}\left[\Psi\left(y_2, b\right) - \Psi\left(y_1, b\right)\right] + \partial_b\langle F^b \rangle
  \end{aligned}
  

Dividing by :math:`\langle \sigma \rangle`, we can derive an equation for :math:`b(\langle z\rangle)`:

.. math::

  \begin{aligned}
  \left.-\partial_t b \right|_{\langle z\rangle} &= -\frac{1}{A_o}\left[\Psi\left(b\left(\langle z \rangle\right), y_2\right) - \Psi\left(b, y_1\right)\right]\partial_{\langle z\rangle}b-\partial_{\langle z\rangle }\langle F^b \rangle \\
    &\equiv -w^\dagger\partial_{\langle z\rangle}b-\partial_{\langle z\rangle }\langle F^b\rangle
  \end{aligned}

where we've defined :math:`w^\dagger` as the residual upwelling due to convergence of volume flux along isopycnals
below :math:`b\left(\langle z\rangle\right)`.

If we re-express the buyoyancy forcing in terms of a diabatic surface forcing :math:`\mathcal{B}_s` and an effective
vertical diffusivity :math:`\kappa_{eff}\equiv\left(\frac{A_1}{A_o}\right)^2\kappa` where :math:`A_1` is the area of the non-incropping 
part of the isopycnal surace, and drop the spatial averages in our notation, the equation can be reduced to:

.. math::

  \begin{aligned}
  \partial_t b \approx -w^\dagger\partial_z b+\partial_z(\kappa_{e\!f\!f}\partial_z b)+\mathcal{B}_s
  \end{aligned}

This final equation, solved by the column model module, is equivalent to the 1-dimensional vertical advection-diffusion equation, except for some re-interpretation of the variables and the introduction of an effective diapycnal diffusivity that accounts for the reduced area of isopycnals incropping into the bottom. 

.. _`Munk (1966)`: https://doi.org/10.1016/0011-7471(66)90602-4
.. _`Jansen and Nadeau (2019)`: https://doi.org/10.1175/JPO-D-18-0187.1
