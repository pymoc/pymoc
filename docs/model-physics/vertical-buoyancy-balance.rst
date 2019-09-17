Vertical Buoyancy Balance - Residual-Mean Advection-Diffusion Equation 
======================================================================

The verical buoyancy profiles within ocean basins (e.g. the Atlantic 
Ocean interior, the North Atlantic Deepwater formation region) is 
maintained by a balance between vertical advection & diffusion. This
balance is calculated in a horizontally averaged sense, an approximation
that is allowable given the flatness of the buoyancy slopwe in the ocean
intetior, where isopycnal outcropping is not a concern. As a note, this
approximation makes this balance unsuitable for modeling of the Southern
Ocean, where distinct :doc:`overturning transport <overturning-transport-closure>` &
:doc:`surface buoyancy <surface-buoyancy-balance>` formulations are utilized.

.. math::

  \begin{aligned}
  \partial_t b \approx -w^\dagger\partial_zb+\partial_z(\kappa_{e\!f\!f}\partial_zb)+\mathcal{B}_s
  \end{aligned}

