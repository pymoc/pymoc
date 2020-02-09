import numpy as np
from pymoc.modules import Psi_SO, Psi_Thermwind, SO_ML, Column, Equi_Column, Neighbor


class ModuleWrapper(object):
  def __init__(self, module, name, neighbors=None):
    self.module = module
    self.name = name
    self.neighbors = neighbors or []

    self.key = name.replace(' ', '_').lower().strip('_')
    self.do_psi_bz = hasattr(module, 'Psibz') and callable(module.Psibz)
    self.b_type = 'b' if isinstance(self.module, Column) or isinstance(
        self.module, Equi_Column
    ) else 'bs'
    self.psi = [0, 0]

  @property
  def module_type(self):
    return self.module.module_type

  def timestep_basin(self, dt=None):
    module = self.module
    wA = 0.0
    for neighbor in self.neighbors:
      if neighbor.direction == 'right':
        wA += neighbor.module_wrapper.psi[0] * 1e6
      else:
        wA -= neighbor.module_wrapper.psi[-1] * 1e6

    module.timestep(wA=wA, dt=dt)

  def update_coupler(self):
    module = self.module
    b1 = None
    b2 = None
    for neighbor in self.neighbors:
      if neighbor.direction == 'left':
        b1 = neighbor.module_wrapper and neighbor.module_wrapper.b
      else:
        b2 = neighbor.module_wrapper and neighbor.module_wrapper.b

    if self.do_psi_bz:
      module.update(b1=b1, b2=b2)
    else:
      module.update(b=b2)
    module.solve()
    self.psi = module.Psibz() if self.do_psi_bz else [module.Psi]

  @property
  def b(self):
    return getattr(self.module, self.b_type)

  def validate_neighbors(self, backlink=False):
    if not backlink:
      for neighbor in self.neighbors:
        linked_neighbor_keys = [
            n.key for n in neighbor.module_wrapper.neighbors
        ]
        if self.key in linked_neighbor_keys:
          raise KeyError(
              'Cannot add module ' + self.name + ' as a neighbor of ' +
              neighbor.module_wrapper.name +
              ' because they are already coupled.'
          )

    if self.module_type == 'coupler':
      # We don't need to explicitly check that a coupler has two or fewer neighbors, as the enforcement of only one left, only one right, and only being able to point left or right implicitly enforces that condition.
      left_neighbor_count = sum(n.direction == 'left' for n in self.neighbors)
      right_neighbor_count = sum(
          n.direction == 'right' for n in self.neighbors
      )
      if left_neighbor_count > 1 or right_neighbor_count > 1:
        raise ValueError(
            'Cannot have a coupler linked in the same direction more than once. Please check your configuration.'
        )

  def backlink_neighbor_links(self):
    for neighbor in self.neighbors:
      neighbor.module_wrapper.neighbors.append(
          Neighbor(
              self.key,
              'left' if neighbor.direction == 'right' else 'right',
              module_wrapper=self
          )
      )
      neighbor.module_wrapper.validate_neighbors(backlink=True)
