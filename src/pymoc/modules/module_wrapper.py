import numpy as np
from pymoc.modules import Psi_SO, Psi_Thermwind, SO_ML, Column, Equi_Column


class ModuleWrapper(object):
  def __init__(self, module, name, neighbors=None):
    self.module = module
    self.name = name
    self.left_neighbors = []
    self.right_neighbors = []

    self.key = name.replace(' ', '_').lower().strip('_')
    self.do_psi_bz = hasattr(module, 'Psibz') and callable(module.Psibz)
    self.b_type = 'b' if isinstance(self.module, Column) or isinstance(
        self.module, Equi_Column
    ) else 'bs'
    self.psi = [0, 0]

    if neighbors:
      self.add_neighbors(neighbors)

  @property
  def module_type(self):
    return self.module.module_type

  def timestep_basin(self, dt=None):
    module = self.module
    wA = 0.0
    for neighbor in self.right_neighbors:
      wA += neighbor.psi[0] * 1e6
    for neighbor in self.left_neighbors:
      wA -= neighbor.psi[-1] * 1e6

    module.timestep(wA=wA, dt=dt)

  def update_coupler(self):
    module = self.module
    b1 = None
    b2 = None
    for neighbor in self.left_neighbors:
      b1 = neighbor.b
    for neighbor in self.right_neighbors:
      b2 = neighbor.b

    if self.do_psi_bz:
      module.update(b1=b1, b2=b2)
    else:
      module.update(b=b2)
    module.solve()
    self.psi = module.Psibz() if self.do_psi_bz else [module.Psi]

  @property
  def b(self):
    return getattr(self.module, self.b_type)

  @property
  def neighbors(self):
    return self.left_neighbors + self.right_neighbors

  def add_neighbor(self, new_neighbor, direction, backlink=False):
    if not direction in ['left', 'right']:
      raise ValueError(
          "Direction of a neighbor must be either 'left' or 'right'."
      )

    if new_neighbor in self.neighbors:
      raise KeyError(
          'Cannot add module ' + new_neighbor.name + ' as a neighbor of ' +
          self.name + ' because they are already coupled.'
      )

    if self.module_type == 'coupler':
      # We don't need to explicitly check that a coupler has two or fewer neighbors, as the enforcement of only one left, only one right, and only being able to point left or right implicitly enforces that condition.
      if len(self.left_neighbors) > 1 or len(self.right_neighbors) > 1:
        raise ValueError(
            'Cannot have a coupler linked in the same direction more than once. Please check your configuration.'
        )

    if direction == 'left':
      self.left_neighbors.append(new_neighbor)
    else:
      self.right_neighbors.append(new_neighbor)

    if not backlink:
      self.backlink_neighbor(new_neighbor)

  def add_neighbors(self, neighbors):
    if len(neighbors) > 0:
      for neighbor in neighbors:
        self.add_neighbor(neighbor['module'], neighbor['direction'])

  def backlink_neighbor(self, neighbor):
    direction = 'left' if neighbor in self.right_neighbors else 'right'
    neighbor.add_neighbor(self, direction, backlink=True)
