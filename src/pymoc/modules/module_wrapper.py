import numpy as np
from pymoc.modules import Column, Equi_Column


class ModuleWrapper(object):
  r"""
  Module Wrapper

  Instances of this class wrap individual phsyical modules 
  (e.g. advective-diffusive columns, thermal wind streamfunctions)
  for use in the model driver. The module wrapper is responsible for 
  timestepping & updating modules, & communication between neighboring
  modules.
  
  By convention geographic north & east are defined as 
  being to the "right" or "rightward", while geographic south & west 
  are defined as being the the "left" or "leftward". Modules are categorized
  as "basins" with time evolving density profiles (e.g. columns, mixed layers),
  or "couplers" with streamfunctions that maintain the dynamical balance between
  neighboring basins (e.g. thermal wind streamfunction, SO transport).
  """
  def __init__(
      self,
      module,
      name,
      left_neighbors=None,
      right_neighbors=None,
  ):
    r"""
    Parameters
    ----------

    module : module class instance
             Physics module being wrapped
    name : string
           Name of the module being wrapped (e.g. Atlantic Basin, AMOC)
    left_neighbors : list; optional
                     List of modules to the "left" of the module being wrapped.
                     Couplers may have at most one left neighbor.
    right_neighbors : list; optional
                      List of modules to the "right" of the module being wrapped.
                      Couplers may have at most one right neighbor.
    """
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

    if left_neighbors or right_neighbors:
      self.add_neighbors(
          left_neighbors=left_neighbors, right_neighbors=right_neighbors
      )

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
    if self.module_type != 'coupler':
      raise TypeError('Cannot use update_coupler on non-coupler modules.')

    module = self.module
    b1 = self.left_neighbors[0].b if len(self.left_neighbors) > 0 else None
    b2 = self.right_neighbors[0].b if len(self.right_neighbors) > 0 else None

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

  def add_left_neighbor(self, new_neighbor, backlinking=False):
    self.validate_neighbor_uniqueness(new_neighbor)
    self.validate_coupler_neighbor_direction('left')
    self.left_neighbors.append(new_neighbor)
    if not backlinking:
      self.backlink_neighbor(new_neighbor)

  def add_right_neighbor(self, new_neighbor, backlinking=False):
    self.validate_neighbor_uniqueness(new_neighbor)
    self.validate_coupler_neighbor_direction('right')
    self.right_neighbors.append(new_neighbor)
    if not backlinking:
      self.backlink_neighbor(new_neighbor)

  def add_neighbors(self, left_neighbors=None, right_neighbors=None):
    if len(left_neighbors) > 0:
      for neighbor in left_neighbors:
        self.add_left_neighbor(neighbor)

    if len(right_neighbors) > 0:
      for neighbor in right_neighbors:
        self.add_right_neighbor(neighbor)

  def validate_neighbor_uniqueness(self, neighbor):
    if neighbor in self.neighbors:
      raise KeyError(
          'Cannot add module ' + neighbor.name + ' as a neighbor of ' +
          self.name + ' because they are already coupled.'
      )

  def validate_coupler_neighbor_direction(self, direction):
    if self.module_type == 'coupler':
      # We don't need to explicitly check that a coupler has two or fewer neighbors, as the enforcement of only one left, only one right, and only being able to point left or right implicitly enforces that condition.
      if len(self.left_neighbors) > 0 and direction == 'left' or len(
          self.right_neighbors
      ) > 0 and direction == 'right':
        raise ValueError(
            'Cannot have a coupler linked in the same direction more than once. Please check your configuration.'
        )

  def backlink_neighbor(self, neighbor):
    if neighbor in self.right_neighbors:
      neighbor.add_left_neighbor(self, backlinking=True)
    else:
      neighbor.add_right_neighbor(self, backlinking=True)
