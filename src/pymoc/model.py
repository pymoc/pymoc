import numpy as np
from pymoc.modules import ModuleWrapper, Neighbor


class Model(object):
  def __init__(self):
    self.basins = []
    self.couplers = []
    self._modules = {}

  def keys(self):
    return self._modules.keys()

  def get_module(self, key):
    if not key:
      return None
    return self._modules[key]

  def validate_neighbors_input(self, neighbors):
    neighbor_keys = [n.key for n in neighbors]
    distinct_neighbor_keys = np.unique(neighbor_keys)
    if len(neighbor_keys) > len(distinct_neighbor_keys):
      raise ValueError(
          'Cannot link basins multiple times. Please check your configuration.'
      )

    for n in neighbors:
      neighbor = self.get_module(n.key)
      if not neighbor:
        raise KeyError(
            'No module present with key ' + n.key +
            ', cannot set as neighbor of ' + name
        )
      n.module_wrapper = neighbor

  def add_module(self, module, name, neighbors=None):
    neighbors = neighbors or []
    self.validate_neighbors_input(neighbors)
    module_wrapper = ModuleWrapper(module, name, neighbors)
    module_wrapper.validate_neighbor_links()
    module_wrapper.backlink_neighbors()

    if hasattr(self, module_wrapper.key):
      raise NameError(
          'Cannot use module name ' + name +
          ' because it would overwrite an existing key or model property.'
      )

    self._modules[module_wrapper.key] = module_wrapper
    if module_wrapper.module_type == 'basin':
      self.basins.append(module_wrapper)
    elif module_wrapper.module_type == 'coupler':
      self.couplers.append(module_wrapper)

    setattr(self, module_wrapper.key, module_wrapper.module)

  def new_module(
      self,
      module_class,
      module_args,
      module_name,
      neighbors=[],
  ):
    self.add_module(
        module_class(**module_args),
        module_name,
        neighbors,
    )

  def run(
      self,
      steps,
      basin_dt,
      coupler_dt,
      snapshot_start=None,
      snapshot_interval=None
  ):
    for coupler in self.couplers:
      coupler.update_coupler()

    for i in range(0, steps):
      self.timestep(i, basin_dt, coupler_dt=coupler_dt)
      if self.snapshot(i, snapshot_start, snapshot_interval):
        yield i
    return

  def timestep(self, step, basin_dt, coupler_dt=0):
    for basin in self.basins:
      basin.timestep_basin(dt=basin_dt)
    if step % coupler_dt == 0:
      for coupler in self.couplers:
        coupler.update_coupler()

  def snapshot(self, step, snapshot_start, snapshot_interval):
    return snapshot_start is not None and snapshot_interval is not None and step >= snapshot_start and (
        step-snapshot_start
    ) % snapshot_interval == 0
