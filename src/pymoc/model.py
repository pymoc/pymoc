import numpy as np
from pymoc.modules import ModuleWrapper, Neighbor


class Model(object):
  def __init__(self):
    self.basins = []
    self.couplers = []
    self._modules = {}

  def keys(self):
    keys = []
    module_wrapper = self.south_module
    while module:
      keys.append(module.key)
      module = module.north
    return keys

  def modules(self):
    modules = []
    module_wrapper = self.south_module
    while module:
      modules.append(module)
      module = module.north
    return modules

  def get_module(self, key):
    if not key:
      return None
    return self._modules[key]

  def add_module(self, module, name, neighbors=None):
    neighbors = neighbors or []
    if len(np.unique([n.key
                      for n in neighbors])) < len([n.key for n in neighbors]):
      raise ValueError(
          'Cannot link basins multiple times. Please check your configuration.'
      )

    for n in neighbors:
      neighbor = self.get_module(n.key)
      if not neighbor:
        raise KeyError('No module present with key ' + n.key)
      n.module_wrapper = neighbor

    module_wrapper = ModuleWrapper(module, name, neighbors)
    if module_wrapper.module_type == 'coupler':
      if len(neighbors) > 2:
        raise ValueError(
            'Streamfunctions cannot connect more than two basins. Please check your configuration.'
        )
      if sum(n.direction == 'left' for n in neighbors
             ) > 1 or sum(n.direction == 'right' for n in neighbors) > 1:
        raise ValueError(
            'Cannot have a coupler linked in the same direction more than once. Please check your configuration.'
        )

    for neighbor in neighbors:
      if module_wrapper.key in [
          n.key for n in neighbor.module_wrapper.neighbors
      ]:
        raise KeyError(
            'Cannot add module ' + module_wrapper.name + ' as a neighbor of ' +
            neighbor.module_wrapper.name + ' because they are already coupled.'
        )

    for neighbor in neighbors:
      neighbor.module_wrapper.neighbors.append(
          Neighbor(
              module_wrapper.key,
              'left' if neighbor.direction == 'right' else 'right',
              module_wrapper=module_wrapper
          )
      )
      if neighbor.module_wrapper.module_type == 'coupler':
        if len(neighbor.module_wrapper.neighbors) > 2:
          raise ValueError(
              'Streamfunctions cannot connect more than two basins. Please check your configuration.'
          )
        if sum(
            n.direction == 'left' for n in neighbor.module_wrapper.neighbors
        ) > 1 or sum(
            n.direction == 'right' for n in neighbor.module_wrapper.neighbors
        ) > 1:
          raise ValueError(
              'Cannot have a coupler linked in the same direction more than once. Please check your configuration.'
          )

    if hasattr(self, module_wrapper.key):
      raise NameError(
          'Cannot use module name ' + name +
          ' because it would overwrite an existing key.'
      )

    if module_wrapper.module_type == 'coupler' and len(
        module_wrapper.neighbors
    ) > 2:
      raise ValueError(
          'Streamfunctions cannot connect more than two basins. Please check your configuration.'
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

  def get_modules_by_type(self, module_type):
    modules = []
    module_wrapper = self.south_module
    while module_wrapper:
      if module_wrapper.module_type == module_type:
        modules.append(module_wrapper)
      module_wrapper = module.north
    return modules

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
