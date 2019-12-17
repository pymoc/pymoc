from pymoc.modules import Module


class Model(object):
  def __init__(self):
    self.south_module = None
    self.north_module = None
    self.basins = []
    self.couplers = []
    self._modules = {}

  def keys(self):
    keys = []
    module = self.south_module
    while module:
      keys.append(module.key)
      module = module.north
    return keys

  def modules(self):
    modules = []
    module = self.south_module
    while module:
      modules.append(module)
      module = module.north
    return modules

  def get_module(self, key):
    if not key:
      return None

    return self._modules[key]
    # module = self.south_module
    # while module and module.key != key:
    #   module = module.north
    # if module and module.key == key:
    #   return module
    # return None

  def add_module(
      self, module, name, north_key=None, south_key=None, do_conv=False
  ):
    north = self.get_module(north_key)
    south = self.get_module(south_key)
    module = Module(
        module,
        name,
        north=north,
        south=south,
        do_conv=do_conv,
    )
    self._modules[module.key] = module

    if north == self.south_module:
      self.south_module = module
    if south == self.north_module:
      self.north_module = module

    if module.module_type == 'basin':
      self.basins.append(module)
    elif module.module_type == 'coupler':
      self.couplers.append(module)

  def get_modules_by_type(self, module_type):
    modules = []
    module = self.south_module
    while module:
      if module.module_type == module_type:
        modules.append(module)
      module = module.north
    return modules

  def run(
      self,
      steps,
      basin_dt,
      coupler_dt,
      snapshot_start=None,
      snapshot_interval=None
  ):
    for i in range(0, steps):
      print(str(i) + '/' + str(steps))
      self.timestep(i, basin_dt, coupler_dt=coupler_dt)
      if self.snapshot(i, snapshot_start, snapshot_interval):
        yield i
    return

  def timestep(self, step, basin_dt, coupler_dt=0):
    for basin in self.basins:
      basin.timestep(dt=basin_dt)
    if step % coupler_dt == 0:
      for coupler in self.couplers:
        coupler.timestep()

  def snapshot(self, step, snapshot_start, snapshot_interval):
    return snapshot_start is not None and snapshot_interval is not None and step >= snapshot_start and (
        step-snapshot_start
    ) % snapshot_interval == 0
