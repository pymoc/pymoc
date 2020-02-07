from pymoc.modules import ModuleWrapper


class Model(object):
  def __init__(self):
    self.south_module = None
    self.north_module = None
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
    module_wrapper = ModuleWrapper(
        module,
        name,
        north=north,
        south=south,
        do_conv=do_conv,
    )

    if hasattr(self, module_wrapper.key):
      raise NameError(
          'Cannot use module name ' + name +
          ' because it would overwrite an existing key.'
      )

    self._modules[module_wrapper.key] = module_wrapper

    if north == self.south_module:
      self.south_module = module_wrapper
    if south == self.north_module:
      self.north_module = module_wrapper

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
      north_key=None,
      south_key=None,
      do_conv=False
  ):
    self.add_module(
        module_class(**module_args),
        module_name,
        north_key=north_key,
        south_key=south_key,
        do_conv=do_conv
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
