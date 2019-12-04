from pymoc.modules import Module


class Model(object):
  def __init__(self):
    self.south_module = None
    self.north_module = None

  def keys(self):
    keys = []
    module = south_module
    while module:
      keys.append(module.key)
      module = module.north
    return keys

  def modules(self):
    modules = []
    module = south_module
    while module:
      modules.append(module)
      module = module.north
    return modules

  def get_module(self, key):
    if not key:
      return None

    module = south_module
    while module and module.key != key:
      module = module.north
    if module.key == key:
      return module
    return None

  def add_module(self, module, name, north_key=None, south_key=None):
    module = Module(
        module,
        name,
        north=self.get_module(north_key),
        south=self.get_module(south_key),
    )

    if module.north == self.south_module:
      self.south_module = module
    if module.south == self.north_module:
      self.north_module = module

  def get_modules_by_type(self, module_type):
    modules = []
    module = south_module
    while module:
      if module.module_type == module_type:
        modules.append(module)
      module = module.north
    return modules

  def run(self, dt, steps):
    for i in range(0, steps):
      self.timestep(dt)

  def timestep(self, dt):
    basins = self.get_modules_by_type('basin')
    couplers = self.get_modules_by_type('coupler')
    for basin in basins:
      wAn = basin.north.module.PsiS() if basin.north else 0
      wAs = basin.south.module.PsiN() if basin.south else 0
      wA = (wAn-wAs) * 1e6
      basin.timestep(wA=wA, dt=dt)
    for coupler in couplers:
      coupler.module.update(b_south=coupler.south, b_north=coupler.north)
      coupler.module.solve()
