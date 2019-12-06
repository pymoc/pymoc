class Module(object):
  def __init__(
      self,
      module,
      name,
      north=None,
      south=None,
  ):
    self.module = module
    self.name = name
    self.north = north
    self.south = south
    self.key = self.name

  @property
  def module_type(self):
    return self.module.module_type

  @property
  def psi_north(self):
    if self.north and self.north.module_type == 'coupler':
      coupler = self.north.module
      if hasattr(coupler, 'Psibz') and callable(coupler, 'Psibz'):
        [Psi, _] = coupler.Psibz()
        return Psi
      else:
        return None
    return None

  @property
  def psi_south(self):
    if self.south and self.south.module_type == 'coupler':
      coupler = self.south.module
      if hasattr(coupler, 'Psibz') and callable(coupler, 'Psibz'):
        [_, Psi] = coupler.Psibz()
        return Psi
      else:
        return coupler.Psi
    return None

  def timestep(self, dt=None):
    if self.module_type == 'basin':
      wAn = self.psi_north
      wAs = self.psi_south
      wA = (wAn-wAs) * 1e6
      self.module.timestep(wA=wA, dt=dt)
      if self.north:
        self.north.south = self
      if self.south:
        self.south.north = self
    elif self.module_type == 'coupler':
      self.module.update(b_south=self.south, b_north=self.north)
      self.module.solve()
      if self.north:
        self.north.south = self
      if self.south:
        self.south.north = self
