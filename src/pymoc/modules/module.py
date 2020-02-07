from datetime import datetime
from pymoc.modules import Psi_SO, Psi_Thermwind, SO_ML, Column, Equi_Column


class Module(object):
  def __init__(
      self,
      module,
      name,
      north=None,
      south=None,
      do_conv=False,
  ):
    if north:
      north.south = self
    if south:
      south.north = self
    self.module = module
    self.name = name
    self.north = north
    self.south = south
    self.do_conv = do_conv
    self.key = name.replace(' ', '_').lower().strip('_')
    self.do_psi_bz = hasattr(module, 'Psibz') and callable(module.Psibz)
    self.b_type = 'b' if isinstance(self.module, Column) or isinstance(
        self.module, Equi_Column
    ) else 'bs'
    self.psi = [0, 0]
    # self.key = ''.join(['_'+i.lower() if i.isupper() else i for i in name]).lstrip('_').replace(' ', '')

  @property
  def module_type(self):
    return self.module.module_type

  def timestep_basin(self, dt=None):
    module = self.module
    wA = 0.0
    if self.north:
      wA += self.north.psi[0] * 1e6
    if self.south:
      wA -= self.south.psi[-1] * 1e6
    module.timestep(wA=wA, dt=dt, do_conv=self.do_conv)

  def update_coupler(self, dt=None):
    module = self.module
    north = self.north
    south = self.south
    if self.do_psi_bz:
      module.update(b2=north and north.b, b1=south and south.b)
    else:
      module.update(b=north and north.b)
    module.solve()
    self.psi = module.Psibz() if self.do_psi_bz else [module.Psi]

  @property
  def b(self):
    return getattr(self.module, self.b_type)
