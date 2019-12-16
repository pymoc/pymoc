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
    # self.key = ''.join(['_'+i.lower() if i.isupper() else i for i in name]).lstrip('_').replace(' ', '')

  @property
  def module_type(self):
    return self.module.module_type

  @property
  def psi_north(self):
    coupler = self.north
    if coupler and coupler.module_type == 'coupler':
      if coupler.do_psi_bz:
        [Psi, _] = coupler.module.Psibz()
        return Psi
      else:
        return 0
    return 0

  @property
  def psi_south(self):
    coupler = self.south
    if coupler and coupler.module_type == 'coupler':
      if coupler.do_psi_bz:
        [_, Psi] = coupler.module.Psibz()
        return Psi
      else:
        return coupler.module.Psi
    return 0

  def timestep(self, dt=None):
    module = self.module
    if self.module_type == 'basin':
      wAn = self.psi_north
      wAs = self.psi_south
      wA = (wAn-wAs) * 1e6
      module.timestep(wA=wA, dt=dt, do_conv=self.do_conv)
    elif self.module_type == 'coupler':
      north = self.north
      south = self.south
      if self.do_psi_bz:
        module.update(b2=north and north.b, b1=south and south.b)
      else:
        module.update(b=north and north.b)
      module.solve()

  @property
  def b(self):
    return self[self.b_type]
    if isinstance(self.module, Column) or isinstance(self.module, Equi_Column):
      return self.module.b
    if isinstance(self.module, SO_ML):
      return self.module.bs
