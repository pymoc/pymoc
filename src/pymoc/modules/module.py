class Module(object):
  def __init__(
    self,
    module,
    name,
    north = None,
    south = None,
  ):
    self.module = module
    self.name = name
    self.north = north
    self.south = south
    self.key = self.name

