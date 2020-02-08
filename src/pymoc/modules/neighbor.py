class Neighbor(object):
  def __init__(self, key, direction, module_wrapper=None):
    self.key = key
    if direction not in ['left', 'right']:
      raise ValueError('Direction must be either left or right.')
    self.direction = direction
    self.module_wrapper = module_wrapper

  def __eq__(self, other):
    return self.key == other.key and self.direction == other.direction
