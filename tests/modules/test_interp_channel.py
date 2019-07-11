import pytest
import sys
sys.path.append('/pymoc/src/modules')
from interp_channel import Interpolate_channel

class TestInterpolate_channel(object):
  def test(self):
    assert 1