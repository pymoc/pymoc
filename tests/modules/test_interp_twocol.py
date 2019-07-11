import pytest
import sys
sys.path.append('/pymoc/src/modules')
from interp_twocol import Interpolate_twocol

class TestInterpolate_twocol(object):
  def test(self):
    assert 1