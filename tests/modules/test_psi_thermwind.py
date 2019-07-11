import pytest
import sys
sys.path.append('/pymoc/src/modules')
from psi_thermwind import Psi_Thermwind

class TestPsi_Thermwind(object):
  def test(self):
    assert 1