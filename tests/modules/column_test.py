import sys
import os
import pytest
src_path = os.path.realpath(os.path.dirname(__file__))
src_path = '/'.join(src_path.split('/')[0:-2]) + '/src/modules'
sys.path.append(src_path)
from column import Column

class TestColumn(object):
  def test_init(self):
    c = Column()
    assert hasattr(c, 'z')