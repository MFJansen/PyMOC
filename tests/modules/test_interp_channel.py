import pytest
import sys
import numpy as np
sys.path.append('/pymoc/src/modules')
from interp_channel import Interpolate_channel

@pytest.fixture(scope="module", params=[
  {
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 81)),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 81)),
    'y': 1e6
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 81)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
  },
  {
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'z': np.asarray(np.linspace(-4.0e3, 0, 81)),
    'bn': np.linspace(0.03, -0.001, 81),
    'bs': np.linspace(0.05, 0.04, 51)
  },
])
def interp_channel_config(request):
  return request.param

class TestInterpolate_channel(object):
  def test_interp_channel_init(self, interp_channel_config):
    if not 'y' in interp_channel_config or not isinstance(interp_channel_config['y'], np.ndarray) or not len(interp_channel_config['y']):
      with pytest.raises(TypeError) as yinfo:
        Interpolate_channel(**interp_channel_config)
      assert(str(yinfo.value) == "y needs to be numpy array providing grid levels")
      return
    if not 'z' in interp_channel_config or not isinstance(interp_channel_config['z'], np.ndarray) or not len(interp_channel_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Interpolate_channel(**interp_channel_config)
      assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
      return
    if not 'bs' in interp_channel_config or not isinstance(interp_channel_config['bs'], np.ndarray):
      with pytest.raises(TypeError) as bsinfo:
        Interpolate_channel(**interp_channel_config)
      assert(str(bsinfo.value) == "('bs', 'needs to be either function, numpy array, or float')")
      return
    if not 'bn' in interp_channel_config or not isinstance(interp_channel_config['bn'], np.ndarray):
      with pytest.raises(TypeError) as bninfo:
        Interpolate_channel(**interp_channel_config)
      assert(str(bninfo.value) == "('bn', 'needs to be either function, numpy array, or float')")
      return

    interp_channel = Interpolate_channel(**interp_channel_config)
    for k in ['z', 'y']:
      assert hasattr(interp_channel, k)

    bs = interp_channel.make_func(interp_channel_config['bs'], 'bs', interp_channel.y)
    for y in interp_channel.y:
      assert(interp_channel.bs(y) == bs(y))

    bn = interp_channel.make_func(interp_channel_config['bn'], 'bn', interp_channel.z)
    for z in interp_channel.z:
      assert(interp_channel.bn(z) == bn(z))
