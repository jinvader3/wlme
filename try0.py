import numpy as np
import scipy.signal as signal
import numpy.fft as fft
from tritopoint import line_intersect_plane
import struct

chirp = np.load('chirp.npy')

rxs = []
for x in range(0, 16):
  rxs.append(np.load('rx%s.npy' % x))

class DistDelayCalc:
  def __init__(self):
    self.ants = []
    self.mdist = 0  

  def add_ant(self, loc):
    self.ants.append(loc)
    locmag = np.linalg.norm(loc)
    if locmag > self.mdist:
      self.mdist = locmag

  def calc_from_as_samples(self, fvector, wave_velocity, digital_sps):
    dist = np.array(self.calc_from(fvector))
    return (dist / wave_velocity) * digital_sps

  def calc_from(self, fvector):
    fvector = -fvector
    fvector = fvector / np.linalg.norm(fvector)
    # The fvector should point inward toward the array.
    p1 = -fvector * (self.mdist + 1.0)
    # Now, p1 and fvector define a plane. Intersect this
    # plane from each of the antennas and calculate the
    # delay from that.
    dists = []

    for ant_ndx in range(0, len(self.ants)):
      ant_loc = self.ants[ant_ndx]
      p = line_intersect_plane(ant_loc, -fvector, p1, fvector)
      dists.append(np.linalg.norm(p - ant_loc))

    base = dists[0]
    for ant_ndx in range(0, len(self.ants)):
      dists[ant_ndx] = dists[ant_ndx] - base

    return dists

def fftseg(data, segsz):
  zeros = segsz - (len(data) % segsz)
  if zeros == 0:
    _data = data
  else:
    # zero extend the data
    _data = np.zeros((len(data) + zeros,), data.dtype)
    _data[0:len(data)] = data
  
  hsegsz = segsz // 2 + 1
  out = np.zeros((_data.shape[0] // segsz * hsegsz,), _data.dtype)


  for x in range(0, len(_data) // segsz):
    segdata = _data[x*segsz:x*segsz+segsz]
    out[x*hsegsz:x*hsegsz+hsegsz] = fft.rfft(segdata)

  return out
 
dc = DistDelayCalc()

from meshcage import WLMEData

data = WLMEData.load_by_path('./100meterpond.wlmadata')

t = np.zeros((3,), np.float64)
for rx in data.rx:
  t += rx.location
t = t / len(data.rx)

for rx in data.rx:
  dc.add_ant(rx.location - t)

theta_steps = 1000
mm = np.zeros((theta_steps, len(chirp)), np.float128)
for theta_ndx in range(0, theta_steps):
  theta = np.pi * 2 / theta_steps * theta_ndx
  fvector = np.array([np.cos(theta), np.sin(theta), 0])
  sdelay = dc.calc_from_as_samples(fvector, 1480, 192000)

  m = np.zeros((mm.shape[1],), np.float128)

  base = 400
  limit = mm.shape[1]

  for x in range(0, len(rxs)):
    off = int(sdelay[x] + base)
    assert(off >= 0)
    rx = rxs[x][off:off + mm.shape[1]]
    m += rx

  out = signal.correlate(m, chirp, mode='full')
  #out = fftseg(out, 128)
  #out = fft.rfft(out)
  mm[theta_ndx, 0:len(out)] = np.abs(out[0:mm.shape[1]])

  #out = rx
  #tmp[int(len(out) * 0.02):] = 0
  #out = fft.irfft(tmp)

  #print(sdelay)
  print('theta=%s len(out)=%s' % (theta, len(out)))
#np.save('test.npy', mm)

#mm = np.load('test.npy')

hmm = np.zeros(mm.shape, np.uint16)

mm_max = np.max(mm)
mm_min = np.min(mm)
mm_delta = mm_max - mm_min
mm_delta *= 0.2
tmp = (mm - mm_min) / mm_delta * 0xffff
tmp[tmp > 0xffff] = 0xffff
hmm[:] = tmp
hmm.tofile('test.data')
print(hmm.shape)
