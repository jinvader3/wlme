import numpy as np
import scipy.signal as signal
import numpy.fft as fft
from tritopoint import line_intersect_plane
import struct

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
    fvector = fvector / np.linalg.norm(fvector)

    dists = []

    p1 = self.ants[0]
    dists.append(0)

    for ant_ndx in range(1, len(self.ants)):
      ant_loc = self.ants[ant_ndx]
      #print('ant_loc=%s fvector=%s p1=%s' % (ant_loc, fvector, p1))
      p, d = line_intersect_plane(ant_loc, fvector, p1, fvector)
      dists.append(np.linalg.norm(p - ant_loc) * d)

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
 
chirp = np.load('chirp.npy')

dc = DistDelayCalc()

from meshcage import WLMEData

data = WLMEData.load_by_path('./testchamber.wlmadata')

#t = np.zeros((3,), np.float64)
#for rx in data.rx:
#  t += rx.location
#t = t / len(data.rx)
#t[0] = data.rx[0].location[0]
#t[1] = data.rx[0].location[1]
#t[2] = data.rx[0].location[2]

rxs = []
for x in range(0, len(data.rx)):
  _data = np.load('rx%s.npy' % x)
  tmp = np.zeros((len(chirp) * 2 + _data.shape[0],), _data.dtype)
  tmp[len(chirp):len(chirp)+len(_data)] = _data
  rxs.append(tmp)
  print('load', np.max(rxs[x]))

for rx in data.rx:
  dc.add_ant(rx.location)

theta_steps = 20

chop_s = 0
chop_l = 1000000
chop_k = 100

psz = len(rxs[0]) - len(chirp) * 2
mm = np.zeros((theta_steps, chop_l // chop_k), np.float128)

#input('chop_s=%s chop_l=%s len(chirp)=%s' % (
#  chop_s, chop_l, len(chirp)
#))

for theta_ndx in range(0, theta_steps):
  theta = np.pi * 2 / theta_steps * theta_ndx
  fvector = np.array([np.cos(theta), np.sin(theta), 0.0])
  sdelay = dc.calc_from_as_samples(fvector, 1480, 192000 * 10)

  #print(sdelay)

  m = np.zeros((psz,), np.float128)

  base = len(chirp)

  for x in range(0, len(rxs)):
    off = int(sdelay[x] + base)
    assert(off >= 0)
    rx = rxs[x][off:off + psz]
    m += rx

  out = m
  #print('out', np.max(out), np.min(out), out.shape)
  out = signal.correlate(m, chirp[::-1], mode='valid')
  #out = fftseg(out, 128)
  #out = fft.rfft(out)
  mm[theta_ndx, :] = np.abs(out[chop_s:chop_s+chop_l:chop_k])

  #out = rx
  #tmp[int(len(out) * 0.02):] = 0
  #out = fft.irfft(tmp)

  #print(sdelay)
  print('theta=%s' % (theta,))
#np.save('test.npy', mm)

#mm = np.load('test.npy')

'''
with open('./plot.txt', 'w') as fd:
  print(x / mm.shape[1])
  for x in range(0, mm.shape[1]):
    fd.write('\n')
    for y in range(0, mm.shape[0]):
      fd.write('%s\n' % mm[y][x])
'''

un = np.unique(mm)
_un = {}

for x in range(0, len(un)):
  print('un', x / len(un))
  _un[un[x]] = x
 
for x in range(0, mm.shape[1]):
  print('un', x / mm.shape[1])
  for y in range(0, mm.shape[0]): 
    mm[y][x] = _un[mm[y][x]]

hmm = np.zeros(mm.shape, np.uint16)
w = mm.shape[1]
for x in range(0, mm.shape[1] // w):
  _mm = mm[:, x*w:x*w+w]
  mm_max = np.max(_mm)
  mm_min = np.min(_mm)
  mm_delta = mm_max - mm_min
  print('mm_delta=%s mm_max=%s mm_min=%s' % (mm_delta, mm_max, mm_min))
  tmp = (_mm - mm_min) / mm_delta * 0xffff
  tmp[tmp > 0xffff] = 0xffff
  mm[:, x*w:x*w+w] = tmp

hmm[:] = mm
hmm.tofile('test.data')
print(hmm.shape)
