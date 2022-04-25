import math
import numpy as np
import scipy.signal as signal
import numpy.fft as fft
from tritopoint import line_intersect_plane
import struct
from meshcage import WLMEData

DIGITAL_SPS = 48000

def beam(ax, ay, az, pad, rxs):
  bdx = ax - rxlocs[0][0]
  bdy = ay - rxlocs[0][1]
  bdz = az - rxlocs[0][2]
  s = rxs[0].shape[0] - pad - 2000
  k = np.zeros((s,), np.float64)
  for i in range(0, len(rxlocs)):
    dx = ax - rxlocs[i][0]
    dy = ay - rxlocs[i][1]
    dz = az  - rxlocs[i][2]
    d0 = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    d1 = math.sqrt(bdx ** 2 + bdy ** 2 + bdz ** 2)
    d = d0 - d1
    ds = (d / 1480) * DIGITAL_SPS
    si = int(pad+ds)
    k += rxs[i][si:si+s]
  return k

def output_image(mm, path='/run/user/1000/test.data'):
  #for x in range(0, mm.shape[1]):
  #  mm[:, x] -= np.average(mm[:, x])
  #mm.swapaxes(1, 0)

  #'''
  un = np.unique(mm)
  _un = {}

  for x in range(0, len(un)):
    print('un', x / len(un))
    _un[un[x]] = x
   
  for x in range(0, mm.shape[1]):
    print('un', x / mm.shape[1])
    for y in range(0, mm.shape[0]): 
      mm[y][x] = _un[mm[y][x]]
  #'''

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
  hmm.tofile(path)
  print(hmm.shape)

chirp = signal.hilbert(np.load('chirp.npy'))
data = WLMEData.load_by_path('./testchamber.wlmadata')

pad = len(chirp) * 2

data.rx.sort(key=lambda x: int(x.name[2:]))

rxs = []
for x in range(0, len(data.rx)):
  rx_name = data.rx[x].name
  _data = np.load('%s.npy' % rx_name)
  tmp = np.zeros((pad * 2 + _data.shape[0],), _data.dtype)
  tmp[pad:pad+len(_data)] = _data
  rxs.append(signal.hilbert(tmp))
  print('load', rx_name)

w = 20
h = 10000
mc = np.ones((w, h), np.float128)
mm = np.zeros((w, h), np.float128)
mm[:, :] *= 0.00001

wave_velocity = 1490
ula_d = 1

cy = 0

thetas = np.linspace(0, np.pi, w)
for y in range(0, len(thetas)):
  theta = thetas[y]

  # The array design frequency.
  ula_f = wave_velocity / ula_d

  print('ula_f=%s y=%s theta=%s' % (ula_f, y / len(thetas), theta))

  ww = np.zeros((len(rxs),), np.complex128)
  for i in range(0, len(rxs)):
    ww[i] = np.exp(-1j * np.pi * math.cos(theta) * i)

  # hann looked interesting.. almost correct

  ww = ww * signal.windows.hann(len(ww))
  #ww *= signal.windows.chebwin(len(ww), at=20, sym=False)

  k = np.zeros((rxs[0].shape[0],), np.complex128)
  for i in range(0, len(rxs)):
    k += rxs[i] * ww[i]

  #k = np.abs(signal.correlate(k, chirp, mode='same'))
  k = np.abs(k)

  for i in range(0, len(k)):
    mm[y, int(i / len(k) * h)] += k[i]
    mc[y, int(i / len(k) * h)] += 1

mm = mm / mc

q = np.zeros((mm.shape[1], 2), np.float64)

qq = np.zeros((400, 400), np.float64)
qqc = np.zeros((400, 400), np.float64)

qqc[:, :] = 0.000001

for x in range(0, mm.shape[1]):
  s = np.array([0, 0], np.float64)
  l = mm.shape[0] // 1
  for y in range(0, l):
    t = (y / l) * np.pi * 2
    v = np.array([np.cos(t), np.sin(t)])
    s += v * np.abs(mm[y, x])
  
  a = math.atan2(s[1], s[0])

  #a += x / mm.shape[1] * np.pi * 0.3

  rx = int(np.cos(a) * (x / mm.shape[1]) * 200) + 200
  ry = int(np.sin(a) * (x / mm.shape[1]) * 200) + 200
  qq[ry, rx] += np.linalg.norm(s)
  qqc[ry, rx] += 1

  q[x, 0] = a
  q[x, 1] = np.linalg.norm(s)

b = 400
q[:, 0] = signal.convolve(q[:, 0], np.ones((b,), q.dtype) / b, mode='same')

with open('test.plot', 'w') as fd:
  for x in range(0, q.shape[0]):
    fd.write('%s %s\n' % (q[x, 0], q[x, 1]))

#for y in range(0, mm.shape[0]):
#  mm[y, :] -= np.average(mm[y, :])

output_image(qq / qqc)
