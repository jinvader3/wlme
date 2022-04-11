import math
import numpy as np
import scipy.signal as signal
import numpy.fft as fft
from tritopoint import line_intersect_plane
import struct
from meshcage import WLMEData

DIGITAL_SPS = 192000 / 10

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

chirp = np.load('chirp.npy')
data = WLMEData.load_by_path('./testchamber.wlmadata')

pad = len(chirp) * 2

rxs = []
for x in range(0, len(data.rx)):
  _data = np.load('rx%s.npy' % x)
  tmp = np.zeros((pad * 2 + _data.shape[0],), _data.dtype)
  tmp[pad:pad+len(_data)] = _data
  rxs.append(tmp)
  print('load', np.max(rxs[x]))

rxlocs = []
#rxlocs.append(np.array([0, 0, 0]))

for i in range(1, len(data.rx)):
  rx = data.rx[i]
  #loc = rx.location - data.rx[0].location
  rxlocs.append(rx.location)

ox = rxlocs[0][0]
oy = rxlocs[0][1]
oz = rxlocs[0][2]
w = 200

'''
for y in range(0, mm.shape[0]):
  ay = y / mm.shape[0] * h - h * 0.5 + oy
  print(ay, y / mm.shape[0])
  for x in range(0, mm.shape[1]):
    ax = x / mm.shape[1] * w - w * 0.5 + ox
    e = np.sum(signal.correlate(k, chirp, mode='same'))
    mm[y][x] = e
'''

mc = np.ones((w, w), np.float128)
mm = np.zeros((w, w), np.float128)

mm[:, :] *= 0.00001

most = []

origin = np.array([ox, oy, oz])

cy = 0
while True:
  urv = np.random.random(3) * 2 - 1
  rv = urv * w - w * 0.5 + origin
  k = beam(rv[0], rv[1], rv[2], pad, rxs)
  #k = k[0:4000]
  e = np.abs(signal.correlate(k, chirp, mode='same'))

  urv = urv / np.linalg.norm(urv)
  
  y = int((math.atan2(urv[2], urv[0]) + np.pi) / (np.pi * 2) * w)
  xp = w / e.shape[0]
  for i in range(0, e.shape[0]):
    v = e[i]
    x = int(xp * i)
    #point = urv * (i / DIGITAL_SPS * 1480)
    #x = int(point[0] * 1 + mm.shape[0] * 0.5)
    #y = int(point[1] * 1 + mm.shape[1] * 0.5)
    try:
      if x >= 0 and y >= 0:
        mm[x, y] += v
        mc[x, y] += 1.0
    except IndexError:
      break
    
  cy += 1
  if cy % 100 == 0:
    print('@', mc)
    output_image(mm / mc)

  '''
  most.append((e, tuple(rv)))
  i += 1
  if i % 100 == 0:
    most.sort(key=lambda x: x[0], reverse=True)
    while len(most) > 1000:
      most.pop()
    print('=')
    with open('test.plot', 'w') as fd:
      for v in most:
        fd.write('%s %s %s\n' % (v[1][0], v[1][1], v[1][2]))
  '''

# Get average.
mm = mm / mc


