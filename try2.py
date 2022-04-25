import math
import numpy as np
import scipy.signal as signal
import numpy.fft as fft
from tritopoint import line_intersect_plane
import struct
from meshcage import WLMEData

def output_image(mm, path='/run/user/1000/test.data'):
  #for x in range(0, mm.shape[1]):
  #  mm[:, x] -= np.average(mm[:, x])
  #mm.swapaxes(1, 0)

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
  '''

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

#data.rx.sort(key=lambda x: int(x.name[2:]))

rxl = []
rxs = []
for x in range(0, len(data.rx)):
  rx_name = data.rx[x].name
  rx_loc = data.rx[x].location
  rxl.append(rx_loc)
  _data = np.load('%s.npy' % rx_name)
  tmp = np.zeros((pad * 2 + _data.shape[0],), _data.dtype)
  tmp[pad:pad+len(_data)] = _data
  rxs.append(signal.hilbert(tmp))
  print('load', rx_name)

tx_loc = data.tx[0].location

wv = 1490
digital_sps = 256000

hh = 100
q = np.zeros((hh, hh), np.float128)

ww = 5

a = data.cage.meshes[0].faces[0][0]
z = 0.0

pp = a - rxl[0]

yspace = np.linspace(pp[0] - ww, pp[0] + ww, q.shape[0])
xspace = np.linspace(pp[1] - ww, pp[1] + ww, q.shape[1])

w = signal.windows.triang(len(rxs))

SA = np.zeros((len(rxs),), np.float64)

#print(a)
#exit()

for yi in range(0, len(yspace)):
  print(yi / len(yspace))
  for xi in range(0, len(xspace)):
    x = xspace[xi]
    y = yspace[yi]
    pos = np.array([x, y, z])
    upos = pos / np.linalg.norm(pos)
    tx_dist = np.linalg.norm((rxl[0] + pos) - tx_loc)
    s = 0 
    for n in range(0, len(rxs)):
      n_rxs = rxs[n]
      rx_dist = np.linalg.norm((rxl[0] + pos) - rxl[n])
      dist = rx_dist #+ tx_dist
      delay = dist / wv * digital_sps
      samp_delay = int(delay)
      #frac_delay = delay - samp_delay
      #theta = math.atan2(upos[1], upos[0])
      #w = np.exp(-1j * np.pi * math.cos(theta) * n)
      #frac_time = frac_delay / digital_sps
      v = n_rxs[samp_delay:samp_delay+len(chirp)]
      v = np.sum((v * np.conjugate(chirp)).real)
      SA[n] = v

    #SA = SA * w

    px = int(xi / len(xspace) * q.shape[1])
    py = int(yi / len(yspace) * q.shape[0])

    #'''
    ss = 0
    for n in range(0, len(rxs) - 1):
      s = 0
      for m in range(n + 1, len(rxs)):
        M_nm = SA[n] * SA[m]
        s += math.copysign(math.sqrt(np.abs(M_nm)), M_nm)
      ss += s
    S_dmas = ss
    #'''

    q[py, px] = S_dmas # np.sum(SA)
    
output_image(q)


