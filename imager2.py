import math
import numpy as np
import scipy.signal as signal
import numpy.fft as fft
from meshcage import WLMEData
import argparse
import os.path
from PIL import Image

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

  Image.fromarray(hmm).save(path)

  #hmm.tofile(path)
  #print(hmm.shape)


def main(args):
  chirp = np.load(args.tx)
  data = WLMEData.load_by_path(args.wlma)

  pad = len(chirp) * 2

  _rxs = np.load(args.rx)
  print(_rxs.shape)
  rxs = np.zeros((_rxs.shape[0], _rxs.shape[1] + pad * 2), _rxs.dtype)

  rxl = []
  #rxs = []
  for x in range(0, len(data.rx)):
    rx_name = data.rx[x].name
    rx_loc = data.rx[x].location
    rxl.append(rx_loc)
    #_data = np.load(os.path.join(args.dir, '%s.npy' % rx_name))
    #tmp = np.zeros((pad * 2 + _data.shape[0],), _data.dtype)
    #tmp[pad:pad+len(_data)] = _data

    _data = _rxs[x, :]
    rxs[x, pad:pad+len(_data)] = _data

    #rxs.append(tmp)
    print('load', rx_name)
    tx_loc = data.tx[0].location

  q = np.zeros((args.r, args.r), np.float128)
  qc = np.zeros((args.r, args.r), np.float128)
  qc[:, :] = 0.00001

  ww = 100

  a = data.cage.meshes[0].faces[0][0]
  z = 0.0

  pp = np.array([args.cx, args.cy, args.cz]) - rxl[0]

  yspace = np.linspace(pp[0] - args.w, pp[0] + args.w, q.shape[0])
  xspace = np.linspace(pp[1] - args.w, pp[1] + args.w, q.shape[1])

  w = signal.windows.triang(len(rxs))

  SA = np.zeros((len(rxs),), np.float64)

  # Determine the chirp response at all positions.
  corr = []
  for n in range(0, len(rxs)):
    corr.append(signal.correlate(rxs[n], chirp, mode='same'))

  # Iterate through spatial positions.
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
        delay = dist / args.wv * args.sps
        samp_delay = int(delay)
        SA[n] = corr[n][samp_delay]

      # From the work of Thomas Kirchner, Franz Sattler, Janek Grohl, and Lena Maier-Hein.
      # "Signed Real-Time Delay Multiply and Sum Beamforming for Multispectral Photoacoustic Imaging"
      ss = 0
      for n in range(0, len(rxs) - 1):
        s = 0
        for m in range(n + 1, len(rxs)):
          M_nm = SA[n] * SA[m]
          s += math.copysign(math.sqrt(np.abs(M_nm)), M_nm)
        ss += s
      S_dmas = ss

      px = int(xi / len(xspace) * q.shape[1])
      py = int(yi / len(yspace) * q.shape[0])
      q[py, px] += S_dmas
      qc[py, px] += 1
      
  output_image(q / qc, args.image)

ap = argparse.ArgumentParser()
ap.add_argument('--wv', type=float, required=True, help='The wave velocity.')
ap.add_argument('--sps', type=int, required=True, help='The samples per second.')
ap.add_argument('--tx', type=str, required=True, help='The input TX chirp for correlation.')
ap.add_argument('--wlma', type=str, required=True, help='The WLMA input file.')
ap.add_argument('--rx', type=str, required=True, help='The RX input streams.')
ap.add_argument('--image', type=str, required=True, help='The output image.')
ap.add_argument('--cx', type=float, required=True, help='The center X value.')
ap.add_argument('--cy', type=float, required=True, help='The center Y value.')
ap.add_argument('--cz', type=float, required=True, help='The center Z value.')
ap.add_argument('--w', type=float, required=True, help='The output spatial width and height.')
ap.add_argument('--r', type=int, required=True, help='The output resolution width and height.')
main(ap.parse_args())
