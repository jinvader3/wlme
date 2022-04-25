'''
  This is a collection of meshes for which intersection tests are made
  using rays. This module also handles loading the mesh into the 
  collection.
'''
import numpy as np
import tritopoint
import json
import scipy.signal as signal
from tritopoint import line_intersect_tri
import math
import numpy.fft as fft
import argparse

def get_zadoff_chu(nzc, q, u):
  out = np.zeros((nzc,), np.complex128)
  cf = nzc % 2
  for n in range(0, nzc):
    samp = np.exp(0-1j * ((np.pi * u * n * (n + cf + 2 * q)) / nzc))
    out[n] = samp
  return out

def get_chirp(sps, time, f0, f1, p2=1.0):
  theta = 0.0
  samps = math.floor(sps * time)
  fd = f1 - f0
  out = np.zeros((samps,), np.float64)
  for x in range(0, samps):
    ifreq = (x ** p2) / samps * fd + f0
    out[x] = math.cos(theta)
    theta += (ifreq * math.pi * 2) / sps
  return out

def main(args):
  digital_sps = 16000
  #chirp = get_chirp(digital_sps, 0.1, 500, 6000)
  chirp = get_zadoff_chu(353, 1, 7).real
  #chirp = np.random.random(1000)

  np.save('chirp.npy', chirp, False)

if __name__ == '__main__':
  ap = argparse.ArgumentParser()
  main(ap.parse_args())
