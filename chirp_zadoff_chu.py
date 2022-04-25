'''
  This is a collection of meshes for which intersection tests are made
  using rays. This module also handles loading the mesh into the 
  collection.
'''
import numpy as np
import scipy.signal as signal
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

def main(args):
  chirp = get_zadoff_chu(args.nzc, args.q, args.u).real
  np.save(args.out, chirp, False)

if __name__ == '__main__':
  ap = argparse.ArgumentParser()
  ap.add_argument('--nzc', type=int, required=True)
  ap.add_argument('--q', type=int, required=True)
  ap.add_argument('--u', type=int, required=True)
  ap.add_argument('--out', type=str, required=True)
  main(ap.parse_args())
