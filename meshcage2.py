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
import random
import argparse

class Ray:
  def __init__(self, origin, vector):
    self.u = origin
    self.v = vector

class RX:
  def __init__(self, name, location, sensitivity=1.0):
    self.name = name
    self.location = location
    self.sensitivity = sensitivity

class TX:
  def __init__(self, name, location, power):
    self.name = name
    self.location = location
    self.power = power

class WLMEData:
  def __init__(self):
    self.rx = []
    self.tx = []
    self.cage = MeshCage()

  def load_by_path(path):
    self = WLMEData()

    with open(path, 'r') as fd:
      data = json.loads(fd.read())
    
    for rx in data['rx']:
      self.rx.append(RX(
        rx['name'],
        np.array(rx['location']), 
        rx['sensitivity']
      ))

    for tx in data['tx']:
      self.tx.append(TX(
        tx['name'],
        np.array(tx['location']),
        tx['power']
      ))

    for mesh in data['meshes']:
      self.cage.add_mesh(Mesh(mesh['polygons'], mesh['meta']))

    return self

class Mesh:
  def __init__(self, faces, meta):
    _faces = []
    for face in faces:
      _face = []
      for ndx in range(3):
        _face.append(
          np.array([face[ndx][0], face[ndx][1], face[ndx][2]])
        )
      _faces.append(_face)
    self.faces = _faces
    self.meta = meta

  def ray_test(self, ray, dbg=False):
    for face in self.faces:
      a = face[0]
      b = face[1]
      c = face[2]
      ip = line_intersect_tri(
        ray.u, ray.v, a, b, c, dbg
      )

      if ip is not None:
        # I can only think of one case where there would be
        # more than one intersection. That is where the
        # intersection is directly upon an edge shared by two
        # or more faces. I'm just going to go with speed over
        # accuracy here and pick the first face.
        return ip, face

    return None

class MeshCage:
  def __init__(self):
    self.meshes = []

  def add_mesh(self, mesh):
    print('added mesh')
    self.meshes.append(mesh)

  def ray_test(self, ray, dbg=False):
    points = []

    for mesh in self.meshes:
      res = mesh.ray_test(ray, dbg)
      if res is not None:
        # Return point of intersection and distance
        # from the ray origin.
        ip = res[0]
        dist = np.linalg.norm(ip - ray.u)
        face = res[1]
        points.append((ip, dist, face, mesh.meta))

    # No intersections.
    if len(points) == 0: return None  

    # Return the first ray to mesh intersection.
    points.sort(key=lambda x: x[1])
    
    return points[0][0], points[0][2], points[0][3]

  def ray_bounces(self, ray, max_bounce_count):
    b = None
    for i in range(max_bounce_count):
      if i > 0:
        yield ray.u, ray.v, b[2]
      b = data.cage.ray_single_bounce(ray)
      if b is None: 
        return
      ray = Ray(b[0] + b[1] * 0.001, b[1])

  def ray_single_bounce(self, ray, dbg=False):
    res = self.ray_test(ray, dbg)
    if res is None:
      # If the mesh cage is not a closed up volume
      # then the ray escapes and no intersection is
      # found.
      return None

    # Let us compute the reflected angle.
    ip, face, meta = res

    a = face[0] 
    b = face[1]
    c = face[2]

    # Get vector representing strike to surface of face.
    d = ip - ray.u
    d = d / np.linalg.norm(d)

    # Get face normal.
    n = np.cross(b - a, c - a)
    # Normalize face normal to unit length one.
    n = n / np.linalg.norm(n)
    #n = Vector3(-n.x, -n.y, -n.z)

    # Compute reflection vector.
    r = d - 2 * np.dot(d, n) * n

    r = r / np.linalg.norm(r)

    return ip, r, meta
    
def main(args):
  txs = np.load(args.tx)
  data = WLMEData.load_by_path(args.wlma)
  wave_velocity = args.wv
  digital_sps = args.sps
  sample_time = args.time

  rx_streams = []
  for rx in data.rx:
    rx_streams.append(np.zeros((int(sample_time * digital_sps),), np.float128))

  # 1. randomly select a triangle from a random mesh
  # 2. randomly pick point on surface of triangle
  # 3. try to intersect tx; if no LOS then goto 1
  # 4. try to intersect each rx
  # 5. if intersect record TX at appropriate position and magnitude
  for cycles in range(0, args.cycles):
    print(cycles / args.cycles)
    # [1]
    mesh = random.choice(data.cage.meshes)
    face = random.choice(mesh.faces)
    # [2]
    a = face[0]
    b = face[1]
    c = face[2]
    ab = b - a
    # Somewhere between a and b.
    d = ab * random.random() + a
    dc = c - d
    # Somewhere between d and e.
    e = dc * random.random() + d
    # [3]
    tx = np.array(data.tx[0].location)
    tx_to_e = e - tx
    tx_to_e = tx_to_e / np.linalg.norm(tx_to_e)
    tx_ray = Ray(tx, tx_to_e)
    hit = data.cage.ray_single_bounce(tx_ray)
    if hit is None:
      # [3]
      continue
    if np.linalg.norm(e - hit[0]) > 0.001:
      continue
    tx_dist = np.linalg.norm(hit[0] - tx_ray.u)
    # [4]
    for rx_ndx in range(0, len(data.rx)):
      rx = data.rx[rx_ndx]
      rxloc = rx.location
      rxs = rx_streams[rx_ndx]
      rx_to_e = e - rxloc
      rx_to_e = rx_to_e / np.linalg.norm(rx_to_e)
      rx_ray = Ray(rxloc, rx_to_e)
      hit = data.cage.ray_single_bounce(rx_ray)
      if hit is None:
        # [4]
        continue
      if np.linalg.norm(e - hit[0]) > 0.001:
        continue
      rx_dist = np.linalg.norm(hit[0] - rx_ray.u)
      dist = rx_dist + tx_dist
      secs = dist / wave_velocity
      samp_index = int(secs * digital_sps)
      sz = min(txs.shape[0], rxs.shape[0] - samp_index)
      # TODO: reflector energy loss needed
      loss_dist = 1.0 / (dist ** 2.0 * np.pi * 4.0)   
      if sz > 0:
        # NOTE: The negative is the phase inversion from a reflection.
        rxs[samp_index:samp_index+sz] += -txs[0:sz] * loss_dist 

  for rx_ndx in range(0, len(rx_streams)):
    print('saving', rx_ndx)
    rxs = rx_streams[rx_ndx]
    rx_name = data.rx[rx_ndx].name
    np.save('%s.npy' % rx_name, rxs, False)

if __name__ == '__main__':
  ap = argparse.ArgumentParser()
  ap.add_argument('--tx', type=str, required=True, help='The NPY (numpy) format transmit file.')
  ap.add_argument('--wlma', type=str, required=True, help='The WLMA (blender export) data.')
  ap.add_argument('--wv', type=float, required=True, help='The wave velocity.')
  ap.add_argument('--sps', type=int, required=True, help='The sampling rate.')
  ap.add_argument('--time', type=float, required=True, help='The sampling time in seconds.')
  ap.add_argument('--cycles', type=int, required=True, help='Number of photons.')
  main(ap.parse_args())
