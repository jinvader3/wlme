'''
  This is a collection of meshes for which intersection tests are made
  using rays. This module also handles loading the mesh into the 
  collection.
'''
import numpy as np
import tritopoint
import json
from tritopoint import line_intersect_tri

class Ray:
  def __init__(self, origin, vector):
    self.u = origin
    self.v = vector

class RX:
  def __init__(self, location, sensitivity=1.0):
    self.location = location
    self.sensitivity = sensitivity

class TX:
  def __init__(self, location, power):
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
        np.array(rx['location']), 
        rx['sensitivity']
      ))

    for tx in data['tx']:
      self.tx.append(TX(
        np.array(tx['location']),
        tx['power']
      ))

    for mesh in data['meshes']:
      self.cage.add_mesh(Mesh(mesh['polygons']))

    return self

class Mesh:
  def __init__(self, faces):
    _faces = []
    for face in faces:
      _face = []
      for ndx in range(3):
        _face.append(
          np.array([face[ndx][0], face[ndx][1], face[ndx][2]])
        )
      _faces.append(_face)
    self.faces = _faces

  def ray_test(self, ray, dbg=False):
    for face in self.faces:
      assert(len(face) == 3)
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
        points.append((ip, dist, face))

    # No intersections.
    if len(points) == 0: return None  

    # Return the first ray to mesh intersection.
    points.sort(key=lambda x: x[1])
    
    return points[0][0], points[0][2]

  def ray_bounces(self, ray, max_bounce_count):
    for i in range(max_bounce_count):
      yield ray.u, ray.v
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
    ip, face = res

    a = face[0] 
    b = face[1]
    c = face[2]

    print('$')
    print('a', a)
    print('b', b)
    print('c', c)
    print('ip', ip)
    print('ray.u', ray.u)

    # Get vector representing strike to surface of face.
    d = ip - ray.u
    d = d / np.linalg.norm(d)

    print('d', d)

    # Get face normal.
    n = np.cross(b - a, c - a)
    # Normalize face normal to unit length one.
    n = n / np.linalg.norm(n)
    #n = Vector3(-n.x, -n.y, -n.z)

    # Compute reflection vector.
    r = d - 2 * np.dot(d, n) * n

    r = r / np.linalg.norm(r)

    print('r', r)

    return ip, r
    

if __name__ == '__main__':
  data = WLMEData.load_by_path('test.wlmadata')
 
  ray = Ray(np.array(data.tx[0].location), np.array([0, 0, 1]))

  for point in data.cage.ray_bounces(ray, 100):
    print(point)

  '''
  with open('scene.txt', 'w') as fd:
    for mesh in data.meshes.meshes:
      for face in mesh.faces:
        fd.write('\n')
        fd.write('%s %s %s\n' % (face[0][0], face[0][1], face[0][2]))
        fd.write('%s %s %s\n' % (face[1][0], face[1][1], face[1][2]))
        fd.write('%s %s %s\n' % (face[2][0], face[2][1], face[2][2]))
        fd.write('%s %s %s\n' % (face[0][0], face[0][1], face[0][2]))
  '''


