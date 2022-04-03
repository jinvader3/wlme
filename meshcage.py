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

# http://resource.npl.co.uk/acoustics/techguides/seaabsorption/
FRESH_WATER_ATT_1HZ_TO_48KHZ_PER_KM = [
  (1, 0.006),
  (10, 0.042),
  (100, 3.462),
  (250, 21.195),
  (500, 84.41),
  (1000, 337.246),
  (2000, 1348.579),
  (3000, 3034.135),
  (5000, 8427.912),
  (10000, 33711.244),
  (20000, 134844.57),
  (48000, 776704.083),
  (96000, 3100990.108),
  (192000, 12403956.538),
]

def atten_fresh_water_8c_shallow_8ph_approx(freq, meters):
  tbl = FRESH_WATER_ATT_1HZ_TO_48KHZ_PER_KM
  for x in range(0, len(tbl)):
    if tbl[x][0] > freq:
      aa = tbl[x-1][1]
      ab = tbl[x][1]
      fa = tbl[x-1][0]
      fb = tbl[x][0]
      db = ((freq - fa) / (fb - fa)) * (ab - aa) + aa
      return meters / 1000 * db
  raise Exception('The frequency is out of upper range.')

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
    

if __name__ == '__main__':
  data = WLMEData.load_by_path('100meterpond.wlmadata')

  #fd = open('points.txt', 'w') 

  # Velocity of sound in water in meters/second at a normal temperature.
  # I'd love to be able to use blender to create a volume of water with
  # a mixture of temperatures just like the real thing because temp will
  # alter the wave velocity and create reflections!!!!
  wave_velocity = 1480.0
  #wave_velocity = 343.0
  # Lets play like we are using sound card sample rates!!! Cheap!!!
  digital_sps = 192000
  # Create a chirp that will be TX'ed into the water.
  chirp = get_chirp(digital_sps, 0.1, 64000, 96000)

  np.save('chirp.npy', chirp, False)

  # The total number of seconds to sample for.
  sample_time = 1
  # Convert the chirp into the frequency domain.
  chirp_fft = fft.rfft(chirp)

  rx_streams = []
  for rx in data.rx:
    rx_streams.append(np.zeros((int(sample_time * digital_sps),), np.float128))

  rc_limit = 100

  for rc in range(0, rc_limit):
    # Create a ray from the TX source with a random trajectory. Let it
    # bounce around and each contact with a surface see if we can ray
    # trace to each RX sensor. If we can ray trace to an RX isotrophic
    # sensor then calculate final attenuate and delay and write that
    # into the RX buffer.
    rv = np.random.random(3) * 2.0 - 1.0
    rv = rv / np.linalg.norm(rv)
    ray = Ray(np.array(data.tx[0].location), rv)

    traveldist = 0.0
    energy = 1.0

    for p in data.cage.ray_bounces(ray, 10):
      # How far was this first run/jump?
      jumpdist = np.linalg.norm(p[0] - ray.u)
      traveldist += jumpdist
      # Verify version 1 meta-data.
      assert(p[2][0])
      # How much energy did we lose?
      energy = energy * p[2][1]
      # How smooth or how much scattering will this bounce have?
      scatter = p[2][2]

      # For each RX sensor.
      for rx_ndx in range(0, len(data.rx)):
        rx = data.rx[rx_ndx]
        rxloc = rx.location

        # See if we can ray trace to the rxloc without intersecting any
        # meshes. To do this first grab a vector from our point to the
        # rxloc.
        rxvec = rxloc - p[0]
        rxdist = np.linalg.norm(rxvec)
        rxvec = rxvec / rxdist
        # Bump the ray forward a bit out of the mesh it hit so it does
        # not re-intersect without going forward.
        # TODO: Corner case.. Need to ignore just this one specific triangular
        #       face it intersected with to begin with.
        rxray = Ray(p[0] + rxvec * 0.00001, rxvec)
        hit = data.cage.ray_single_bounce(rxray)
        # If the cage is closed we WILL hit something. Just make sure what
        # we hit was further than the rxloc.
        if hit is None:
          hitdist = rxdist + 1.0
        else:
          hitdist = np.linalg.norm(hit[0] - ray.u)

        if hitdist > rxdist:
          # Okay, now it did reach the rxloc without intersection. Calculate
          # the total attenuation and delay.
          totaldist = traveldist + rxdist
          # TODO: I think I got this right for a sphere expanding.
          distloss = 1 / (totaldist ** 2.0 * math.pi * 4)

          binfreqs = fft.rfftfreq(len(chirp), 1.0 / digital_sps)
          #'''
          tmp_fft = np.zeros((len(chirp_fft),), chirp_fft.dtype)
          for x in range(0, len(tmp_fft)):
            fldb = atten_fresh_water_8c_shallow_8ph_approx(binfreqs[x], totaldist)
            fl = 1 / (10 ** (fldb / 20))
            tmp_fft[x] = chirp_fft[x] * fl
          tmp = fft.irfft(tmp_fft) 
          #'''
          #tmp = chirp
          #print('totaldist=%s energy=%s distloss=%s fl=%s' % (totaldist, energy, distloss, fl))

          # Convert the delay distance into seconds but add the final distance.
          t = totaldist / wave_velocity
          # Convert the delay into a digital sample delay.
          # TODO: Do fractional sample shifting.
          sample_index = int(t * digital_sps)
          #print('@@@', rx_ndx, sample_index, t, totaldist, rxdist)
          # TODO: Give each sensor a different startup and add noise to its sample rate.
          # This is a big element I want to investigate. And, different types of noise
          # could show me just how bad things were or were not during real testing.

          # Calculate maximum size. Don't overwrite buffer or underwrite it.
          sz = min(chirp.shape[0], rx_streams[0].shape[0] - sample_index)
          if sz > 0:
            # The heart of the simulator is here.
            rx_streams[rx_ndx][sample_index:sample_index+sz] = \
              tmp[0:sz] * energy * distloss
          print('.', rc / rc_limit)
        else:
          # It hit something before it got there.
          pass

  for rx_ndx in range(0, len(rx_streams)):
    print('saving', rx_ndx)
    rxs = rx_streams[rx_ndx]
    np.save('rx%s.npy' % rx_ndx, rxs, False)

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
