'''
  https://math.stackexchange.com/questions/2934157/how-to-calculate-the-intersection-point-between-a-straight-and-a-triangle-in-py/2934219#2934219

  Modifed the accepted answer to detect and not backwards intersect. Fixed division by zero corner case.
'''
import numpy as np
import math
from math import copysign

def dif(u,v):
  s = []
  for i in range(len(u)):
    s.append(u[i]-v[i])
  return s

def mag(u):
  s = 0
  for i in range(len(u)):
    s = s + u[i] * u[i]
  return math.sqrt(s)

def sub(u,v):
  s = []
  for i in range(len(u)):
    s.append(u[i]-v[i])
  return s

def add(u,v):
  s = []
  for i in range(len(u)):
    s.append(u[i]+v[i])
  return s

def divs(u, v):
  s = []
  for i in range(len(u)):
    s.append(u[i] / v)
  return s

def prod(u,c):
  v = []
  for i in range(len(u)):
    v.append(u[i]*c)
  return v

def isin(p, a, b, c):
  def angle(p, a, b):
    a = np.subtract(a, p)
    b = np.subtract(b, p)
    ad = np.linalg.norm(a)
    bd = np.linalg.norm(b)
    if ad == 0 or bd == 0:
      return 0
    a = a / ad
    b = b / bd
    try:
      return math.acos(np.dot(a, b))
    except ValueError as err:
      print('math.acos error: %s' % err)
      return 0
  # We should have a full 360-degrees with no less
  # and no more. The only difference being error
  # from finite precision and PI.
  aa = angle(p, a, b)
  ab = angle(p, b, c)
  ac = angle(p, c, a)

  #print('angle check', aa, ab, ac, aa + ab + ac, math.pi * 2 - (aa + ab + ac))

  return (math.pi * 2 - (aa + ab + ac)) < 0.1

def eq(u, v):
  if len(u) != len(v):
    return False
  for i in range(len(u)):
    if u[i] != v[i]:
      return False
  return True

def debug_vector(path, a, b):
  with open(path, 'w') as fd:
    fd.write('%s %s %s %s %s %s\n' % (a[0], a[1], a[2], b[0], b[1], b[2]))

def debug_point(path, a):
  with open(path, 'w') as fd:
    fd.write('%s %s %s\n' % (a[0], a[1], a[2]))

def debug_dot(path, a):
  with open(path, 'w') as fd:
    fd.write('%s %s %s\n' % (a[0], a[1], a[2]))

def debug_triangle(path, a, b, c):
  with open(path, 'w') as fd:
    fd.write('%s %s %s\n' % (a[0], a[1], a[2]))
    fd.write('%s %s %s\n' % (b[0], b[1], b[2]))
    fd.write('%s %s %s\n' % (c[0], c[1], c[2]))
    fd.write('%s %s %s\n' % (a[0], a[1], a[2]))

def line_intersect_tri(q0, w, p1, p2, p3, dbg=False):
  # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersectioni

  # The equation for a point on a plane. Where:
  #   - n is the plane normal
  #   - p is a point to be tested
  #   - p_0 is a known valid point on the plane
  #   - (p - p_0) * n = 0

  # The equation for a line is:
  #   - p = l_0 + l * d
  #   - d is a member of the set of real numbers
  #   - l_0 is any point on the line
  #   - l is the vector representing the direction of the line
  #   - d is the scalar representing the distance along the line from point l_0
  
  # Using basic algebra replace the point to be tested with the equation of a
  # line. Then, reorder the equation so that we solve for `d` the distance along
  # the line respective of the line vector.

  # Using p1 (one point from the triangle) as a point on the plane with n being
  # the triangle/plane normal. The `q0` is a point on the line, `w` is the line
  # vector.
 
  n = np.cross(p2 - p1, p3 - p1)
  n = n / np.linalg.norm(n)  

  denom = np.dot(w, n)
  if np.linalg.norm(denom) == 0:
    # The line and plane are parallel.
    #print('paralle')
    return None

  numer = np.dot((p1 - q0), n)

  if np.linalg.norm(numer) == 0:
    # The line is contained completely in the plane.
    #print('line completed in plane')
    return None

  d = numer / denom

  w = w / np.linalg.norm(w)

  p = q0 + w * d

  # See if the intersection point is within the triangle.
  if isin(p, p1, p2, p3):
    wa = p - q0
    wa = wa / np.linalg.norm(wa)

    if 1.0 - np.dot(w, wa) < 0.001:
      return p

    return None
  else:
    return None

def line_intersect_plane(q0, w, p1, n):
  # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersectioni

  # The equation for a point on a plane. Where:
  #   - n is the plane normal
  #   - p is a point to be tested
  #   - p_0 is a known valid point on the plane
  #   - (p - p_0) * n = 0

  # The equation for a line is:
  #   - p = l_0 + l * d
  #   - d is a member of the set of real numbers
  #   - l_0 is any point on the line
  #   - l is the vector representing the direction of the line
  #   - d is the scalar representing the distance along the line from point l_0
  
  # Using basic algebra replace the point to be tested with the equation of a
  # line. Then, reorder the equation so that we solve for `d` the distance along
  # the line respective of the line vector.

  # Using p1 (one point from the triangle) as a point on the plane with n being
  # the triangle/plane normal. The `q0` is a point on the line, `w` is the line
  # vector.
 
  denom = np.dot(w, n)
  if np.linalg.norm(denom) == 0:
    # The line and plane are parallel.
    raise Exception('Line is parallel to the plane.')

  numer = np.dot((p1 - q0), n)

  if np.linalg.norm(numer) == 0:
    # The line is contained completely in the plane.
    #print('line completed in plane')
    return q0, 1

  d = numer / denom

  w = w / np.linalg.norm(w)

  p = q0 + w * d

  # Don't backfire into the plane.
  wa = p - q0
  wa = wa / np.linalg.norm(wa)
  if 1.0 - np.dot(w, wa) < 0.001:
    return p, 1
  return p, -1

if __name__ == '__main__':
  p1 = np.array([0,0,0])
  p2 = np.array([10,0,0])
  p3 = np.array([0,10,0])
  q0 = np.array([5,5,-5])
  w  = np.array([0,0,1])

  # [5, 5, 0]

  print(line_intersect_tri(q0, w, p1, p2, p3))

