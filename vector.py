import math
import numpy

class Vector3:
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z

  def __repr__(self):
    return '<%s %s %s>' % (self.x, self.y, self.z)

  def cross(a, b):
    #i = a.y * b.z - a.z * b.y
    #j = a.x * b.z - a.z * b.x
    #k = a.x * b.y - a.y * b.x
    #return Vector3(i, j, k)
    _a = [a.x, a.y, a.z]
    _b = [b.x, b.y, b.z]
    c = numpy.cross(_a, _b)
    return Vector3(c[0], c[1], c[2])  


  def __add__(a, b):
    return Vector3(a.x + b.x, a.y + b.y, a.z + b.z)

  def __sub__(a, b):
    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z)

  def __mul__(a, b):
    if type(b) is Vector3:
      #return Vector3(a.x * b.x, a.y * b.y, a.z * b.z)
      return a.cross(b)
    else:
      return Vector3(a.x * b, a.y * b, a.z * b)

  def __truediv__(a, b):
    if type(b) is Vector3:
      return Vector3(a.x / b.x, a.y / b.y, a.z / b.z)
    else:
      return Vector3(a.x / b, a.y / b, a.z / b)

  def __abs__(a):
    return math.sqrt(a.x ** 2 + a.y ** 2 + a.z ** 2)
