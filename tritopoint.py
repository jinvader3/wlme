'''
  https://math.stackexchange.com/questions/2934157/how-to-calculate-the-intersection-point-between-a-straight-and-a-triangle-in-py/2934219#2934219

  Modifed the accepted answer to detect and not backwards intersect. Fixed division by zero corner case.
'''
import numpy as np
import math

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

def isin(p,p1,p2,p3):
  u = dif(p1,p2)
  v = dif(p3,p2)
  w = dif(p,p2)
  l1 = np.dot(u,w)
  if l1 < 0: 
    return False
  l2 = np.dot(v,w)
  if l2 < 0: 
    return False
  u = dif(p1,p3)
  v = dif(p2,p3)
  w = dif(p,p3)
  l1 = np.dot(u,w)
  if l1 < 0: 
    return False
  l2 = np.dot(v,w)
  if l2 < 0: 
    return False
  return True

def eq(u, v):
  if len(u) != len(v):
    return False
  for i in range(len(u)):
    if u[i] != v[i]:
      return False
  return True

def line_intersect_tri(q0, w, p1, p2, p3):
  u  = dif(p1,p2)
  v  = dif(p3,p2)
  n  = list(np.cross(u,v))
  p0 = prod(add(add(p1,p2),p3),1/3.)
  a = np.dot(w, n)
  if a == 0:
    print('external point via a')
    return None

  #pi = add(q0,prod(w,-np.dot(dif(q0,p0),n)/np.dot(w,n)))
  pi = add(q0,prod(w, -np.dot(dif(q0,p0),n) / a))

  if isin(pi,p1,p2,p3):
    w = divs(w, mag(w))
    wa = sub(pi, q0)
    wa = divs(wa, mag(wa))

    print('w=%s wa=%s' % (w, wa))    

    if eq(w, wa):
      print(pi, " internal point")
      return pi

    print(pi, 'backwards intersection')
    return None
  else:
    print(pi, " external point")
    return None

if __name__ == '__main__':
  p1 = [0,0,0]
  p2 = [10,0,0]
  p3 = [0,10,0]
  q0 = [5,5,-5]
  w  = [0,0,1]

  line_intersect_tri(q0, w, p1, p2, p3)

