from math import acos, degrees, trunc

#------------------------------------------#
class Vec:
  def __init__(self, x,y,z):
    self.x= x
    self.y= y
    self.z= z
  
  def ggT(self,g,k):
    g= abs(g)
    k= abs(k)
    while k>0:
      g,k = k,g%k
    return g
  
  def __add__(self, other):
    return Vec(self.x+other.x,self.y+other.y,self.z+other.z)

  def __sub__(self, other):
    return Vec(self.x-other.x,self.y-other.y,self.z-other.z)
  
  def __repr__(self):
    return f"{'{'}{self.x}\\{self.y}\\{self.z}{'}'}"
  
  def __bool__(self):
    return abs(self.x)>1e-6 or abs(self.y)>1e-6 or abs(self.z)>1e-6
  
  def __trunc__(self):
    m= abs(self)
    return Vec(self.x/m, self.y/m, self.z/m)
  
  def __abs__(self):
    return ((self.x)**2+(self.y)**2+(self.z)**2)**0.5
  
  def __eq__(self, other):
    return self.x==other.x & self.y==other.y & self.z==other.z
  
  def __mul__(self, t):
    if isinstance(t,int) or isinstance(t,float):
      return Vec(self.x*t, self.y*t, self.z*t)
    else:
      return self.x*t.x + self.y*t.y + self.z*t.z
  
  def __rmul__(self, t):
    if isinstance(t,int) or isinstance(t,float):
      return Vec(self.x*t, self.y*t, self.z*t)
    else:
      return self.x*t.x + self.y*t.y + self.z*t.z
  
  def __matmul__(self, other):
    r1= self.y*other.z - self.z*other.y
    r2= self.z*other.x - self.x*other.z
    r3= self.x*other.y - self.y*other.x
    return Vec(r1,r2,r3)
    
  def __contains__(self, other):
    return not self@other
  
  def __invert__(self):
    r1= self.x
    r2= self.y
    r3= self.z
    if isinstance(r1,int) and isinstance(r2,int) and isinstance(r3,int):
      t= self.ggT(self.ggT(r1,r2),r3)
      return Vec(r1//t, r2//t, r3//t)
    else:
      return self
  
  def __or__(self, other):
    return Lin(self, other-self)
  
  def __xor__(self, other):
    cosphi= self*other/abs(self)/abs(other)
    return degrees(acos(cosphi))

#------------------------------------------#
class Lin:
  def __init__(self, a,v):
    self.a= a
    self.v= v
  
  def __repr__(self):
    return f"{self.a} +t{self.v}"
  
  def __add__(self, t):     # Einsetzen von t>=0 in Param'gl
    return self.a+t*self.v
  
  def __sub__(self, t):     # Einsetzen von t<0  in Param'gl
    return self.a-t*self.v
    
  def __contains__(self, point):
    return (point-self.a) in self.v
    
  def __eq__(self, other):
    return (self.v in other.v) and (self.a in other)
      
  def __and__(self, other):
    if self.v in other.v:
      if self.a in other:
        raise RuntimeError('identische Geraden!')
      else:
        raise RuntimeError('parallele Geraden!')
    else:
      trans= other.a - self.a # verbindet Referenzpunkte
      norml= self.v @ other.v # steht senkrecht zu beiden Ger
      if abs(trans*norml) > 1e-6:
        raise RuntimeError('windschiefe Geraden!')
      else:
        a11= self.v.x; a12= -other.v.x; a13= other.a.x-self.a.x
        a21= self.v.y; a22= -other.v.y; a23= other.a.y-self.a.y
        a31= self.v.z; a32= -other.v.z; a33= other.a.z-self.a.z
        detlinks= a11*a22-a21*a12
        if detlinks != 0:
          dets= a12*a23-a22*a13 # parameter fuer self
          dett= a11*a23-a21*a13 # parameter fuer other
          s= -dets/detlinks
          t= dett/detlinks
        else:
          detlinks= a21*a32-a31*a22
          dets= a22*a33-a32*a23
          dett= a21*a33-a31*a23
          s= -dets/detlinks
          t= dett/detlinks
        return self+s

#------------------------------------------#
class Pla:
  def __init__(self, a,b,c,d):
    self.a= a
    self.b= b
    self.c= c
    self.d= d
    self.n= Vec(a,b,c)
  
  def pl_str(self, x):
    if x>=0:
      return '+'+str(x)
    else:
      return str(x)
      
  def __repr__(self):
    return f"{str(self.a)}x{self.pl_str(self.b)}y{self.pl_str(self.c)}z={str(self.d)}"
  
  def __contains__(self, other):
    if isinstance(other, Vec):
      return self.n*other == self.d
    if isinstance(other, Lin):
      return (self.n*other.v == 0) and (other.a in self)

#------------------------------------------#  
a= Vec(1,2,3)
E= Pla(2,2,1,9)
print(E)
print(a in E)