import numpy as np

def Delta(Theta1Rad, Phi1Rad, Theta2Rad, Phi2Rad):
  r''' Computes angular difference between two points in spherical coordinates.
       All angles are in radians with Theta defined from North pole of sphere.

  Parameters
  ----------
  Theta1Rad : array
    list of theta values for location 1 [radians]
  Phi1Rad : array
    list of phi values for location 1   [radians]
  Theta2Rad : array
    list of theta values for location 2 [radians]
  Phi2Rad : array
    list of phi values for location 2   [radians]

  Returns
  -------
  delta : array
    array with space angle difference between locations 1 & 2 [radians]
  '''

  DeltaPhiRad = Phi1Rad - Phi2Rad

  Vector = np.cos(Theta1Rad)*np.cos(Theta2Rad) + \
                np.sin(Theta1Rad)*np.sin(Theta2Rad)*np.cos(DeltaPhiRad)

  delta = np.zeros( len(Vector) )
  delta[ Vector >=  1.0 ] = 0.0
  delta[ Vector <= -1.0 ] = np.pi

  mask = (-1.0 < Vector) & (Vector < 1.0)
  delta[ mask ] = np.arccos( Vector[mask] )

  return delta


def stats(x, verbose = True):

  if verbose:
    print(" stats for %d entries" % x.size)

  mean = x.sum()/x.size
  median = np.percentile(x, 50)

  return mean, median

def weighted_stats(x, xmin, xmax, w, nbin = 5000):
  r''' Get mean and median using event weights

  Parameters
  ----------
  x : array
    array of values
  xmin : float
    minimum value of x
  xmax : float
    maximum value of x
  w : array
    array of weights

  Returns
  -------
  weighted mean, median, and p-value of median
  '''

  # make sure weight is normalized to 1
  w /= w.sum()  

  # create weighted histogram
  h, bins = np.histogram(x, np.linspace(xmin,xmax,nbin+1), weights = w)

  # search for median
  p = 0.
  for i in range(h.size):
    p += h[i]
    median = bins[i+1]
    if p > 0.5: break

  # calculate weighed mean
  mean = (x*w).sum() / w.sum()

  return mean, median, p
