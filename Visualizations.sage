import matplotlib.pyplot as plt
import pdb
load("Bounds.sage")

def plot_singleton_bound(n=1000):
  X = [k for k in xrange(1, n+1)]
  Y = [n-k+1 for k in X]
  return (("k", X), ("max d", Y))

def plot_plotkin_bound(n=1000):
  X = [k for k in xrange(1, n+1)]
  Y = [plotkin_bound(n=n, k=k, pp=False) for k in X]
  return (("k", X), ("max d", Y))

def plot_griesmer_bound(k=500, max_d=1000):
  X = [d for d in xrange(1, max_d+1)]
  Y = [griesmer_bound(d=d, k=k, pp=False) for d in X]
  return (("d", X), ("min n", Y))

def plot_hamming_bound(q=2, k=500, max_d=50):
  X = [d for d in xrange(1, max_d+1)]
  f = lambda d: lambda n: int(q^(n-k)) - sum([int(binomial(n, j)*(q-1)^j) for j in xrange(int((d-1)/2)+1)])
  Y = [bisect(f(d), 0, k, (k+d)*2) for d in X]
  return(("d", X), ("min n", Y))

def plot_gilbert_varshamov_bound(q=2, k=500, max_d=50):
  X = [d for d in xrange(1, max_d+1)]
  f = lambda d: lambda n: int(q^(n-k)) - int(sum([binomial(n-1, j)*(q-1)^j for j in xrange(d-1)]))
  Y = [bisect(f(d), 0, k, (k+d)*2) for d in X]
  return(("d", X), ("min n", Y))


def bisect(fun, x, lo, hi):
  if lo > hi:
    return -1
  mid = int((lo+hi)/2)
  if mid == lo:
    return mid
  if fun(mid) > x:
    return bisect(fun, x, lo, mid)
  if fun(mid) < x:
    return bisect(fun, x, mid, hi)
