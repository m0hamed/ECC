import matplotlib.pyplot as plt

load("Bounds.sage")

def plot_singleton_bound(n=1000):
  X = [k for k in xrange(1, n+1)]
  Y = [n-k+1 for k in X]
  draw_graph("Singleton Bound", ("k", X), ("max d", Y), "n=%i"%n)

def plot_plotkin_bound(n=50):
  X = [k for k in xrange(1, n+1)]
  Y = [plotkin_bound(n=n, k=k, pp=False) for k in X]
  draw_graph("Plotkin Bound", ("k", X), ("max d", Y), "n=%i"%n)

def plot_griesmer_bound(k=500, max_d=1000):
  X = [d for d in xrange(1, max_d+1)]
  Y = [griesmer_bound(d=d, k=k, pp=False) for d in X]
  draw_graph("Griesmer Bound", ("d", X), ("min n", Y), "k=%i"%k)

def plot_hamming_bound(q=2, k=500, max_d=50):
  X = [d for d in xrange(1, max_d+1)]
  f = lambda d: lambda n: int(q^(n-k)) - sum([int(binomial(n, j)*(q-1)^j) for j in xrange(int((d-1)/2)+1)])
  Y = [bisect(f(d), 0, k, (k+d)*2) for d in X]
  draw_graph("Hamming Bound", ("d", X), ("min n", Y), "k=%i, q=%i"%(k, q))

def plot_gilbert_varshamov_bound(q=2, k=500, max_d=50):
  X = [d for d in xrange(1, max_d+1)]
  f = lambda d: lambda n: int(q^(n-k)) - int(sum([binomial(n-1, j)*(q-1)^j for j in xrange(d-1)]))
  Y = [bisect(f(d), 0, k, (k+d)*2) for d in X]
  draw_graph("Gilbert Varshamov Bound", ("d", X), ("min n", Y), "k=%i, q=%i"%(k, q))

def draw_graph(name, X, Y, text=None):
  (xlabel, Xs), (ylabel, Ys) = X, Y
  plt.plot(Xs, Ys, 'r')
  plt.title(name)
  if text: plt.text(max(Xs)*0.8, max(Ys)*0.8, text, fontdict={"size": 16})
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.show()

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


if __name__ == "__main__":
  plot_singleton_bound()
  plot_plotkin_bound()
  plot_griesmer_bound()
  plot_hamming_bound()
  plot_gilbert_varshamov_bound()
  pass
