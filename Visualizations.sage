import matplotlib.pyplot as plt
from itertools import cycle

load("Bounds.sage")

def plot_singleton_bound(n=1000, plot=True):
  X = [(n-k+1)/n for k in xrange(1, n+1)]
  Y = [k/n for k in xrange(1, n+1)]
  if plot:
    draw_graph("Singleton Bound", ("delta=d/n", X), ("R=k/n", Y))
  else:
    return ("Singleton Bound", X, Y)

def plot_plotkin_bound(n=10, plot=True):
  Ks = range(1, n+1)
  X = [k/n for k in Ks]
  Y = [plotkin_bound(n=n, k=k, pp=False)/n for k in Ks]
  if plot:
    draw_graph("Plotkin Bound", ("delta=d/n", Y), ("R=k/n", X))
  else:
    return ("Plotkin Bound", X, Y)

def plot_griesmer_bound(k=500, max_d=1000, plot=True):
  Ds = range(1, max_d+1, 2)
  Y = [griesmer_bound(d=d, k=k, pp=False) for d in Ds]
  X = [d/n for d,n in zip(Ds, Y)]
  Y = map(lambda n: k/n, Y)
  if plot:
    draw_graph("Griesmer Bound", ("delta=d/n", X), ("R=k/n", Y))
  else:
    return ("Griesmer Bound", X, Y)

def plot_hamming_bound(q=2, k=300, max_d=300, plot=True):
  Ds = range(1, max_d+1, 2)
  f = lambda d: lambda n: int(q^(n-k)) - sum([int(binomial(n, j)*(q-1)^j) for j in xrange(int((d-1)/2)+1)])
  Ns = [bisect(f(d), 0, k, (k+d)*2) for d in Ds]
  X = [d/n for d,n in zip(Ds, Ns)]
  Y = [k/n for n in Ns]
  if plot:
    draw_graph("Hamming Bound", ("delta=d/n", X), ("R=k/n", Y))
  else:
    return ("Hamming Bound", X, Y)

def plot_gilbert_varshamov_bound(q=2, k=300, max_d=300, plot=True):
  Ds = range(1, max_d+1, 2)
  f = lambda d: lambda n: int(q^(n-k)) - int(sum([binomial(n-1, j)*(q-1)^j for j in xrange(d-1)]))
  Ns = [bisect(f(d), 0, k, (k+d)*2) for d in Ds]
  X = [d/n for d,n in zip(Ds, Ns)]
  Y = [k/n for n in Ns]
  #print zip(map(float,X),map(float,Y))
  if plot:
    draw_graph("Gilbert Varshamov Bound", ("delta=d/n", X), ("R=k/n", Y))
  else:
    return ("Gilbert Varshamov Bound", X, Y)

def draw_graphs(graph_array):
  lines = ["-","--","-.",":"]
  linecycler = cycle(lines)
  plt.xlabel("delta=d/n", fontdict={"size": 28})
  plt.ylabel("R=k/n", fontdict={"size": 28})
  plt.axis([0,0.5,0,1])
  plt.hold(True)
  for label, X, Y in graph_array:
    plt.plot(X, Y, next(linecycler), label=label)
    plt.tick_params(axis='both', labelsize=22)
  plt.legend(prop={'size':22})
  plt.show()

def draw_graph(name, X, Y, text=None):
  (xlabel, Xs), (ylabel, Ys) = X, Y
  plt.plot(Xs, Ys, 'r')
  plt.title(name)
  if text: plt.text(max(Xs)*0.8, max(Ys)*0.8, text, fontdict={"size": 16})
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.show()

def bisect(fun, x, lo, hi, eps=0.5):
  if lo > hi:
    return -1
  mid = (lo+hi)/2
  #print float(lo), float(mid), float(hi), float(mid-lo), fun(mid)
  if mid-lo < eps:
    return mid
  if fun(mid) > x:
    return bisect(fun, x, lo, mid)
  if fun(mid) < x:
    return bisect(fun, x, mid, hi)
  if fun(mid) == x:
    return mid


def plot_linear_code_test(F,n,k,tests):
  results = test_decode(F,k, n, tests)
  X = range(1,n+1)

  fig = plt.figure()

  ax = plt.gca()
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)

  plt.title("Random (%d,%d,%d) code in F_%d." % (n,k,results[0].minimum_distance(),F.order()))
  plt.plot(X, results[1], label='Nearest Neighbor')
  plt.plot(X, results[2], label='Syndrome')
  plt.legend()
  plt.ylabel("Successful Decoding (%)")
  plt.xlabel("Number of errors")
  fig.savefig('decodeTester%d-%d-%d-F_%d.pdf'% (n,k,results[0].minimum_distance(),F.order()))
  plt.close(fig)

def plot_hamming_test(F,hamming_m,tests,crossovers):
  HC = HammingCode(hamming_m, F)
  NC = NoCode(F, HC._rank)

  NC_bit_error_rates = []
  HC_bit_error_rates = []
  NC_block_error_rates = []
  HC_block_error_rates = []

  for crossover in crossovers:
    Ch1 = QarySymmetricChannel(NC,NC.decode,NC.unencode,crossover/float(100))
    NC_results = Ch1.transmit(tests)
    NC_bit_error_rates += [NC_results[0]*100]
    NC_block_error_rates += [NC_results[1]*100]
    Ch2 = QarySymmetricChannel(HC,HC.decode,HC.unencode,crossover/float(100))
    HC_results = Ch2.transmit(tests)
    HC_bit_error_rates += [HC_results[0]*100]
    HC_block_error_rates += [HC_results[1]*100]



  fig = plt.figure()

  ax = plt.gca()
  plt.axis('equal')
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)

  plt.title("%dth Hamming code (%d,%d) in F_%d vs No code" % (hamming_m,HC._length,HC._rank,F.order()))
  plt.plot(crossovers, NC_bit_error_rates, label='No Code')
  plt.plot(crossovers, HC_bit_error_rates, label='4th Hamming Code')
  plt.legend(loc=4)
  plt.ylabel("Bit error rate (%)")
  plt.xlabel("Qary channel crossover probability (%)")
  fig.savefig('%dhammingVsNone-F_%d-bit.pdf'% (hamming_m,F.order()))
  plt.close(fig)

  fig = plt.figure()

  ax = plt.gca()
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)

  plt.title("%dth Hamming code (%d,%d) in F_%d vs No code" % (hamming_m,HC._length,HC._rank,F.order()))
  plt.plot(crossovers, NC_block_error_rates, label='No Code')
  plt.plot(crossovers, HC_block_error_rates, label='4th Hamming Code')
  plt.legend(loc=4)
  plt.ylabel("Block error rate (%)")
  plt.xlabel("Qary channel crossover probability (%)")
  fig.savefig('%dhammingVsNone-F_%d-block.pdf'% (hamming_m,F.order()))
  plt.close(fig)

if __name__ == "__main__":
  #plot_singleton_bound()
  #plot_plotkin_bound()
  #plot_griesmer_bound()
  #plot_hamming_bound()
  #plot_gilbert_varshamov_bound()

  # F.<a> = GF(2)
  # plot_linear_code_test(F,10,3,1000)

  # F.<a> = GF(7)
  #crossovers = range(1,16,1)
  #plot_hamming_test(F,2000,2)

  pass