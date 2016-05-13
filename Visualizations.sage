import matplotlib.pyplot as plt
from itertools import cycle

load("Bounds.sage")
load("DecodeTester.sage")
load("QarySymmetricChannel.sage")

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

def transmit_over_codes(tests,crossovers,codes):
  for code in codes:
    code["bit_error_rates"] = []
    code["block_error_rates"] = []
    code["decoding_time"] = 0
    code["interpolation"] = 0
    code["root_finding"] = 0

  for crossover in crossovers:
    for code in codes:

      results = QarySymmetricChannel(code["code"],code["decode"],code["unencode"],crossover/float(100)).transmit(tests)
      code["bit_error_rates"] += [results[0]*100]
      code["block_error_rates"] += [results[1]*100]
      if type(results[2]) == tuple:
        code["decoding_time"] += results[2][0]
        code["interpolation"] += results[2][1]
        code["root_finding"] += results[2][2]
      else:
        code["decoding_time"] += results[2]

  for code in codes:
    print "Average decoding time for %s: %f sec"%(code["name"],code["decoding_time"]/float(len(crossovers)))
    if code["interpolation"] > 0:
      print "Average interpolation time for %s: %f sec"%(code["name"],code["interpolation"]/float(len(crossovers)))
    if code["root_finding"] > 0:
      print "Average root finding time for %s: %f sec"%(code["name"],code["root_finding"]/float(len(crossovers)))

def plot_codes(tests,crossovers,codes,title,filename,y_axis="bit_error_rates"):

  transmit_over_codes(tests,crossovers,codes)

  fig = plt.figure()

  ax = plt.gca()
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)

  plt.title(title)
  for code in codes:
    plt.plot(crossovers, code[y_axis], label=code["name"])
  plt.legend(prop={'size':10},loc=2)
  if y_axis=="bit_error_rates":
    plt.ylabel("Bit error rate (%)")
  else:
    plt.ylabel("Block error rate (%)")
  plt.xlabel("Qary channel crossover probability (%)")
  fig.savefig(filename)
  plt.close(fig)

def plot_hamming_vs_none():
  hamming_m = 4
  F.<a> = GF(2)
  HC = HammingCode(hamming_m, F)
  HC_desc = {"code":HC, "name":"4th binary Hamming Code", "decode":HC.decode, "unencode":HC.unencode}
  NC = NoCode(GF(2), HC._rank)
  NC_desc = {"code":NC, "name":"No Code", "decode":NC.decode, "unencode":NC.unencode}
  crossovers = range(1,5,1)

  title = "%dth Hamming code (%d,%d) in F_%d vs No code" % (hamming_m,HC._length,HC._rank,F.order())
  filename = '%dhammingVsNone-F_%d-bit.pdf'% (hamming_m,F.order())
  plot_codes(1000,crossovers,[HC_desc,NC_desc],title,filename,y_axis="bit_error_rates")

  title = "%dth Hamming code (%d,%d) in F_%d vs No code" % (hamming_m,HC._length,HC._rank,F.order())
  filename = '%dhammingVsNone-F_%d-block.pdf'% (hamming_m,F.order())
  plot_codes(1000,crossovers,[HC_desc,NC_desc],title,filename,y_axis="block_error_rates")

def plot_reed_solomon_different_rates():
  F.<a> = GF(16)
  length = 12
  rank = 8
  RS1 = ReedSolomonCode(length, rank, F)
  RS1_desc = {"code":RS1, "name":("(%d,%d,%d) RS code in F_%d" % (RS1._length,RS1._rank,RS1._min_dist,F.order())), "decode":RS1.bw_decode, "unencode":None}

  RS2 = ReedSolomonCode(length+2, rank, F)
  RS2_desc = {"code":RS2, "name":("(%d,%d,%d) RS code in F_%d" % (RS2._length,RS2._rank,RS2._min_dist,F.order())), "decode":RS2.bw_decode, "unencode":None}

  RS3 = ReedSolomonCode(length+4, rank, F)
  RS3_desc = {"code":RS3, "name":("(%d,%d,%d) RS code in F_%d" % (RS3._length,RS3._rank,RS3._min_dist,F.order())), "decode":RS3.bw_decode, "unencode":None}

  crossovers = range(1,21,2)

  plot_codes(1000,crossovers,[RS1_desc,RS2_desc,RS3_desc],"Reed Solomon Simulation","rssim_3.pdf",y_axis="block_error_rates")

def guruswami_sudan_timing_simulation():
  F.<a> = GF(13)
  length = 11
  rank = 5
  crossovers = [10]
  RS1 = ReedSolomonCode(length, rank, F)
  codes = []
  codes.append({"code":RS1, "name":("(%d,%d,%d) RS code in F_%d(EEA)" % (RS1._length,RS1._rank,RS1._min_dist,F.order())), "decode":RS1.eea_decode, "unencode":None})
  codes.append({"code":RS1, "name":("(%d,%d,%d) RS code in F_%d(BW)" % (RS1._length,RS1._rank,RS1._min_dist,F.order())), "decode":RS1.bw_decode, "unencode":None})
  for ip in ["knh","book"]:
    for root in ["simple","book","roth"]:
      codes.append({"code":RS1, "name":("(%d,%d,%d) RS code in F_%d(GS,%s,%s)" % (RS1._length,RS1._rank,RS1._min_dist,F.order(),ip,root)), "decode":lambda x: RS1.gs_decode(x, interpolate=ip, roots=root, timing=True), "unencode":None})
  transmit_over_codes(10,crossovers,codes)

if __name__ == "__main__":
  #plot_singleton_bound()
  #plot_plotkin_bound()
  #plot_griesmer_bound()
  #plot_hamming_bound()
  #plot_gilbert_varshamov_bound()

  # F.<a> = GF(2)
  # plot_linear_code_test(F,10,3,1000)

  # plot_hamming_vs_none()


  pass

guruswami_sudan_timing_simulation()
