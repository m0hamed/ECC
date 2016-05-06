load("BasicLinearCode.sage")

class HammingCode(BasicLinearCode):
  #Todo: commenting
  def __init__(self, m):
    length = 2^m-1
    rank = length-m
    Ir = matrix.identity(rank)
    gen_rows = []

    #untested, unbased, completely derivative determination of table (-A^T)
    # its tested now, but the rest is still true :)
    for i in range(0,m-1):
      base = ([0]*(i))+([1]*(length-rank-i))
      while True:
        gen_rows.insert(0,copy(base))
        if not next_permutation(base):
          break
    A = matrix(GF(2),(gen_rows));

    self._generator_matrix = matrix(GF(2),matrix.identity(rank).augment(A))
    self._parity_check_matrix = A.transpose().augment(matrix.identity(length-rank))
    self._base_ring = GF(2)
    self._rank, self._length = self._generator_matrix.dimensions()

  def parity_check_matrix(self):
    return self._parity_check_matrix

  def minimum_distance(self):
    return 3

#test of hamming generation
#m = 4
#h = hamming(m)
#print "m=", m, ", Hamming(", h[0], ",", h[1], ")"
#print "Generator:\n", h[2]
#print "Parity:\n", h[3]
#print "H*G^T:\n", h[3]*h[2].transpose()
