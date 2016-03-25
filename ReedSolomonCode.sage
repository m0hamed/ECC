load("BasicLinearCode.sage")
import pdb

class ReedSolomonCode(BasicLinearCode):
  def __init__(self, n, k, F, alphas=None):
    self._base_ring = F
    self._rank = k
    self._length = n

    if n<k:
      raise Exception, "Length cannot be smaller than dimension"
    if n>len(F.list()):
      raise Exception, "Length cannot be larger than size of field"

    if alphas is not None:
      if len(alphas) > len(set(alphas)):
        raise Exception, "Alphas must be distinct"
      if len(alphas) != self._length:
        raise Exception, "Number of alphas must be equal to the length"
      self.alphas = alphas
    else:
      self.alphas = self._base_ring.list()[:n]

    # The gen matrix is the vandermonde matrix
    self._generator_matrix = matrix(F, V(self.alphas, k))

  def _get_deltas(self):
    return [
        reduce(
          lambda x,y: x*y,
          [self._base_ring(alphai - alphaj)^-1 for alphaj in self.alphas
            if alphai != alphaj]
        ) for alphai in self.alphas]

  def parity_check_matrix(self):
    return matrix(self._base_ring, V(self.alphas, self._length - self._rank)) \
        *diag(self._get_deltas())

  # Berlekamp-Welsh decoding
  # No failure checking!!
  def bw_decode(self, received_word):
    t = floor((self._length-self._rank)/2.0)
    l0 = self._length - 1 - t
    l1 = self._length - 1 - t - (self._rank - 1)
    m = [[alpha^i for i in range(l0+1)] + [r*alpha^i for i in range(l1+1)]
        for alpha, r in zip(self.alphas, received_word)]
    m = matrix(self._base_ring, m)
    # get the last row of the right kernel since this has the lowest weight
    Qs = m.right_kernel_matrix()[-1]
    # split them into Q0 and Q1
    Q0 = Qs[:l0+1]
    Q1 = Qs[l0+1:]
    # convert lists into polinomials for division
    R.<x> = self._base_ring['x']
    Q0x = R(list(Q0))
    Q1x = R(list(Q1))
    # get the quotient from the division, (should do some failure checking here)
    g, _ = (-Q0x).quo_rem(Q1x)
    # return the coefficients padded with 0s when needed
    return g.padded_list(self._rank)


# returns the diagonalization of a vector
def diag(vector):
  return matrix(len(vector), lambda i,j: vector[i] if i==j else 0)

# returns the vandermonde matrix coresponding to the alphas and rank
def V(alphas, rank):
  return [[x^i for x in alphas] for i in range(rank)]

