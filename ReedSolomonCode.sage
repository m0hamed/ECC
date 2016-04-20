load("BasicLinearCode.sage")
import pdb

class ReedSolomonCode(BasicLinearCode):
  def __init__(self, n, k, F, alphas=None):
    self._base_ring = F
    self._rank = k
    self._length = n
    self._min_dist = n-k+1
    self._t = floor((self._min_dist-1)/2)

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
    self._parity_check_matrix = None
    self._gammas = None

  def _get_gammas(self):
    if (self._gammas is None):
      self._gammas = [
          reduce(
            lambda x,y: x*y,
            [self._base_ring(alphai - alphaj)^-1 for alphaj in self.alphas
              if alphai != alphaj]
          ) for alphai in self.alphas]
    return self._gammas

  def parity_check_matrix(self):
    if (self._parity_check_matrix is None):
      self._parity_check_matrix = matrix(self._base_ring, V(self.alphas, self._length - self._rank)) \
          *diag(self._get_gammas())
    return self._parity_check_matrix

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

  # Peterson algorithm decoding through the Extended Euclidean Algorithm
  # works up to lagrange interpolation, fails horribly there
  def eea_decode(self, received_word):
    #calculate syndromes for the received word
    syndromes = self.get_syndromes(received_word)

    #generate the syndrome polynomial
    PF.<x> = self._base_ring[]
    syndrome_polynomial = reduce(lambda z,y: z+y, [int(syndromes[self._t*2-1-i])*x^i for i in range(self._t*2)])
    #run the extended euclidean algorithm between x^2t and the syndrome polynomial, up to the first iteration where the degree of the remainder is less than t
    eea_results = EEA(x^(2*self._t),syndrome_polynomial,self._t)

    #the alphas at the error positions are the roots of the lambda function, which is equal to g(x) (scaled)
    error_alphas = [x[0] for x in eea_results.roots()]
    correct_alphas = set(self.alphas)-set(error_alphas)
    #calculate the lagrange interpolation

    lagrange_interpolation = self._lagrange(correct_alphas,received_word)
    # print "lagrange:", lagrange_interpolation

    if lagrange_interpolation.degree() >= self._rank:
      return None
    return vector(self._base_ring,lagrange_interpolation.list()+[0]*(self._rank-len(lagrange_interpolation.list())))

  # returns the lagrange interpolation for the current received word and its correct positions 
  def _lagrange(self, correct_alphas, received_word):
    PF.<x> = self._base_ring[]

    alphamap = {a_i:i for i, a_i in enumerate(self.alphas)}
    return sum([
      received_word[alphamap[a_i]]*reduce(
        lambda z,y: z*y, [(x-a_j)/(a_i-a_j) 
        for a_j in correct_alphas if a_j!=a_i]
        ) for a_i in correct_alphas])

  # returns the syndromes of the received word
  def get_syndromes(self, received_word):
    return [
      reduce(
        lambda x,y: x+y, [self.alphas[h]^(m)*self._get_gammas()[h]*received_word[h] 
        for h in range(self._length)]
        )for m in range(self._t*2)]

  # returns the message encoded through evaluation
  def eval_encode(self, message):
    return [
      reduce(
        lambda x,y: x+y, [message[i]*alpha^i for i in range(self._rank)]
        )for alpha in self.alphas]

# returns the diagonalization of a vector
def diag(vector):
  return matrix(len(vector), lambda i,j: vector[i] if i==j else 0)

# returns the vandermonde matrix coresponding to the alphas and rank
def V(alphas, rank):
  return [[x^i for x in alphas] for i in range(rank)]

# returns the extended euclidean algorithm results ran up to the first iteration where the degree of the remainder is less than t
# return value is gi
def EEA(a,b,t):
  if a.degree()<b.degree():
    a,b = b,a
  M = matrix([[a,0],[b,1]])

  while M[1,0] != 0 and t <= M[1,0].degree():
    quo = M[0,0].quo_rem(M[1,0])[0]
    M = matrix([[0,1],[1,-quo]])*M

  return M[1][1]

# RS1 = ReedSolomonCode(7, 3, GF(13), alphas=[2,3,4,5,6,7,8])
# print "n:",RS1._length,", k:",RS1._rank,", t:",RS1._t
# msg = vector(GF(13),[5,0,1])
# cw = RS1.encode(msg)
# error = vector(GF(13),[1,0,0,0,6,0,0])
# received = cw + error
# print received
# print RS1.eval_encode(msg)

# print "bw decode, no error:\n",RS1.bw_decode(cw)
# print "bw decode, error(s):\n",RS1.bw_decode(received)
# print "eea decode, no error:\n",RS1.eea_decode(cw)
# print "eea decode, error(s):\n",RS1.eea_decode(received)
