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
    self._lagrange_functions = None

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

    #the positions of the errors are the roots of the lambda function, which is equal to g(x) (scaled)
    error_positions = [x[0] for x in eea_results[2].roots()]

    #calculate the lagrange interpolation
    lagrange_interpolation = 0
    k = 0
    # print self._get_lagrange_functions()
    for i in range(self._length):
      if i not in error_positions:
        # print i,received_word[i], received_word[i]*self._get_lagrange_functions()[i]
        lagrange_interpolation += received_word[i]*self._get_lagrange_functions()[i]
        k += 1
      if k == self._rank:
        break
  
    return lagrange_interpolation

  # returns the constituents of the lagrange interpolation
  def _get_lagrange_functions(self):
    if self._lagrange_functions is None:
      PF.<x> = self._base_ring[]
      self._lagrange_functions = [
        self._get_gammas()[i]*reduce(
          lambda z,y: z*y, [(x-self.alphas[j]) 
          for j in range(self._length) if j!=i]
          ) for i in range(self._length)]
    return self._lagrange_functions

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
# return value is a tuple: (ri, fi, gi)
def EEA(a,b,t):
  if a.degree()<b.degree():
    a,b = b,a
  M = matrix([[a,1,0],[b,0,1]])

  while M[1,0] != 0 and t <= M[0,0].degree():
    quo = M[0,0].quo_rem(M[1,0])[0]
    M = matrix([[0,1],[1,-quo]])*M

  return M[0]

# RS1 = ReedSolomonCode(7, 3, GF(7))
# msg = vector(GF(7),[0,1,1])
# cw = RS1.encode(msg)
# error = vector(GF(7),[1,0,0,0,0,3,0])
# received = cw + error
# print received
# print RS1.eval_encode(msg)

# print RS1.bw_decode(cw)
# print "decode, no error:\n",RS1.eea_decode(cw)
# print "decode, 1 error:\n",RS1.eea_decode(received)


# F.<x> = GF(13)[]
# a = x^3+x^2-x+4
# b = x^2+x+1
# c = 6*x^3 + 5*x^2 + 9*x + 7

# print "roots:",c.roots()
# print "polynomials:",a,",",b
# d,u,v = xgcd(a,b)
# print "EEA:",d,",",u,",",v
# print "custom EEA:\n",EEA(a,b,0)
