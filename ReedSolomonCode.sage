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

  # Guruswami-sudan decoding algorithm
  def gs_decode(self, received_word, tau, s, l):
    # get polynomial coefficients
    Qlist = self._interpolate(received_word, tau, s, l)
    # convert coefficients into sage polynomial
    Qpoly = self._make_poly(Qlist, tau, s, l)
    # TODO: Find roots
    return self._find_roots(Qpoly)

  # takes a list of coefficients and returns a sage polynomial
  def _make_poly(self, Qlist, tau, s, l):
    # create polynomial ring
    R.<x,y> = self._base_ring[]
    # final polynomial Q(x,y)
    Qpoly = 0
    # a goes up to and including l
    for a in xrange(l+1):
      # polynomial Qa
      Qa = 0
      # b goes up to and including l_a
      for b in xrange(s*(self._length - tau) - a*(self._rank - 1) - 1 + 1):
        Qa += Qlist.pop(0)*(x^b)
      Qpoly += Qa*(y^a)
    return Qpoly

  # returns a long list of coefficients of Q(x,y) as a list
  def _interpolate(self, received_word, tau, s, l):
    # matrix for the linear equations
    m = []
    for alpha, rec in zip(self.alphas, received_word):
      # row of linear equations for each alpha and received word element
      row = []
      # a goes up to and including l
      for a in xrange(l+1):
        # b goes up to and including l_a
        for b in xrange(s*(self._length - tau) - a*(self._rank - 1) - 1 + 1):
          # the coefficent of Q_ab
          coeff = 0
          # h goes up to and including a but excluding s
          for h in xrange(min(s, a + 1)):
            # r goes up to and including b but excluding s-h
            for r in xrange(min(s-h, b + 1)):
              coeff += binomial(a,h)*binomial(b,r)*rec^(a-h)*alpha^(b-r)
          row.append(coeff)
      m.append(row)
    # make a sage matrix out of it
    m = matrix(self._base_ring, m)
    # get one of the rows of the right kernel matrix as a list
    return list(m.right_kernel_matrix()[-1])

  #factorizes the Q polynomial, finds factors of type "y-f(x)", and returns the list of f(x)
  def _find_roots(self, Qpoly):
    R.<x,y> = self._base_ring[]
    results = []
    #factorize the multivariate polynomial, and iterate over the factors
    #factor[0] holds the factor, factor[1] the number of occurences of said factor
    for factor in Qpoly.factor():
      #if the current polynomial has "1*y" as a monomial, and
      #it has a degree for y of 1, and
      #it has a degree for x less than the rank
      if (1,y) in factor[0] and factor[0].degree(y) == 1 and factor[0].degree(x) < self._rank:
        #factor[0] is of the form y - f(x), get poly = f(x)
        poly = y-factor[0]
        #get the list of coefficients for x, cannot use .list() since this still has a type of multivariate polynomial
        coeff = [poly.monomial_coefficient(x^d) for d in range(poly.degree(x)+1)]
        #pad the coefficients with zeroes up to rank, make it a vector, and add to result list
        results.append(vector(self._base_ring,coeff+[0]*(self._rank-len(coeff))))
    
    if len(results) == 0:
      return None
    elif len(results) == 1:
      return results[0]
    else:
      return results

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
# # print "n:",RS1._length,", k:",RS1._rank,", t:",RS1._t
# msg = vector(GF(13),[5,1,0])
# cw = RS1.encode(msg)
# error = vector(GF(13),[1,0,0,0,6,0,0])
# received = cw + error
# print received
# print RS1.eval_encode(msg)

# print "bw decode, no error:\n",RS1.bw_decode(cw)
# print "bw decode, error(s):\n",RS1.bw_decode(received)
# print "eea decode, no error:\n",RS1.eea_decode(cw)
# print "eea decode, error(s):\n",RS1.eea_decode(received)
# print "gs decode, error(s):\n",RS1.gs_decode(received, RS1._t, 1, 1)


