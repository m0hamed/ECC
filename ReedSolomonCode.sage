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

  def minimum_distance(self):
    return self._min_dist

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

  def _get_decode_params(self, tau=None):
    # get n,k for easy coding
    n, k = self._length, self._rank
    # caluculate johnson radius
    johnson = n - sqrt(n*(k-1))

    #check if specified tau is too large
    if tau is not None and tau >= floor(johnson):
      print "Specified tau too large!"

    # set best tau possibly overriding incorrect tau
    if tau is None or tau >= floor(johnson):
      tau = floor(johnson)

    # if tau is still too big (integer value of johnson)
    if tau >= johnson:
      tau -= 1

    # flag to break out of loop
    stop = False
    # l up to 30, just some arbitrary number
    for l in xrange(1, 30):
      if stop:
        break
      # s should be less than l
      for s in xrange(1, l+1):
        bound = n*(2*l-s+1)/(2*(l+1))-l*(k-1)/(2*s)
        if tau < bound:
          stop = True
          break
    # if no good params found
    else:
      # return half min dist and 1,1 for normal non-list decoding
      return int((n-k+1)/2), 1, 1

    return tau, s, l

  # Guruswami-sudan decoding algorithm
  def gs_decode(self, received_word, tau=None, s=None, l=None, interpolate="knh", roots="simple", timing=False):

    if tau is None or s is None or l is None:
      tau, s, l = self._get_decode_params(tau)

    # print "using tau=%i, s=%i, l=%i" % (tau, s, l)
    # get interpolation polynomial
    if timing:
      start = time.clock()
    if interpolate == "knh":
      Qpoly = self._knh_interpolate(received_word, tau, s, l)
    elif interpolate == "book":
      Qpoly = self._interpolate(received_word, tau, s, l)
    else:
      return None
    if timing:
      end = time.clock()
      interpolation_time = end-start

    # print Qpoly

    # find roots
    if timing:
      start = time.clock()
    if roots == "simple":
      words = self._find_roots(Qpoly)
    elif roots == "book":
      words = self._find_roots_book(Qpoly)
    elif roots == "roth":
      words = self._find_roots_roth(Qpoly)
    else:
      return None
    if timing:
      end = time.clock()
      root_finding_time = end-start

    # print "decode results:", words

    #return codewords ...
    # codewords = [self.eval_encode(word) for word in words]
    # print [codeword for codeword in codewords if self.distance(codeword,received_word) <= tau]

    #.. or return unencoded words
    # return only the words whose codewords are in the required distance range

    # print "distances:", [self.distance(self.eval_encode(word),received_word) for word in words]
    results = [word for word in words if self.distance(self.eval_encode(word),received_word) <= tau]

    if timing:
      return (results, interpolation_time, root_finding_time)
    return results

  def _check_multiplicity(self, poly, encoded_word, s):
    for alpha, enc in zip(self.alphas, encoded_word):
      for bound in xrange(s):
        for dx in xrange(bound+1):
          dy = bound-dx
          if poly(x+alpha, y+enc).monomial_coefficient(x^dx*y^dy) != 0:
            print "Multiplicity test failed"
            print "dx=%i, dy=%i, poly=%s" % (dx, dy, str(poly(x+alpha, y+enc)))
            return False

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

  # takes a sage polynomial and returns a list of coefficients
  def _make_list(self, Qpoly, tau, s, l):
    R.<x,y> = self._base_ring[]
    Qlist = []
    for a in xrange(l+1):
      # b goes up to and including l_a
      for b in xrange(s*(self._length - tau) - a*(self._rank - 1) - 1 + 1):
        Qlist.append(Qpoly.monomial_coefficient(y^a*x^b))
    return Qlist

  # Interpolates Q(x,y) by solving a system of linear equations
  def _interpolate(self, received_word, tau, s, l):
    # matrix for the linear equations
    m = []
    for alpha, rec in zip(self.alphas, received_word):
      # h goes up to but excluding s
      for h in xrange(s):
        # r goes up to but excluding s-h
        for r in xrange(s-h):
          # row of coefficients representing a linear equation
          row = []
          # a goes up to and including l
          for a in xrange(l+1):
            # b goes up to and including l_a
            for b in xrange(s*(self._length - tau) - a*(self._rank - 1) - 1 + 1):
              # the coefficent of Q_ab
              coeff = 0
              # coeff is not 0 if condition met
              if a >= h and b >= r:
                coeff = binomial(a,h)*binomial(b,r)*(rec^(a-h))*(alpha^(b-r))
              row.append(coeff)
          m.append(row)
    # make a sage matrix out of it
    m = matrix(self._base_ring, m)
    # get one of the rows of the right kernel matrix as a list (polynomial coefficients)
    Qlist = list(m.right_kernel_matrix()[-1])
    #print Qlist, len(Qlist)
    # convert coefficients into sage polynomial
    return self._make_poly(Qlist, tau, s, l)

  def _knh_interpolate(self, received_word, tau, s, l):
    R.<x,y> = self._base_ring[]
    B = [y^i for i in xrange(l+1)]
    for alpha, rec in zip(self.alphas, received_word):
      Hs = [self._build_H_matrix(b, alpha, rec, s) for b in B]
      for bound in xrange(s):
        for dx in xrange(bound+1):
          dy = bound-dx
          minwdeg = 10000
          minbj = None
          minHj = None
          for bj, Hj in zip(B, Hs):
            if Hj[dx][dy] != 0:
              wdeg = bj.weighted_degree({x:1, y:self._rank-1})
              if wdeg < minwdeg:
                minwdeg = wdeg
                minbj = bj
                minHj = Hj
          Bnew = []
          Hsnew = []
          for bj, Hj in zip(B, Hs):
            if bj == minbj:
              continue
            Bnew.append(bj - (Hj[dx][dy]/minHj[dx][dy])*minbj)
            Hsnew.append(Hj - (Hj[dx][dy]/minHj[dx][dy])*minHj)
          Bnew.append((x-alpha)*minbj)
          Hsnew.append(self._shift_down_H_matrix(minHj, s))
          B = Bnew
          Hs = Hsnew
    return min(B, key=lambda bj: bj.weighted_degree({x:1, y:self._rank-1}))

  # shifts the hasse derivative matrix down a row
  def _shift_down_H_matrix(self, H, s):
    # transform back into mutuable list
    H = list(H)
    # add new row
    H.insert(0, [0]*s)
    # remove last row
    H.pop()
    # change back into matrix and return
    return matrix(self._base_ring, H)

  # builds a matrix of Hasse derivatives
  def _build_H_matrix(self, poly, x0, y0, s):
    # set up x and y
    R.<x,y> = self._base_ring[]
    # initial H matrix of all 0s
    H = [[0]*s for _ in xrange(s)]
    for dx in xrange(s):
      for dy in xrange(s-dx):
        # set correct hasse derivative value
        H[dx][dy] = poly(x+x0, y+y0).monomial_coefficient(x^dx*y^dy)
    return matrix(self._base_ring, H)


  #factorizes the Q polynomial, finds factors of type "y-f(x)", and returns the list of f(x)
  def _find_roots(self, Qpoly):
    R.<x,y> = self._base_ring[]
    results = []
    #factorize the multivariate polynomial, and iterate over the factors
    #factor[0] holds the factor, factor[1] the number of occurences of said factor
    for factor in Qpoly.factor():
      #it has a degree for y of 1, and
      #it has a degree for x less than the rank
      if factor[0].degree(y) == 1 and factor[0].degree(x) < self._rank:
        #divide the polynomial with the monomial coefficient of y, getting factor[0] to the form y - f(x), then get poly = f(x)
        poly = y-factor[0]/factor[0].monomial_coefficient(y)
        #get the list of coefficients for x, cannot use .list() since this still has a type of multivariate polynomial
        coeff = [poly.monomial_coefficient(x^d) for d in range(poly.degree(x)+1)]
        #pad the coefficients with zeroes up to rank, make it a vector, and add to result list
        results.append(vector(self._base_ring,coeff+[0]*(self._rank-len(coeff))))

    return results

  #factorizes the Q polynomial, finds factors of type "y-f(x)", and returns the list of f(x)
  #uses book algorithm
  #does not work if q is extension field, factorization of extension of extension field not implemented in sage
  def _find_roots_book(self, Qpoly, debug=False):
    #multivariate polynomial ring over GF(q), for the Qpoly
    R.<x,y> = self._base_ring[]
    #univariate polynomial ring over GF(q), for operation of coefficients of y^i in the multivariate polynomial

    # Z.<z> = self._base_ring[]
    Z = PolynomialRing(self._base_ring, 'x')
    z = Z.gen()

    #book example recreation irreducible polynomial, for debug purposes
    #h = 1 + a^11*z + a^6*z^2 + a^8*z^3 + a^11*z^4 + a^12*z^5 + a^13*z^6 + z^7

    #create a polynomial of rank degree, that is irreducible over Fq
    h = Z.irreducible_element(self._rank)

    #quotient ring in GF(q^k), to map the coefficients of y^i to field elements instead of univariate polynomials of x
    U.<u> = self._base_ring.extension(h)
    #univariate polynomial ring with coefficients in GF(q^k), to hold the polynomial of y to factorize
    P.<p> = U[]

    if debug:
      print "Irreducible polynomial:"
      extension_poly_print( self._base_ring, h)
      extension_poly_print( self._base_ring, U.modulus())

    #get a list of coefficients of powers of y, each coefficient is a polynomial of x. Also converts the coefficients to univariate
    y_coefficients = [Z(Qpoly.coefficient({y:d})) for d in range(Qpoly.degree(y)+1)]

    if debug:
      print "Qpoly coefficients of y^i:"
      for i,xc in enumerate(y_coefficients):
        if xc != 0:
          print "  (", xc, ") * y^", str(i)

    #convert each list element from a multivariate polynomial to a univariate polynomial
    # y_coefficients = [sum([poly.monomial_coefficient(x^d)*z^d for d in range(poly.degree(x)+1)]) for poly in y_coefficients]

    # if debug:
    #   print "Qpoly coefficients of y^i, converted to univariate polynomials:"
    #   for i,xc in enumerate(y_coefficients):
    #     if xc != 0:
    #       print "  (",
    #       extension_poly_print( self._base_ring, xc, False)
    #       print ") * y^", str(i)


    #get the Φ mapped elements from the coefficients, which are the coefficients modulo the polynomial generating the extension field
    y_coefficients_modulo = [coeff.quo_rem(h)[1] if coeff !=0 else 0 for coeff in y_coefficients]

    #convert them to the proper type (now each coefficient is a member of a ring, not a univariate polynomial)
    y_coefficients_modulo = [U(coeff) for coeff in y_coefficients_modulo]

    if debug:
      print "Qpoly coefficients of y^i, set in extension field:"
      for i,xc in enumerate(y_coefficients_modulo):
        if xc != 0:
          print "  (",
          extension_poly_print( self._base_ring, xc, False)
          print ") ...* y^", str(i)

    # get the univariate polynomial, with coefficients in the larger extension field
    Qpoly_univariate = sum([y_coefficients_modulo[d]*p^d for d in range(Qpoly.degree(y)+1)])

    if debug:
      print "Qpoly new normal:"
      print Qpoly_univariate

    #attempt to create a direct extension field with one primitive element, to map to and from the field with two primitive elements, and get over factorization not being implemented
    # S.<s> = self._base_ring.extension(self._rank)
    # L.<l> = S[]

    #attempt to get a polynomial in the extension field of p^(i+k), where k rank and q = p^i
    #instead of a polynomial in the extension q^k of the original extension q = p^i
    #this is done to try and circumvent the non-implementation of factorization in the "double" extension field
    #fails due to no possible coercion

    # Qpoly_univariate = sum([S(y_coefficients_modulo[d])*l^d for d in range(Qpoly.degree(y)+1)])

    results = []
    for (factor,_) in Qpoly_univariate.factor():
      #if it has a degree for y of 1
      if factor.degree() == 1:
        #divide the polynomial with the monomial coefficient of y, getting factor to the form y - f(x), then get poly = f(x)
        poly = p-factor/factor.list()[1]

        #convert poly from univariate polynomial in extension field, to an element in the extension field
        if poly != 0:
          poly = poly.list()[0]
        #convert poly from an element in the extension field to a univariate polynomial in the base field
        poly = Z(poly)

        #if it has a degree for x less than the rank, continue
        if poly.degree() >= self._rank:
          continue;

        #get the coefficients of the polynomial
        coefficients = poly.list()

        #pad the coefficients with zeroes up to rank, make it a vector, and add to result list
        results.append(vector(self._base_ring,coefficients+[0]*(self._rank-len(coefficients))))

    return results

#factorizes the Q polynomial, finds factors of type "y-f(x)", and returns the list of f(x)
  def _find_roots_roth(self, Qpoly):
    words = self._reconstruct(Qpoly,0,[[]])
    return [vector(self._base_ring,word) for word in words]


  #
  def _reconstruct(self,Qpoly,i,prev, debug=False):
    R.<x,y> = self._base_ring[]

    Z = PolynomialRing(self._base_ring, 'y')
    z = Z.gen()

    ret = []

    d = Qpoly.degree(y)
    r = 0
    while true:
      r += 1
      division_results = Qpoly.quo_rem(x)

      #what does "find the largest integer r for which Q(x,y)/x^r is still a (bivariate) polynomial" mean?
      #following ifs are all interpretations, none of them work exactly correct

      # if Qtest[0] != 0 and Qtest[1].degree(y) < 1:
      if division_results[0] != 0 and division_results[1] == 0:
      # if Qtest[0] != 0 and Qtest[0].degree(y) == d:
        Qpoly = division_results[0]
        if debug:
          print " division with x^"+str(r),division_results
      else:
        if debug:
          print "division end at r="+str(r)+": ", division_results, division_results[1].degree(y)
        break

    #evaluate polynomial with x=0. and convert to univariate
    M = Z(Qpoly(x=0))

    if debug:
      print "poly to find roots in:", M

    #iterate over roots
    if M == 0:
      return ret
    for (root,_) in M.roots():
      #add root to list of roots that represent a polynomial
      tmpprev = [j+[root] for j in prev]

      if debug:
        print "found root at pos "+str(i)+":", root, " :", tmpprev

      if i < self._rank-1:
        M_new = Qpoly(y=x*y+root)
        #recursive calling into rank depth
        ret += self._reconstruct(M_new,i+1,tmpprev)
      else:
        ret += tmpprev

    return ret

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
    return self.get_vector(g.padded_list(self._rank))

  # Peterson algorithm decoding through the Extended Euclidean Algorithm
  # works up to lagrange interpolation, fails horribly there
  def eea_decode(self, received_word):
    #calculate syndromes for the received word
    syndromes = self.get_syndromes(received_word)

    #generate the syndrome polynomial
    PF.<x> = self._base_ring[]

    syndrome_polynomial = sum([syndromes[self._t*2-1-i]*x^i for i in range(self._t*2)])
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
    #essentially reverse of self.alphas
    alphamap = {a_i:i for i, a_i in enumerate(self.alphas)}
    # creates a list of largange functions, which are the product of (x-a_j)/(a_i-a_j) for a_i != a_j
    #return the sum of the list
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
    return vector(self._base_ring, [
      sum([message[i]*alpha^i for i in range(self._rank)])
      for alpha in self.alphas])

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

#temp function, to help with debugging extension field decoding
def extension_poly_print(base_ring, poly, nl=True):
  P.<a> = base_ring
  PF.<x> = base_ring[]

  if base_ring.is_prime_field():
    print poly,
  else:
    rev = {}
    for i in range(0,base_ring.order()):
      print a^i
      rev[a^i]="a^"+str(i)
    rev[0]="0"
    rev[1]="1"

    coeffs = poly.list();
    for i,j in enumerate(coeffs):

      coeff = j
      if coeff != 0:
        print rev[coeff],
        if i == 1:
          print "*x",
        elif i > 1:
          print "*x^"+str(i),
        if i < len(coeffs)-1:
          print " + ",
  if nl:
    print


##### simple example
# RS1 = ReedSolomonCode(7, 3, GF(13))
# print "n:",RS1._length,", k:",RS1._rank,", t:",RS1._t
# msg = vector(GF(13),[5,6,1])
# cw = RS1.encode(msg)
# print cw
# error = vector(GF(13),[1,0,0,1,0,1,0])
# received = cw + error
# print received
# # print RS1.eval_encode(msg)

# # print "bw decode, no error:\n",RS1.bw_decode(cw)
# print "bw decode, error(s):\n",RS1.bw_decode(received)
# # print "eea decode, no error:\n",RS1.eea_decode(cw)
# print "eea decode, error(s):\n",RS1.eea_decode(received)
# for ip in ["knh","book"]:
#   for root in ["simple","book","roth"]:
#     print "gs decode("+ip+","+root+"), error(s):\n",RS1.gs_decode(received, interpolate=ip, roots=root)

#### book example, p.134

# F.<a> = GF(16)
# alphas = [a^i for i in range(0,15)]
# for i,alpha in enumerate(alphas):
#   print "a^"+str(i), "=", alpha

# RS2 = ReedSolomonCode(int(15),int(7),F,alphas)
# # msg = vector(F,[a^5,a^2,0])
# # cw = RS2.eval_encode(msg)
# # error = vector(F,[a^1,0,0,0,0,0,0])
# # received = vector(F,[0,0,a^11,0,a^12,a^11,0,0,0,0,0,0,a^3,0,a^7])
# received = vector(F,[1,0,0,1,0,0,1,0,0,1,0,0,1,0,0])
# print msg
# # print "bw decode, error(s):\n",RS2.bw_decode(received)
# # print "eea decode, error(s):\n",RS2.eea_decode(received)

# for ip in ["book","knh"]:
#   for root in ["simple","roth"]:
#     decoded = RS2.gs_decode(received, tau=5, s=4, l=6, interpolate=ip, roots=root)
#     print "gs decode("+ip+","+root+"), error(s):\n", decoded

# Roth paper example: http://www.cs.technion.ac.il/~ronny/PUB/rs.pdf

# RS3 = ReedSolomonCode(18, 2, GF(19), alphas=range(1,19))
# print "n:",RS3._length,", k:",RS3._rank,", t:",RS3._t, "alphas:",RS3.alphas
# msg = vector(GF(19),[18,14])
# cw = RS3.encode(msg)
# error = vector(GF(19),[11, 16, 17, 12, 17, 0, 0, 2, 14, 0, 0, 0, 3, 0, 14, 8, 11, 15])
# received = cw + error

# for ip in ["book","knh"]:
#   for root in ["simple","book","roth"]:
#     print "gs decode("+ip+","+root+"), error(s):\n",RS3.gs_decode(received, tau=12, s=1, l=4, interpolate=ip, roots=root)

#### Qpolys to find roots in from Roth paper example, something goes wrong

# R.<x,y> = GF(19)[]
# Qpoly =  4 + 12*x + 5*x^2 + 11*x^3 + 8*x^4 + 13*x^5
# Qpoly += (14 + 14*x + 9*x^2 + 16*x^3 + 8*x^4)*y
# Qpoly += (14 + 13*x + x*2)*y^2
# Qpoly += (2 + 11*x + x^2)*y^3
# Qpoly += 17*y^4
# print Qpoly
# print "root finding for Q polynomial (simple), :\n",RS3._find_roots(Qpoly)
# print "root finding for Q polynomial (book), :\n",RS3._find_roots_book(Qpoly)
# print "root finding for Q polynomial (roth), :\n",RS3._find_roots_roth(Qpoly)


# Qpoly =  8 +        12*x^2 + 9*x^3 + 8*x^4
# Qpoly += (5 + 14*x + 7*x^2 + 15*x^3 + 4*x^4)*y
# Qpoly += (12 + 12*x + 15*x*2 + 4*x^3)*y^2
# Qpoly += (9 + 10*x + 14*x^2)*y^3
# Qpoly += (13+x)*y^4
# print Qpoly
# print "root finding for Q polynomial (simple), :\n",RS3._find_roots(Qpoly)
# print "root finding for Q polynomial (book), :\n",RS3._find_roots_book(Qpoly)
# print "root finding for Q polynomial (roth), :\n",RS3._find_roots_roth(Qpoly)



