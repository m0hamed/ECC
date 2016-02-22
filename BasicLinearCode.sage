import itertools

class BasicLinearCode:
  # generator_matrix is of type sage matrix in a certain field
  # a random generator matrix can be generated, by giving base_ring, rank, and length
  def __init__(self, generator_matrix=None, base_ring=GF(2), rank=1, length=2):
    if generator_matrix is not None:
       self._generator_matrix = generator_matrix
    else:
       self._generator_matrix = matrix.random_echelonizable(MatrixSpace(base_ring, rank, length, sparse=False), rank, length)
    # this is the field of the matrix
    self._base_ring = self._generator_matrix.base_ring()
    self._rank, self._length = self._generator_matrix.dimensions()
    self._syndrome_dict = None
    self._codewords = None

  # return generator matrix
  def generator_matrix(self):
    return self._generator_matrix

  # simple encode by multiplying with the generator matrix.
  # Assumes msg is a vector in the same base field
  def encode(self, msg):
    msg = self.get_vector(msg)
    return msg*self._generator_matrix

  # uses the reduced form of the generator matrix to encode
  # Assumes msg is a vector in the same base field
  def encode_systematic(self, msg):
    msg = self.get_vector(msg)
    return msg*self._generator_matrix.echelon_form()

  # returns the first 1 in each row of the reduced form of the generator matrix.
  # This basically corresponds to the indices in the code word that are the same as the message
  def systematic_positions(self):
    positions = []
    for row in self._generator_matrix.echelon_form():
      positions.append(list(row).index(1))
    return positions

  # just retrieves the parts of the code word corresponding to the systematic indicies
  def unencode_systematic(self, codeword):
    return vector(self._base_ring, map(codeword.get, self.systematic_positions()))

  # returns a list of all possible code words by generating a list of all possible messages and encoding them
  def codewords(self):
    if self._codewords is None:
      self._codewords = []
      # messages are generated by the cartesian product of the set possible elements of the base field n times
      # where n is the dimension of the code (the number of rows in the generator matrix / the size of the messages).
      for msg in itertools.product(*[self._base_ring.list() for _ in xrange(self._generator_matrix.nrows())]):
        self._codewords.append(self.encode(vector(self._base_ring, msg)))
    return self._codewords

  # returns the minimum hamming weight of all the codewords except the 0 one.
  def minimum_distance(self):
    return min([codeword.hamming_weight() for codeword in self.codewords() if codeword.hamming_weight() > 0])

  # the parity check matrix is basically the right kernel of the generator matrix
  def parity_check_matrix(self):
    return self._generator_matrix.right_kernel_matrix()

  #check message to be encoded: if simple list, create a matrix
  def get_vector(self, msg):
    return vector(self._base_ring,msg)

  # returns the hamming distance between a received word and a codeword.
  def distance(self,word,codeword):
    return sum([0 if (word[i]==codeword[i]) else 1 for i in range(0,len(word))]);

  # returns the decoded word using either mode="nearest" or mode="syndrome"
  def decode(self,word,mode="nearest"):
    if mode == "nearest":
      return self.unencode_systematic(self.nearest_neighbour(word))
    elif mode == "syndrome":
      return self.unencode_systematic(self.syndrome_coset_leader_correction(word))
    else:
      return None


  # returns the nearest neighbour codeword to the received word
  def nearest_neighbour(self,word):
    word = self.get_vector(word)
    dict = {self.distance(word,codeword):codeword for codeword in self.codewords()}
    return dict[min(dict)]
    #would be more efficient if min is calculated during the distance calculation loop, and no dictionary is kept, but we only care about asymptotic complexity :P

  # returns the codeword after the coset leader error has been subtracted from the received word
  def syndrome_coset_leader_correction(self,word):
    word = self.get_vector(word)
    syndrome = self.get_syndrome(word)
    if syndrome.is_zero():
      return word
    else:
      self.generate_syndromes()
      if syndrome in self._syndrome_dict:
        return word-self._syndrome_dict[syndrome]
      else:
        return None

  # calculate syndrome for received message
  def get_syndrome(self,word):
    word = self.get_vector(word)
    ret = word*self.parity_check_matrix().transpose()
    ret.set_immutable() #using this as a dictionary key requires it to be hashable, and therefore immutable
    return ret


  # generate syndrome dictionary
  def generate_syndromes(self):
    if self._syndrome_dict is not None:
      return
    self._syndrome_dict = {}
    self._syndromes_generated = 0
    #calculate number of unique symdromes
    syndrome_dict_full_size = (int(self._base_ring.last())+1)^(self._length-self._rank)
    base_errors = [[]]
    for nerrors in range(1,self._length+1):
      #create generator of unique ordered error permutations of length nerrors
      base_errors_gen = self.error_permutations(nerrors,base_errors)
      base_errors = []
      for base_set in base_errors_gen:
        #hold generated unique ordered error permutations, as base for for next generator
        base_errors.append(base_set)
        #append zeroes up to the code's length for a complete error vector
        permutation = ([0] * (self._length - len(base_set)))
        permutation.extend(base_set)
        while True:
          self._syndromes_generated += 1
          error = self.get_vector(permutation);
          syndrome = self.get_syndrome(error)
          if syndrome not in self._syndrome_dict:
            print "syndrome table:", "{0:.2f}".format(float(len(self._syndrome_dict)/syndrome_dict_full_size)*100)," percent full\r",
            self._syndrome_dict[syndrome] = error
          if len(self._syndrome_dict) == syndrome_dict_full_size:
            return
          #run through all lexicographical permutations of the base error
          if not next_permutation(permutation):
            break

  #generator for all unique lists of nerrors length in finite field, without zero values and ordered
  def error_permutations(self, nerrors, prevlists):
    if nerrors > 0:
      for prevlist in prevlists:
        limit = (int(self._base_ring.last()) if len(prevlist) == 0 else prevlist[0])
        for i in range(1, limit+1):
          yield [i] + prevlist
    else:
      yield []

  # Querying of the singleton bound, provide d,n,k and get a true or false
  # answer if it is possible or provide just two of them and get the value of
  # the third
  @staticmethod
  def singleton_bound(**kwargs):
    variables = ["d", "n", "k"]
    values = map(kwargs.get, variables)
    count_nones = sum([1 for x in values if x is None])
    if count_nones > 1:
      raise Exception, "Please leave at most only one variable undefined"
    d, n, k = values
    if count_nones == 0:
      if d <= n-k+1:
        print "A linear Code with length %i and dimension %i can have a distance %i within the singleton bound" % (n,k,d)
        return True
      else:
        print "A linear Code with length %i and dimension %i can have a distance of at most %i within the singleton bound" % (n,k,n-k+1)
        return False
    if d is None:
      print "A linear Code with length %i and dimension %i can have a distance of at most %i within the singleton bound" % (n,k,n-k+1)
      return n-k+1
    elif n is None:
      print "A linear Code with dimension %i and minimum distance %i must have a length of at least %i within the singleton bound" % (k,d,d+k-1)
      return d+k-1
    elif k is None:
      print "A linear Code with length %i and minimum distance %i can have a dimension of at most %i within the singleton bound" % (n,d,n-d+1)
      return n-d+1

#lexicographic permutation generation
#https://www.nayuki.io/page/next-lexicographical-permutation-algorithm
def next_permutation(arr):
  # Find non-increasing suffix
  i = len(arr) - 1
  while i > 0 and arr[i - 1] >= arr[i]:
      i -= 1
  if i <= 0:
      return False

  # Find successor to pivot
  j = len(arr) - 1
  while arr[j] <= arr[i - 1]:
      j -= 1
  arr[i - 1], arr[j] = arr[j], arr[i - 1]

  # Reverse suffix
  arr[i : ] = arr[len(arr) - 1 : i - 1 : -1]
  return True

#temporary (hopefully correct) construction of m'th Hamming code, generator and parity matrices in systematic form
def hamming(m):
  length = 2^m-1
  rank = length-m
  Ir = matrix.identity(rank)
  gen_rows = []

  #untested, unbased, completely derivative determination of table (-A^T)
  for i in range(0,m-1):
    base = ([0]*(i))+([1]*(length-rank-i))
    while True:
      #why do I have to copy????
      gen_rows.insert(0,copy(base))
      if not next_permutation(base):
        break
  A = matrix(GF(2),(gen_rows));

  generator = matrix(GF(2),matrix.identity(rank).augment(A))

  parity = A.transpose().augment(matrix.identity(length-rank))

  return (length,rank,generator,parity)

#test of hamming generation
#m = 4
#h = hamming(m)
#print "m=", m, ", Hamming(", h[0], ",", h[1], ")"
#print "Generator:\n", h[2]
#print "Parity:\n", h[3]
#print "H*G^T:\n", h[3]*h[2].transpose()
