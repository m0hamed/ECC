import itertools

class BasicLinearCode:
  # generator_matrix is of type sage matrix in a certain field
  # a random generator matrix can be generated, by giving base_ring, rank, and dimension
  def __init__(self, generator_matrix=None, base_ring=GF(2), rank=1, dimension=2):
    if generator_matrix is not None:
       self._generator_matrix = generator_matrix
    else:
       self._generator_matrix = matrix.random_echelonizable(MatrixSpace(base_ring, rank, dimension, sparse=False), rank, dimension)
    # this is the field of the matrix
    self._base_ring = self._generator_matrix.base_ring()
    self._rank, self._dimension = self._generator_matrix.dimensions()
    self._syndrome_dict = None

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
    codewords = []
    # messages are generated by the cartesian product of the set possible elements of the base field n times
    # where n is the dimension of the code (the number of rows in the generator matrix / the size of the messages).
    for msg in itertools.product(*[self._base_ring.list() for _ in xrange(self._generator_matrix.nrows())]):
      codewords.append(self.encode(vector(self._base_ring, msg)))
    return codewords

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
    for error in self.full_error_permutations():
      error = self.get_vector(error)
      syndrome = self.get_syndrome(error)
      if syndrome in self._syndrome_dict:
        if self._syndrome_dict.hamming_weight()>error.hamming_weight():
          self._syndrome_dict[syndrome] = error
      else:
        self._syndrome_dict[syndrome] = error

  #generate all errors in 
  #error model for Finite Fields over GF(2): each error is an 
  def full_error_permutations(self):
    error_correcting_capacity = int((self.minimum_distance()-1)/2)
    generators = [self.error_permutations(dist) for dist in range(1,error_correcting_capacity+1)]
    ret = []
    for generator in generators:
      for base_set in generator:
        base_set.extend([0] * (self._dimension - len(base_set)))
        ret += set(itertools.permutations(base_set))
    return ret


  def error_permutations(self, sum):
    if sum > 0:
      for sublist in self.error_permutations(sum-1):
        if len(sublist) < self._dimension: #constrict errors to contain up to dimension elements
          yield [1] + sublist
        if sublist and (len(sublist) < 2 or sublist[1] > sublist[0]):
          if sublist[0] < self._base_ring.last(): #constrict errors to values of the finite field
            yield [sublist[0] + 1] + sublist[1:]
    else: #stop recursion
      yield []