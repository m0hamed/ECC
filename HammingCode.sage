load("BasicLinearCode.sage")

class HammingCode(BasicLinearCode):
  # creates the mth binary Hamming code, calculates one possible generator matrix, and its corresponding parity check matrix
  # also calculates the bit positions of the message that generate each parity bit
  def __init__(self, m, base_ring=GF(2)):
    self._base_ring = base_ring
    self._m = m;
    self._length = (base_ring.order()^m-1)/(base_ring.order()-1)
    self._rank = self._length-m

    gen_rows = []
    #lexicographically generate all vectors of length m and Hamming distance i>=2 for binary codes
    # if False and base_ring == GF(2):
    #   for i in range(2,m+1):  #iterate over all Hamming distances
    #     base = [0]*(m-i)+[1]*i #create a list of field elements
    #     while True:
    #       new_row = copy(base)
    #       gen_rows.insert(0,new_row)
    #       if not next_permutation(base):
    #         break

    #generate all pairwise linearly independent vectors of length m
    for vector in VectorSpace(base_ring,m):
      #skip all vectors of Hamming distance < 2 for binary codes, those will be added as an identity matrix to the parity check matrix
      #for qary codes, skip just the zero vector
      if self.distance(vector,[0]*m)<2:
        continue

      #check for linear dependence with all previously selected vectors, skip those that are dependant
      #this check is not necessary for binary vectors, as different binary vectors have zeroes in different positions
      linear_dependence = False

      if (base_ring != GF(2)):
        for row in gen_rows:
          if self._are_linearly_dependant(row,vector):
            linear_dependence = True
            break

      if not linear_dependence:
        gen_rows.append(vector)

    # put all the vectors in a matrix
    A = matrix(self._base_ring,(gen_rows));

    #create the generator and parity check matrices
    #generator matrix will be in echelon form by adding an identity matrix an the beginning
    self._generator_matrix = matrix.identity(self._base_ring,self._rank).augment(A)
    self._parity_check_matrix = (-A.transpose()).augment(matrix.identity(self._base_ring,m))

    #for binary Hamming codes, create a dictionary of sets, each set containing the message bits that generate the parity bits
    if (base_ring == GF(2)):
      self._parity_positions = {}
      for parpos in range(self._m):
        self._parity_positions[parpos] = set()
        for i  in range(self._rank):
          if self._generator_matrix[i,self._rank+parpos] != 0:
            self._parity_positions[parpos].add(i)

  #check if two vectors are linearly dependant, by dividing them one by one and checking the quo
  def _are_linearly_dependant(self,x,y):
    #return matrix([x,y]).echelon_form()[1] == 0
    div = None
    for i in range(self._m):
      if x[i] == 0 and y[i] == 0:
        continue
      elif x[i] == 0 or y[i] == 0:
        return False
      divNew = x[i].quo_rem(y[i])
      if divNew[1] != 0:
        return False;
      divNew = divNew[0]
      if div is not None and div != divNew:
        return False
      div = divNew
    return True

  def parity_check_matrix(self):
    return self._parity_check_matrix

  def minimum_distance(self):
    return 3

  def encode(self,msg):
    #for binary Hamming codes, we can calculate parity bits directly by adding the correct message bits
    if self._base_ring == GF(2):
      msg = list(msg)
      paritybits = self._calculate_parity(msg)

      return self.get_vector(msg+paritybits)
    #for higher fields, the codeword can be generated by multiplying with the generator matrix, which is in systematic form
    return self.get_vector(msg)*self._generator_matrix

  #calculate the parity bits for binary codes, by adding the correct message bits
  def _calculate_parity(self,msg):
    paritybits = [0]*self._m

    for parbit in range(self._m):
      for parity_constituent in self._parity_positions[parbit]:
        paritybits[parbit] += msg[parity_constituent]
    return paritybits

  #unencoding involves simply returning the first rank symbols of the codeword
  def unencode(self,msg):
    return self.get_vector(msg[:self._rank])

  #decodes only binary Hamming codes, through the syndrome of the received word
  def decode(self,word):
    word = self.get_vector(word)
    #for binary codes, we can check the message correctness faster by generating the parity bits, and comparing them to the received ones
    if self._base_ring == GF(2) and self.is_codeword_binary(word):
      return word

    parity_matrix_transposed = self._parity_check_matrix.transpose()

    syndrome = word*parity_matrix_transposed
    #if the syndrome is the zero vector, then it is a correct codeword (check applies only to qary Hamming codes)
    if syndrome == 0:
      return word

    #the syndrome vector for binary hamming codes is a column of the parity check matrix, the position of which signified the error position
    #by flipping that bit, we get the correct decoding (only applies for 1 error) 
    if (self._base_ring == GF(2)):
      [incorrect_message_pos] = [i for i,x in enumerate(parity_matrix_transposed) if x == syndrome]
      word[incorrect_message_pos] += 1
      return word
    else:
      [incorrect_message_pos] = [i for i,x in enumerate(parity_matrix_transposed) if self._are_linearly_dependant(x,syndrome)]
      # get position of first non zero element in list
      i = next((i for i, x in enumerate(syndrome) if x), None)
      word[incorrect_message_pos] -= syndrome[i].quo_rem(parity_matrix_transposed[incorrect_message_pos,i])[0]
      
      return word

  #check if a received word for a binary Hamming code is a codeword
  #for binary codes, we can check the message correctness faster by generating the parity bits, and comparing them to the received ones
  def is_codeword_binary(self,word):
    word_paritybits = self.get_vector(word[-self._m:])
    gen_paritybits = self.get_vector(self._calculate_parity(self.unencode(word)))

    if word_paritybits == gen_paritybits:
      return True
    return False
    
#test of hamming generation
# m = 2

# F.<a> = GF(9)
# h = HammingCode(m, F)
# print "m=", m, ", Hamming(", h._rank,",", h._length, ")"
# print F.list()
# print "Generator:\n", h._generator_matrix.str()
# print "Parity:\n", h._parity_check_matrix

# binaries = list(itertools.product([0, 1], repeat=h._rank))

# for msg in binaries:
# msg = binaries[-3]
# msg = [0,1]
# msg = [0, a, a + 1, 2*a + 1, 2, 2*a, 2*a + 2, a + 2]
# m = vector(F,msg)
# cw = h.encode(m)
# err = vector(F,[0,0,0,0,0,0,0,0,1,0])#,0,0,0,0,0])
# rw = cw + err
# # # # if (h.encode_systematic(m)!=h.encode(m) or h.unencode_systematic(cw)!=h.unencode(cw)):
# print "message",m
# print "encoded",cw
# print "error", err
# print "received", rw
# print h.encode_systematic(m)!=h.encode(m) or h.unencode_systematic(cw)!=h.unencode(cw)
# decoded = h.decode(rw)
# print "Decode",decoded
# print "unencode",h.unencode(decoded)

# # print "Generator Echelon:\n", h._generator_matrix.echelon_form()
# print "H*G^T:\n", h._parity_check_matrix*h._generator_matrix.transpose()
