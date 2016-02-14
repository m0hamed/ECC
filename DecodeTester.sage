#tests the decoding capability of both syndrome and nearest neighbors algorithms
#for a number of errors introduced to a random codeword from one up to its dimension
def test_decode(finite_field_order,rank,dimension,tests_per_number_of_errors):
  lc = BasicLinearCode(base_ring=GF(finite_field_order),rank=rank,dimension=dimension)
  codewords = lc.codewords()
  for nerrors in range(1,dimension+1):
    unsuccessful_nearest = 0
    successful_nearest = 0
    unsuccessful_syndrome = 0
    successful_syndrome = 0
    for test in range(0,tests_per_number_of_errors):
      codeword = codewords[randint(0,len(codewords)-1)]
      received_word = codeword + random_error(finite_field_order,dimension,nerrors)
      syndrome_decoded = lc.syndrome_coset_leader_correction(received_word)
      if syndrome_decoded == codeword:
        successful_syndrome += 1
      else:
        unsuccessful_syndrome += 1
      nn_decoded = lc.nearest_neighbour(received_word)
      if nn_decoded == codeword:
        successful_nearest += 1
      else:
        unsuccessful_nearest += 1
    print "Decoding", tests_per_number_of_errors, "received words with", nerrors, "errors:"
    print "Syndrome decoding:", successful_syndrome, "successful,", unsuccessful_syndrome, "unsuccessful"
    print "NearestN decoding:", successful_nearest, "successful,", unsuccessful_nearest, "unsuccessful\n"

#creates a vector in the finite field with the defined order, of specified dimension with nerrors errors at random positions with random non_zero values in the finite field
def random_error(finite_field_order,dimension,nerrors):
  field = GF(finite_field_order);
  errors = 0
  ret = vector(field,dimension)
  while errors < nerrors:
    errorpos = randint(0,dimension-1)
    if ret[errorpos] == 0:
      ret[errorpos] = randint(1,int(field.last()))
      errors += 1
  return ret

#test_decode(2,3,10,100)