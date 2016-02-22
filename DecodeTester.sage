load("BasicLinearCode.sage")
#tests the decoding capability of both syndrome and nearest neighbors algorithms
#for a number of errors introduced to a random codeword from one up to its length
def test_decode(finite_field_order,rank,length,tests_per_number_of_errors):
  print "Generating linear code of order:", finite_field_order, "rank:", rank, "length:", length
  lc = BasicLinearCode(base_ring=GF(finite_field_order),rank=rank,length=length)
  print "Minimum Distance: ", lc.minimum_distance(), "Guaranteed error correction:", int((lc.minimum_distance()-1)/2)
  codewords = lc.codewords()
  print "codewords generated:", len(codewords)
  lc.generate_syndromes()
  print "syndromes generated:", lc._syndromes_generated, "syndrome space size:", len(lc._syndrome_dict) 

  for nerrors in range(1,length+1):
    unsuccessful_nearest = 0
    successful_nearest = 0
    unsuccessful_syndrome = 0
    successful_syndrome = 0
    for test in range(0,tests_per_number_of_errors):
      codeword = codewords[randint(0,len(codewords)-1)]
      received_word = codeword + random_error(finite_field_order,length,nerrors)
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

#creates a vector in the finite field with the defined order, of specified length with nerrors errors at random positions with random non_zero values in the finite field
def random_error(finite_field_order,length,nerrors):
  field = GF(finite_field_order);
  errors = 0
  ret = vector(field,length)
  while errors < nerrors:
    errorpos = randint(0,length-1)
    if ret[errorpos] == 0:
      ret[errorpos] = randint(1,int(field.last()))
      errors += 1
  return ret

#test_decode(3,3,10,1000)