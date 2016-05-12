import matplotlib.pyplot as plt
import time
load("BasicLinearCode.sage")

#tests the decoding capability of both syndrome and nearest neighbors algorithms
#for a number of errors introduced to a random codeword from one up to its length
def test_decode(base_ring,rank,length,tests_per_number_of_errors):
  print "Generating linear code of order:", base_ring.order(), "rank:", rank, "length:", length
  while True:#get a code with at least t, for a nicer graph :P
    lc = BasicLinearCode(base_ring=base_ring,rank=rank,length=length)
    if (int((lc.minimum_distance()-1)/2) >= 1):
      break

  print "Minimum Distance: ", lc.minimum_distance(), "Guaranteed error correction:", int((lc.minimum_distance()-1)/2)
  codewords = lc.codewords()
  print "codewords generated:", len(codewords)
  #generate syndromes, and measure time
  start = time.clock()
  lc.generate_syndromes()
  end = time.clock()
  syndrome_generation_time = end - start

  print "syndromes generated:", lc._syndromes_generated, "syndrome space size:", len(lc._syndrome_dict) # base_ring.order()^(length-rank)-1

  syndrome_decoding_time = 0
  nn_decoding_time = 0
  results_nearest = []
  results_syndrome = []
  #for every possible number of errors
  for nerrors in range(1,length+1):
    unsuccessful_nearest = 0
    successful_nearest = 0
    unsuccessful_syndrome = 0
    successful_syndrome = 0

    #run the tests
    for test in range(0,tests_per_number_of_errors):
      codeword = codewords[randint(0,len(codewords)-1)]
      received_word = codeword + random_error(base_ring,length,nerrors)

      #measure time elapsed for syndrome decoding
      start = time.clock()
      syndrome_decoded = lc.decode(received_word, mode="syndrome")
      end = time.clock()
      syndrome_decoding_time += end-start

      if syndrome_decoded == codeword:
        successful_syndrome += 1
      else:
        unsuccessful_syndrome += 1

      #measure time elapsed for nearest neighbour decoding
      start = time.clock()
      nn_decoded = lc.decode(received_word, mode="nearest")
      end = time.clock()
      nn_decoding_time += end-start

      if nn_decoded == codeword:
        successful_nearest += 1
      else:
        unsuccessful_nearest += 1

    #calculate successful decoding rate
    results_nearest += [float(100*successful_nearest/tests_per_number_of_errors)]
    results_syndrome += [float(100*successful_syndrome/tests_per_number_of_errors)]

    #test results, comprehensive
    print "Decoding", tests_per_number_of_errors, "received words with", nerrors, "errors:"
    print "Syndrome decoding:", successful_syndrome, "successful,", unsuccessful_syndrome, "unsuccessful"
    print "NearestN decoding:", successful_nearest, "successful,", unsuccessful_nearest, "unsuccessful\n"

  #times elapsed
  print "Syndrome generation time:",syndrome_generation_time,"sec"
  print "Average Syndrome decoding time:",syndrome_decoding_time/(tests_per_number_of_errors*length),"sec"
  print "Average Nearest Neighbor decoding time:",nn_decoding_time/(tests_per_number_of_errors*length),"sec"
  return (lc,results_nearest,results_syndrome)

#creates a vector in the finite field with the defined order, of specified length with nerrors errors at random positions with random non_zero values in the finite field
def random_error(base_ring,length,nerrors):
  errors = 0
  ret = vector(base_ring,length)
  while errors < nerrors:
    errorpos = randint(0,length-1)
    if ret[errorpos] == 0:
      ret[errorpos] = randint(1,int(base_ring.last()))
      errors += 1
  return ret

##example test
# F.<a> = GF(2)
# n = 10
# k = 3
# tests = 1000
# test_decode(F,k, n, tests)