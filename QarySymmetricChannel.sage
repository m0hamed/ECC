load("ReedSolomonCode.sage")
load("HammingCode.sage")

#linear code with rank = length, ie with no error correction
class NoCode(BasicLinearCode):
  def __init__(self, base_ring, length=2):

    self._generator_matrix = matrix.identity(base_ring,length)
    self._base_ring = base_ring
    self._rank = length
    self._length = length

  def encode(self, msg):
    return self.get_vector(msg)

  def unencode(self, msg):
    return self.get_vector(msg)

  def decode(self, msg):
    return self.get_vector(msg)

  def minimum_distance(self):
    return 1


class QarySymmetricChannel:
  def __init__(self, code, decoding_algorithm, unencoding_algorithm, crossover):
    self._code = code
    self._decoding_algorithm = decoding_algorithm
    self._unencoding_algorithm = unencoding_algorithm
    self._base_ring = code._base_ring
    self._message_length = code._rank
    self._codeword_length = code._length
    self._crossover = crossover

  def random_message(self):
    return vector(self._base_ring,[self._base_ring.list()[randint(0,self._base_ring.order()-1)] for _ in range(self._message_length)])

  def random_error(self):
    return vector(self._base_ring,[self._base_ring.list()[randint(1,self._base_ring.order()-1)]  if random()<self._crossover else 0 for _ in range(self._codeword_length)])

  def transmit(self, n):
    print "Transmitting %d messages over (%d,%d,%d) code in F_%d, with crossover probability %f." % (n,self._codeword_length,self._message_length,self._code.minimum_distance(),self._base_ring.order(),self._crossover)
    bit_errors = 0
    block_errors = 0
    for _ in range(n):
      msg = self.random_message()
      codeword = self._code.encode(msg)
      error = self.random_error()
      received = codeword + error
      decoded = self._decoding_algorithm(received)
      if (self._unencoding_algorithm is not None):
        decoded = self._unencoding_algorithm(received)

      if msg != decoded:
        block_errors += 1
        if decoded is None:
          bit_errors += self._message_length
        else:
          bit_errors += sum([msg[i] != decoded[i] for i in range(self._message_length)])

    print "bit errors:",bit_errors
    print "block errors:",block_errors
    bit_error_rate = float(bit_errors)/(float(n)*self._message_length)
    block_error_rate = float(block_errors)/float(n)
    print "bit error rate:",bit_error_rate
    print "block error rate:",block_error_rate
    return (bit_error_rate, block_error_rate)

# NC = NoCode(F, 6)
# Ch1 = QarySymmetricChannel(NC,NC.decode,0.1)
# Ch1.transmit(1000)

# RS1 = ReedSolomonCode(14, 6, F)
# Ch1 = QarySymmetricChannel(RS1,RS1.eea_decode,0.1)
# Ch1.transmit(1000)
