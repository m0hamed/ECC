load("ReedSolomonCode.sage")

class QarySymmetricChannel:
  def __init__(self, code, error_chance):
    self._code = code
    self._base_ring = code._base_ring
    self._message_length = code._rank
    self._codeword_length = code._length
    self._error_chance = error_chance

  def random_message(self):
    return vector(self._base_ring,[randint(0,int(self._base_ring.last()))for _ in range(self._message_length)])

  def random_error(self):
    return vector(self._base_ring,[randint(1,int(self._base_ring.last())) if random()<self._error_chance else 0 for _ in range(self._codeword_length)])

  def transmit(self, n):
    bit_errors = 0
    block_errors = 0
    for _ in range(n):
      msg = self.random_message()
      codeword = self._code.encode(msg)
      error = self.random_error()
      received = codeword + error
      decoded = self._code.eea_decode(received)
      if msg != decoded:
        block_errors += 1
        if decoded is None:
          bit_errors += self._message_length
        else:
          bit_errors += sum([msg[i] != decoded[i] for i in range(self._message_length)])

    print "bit errors:",bit_errors
    print "block errors:",block_errors
    print "bit error rate:",float(bit_errors)/(float(n)*self._message_length)
    print "block error rate:",float(block_errors)/float(n)

RS1 = ReedSolomonCode(7, 3, GF(13))
Ch1 = QarySymmetricChannel(RS1,0.5)
Ch1.transmit(10000)
