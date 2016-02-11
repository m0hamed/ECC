
class BasicLinearCode:
  def __init__(self, generator_matrix):
    self._generator_matrix = generator_matrix
    self._base_ring = generator_matrix.base_ring()

  def generator_matrix(self):
    return self._generator_matrix

  def encode(self, msg):
    return msg*self._generator_matrix

  def encode_systematic(self, msg):
    return msg*self._generator_matrix.echelon_form()

  def systematic_positions(self):
    positions = []
    for row in self._generator_matrix.echelon_form():
      positions.append(list(row).index(1))
    return positions

  def unencode_systematic(self, codeword):
    return vector(self._base_ring, map(codeword.get, self.systematic_positions()))
