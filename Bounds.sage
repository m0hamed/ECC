# extract the n,d,k variables out of dictionary vdict and check the number of unknowns
def get_ndk(vdict, max_nones=1):
  variables = ["n", "d", "k"]
  values = map(vdict.get, variables)
  count_nones = sum([1 for x in values if x is None])
  if count_nones > max_nones:
    raise Exception, "Please leave at most only %i variable undefined" % max_nones
  n, d, k = values
  return n, d, k, count_nones

# Querying of the singleton bound, provide d,n,k and get a true or false
# answer if it is possible or provide just two of them and get the value of
# the third
def singleton_bound(pp=True, **kwargs):
  n, d, k, count_nones = get_ndk(kwargs)
  if count_nones == 0:
    if d <= n-k+1:
      if pp: print "A linear Code with length %i, dimension %i and maximum distance %i can exist within the singleton bound" % (n,k,d)
      return True
    else:
      if pp: print "A linear Code with length %i and dimension %i can have a distance of at most %i within the singleton bound" % (n,k,n-k+1)
      return False
  if d is None:
    if pp: print "A linear Code with length %i and dimension %i can have a distance of at most %i within the singleton bound" % (n,k,n-k+1)
    return n-k+1
  elif n is None:
    if pp: print "A linear Code with dimension %i and minimum distance %i must have a length of at least %i within the singleton bound" % (k,d,d+k-1)
    return d+k-1
  elif k is None:
    if pp: print "A linear Code with length %i and minimum distance %i can have a dimension of at most %i within the singleton bound" % (n,d,n-d+1)
    return n-d+1

# Querying of the plotkin bound, provide d,n,k and get a true or false
# answer if it is possible or provide just two of them and get the value of
# the third. This is not available for the dimension since it is hard to isolate
def plotkin_bound(pp=True, **kwargs):
  n, d, k, count_nones = get_ndk(kwargs)
  if count_nones == 0:
    if d <= n*2^(k-1)/(2^k-1):
      if pp: print "A linear Code with length %i, dimension %i and a maximum distance %i can exist within the plotkin bound" % (n,k,d)
      return True
    else:
      if pp: print "A linear Code with length %i and dimension %i can have a distance of at most %i within the plotkin bound" % (n,k,n*2^(k-1)/(2^k-1))
      return False
  if d is None:
    if pp: print "A linear Code with length %i and dimension %i can have a distance of at most %i within the plotkin bound" % (n,k,n*2^(k-1)/(2^k-1))
    return int(n*2^(k-1)/(2^k-1))
  elif n is None:
    if pp: print "A linear Code with dimension %i and minimum distance %i must have a length of at least %i within the plotkin bound" % (k,d,d*2^(k-1)/(2^k-1))
    return int(d*(2^k-1)/2^(k-1))
  elif k is None:
    raise Exception, "You can't query the plotkin bound on the parameter k"

# Querying of the griesmer bound, provide d,n,k and get a true or false answer
# if it is possible or provide just two of them and get the value of the third.
# This is not available for the dimension or minimum distance since they are
# hard to isolate
def griesmer_bound(pp=True, **kwargs):
  n, d, k, count_nones = get_ndk(kwargs)
  if count_nones == 0:
    if n >= sum([ceil(d/2**i) for i in xrange(k)]):
      if pp: print "A linear Code with length %i, dimension %i and maximum distance %i can exist within the griesmer bound" % (n,k,d)
      return True
    else:
      if pp: print "A linear Code with length %i and dimension %i can have a distance of at most %i within the griesmer bound" % (n,k,n*2^(k-1)/(2^k-1))
      return False
  if d is None:
    raise Exception, "You can't query the griesmer bound on the parameter d"
  elif k is None:
    raise Exception, "You can't query the griesmer bound on the parameter k"
  elif n is None:
    min_n = sum([ceil(d/2**i) for i in xrange(k)])
    if pp: print "A linear Code with dimension %i and minimum distance %i must have a length of at least %i within the griesmer bound" % (k,d,min_n)
    return min_n

# Querying of the hamming bound, provide d,n,k,q and get a true or false answer
# if it is possible within the bound.
def hamming_bound(pp=True, **kwargs):
  variables = ["n", "d", "k", "q"]
  values = map(kwargs.get, variables)
  count_nones = sum([1 for x in values if x is None])
  n, d, k, q = values

  if count_nones > 0:
    raise Exception, "Cannot query hamming bound with an unknown"

  if sum([binomial(n, j)*(q-1)^j for j in xrange(int((d-1)/2)+1)]) <= q^(n-k):
    if pp: print "A linear Code with length %i, dimension %i, field size %i and maximum distance %i can exist within the hamming bound" % (n,k,q,d)
    return True
  else:
    if pp: print "A linear Code with length %i, dimension %i, field size %i and maximum distance %i cannot exist within the hamming bound" % (n,k,q,d)
    return False

# Querying of the gilbert varshamov bound, provide d,n,k,q and get a true or false answer
# if it is possible within the bound.
def gilbert_varshamov_bound(pp=True, **kwargs):
  variables = ["n", "d", "k", "q"]
  values = map(kwargs.get, variables)
  count_nones = sum([1 for x in values if x is None])
  n, d, k, q = values

  if count_nones > 0:
    raise Exception, "Cannot query Gilbert Varshamov bound with an unknown"

  if q^(n-k) > sum([binomial(n-1, j)*(q-1)^j for j in xrange(d-1)]):
    if pp: print "A linear Code with length %i, dimension %i, field size %i and maximum distance %i can exist within the Gilbert Varshamov bound" % (n,k,q,d)
    return True
  else:
    if pp: print "A linear Code with length %i, dimension %i, field size %i and maximum distance %i cannot exist within the Gilbert Varshamov bound" % (n,k,q,d)
    return False

