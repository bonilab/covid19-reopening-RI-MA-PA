from scipy.stats import multinomial
multinom_logpmf = multinomial

def py_ll(means, death_count):
  x = death_count.astype(int)
  n = x.sum(axis=-1)
  p = means / means.sum(axis=-1,keepdims=True)
  return multinomial.logpmf(x = x, n = n, p = p).sum()

def py_ll_dangeous(means, death_count):
  "Same as above, but doens't check for valid inputs"
  x = death_count.astype(int)
  n = x.sum(axis=-1)
  p = means / means.sum(axis=-1,keepdims=True)
  return multinomial._logpmf(x = x, n = n, p = p).sum()
