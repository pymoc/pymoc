import numpy as np


def gridit(x1, x2, f):
  n1 = len(x1)
  n2 = len(x2)
  array = np.zeros((n1, n2))
  for i in range(0, n1):
    for j in range(0, n2):
      array[i, j] = f(x1[i], x2[j])
  return array
