# Proof of concept - calculate epsilon (error bound)

from math import sqrt
from scipy.stats import norm

from kmerlight_theory import variance_of_absolute_estimate

delta = 0.05

r = 2**18
t = 7

F0 = 1.3 * 10 ** 7
m = 50
fi = 40000

sigma2 = variance_of_absolute_estimate(F0, fi, r, t)
alpha = delta / m

epsilon = norm.isf(alpha/2) * sqrt(sigma2) / fi
print(epsilon)
