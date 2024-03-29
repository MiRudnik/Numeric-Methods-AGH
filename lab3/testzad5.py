import numpy as np

A = np.matrix([[0.0001, -5.0300, 5.8090, 7.8320],
               [2.2660, 1.9950,  1.2120, 8.0080],
               [8.8500, 5.6810,  4.5520, 1.3020],
               [6.7750, -2.253,  2.9080, 3.9700]])

b = np.matrix([9.5740, 7.2190, 5.7300, 6.2910]).transpose()

x = np.linalg.solve(A, b)

# Checking
print(np.allclose(np.dot(A, x), b))

print(x)
