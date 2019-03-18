import numpy as np
import matplotlib.pyplot as plt

# for matrix: [16 , 3 | 11
#              7, -11 | 13]

errorJacobi = [1.545454545, 0.3181818182, 0.1844008264, 0.03796487603, 0.02200237134, 0.004529899981, 0.002625282944,
               0.0005404994296, 0.0003132439876, 6.449140921e-005, 3.737570307e-005, 7.69499769e-006, 4.459600934e-006,
               9.181531333e-007, 5.32111475e-007, 1.095523625e-007, 6.34905738e-008, 1.30715887e-008, 7.575579941e-009,
               1.559678164e-009, 9.039043958e-010, 1.860980259e-010, 1.078521716e-010, 2.22050156e-011]

errorGauss = [1.863636364, 0.3494318182, 0.04169356921, 0.004974800872, 0.000593584195, 7.08253869e-005,
              8.450756392e-006, 1.008328888e-006, 1.203119696e-007, 1.435540553e-008, 1.712860964e-009,
              2.043755165e-010, 2.438571567e-011]

errorSOR = [1.755340909, 0.2889200994, 0.003720568559, 0.0007143638979, 8.071524194e-006, 1.723882833e-006,
            3.342626231e-008, 4.052837621e-009, 1.147102413e-010, 9.250711308e-012]


plt.plot(range(1, len(errorJacobi)+1), errorJacobi, 'r', label = 'Jacobi')
plt.plot(range(1, len(errorGauss) + 1), errorGauss, 'b', label = 'Gauss-Seidel')
plt.plot(range(1, len(errorSOR) + 1), errorSOR, 'g', label = 'SOR')

plt.title('Comparison of methods solving linear equations')
plt.legend()
plt.xlabel('Iteration')
plt.xticks(np.arange(0, 25, step=2))
plt.ylabel('Max difference in values between iterations')
plt.yscale('log')
plt.grid(True)
plt.show()
