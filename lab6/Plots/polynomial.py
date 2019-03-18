import numpy as np
import matplotlib.pyplot as plt

# pow(x,5) + 3.5*pow(x,4) - 2.5*pow(x,3) - 12.5*pow(x,2) + 1.5*x + 9

intervals = [5, 10, 20, 30, 40, 50]

rectangle = [39, 26.8125, 21.8613, 20.3876, 19.6842, 19.2729]

trapezoidal = [24, 19.3125, 18.1113, 17.8876, 17.8092, 17.7729]

simpson = [17.75, 17.7109, 17.7085, 17.7084, 17.7083, 17.7083]

value = 17.7083


plt.plot(intervals, rectangle, 'ro-', label='Rectangle rule')
plt.plot(intervals, trapezoidal, 'bo-', label='Trapezoidal rule')
plt.plot(intervals, simpson, 'go-', label='Simpson\'s rule')
plt.axhline(y=value, color='k', linestyle='-', label='Exact value')

plt.title('Comparison of quadratures for polynomial function')
plt.legend()
plt.xlabel('Number of intervals')
plt.xticks(np.arange(0, 51, step=5))
plt.ylabel('Integral Value')
#plt.yticks(np.arange(116, 136, step=2))
#plt.ylim(17.70, 17.82)
plt.grid(True)
plt.savefig('polynomial.png')
plt.show()

