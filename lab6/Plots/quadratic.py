import numpy as np
import matplotlib.pyplot as plt

# 3*pow(x,2) + 2*x - 6

intervals = [5, 10, 20, 30, 40, 50]

rectangle = [165, 141.875, 130.781, 127.153, 125.352, 124.275]

trapezoidal = [122.5, 120.625, 120.156, 120.069, 120.039, 120.025]

simpson = [120, 120, 120, 120, 120, 120]

value = 120


plt.plot(intervals, rectangle, 'ro-', label='Rectangle rule')
plt.plot(intervals, trapezoidal, 'bo-', label='Trapezoidal rule')
plt.plot(intervals, simpson, 'go-', label='Simpson\'s rule')
plt.axhline(y=value, color='k', linestyle='-', label='Exact value')

plt.title('Comparison of quadratures for quadratic function')
plt.legend()
plt.xlabel('Number of intervals')
plt.xticks(np.arange(0, 51, step=5))
plt.ylabel('Integral Value')
plt.yticks(np.arange(116, 136, step=2))
plt.ylim(115, 135)
plt.grid(True)
plt.savefig('quadratic.png')
plt.show()

