import numpy as np
import matplotlib.pyplot as plt

# 2.5*x + 3

intervals = [5, 10, 20, 30, 40, 50]

rectangle = [52.5, 49.375, 47.8125, 47.2917, 47.0312, 46.875]

trapezoidal = [46.25, 46.25, 46.25, 46.25, 46.25, 46.25]

simpson = [46.25, 46.25, 46.25, 46.25, 46.25, 46.25]

value = 46.25


plt.plot(intervals, rectangle, 'ro-', label='Rectangle rule')
plt.plot(intervals, trapezoidal, 'bo-', label='Trapezoidal rule')
plt.plot(intervals, simpson, 'go-', label='Simpson\'s rule')
plt.axhline(y=value, color='k', linestyle='-', label='Exact value')

plt.title('Comparison of quadratures for linear function')
plt.legend()
plt.xlabel('Number of intervals')
plt.xticks(np.arange(0, 51, step=5))
plt.ylabel('Integral Value')
plt.grid(True)
plt.savefig('linear.png')
plt.show()
