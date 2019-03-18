import matplotlib.pyplot as plt
import numpy as np


px = [10.0, 20.0, 40.0, 45.0, 60.0, 65.0, 75.0, 80.0]
py = [40.0, 40.0, 60.0, 80.0, 90.0, 110.0, 100.0, 130.0]

x = np.arange(0.0, 100.0, 1.0)
y = 1.24249*x + 19.9022

plt.plot(x, y, 'b-')
plt.plot(px, py, 'ro')
plt.xlabel('X')
plt.ylabel('Y')
plt.title("Linear regression")
plt.grid()

plt.text(85, 120, "y = 1.24249*x + 19.9022", size=8, rotation=37.,
         ha="center", va="center",
         bbox=dict(boxstyle="round",
                   ec=(0.0, 0.0, 1.0),
                   fc=(0.8, 0.8, 0.8),
                   )
         )

plt.savefig("zad5.png")
plt.show()
