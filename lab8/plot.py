import csv
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d # wcale nie unused!

x = []
y = []
z = []

with open('explicit_euler.csv', 'rt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        x.append(float(row[0]))
        y.append(float(row[1]))
        z.append(float(row[2]))


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z, 'r-')
plt.title("Explicit Euler")
plt.savefig("plot.png")
plt.show()
