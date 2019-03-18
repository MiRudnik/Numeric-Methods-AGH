import csv
import matplotlib.pyplot as plt
import numpy as np
import math


def analytical(time):
    return math.exp(-time*0.2)


x = []
y = []

with open('runge_kutta1.csv', 'rt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        x.append(float(row[0]))


with open('runge_kutta.csv', 'rt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        y.append(float(row[0]))

t = np.arange(0.01, 25.01, 0.01)
plt.plot(t, x, 'r-', label='λ=0.2')
plt.plot(t, y, 'b-', label='λ=1.0')
plt.title("Rozpad promieniotwórczy")
plt.xlabel("czas [s]")
plt.ylabel("u(t)")
plt.legend()
plt.grid()
plt.savefig("Exponential Decay.png")
plt.show()
