import csv
import math as m
import matplotlib.pyplot as plt
import scipy.fftpack as fft


data = []
with open('computed_data.csv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        re = float(row[0])
        im = float(row[1])
        data.append(m.sqrt(re*re+im*im))

for e in data:
    print(e)

f_s = 1
freqs = fft.fftfreq(len(data)) * f_s
fig, ax = plt.subplots()
ax.stem(freqs[:64], data[:64])
ax.set_xlabel('Frequency in hits/month')
ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
ax.set_ylim(-1, 15000)
plt.show()
