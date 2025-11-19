import matplotlib.pyplot as plt
import numpy as np
import csv

# Datei einlesen
time1 = []
max_pressure1 = []
time2 = []
max_pressure2 = []
time3 = []
max_pressure3 = []

with open("../build/data_0.074.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        time1.append(float(row["Time"]))
        max_pressure1.append(float(row["MaxPressure"]))

with open("../build/data_0.076.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        time2.append(float(row["Time"]))
        max_pressure2.append(float(row["MaxPressure"]))

with open("../build/data_0.1.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        time3.append(float(row["Time"]))
        max_pressure3.append(float(row["MaxPressure"]))

# Plot
l = 0.01
time1 = np.exp(l*np.array(time1))
time2 = np.exp(l*np.array(time2))
time3 = np.exp(l*np.array(time3))
plt.plot(time2, max_pressure2, label='dt=0.076',linestyle='dotted')
plt.plot(time3, max_pressure3, label='dt=0.1',linestyle='--')
plt.plot(time1, max_pressure1, label='dt=0.074')

plt.xlabel("Time")
plt.ylabel("MaxPressure")
plt.title("MaxPressure")
plt.grid(True)
plt.legend()
plt.yscale('log')
plt.xscale('linear')
plt.ylim(1e-1, 5)

plt.show()