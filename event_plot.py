import matplotlib.pyplot as plt

fs = open("event_sequence_hawkes_univ1")

x_axis = []
y_axis = []

for line in fs:
	a = float(line.strip())
	x_axis.append(a)
	y_axis.append(1)

fs.close()

plt.stem(x_axis,y_axis)
plt.ylim([0,2])
plt.xlim([0,15])
plt.show() 
