import matplotlib.pyplot as plt
x = []
b = []
t = []
D = []

fs = open("out_put")
for line in fs:
	s,b1,n = map(float,line.strip().split(" "))
	x.append(s)
	b.append(b1)

fs.close()

fs = open("lambda_t")
for line in fs:
	s,t1 = map(float,line.strip().split(" "))
	D.append(s)
	t.append(t1)

fs.close()

plt.step(x,b,color='r')
plt.step(D,t,color='b')
plt.show()
