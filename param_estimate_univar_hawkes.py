import scipy.optimize as opt
import math

#event = []

#fs = open("event_sequence_hawkes_univ1")

#for line in fs:

#def t_1(x1,x2):
#	return ((1.5-x1+x1*x2)**2)

#def t_2(x1,x2):
#	return ((2.25-x1+x1*(x2**2))**2)

#def t_3(x1,x2):
#	return ((2.625-x1+x1*(x2**3))**2)


#def t_1(x1,x2):
#	return (x1+2*x2-7)**2

#def t_2(x1,x2):
#	return (2*x1+x2-5)**2

#x = [0,0]

#fun = lambda x: t_1(x[0],x[1]) + t_2(x[0],x[1]) + t_3(x[0],x[1])

#fun = lambda x: t_1(x[0],x[1]) + t_2(x[0],x[1])
#sol = opt.minimize(fun,x,method='Nelder-Mead')

#print sol.x

event = []


def t_1(x1):
	global event
	return event[-1]*x1

def t_2(x1,x2):
	global event
	t_n = event[-1]
	s = 0.0
	for i in xrange(len(event)-1):
		a = math.exp(-x2*(t_n-event[i]))
		s+=a-1

	s=s*(x1/x2)
	return s

def t_3(x1,x2,x3):
	global event
	R = [0.0 for i in xrange(len(event))]
	for i in xrange(1,len(R)):
		R[i] = (math.exp(-x3*(event[i]-event[i-1])))*(1+R[i-1])
	s=0.0	
	for i in xrange(len(R)):
		s+=math.log(x1+x2*R[i])
	return s

if __name__=="__main__":
	fs = open("event_sequence_hawkes_univ1")
	for line in fs:
		event.append(float(line.strip()))
	fs.close()
	x = [0.5,0.3,0.7]
	fun = lambda x: t_1(x[0]) - t_2(x[1],x[2]) - t_3(x[0],x[1],x[2])
	soln = opt.minimize(fun,x,method='Nelder-Mead')
	print soln.x


