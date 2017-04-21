import scipy.optimize as opt
import math
import random as re
import os
import matplotlib.pyplot as plt


def t_1(x1,event):
	return event[-1]*x1

def t_2(x1,x2,event):
	t_n = event[-1]
	s = 0.0
	for i in xrange(len(event)-1):
		a = math.exp(-x2*(t_n-event[i]))
		s+=a-1

	s=s*(x1/x2)
	return s

def t_3(x1,x2,x3,event):
	R = [0.0 for i in xrange(len(event))]
	for i in xrange(1,len(R)):
		R[i] = (math.exp(-x3*(event[i]-event[i-1])))*(1+R[i-1])
	s=0.0	
	for i in xrange(len(R)):
		s+=math.log(x1+x2*R[i])
	return s

def plot_event(x):
	y = [1 for i in xrange(len(x))]
	return x,y

def compare(event_or,l_p,s):
	event_gen = []
	fs = open("event_sequence_hawkes_univ")
	for line in fs:
		event_gen.append(l_p+float(line.strip()))
	event_real = event_or[s:s+50]
	#print event_or[:s+50]
	#print len(event_gen),len(event_real)
	#plt.xlim([l_p-1,l_p+30])
	plt.ylim([0,2])
	#plt.xscale('log')
	x,y = plot_event(event_gen)
	plt.stem(x,y,'b',markerfmt='bo',label='synt')
	x,y = plot_event(event_real)
	plt.stem(x,y,'r',markerfmt='ro',label='real')
	plt.show()
	
def obtain_node_arrivals(fname):
	node_arr = {}
	node = []
	fs = open(fname)
	for line in fs:
		t,u,v = map(int,line.strip().split("\t"))
		if u not in node_arr:
			node_arr[u] = []
			node_arr[u].append(t)
		else:
			if node_arr[u][-1]!=t:
				node_arr[u].append(t)
		if v not in node_arr:
			node_arr[v] = []
			node_arr[v].append(t)
		else:
			if node_arr[v][-1]!=t:
				node_arr[v].append(t)	

	fs.close() 

	for key in node_arr:
		node.append(key)

	return node_arr,node

def select_node_arrivals(node_arr,node):
	event = [0.0]
	event_or = [0.0]
	u = node[1]
	t = int(len(node_arr[u])*0.8)
	#for i in xrange(t):
	#	event.append(node_arr[u][i])
	i_a_t = []
	for i in xrange(1,len(node_arr[u])):
		a = float(node_arr[u][i] - node_arr[u][i-1])/3600
		i_a_t.append(a)

	for i in xrange(1,len(node_arr[u])):
		a = event_or[i-1] + i_a_t[i-1] 
		event_or.append(a)
		if i<t:
			event.append(a)
	return event,event_or

def read_arrivals(fname):
	event = [0.0]
	event_or = [0.0]
	fs = open(fname)
	intm = [] 

	for line in fs:
		intm.append(int(line.strip()))
	fs.close()

	t = int(len(intm)*0.8)

	i_a_t = []
	for i in xrange(1,len(intm)):
		a = float(intm[i] - intm[i-1])/3600
		i_a_t.append(a)

	for i in xrange(1,len(intm)):
		a = event_or[i-1] + i_a_t[i-1] 
		event_or.append(a)
		if i<t:
			event.append(a)
	return event,event_or

if __name__=="__main__":
	
	#node_arr,node = obtain_node_arrivals("combined_data_mathoverflow")  # this routine is to be used when the input is a time stamped edge list	
		
	#event,event_or = select_node_arrivals(node_arr,node) # to be used in conjunction with obtain_node_arrivals

	event,event_or = read_arrivals("node_arrivals")	# to be used if the input is arrival times only

	#print event
	#print event_or
	#p_t_g = len(event_or) - len(event)  # number of points to generate
	p_t_g = 50
	l_p = event[-1]	 # time of the last event from which to generate
	
	print p_t_g,l_p

	x = [0.5,0.3,0.5]
	bnds = ((0.001,None),(0.001,None),(0.001,None))
	event.sort()	
	fun = lambda x: t_1(x[0],event) - t_2(x[1],x[2],event) - t_3(x[0],x[1],x[2],event)
	soln = opt.minimize(fun,x,method='L-BFGS-B',bounds=bnds)
	print soln.x[0], soln.x[1], soln.x[2]
	
	os.system("./univar_hawkes "+str(p_t_g)+" "+str(soln.x[0])+" "+str(soln.x[1])+" "+str(soln.x[2]))

	print "generated events..."
	
	compare(event_or,l_p,len(event))			
