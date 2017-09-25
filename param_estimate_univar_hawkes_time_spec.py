import scipy.optimize as opt
import math
import random as rand
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import time

__author__ = "Sandipan Sikdar"

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

def compare(event_or,l_p,s,v,fname):
	event_gen = []
	fs = open(fname)
	for line in fs:
		event_gen.append(l_p+float(line.strip()))
	fs.close()
	if v==1:
		event_real = event_or[:s]
	else:
		event_real = event_or[s:]
	#print len(event_real),len(event_gen)
	'''
	plt.ylim([0,2])	
	x,y = plot_event(event_gen)
	plt.stem(x,y,'b',markerfmt='bo',label="synt")
	x,y = plot_event(event_real)
	plt.stem(x,y,'r',markerfmt='ro',label="real")
	plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
	plt.show()
	'''
	return RMSE(event_real,event_gen)
	
def obtain_node_arrivals(fname):
	node_arr = {}
	node = []
	fs = open(fname)
	for line in fs:
		if re.match("%.*",line):
			continue
		temp = map(int,line.strip().split(" "))
		u = temp[0]
		v = temp[1]
		t = temp[-1]
		#print u,v,t
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
	for n in node_arr:
		node_arr[n].sort()
		node.append(n)

	return node_arr,node

def find_time_span(node_arr):
	time = []
	for n in node_arr:
		time.append(node_arr[n][-1])
	span = max(time) - min(time)
	span*=0.8
	return min(time)+span	


def select_node_arrivals(node_arr,node,n,t):
	event = [0.0]
	event_or = [0.0]
	u = node[n]
	#t = int(len(node_arr[u])*0.8)
	#t /= 3600.0
	i_a_t = []
	for i in xrange(1,len(node_arr[u])):
		a = float(node_arr[u][i] - node_arr[u][i-1])
		i_a_t.append(a)
	i = 1
	#print "time",t
	for i in xrange(1,len(node_arr[u])):
		a = event_or[i-1] + i_a_t[i-1] 
		event_or.append(a)
		#print node_arr[u][i],t
		if node_arr[u][i]<t:
			event.append(a)
	for i in xrange(len(event_or)):
		event_or[i]/=3600.0
		if i<len(event):
			event[i]/=3600.0
	#print len(event),len(event_or)
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
		a = float(intm[i] - intm[i-1])/(3600)
		i_a_t.append(a)

	for i in xrange(1,len(intm)):
		a = event_or[i-1] + i_a_t[i-1] 
		event_or.append(a)
		if i<t:
			event.append(a)
	return event,event_or

def RMSE(event_synt,event_or):   #calculates root mean squared error
	i_a_t_synt = []
	i_a_t_or = []

	for i in xrange(1,len(event_synt)):
		i_a_t_synt.append(event_synt[i] - event_synt[i-1])
		i_a_t_or.append(event_or[i] - event_or[i-1])	
	
	s = 0.0
	for i in xrange(len(i_a_t_synt)-1):
		s += (i_a_t_synt[i] - i_a_t_or[i])**2
	s = s/len(i_a_t_synt)
	return math.sqrt(s)

def MAE(event_synt,event_or):  #calculates mean absolute error
	i_a_t_synt = []
	i_a_t_or = []

	for i in xrange(1,len(event_synt)):
		i_a_t_synt.append(event_synt[i] - event_synt[i-1])
		i_a_t_or.append(event_or[i] - event_or[i-1])	

	s = 0.0
	for i in xrange(len(i_a_t_synt)-1):
		s += math.fabs(i_a_t_synt[i] - i_a_t_or[i])
	s = s/len(i_a_t_synt)
	return s

def minimum(a,b):
	if a>b:
		return b
	else:
		return a

def PRMSE(event_synt,event_or):
	k = minimum(len(event_synt),len(event_or))
	s = 0.0
	for i in xrange(k):
		s += (event_synt[i] - event_or[i])**2
	T = event_or[-1]

	if len(event_synt)>k:
		for i in xrange(k,len(event_synt)):
			s += (T - event_synt[i])**2
	else:
		for i in xrange(k,len(event_or)):
			s += (T - event_or[i])**2

	return math.sqrt(s)

def estimate_univar_hawkes(ver):

	node_arr,node = obtain_node_arrivals("CollegeMsg.txt")  # this routine is to be used when the input is a time stamped edge list

	rmse_all = []

	time_f = find_time_span(node_arr)

	for i in xrange(len(node)):

		print "considering node",node[i]	
		
		event,event_or = select_node_arrivals(node_arr,node,i,time_f) # to be used in conjunction with obtain_node_arrivals

		if len(event_or)<100:
			print "not enough events to learn parameter"
			continue

		iterations = 1

		rmse_sim = []

		x = [0.5,0.3,0.5]
		bnds = ((0.001,None),(0.001,None),(0.001,None))
		event.sort()	
		fun = lambda x: t_1(x[0],event) - t_2(x[1],x[2],event) - t_3(x[0],x[1],x[2],event)
		soln = opt.minimize(fun,x,method='L-BFGS-B',bounds=bnds)
		#soln = opt.minimize(fun,x,method='Nelder-Mead')
		print "solution", soln.x[0], soln.x[1], soln.x[2]

		fs = open("intm","w")
		for obj in event:
			fs.write(str(obj))
			fs.write("\n")
		fs.close()

		while iterations>0:
		
			if ver==1:  # simulations matching the training set...................

				time.sleep(1)
			
				os.system("./univar_hawkes "+str(len(event))+" "+str(soln.x[0])+" "+str(soln.x[1])+" "+str(soln.x[2])+" 0 intm")

				#print "generated events..."
				
				rmse = compare(event_or,0,len(event),1,"event_sequence_hawkes_univ") # compares on the training data

			else:      # simulations matching the test set.........................

				p_t_g = len(event_or) - len(event)  # number of points to generate
		
				l_p = event[-1]	 # time of the last event from which to generate

				time.sleep(1)
				
				os.system("./univar_hawkes "+str(p_t_g)+" "+str(soln.x[0])+" "+str(soln.x[1])+" "+str(soln.x[2])+" 0 intm")

				#print "generated events..."
				try:	
		
					rmse = compare(event_or,l_p,len(event),0,"event_sequence_hawkes_univ") # compares on the training data
					rmse_sim.append(rmse)
				except:
					print "not enough events"
	
			iterations-=1

		if len(rmse_sim)>0:
			print "mean_error",np.mean(rmse_sim)
			rmse_all.append(np.mean(rmse_sim))
		print "--------------------"

	return rmse_all


if __name__=="__main__":
	
	ver = int(sys.argv[1])

	rmse_all = estimate_univar_hawkes(ver)

	print np.mean(rmse_all),np.std(rmse_all)
		
