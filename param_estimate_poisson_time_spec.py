import os
import sys
from param_estimate_univar_hawkes import RMSE,MAE,obtain_node_arrivals,select_node_arrivals,compare,find_time_span
import numpy as np
import time

__author__ = "Sandipan Sikdar"

def estimate_poisson(ver):

	node_arr,node = obtain_node_arrivals("CollegeMsg.txt")  # this routine is to be used when the input is a time stamped edge list

	rmse_all = []

	time_f = find_time_span(node_arr)

	for i in xrange(len(node)):

		print "considering node",node[i]	
		
		event,event_or = select_node_arrivals(node_arr,node,i,time_f)

		if len(event)<100:
			print "not enough events to learn parameter"
			continue

		intensity = len(event)/float(event[-1] - event[0])
	
		iterations = 10
		rmse_sim = []

		l_p = len(event_or) - len(event)

		while(iterations>0):
			time.sleep(1)
			if int(ver)==1:
				os.system("./poisson "+str(l_p)+" "+str(intensity)+" "+sys.argv[1])

			else:
				os.system("./poisson "+str(event_or[-1])+" "+str(intensity)+" "+sys.argv[1])

			try:
				rmse = compare(event_or,l_p,len(event),0,"event_sequence_poisson")
				rmse_sim.append(rmse)
			except:
				print "not enough events"
			#print iterations,rmse
			iterations-=1
		if len(rmse_sim)>0:
			print np.mean(rmse_sim)
			rmse_all.append(np.mean(rmse_sim))
		
	return rmse_all


if __name__=="__main__":

	ver = sys.argv[1]	
	rmse_all = estimate_poisson(ver)
	print np.mean(rmse_all),np.std(rmse_all)

