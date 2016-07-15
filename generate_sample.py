#!/usr/bin/env python

import math
import random
import numpy as np
import os
import sys

if __name__=='__main__':


   trappe=open(sys.argv[1],"w")

   mu=0.5
   sigma=0.1
   nb_events=10000

   counter=0
   for i in range(1000000):
     if i%10000==0: print " %s generated so far" % i


     if 0:
       r=np.random.normal(loc=mu, scale=sigma)
       phi=random.random()*2*math.pi
       x=r*math.cos(phi)
       y=r*math.sin(phi)
       if abs(x)<1 and abs(y)<1:
#         trappe.write(" %s  %s  \n" % (x,y))
         trappe.write(" %s  %s  \n" % (r,phi))

     x=random.random()
     y=random.random()

     eps=0.1
     fct=1.0+eps*(2.0*x)
     max_weight=1.0+2*eps

     if fct>max_weight*random.random():
#     if 1:
       trappe.write(" %s  %s  \n" % (x,y))
       counter+=1
     if counter==nb_events: break

   trappe.close()
