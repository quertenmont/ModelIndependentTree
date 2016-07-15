#!/usr/bin/env python

import math
import random
import numpy as np
import os
import itertools
from Tkinter import *
import sys

def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x-r, y-r, x+r, y+r, **kwargs)
Canvas.create_circle = _create_circle

def _create_circle_arc(self, x, y, r, **kwargs):
    if "start" in kwargs and "end" in kwargs:
        kwargs["extent"] = kwargs["end"] - kwargs["start"]
        del kwargs["end"]
    return self.create_arc(x-r, y-r, x+r, y+r, **kwargs)
Canvas.create_circle_arc = _create_circle_arc

class histo:

  def __init__(self,nbbins, output, lower=0.0, upper=1.0):

    self.lower=lower
    self.upper=upper
    self.bars={}
    self.out=output
    self.nb_bin=nbbins
    for bin in range(self.nb_bin):self.bars[bin]=0

  def put_point(self, x_value,yvalue):

    if x_value<self.lower or x_value > self.upper: return
    bin=int((x_value-self.lower)*self.nb_bin/(self.upper-self.lower) )

#    print bin
#    print x_value 
#    print " "

    if bin > self.nb_bin-1:
       print 'got an overflow'
       return
    if bin < 0:
       print 'got an underflow'
       return

    self.bars[bin]=self.bars[bin]+yvalue

  def print_histo(self, rescale=1.0):

    trappe=open(self.out, "w")
    #firstline=str('# ')+str(median)+ "  "+ str(pvalue)+"  "+str(pvalue2)+"\n"
    for bin in range(self.nb_bin):
      xvalue=self.lower+(float(bin)+0.5)/float(self.nb_bin)*(self.upper-self.lower)
      trappe.write(str(xvalue)+"  "+str(self.bars[bin]*rescale)+"\n")
    trappe.close()

class misc:

    @staticmethod

    def get_d(tag1,tag2):
      # return distance of two binary code

      if len(tag1)!=len(tag2): print "Problem with the binary codes"
      d=len(tag1)
      for i in range(len(tag1)):
        if tag1[i]==tag2[i]: 
           d=d-1
        else :
           break # binary code differ, d should no be further decreased

#      print "tag1, tag2, d: %s  %s  %s" % (tag1,tag2,str(d))
      return d 

    def get_spTy(momenta):

       ptot=Momentum(0.0,0.0,0.0,0.0)
       for part in momenta.keys():
          if part>-1:
              ptot=ptot.add(momenta[part])
       s=ptot.dot(ptot)
       pT=math.sqrt(ptot.px**2+ptot.py**2)
       y = math.log((ptot.E+ptot.pz)/(ptot.E-ptot.pz))*0.5
       return s, pT,y

class invariants:

    def __init__(self,num_legs):

       self.num_legs=num_legs
       self.get_all_s_invariants()
       self.init_boundaries()

    def print_boundaries(self):
       print "Phase-space boundaries: "
       print "Invariant     min value        max value"
       for inv in self.list_invtuple:
          print " %s        %s            %s  " % (inv, self.min_inv[inv], self.max_inv[inv])
       print '  '
       print " %s        %s            %s  " % ('y  ', self.ymin,self.ymax)
       print " %s        %s            %s  " % ('pT  ', self.pTmin,self.pTmax)
       print " %s        %s            %s  " % ('s   ', self.stotmin,self.stotmax)
       print '  '
       for i in range(self.num_legs):
           print 'Particle %s :' % str(i)
           print " %s        %s            %s  " % ('y  ', self.ypartmin[i],self.ypartmax[i])
           print " %s        %s            %s  " % ('pT  ', self.pTpartmin[i],self.pTpartmax[i])
           print " %s        %s            %s  " % ('E   ', self.Epartmin[i],self.Epartmax[i])



    def set_baseinvariants(self,momenta):
       """ return a dictionary 
             keys = tuple     value = invariant
                    (i,j)             pi.pj
       """
       self.sij={}
       for i in range(self.num_legs-1):
           for j in range(i+1,self.num_legs):
               self.sij[(i,j)]=momenta[i].dot(momenta[j])

       s,pT, y = misc.get_spTy(momenta)
       self.s=s
       self.pT=pT
       self.y=y

    def get_all_s_invariants( self, initlist=0,length=0):
       """ return the list of all the tuples defining the invariants 
           Example with num_legs=4:
           (0,1)      <->   p0.p1         
           (0,2)      <->   p0.p2
           (0,3)      <->   p0.p3
           (1,2)      <->   p1.p2
           (1,3)      <->   p1.p3
           (2,3)      <->   p2.p3
           (1,2,3)    <->   p1.p2 + p2.p3
           (0,2,3)    <->   p0.p2 + p2.p3
           (0,1,3)    <->   p0.p1 + p1.p3
           (0,1,2)    <->   p0.p1 + p1.p2
           (0,1,2,3)  <->   p0.p1+p0.p2+p0.p3+p1.p2+p1.p3+p2.p3

       """

       # 
       if initlist==0:
           initlist=[]
           dico_length2tuple={}  # dico to record all tuples of a given length (= key) 
           list_of_inv=[]  # list to be filled with all the tuples
           for i in range (self.num_legs-1):
               initlist.append([i])
               for element in initlist:
                  fixedlength_inv=self.get_all_s_invariants( initlist=[element],length=self.num_legs-i)
                  for item2 in fixedlength_inv:
                      list_of_inv.append(tuple(item2))
           self.list_invtuple=list_of_inv

       else:
         outerlist=[]
         for item in initlist:
           if len(item)==length: continue
           # need to append one extra element 
           init_index=item[-1]+1   # smallest allowed integer for the next element
           fin_index=self.num_legs-(length-len(item))  # largest allowed integer for the next element
           newlist=[]  # to be filled with all new list obtained by appending one extra digit to item
           for i in range(init_index, fin_index+1):
               newlist.append(item+[i])

           if len(item)<length-1: # this is where the recusrion takes place
               newlist=self.get_all_s_invariants( initlist=newlist,length=length)
           outerlist+=newlist
         return outerlist

    def evaluate_inv(self,invtuple):

       """
        Evaluate the invariant associated with invtuple
        E.g. if invtuple = (i,j,k), evaluates pi.pj+pi.pk+pj.k
       """

       if not hasattr(self, 'sij'): print "Error: list_invtuple not defined"
       length_tuple=len(invtuple)
       inv=0.0
       for i in range(0,length_tuple-1):
           for j in range(i+1,length_tuple):
               inv+=self.sij[(invtuple[i],invtuple[j])]
       return inv


    def set_curr_inv(self):
       """ 
         Evaluate all the invariants defined in list_invtuple
       """

       if not hasattr(self, 'list_invtuple'): print "Error: list_invtuple not defined"

       self.curr_inv={}
       for invtuple in self.list_invtuple:
           self.curr_inv[invtuple]=self.evaluate_inv(invtuple)

    def update_boundaries(self, momenta):

       for inv in self.list_invtuple:
           if self.curr_inv[inv]>self.max_inv[inv]:
               self.max_inv[inv]=self.curr_inv[inv]
           if self.curr_inv[inv]<self.min_inv[inv]:
               self.min_inv[inv]=self.curr_inv[inv]

       if self.s<self.stotmin: self.stotmin=self.s
       if self.s>self.stotmax: self.stotmax=self.s
       if self.y<self.ymin: self.ymin=self.y
       if self.y>self.ymax: self.ymax=self.y
       if self.pT<self.pTmin: self.pTmin=self.pT
       if self.pT>self.pTmax: self.pTmax=self.pT

       for i in range(self.num_legs):
           curr_pT=momenta[i].get_pT()
           curr_y=momenta[i].get_y()

           if momenta[i].E < self.Epartmin[i]:self.Epartmin[i]=momenta[i].E
           if momenta[i].E > self.Epartmax[i]:self.Epartmax[i]=momenta[i].E

           if curr_pT < self.pTpartmin[i]:self.pTpartmin[i]=curr_pT
           if curr_pT > self.pTpartmax[i]:self.pTpartmax[i]=curr_pT

           if curr_y < self.ypartmin[i]:self.ypartmin[i]=curr_y
           if curr_y > self.ypartmax[i]:self.ypartmax[i]=curr_y


    def pass_cut(self, momenta):

       #if not hasattr(self, 'list_invtuple'): print "Error: list_invtuple not defined"
       #if not hasattr(self, 'min_inv'): print "Error: min_inv not defined"
       #if not hasattr(self, 'max_inv'): print "Error: max_inv not defined"
       #if not hasattr(self, 'curr_inv'): print "Error: curr_inv not defined"

       # cuts on the invariants: 
       for inv in self.list_invtuple:
           if self.curr_inv[inv]>self.max_inv[inv] or self.curr_inv[inv]<self.min_inv[inv]:
               #print 'fails with inv  ' + str(inv) 
               return 0

       # cuts on E, pT, y of each final state particle:
       for i in range(self.num_legs):
           curr_pT=momenta[i].get_pT()
           curr_y=momenta[i].get_y()

           if momenta[i].E < self.Epartmin[i]: return 0
           if momenta[i].E > self.Epartmax[i]: return 0

           if curr_pT < self.pTpartmin[i]: return 0
           if curr_pT > self.pTpartmax[i]: return 0

           if curr_y < self.ypartmin[i]: return 0
           if curr_y > self.ypartmax[i]: return 0

       return 1

    def init_boundaries(self):
       """ initialize the maximum value and minimum value of the invariants """

       if not hasattr(self, 'list_invtuple'): print "Error: list_invtuple not defined"
       self.min_inv={}
       self.max_inv={}

       self.stotmin=10E20
       self.stotmax=0.0
       self.pTmin=10E20
       self.pTmax=0.0
       self.ymin=10E20
       self.ymax=-10E20
#
       self.Epartmin={}
       self.Epartmax={}
       self.pTpartmin={}
       self.pTpartmax={}
       self.ypartmin={}
       self.ypartmax={}

       for inv in self.list_invtuple:
           self.min_inv[inv]=10E20
           self.max_inv[inv]=0

       for i in range(self.num_legs):
           self.Epartmin[i]=10E20
           self.Epartmax[i]=0.0
           self.pTpartmin[i]=10E20
           self.pTpartmax[i]=0.0
           self.ypartmin[i]=10E20
           self.ypartmax[i]=-10E20

class lhco_events:

    def __init__(self,path):
        """ routines to read event kinematics from the lhco event file  """

        self.path=path

    def open_file(self):
        # open the event file
        self.trappe=open(self.path, 'r')

    def close_file(self):
        # open the event file
        self.trappe.close()

    def get_next_event(self, do_print='0'):
        """ read the first event to determine the final state content """

        # fix a conventional ordering
        particle_names=['electron','positron', 'muon', 'antimuon','tau', 'antitau', 'jet', 'bjet' ]

        first_event={}
        # read the first event
        while 1:
            line=self.trappe.readline()
            if (line==""):
                return 0
            if (line[0]=="#"): continue # skip the first comment line in the event
            line=line.replace('\n','') # remove end of line statement
            line=line.split()
            if len(line)==3:   # means this is the second line in the event, contain the ID of the event
              id_event=line[1]
              self.curr_event=id_event
            else:
              curr_particle=Particle(line[1],line[2],line[3],line[4],line[5],line[6], line[7])
              if curr_particle["name"] in first_event.keys(): # here we need to use a certain ordering
                 first_event[curr_particle["name"]]=misc.ptordering(first_event[curr_particle["name"]],curr_particle)
              else:
                 first_event[curr_particle["name"]]=[curr_particle]
              if curr_particle["id"]==6: break # we have read the first event entirely -> quit the loop 

        external=-1 # start label from 0
        momenta={}
        for name in particle_names:
            if name in first_event.keys():
                counter=-1
                for part in first_event[name]:
                    counter+=1
                    external+=1
                    momenta[external]=first_event[name][counter]['p'].copy()
                    if do_print:
                        print 'particle %s : %s ' % (external, name)
#       azimuthal rotation such that momenta[0] has no py component
        phi = momenta[0].get_phi()
        for mom in momenta:
           momenta[mom]=momenta[mom].Trotate(phi)
        return momenta


class Particle(dict):
    def __init__(self, id, eta,phi,pT, mass,charge,btag):
        self["id"]=int(id)
        self["eta"]=float(eta)
        self["phi"]=float(phi)
        self["pT"]=float(pT)
        self["m"]=float(mass)
        self["charge"]=float(charge)
        self["btag"]=int(btag)
        self.get4momentum()

        self["name"]="NAN"
        if self["id"]==3 and self["charge"] >0:
            self["name"]="antitau"
        if self["id"]==3 and self["charge"] <0:
            self["name"]="tau"
        if self["id"]==2 and self["charge"] >0:
            self["name"]="antimuon"
        if self["id"]==2 and self["charge"] <0:
            self["name"]="muon"
        if self["id"]==1 and self["charge"] >0:
            self["name"]="positron"
        if self["id"]==1 and self["charge"] <0:
            self["name"]="electron"
        if self["id"]==4 and abs(btag)<1:
            self["name"]="jet"
        if self["id"]==4 and abs(btag)>0:
            self["name"]="bjet"
        if self["id"]==6 :
            self["name"]="missing"

    def printout(self):
        print "final state particle: "
        print str(self["id"])+"  "+ self["p"].nice_string()

    def get4momentum(self):
        px=self["pT"]*math.cos(self["phi"])
        py=self["pT"]*math.sin(self["phi"])
        e2eta=math.exp(2.0*self["eta"])
        costh=(e2eta-1.0)/(e2eta+1.0)
        sinth=math.sqrt(1-costh**2)
        mom=self["pT"]/sinth
        pz=mom*costh
        energy=math.sqrt(mom**2+self["m"]**2)
        self["p"]=Momentum(energy,px,py,pz)

class Momentum:
    """A class to handel 4-vectors and the associated operations """
    def __init__(self,E,px,py,pz, m=-1):
        self.px=px
        self.py=py
        self.pz=pz
        self.E=E
        self.mod2=px**2+py**2+pz**2
        self.sq=E**2-self.mod2

        if m>-1: self.m=m
#        control = E**2+self.mod2
#        if m==-1:
#          if not control:
#              self.m = 0
#          elif self.sq/control < 1e-8:
#              self.m=0.0
#          else:
#            self.m=math.sqrt(self.sq)

    def Trotate(self,phi):
       """ rotation in the transverse plane"""

       px=self.px*math.cos(phi)+self.py*math.sin(phi)
       py=-self.px*math.sin(phi)+self.py*math.cos(phi)
       E=0.0+self.E
       pz=0.0+self.pz

       return Momentum(E,px,py,pz)


    def get_pT(self):
        return math.sqrt(self.px**2+self.py**2)

    def get_phi(self):

      if self.px >0.0:
          phi= math.atan(self.py/self.px)
      elif self.px <0.0 :
          phi= math.atan(self.py/self.px) + math.pi
      elif self.py >= 0.0:
          phi= math.pi/2.0
      elif self.py<0.0 :
          phi= -math.pi/2.0

      if (phi<0.0): phi+=2.0*math.pi
      return phi

    def get_delta_R2(self, mom):

      p1x=self.px
      p1y=self.py

      p2x=mom.px
      p2y=mom.py

      denom=math.sqrt(p1x**2+p1y**2)*math.sqrt(p2x**2+p2y**2)
      temp=max(-0.999999999, (p1x*p2x+p1y*p2y)/denom)
      temp=min(0.9999999999, temp)
      delta_phi=math.acos(temp)

      return delta_phi**2+(self.get_y()-mom.get_y())**2

    def get_theta(self):

      theta=math.acos(self.pz/math.sqrt(self.px**2+self.py**2+self.pz**2))
      return theta

    def get_y(self):
        try:
           return 0.5*math.log((self.E+self.pz)/(self.E-self.pz))
        except:
           print "Problem in get_y"
           print self.nice_string()

    def norm3(self):
        """ retrun the spacial norm of the vector"""
        return math.sqrt(self.px**2+self.py**2+self.pz**2)

    def dot3(self,q):
        """ return |p|^2 (spatial components only) """
        return self.px*q.px+self.py*q.py+self.pz*q.pz

    def dot(self,q):
        """ Minkowski inner product """
        return self.E*q.E-self.px*q.px-self.py*q.py-self.pz*q.pz

    def deltaR(self,q):
        """ Minkowski inner product """
        phi1=self.get_phi()
        phi2=q.get_phi()
        y1=self.get_y()
        y2=q.get_y()
        return math.sqrt((phi1-phi2)**2+ (y1-y2)**2)


    def subtract(self,q):
        tot=Momentum(self.E-q.E,self.px-q.px,self.py-q.py,self.pz-q.pz)
        return tot

    def add(self,q):
        tot=Momentum(self.E+q.E,self.px+q.px,self.py+q.py,self.pz+q.pz)
        return tot

    __add__ = add

    def nice_string(self):
        return str(self.E)+" "+str(self.px)+" "+str(self.py)+" "+str(self.pz)

    __str__ = nice_string

    def boost(self, q):
        """ boost a vector from a frame where q is at rest to a frame where q is given 
                This routine has been taken from HELAS
        """
        qq = q.mod2

        if (qq > 1E-10*abs(q.E)):
            pq=self.dot3(q)
            m=q.m
            #if (abs(m-self.mA)>1e-6): print "warning: wrong mass"
            lf=((q.E-m)*pq/qq+self.E)/m
            if hasattr(self, 'm'):
              pboost=Momentum((self.E*q.E+pq)/m, self.px+q.px*lf,\
                              self.py+q.py*lf,self.pz+q.pz*lf,self.m)
            else:
              pboost=Momentum((self.E*q.E+pq)/m, self.px+q.px*lf,\
                              self.py+q.py*lf,self.pz+q.pz*lf)
        else:
            if hasattr(self, 'm'):
              pboost=Momentum(self.E,self.px,self.py,self.pz, self.m)
            else:
              pboost=Momentum(self.E,self.px,self.py,self.pz)

        return pboost

    def copy(self):
        copy_mom=Momentum(self.E,self.px,self.py,self.pz)
        return copy_mom

class superbox:
  """
    a set of classification trees
  """

  def __init__(self, all_events=None, nb_trees=None, nb_events=10000,nb_layers=5):
       
      if all_events and nb_trees and nb_layers and nb_events :
         self.nb_trees=nb_trees
         self.all_trees=[]
         self.nb_layers=nb_layers
         for i in range(nb_trees):
            events=[evt for evt in all_events[i*nb_events:(i+1)*nb_events]]
            self.all_trees.append(box(events))
            self.all_trees[i].generate_all_daughters(nb_layers)
            print "Tree nb %s done" % str(i+1)

  def get_boxes(self,feat1=0,feat2=1,tree=0):

      curr_tree= self.all_trees[tree]
      curr_tree.get_boxes(feat1,feat2,curr_tree.boundaries[feat1],curr_tree.boundaries[feat2])

  def read_super_tree(self,inputfile):

      trappe=open(inputfile, "r")
      tree_def=[]
      self.all_trees=[]
      count_trees=0
      line=trappe.readline()
      line=trappe.readline()
      while 1:
         line=trappe.readline()
         if line=="": 
            self.all_trees.append(box_from_file(tree_def))
            count_trees+=1
            break
         line=line.replace("\n", "")
         line=line.split()
         if line[0]=="T":
            self.all_trees.append(box_from_file(tree_def))
            count_trees+=1
            tree_def=[]
         else:
            tree_def.append(line)
      self.nb_trees=count_trees
      print "Number of trees: %s" % self.nb_trees

  def record_super_tree(self, outputfile):


      trappe=open(outputfile, "w")

      trappe.write("Number of trees: %s \n" % self.nb_trees)
      for tree in self.all_trees:
         trappe.write("T \n")
         record=tree.record_tree()
         for line in record:
            trappe.write(line)
      trappe.close

  def get_weight(self, one_point):
      
      weight=0.0
      for tree in self.all_trees:
         weight+=tree.get_weight(one_point)
      return weight/float(self.nb_trees)

  def set_deviation(self,experiment):
      for tree in self.all_trees:
         tree.set_deviation(experiment)
 
  def get_deviation(self,one_point):

      dev=0.0
      for tree in self.all_trees:
         dev+=tree.get_deviation(one_point)
      return dev/float(self.nb_trees)

  def set_cumul_deviation(self, experiment):
      # initialize the function \bar{F}
      for event in experiment:
         dev=self.get_deviation(event)
         for tree in self.all_trees:
             tree.update_cumul_deviation(event,dev) 

      for tree in self.all_trees:
         all_cumul=[]
         tag=""
         tree.get_cumul(all_cumul,tag)
         tree.all_cumul=all_cumul

  def get_event_tag(self,experiment):
     """ build the table
         evt1  tag1  tag2 ...
         evt2  tag1  tag2 ...

         where tag1, tag2, ... are the indices for the leaf with evt i in tree 1, tree2, ...
     """

     event_tags=[]

     for evt in experiment:
       dev=self.get_deviation(evt)
       tags=[]
       for i,tree in enumerate(self.all_trees):
          tag=""
          tag=tree.get_tag_event(evt,tag)
          nb_layers=len(tag)
          tags.append(int(tag,2))
       event_tags.append((dev,tags))

     self.event_tags=event_tags
     return nb_layers

  def get_two_point_fct_long(self,experiment):

     nb_layers=self.get_event_tag(experiment)
     event_tags=self.event_tags
 
     histo_d=histo(nb_layers, "histo.out", lower=0.5, upper=float(nb_layers)+0.5)

#    load an array with all distances:
     print "Loading the array with all distances:"
     all_distances=[]
     for i in range(2**nb_layers):
       tagi=(nb_layers-len(bin(i)[2:]))*'0'+bin(i)[2:]
       all_distances.append((2**nb_layers)*[0])
       for j in range(2**nb_layers):
          tagj=(nb_layers-len(bin(j)[2:]))*'0'+bin(j)[2:]
          all_distances[i][j]=misc.get_d(tagi,tagj)
     all_distances=tuple(all_distances)

  
     corr=0
     for i in range(len(event_tags)-1):
        for j in range(i+1,len(event_tags),1):
            distance=0
            for k in range(self.nb_trees):
#              test1=misc.get_d(event_tags[i][1][k],event_tags[j][1][k])
               distance+=all_distances[event_tags[i][1][k]][event_tags[j][1][k]]
            FF=event_tags[i][0]*event_tags[j][0]
            dist=float(distance)/float(self.nb_trees)
            histo_d.put_point(dist,FF)
            #if all_distances[event_tags[i][1][k]][event_tags[j][1][k]]==1:
            corr+=dist*FF
     
     del all_distances
     factor=2.0/float(len(event_tags))**2
     histo_d.print_histo(rescale=factor) 
     return corr*factor

  def get_two_point_fct(self):
 
     corr=0
     count=0
     for tree in self.all_trees:
        curr_corr=tree.get_two_point_fct(count)
        corr+=curr_corr
        print "corr for tree %s : %s " % (count,curr_corr)
        count+=1
     corr=corr/float(len(self.all_trees))
     return corr

  def get_info_tree(self, index_tree):

      self.all_trees[index_tree].print_info_tree()



class box_from_file:
  """ A tree that is extracted from a text file
  """

  def __init__(self, record):

     line=record[0]
     if line[0]=="C":
        self.cut_variable=int(line[1])
        self.cut_value=float(line[2])
        del record[0]
        self.d0=box_from_file(record)        
        self.d1=box_from_file(record)        
     elif line[0]=="L": #we have reached a leaf
        self.size=int(line[1]) 
        self.vol=float(line[2])
        self.obs=0
        self.cumul=0
        del record[0] 
     else:
        print "Unknown info in tree"
        print line

  def get_info_tree(self,nb_event_in_box, deviations,cumul):

     if not hasattr(self,"d0") and not hasattr(self,"d1"):
       nb_event_in_box.append(self.size)
       deviations.append(self.deviation)
       cumul.append(self.cumul)

     elif hasattr(self,"d0") and  hasattr(self,"d1"):
#    recursion:
       self.d0.get_info_tree(nb_event_in_box,deviations,cumul)
       self.d1.get_info_tree(nb_event_in_box,deviations,cumul)

  def get_tag_event(self,event,tag):

     if not hasattr(self,"d0") and not hasattr(self,"d1"):
       return tag
     else:
       if event[self.cut_variable]<self.cut_value:
          return self.d0.get_tag_event(event,tag+"0")
       else:
          return self.d1.get_tag_event(event,tag+"1")


  def get_cumul(self, all_cumul,tag):
     # fill the list all_cumul =[(bar{F},tag), ...]
     # where tag is the binary code of the leaf 

     if not hasattr(self,"d0") and not hasattr(self,"d1"):
       all_cumul.append((self.cumul,tag))

     elif hasattr(self,"d0") and  hasattr(self,"d1"):
#    recursion:
       self.d0.get_cumul(all_cumul,tag+"0")
       self.d1.get_cumul(all_cumul,tag+"1")

  def print_info_tree(self):

     print "Information about this tree: "
     nb_event_in_box=[]
     deviations=[]
     cumul=[]
     self.get_info_tree(nb_event_in_box, deviations, cumul)
     # compute sum of \bar F in each leaf:
     tot_cumul=0
     for cumul in self.all_cumul:
       tot_cumul+=cumul[0]
     tot_dev=0
     for dev in deviations:
       tot_dev+=dev

     for i in range(len(nb_event_in_box)) :
        print "BOX %s" %i
        print "nb events: %s " % nb_event_in_box[i]
        print "deviation :  %s " % deviations[i]
        print "cumul     :  %s " % self.all_cumul[i][0]
        print "tag       :  %s " % self.all_cumul[i][1]
        print "  "
     print "tot cumul: %s" % tot_cumul 
     print "tot dev: %s" % tot_dev 

  def get_two_point_fct(self,count):
     # compute the correlation function for the tree
    
     corr=0
     counter=0
     for i in range(len(self.all_cumul)-1):
       for j in range(i+1,len(self.all_cumul),1):
          counter+=1
          dist=misc.get_d(self.all_cumul[i][1],self.all_cumul[j][1])
#          print i, j, dist
          #if dist==1:
          corr+=float(dist)*self.all_cumul[i][0]*self.all_cumul[j][0]
          Fsq=self.all_cumul[i][0]*self.all_cumul[j][0]
     if counter!=len(self.all_cumul)*(len(self.all_cumul)-1)/2:print "warning: wrong nb of pairs of leaves"
     corr=corr/float(self.nb_exp_evts)**2*2.0
     return corr

  def update_cumul_deviation(self,one_point,dev):
     # update the cumul function \bar{F} 
     # that-is-to-say: 1) go into the leaf one_point is falling into
     #                 2) in that leaf, increment the variable "cumul" by the quanity "dev"
     if not hasattr(self, "cut_value"):
          self.cumul+=dev
          return
     else:
       if one_point[self.cut_variable]<self.cut_value:
          return self.d0.update_cumul_deviation(one_point,dev)
       else:
          return self.d1.update_cumul_deviation(one_point,dev)

  def increment_box(self,one_point):

     if not hasattr(self, "cut_value"):
          self.obs+=1
          return
     else:
       if one_point[self.cut_variable]<self.cut_value:
          return self.d0.increment_box(one_point)
       else:
          return self.d1.increment_box(one_point)

  def get_nb_points_in_tree(self,sizes):

     if not hasattr(self, "cut_value"):
          sizes.append(self.size)
          return
     else:
          self.d0.get_nb_points_in_tree(sizes)
          self.d1.get_nb_points_in_tree(sizes)

  def get_deviation(self, one_point):

     if not hasattr(self, "cut_value"):
       return self.deviation
     else:
       if one_point[self.cut_variable]<self.cut_value:
          return self.d0.get_deviation(one_point)
       else:
          return self.d1.get_deviation(one_point)

  def get_weight(self, one_point):

     if not hasattr(self, "cut_value"):
       return float(self.size)/self.vol
     else:
       if one_point[self.cut_variable]<self.cut_value:
          return self.d0.get_weight(one_point)
       else:
          return self.d1.get_weight(one_point)

  def set_deviation2(self, nb_exp_events, nb_points_in_tree):

     if not hasattr(self, "cut_value"):
       self.deviation=1.0-float(self.size)/float(self.obs)/float(nb_points_in_tree)*float(nb_exp_events)
       #print "dev in box %s" % self.deviation
       #print "obs in box %s" % self.obs
     else:
          self.d0.set_deviation2(nb_exp_events,nb_points_in_tree)
          self.d1.set_deviation2(nb_exp_events, nb_points_in_tree)

  def set_deviation(self,experiment):
 
      sizes=[]
      self.get_nb_points_in_tree(sizes)
      nb_points=0
      for size in sizes: nb_points+=size 
      self.nb_points_in_tree=nb_points

      count_evt=0
      for event in experiment:
          count_evt+=1
          self.increment_box(event)
      self.nb_exp_evts=count_evt
      self.set_deviation2(count_evt, nb_points)   


class box:

   """   
    A box is characterized by 
      - a list of events in the box (= list of integers)
      - a set of boudaries (up and down values for each feature)
          [[down_0,up_0], [down_1, up_1], ...]
      - perhaphs two daughter boxes
   """
 
   def __init__(self,all_events,list_events=None,boundaries=None):
     """
     events are the events inside the box (list of tuple)
     [ (id, feat1, feat2, ...), ... ]
     """
     self.all_events=all_events   # essentially a global variables with all the evnts
     if (list_events):
       self.events=list_events
     else:
       self.events=list(xrange(len(all_events)))
     if (boundaries):
       self.boundaries=boundaries
     else:
       self.get_boundaries()
     self.size=len(self.events)
#     self.gradients=[]

   def get_vol(self):
     """ return vol of the box""" 

     vol=1.0
     for bound in self.boundaries:
        vol=vol*(bound[1]-bound[0])
     return vol

   def get_weight(self, one_point):

     if not hasattr(self, "cut_value"):
       return float(self.size)/self.get_vol()
     else: 
       if one_point[self.cut_variable]<self.cut_value:
          return self.d0.get_weight(one_point)
       else: 
          return self.d1.get_weight(one_point)

   def record_tree(self):

     tree_record=[]
     self.record_sub_tree(tree_record)
     return tree_record

   def record_sub_tree(self,tree_record):

     if hasattr(self,"d0") and  hasattr(self,"d1"):
#    recursion:  
       tree_record.append("C %s %s \n" % (self.cut_variable, self.cut_value))
       self.d0.record_sub_tree(tree_record)
       self.d1.record_sub_tree(tree_record)
     else: 
       tree_record.append("L %s %s \n" % (self.size, self.get_vol()))

   def print_tree(self):

     print "Information about this tree: "
     nb_event_in_box=[]
     vol_boxes=[]
     self.print_sub_tree(nb_event_in_box,vol_boxes)
     for i in range(len(nb_event_in_box)) :
        print "BOX %s" %i
        print "nb events: %s " % nb_event_in_box[i]
        print "vol box :  %s " % vol_boxes[i]
        print "  "
     
   def print_sub_tree(self,nb_event_in_box,vol_boxes):
       

     if not hasattr(self,"d0") and not hasattr(self,"d1"):
       nb_event_in_box.append(self.size)
       vol_boxes.append(self.get_vol())

     elif hasattr(self,"d0") and  hasattr(self,"d1"):
#    recursion:
       self.d0.print_sub_tree(nb_event_in_box,vol_boxes)
       self.d1.print_sub_tree(nb_event_in_box,vol_boxes)
     else:
       print "Warning; problem i the classification tree (a)"

   def print_boundaries(self):
 
     nb_features=len(self.all_events[0])
     for i in range(nb_features):
        print "Lower value of feaure %s : %s " % (i,self.boundaries[i][0])
        print "Upper value of feaure %s : %s " % (i,self.boundaries[i][1])

   def get_boundaries(self):

     nb_features=len(self.all_events[0])
     self.boundaries=[]
     for i in range(nb_features):
         down=self.all_events[0][i]
         up=down
         for evt in self.events:
            if self.all_events[evt][i]>up:
               up=self.all_events[evt][i]
            if self.all_events[evt][i]<down:
               down=self.all_events[evt][i]
         up=up+(up-down)*1e-6
         down=down-(up-down)*1e-6
         self.boundaries.append([down,up])

   def get_feature(self):
     """
      Identify the feature for which the gradient is the largest
     """

     nb_features=len(self.all_events[0])
     grad_max=0
     for i in range(nb_features):
         down=self.boundaries[i][0]
         up=self.boundaries[i][1]
         medium_value=down+(up-down)/2.0
         step=(up-down)/5.0
         a=[0,0,0,0,0]
         
         for evt in self.events:
            if self.all_events[evt][i]<down+step:
               a[0]+=1
            elif self.all_events[evt][i]<down+2.0*step:
               a[1]+=1
            elif self.all_events[evt][i]<down+3.0*step:
               a[2]+=1
            elif self.all_events[evt][i]<down+4.0*step:
               a[3]+=1
            else:
               a[4]+=1

#         print "a %s , n %s " % (a, self.size)
         grad=1e-16
         if a[0]>0 and a[1]>0 :
           grad+=abs(float(a[1]-a[0])/float(a[0]+a[1]))
         if a[1]>0 and a[2]>0 :
           grad+=abs(float(a[2]-a[1])/float(a[1]+a[2]))
         if a[2]>0 and a[3]>0 :
           grad+=abs(float(a[3]-a[2])/float(a[2]+a[3]))
         if a[3]>0 and a[4]>0 :
           grad+=abs(float(a[4]-a[3])/float(a[3]+a[4]))
#         print "grad", grad
         if grad_max<grad:
            grad_max=grad
            good_var=i
#         self.gradients.append(grad)
     return good_var

   def generate_all_daughters(self,nb_layers,min_cut=0):
     """
       build up the whole classification tree
     """

     var=self.get_feature()
     self.generate_daughters(var)
     #print "Two daughters have been generated"
#      recurrence:

#     if self.d0.size>min_cut:  #criteria = min number of events
     if nb_layers>1:            #criteria = number of layers
       self.d0.generate_all_daughters(nb_layers-1)

#     if self.d1.size>min_cut:
     if nb_layers>1:            #criteria = number of layers
       self.d1.generate_all_daughters(nb_layers-1)


   def generate_daughters(self,var):
     """
       define two equal size boxes by dividing the current box along the variable "var"
     """

#     first construct a dico {feat_var : id_evt}

     dico_var2id={}
     features=[]
     for evt in self.events: 
        dico_var2id[self.all_events[evt][var]]=evt
        features.append(self.all_events[evt][var])

#    sort the feature values
     features.sort()

#    determine the cut value
     mu=float(self.size)/2.0-0.5 # since python list starts at position 0

#     std=math.sqrt(float(self.size))
     std=0.001
     cut_evt=int(round(np.random.normal(mu,std)))
     cut_val=(features[cut_evt-1]+features[cut_evt+1])/2.0

     self.cut_variable=var
     self.cut_value=cut_val


#    old
#     if self.size%2==0:
#        median_val=(features[self.size/2-1]+features[self.size/2])/2.0
#     else:
#        if random.random()<0.5:
#          median_val=(features[self.size/2-1]+features[self.size/2])/2.0
#        else: 
#          median_val=(features[self.size/2+1]+features[self.size/2])/2.0

#    create the list of events associated with daughters a and b
     list_a=[]
     list_b=[]

     for evt in self.events:
       if self.all_events[evt][var]<cut_val:
          list_a.append(evt)
       else:
          list_b.append(evt)

#    create the boundaries associated with daughters a and b
#    PAY ATTENTION THAT THE BOUNDARIES IN self.boundaries, boundaries_a and 
#    boundaries_b DO NOT HAVE THE SAME ADDRESS
     boundaries_a=[]
     boundaries_b=[]
     for item in self.boundaries:
         boundaries_a.append(list(item))
         boundaries_b.append(list(item))
     boundaries_a[var][1]=cut_val
     boundaries_b[var][0]=cut_val

#     print boundaries_a
#     print boundaries_b
#     print self.boundaries
#     print id(self.boundaries)
#     print id(boundaries_a)
#     print id(boundaries_b)
#     stop

     self.d0=box(self.all_events,list_a,boundaries_a)
     self.d1=box(self.all_events,list_b,boundaries_b)

#     print "Boundaries of d0"
#     self.d0.print_boundaries()
#     print "Boundaries of d1"
#     self.d1.print_boundaries()

     del dico_var2id
     del features

   def get_boxes(self, var1, var2,boundvar1,boundvar2):

     master = Tk()     
     totwidth=1000
     totheight=720
     w = Canvas(master, width=totwidth, height=totheight)
     w.pack()
     self.draw_boxes( var1, var2,boundvar1,boundvar2,w,totwidth,totheight)     
     mainloop() 

   def draw_boxes(self, var1, var2,boundvar1,boundvar2,w,totwidth,totheight):
     """ draw the boxes using Tkinter"""
  
     xlow=10+(totwidth-20)*(self.boundaries[var1][0]-boundvar1[0])/(boundvar1[1]-boundvar1[0])
     ylow=10+(totheight-20)*(self.boundaries[var2][0]-boundvar2[0])/(boundvar2[1]-boundvar2[0])

     xup=10+(totwidth-20)*(self.boundaries[var1][1]-boundvar1[0])/(boundvar1[1]-boundvar1[0])
     yup=10+(totheight-20)*(self.boundaries[var2][1]-boundvar2[0])/(boundvar2[1]-boundvar2[0])

#w.create_line(0, 0, 200, 100)
#w.create_line(0, 100, 200, 0, fill="red", dash=(4, 4))
#w.create_rectangle(50, 25, 150, 75, fill="blue")
#     print "Drawing one box with coordinates %s %s %s %s " % (xlow, ylow, xup, yup)
     w.create_rectangle(xlow, ylow, xup, yup)

#    recurrence
     if hasattr(self,"d0"): 
        self.d0.draw_boxes(var1,var2,boundvar1,boundvar2,w,totwidth,totheight)
        self.d1.draw_boxes(var1,var2,boundvar1,boundvar2,w,totwidth,totheight)
     else:
        colorset=["blue", "green", "red", "magenta", "black", "cyan"]
        rand=random.randint(0,5)
        color=colorset[rand]
        for evt in self.events:
           x=self.all_events[evt][var1] 
           y=self.all_events[evt][var2]
           x=10+(totwidth-20)*(x-boundvar1[0])/(boundvar1[1]-boundvar1[0]) 
           y=10+(totheight-20)*(y-boundvar2[0])/(boundvar2[1]-boundvar2[0])
#           w.create_circle(x, y, 2, fill=color)
           

def get_all_features(filename):

    print "Scanning the event file ..."
    counter_event=0
    first=1
    features_tree=[]

    if filename[-4:]=="lhco":
       events=lhco_events(filename)
       events.open_file()
       while 1:
         # read event
         momenta=events.get_next_event(do_print=0)
         if momenta==0: break
         nfinal=len(momenta)

         if  counter_event==max_nb_events: break

         counter_event+=1

         if first:
             invariants_info=invariants(nfinal)
#             print "weight: "+ str(weight_value)

         invariants_info.set_baseinvariants(momenta)
         invariants_info.set_curr_inv()
         invariants_info.update_boundaries(momenta)

#       features_list=[  ]
#       list_2inv=[]
#       for i in range(invariants_info.num_legs-1):
#           for j in range(i+1,invariants_info.num_legs):
#               list_2inv.append(invariants_info.sij[(i,j)])
#       list_2inv.sort()
#       curr_features_list+=list_2inv
#
#       for i in range(len(momenta)):
#           curr_features_list.append(momenta[i].get_phi())
#           curr_features_list.append(momenta[i].get_pT())
#           curr_features_list.append(momenta[i].get_y())

         features_tree.append((momenta[0].get_pT(),momenta[0].get_y()))

#       features.append(misc.get_feature_one_event(invariants_info,weight_value, momenta))

         first=0

    else:
       nb_features=-1
       trappe=open(filename,'r')
       while 1:
         line=trappe.readline()
         if line=="": break
         line=line.replace("\n", "")
         feature_list=line.split()
         if first: 
            nb_features=len(feature_list)
            first=0
         else:
            if nb_features!=len(feature_list):
              print "wraning: problem with the features"
         features_tree.append(tuple([float(el) for el in feature_list]))

    features_tree=tuple(features_tree)
    return features_tree

if __name__=='__main__':

#    max_nb_events=100000

#    run_name = 'pp_w_May14'
#    events=lhco_events('pp_w_May14/pp_w.lhco')
#    events.open_file()

#    PSname=run_name+'PS'

#    curr_dir=os.getcwd()
#    res_PS=os.path.join(curr_dir, PSname)

#    if not os.path.isdir(res_PS): os.mkdir(res_PS)
    

#    print 'number of events : %s' % counter_event
#    invariants_info.print_boundaries()
#    print "  "
#    print 'feature tree evt 0: '
#    print features_tree[0]
#    print 'feature tree evt 1: '
#    print features_tree[1]

#    for i in range(1):
#      print "classification tree %i" % i
#      classi=box(features_tree)
#      classi.print_boundaries()
#      classi.generate_all_daughters(5)
#      classi.print_tree()
#      tree_text=classi.record_tree()
#      print tree_text
#      classi.get_boxes(0,1,classi.boundaries[0],classi.boundaries[1])

#    stop 

    mincut=100

 
#    print "40 trees:"

    if 1:
      features_tree=get_all_features(sys.argv[1])
      supertree=superbox(all_events=features_tree,nb_trees=100,nb_events=10000,nb_layers=6)
      supertree.get_boxes(feat1=0,feat2=1,tree=0)
     
#      supertree.record_super_tree("SUPERTREES/Supertree_oct21_lin.out")
      supertree.record_super_tree("SUPERTREES/Supertree_test.out")
#    stop

    if 0:
  
      test_supertree=superbox()

#      test_supertree.read_super_tree("SUPERTREES/Supertree_oct21_lin.out")

      test_supertree.read_super_tree("SUPERTREES/Supertree_test.out")

      exp_sample=sys.argv[1]
      features_exp=get_all_features(exp_sample)
      test_supertree.set_deviation(features_exp)
      test_supertree.set_cumul_deviation(features_exp) 
      test_supertree.get_info_tree(0)
      print "start two-point fct"
      corr=test_supertree.get_two_point_fct()
      print "corr: %s" % corr
      corr=test_supertree.get_two_point_fct_long(features_exp)
      print "corr: %s" % corr


#    for i in range(200):
#      x1=-1.0+(i+1)/101.0
#      x2=0.0
#      one_point=[x1, x2]
##      weight1=supertree.get_weight(one_point)
#      weight2=test_supertree.get_weight(one_point)
#      dev=test_supertree.get_deviation(one_point)
##      print "%s  %s %s " % (x1,weight1, weight2)
#      print "%s  %s %s " % (x1,weight2, dev)
   

