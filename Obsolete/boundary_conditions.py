#!/usr/bin/python
from numpy import piecewise
top = (('transmissive',(),(0.,1.)))
bottom = (('transmissive',(),(0.,1.)))
left = (('dirichlet',(...),(0.,.5)),('dirichlet',(...),(.5,1.)))
right = (('transmissive',(...),(0.,1.)))



if __name__ == "__main__":


    
top=[]
bottom=[]
left=[]
right=[]
front=[]
back=[]

top.append(('transmissive',(),(0.,1.)))
bottom.append(('transmissive',(),(0.,1.)))
front.append(('transmissive',(),(0.,1.)))
back.append(('transmissive',(),(0.,1.)))
left.append(('dirichlet',(...),(0.,.5)))
left.append(('dirichlet',(...),(.5,1.)))
right.append(('transmissive',(),(0.,1.)))

for i in bound:
    
