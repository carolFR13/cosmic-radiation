#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eficiencias

@author: carol
"""
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from numba import jit 



#DEFINIMOS METODO MC
@jit(nopython=True) #non utiliza o interprete de python para mellorar o rendemento
def MC_A(N, L, d1,d2):
    contas=[0,0,0,0]
    for i in range(N):
        phi=2*np.pi*rd.random() 
        theta=np.arccos((1-rd.random())**(1/4))
        
        x1=rd.uniform(0,L) #punto de impacto na placa 1
        y1=rd.uniform(0,L)
        
        x2=x1-d1*np.tan(theta)*np.cos(phi) 
        y2=y1-d1*np.tan(theta)*np.sin(phi)
        
        x3=x1-2*d1*np.tan(theta)*np.cos(phi)
        y3=y1-2*d1*np.tan(theta)*np.sin(phi)
        
        x4=x1-(2*d1+d2)*np.tan(theta)*np.cos(phi)
        y4=y1-(2*d1+d2)*np.tan(theta)*np.sin(phi)
        
        x=[x1,x2,x3,x4]
        y=[y1,y2,y3,y4]
        
        for i in range(len(x)):
            if 0<x[i]<L and 0<y[i]<L : #se pasa polo detector rexistrase a conta
                contas[i]+=1 
    return contas

def MC_B(N, L, d1,d2):
    contas=[0,0,0,0]
    for i in range(N):
        phi=2*np.pi*rd.random() 
        theta=np.arccos((1-rd.random())**(1/4))
        
        x2=rd.uniform(0,L) #punto de impacto na placa 1
        y2=rd.uniform(0,L)
        
        x1=x2+d1*np.tan(theta)*np.cos(phi) 
        y1=y2+d1*np.tan(theta)*np.sin(phi)
        
        x3=x2-d1*np.tan(theta)*np.cos(phi)
        y3=y2-d1*np.tan(theta)*np.sin(phi)
        
        x4=x2-(d1+d2)*np.tan(theta)*np.cos(phi)
        y4=y2-(d1+d2)*np.tan(theta)*np.sin(phi)
        
        x=[x1,x2,x3,x4]
        y=[y1,y2,y3,y4]
        
        for i in range(len(x)):
            if 0<x[i]<L and 0<y[i]<L : #se pasa polo detector rexistrase a conta
                contas[i]+=1 
    return contas


def MC2_A(N, L, d1,d2,trigger):
    contas_totais=0
    contas=[0,0,0]
    for i in range(N):
        phi=2*np.pi*rd.random() 
        theta=np.arccos((1-rd.random())**(1/4))
        
        x1=rd.uniform(0,L) #punto de impacto na placa 1
        y1=rd.uniform(0,L)
        
        x2=x1-d1*np.tan(theta)*np.cos(phi) 
        y2=y1-d1*np.tan(theta)*np.sin(phi)
        
        x3=x1-2*d1*np.tan(theta)*np.cos(phi)
        y3=y1-2*d1*np.tan(theta)*np.sin(phi)
        
        x4=x1-(2*d1+d2)*np.tan(theta)*np.cos(phi)
        y4=y1-(2*d1+d2)*np.tan(theta)*np.sin(phi)
        
        x=[x1,x2,x3,x4]
        y=[y1,y2,y3,y4]
        
        if trigger==1: #trigger en A
            if 0<x1<L and 0<y1<L:
                contas_totais+=1
                if 0<x2<L and 0<y2<L:
                    contas[0]+=1
                if 0<x3<L and 0<y3<L:
                    contas[1]+=1
                if 0<x4<L and 0<y4<L:
                    contas[2]+=1
        if trigger==2: #trigger AB
            if 0<x1<L and 0<y1<L:
                if 0<x2<L and 0<y2<L:
                    contas_totais+=1
                    if 0<x3<L and 0<y3<L:
                        contas[0]+=1
                    if 0<x4<L and 0<y4<L:
                        contas[1]+=1
        if trigger==3:
            if 0<x1<L and 0<y1<L:
                if 0<x2<L and 0<y2<L:
                    if 0<x3<L and 0<y3<L:
                        contas_totais+=1
                        if 0<x4<L and 0<y4<L:
                            contas[0]+=1     
    return contas,contas_totais



def MC2_B(N, L, d1,d2,trigger):
    contas_totais=0
    contas=[0,0]
    for i in range(N):
        phi=2*np.pi*rd.random() 
        theta=np.arccos((1-rd.random())**(1/4))
        
        x2=rd.uniform(0,L) #punto de impacto na placa 1
        y2=rd.uniform(0,L)
        
        x1=x2+d1*np.tan(theta)*np.cos(phi) 
        y1=y2+d1*np.tan(theta)*np.sin(phi)
        
        x3=x2-d1*np.tan(theta)*np.cos(phi)
        y3=y2-d1*np.tan(theta)*np.sin(phi)
        
        x4=x2-(d1+d2)*np.tan(theta)*np.cos(phi)
        y4=y2-(d1+d2)*np.tan(theta)*np.sin(phi)
        
        x=[x1,x2,x3,x4]
        y=[y1,y2,y3,y4]
        
        if trigger==4:
            if 0<x2<L and 0<y2<L: #impoñemos que pase por B
                if 0<x3<L and 0<y3<L: #impoñemos que pase por C
                    contas_totais+=1
                    if 0<x1<L and 0<y1<L:
                        contas[0]+=1
                    if 0<x4<L and 0<y4<L:
                        contas[1]+=1
        if trigger==5:
            if  0<x2<L and 0<y2<L:
                if 0<x3<L and 0<y3<L:
                    if 0<x4<L and 0<y4<L: #pasa por BCD
                        contas_totais+=1
                        if 0<x1<L and 0<y1<L:
                            contas[0]+=1
    return contas,contas_totais




def MC_ABCnD(N, L, d1,d3):
    contas=[0,0,0,0]
    for i in range(N):
        phi=2*np.pi*rd.random() 
        theta=np.arccos((1-rd.random())**(1/4))
        
        x1=rd.uniform(0,L) #punto de impacto na placa 1
        y1=rd.uniform(0,L)
        
        x2=x1-d1*np.tan(theta)*np.cos(phi) 
        y2=y1-d1*np.tan(theta)*np.sin(phi)
        
        x3=x1-2*d1*np.tan(theta)*np.cos(phi)
        y3=y1-2*d1*np.tan(theta)*np.sin(phi)
        
        x4=x1-(2*d1+d3)*np.tan(theta)*np.cos(phi)
        y4=y1-(2*d1+d3)*np.tan(theta)*np.sin(phi)

        x=[x1,x2,x3,x4]
        y=[y1,y2,y3,y4]

        for i in range(len(x)):
            if 0<x[i]<L and 0<y[i]<L : #se pasa polo detector rexistrase a conta
                contas[i]+=1 
    return contas

N=1000000; L=30.1
d1=2.5 ; d2=2.7; d3=5.4 #(distancia de C a D + distancia laminas plomo)

#1: trigger A ; 2: trigger  AB ; 3: trigger ABC 
#4: trigger BC ; 5: trigger BCD


#contas_A=MC_A(N,L,d1,d2)
#contas_B=MC_B(N,L,d1,d2)

#contas,total=MC2_A(N,L,d1,d2,3)
#contas,total=MC2_B(N,L,d1,d2,5)

contas=MC_ABCnD(N,L,d1,d3)
