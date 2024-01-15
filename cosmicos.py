#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 20:11:42 2023

@author: carol
"""

import numpy as np
import matplotlib.pyplot as plt
from modulos_cosmicos import *
import math
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from uncertainties import ufloat


plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams["figure.figsize"] = [8.0, 5.50]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = "serif"

#plt.savefig('grafica1.png', dpi=300)


files1=['medidas1000V_410.txt','medidas1000V_1000.txt','medidas1000V_581.txt','coincidenciasaccidentais.txt','medida_alex_abcnod.txt']
files2=['v/800V_sentrapo.txt','v/800V_contrapo.txt','v/825V.txt','v/850V.txt','v/875V.txt','v/900V.txt','v/925V.txt','v/950V.txt'] #todas as medidas con trigger en AD e variando a voltaxe de B e C
files3=['xeo/distancia1_300.txt','xeo/distancia1_400.txt','xeo/distancia2_133.txt'] #datos para obter a eficiencia xeometrica
files4=['intr/triggerAD_600.txt','intr/triggerBC_1200_malalineado.txt','intr/triggerBC_benalineado.txt','intr/triggerBD_1200_ABcambiados.txt']
files5=['noites/noite1.txt','noites/noite2.txt','noites/noite3.txt']
files6=['nico/BC_900v.txt','nico/BC_850v.txt','nico/medidasABCD1000V_mercores_500contas.txt','nico/umbral200.txt','nico/umbral500.txt','nico/umbral1000.txt']
files7=['alex/A.txt','alex/ABC.txt','alex/ABCD.txt','alex/abcnod.txt']



m, t = read(files5[0])
m=supresionCeros(m,t, -80)

#------Tasas umbrais---------

'''
Neste apartado obtemos os valores das tasas para 
os distintos umbras con valores medidos e posteriormente 
filtrando os datos co umbral desexado
'''



# u1=-200;u2=-500;u3=-1000

# m_n1,t_n1=read(files6[-3])
# m_n1_no,t1_no=novoUmbral(m_n1,t_n1, u1,0,0,0)
# m_n2,t_n2=read(files6[-2])
# m_n2_no,t2_no=novoUmbral(m_n2,t_n2, u2,0,0,0)
# m_n3,t_n3=read(files6[-1])
# m_n3_no,t3_no=novoUmbral(m_n3,t_n3, u3,0,0,0)

# m,t=read(files1[1])
# m=supresionCeros(m,t, -80)
# m1,t1=novoUmbral(m,t,u1,0,0,0)
# m2,t2=novoUmbral(m,t,u2,0,0,0)
# m3,t3=novoUmbral(m,t,u3,0,0,0)

# tasa1,u_t1=tasa(len(m_n1),t_n1[-1],t_n1[0])
# tasa2,u_t2=tasa(len(m_n2),t_n2[-1],t_n2[0])
# tasa3,u_t3=tasa(len(m_n3),t_n3[-1],t_n3[0])

# tasa1_,u_t1_=tasa(len(m1),t1[-1],t1[0])
# tasa2_,u_t2_=tasa(len(m2),t2[-1],t2[0])
# tasa3_,u_t3_=tasa(len(m3),t3[-1],t3[0])




#-------- Eficiencias e fluxo dos muons incidentes -------

'''
Calculamos os valores da eficiencia intrínseca e os valores das nicertezas
das eficiencias obtidas co montecarlo

correximos os valores da ef_intr de D pero da >1
obtemos os valores expetimentais das ef_xeo para d=30 e d=80 cm
calculamos o fluxo de particulas incidentes para distintos archivos

'''

# umbral=-250

# #trigger AD
# m, t = read(files4[0])
# m=supresionCeros(m,t, -80)
# m,t=novoUmbral(m,t, -150,-150,0,0)
# min_ampl,min_t=minimos_amplitude(m)

# ef_i_1=ef_intr(min_ampl,umbral)
# ef_i_B=ef_i_1[3]
# ef_i_C=ef_i_1[2]

# # N=1000000
# # contas=[931847, 866570, 805275, 349867, 84482]

# # ef_xeo_D = ufloat(contas[2], np.sqrt(contas[2])) / ufloat(N, np.sqrt(N))
# # ef_xeo_C=ufloat(contas[1], np.sqrt(contas[1])) / ufloat(N, np.sqrt(N))
# # ef_xeo_B=ufloat(contas[0], np.sqrt(contas[0])) / ufloat(N, np.sqrt(N))


# #trigger BC
# m2, t2 = read(files4[2])
# m2=supresionCeros(m2,t2, -80)
# m2_,t2_=novoUmbral(m2,t2, -120,-120,0,0)
# min_ampl_2,min_t_2=minimos_amplitude(m2_)

# ef_i_2=ef_intr(min_ampl_2,umbral)
# ef_i_A=ef_i_2[2]
# ef_i_D=ef_i_2[3]/ef_xeo_D

# chd=min_ampl_2[3]
# n_contas=len(chd);ef=0
# for i in range(n_contas):
#     if chd[i]<=umbral:
#         ef+=1

# #print('ef_xeo_D',ef_xeo_D)
# #ef_i_D=(ef_xeo_D*ufloat(ef, np.sqrt(ef))) / ufloat(n_contas, np.sqrt(n_contas))
# #print('ef_i_D',ef_i_D)

# m2,t2=read(files3[1])
# m3,t3=read(files3[2])
# m5,t5=read(files1[1])


# tasa_d1,u1=tasa(len(m2),t2[-1],t2[0])
# tasa_d2,u2=tasa(len(m3), t3[-1], t3[0])
# tasa_3,u3=tasa(len(m5), t5[-1], t5[0])

# t_d1=ufloat(tasa_d1,u1)
# t_d2=ufloat(tasa_d2,u2)
# t_0=ufloat(tasa_3,u3)

# e_a1=ufloat(0.931,0.039)
# e_b1=ufloat(0.992,0.057)
# e_c1=ufloat(0.988,0.057)
# e_d1=ufloat(0.858,0.037)
# e_xeo1=ufloat(0.8053,0.0012)
# e_a2=ufloat(0.929,0.039)
# e_b2=ufloat(0.988,0.057)

# N=1000000
# contas=[931847, 866570, 805275, 349867, 84482]

# ef_xeo_1 = ufloat(contas[-2], np.sqrt(contas[2])) / ufloat(N, np.sqrt(N))
# ef_xeo_2=ufloat(contas[-1], np.sqrt(contas[1])) / ufloat(N, np.sqrt(N))

# m, t = read(files1[1])
# m=supresionCeros(m,t, -80)
# m,t=novoUmbral(m,t, -250,0,0,0)
# m2, t2 = read(files7[1])
# m2=supresionCeros(m2,t2, -80)
# m2,t2=novoUmbral(m2,t2, -250,-250,-250,0)
# m3, t3 = read(files7[2])
# m3=supresionCeros(m3,t3, -80)
# m3,t3=novoUmbral(m3,t3, -250,-250,-250,-250)
# m4, t4 = read(files7[3])
# m4=supresionCeros(m4,t4, -80)
# m4,t4=novoUmbral(m4,t4, -250,-250,-250,None)

# fluxo1=fluxo_muons(m, t)
# fluxo2=fluxo_muons(m2, t2)
# fluxo3=fluxo_muons(m3, t3)
# fluxo4=fluxo_muons(m4, t4)

# time=(t[-1]-t[0])/1000000
# time2=(t2[-1]-t2[0])/1000000
# time3=(t3[-1]-t3[0])/1000000
# time4=(t4[-1]-t4[0])/1000000
# s_t=np.sqrt(2)*10**-6

# triggers=ufloat(len(m), np.sqrt(len(m)))
# triggers2=ufloat(len(m2), np.sqrt(len(m2)))
# triggers3=ufloat(len(m3), np.sqrt(len(m3)))
# triggers4=ufloat(len(m4), np.sqrt(len(m4)))


# m, t = read(files1[1])
# m=supresionCeros(m,t, -80)
# m,t=novoUmbral(m,t, -120,-120,-120,-120)
# m2, t2 = read(files5[0])
# m2=supresionCeros(m2,t2, -80)
# m2,t2=novoUmbral(m2,t2, -150,-150,-150,None)

# fluxo1=fluxo_muons(m, t)
# fluxo2=fluxo_muons(m2, t2)

# time=(t[-1]-t[0])/1000000
# time2=(t2[-1]-t2[0])/1000000

# triggers=ufloat(len(m), np.sqrt(len(m)))
# triggers2=ufloat(len(m2), np.sqrt(len(m2)))




#-----------Representacion de un trigger tipico----------
'''
representamos un dos triggeres que tivesen amplitude de -1V aprox
'''

# m, t = read(files1[2])
# m_nova=supresionCeros(m,t, -80)
# m,t=novoUmbral(m_nova,t, -250,-250,-250,-250)

# ns=m[5][0];cha=m[5][1];chb=m[5][2];chc=m[5][3];chd=m[5][4] #escollemos o sexto trigger, por exemplo
# plt.plot()
# plt.plot(ns[:-1],cha[:-1],label='Canal A')
# plt.plot(ns[:-1],chb[:-1],label='Canal B')
# plt.plot(ns[:-1],chc[:-1],label='Canal C')
# plt.plot(ns[:-1],chd[:-1],label='Canal D')
# plt.xlabel('t (ns)')
# plt.ylabel('Amplitude (mV)')
# plt.grid(alpha=0.5)
# plt.legend(loc='lower right')
# plt.show()

#-----------Graficas amplitudes distintos trigger
'''
representamos as graficas dos minimos das amplitudes 
para distintas condicions de trigger
'''

# files_AD=['v/850V.txt','v/900V.txt','v/950V.txt']
# files_BC=['nico/BC_900v.txt','nico/BC_850v.txt']
# files_AD.reverse()
# labels=['HV=850 V','HV=900 V','HV=950 V']
# labels.reverse()
# hist_amplitudes(files_AD,labels)


#----------Histogramas intervalos de tempo fixos------------

'''
obtemos os histogramas de poisson para intervalos de tempo fixos. 
cambiamos o delta e o numero de bins en cada caso 
para adecuar o noso histograma (buscando o chi2 mais pequeno)
'''

# m, t = read(files1[1])
# m2, t2 = read(files5[1])
# m3, t3 = read(files5[2])
# m=supresionCeros(m,t, -80)
# m2=supresionCeros(m2,t2, -80)
# m3=supresionCeros(m3,t3, -80)
# umbral=-120
# m,t=novoUmbral(m,t, umbral,umbral,umbral,None)
# m2,t2=novoUmbral(m2,t2, umbral,umbral,umbral,None)
# m3,t3=novoUmbral(m3,t3, umbral,umbral,umbral,None)

# delta=15

# Nrepeat=n_contas(t,delta)
# Nrepeat2=n_contas(t2,delta)
# #Nrepeat=n_contas(t3,delta)
# #Nrepeat=np.concatenate((Nrepeat,Nrepeat2), axis=None)

# print('\n')
# print('Arquivo ABCD carol')
# print('Intervalo de tempo', delta)
# print('\n')

# Nmax=int(max(Nrepeat))

# n_bins = math.ceil(math.log2(sum(Nrepeat))) + 1
# n_bins=18
# hist, bordes = np.histogram(Nrepeat, bins=n_bins,range=(0,Nmax))
# media=media_(hist)

# #print('A media do histograma é',media)
# print('O numero de bins é', n_bins)
# chiv,chiv_g=test_chi(hist,media)

# ancho_bin=Nmax/n_bins

# x=bordes+ancho_bin/2
# x=x[:-1]
# #x=np.linspace(ancho_bin/2,bordes[-1]+ancho_bin/2,len(chiv))
# f = interp1d(x, chiv, kind='cubic')
# f2= interp1d(x, chiv_g, kind='cubic')
# x_new = np.linspace(x[0], x[-1], 100)
# y_new = f(x_new)
# y_new_2=f2(x_new)

# plt.figure()
# plt.hist(Nrepeat,bins=n_bins,range=(0,Nmax), alpha=0.99, rwidth=0.85,color='lightskyblue')
# plt.plot(x_new,y_new,'--',color='darkslateblue')
# plt.plot(x_new,y_new_2,'--',color='mediumpurple')
# plt.plot(x,chiv,'.',label='Poisson',color='darkslateblue')
# plt.plot(x,chiv_g,'.',label='Gauss',color='mediumpurple')

# plt.title('dt ='+str(delta)+'s')
# plt.grid(axis='y',alpha=0.75)
# plt.xlabel('#contas')
# plt.ylabel('Frecuencia')
# plt.legend(loc='upper right')
# plt.show()


#---------Histograma diferencias tempo entre disparos--------------
'''
representamos as distribucions tomando os intervalos de tempo entre sucesivos disparos.
aqui axustamos a unha exponencial da forma a*exp(-b*t)+c
se usamos só un parametro no axuste (r*exp(-r*t)) non dá nada... 
a exponencial son todos valores de aproximadamente 0 non sei por que
'''

# m, t = read(files5[0])
# #m2, t2 = read(files5[1])
# #m3, t3 = read(files5[2])
# m=supresionCeros(m,t, -80)
# #m2=supresionCeros(m2,t2, -80)
# #m3=supresionCeros(m3,t3, -80) 
# m,t=novoUmbral(m,t, -250,-250,-250,None)
# #m2,t2=novoUmbral(m2,t2, -150,-150,-150,None)
# #m3,t3=novoUmbral(m3,t3, -150,-150,-150,None)

# dt=np.diff(t)/1e6
# #dt2=np.diff(t2)/1e6
# #dt3=np.diff(t3)/1e6

# #dt=np.concatenate((dt1,dt2,dt3))

# n_bins = math.ceil(math.log2(len(dt))) + 1
# n_bins=46
# lim=300
# bin_range = (0, lim)
# #n_bins=30

#--------Axuste a unha exponencial
# hist, bordes = np.histogram(dt, bins=n_bins,range=bin_range)

# ancho_bins = bordes[1] - bordes[0]
# x = bordes[:-1] + ancho_bins/2
# sigma = 1 / np.sqrt(hist + 1e-6)
# bounds = ([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf])
# popt, pcov = curve_fit(exponencial, x, hist, sigma=sigma)
# popt2, pcov2 = curve_fit(exponencial_sola, x, hist, sigma=sigma)
# axuste = exponencial(x, *popt)
# axuste2 = exponencial_sola(x, *popt2)
# print('\n')
# print('Empregando tres parametros para o axuste:')
# chi2,g_l,chi2_r=chi_axuste(hist,axuste)

# print('\n')
# print('Empregando un so parametro para o axuste:')
# chi2_,g_l_,chi2_r_=chi_axuste(hist,axuste2)

# a,b,c=popt
# perr = np.sqrt(np.diag(pcov))
# a_err, b_err, c_err = perr
# r=popt2
# r_err=np.sqrt(np.diag(pcov2))
# b_float=ufloat(b,b_err)

# print('\n')
# print('O número de bins é',n_bins)
# print('\n')
# print('a=',a,a_err)
# print('b=',b,b_err)
# print('c=',c,c_err)
# print('\n')
# print('r=',r,r_err)


# # Graficar el histograma y la función exponencial ajustada
# x2=np.linspace(x[0],x[-1],100)
# plt.grid(alpha=0.75)
# plt.hist(dt, bins=n_bins,range=bin_range,alpha=0.99, rwidth=0.85,color='lightskyblue')
# plt.plot(x2, exponencial(x2, *popt), '--',color='slateblue', label='$a e^{-bt}+c$')
# #plt.plot(x2, exponencial_sola(x2, *popt2), '--',color='slateblue', label='$r e^{-rt}$')
# plt.xlabel('t(s)')
# plt.ylabel('Frecuencia')

# plt.legend()
# plt.show()

#---------Estudo da vida media do muon
'''
Aqui volvemos a estudar distribucion de intervalos de tempo variables pero
agora con datos dos tempos dos segundos picos

obtemos o axuste a unha exponencial usando distintos tipos de datos... picos en
B,C,D ou picos en C ou picos en D
'''

umbral=150
m, t = read(files5[0])
m2, t2 = read(files5[1])
m3, t3 = read(files5[2])
m=supresionCeros(m,t, -80)
m2=supresionCeros(m2,t2, -80)
m3=supresionCeros(m3,t3, -80) 
m,t=novoUmbral(m,t, -umbral,-umbral,-umbral,None)
m2,t2=novoUmbral(m2,t2, -umbral,-umbral,-umbral,None)
m3,t3=novoUmbral(m3,t3, -umbral,-umbral,-umbral,None)


channel=4
min_t=100
print('\n')
print('canal',channel)
print('umbral',umbral)
print('\n')
ele_t,secPeaks=searchSecondMinimum(m,channel,min_t,-umbral)
print('O número de segundos picos no arquivo 1 é',secPeaks)
ele_t_2,secPeaks2=searchSecondMinimum(m2,channel,min_t,-umbral)
print('O número de segundos picos no arquivo 2 é',secPeaks2)
ele_t_3,secPeaks3=searchSecondMinimum(m3,channel,min_t,-umbral)
print('O número de segundos picos no arquivo 3 é',secPeaks3)


ele_t_tot=[]
ele_t_tot_3=[]
ele_t_tot_4=[]

channels=[2,3,4]

for i in range(len(channels)):
    print('channel',channels[i])
    print('\n')
    ele_t,secPeaks=searchSecondMinimum(m,channels[i],min_t,-umbral)
    print('O número de segundos picos no arquivo 1 é',secPeaks)
    for j in range(len(ele_t)):
        ele_t_tot.append(ele_t[j])
        if i==1: #se o canal é o 3
            ele_t_tot_3.append(ele_t[j])
        if i==2: #se o canal é o 4
            ele_t_tot_4.append(ele_t[j])
    ele_t_2,secPeaks2=searchSecondMinimum(m2,channels[i],min_t,-umbral)
    for k in range(len(ele_t_2)):
        ele_t_tot.append(ele_t_2[k])
        if i==1:
            ele_t_tot_3.append(ele_t_2[k])
        if i==2:
            ele_t_tot_4.append(ele_t_2[k])
    print('O número de segundos picos no arquivo 2 é',secPeaks2)
    ele_t_3,secPeaks3=searchSecondMinimum(m3,channels[i],min_t,-umbral)
    for l in range(len(ele_t_3)):
        ele_t_tot.append(ele_t_3[l])
        if i==1:
            ele_t_tot_3.append(ele_t_3[l])
        if i==2:
            ele_t_tot_4.append(ele_t_3[l])
    print('O número de segundos picos no arquivo 3 é',secPeaks3)

ele_t_tot=[6024, 6028, 7088, 1956, 1452, 1332, 1336, 1720, 2904, 4952, 4956, 356, 360, 
            3152, 1680, 308, 1192, 696, 532, 4044, 4048, 120, 10120, 3128, 3132, 6020, 
            2932, 2936, 4492, 2184, 436, 4952, 932, 936, 264, 1548, 1552, 116, 120, 1480, 
            14260, 300, 304, 268, 11088, 11092, 4768, 3488, 3492, 1852, 220, 104, 956, 3808, 
            3812, 2216, 2220, 5308, 736, 1484, 504, 216, 116,120] #len=64

ele_t_tot_3=[10120, 3128, 3132, 6020, 2932, 2936, 4492, 2184, 436, 4952, 932, 936, 264, 1548, 
              1552, 116, 120, 1480, 14260, 300, 304, 268, 11088, 11092, 4768]

ele_t_tot_4=[3488, 3492, 1852, 220, 104, 956, 3808, 3812, 2216, 2220, 5308, 736, 1484, 504, 216, 116, 120]

n_bins = math.ceil(math.log2(len(ele_t_tot))) + 1
n_bins=7
lim=max(ele_t_tot)
bin_range=(0,lim)

hist, bordes = np.histogram(ele_t_tot, bins=n_bins,range=bin_range)

ancho_bins = bordes[1] - bordes[0]
x = bordes[:-1] + ancho_bins/2
sigma = 1 / np.sqrt(hist)

popt, pcov = curve_fit(exponencial_tau, x, hist, sigma=sigma,bounds=([0,1000],[250,4000]))
axuste = exponencial_tau(x, *popt)

chi2,g_l,chi2_r=chi_axuste(hist,axuste,2)

perr = np.sqrt(np.diag(pcov))
s_N, s_tau = perr

print('\n')
print('n_bins',n_bins)
print('\n')
print('tau=',popt[1])
print('s_tau=',s_tau)

# Graficar el histograma y la función exponencial ajustada
x2=np.linspace(x[0],x[-1],100)
plt.grid(alpha=0.75)
plt.hist(ele_t_tot, bins=n_bins,range=bin_range,alpha=0.7, rwidth=0.85,color='cornflowerblue')
plt.plot(x2, exponencial_tau(x2, *popt), '--',color='rebeccapurple', label=r'$N_0 exp(-t/ \tau)$')
#plt.plot(x2, exponencial_sola(x2, *popt2), '--',color='slateblue', label='$r e^{-rt}$')
plt.xlabel('t(s)')
plt.ylabel('Frecuencia')

plt.legend()
plt.show()

n_bins_1=15
n_bins_2=8
bin_range=(0,14260)

hist4, bordes4 = np.histogram(ele_t_tot_4, bins=n_bins_1,range=bin_range)
hist5, bordes5 = np.histogram(ele_t_tot_4, bins=n_bins_2,range=bin_range)


plt.hist(ele_t_tot_4, bins=n_bins_2,range=(0,14260),alpha=0.9, rwidth=0.85,color='cornflowerblue',label='Detector D')
plt.hist(ele_t_tot_3, bins=n_bins_1,range=bin_range,alpha=0.7, rwidth=0.85,color='lightskyblue',label='Detector B')
plt.grid(axis='y',alpha=0.75)
plt.xlabel('t(s)')
plt.ylabel('Frecuencia')
plt.legend()
plt.show()



m,t=read(files_BC[2])
m1,t1=minimos_amplitude(m)
import numpy as np
import matplotlib.pyplot as plt
covarianza=np.cov(m1[0],m1[1])
correlation_matrix = np.corrcoef(m1[0],m1[3])
print(covarianza, correlation_matrix)
plt.scatter(m1[0], m1[3], color='darkblue',alpha=0.25,marker=".")
plt.xlim(-250,-2150)
plt.xlabel('Amplitude A(mV)')
plt.ylabel('Amplitude D(mV)')
plt.title('Resultados da correlación entre detectores')
plt.ylim(0,-2150)