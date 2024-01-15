#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:16:19 2023

@author: carol
"""
import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats as stats
from uncertainties import ufloat
import pandas as pd

def read(file):
    #m é a matriz con todos os triggers, cada elemento do array son 5 arrays que se corresponden con ns, cha, chb e chc
    file = open(file, "r") 
    triggers=0
    t=[]
    m=[]
    ns=[];cha=[];chb=[];chc=[];chd=[]; it=[]
    
    for l in file:
        line = l.strip()
        if len(line) == 0: 
            continue #continue if empty line
        if len(line.split()) == 1: #if only one word: time in microseconds
            if triggers==0: 
                time_init = int(line) #definimos o t inicial
            else:
                it=np.stack((ns,cha,chb,chc,chd))
                m.append(it) #cada un dos elementos de m é un array con ns,cha,chb,chc e chd
                ns=[];cha=[];chb=[];chc=[];chd=[]; it=[] 
            triggers=triggers+1 
            time = int(line) 
            t.append(time) 
            
        else:
            nsp, chap, chbp, chcp, chdp = line.split(" ")
            ns.append(int(nsp))
            cha.append(int(chap))
            chb.append(int(chbp))
            chc.append(int(chcp))
            chd.append(int(chdp))
    it=np.stack((ns,cha,chb,chc,chd))
    m.append(it) #incluímos o ultimo elemento
    return m, t

def read_hector(file):
    file = open(file, "r") 
    triggers=0
    m = []
    t = []
    for l in file:
        it = {}
        line = l.strip()
        if len(line) == 0: 
            continue 
        if len(line.split()) == 1:
            if triggers==0: 
                time_init = int(line)
            triggers=triggers+1 
            time = int(line) 
            t.append(time) 
            continue
        else:
            ns, cha, chb, chc, chd = line.split(" ")
            it['trigger'] = triggers 
            it['time_ns'] = ns 
            it['channelA'] = cha 
            it['channelB'] = chb 
            it['channelC'] = chc 
            it['channelD'] = chd 
            m.append(it)
    return m, t, time_init, time, triggers

def print_info(file):
    m, t, time_init, time, triggers=read_hector(file)
    ttime = (time-time_init)/1000000
    print('Arquivo:',file)
    print('Initial time: {}' .format(time_init))
    print('Final time: {}' .format(time))
    print('Total time in seconds: {}' .format(ttime))
    print('Total numer of triggers: {}' .format(triggers))
    print('Total numer of triggers per second: {}' .format(triggers/ttime))
    
def supresionCeros(m,t,supresion): 
    #devolveche a mesma matriz m pero eliminando todas as filas cuxos valores dos voltaxes sexan todos >supresion
    n=len(m)
    m_nova=[]
    index_delete=[]
    with open("temp.txt", "w") as file: #crea un archivo coas coincidencias sen o ruido
        for i in range(n): #para cada trigger
            ns=m[i][0] ; cha=m[i][1] ; chb=m[i][2] ; chc=m[i][3] ; chd=m[i][4]
            index_delete=[]
            for j in range(len(cha)): #suprimimos os ceros, non cambiamos o num de triggers
                if (cha[j]<supresion) or (chb[j]<supresion) or (chc[j]<supresion) or (chd[j]<supresion): 
                    continue
                else:
                    index_delete.append(j) #indeces dos vectores que hai que borrar 
            #print(len(index_delete))
            ns=np.delete(ns,index_delete)
            cha=np.delete(cha,index_delete)
            chb=np.delete(chb,index_delete)
            chc=np.delete(chc,index_delete)
            chd=np.delete(chd,index_delete)
            it=np.stack((ns,cha,chb,chc,chd))
            file.write(str(t[i]))
            file.write('\n')
            np.savetxt(file, np.transpose(it),fmt='%d')
            file.write('\n')
            m_nova.append(it)
    return m_nova

def novoUmbral(m,t,A,B,C,D): #A:umbral en A; B:umbral en B; ...;D:None/valor do umbral
    #devolveche as matrices m e t pero eliminando todos os triggers con valores > umbral
    m_nova=[]
    t_novo=[]
    eliminados=0
    for i in range(len(m)): #para cada trigger
        ns=m[i][0];cha=m[i][1];chb=m[i][2];chc=m[i][3];chd=m[i][4]
        for j in range(len(ns)):
            if ns[j]==40:
                if D==None:
                    if (cha[j]>A) or (chb[j]>B) or (chc[j]>C):
                        eliminados+=1
                        #print('Trigger',i,'eliminado')
                    else:
                        m_nova.append(m[i])
                        t_novo.append(t[i])
                        #print('Trigger', i, 'supera o umbral')
                else:
                    if (cha[j]>A) or (chb[j]>B) or (chc[j]>C) or (chd[j]>D):
                        eliminados+=1
                        #print('Trigger',i,'eliminado')
                    else:
                        m_nova.append(m[i])
                        t_novo.append(t[i])
                        #print('Trigger', i, 'supera o umbral')
    porc_elim=(eliminados/len(m))*100
    print('Os triggers eliminados co novo umbral foron', eliminados,'de',len(m),'triggers,', porc_elim,'%')
    return m_nova, t_novo 

def get_parabola_min(x1, y1, x2, y2, x3, y3):
    # Despejando los coeficientes de la parábola
    a = ((y2 - y1) * (x3 - x1) - (y3 - y1) * (x2 - x1)) / ((x1 - x2) * (x3 ** 2 - x1 ** 2) - (x1 - x3) * (x2 ** 2 - x1 ** 2))
    b = ((y2 - y1) - a * (x2 ** 2 - x1 ** 2)) / (x2 - x1)
    c = y1 - a * x1 ** 2 - b * x1
    # Encontrando las coordenadas del vértice
    xv = -b / (2 * a)
    yv = a * xv ** 2 + b * xv + c
    return xv, yv
            
def minimos_amplitude(m):

    minA=[]; minB=[]; minC=[]; minD=[]
    minA_t=[];minB_t=[];minC_t=[];minD_t=[]

    for i in range(len(m)):
        ns=m[i][0] ; cha=m[i][1] ; chb=m[i][2] ; chc=m[i][3] ; chd=m[i][4]
        #print(i)
        
        for j, channel in enumerate([cha,chb,chc,chd]):
            index_min=np.argmin(channel)
            if index_min==0:
                t_min,ampl_min=ns[index_min],channel[index_min]
                #print(t_min,ampl_min)
            elif index_min==len(channel)-1:
                t_min,ampl_min=ns[index_min],channel[index_min]
                #print(t_min,ampl_min)
            else:
                t_min,ampl_min=get_parabola_min(ns[index_min-1],channel[index_min-1],ns[index_min],channel[index_min],ns[index_min+1],channel[index_min+1])
                #print(t_min,ampl_min)
            if j == 0:
                minA.append(ampl_min)
                minA_t.append(t_min)
            elif j == 1:
                minB.append(ampl_min)
                minB_t.append(t_min)
            elif j == 2:
                minC.append(ampl_min)
                minC_t.append(t_min)
            elif j == 3:
                minD.append(ampl_min)
                minD_t.append(t_min)
        
    minimos_ampl=np.stack((minA,minB,minC,minD))
    minimos_t=np.stack((minA_t,minB_t,minC_t,minD_t))
    return minimos_ampl,minimos_t

def hist_amplitudes(files,labels,colors,alphas):
    m_amp=[]
    for file in files:
        m, t = read(file)
        print('Número de triggers',len(m))
        m_nova=supresionCeros(m,t, -80)
        m_nova_2,t_novo_2=novoUmbral(m_nova,t, -200,-200,0,0)
        min_amp, minimos_t=minimos_amplitude(m_nova_2)
        print('O número de mínimos é',(len(min_amp[0])))
        m_amp.append(min_amp)
        n_canales=len(m_amp[0])
    
    plt.clf()
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,sharex=True,figsize=(10, 8))
    
    ax1.grid(axis='y',alpha=0.75)
    ax2.grid(axis='y',alpha=0.75)
    ax3.grid(axis='y',alpha=0.75)
    ax4.grid(axis='y',alpha=0.75)
    
    for i in range(len(m_amp)):
        chb=m_amp[i][0];chc=m_amp[i][1];chd=m_amp[i][2];cha=m_amp[i][3]
        n_canales=30
        bin_range = (-2000, 0)      
        hist1, bins1, _ = ax1.hist(cha,bins=n_canales,range=bin_range,color=colors[i],alpha=alphas[i],rwidth=1,label=labels[i],edgecolor='black')
        hist2, bins2, _ = ax2.hist(chb,bins=n_canales,range=bin_range,color=colors[i],alpha=alphas[i],rwidth=1,label=labels[i],edgecolor='black')
        hist3, bins3, _ = ax3.hist(chc,bins=n_canales,range=bin_range,color=colors[i],alpha=alphas[i],rwidth=1,label=labels[i],edgecolor='black')
        hist4, bins4, _ = ax4.hist(chd,bins=n_canales,range=bin_range,color=colors[i],alpha=alphas[i],rwidth=1,label=labels[i],edgecolor='black')
        
        bin_centers = (bins1[:-1] + bins1[1:]) / 2
        bin_width = bins1[1] - bins1[0]
        dot_x = bin_centers
        
        errors1 = np.sqrt(hist1)
        errors2 = np.sqrt(hist2)
        errors3 = np.sqrt(hist3)
        errors4 = np.sqrt(hist4)
        
        # ax1.errorbar(bin_centers, hist1, yerr=errors1, fmt='none', color=colors[i], linewidth=1.5)
        # ax2.errorbar(bin_centers, hist2, yerr=errors2, fmt='none', color=colors[i], linewidth=1.5)
        # ax3.errorbar(bin_centers, hist3, yerr=errors3, fmt='none', color=colors[i], linewidth=1.5)
        # ax4.errorbar(bin_centers, hist4, yerr=errors4, fmt='none', color=colors[i], linewidth=1.5)
        
        
        # ax1.plot(dot_x, hist1, marker='.', color=colors[i], linestyle='none')
        # ax2.plot(dot_x, hist2, marker='.', color=colors[i], linestyle='none')
        # ax3.plot(dot_x, hist3, marker='.', color=colors[i], linestyle='none')
        # ax4.plot(dot_x, hist4, marker='.', color=colors[i], linestyle='none')

    #,bins=np.linspace(-2000,0,n_canales+1),width=(2000/n_canales)
    ax1.set_xlabel('Amplitude (mV)')
    ax1.set_ylabel('Frecuencia')
    ax2.set_xlabel('Amplitude (mV)')
    ax2.set_ylabel('Frecuencia')
    ax3.set_xlabel('Amplitude (mV)')
    ax3.set_ylabel('Frecuencia')
    ax4.set_xlabel('Amplitude (mV)')
    ax4.set_ylabel('Frecuencia')

    ax1.legend(loc='best',fontsize='x-small')
    ax2.legend(loc='best',fontsize='x-small')
    ax3.legend(loc='best',fontsize='x-small')
    ax4.legend(loc='best',fontsize='x-small')

    ax1.set_title('Canal A')
    ax2.set_title('Canal B')
    ax3.set_title('Canal C')
    ax4.set_title('Canal D')

    plt.xlim(-2000,0)
    plt.savefig('hist_amplitudes.png', dpi=300)
    plt.show()
   

def minimos_amplitude_mal(m): 
    #devolve o valor minimo do vector sen facer o calculo analitico. no vale
    minA=[]; minB=[]; minC=[]; minD=[]
    for i in range(len(m)):
        ns=m[i][0] ; cha=m[i][1] ; chb=m[i][2] ; chc=m[i][3] ; chd=m[i][4]
        minA.append(np.min(cha))
        minB.append(np.min(chb))
        minC.append(np.min(chc))
        minD.append(np.min(chd))
    minimos = np.array([minA, minB, minC, minD])
    return minimos

def fluxo_muons(m,t):
    L=ufloat(30.1, 0.1)
    A=L**2
    triggers=ufloat(len(m), np.sqrt(len(m)))
    time=(t[-1]-t[0])/1000000 #tempo total en segundos
    tasa=(triggers/(A*time))*60 #tasa en cm^-2 min^-1 
    print('O número de contas é', triggers)
    print('A diferenza de tempos en s {:.8f}'.format(time))
    print('A área é', A)
    print('A tasa é', tasa)
    return tasa #teoricamente este valor deberia ser de 1

def tasa(contas,t,t0):
    ttime=(t-t0)/1000000
    tasa=contas/ttime
    u_contas=(np.sqrt(contas)) #asumindo que temos distribucion de poisson
    u_t=0
    u_tasa=np.sqrt((contas*u_t/ttime)**2+(u_contas/ttime)**2)
    return tasa,u_tasa

def n_contas(t,delta):
    ttime = (t[-1]-t[0])/1000000
    num_int=int(ttime/delta)+1
    Nrepeat=np.zeros(num_int)
    for coin_time in t:
        int_coin1=int((coin_time-t[0])/1e6/delta)
        Nrepeat[int_coin1]+=1
    return Nrepeat


def ef_intr(min_amp,umbral): #min_amp:vector cos mínimos da amplitude
    ef_intr=[]
    ef=0
    for channel in min_amp:
        n_contas=len(channel)
        for i in range(n_contas):
            if channel[i]<=umbral:
                ef+=1
        ef_intr_v= ufloat(ef, np.sqrt(ef)) / ufloat(n_contas, np.sqrt(n_contas))
        ef=0
        ef_intr.append(ef_intr_v)
    return ef_intr

def searchSecondMinimum(m,channel,min_time,ampl):
    ele_t=[]
    time=0 ; secPeaks=0
    n_triggers=len(m)
    for i in range(len(m)):
        ns=m[i][0] 
        for j in range(len(ns)):
            if ns[j]>min_time:
                if m[i][channel][j]<ampl:
                    time=ns[j]
                    ampl=m[i][channel][j]
                    if time>0:
                        ele_t.append(time)
                        secPeaks=secPeaks+1
                    time=0
    return ele_t,secPeaks
            
def searchSecondMinimum_hector(m,n_triggers,channel,min_time,min_ampl): #min_time=100;min_amp=400
    ele_t = [] 
    time=0
    counter=1
    secPeaks=0
    while counter < n_triggers:
        ampl=min_ampl #reset
        for item in m:
            if item['trigger'] == counter:
                if int(item['time_ns']) > min_time : # minimum time
                    if int(item[channel]) < ampl : 
                        time = int(item['time_ns']) 
                        ampl = int(item[channel]) 
            else: 
                if time>0 :
                    ele_t.append(time) #take time at minimum
                secPeaks=secPeaks+1 
                time=0 #reset
                counter=counter+1
    return ele_t,secPeaks #returns times (vector) and number of secondary peaks
        
def Poisson(n,l):
    n = int(n)
    P=((l**n)/math.factorial(n))*np.exp(-l)
    return P

def Gauss(n,l):
    G=(1/(2*np.pi*l)**0.5)*np.exp((-(n-l)**2)/(2*l))
    return G

def mejorartest(x): #funcion para mellorar o test de chi cuadrado agrupando os bins con numero de contas menores a 5
    it=0 ; it2=0
    while x[0]<5:
        it+=1 #num de valores fusionados ao principio do vector
        x[1]+=x[0] ; x=np.delete(x,0)
    while x[len(x)-1]<5:
        it2+=1 #numero de valores fusionados ao final do vector
        x[len(x)-2]+=x[len(x)-1] ; x=np.delete(x,len(x)-1)
    return(x,it,it2)

def media_(x,frecuencia):
    frecuencia=list(frecuencia)
    numerador=0;denominador=0
    for i in range (len(frecuencia)):
        numerador+=(frecuencia[i]*x[i]) ; denominador+=frecuencia[i]
    return numerador/denominador      
    

def tabla_chi(contas,frecuencia,nome,it1,it2): #non pode estar normalizado !! 
    #para estudar a veracidade dunha hipotese estadistica non para o axuste dunha función !
    media=media_(contas,frecuencia)
    k=len(frecuencia) #numero de bins
    n=sum(frecuencia) #numero de contas
    xiOi=[];Eiv_Poiss=[];Eiv_Gauss=[]
    chiv_poiss=[];chiv_gauss=[]
    g_l=k-2
    
    for i in range(len(frecuencia)):
        xiOi.append(contas[i]*frecuencia[i]) 
        Ei_Poiss=n*Poisson(contas[i],media) ; Ei_Gauss=n*Gauss(contas[i],media)
        Eiv_Poiss.append(Ei_Poiss) ; Eiv_Gauss.append(Ei_Gauss)
    
    for i in range(it1):
        contas=np.delete(contas,0)
        frecuencia[1]+=frecuencia[0] ; frecuencia=np.delete(frecuencia,0)
        Eiv_Poiss[1]+=Eiv_Poiss[0] ; Eiv_Poiss=np.delete(Eiv_Poiss,0)
        xiOi[1]+=xiOi[0] ;xiOi=np.delete(xiOi,0)
        Eiv_Gauss[1]+=Eiv_Gauss[0] ; Eiv_Gauss=np.delete(Eiv_Gauss,0)
    
    for i in range(it2):
        contas=np.delete(contas,len(contas)-1)
        xiOi[len(xiOi)-2]+=xiOi[len(xiOi)-1] ;xiOi=np.delete(xiOi,len(xiOi)-1)
        frecuencia[len(frecuencia)-2]+=frecuencia[len(frecuencia)-1] ; frecuencia=np.delete(frecuencia,len(frecuencia)-1)
        Eiv_Poiss[len(Eiv_Poiss)-2]+=Eiv_Poiss[len(Eiv_Poiss)-1] ; Eiv_Poiss=np.delete(Eiv_Poiss,len(Eiv_Poiss)-1)
        Eiv_Gauss[len(Eiv_Gauss)-2]+=Eiv_Gauss[len(Eiv_Gauss)-1] ; Eiv_Gauss=np.delete(Eiv_Gauss,len(Eiv_Gauss)-1)
 
    
    for i in range(len(frecuencia)):
        chi_poiss=((frecuencia[i]-Eiv_Poiss[i])**2)/Eiv_Poiss[i]
        chi_gauss=((frecuencia[i]-Eiv_Gauss[i])**2)/Eiv_Gauss[i]
        chiv_poiss.append(chi_poiss)
        chiv_gauss.append(chi_gauss)
       
    media2=media_(contas,frecuencia)
    
    data = {'contas': contas,
            'Oi': frecuencia,
            'xiOi':xiOi,
            'Ei_Poiss':Eiv_Poiss,
            'Ei_Gauss':Eiv_Gauss,
            'chi_poiss':chiv_poiss,
            'chi_gauss':chiv_gauss
            }
    
    df = pd.DataFrame(data)
    df.to_excel(nome, index=False)  
    chi_poiss=sum(chiv_poiss)
    chi_gauss=sum(chiv_gauss)
    p_value=1-stats.chi2.cdf(chi_poiss, g_l)
    p_value_gauss=1-stats.chi2.cdf(chi_gauss, g_l)
    
    print('\n')
    print('A media do histograma é', media)
    print('\n')
    
    print('A media despois de reagrupar do histograma é', media2)
    print('\n')
        
    print('Valores para a distribución de Poisson')
    print('chi_squared =', chi_poiss)
    print('Graos de liberdade:',g_l)
    print('chi_reducido=',chi_poiss/g_l)
    print('1-alpha = ', p_value)

    print('Valores para a distribución de Gauss')
    print('chi_squared =', chi_gauss)
    print('Graos de liberdade:',g_l)
    print('chi_reducido=',chi_gauss/g_l)
    print('1-alpha = ', p_value_gauss)
    
    return g_l,chi_poiss,chi_gauss,p_value,p_value_gauss,Eiv_Poiss,Eiv_Gauss


def Exp(r,t):
    E=r*np.exp(-r*t)
    return E


def tabla_chi_exp(x,hist,media2,nome): #x: punto medio do bin
    #para estudar a veracidade dunha hipotese estadistica non para o axuste dunha función !
    #media=media_(x,hist)
    media=media2
    k=len(hist) #numero de bins
    n=sum(hist) #numero de contas
    xiOi=[];Eiv_exp=[]
    chiv_exp=[]
    g_l=k-2
    r=1/media
    for i in range(len(x)):
        xiOi.append(x[i]*hist[i]) 
        Ei_exp=n*Exp(r,x[i]) 
        Eiv_exp.append(Ei_exp) 
        chi_exp=((hist[i]-Ei_exp)**2)/Ei_exp
        chiv_exp.append(chi_exp)
    
    data = {'xi': x,
            'Oi': hist,
            'xiOi':xiOi,
            'Ei_exp':Eiv_exp,
            'chi_exp':chiv_exp,
            }
    
    df = pd.DataFrame(data)
    df.to_excel(nome, index=False)  
    chi_exp=sum(chiv_exp)
    xiOi_tot=sum(xiOi)
    p_value=1-stats.chi2.cdf(chi_exp, g_l)
    
    print('A media do histograma é', media)
    print('\n')
        
    print('Valores para a distribución exponencial')
    print('chi_squared =', chi_exp)
    print('Graos de liberdade:',g_l)
    print('chi_reducido=',chi_exp/g_l)
    print('1-alpha = ', p_value)
    
    return g_l,chi_exp,p_value,xiOi_tot,Eiv_exp


def trat(x): #media e incerteza da media dunha serie de datos
    mx = np.mean(x)
    n=len(x)
    sa = np.sqrt((sum((x-mx)**2)/(n-1)))
    smx = sa/np.sqrt(n)
    return mx, smx

def media_ponderada(datos, incertidumbres):
    numerador = np.sum(datos / incertidumbres**2)
    denominador = np.sum(1 / incertidumbres**2)
    media_ponderada = numerador / denominador
    error_media_ponderada = np.sqrt(1 / denominador)
    return media_ponderada, error_media_ponderada

def exponencial(x, a, b, c):
    return a * np.exp(-b * x) + c

def exponencial_sola(x,r):
    return r * np.exp(-r * x)

def exponencial_tau(x,N,tau):
    return N*np.exp(-x/tau)

def chi_axuste(hist,axuste,ligaduras):
    diferencia = hist - axuste
    chi2 = np.sum(diferencia**2 / np.sqrt(hist + 1e-6))
    n_bins = len(hist)
    g_l = n_bins - ligaduras
    chi2_r=chi2/g_l
    print('chi_squared =', chi2)
    print('Graos de liberdade:',g_l)
    print('chi_squared reducido=', chi2_r)
    print('Percentil do text X2 = ', 1-stats.chi2.cdf(chi2, g_l))
    return chi2,g_l,chi2_r

