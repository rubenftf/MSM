# -*- coding: utf-8 -*-
import numpy as np
from pylab import plt
import fit_class as f
import datetime

start_at=datetime.datetime.now()

Em=1300

data_part=10
time_p=np.loadtxt('time_P.txt', unpack=True)[::data_part]
polarization_data=np.loadtxt('{}P.txt'.format(Em), unpack=True)[::data_part]
time=np.loadtxt('time_S.txt', unpack=True)[::data_part]
strain_data=np.loadtxt('{}S.txt'.format(Em), unpack=True)[::data_part]
data_len=np.shape(strain_data)[0]
time=np.array([time[i] for i in range(data_len)])

p_opt_g=[0.42,2.9e4,3.28e4,3.25e4,0.28,2.0,2.0,1.0,3.9,2.8,0.012]

ittera=1
d_a0=0
d_Ea1=0
d_Ea2=0
d_Ea3=0
d_alpha=0
d_beta=0.5
d_gamma=0
d_s=0
d_x0=0.1
d_eps=0
d_k=-0

p=1
s=1

if p==1:

    fig, a = plt.subplots(figsize=(12.5, 10))
    a.plot(time_p,polarization_data/10**2,'g^')
    for i in range(ittera):

        a0=p_opt_g[0]+d_a0*i
        Ea1=p_opt_g[1]+d_Ea1*i
        Ea2=p_opt_g[2]+d_Ea2*i
        Ea3=p_opt_g[3]+d_Ea3*i
        alpha=p_opt_g[4]+d_alpha*i
        betta=p_opt_g[5]+d_beta*i
        gamma=p_opt_g[6]+d_gamma*i
        s=p_opt_g[7]+d_s*i
        x0=p_opt_g[8]+d_x0*i
        eps=p_opt_g[9]+d_eps*i
        k=p_opt_g[10]+d_k*i
        p_opt=[a0,Ea1,Ea2,Ea3,alpha,betta,gamma,s,x0,eps,k]
        P=f.fit(Em,p_opt)
        polarization_arr=[P.polarization(i) for i in time_p]
        a.plot(time_p,polarization_arr,linewidth=3.0)

    a.set_xscale('log')
    a.set_xlabel(r'$t$, seconds',fontsize=30)
    a.tick_params(axis='x', pad=10)
    a.set_ylabel(r'$ \vartriangle$$P$,$\mu$C/cm$^2$',fontsize=30)
    a.tick_params(axis='y', pad=10)
    a.xaxis.labelpad = 10
    a.yaxis.labelpad = 10
    plt.ylim(-0.1,1.1)
    plt.tick_params(labelsize=30)
    plt.show()
    end_at=datetime.datetime.now()
    print ("time consumed: {}".format(end_at-start_at))

start_at=datetime.datetime.now()
if s==1:

    fig, a = plt.subplots(figsize=(12.5, 10))
    a.plot(time,strain_data,'g^')
    for i in range(ittera):

        a0=p_opt_g[0]+d_a0*i
        Ea1=p_opt_g[1]+d_Ea1*i
        Ea2=p_opt_g[2]+d_Ea2*i
        Ea3=p_opt_g[3]+d_Ea3*i
        alpha=p_opt_g[4]+d_alpha*i
        betta=p_opt_g[5]+d_beta*i
        gamma=p_opt_g[6]+d_gamma*i
        s=p_opt_g[7]+d_s*i
        x0=p_opt_g[8]+d_x0*i
        eps=p_opt_g[9]+d_eps*i
        k=p_opt_g[10]+d_k*i

        p_opt=[a0,Ea1,Ea2,Ea3,alpha,betta,gamma,s,x0,eps,k]
        S=f.fit(Em,p_opt)
        strain_arr=[S.strain(i) for i in time]
        a.plot(time,strain_arr,linewidth=3.0)

    a.set_xscale('log')
    a.set_xlabel(r'$t$, seconds',fontsize=30)
    a.tick_params(axis='x', pad=10)
    a.set_ylabel(r'$ \vartriangle$$S$,$\mu$C/cm$^2$',fontsize=30)
    a.tick_params(axis='y', pad=10)
    a.xaxis.labelpad = 10
    a.yaxis.labelpad = 10
    plt.ylim(-0.4,0.2)
    plt.tick_params(labelsize=30)
    plt.show()

end_at=datetime.datetime.now()

print (end_at-start_at)