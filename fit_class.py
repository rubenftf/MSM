# -*- coding: utf-8 -*-
import numpy as np
import scipy.optimize
import scipy.integrate
from pylab import plt

class fit():

    def __init__(self,Em,a):
        self.tau0=0.8e-11
        self.Em=float(Em)
        self.int_r0=0.5e3
        self.int_r1=3e3
        self.Ps=0.38
        self.a=a

    def tau1(self,E):
        return self.tau0*np.exp(self.a[1]/E)
    def tau2(self,E):
        return self.tau0*np.exp(self.a[2]/E)
    def tau3(self,E):
        return self.tau0*np.exp(self.a[3]/E)
    def G(self,E):
        return 1/self.Em*1/np.pi*self.a[10]/((E/self.Em-1)**2+self.a[10]**2)

    def Q(self,t,E2):
        a=self.a
        Em=self.Em
        tau1=self.tau1
        tau2=self.tau2
        return lambda t1: a[4]/tau1(Em)*(t1/tau1(Em))**(a[4]-1)*np.exp(-(t1/tau1(Em))**a[4]-((t-t1)/tau2(E2))**a[5])

    def int_Q(self,t,E2):
        Q=self.Q
        return scipy.integrate.quad(Q(t, E2),0,t,epsabs=0.01)[0]

    def Q_integ(self,t):
        int_r0=self.int_r0
        int_r1=self.int_r1
        return scipy.integrate.quad(lambda E2: self.G(E2)*self.int_Q(t, E2),int_r0,int_r1)[0]

    def dp90z(self,t):
        a=self.a
        return a[0]*self.Ps*(2*(1-np.exp(-(t/self.tau1(self.Em))**a[4]))-self.Q_integ(t))

    def dp180z(self,t):
        int_r0=self.int_r0
        int_r1=self.int_r1
        a=self.a
        G=self.G
        dp180z_0=lambda E: 2*self.Ps*(1-a[0])*(1-np.exp(-(t/self.tau3(E))**a[6]))
        return scipy.integrate.quad(lambda E: G(E)*dp180z_0(E),int_r0,int_r1)[0]

    def polarization(self,t):
        return self.dp90z(t)+self.dp180z(t)

    def strain(self,t):
        eps=self.a[9]*1e3
        qzz=self.a[8]
        a=self.a
        return -a[7]*(a[0]*self.Q_integ(t))+2*8.8*1e-12*eps*qzz*self.Em*1e3*(self.dp90z(t)+self.dp180z(t)-self.Ps)

    def plot_G(self):
        int_r0=self.int_r0
        int_r1=self.int_r1
        x=[i*10 for i in range (1,250)]
        y=[self.G(i) for i in x]
        plt.plot(x,y)
        plt.show()
        print (scipy.integrate.quad(lambda E: self.G(E),int_r0,int_r1)[0])