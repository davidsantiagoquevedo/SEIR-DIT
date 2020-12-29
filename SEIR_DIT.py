#!/usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
from scipy.integrate import odeint
import sympy as sp
import math as mt

########################################################################
################ ########## SEIR_DIT CLASS ########## ##################
########################################################################
class SEIR_DIT():
	def __init__(self):
		self.dummy=True

	def set_params(self,beta,theta,phi,phi_,
				  beta_0=1.3743,beta_1=0.019899665,beta_H=0.01,r=0.65,
				  omega=1/4.6, gamma_M=1/2.1, sigma_C=1/3., sigma_CSI=1/4.1, sigma_CT=1/7.1,
				  gamma_HR=1/9., nu=1/14.8, gamma_R=1/3., sigma_HD=1/9., sigma_UD=1/11.1,
				  delta_M=0.9790248, delta_HR=0.6795883, delta_HD=0.1068193, delta_UR=0.1388350,
				  xi_PCR=0.85, xi_AG=0.75, T_PCR=2., T_AG=1., psi=3.,
				  alpha=1/1., t0_DIT=184, n=15.2, nt=0, rho=1/12., q=0.5):
					  
		self.beta_0=beta_0			#Transmission rate outside household
		self.beta_1=beta_1			#Transmission rate within household
		self.beta_H=beta_H			#Transmission rate within hospitals
		self.r=r					#Rate of self-isolation for case infections
		self.beta=beta				#Transmission rate
		
		self.omega=omega			#Mean incubation period
		self.gamma_M=gamma_M		#Mean duration of mild infection
		self.sigma_C=sigma_C		#Mean time upon self-isolation for severe infections
		self.sigma_CSI=sigma_CSI	#Mean duration of isolation for severe infections prior hospitalization
		self.sigma_CT=sigma_CT		#Mean duration of isolation for traced contacts with severe infection prior hospitalization
		self.gamma_HR=gamma_HR		#Mean duration of hospitalization for non-critical cases if survive
		self.nu=nu					#Mean duration in ICU if survive
		self.gamma_R=gamma_R		#Mean duration of stepdown post ICU
		self.sigma_HD=sigma_HD		#Mean duration of hospitalization for non-critical cases if die
		self.sigma_UD=sigma_UD		#Mean duration in ICU if die

		self.delta_M=delta_M		#Probability of mild infections 
		self.delta_HR=delta_HR		#Probability of recovery for hospitalized infections requiring a general hospital bed
		self.delta_HD=delta_HD		#Probability of dying for hospitalized infections requiring a general hospital bed 
		self.delta_UR=delta_UR		#Probability of recovery for hospitalized infections requiring an ICU bed
		self.delta_UD=(1.			#Probability of dying for hospitalized infections requiring an ICU bed 
					 -delta_HR-delta_HD-delta_UR)

		self.xi_PCR=xi_PCR	#Sensitivity of RT-PCR test	
		self.xi_AG=xi_AG	#Sensitivity of Antigen test	
		self.T_PCR=T_PCR	#Mean time to results of RT-PCR test
		self.T_AG=T_AG		#Mean time to results of Antigen test
		self.psi=psi		#Proportion of non-infected suspected cases

		self.alpha=alpha	#Mean time between onset of symptoms and detection
		self.n=float(n)		#Averge number of contacts
		self.nt=float(nt)	#Traced contacts of suspected index cases
		self.rho=rho		#Mean duration of isolation period
		self.q=q			#Rate of accomplishment of isolation DIT strategy
		self.t0_DIT=t0_DIT	#Initial time of DIT strategy

		self.theta=theta	#Detection rate
		self.phi=phi		#Tracing rate of exposed contacts
		self.phi_=phi_		#Tracing rate of non-exposed contacts

	def set_initial(self, N0, E0=0, ET0=0, IM0=0, IMD0=0, IMT0=0, 
					IC0=0, ICT0=0, ICSI0=0, IHR0=0, IUR0=0, IHD0=0, IUD0=0, 
					IR0=0, R0=0, D0=0,
					QIMD10=0, QET0=0, QIMT0=0, QIMT10=0, QICT0=0, QS10=0, QS20=0):

		self.N0=N0
		self.E0=E0
		self.ET0=ET0

		self.IM0=IM0
		self.IMD0=IMD0
		self.IMT0=IMT0

		self.IC0=IC0
		self.ICT0=ICT0
		self.ICSI0=ICSI0

		self.IHR0=IHR0
		self.IUR0=IUR0
		self.IHD0=IHD0
		self.IUD0=IUD0
		self.IR0=IR0
		self.R0=R0
		self.D0=D0

		self.QIMD10=QIMD10
		self.QET0=QET0
		self.QIMT0=QIMT0
		self.QIMT10=QIMT10
		self.QICT0=QICT0
		self.QS10=QS10
		self.QS20=QS20

		self.S0= (self.N0 - self.E0 - self.ET0 - self.IM0 - self.IMD0  - self.IMT0 
				  - self.IC0 - self.ICT0 - self.ICSI0 - self.IHR0 - self.IUR0 - self.IHD0 
				  - self.IUD0 - self.IR0 - self.R0 - self.D0 
				  - self. QIMD10 - self.QET0 - self. QIMT0 - self.QIMT10 - self. QICT0 - self.QS10 - self.QS20)

	def ODES(self,y,t):
		S, E, ET, IM, IMD, IMT, IC, ICT, ICSI, IHR, IUR, IHD, IUD, IR, R, D, N, QIMD1, QET, QIMT, QIMT1, QICT, QS1, QS2 = y
		beta=self.beta(t)
		theta=self.theta(t,IM)
		b=beta/self.n
		
		dSdt = (-beta*S/float(N)*(IM+IC+(1.-self.q)*(IMD+IMT+ICT)+(1.-self.r)*ICSI)
				-self.beta_H*S/float(N)*(IHR+IUR+IHD+IUD+IR)
				-(self.nt+1.)*(self.psi-1.)*theta*self.alpha*S/float(N)*IM
				-self.n*(1.-b)*self.phi_(t)*theta*S/float(N)*IM
				+ 1./self.T_AG*QS1
				+1./(1./self.omega + self.T_PCR)*QS2
				) 
				

		dEdt =(beta*S/float(N)*((1.-self.phi(t)*theta)*IM+IC
							  +(1.-self.r)*ICSI
							  +(1-self.q)*(IMD+IMT+ICT))
			  +self.beta_H*S/float(N)*(IHR+IUR+IHD+IUD+IR)
			  -self.omega*E
			  )
			   
		dETdt = beta*self.phi(t)*theta*IM*S/float(N) - self.omega*ET

		dIMdt = self.delta_M*self.omega*E - theta*self.alpha*IM - (1.-theta)*self.gamma_M*IM
		
		dIMDdt = theta*self.alpha*IM - 1./(1./self.gamma_M - 1./self.alpha)*IMD
		
		dIMTdt = self.delta_M*self.omega*ET  - self.gamma_M*IMT

		dICdt = (1.-self.delta_M)*self.omega*E - self.sigma_C*IC
		
		dICSIdt = self.sigma_C*IC - self.sigma_CSI*ICSI
		
		dICTdt = (1.-self.delta_M)*self.omega*ET - self.sigma_CT*ICT

		dIHRdt = self.delta_HR*(self.sigma_CSI*ICSI+self.sigma_CT*ICT) - self.gamma_HR*IHR
		
		dIURdt = self.delta_UR*(self.sigma_CSI*ICSI+self.sigma_CT*ICT) - self.nu*IUR
		
		dIHDdt = self.delta_HD*(self.sigma_CSI*ICSI+self.sigma_CT*ICT) - self.sigma_HD*IHD
		
		dIUDdt = self.delta_UD*(self.sigma_CSI*ICSI+self.sigma_CT*ICT) - self.sigma_UD*IUD

		dIRdt = self.nu*IUR - self.gamma_R*IR
		
		dRdt = (self.gamma_R*IR
				+self.gamma_HR*IHR
				+(1.-theta)*self.gamma_M*IM
				+1./(1./self.gamma_M-1./self.alpha)*IMD 
				+self.gamma_M*IMT
				)

		dDdt = self.sigma_HD*IHD + self.sigma_UD*IUD
		
		dNdt = -self.sigma_HD*IHD - self.sigma_UD*IUD
		
		#Isolation of index cases and contacts
		dQIMD1dt = self.xi_AG*(1./self.T_AG)*IMD - self.rho*QIMD1
		
		dQETdt = beta*self.phi(t)*theta*IM*S/float(N) - self.omega*QET
		
		dQIMTdt = self.delta_M*self.omega*QET - 1./float(self.T_PCR)*QIMT
		
		dQIMT1dt = (self.xi_PCR*1./float(self.T_PCR)*QIMT 
					- 1./((1./self.rho)
					-(1./self.omega)
					-self.T_PCR)*QIMT1)
		
		dQICTdt = (1.-self.delta_M)*self.omega*QET - self.sigma_CT*QICT
		
		dQS1dt= ((self.psi-1.)*(self.nt+1.)*self.alpha*theta*IM*S/float(N) 
				- 1./self.T_AG*QS1)
		
		dQS2dt = (self.n*(1.-b)*self.phi_(t)*theta*IM*S/float(N) 
				- 1./(1./self.omega + self.T_PCR)*QS2) 
						
		return [dSdt, dEdt, dETdt, dIMdt, dIMDdt, dIMTdt, dICdt, dICTdt, dICSIdt, 
				dIHRdt, dIURdt, dIHDdt, dIUDdt, dIRdt, dRdt, dDdt, dNdt,
				dQIMD1dt ,dQETdt, dQIMTdt, dQIMT1dt, dQICTdt, dQS1dt, dQS2dt]

	def solve(self,t0,tf,dt):
		self.t0=t0
		self.tf=tf
		self.dt_=1/dt
		y0= [self.S0, self.E0, self.ET0, self.IM0, self.IMD0, self.IMT0, 
			 self.IC0, self.ICT0, self.ICSI0, self.IHR0, self.IUR0, 
			 self.IHD0, self.IUD0, self.IR0, self.R0, self.D0, self.N0,
			 self.QIMD10, self.QET0, self.QIMT0, self.QIMT10, self.QICT0, self.QS10, self.QS20]

		t= np.linspace(self.t0, self.tf, (self.tf-self.t0)*self.dt_+1)
		self.t_=t
		solution= odeint(self.ODES,y0,t)

		self.S=solution.T[0]
		self.E=solution.T[1]
		self.ET=solution.T[2]

		self.IM=solution.T[3]
		self.IMD=solution.T[4]
		self.IMT=solution.T[5]

		self.IC=solution.T[6]
		self.ICT=solution.T[7]
		self.ICSI=solution.T[8]

		self.IHR=solution.T[9]
		self.IUR=solution.T[10]
		self.IHD=solution.T[11]
		self.IUD=solution.T[12]

		self.IR=solution.T[13]
		self.R=solution.T[14]
		self.D=solution.T[15]
		self.N=solution.T[16]

		self.QIMD1=solution.T[17]
		self.QET=solution.T[18]
		self.QIMT=solution.T[19]
		self.QIMT1=solution.T[20]
		self.QICT=solution.T[21]
		self.QS1=solution.T[22]
		self.QS2=solution.T[23]

	def count_tracing_isolation(self):
		self.index_cases=[]
		self.isolated_total=[]
		self.isolated_contacts=[]
		self.traced_contacts=[]
		for i in range(len(self.t_)):
			theta=self.theta(self.t_[i],self.IM[i])	
			beta=self.beta(self.t_[i])
			b=beta/self.n
			
			self.index_cases.append(theta*self.IM[i]*self.psi)
			self.isolated_total.append((self.QIMD1[i]+self.QET[i]+
								  self.QIMT[i]+self.QIMT1[i]+
								  self.QICT[i]+
								  self.QS1[i]+self.QS2[i])*self.q)
			self.isolated_contacts.append((self.n*(1-b)*self.phi_(self.t_[i])*theta
										  +(self.psi-1.)*(self.nt)*self.alpha*theta 
										  + beta*self.phi(self.t_[i])*theta)*self.S[i]*self.IM[i]/float(self.N[i])*self.q)
			self.traced_contacts.append((self.n*(1-b)*self.phi_(self.t_[i])*theta
										  +(self.psi-1.)*(self.nt)*self.alpha*theta 
										  + beta*self.phi(self.t_[i])*theta)*self.S[i]*self.IM[i]/float(self.N[i]))
	def count_tests(self):
		self.N_PCR=[]
		self.N_AG=[]
		self.pos_PCR=[]
		self.pos_AG=[]
		for i in range(len(self.t_)):
			theta=self.theta(self.t_[i],self.IM[i])		
			#Before DIT, the level of testing is unknown
			if self.t_[i]<self.t0_DIT:
				self.N_PCR.append(0)
				self.pos_PCR.append(0)
				self.N_AG.append(0)
				self.pos_AG.append(0)
			#Testing starts with DIT strategy	
			if self.t_[i]>=self.t0_DIT:
				beta=self.beta(self.t_[i])
				b=beta/self.n
				
				self.N_PCR.append(self.omega*self.QET[i]
								 +self.n*(1-b)*self.phi_(self.t_[i])*theta*self.S[i]*self.IM[i]/float(self.N[i]))
				self.pos_PCR.append(self.omega*self.QET[i]*self.xi_PCR)
				self.N_AG.append(self.sigma_CSI*self.ICSI[i]*self.psi+theta*self.alpha*self.IM[i]*self.psi)
				self.pos_AG.append((self.sigma_CSI*self.ICSI[i]+theta*self.alpha*self.IM[i])*self.xi_AG)
		self.R_pos_PCR=np.array(self.pos_PCR)/np.array(self.N_PCR)
		self.R_pos_AG=np.array(self.pos_AG)/np.array(self.N_AG)

	def count_icu(self):
		self.icu_occupancy=np.array(self.IUR)+np.array(self.IUD)
		
	def count_daily_deaths(self):
		self.daily_deaths=[]
		for i in range(len(self.t_)):
			self.daily_deaths.append(self.sigma_HD*self.IHD[i] + self.sigma_UD*self.IUD[i])
	
	def calculate_incidence(self):
		self.incidence=[]
		for i in range(len(self.t_)):
			beta=self.beta(self.t_[i])
			self.incidence.append(
				beta*self.S[i]/float(self.N[i])*(self.IM[i]+self.IC[i]
								+(1.-self.q)*(self.IMD[i]+self.IMT[i]+self.ICT[i])
								+(1.-self.r)*self.ICSI[i])
				+self.beta_H*self.S[i]/float(self.N[i])*(self.IHR[i]+self.IUR[i]+self.IHD[i]+self.IUD[i]+self.IR[i]))
	
	def calculate_prevalence(self):
		self.prevalence= (np.array(self.IM)
						 +np.array(self.IMD)
						 +np.array(self.IMT)
						 +np.array(self.IC)
						 +np.array(self.ICT)
						 +np.array(self.ICSI)
						 +np.array(self.IHR)
						 +np.array(self.IUR)
						 +np.array(self.IHD)
						 +np.array(self.IUD)
						 +np.array(self.IR))
						 
	def calculate_attack_rate(self):
		self.attack_rate = np.array(self.R)/self.N0
		
	def calculate_rt(self):
		omega=self.omega
		gamma_M=self.gamma_M
		sigma_C=self.sigma_C
		sigma_CSI=self.sigma_CSI
		sigma_CT=self.sigma_CT
		gamma_HR=self.gamma_HR
		nu=self.nu
		gamma_R=self.gamma_R
		sigma_HD=self.sigma_HD
		sigma_UD=self.sigma_UD

		delta_M = self.delta_M
		delta_HR = self.delta_HR
		delta_UR = self.delta_UR
		delta_HD = self.delta_HD
		delta_UD = self.delta_UD
		
		alpha=self.alpha
		
		beta_H= self.beta_H
		q=self.q
		r=self.r
		self.rt=[]
		for i in range(len(self.t_)):
			beta=self.beta(self.t_[i])
			theta=self.theta(self.t_[i],self.IM[i])	
			phi=self.phi(self.t_[i])	
								
			K11= (beta*((1.-phi*theta)*delta_M/(alpha*theta+gamma_M*(1.-theta))
						+(1.-q)*alpha*delta_M*theta*(1./gamma_M-1./alpha)/(alpha*theta+gamma_M*(1.-theta))
						+(1.-delta_M)/sigma_C
						+(1.-r)*(1.-delta_M)/sigma_CSI)
				  +beta_H*(1-delta_M)*(delta_HR/gamma_HR
						  +delta_UR/nu
						  +delta_HD/sigma_HD
						  +delta_UD/sigma_UD
						  +delta_UR/gamma_R))			
			K12= (beta*((1.-q)*delta_M/gamma_M
						+(1.-q)*(1.-delta_M)/sigma_CT)
				  +beta_H*(1-delta_M)*(delta_HR/gamma_HR
						  +delta_UR/nu
						  +delta_HD/sigma_HD
						  +delta_UD/sigma_UD
						  +delta_UR/gamma_R))
			K21= beta*phi*theta*delta_M/(alpha*theta+gamma_M*(1.-theta))
			lmbd1=(K11 + mt.sqrt(K11**2 + 4*K12*K21))/2.
			reff =lmbd1*self.S[i]/float(self.N[i])
			self.rt.append(reff)
		
	def calculate_beta(self,Rt,S,N):
		omega=self.omega
		gamma_M=self.gamma_M
		sigma_C=self.sigma_C
		sigma_CSI=self.sigma_CSI
		sigma_CT=self.sigma_CT
		gamma_HR=self.gamma_HR
		nu=self.nu
		gamma_R=self.gamma_R
		sigma_HD=self.sigma_HD
		sigma_UD=self.sigma_UD

		delta_M = self.delta_M
		delta_HR = self.delta_HR
		delta_UR = self.delta_UR
		delta_HD = self.delta_HD
		delta_UD = self.delta_UD
		
		alpha=self.alpha
		theta=self.theta
		phi=self.phi
		
		beta_H= self.beta_H
		beta_ = sp.symbols('beta_')
		q=self.q
		r=self.r

		self.beta_solution=[]
		
		#Dominant eigen-value
		K11_= (beta_*((1.-phi*theta)*delta_M/(alpha*theta+gamma_M*(1.-theta))
					+(1.-q)*alpha*delta_M*theta*(1./gamma_M-1./alpha)/(alpha*theta+gamma_M*(1.-theta))
					+(1.-delta_M)/sigma_C
					+(1.-r)*(1.-delta_M)/sigma_CSI)
			 +beta_H*(1-delta_M)*(delta_HR/gamma_HR
					  +delta_UR/nu
					  +delta_HD/sigma_HD
					  +delta_UD/sigma_UD
					  +delta_UR/gamma_R))			
		K12_= (beta_*((1.-q)*delta_M/gamma_M
					+(1.-q)*(1.-delta_M)/sigma_CT)
			  +beta_H*(1-delta_M)*(delta_HR/gamma_HR
					  +delta_UR/nu
					  +delta_HD/sigma_HD
					  +delta_UD/sigma_UD
					  +delta_UR/gamma_R))
		K21_= beta_*phi*theta*delta_M/(alpha*theta+gamma_M*(1.-theta))
		
		lmbd1_=(K11_ + sp.sqrt(K11_**2 + 4*K12_*K21_))/2.
		
		if len(Rt)!=len(S) or len(Rt)!=len(N) or len(S)!=len(N):
			print 'Todos los vectores deben tener el mismo largo'
		else:
			for i in range(len(Rt)):
				Rt_beta=lmbd1_*S[i]/N[i]
				expr= Rt_beta - Rt[i]
				sol= sp.solve(expr)
				self.beta_solution.append(max(sol))
		del(beta_)
			
					
				
