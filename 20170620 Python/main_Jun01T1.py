# main_Jun01T1.py

# Python program for Task 1 of Jun 01,2017
# June 1,2017
# May 22-24, 2017
# 20170411 PM
# Mar18,2017
# Task #20170221
# Related to May16/Oct17,2016
# edited in Feb21,2017
import numpy as np
import numpy.linalg as LA
import scipy.special as special
import scipy.integrate as integrate
from scipy.optimize import fsolve
import os

#
class Epsilon_u:
	"""class of epsilon_u"""

	def __init__(self, typen):
		# self.type is object variable
		self.typen = typen

	def Assign(self):
		if self.typen == 'I':
			return np.sqrt(3)
		elif self.typen == 'II':
			return np.sqrt(2)/(np.sqrt(np.pi-2))
		elif self.typen == 'III':
			return 1
		elif self.typen == 'IV':
			return (3-np.sqrt(3))/2
		elif self.typen == 'O':
			return (np.pi/4)/(np.sin(np.pi/4))

#
class CaseCode:
	"""docstring for ClassName"""
	def __init__(self, case):
		self.case = case

	def VT(self):
		# 1
		if self.case == 'pstar':
			return [0.01*x for x in range(18,251,1)]
		# 2
		elif self.case == 'b':
			return [0.01*x for x in range(1,101,1)]
		# 3
		elif self.case == 'phi':
			return [0.01*x for x in range(30,96,1)]
		# 4
		elif self.case == 'sigma':
			return [0.01*x for x in range(1,121,1)]
		# 5
		elif self.case == 'beta':
			return [0.01*x for x in range(40,100,1)]
		# 6
		elif self.case == 'lambda':
			return [0.01*x for x in range(2,64,1)]
		# 7
		elif self.case == 'c_f':
			return [0.01*x for x in range(12,131,1)]
		# 8
		elif self.case == 'r':
			return [0.001*x for x in range(3,151,1)]
		# 9
		elif self.case == 'c_p':
			return [0.001*x for x in range(10,201,1)]
		# 10
		elif self.case == 'delta':
			return [0.01*x for x in range(21,71,1)]

	def Str_Var(self):
		"""LaTex string"""
		# 1
		if self.case == 'pstar':
			return r'p^*'
		# 2
		elif self.case == 'b':
			return r'b'
		# 3
		elif self.case == 'phi':
			return r'\phi'
		# 4
		elif self.case == 'sigma':
			return r'\sigma'
		# 5
		elif self.case == 'beta':
			return r'\beta'
		# 6
		elif self.case == 'lambda':
			return r'\lambda'
		# 7
		elif self.case == 'c_f':
			return r'c^F'
		# 8
		elif self.case == 'r':
			return r'r'
		# 9
		elif self.case == 'c_p':
			return r'c^P'
		# 10
		elif self.case == 'delta':
			return r'\delta'

	def Tell(self):
		print('Case is: {}'.format(self.case))
		print('The array of variable is:{}'.format(CaseCode.VT(self)))
		print('String of variable is:{}'.format(CaseCode.Str_Var(self)))

#
class Operator(CaseCode):
	"""Main operator to get the solution

	epsilon_d,epsilon_c,theta,alpha"""

	# para is the dictionary updated
	def __init__(self, case, typen, para):
		CaseCode.__init__(self, case)
		self.typen = typen
		self.para = para

	def fun_q_theta(self, x):
		return self.para['A']* x **(-self.para['B1']/self.para['B2'])

	def fun_F_x(self, x):
		if self.typen == 'O':
			temp1 = (np.pi/4.0)/np.sin(np.pi/4.0) - x
			temp2 = (np.pi/2.0)/np.sin(np.pi/2.0) - ((np.pi/4.0)/np.sin(np.pi/4.0))**2.0
			temp3 = 1.0 + (temp1/np.sqrt(temp2))**(-4.0)
			return 1.0-(temp3)**(-1.0)
		elif self.typen == 'I':
			if x < -1.0*self.para['epsilon_u']:
				return 0
			elif x > self.para['epsilon_u']:
				return 1.0
			else:
				return (x + self.para['epsilon_u'])/(2.0*self.para['epsilon_u'])

		elif self.typen == 'II':
			temp=np.sqrt(1.0/np.pi)-x*np.sqrt(.5-1.0/np.pi)
			return 1.0-special.erf(temp)

		elif self.typen == 'III':
			return (.5-.5 * special.erf((np.log(-x+1.0)+ \
				.5 * np.log(2.0))/np.sqrt(2.0 * np.log(2.0))))

		elif self.typen == 'IV':
			return ((3.0-2.0*x)/np.sqrt(3.0))**(-3.0)

	def integrand(self,x):
		return (1-self.fun_F_x(x))

	def fun_int_F(self,a,b):
		temp = integrate.quad(self.integrand,a,b)
		return temp[0]

	def Equation(self,x):
		epsilon_d = x[0]
		epsilon_c = x[1]
		theta = x[2]
		alpha = x[3]

		q_theta = self.fun_q_theta(theta)

		int_Fdu = self.fun_int_F(epsilon_d,self.para['epsilon_u'])
		int_Fcu = self.fun_int_F(epsilon_c,self.para['epsilon_u'])
		####
		part11 = epsilon_d + self.para['lambda'] * \
			int_Fdu / (self.para['r'] + self.para['lambda'] + theta * q_theta)
		part12 = (self.para['b'] - self.para['pstar'])/self.para['sigma']
		####
		temp211 = self.para['delta'] * (self.para['r'] + self.para['lambda']\
			-self.para['phi'] * (1.0-self.fun_F_x(epsilon_c)))\
				 * (epsilon_c-epsilon_d)
		temp212 = self.para['r'] + self.para['lambda'] + theta * q_theta
		temp213 = (self.para['lambda']/(self.para['r']+self.para['lambda'])\
			-self.para['lambda'] * self.para['delta']/temp212) * int_Fcu
		part21 = epsilon_c-temp211/((1.0-self.para['phi']) * temp212)+temp213

		temp2211 = ((self.para['beta']+self.para['phi'] * \
			(1.0-self.para['beta'])) * self.para['c_f'] * (1.0-alpha))/((1.0-self.para['beta']) * \
			(1.0-self.para['phi']))
		temp2212 = (self.para['beta'] * self.para['c_p'] * alpha)/(1.0-self.para['beta'])
		temp221 = theta * (temp2211 + temp2212)
		part22 = self.para['delta'] * epsilon_d+((1.0-self.para['delta']) * \
			(self.para['b']-self.para['pstar'])+temp221)/self.para['delta']
		####
		part31 = 1.0/q_theta

		temp321 = (1.0-self.para['beta']) * (1.0-self.para['phi'])/self.para['c_f']
		temp322 = self.para['sigma'] * (self.para['epsilon_u']-epsilon_c)/(\
			self.para['r']+self.para['lambda'])+\
			self.para['delta'] * self.para['sigma'] * (epsilon_c-epsilon_d)/(\
				(1.0-self.para['phi']) * (self.para['r']+self.para['lambda']+theta * q_theta))
		part32 = temp321 * temp322
		####
		part41 = 1.0/q_theta

		temp421 = (1.0-self.para['beta']) * self.para['delta'] * self.para['sigma'] * (\
			self.para['epsilon_u']-epsilon_d)
		temp422 = self.para['c_p'] * (self.para['r']+self.para['lambda']+theta * q_theta)
		part42 = temp421/temp422

		return [
			part11-part12,
			part21-part22,
			part31-part32,
			part41-part42
		]
	def Solution(self):
		# Tuple
		# Initial value guess
		guess = (-5.0,-1.0,1.0,0.5)

		# fsolve is a wrapper around MINPACK's hybrd and hybrj algorithms
		# [epsilon_d,epsilon_c,theta,alpha]
		# xtol : The calculation will terminate if the relative error between
		# two consecutive iterates is at most xtol.
		# maxfev : The maximum number of calls to the function.
		sol, infodict, ier, mesg = fsolve(self.Equation,guess,full_output=True)


		if ier != 1:
			print("ier = {}".format(ier))
			print("sol = {}".format(sol))
			print("infodict = {}".format(infodict))
			return [sol,ier]
		else:
		 	return [sol,ier]


	def Tell(self):
		CaseCode.Tell(self)
		print('Type of F(x) is:{}'.format(self.typen))
		print('Parameters are:{}'.format(self.para))
		print('Solution is:{}'.format(self.Solution))
		print('Error is:{}'.format(self.Error))

# main part
taskcode = 'Jun01T1'
str_para = 'paraJun01V27'
# Example on Windows:
# source = ['"C:\\My Documents"','C:\\Code']
# for Mac OS or Linux
path_data = '/Users/huwei/Dropbox/On_Local/Python/20170620 Python/'
if not os.path.exists(path_data):
	print('path_data does not exist!')
# casecode: a list
cell_typen = ['I','II','III','IV','O']
cell_case = ['pstar','b','phi','sigma','beta','lambda',\
'c_f','r','c_p','delta']

v_tt=[3]
v_cc=range(1,11)

# loop for type of function
for tt in v_tt:
	typen = cell_typen[tt-1]
	print(typen)
	# loop for case
	for cc in v_cc:
		# case code
		case = cell_case[cc-1]

		myParameter_Jun01V27 = {
			'pstar': 1,
			'b': 0.2,
			'phi': 0.6,
			'sigma': 0.06,
			'beta': 0.72,
			'lambda': 0.1,
			'c_f': 0.36,
			'r': 0.012,
			'c_p': 0.06,
			'delta': 0.53,

			'A': 1.355,
			'B1': 18.0,
			'B2': 25.0,
			'mu': 0.08,
			'epsilon_u':Epsilon_u(typen).Assign()
		}


		# array of variable
		vt=CaseCode(case).VT()
		str_fig = taskcode+'_'+typen+'_'+case

		vepsilon_d=[]
		vepsilon_c=[]
		vtheta=[]
		valpha=[]
		vflag=[]

		for i in range(0,len(vt)):
			# initialize parameters in each loop
			Paras = myParameter_Jun01V27
			Paras[case] = vt[i]
			tresult = Operator(case,typen,Paras).Solution()
			vepsilon_d.append(tresult[0][0])
			vepsilon_c.append(tresult[0][1])
			vtheta.append(tresult[0][2])
			valpha.append(tresult[0][3])
			vflag.append(tresult[1])

		np.savetxt(path_data+os.sep+str_fig+'_'+str_para+'_epsilon_d.txt',vepsilon_d)
		np.savetxt(path_data+os.sep+str_fig+'_'+str_para+'_epsilon_c.txt',vepsilon_c)
		np.savetxt(path_data+os.sep+str_fig+'_'+str_para+'_theta.txt',vtheta)
		np.savetxt(path_data+os.sep+str_fig+'_'+str_para+'_alpha.txt',valpha)
		np.savetxt(path_data+os.sep+str_fig+'_'+str_para+'_flag.txt',vflag)

		################################################
		#print typen
		#print CaseCode(case).VT()[1]
		#print Epsilon_u(typen).Assign()
		#print myParameter_Jun01V27['pstar']
		#basic
