# main_Jun01T1_Plot.py

# Python PLOT program for Task 1 of Jun 01,2017
# June 1,2017
# May 22-24, 2017
# 20170411 PM
# Mar18,2017
# Task #20170221
# Related to May16/Oct17,2016
# edited in Feb21,2017

import numpy as np
import matplotlib.pyplot as plt
from itertools import compress
import os

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

# main part
taskcode = 'Jun01T1'
str_para = 'paraJun01V27'
# Example on Windows:
# source = ['"C:\\My Documents"','C:\\Code']
# for Mac OS or Linux
path_fig = '/Users/huwei/Dropbox/On_Local/Python/20170620 Python'
if not os.path.exists(path_fig):
	print('path_fig does not exist!')

# Example on Windows:
# source = ['"C:\\My Documents"','C:\\Code']
# for Mac OS or Linux
path_data = '/Users/huwei/Dropbox/On_Local/Python/20170620 Python'
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
	# loop for case
	for cc in v_cc:
		# case code
		case = cell_case[cc-1]

		# array of variable
		vt=CaseCode(case).VT()
		str_fig = taskcode+'_'+typen+'_'+case

		#########################
		vepsilon_d = np.loadtxt(path_data+os.sep+str_fig+'_'+str_para+'_epsilon_d.txt')
		vepsilon_c = np.loadtxt(path_data+os.sep+str_fig+'_'+str_para+'_epsilon_c.txt')
		vtheta = np.loadtxt(path_data+os.sep+str_fig+'_'+str_para+'_theta.txt')
		valpha = np.loadtxt(path_data+os.sep+str_fig+'_'+str_para+'_alpha.txt')
		vflag = np.loadtxt(path_data+os.sep+str_fig+'_'+str_para+'_flag.txt')
		##########################
		alpha_flag=np.logical_and(\
			np.array(valpha)>0,\
			np.array(valpha)<1)

		theta_flag=np.array(vtheta)>0

		flag = np.logical_and(\
			np.logical_and(\
			np.array(vflag)==1,\
			alpha_flag),theta_flag)

		# Set the font dictionaries (for plot title and axis titles)
		title_font = {'fontname':'Arial', 'size':'16', 'color':'black', 'weight':'normal',\
			'verticalalignment':'bottom'}
			# Bottom vertical alignment for more space
		axis_font = {'fontname':'Arial', 'size':'18'}
		# print flag
		f1 = plt.figure(1,figsize=(16,8),dpi=100)
		plt.subplot(221)
		plt.plot(list(compress(vt, flag)),\
			list(compress(vepsilon_d,flag)),\
			linewidth=2)
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		plt.xlabel(r'$ '+CaseCode(case).Str_Var()+r' $',**axis_font)
		plt.ylabel(r'$ \epsilon^{d} $',**axis_font)
		str_title = r'$ '+taskcode+','+typen+','+str_para+r' $'
		plt.title(str_title, **title_font)
		plt.xlim((min(list(compress(vt, flag))),\
			max(list(compress(vt, flag)))))

		plt.subplot(222)
		plt.plot(list(compress(vt, flag)),\
			list(compress(vepsilon_c,flag)),\
			linewidth=2)
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		plt.xlabel(r'$ '+CaseCode(case).Str_Var()+r' $',**axis_font)
		plt.ylabel(r'$ \epsilon^{c} $',**axis_font)
		str_title = r'$'+taskcode+','+typen+','+str_para+'$'
		plt.title(str_title, **title_font)
		plt.xlim((min(list(compress(vt, flag))),\
			max(list(compress(vt, flag)))))

		plt.subplot(223)
		plt.plot(list(compress(vt, flag)),\
			list(compress(vtheta,flag)),\
			linewidth=2)
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		plt.xlabel(r'$ '+CaseCode(case).Str_Var()+r' $',**axis_font)
		plt.ylabel(r'$ \theta $',**axis_font)
		plt.xlim((min(list(compress(vt, flag))),\
			max(list(compress(vt, flag)))))

		plt.subplot(224)
		plt.plot(list(compress(vt, flag)),\
			list(compress(valpha,flag)),\
			linewidth=2)
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		plt.xlabel(r'$ '+CaseCode(case).Str_Var()+r' $',**axis_font)
		plt.ylabel(r'$ \alpha $',**axis_font)
		plt.xlim((min(list(compress(vt, flag))),\
			max(list(compress(vt, flag)))))

		str_fig2 = path_fig+os.sep+str_fig+'_'+str_para+'_sol.png'

		f1.savefig(str_fig2);
		plt.close()
		# for Windows
		# mng = plt.get_current_fig_manager()
		# mng.frame.Maximize(True)
		#plt.show()
