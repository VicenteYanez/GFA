#! /usr/bin/python3

"""
@author: Francisco Garcia

Command line user interface for building trajectory models
"""

import sys
import os
import pdb

import numpy as np

from loadGPS import load1stations
from model_accions import *
import matplotlib.pyplot as plt


class color():
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


# parametros
codigo = sys.argv[1]
estation = sys.argv[2]

path = os.path.dirname(os.path.abspath(__file__))
directory = '{}/../example_files/{}/'.format(path, codigo)
model_file = '{}modelo/{}.txt'.format(directory, estation)

# parametros del nuevo modelo
ok = "n"
while ok in ['n', 'N']:
	stndr = input(color.BOLD+"Do you want to use the standar values? y/n ")
	print(" ")
	if stndr in ['y', 'Y', 's', 'S', 'O', 'o']:
		parametros = {"polinomio": 1,
		              "saltos": [2010.15, 2014.24832],
		              "Escala curva log": [1],
		              "Periodos Fourier": [1],
		              "Inicio log": [2010.15]}
		tsdata, modeldata = load1stations(estation, directory, True, model=False)
		modelo, residual = build_model(tsdata, parametros)
		save_model(modelo, model_file)
		upgrade_list(estation, parametros, residual, directory)
		pl = 1
		hs = [2010.15, 2014.24832]
		ls = [1]
		fp = [1]
		lb = [2010.15]

	else:
		"""
		polinomial part
		"""
		ok_pl = input(color.BOLD+"The stndr value for polinomial function is 1, do you want to change it? y/n "+color.END)
		print(" ")
		if ok_pl in ['y', 'Y', 's', 'S', 'O', 'o']:
			pl = int(input(color.BOLD+color.PURPLE+"Order for the polimolial function: "+color.END))
			print(" ")
		else:
			pl = 1
		"""
		jumps
		"""
		ok_hs = input(color.BOLD+"the stndr number of jumps is 2 [2010.01, 2014.24832], do you want to change it? y/n "+color.END)
		print(" ")
		if ok_hs in ['y', 'Y', 's', 'S', 'O', 'o']:

			ok_hs2 = input(color.BOLD+"do you want to select the dates? y/n"+color.END+", if not, last dates will be used ")
			print(" ")
			if ok_hs2 in ['y', 'Y', 's', 'S', 'O', 'o']:
				hs = []
				n_hs = int(input(color.BOLD+color.PURPLE+"Number of jumps: "+color.END))
				print(" ")
				for g in range(n_hs):
					hsd = float(input(color.BOLD+color.PURPLE+"Date for the "+str(g+1)+" heaviside jump: "+color.END))
					print(" ")
					hs.append(hsd)
			else:
				print(color.BOLD+color.RED+"The last values will be used for jumps "+str(hs)+color.END)
				print(" ")
		else:
			hs =[2010.15, 2014.24832]
		"""
		logs
		"""
		ok_lb = input(color.BOLD+"the stndr number of log curves is 2 [2010.01, 2014.24832] [1, 1], do you want to change it? y/n "+color.END)
		print(" ")
		if ok_lb in ['y', 'Y', 's', 'S', 'O', 'o']:

			ok_lb2 = input(color.BOLD+"do you want to select the dates? y/n, if n, last dates will be used "+color.END)
			print(" ")

			if ok_lb2 in ['y', 'Y', 's', 'S', 'O', 'o']:
				lb = []
				ls = []
				n_lb = int(input(color.BOLD+color.PURPLE+"Number of log curves: "+color.END))
				print(" ")
				for g in range(n_lb):
					lbd = float(input(color.BOLD+color.PURPLE+"Date for the "+str(g+1)+" log curve: "+color.END))
					print(" ")
					lb.append(lbd)

					lsd = float(input(color.BOLD+color.PURPLE+"Scale for the "+str(g+1)+" log curve: "+color.END))
					print(" ")
					ls.append(lsd)
			else:
				print(color.BOLD+color.RED+"The last values will be used for log curves"+str(lb)+color.END)
				print(" ")
				print(color.BOLD+color.RED+"with theese scales"+str(ls)+color.END)
				print(" ")
		else:
			lb =[2010.01, 2014.24832]
			ls =[1,1]
		"""
		fourier
		"""
		ok_fp = input(color.BOLD+"the stndr number of Fourier periods is 1 with a period of 1 year, do you want to change it? y/n "+color.END)
		print(" ")
		if ok_fp in ['y', 'Y', 's', 'S', 'O', 'o']:
			fp = []
			n_fp = int(input(color.BOLD+color.PURPLE+"Number of fourier periods: "+color.END))
			print(" ")
			for f in range(n_fp):
				fpd = float(input(color.BOLD+color.PURPLE+"period for the "+str(f+1)+" fourier serie: "+color.END))
				print(" ")
				fp.append(fpd)

		else:
			fp = [1]
		print(pl)
		print(hs)
		print(ls)
		print(fp)
		print(lb)

		parametros = {"polinomio": pl,
		              "saltos": hs,
		              "Escala curva log": ls,
		              "Periodos Fourier":  fp,
		              "Inicio log": lb}
		tsdata, modeldata = load1stations(estation, directory, True, model=False)
		modelo, residual = build_model(tsdata, parametros)
		save_model(modelo, model_file)
		upgrade_list(estation, parametros, residual, directory)

	ts_dates=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[0])
	ts_e    =np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[1])
	ts_n    =np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[2])
	ts_z    =np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[3])
	ts_e_err=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[4])
	ts_n_err=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[5])
	ts_z_err=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[6])
	tm_dates=np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+".txt", usecols=[0])
	tm_e    =np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+".txt", usecols=[1])
	tm_n    =np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+".txt", usecols=[2])
	tm_z    =np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+".txt", usecols=[3])
	"""
	East component
	"""
	plt.subplot(311)
	plt.errorbar(ts_dates, ts_e, yerr=ts_e_err, markersize=5., fmt='bo', ecolor='b', elinewidth=1.)
	plt.plot(tm_dates, tm_e, '-', c='g', linewidth=2., alpha=1.)
	plt.grid(True)
	plt.ylabel('Displacement mm')
	plt.title('Time serie and trajectory model for '+estation)
	"""
	North component
	"""
	plt.subplot(312)
	plt.errorbar(ts_dates, ts_n, yerr=ts_n_err, markersize=5., fmt='bo', ecolor='b', elinewidth=1.)
	plt.plot(tm_dates, tm_n, '-', c='g', linewidth=2., alpha=1.)
	plt.grid(True)

	"""
	Vertical component
	"""
	plt.subplot(313)
	plt.errorbar(ts_dates, ts_z, yerr=ts_z_err, markersize=5., fmt='bo', ecolor='b', elinewidth=1.)
	plt.plot(tm_dates, tm_z, '-', c='g', linewidth=2., alpha=1.)
	plt.grid(True)
	plt.show()

	ok = input(color.BOLD+"Do you like the result? y/n "+color.END)
	print(" ")

if 'parametros' in locals():
    ok_no_f = input(color.BOLD+"Do you want to compute a solution with no seasonal component? y/n "+color.END)
    print(" ")
    if ok_no_f in ['y', 'Y', 's', 'S', 'O', 'o']:
    	model_file = '{}modelo/{}_nof.txt'.format(directory, estation)
    	parametros = {"polinomio": pl,
    			      "saltos": hs,
    			      "Escala curva log": ls,
    			      "Periodos Fourier":  fp,
    			      "Inicio log": lb}
    	tsdata, modeldata = load1stations(estation, directory, True, model=False)
    	modelo, residual = build_model_sf(tsdata, parametros)
    	save_model(modelo, model_file)
    	upgrade_list(estation, parametros, residual, directory)

    	ts_dates=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[0])
    	ts_e    =np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[1])
    	ts_n    =np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[2])
    	ts_z    =np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[3])
    	ts_e_err=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[4])
    	ts_n_err=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[5])
    	ts_z_err=np.loadtxt("../example_files/"+codigo+"/series/"+estation+".txt", usecols=[6])
    	tm_dates=np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+"_nof.txt", usecols=[0])
    	tm_e    =np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+"_nof.txt", usecols=[1])
    	tm_n    =np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+"_nof.txt", usecols=[2])
    	tm_z    =np.loadtxt("../example_files/"+codigo+"/modelo/"+estation+"_nof.txt", usecols=[3])
    	"""
    	East component
    	"""
    	plt.subplot(311)
    	plt.errorbar(ts_dates, ts_e, yerr=ts_e_err, markersize=5., fmt='bo', ecolor='b', elinewidth=1.)
    	plt.plot(tm_dates, tm_e, '-', c='g', linewidth=2., alpha=1.)
    	plt.grid(True)
    	plt.ylabel('Displacement mm')
    	plt.title('Time serie and trajectory model for '+estation)
    	"""
    	North component
    	"""
    	plt.subplot(312)
    	plt.errorbar(ts_dates, ts_n, yerr=ts_n_err, markersize=5., fmt='bo', ecolor='b', elinewidth=1.)
    	plt.plot(tm_dates, tm_n, '-', c='g', linewidth=2., alpha=1.)
    	plt.grid(True)

    	"""
    	Vertical component
    	"""
    	plt.subplot(313)
    	plt.errorbar(ts_dates, ts_z, yerr=ts_z_err, markersize=5., fmt='bo', ecolor='b', elinewidth=1.)
    	plt.plot(tm_dates, tm_z, '-', c='g', linewidth=2., alpha=1.)
    	plt.grid(True)
    	plt.show()
    else:
    	print(color.BOLD+color.RED+"OK, good luck!"+color.END)
    	print(" ")
else:
    print('Bye')
