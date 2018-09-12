import os, sys, numpy, time, ttk, traceback, tkMessageBox, matplotlib, itertools, collections, subprocess, operator, math
import Tkinter as Tk
matplotlib.use('TkAgg')
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
import matplotlib.pyplot as plt
from tkFileDialog import askopenfilename
from optparse import OptionParser
from scipy import stats

parser = OptionParser()
parser.add_option("--v", action='store_true',help="Show Spinnaker output. Default = false",default=False,dest='show_spin_out')
parser.add_option("--d", action='store_true',help="Show debug output. Default = false",default=False,dest='debug_out')
(options, args) = parser.parse_args()

plotlabelsize=20

show_spin_out = options.show_spin_out
debug_out = options.debug_out
#define root window
root = Tk.Tk()
bgcolor = '#BFBFBF'
Tk.Grid.columnconfigure(root,0,weight=1)
Tk.Grid.rowconfigure(root,0,weight=1)
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("%dx%d+0+0" % (w*0.4, h*0.75))
root.wm_title("Spinteractive v3")
root.config(bg=bgcolor)
root.resizable(width=True,height=True)
#root.geometry("%dx%d+0+0" % (w, h))

mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['ytick.major.pad'] = 8
mpl.rcParams['xtick.major.pad'] = 8

global variables
variables = ['def_cri','def_v0','def_hb1','def_hb2','def_b0','def_b1','def_Z1','def_hv','def_gmode','def_D0','def_mu','def_modestr','def_ablossesstr','def_velfield','def_nu1','def_nu2','def_nu3','def_nu4','def_norm','def_epsilon','def_fwhm']

juststarted = True

def quitit():
	pars = [def_cri,def_v0,def_hb1,def_hb2,def_b0,def_b1,def_Z1,def_hv,def_gmode.get(),def_D0,def_mu,def_modestr.get(),def_ablossesstr.get(),def_velfield.get(),frequency1field.get(),frequency2field.get(),frequency3field.get(),frequency4field.get(),def_norm.get(),def_epsilon.get(),FWHMfield.get()]
	numpy.savetxt('pars.last',pars,fmt='%s')
	sys.exit()

def saveplotpars(fname):
	pars = [["Gamma:",def_cri],["b0:",def_b0],["b1:",def_b1],["hb1:",def_hb1],["hb2:",def_hb2],["Z1:",def_Z1],["galaxy mode:",def_gmode.get()],["mode:",def_modestr.get()],["v0:",def_v0],["hv:",def_hv],["abiatic losses:",def_ablossesstr.get()],["velocity field:",def_velfield.get()],["D0:",def_D0],["mu:",def_mu],["freq 1:",frequency1field.get()],["freq 2:",frequency2field.get()],["freq 3:",frequency3field.get()],["freq 4:",frequency4field.get()],['normalize:',def_norm.get()],['epsilon:',def_epsilon.get()],["rChi^2_1: ",chi_array[0]],["rChi^2_2: ",chi_array[1]],["rChi^2_spi1: ",chi_array[2]],["RMS chi: ",get_global_chi()],["fwhm: ",def_fwhm]]
	numpy.savetxt(fname,pars,fmt='%s')

#define needed variables or load last used values
try:
	lastpars = numpy.loadtxt('pars.last',dtype=object)
	def_cri= float(lastpars[0])
	def_v0 = float(lastpars[1])
	def_hb1= float(lastpars[2])
	def_hb2= float(lastpars[3])
	def_b0 = float(lastpars[4])
	def_b1 = float(lastpars[5])
	def_Z1 = float(lastpars[6])
	def_hv = float(lastpars[7])
	def_gmode = Tk.IntVar()
	def_gmode.set(int(lastpars[8]))


	def_D0 = float(lastpars[9])
	def_mu = float(lastpars[10])

	def_modestr = Tk.StringVar()
	def_modestr.set(lastpars[11])

	if def_modestr.get()=="Advection":
		def_mode = 1
	else:
		def_mode = 2

	def_ablossesstr = Tk.StringVar()
	def_ablossesstr.set(lastpars[12])

	if def_ablossesstr.get()=="Yes":
		def_ablosses = 1
	else:
		def_ablosses = -1

	def_velfield = Tk.IntVar()
	def_velfield.set(int(lastpars[13]))
	def_nu1 = lastpars[14]
	def_nu2 = lastpars[15]
	def_nu3 = lastpars[16]
	def_nu4 = lastpars[17]
	def_norm = Tk.IntVar()
	def_norm.set(int(lastpars[18]))
	def_epsilon = Tk.IntVar()
	def_epsilon.set(int(lastpars[19]))
	def_fwhm = float(lastpars[20])

	print 'Parameters from last file loaded succesfully'
except:
	if debug_out:
			traceback.print_exc()
	print 'No last file found, loading defaults'
	def_modestr = Tk.StringVar()
	def_modestr.set("Advection")
	def_gmode = Tk.IntVar()
	def_gmode.set(1)
	def_cri=2.3
	def_v0=150
	def_b0=10
	def_hb1=0.5
	def_hb2=4
	def_b1 = 7
	def_hv = 1
	def_Z1 = 10
	def_D0 = 2.0
	def_mu = 0.5
	def_velfield = Tk.IntVar()
	def_velfield.set(0)
	def_ablossesstr = Tk.StringVar()
	def_ablossesstr.set("Yes")
	if def_ablossesstr.get()=="Yes":
		def_ablosses = 1
	else:
		def_ablosses = -1

	if def_modestr.get()=="Advection":
		def_mode = 1
	else:
		def_mode = 2
	def_nu1 = "1.5e9"
	def_nu2 = "6.0e9"
	def_nu3 = "6.20e9"
	def_nu4 = "8.40e9"
	def_epsilon = Tk.IntVar()
	def_epsilon.set(-1)
	def_norm = Tk.IntVar()
	def_norm.set(1)
	def_fwhm = 1.2

	##########################################################

def plot_dummylegend():
#	global def_velfield
	if def_ablossesstr.get()=="Yes":
		def_ablosses = 1
	else:
		def_ablosses = -1
	
	stuff = range(len(variables))
	for i in range(len(stuff)):
		stuff[i]=''
	stuff[0]=r'$\gamma$: '+str(round(def_cri,3))
	stuff[1]=r'$b_0$: '+str(round(def_b0,3))+'$\mu$G'
	stuff[2]=r'$b_1$: '+str(round(def_b1,3))+'$\mu$G'
	stuff[3]=r'$hb_1$: '+str(round(def_hb1,3))+' kpc'
	stuff[4]=r'$hb_2$: '+str(round(def_hb2,3))+' kpc'
	stuff[5]=r'$Z_1$: '+str(round(def_Z1,3))+' kpc'
	stuff[6]='gmode: '+str(def_gmode.get())
	stuff[8]='mode: '+str(def_modestr.get())
	if def_modestr.get()=='Advection':
		stuff[9]=r'$v_0$: '+str(round(def_v0,3))+' km/s'
		stuff[10]=r'$h_V$: '+str(round(def_hv,3))+' kpc'
		stuff[11]='Ad. Losses: '+str(def_ablossesstr.get())
		stuff[12]='Vel. Field: '+str(def_velfield.get())
	else:
		stuff[9]=r'$D_0$: '+str(round(def_D0,3))
		stuff[10]=r'$\mu_{diff}$: '+str(round(def_mu,3))
	stuff[13] = 'Epsilon: '+str(def_epsilon.get())
	stuff[14] = 'Normalization: '+str(def_norm.get())
	stuff[15] = 'RMS '+r'$\chi^2$: ' + str(round(get_global_chi(),2))
	l1 = ax0.legend(stuff,loc=10,handlelength=0, handletextpad=0,markerscale=0,fontsize=22,scatterpoints=1)
	for i in range(len(variables)):
		l1.legendHandles[i].set_visible(False)
	parametercanvas.draw()

def get_global_chi():
	global chi_array
	chis = numpy.array(chi_array)
	if len(chi_array)!=0:
		totalchi = (numpy.nanmean(chis**2))**0.5
		totalchi = numpy.nanmean(chis)
		return totalchi
	else:
		return 0

frequencies = []
simplot_i1 = None
simplot_i2 = None
simplot_i3 = None
simplot_i4 = None

simplot_a12 = None
simplot_a23 = None
simplot_a34 = None

chi_array = []

spi1 = None
spi2 = None
spi3 = None

chiplot_1 = None
chiplot_2 = None
chiplot_3 = None
chiplot_4 = None
chiplot_5 = None
chiplot_6 = None
chiplot_7 = None

haserror = False

############# setup parameter frame and plot
parameterframe = Tk.Frame(master=root)	
parameterframe.place(relx=0,rely=0,relheight=1,relwidth=0.15)

parameterfig = plt.figure()
parametercanvas = FigureCanvasTkAgg(parameterfig, master=parameterframe)


ax0 = parameterfig.add_axes([0.0,0.0,1,1])
ax0.axis('off')
for i in range(len(variables)+4):
	ax0.plot(0,0)
plot_dummylegend()
parametercanvas.show()
parametercanvas.get_tk_widget().place(relx=0,rely=0,relwidth=1,relheight=1)
##########################################################

############# setup main plot frame

plotframe = Tk.Frame(master=root)
plotframe.place(relx=0.15,rely=0,relheight=1,relwidth=0.705)
plotframe.config(bg=bgcolor)
##########################################################

plots = []
axes = []
for i in range(7):
	plots.append(None)


def onclick(event):
		distlabelvar.set(str(round(event.xdata,2))+' kpc: '+str(round(event.ydata,2)))

plotfig = plt.figure()

plotcanvas = FigureCanvasTkAgg(plotfig, master=plotframe)
plotcanvas.get_tk_widget().place(relx=0,rely=0,relwidth=1,relheight=1)
plotcanvas.show()
cid = plotfig.canvas.mpl_connect('button_press_event', onclick)
simplot_a12 = None

def onclick(event):
		print event.xdata
		print event.ydata

def calc_spi(freq1,freq2,data1,data2):
	global show_spin_out, debug_out
	try:
		spi = numpy.log(data1/data2)/numpy.log(freq1/freq2)
		return spi
	except:
		if debug_out:
			traceback.print_exc()
		return [0]


def initial_plot():
	global simplot_i1, simplot_i2, simplot_i3, simplot_i4, simplot_a12, simplot_a23, simplot_a34, chi_array, spi1, spi2, spi3, data, plotfig, show_spin_out, debug_out, chi_array, chiplot_1,  chiplot_2, chiplot_3, chiplot_4, chiplot_5, chiplot_6, chiplot_7,def_norm, data, spi1err, spi2err, spi3err, spie4err
	#hack to have frequencies in an array for plotting purpouses
	frequencies = []
	frequencies.append(float(frequency1field.get()))
	frequencies.append(float(frequency2field.get()))	
	frequencies.append(float(frequency3field.get()))
	frequencies.append(float(frequency4field.get()))

	numfig = sum(filenames)+sum(filenames)-1
	maxheight = 0.85/numfig
	try:
		d2 = numpy.array(data)
	except ValueError:
		tkMessageBox.showinfo("Error", "Array shape mismatch. Does one file have errorbars and others do not?")
		if debug_out:
			traceback.print_exc()
		return
	minx = numpy.nanmin(d2[:,:,0])
	maxx = numpy.nanmax(d2[:,:,0])
	plotfig.clf()
	for i in range(numfig-1,-1,-1):
		plotfig.add_axes([0.1,(i*maxheight)+0.1,0.85,maxheight])

	for i in range(len(plotfig.axes)):
		plotfig.axes[i].grid()
		pad = (maxx-minx)*0.05/1000.
#		plotfig.axes[i].set_xlim(minx/1000-pad,maxx/1000+pad)
		plotfig.axes[i].set_xlim(0.-pad,15.+pad)
		if i%2==0:
			plotfig.axes[i].yaxis.tick_right()


	for i in range(0,numfig-1):
		plotfig.axes[i].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')

	for i in range(d2.shape[0]):
		if int(def_norm.get())==1:
			normfac = numpy.max(d2[i,:,1])
		else:
			normfac = 1
		global haserror
		if d2[i,:].shape[1]==3:
			plotfig.axes[i].errorbar(d2[i,:,0]/1000,d2[i,:,1]/normfac,d2[i,:,2]/normfac,ls='none')
			plotfig.axes[i].scatter(d2[i,:,0]/1000,d2[i,:,1]/normfac,s=50,label=str(frequencies[i]/1e9)+' GHz')
			haserror = True
		else:
			plotfig.axes[i].plot(d2[i,:,0]/1000,d2[i,:,1]/normfac,label=str(frequencies[i]/1e9)+' GHz')
			haserror = False

	if d2.shape[0]==1:
		try:
			simplot_i1, = plotfig.axes[0].plot(d2[0,:,0]/1000,numpy.zeros(shape=len(d2[0,:,0])),color='red',label=r'Sim')
			xlim = plotfig.axes[0].get_xlim()
			ylim = plotfig.axes[0].get_ylim()
			chiplot_1, = plotfig.axes[0].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			plotfig.axes[0].set_ylabel('Normalized intensity',fontsize=18)
			plotfig.axes[0].set_xlabel('Distance from disk [kpc]',fontsize=18)
		except:
			if debug_out:
				traceback.print_exc()
			pass

	if d2.shape[0]==2:
		try:
			spi1 = calc_spi(numpy.float(frequency1field.get()),numpy.float(frequency2field.get()),d2[0,:,1],d2[1,:,1])
			if haserror:
				spi1err = (1/numpy.log(numpy.float(frequency1field.get())/numpy.float(frequency2field.get())))*((d2[0,:,2]/d2[0,:,1])**2+(d2[1,:,2]/d2[1,:,1])**2)**0.5
				plotfig.axes[2].errorbar(d2[0,:,0]/1000.,spi1,spi1err,color='green',ls='none')
#				plotfig.axes[2].scatter(d2[0,:,0]/1000.,spi1,s=50,color='green',label=r'$\alpha$ 1-2')
				plotfig.axes[2].scatter(d2[0,:,0]/1000.,spi1,s=50,color='green',label=r'$\alpha$')
			else:
				plotfig.axes[2].scatter(d2[0,:,0]/1000.,spi1,color='green',label=r'$\alpha$ 1-2')
		except:
			if debug_out:
				traceback.print_exc()
			pass

		try:
			simplot_i1, = plotfig.axes[0].plot(d2[0,:,0]/1000.,spi1*0,color='red',label=r'Sim')
			simplot_i2, = plotfig.axes[1].plot(d2[0,:,0]/1000.,spi1*0,color='red',label=r'Sim')
			simplot_a12, = plotfig.axes[2].plot(d2[0,:,0]/1000.,spi1,color='red',label=r'Sim')
			xlim = plotfig.axes[0].get_xlim()
			ylim = plotfig.axes[0].get_ylim()
			chiplot_1, = plotfig.axes[0].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[1].get_xlim()
			ylim = plotfig.axes[1].get_ylim()
			chiplot_2, = plotfig.axes[1].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[2].get_xlim()
			ylim = plotfig.axes[2].get_ylim()
			chiplot_3, = plotfig.axes[2].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			plotfig.axes[2].set_xlabel('Distance from disk [kpc]',fontsize=18)
#			plotfig.axes[0].set_ylabel('Normalized Intensity',fontsize=18)
			plotfig.axes[1].set_ylabel('Normalized intensity',fontsize=18,y=1)
			plotfig.axes[2].set_ylabel('Spectral index',fontsize=18)
		except:
			if debug_out:
				traceback.print_exc()
			pass

	if d2.shape[0]==3:
		try:
			spi1 = calc_spi(numpy.float(frequency1field.get()),numpy.float(frequency2field.get()),d2[0,:,1],d2[1,:,1])
			if haserror:
				spi1err = (1/numpy.log(numpy.float(frequency1field.get())/numpy.float(frequency2field.get())))*((d2[0,:,2]/d2[0,:,1])**2+(d2[1,:,2]/d2[1,:,1])**2)**0.5
				plotfig.axes[3].errorbar(d2[0,:,0],spi1,spi1err,color='green',label=r'$\alpha$ 1-2',ls='none')
				plotfig.axes[3].scatter(d2[0,:,0],spi1,s=50)
			else:
				plotfig.axes[3].scatter(d2[0,:,0],spi1,color='green',label=r'$\alpha$')
			plotfig.axes[3].legend(scatterpoints=1,fontsize=22)
		except:
			if debug_out:
				traceback.print_exc()
			pass

		try:
			spi2 = calc_spi(numpy.float(frequency2field.get()),numpy.float(frequency3field.get()),d2[1,:,1],d2[2,:,1])
			if haserror:
				spi2err = (1/numpy.log(numpy.float(frequency1field.get())/numpy.float(frequency2field.get())))*((d2[1,:,2]/d2[1,:,1])**2+(d2[2,:,2]/d2[2,:,1])**2)**0.5
				plotfig.axes[4].errorbar(d2[0,:,0],spi2,spi2err,color='green',label=r'$\alpha$ 2-3',ls='none')
				plotfig.axes[4].scatter(d2[0,:,0],spi2,s=50)
			else:
				plotfig.axes[4].scatter(d2[0,:,0],spi2,color='green',label=r'$\alpha$ 2-3')
			plotfig.axes[4].legend(scatterpoints=1,fontsize=22)
		except:
			if debug_out:
				traceback.print_exc()
			pass

		try:
			simplot_i1, = plotfig.axes[0].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_i2, = plotfig.axes[1].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_i3, = plotfig.axes[2].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_a12, = plotfig.axes[3].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_a23, = plotfig.axes[4].plot(d2[0,:,0],spi2,color='red',label=r'Sim')
			xlim = plotfig.axes[0].get_xlim()
			ylim = plotfig.axes[0].get_ylim()
			chiplot_1, = plotfig.axes[0].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[1].get_xlim()
			ylim = plotfig.axes[1].get_ylim()
			chiplot_2, = plotfig.axes[1].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[2].get_xlim()
			ylim = plotfig.axes[2].get_ylim()
			chiplot_3, = plotfig.axes[2].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[3].get_xlim()
			ylim = plotfig.axes[3].get_ylim()
			chiplot_4, = plotfig.axes[3].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[4].get_xlim()
			ylim = plotfig.axes[4].get_ylim()
			chiplot_5, = plotfig.axes[4].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)

		except:
			if debug_out:
				traceback.print_exc()
			pass


	if d2.shape[0]==4:
		try:
			spi1 = calc_spi(numpy.float(frequency1field.get()),numpy.float(frequency2field.get()),d2[0,:,1],d2[1,:,1])
			if haserror:
				spi1err = (1/numpy.log(numpy.float(frequency1field.get())/numpy.float(frequency2field.get())))*((d2[0,:,2]/d2[0,:,1])**2+(d2[1,:,2]/d2[1,:,1])**2)**0.5
				plotfig.axes[4].errorbar(d2[0,:,0],spi1,spi1err,color='green',label=r'$\alpha$ 1-2',ls='none')
				plotfig.axes[4].scatter(d2[0,:,0],spi1,s=50,color='green')
			else:
				plotfig.axes[4].scatter(d2[0,:,0],spi1,color='green',label=r'$\alpha$ 1-2')
			plotfig.axes[4].legend(scatterpoints=1)
		except:
			if debug_out:
				traceback.print_exc()
			pass

		try:
			spi2 = calc_spi(numpy.float(frequency2field.get()),numpy.float(frequency3field.get()),d2[1,:,1],d2[2,:,1])
			if haserror:
				spi2err = (1/numpy.log(numpy.float(frequency1field.get())/numpy.float(frequency2field.get())))*((d2[1,:,2]/d2[1,:,1])**2+(d2[2,:,2]/d2[2,:,1])**2)**0.5
				plotfig.axes[5].errorbar(d2[0,:,0],spi2,spi2err,color='green',label=r'$\alpha$ 2-3')
				plotfig.axes[5].scatter(d2[0,:,0],spi2,s=50,color='green')
			else:
				plotfig.axes[5].scatter(d2[0,:,0],spi2,color='green',label=r'$\alpha$ 2-3',ls='none')
			plotfig.axes[5].legend(scatterpoints=1)
		except:
			if debug_out:	
				traceback.print_exc()
			pass

		try:
			spi3 = calc_spi(numpy.float(frequency3field.get()),numpy.float(frequency4field.get()),d2[2,:,1],d2[3,:,1])
			if haserror:
				spi3err = (1/numpy.log(numpy.float(frequency1field.get())/numpy.float(frequency2field.get())))*((d2[2,:,2]/d2[2,:,1])**2+(d2[3,:,2]/d2[3,:,1])**2)**0.5
				plotfig.axes[6].errorbar(d2[0,:,0],spi3,spi3err,color='green',label=r'$\alpha$ 2-4')
				plotfig.axes[6].scatter(d2[0,:,0],spi3,s=50,color='green')
			else:
				plotfig.axes[6].plot(d2[0,:,0],spi3,color='green',label=r'$\alpha$ 2-4',ls='none')
			plotfig.axes[6].legend(scatterpoints=1)
		except:
			if debug_out:
				traceback.print_exc()
			pass

		try:
			simplot_i1, = plotfig.axes[0].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_i2, = plotfig.axes[1].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_i3, = plotfig.axes[2].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_i4, = plotfig.axes[3].plot(d2[0,:,0],spi1,color='red',label=r'Sim')
			simplot_a12, = plotfig.axes[4].plot(d2[0,:,0],spi1,color='red',label=r'Sim 1-2')
			simplot_a23, = plotfig.axes[5].plot(d2[0,:,0],spi2,color='red',label=r'Sim 2-3')
			simplot_a34, = plotfig.axes[6].plot(d2[0,:,0],spi3,color='red',label=r'Sim 3-4')
			xlim = plotfig.axes[0].get_xlim()
			ylim = plotfig.axes[0].get_ylim()
			chiplot_1, = plotfig.axes[0].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[1].get_xlim()
			ylim = plotfig.axes[1].get_ylim()
			chiplot_2, = plotfig.axes[1].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[2].get_xlim()
			ylim = plotfig.axes[2].get_ylim()
			chiplot_3, = plotfig.axes[2].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[3].get_xlim()
			ylim = plotfig.axes[3].get_ylim()
			chiplot_4, = plotfig.axes[3].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[4].get_xlim()
			ylim = plotfig.axes[4].get_ylim()
			chiplot_5, = plotfig.axes[4].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[5].get_xlim()
			ylim = plotfig.axes[5].get_ylim()
			chiplot_6, = plotfig.axes[5].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)
			xlim = plotfig.axes[6].get_xlim()
			ylim = plotfig.axes[6].get_ylim()
			chiplot_7, = plotfig.axes[6].plot(numpy.mean(xlim),numpy.mean(ylim),label=r"$\chi$="+str(0),alpha=0)

		except:
			if debug_out:
				traceback.print_exc()
			pass

	for i in range(len(plotfig.axes)):
		plotfig.axes[i].legend(prop={'size': 12},scatterpoints=1,fontsize=22)
		
	global def_scale
	if def_scale.get()=='Log':
		for i in range(sum(filenames)):
			plotfig.axes[i].set_yscale('log')
	else:
		for i in range(sum(filenames)):
			plotfig.axes[i].set_yscale('linear')
	
	run()
	
	plotcanvas.draw()
	plot_dummylegend()


def chisq(x,y,error):
   	x = numpy.array(x)
	y = numpy.array(y)
	#x is measurement, y is model
	result = numpy.nansum( ((x-y)/error)**2 )
	deg = 1 #only one parameter is changed at a time
	reduced = result/(len(x)-1-deg) 
	return reduced

def run():
	global filenames, def_mode, def_v0, def_cri, def_b0, def_b1, def_hb1, def_hb2, def_D0, def_gmode, def_Z1, def_hv, def_mode, def_ablosses, def_mu, def_velfield, def_norm, def_epsilon, data, plotfig, simplot_i1, simplot_i2, simplot_i3, simplot_i4, simplot_a12, simplot_a23, simplot_a34, chi_array, show_spin_out, debug_out, chiplot_1,  chiplot_2, chiplot_3, chiplot_4, chiplot_5, chiplot_6, chiplot_7, tosearch, best_results_array, def_fwhm
	global applyvalues
	try:
		if (applyvalues==1):
			for i in range(len(tosearch)):
				exec('{0}={1}'.format(tosearch[i],best_results_array[i]),globals())
			del tosearch,best_results_array
			applyvalues=0
	except:
		pass
	#define parameters to run spinnaker
	parfile = open('parameters-template','r')
	parameters = parfile.read()
	parfile.close()
	pars2 = parameters
	d2 = numpy.array(data)
	halosize = numpy.max(d2[0,:,0])/1000.0
	gridsize = len(d2[0,:,0])-1

	###new gridsize
	# length of datapoints -1
	data_length = len(d2[0,:,0])-1
	lim = 200/data_length
	gridsize = lim*data_length
	steps = 1.0*max(d2[0,:,0])/gridsize

	try:
		nu1 = numpy.float(frequency1field.get())
		pars2 = pars2.replace('nu_1 = 1.5e9# [Hz]','nu_1 = '+str(nu1)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()
		pass

	try:
		nu2 = numpy.float(frequency2field.get())
		pars2 = pars2.replace('nu_2 = 6.0e9# [Hz]','nu_2 = '+str(nu2)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()
		pass

	try:
		nu3 = numpy.float(frequency3field.get())
		pars2 = pars2.replace('nu_3 = 6.20e9# [Hz]','nu_3 = '+str(nu3)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()
		pass

	try:
		nu4 = numpy.float(frequency4field.get())
		pars2 = pars2.replace('nu_4 = 8.40e9# [Hz]','nu_4 = '+str(nu4)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()	
		pass
	pars2 = pars2.replace('z_halo = 8.','z_halo = '+str(halosize))	#enter here the parameter you wish to change
	pars2 = pars2.replace('grid_size = 200','grid_size = '+str(gridsize))
	pars2 = pars2.replace('V0 = 300.0e5','V0 = '+str(def_v0)+'e5')
	pars2 = pars2.replace('gamma_in = 2.7','gamma_in = '+str(def_cri))
	pars2 = pars2.replace('B0 = 13.5e-6','B0 = '+str(def_b0)+'e-6')
	pars2 = pars2.replace('B1 = 8.6e-6','B1 = '+str(def_b1)+'e-6')
	pars2 = pars2.replace('h_B1 = 0.6','h_B1 = '+str(def_hb1))
	pars2 = pars2.replace('h_B2 = 5.0','h_B2 = '+str(def_hb2))
	pars2 = pars2.replace('D0 = 2.0e28','D0 = '+str(def_D0)+'e28')
	pars2 = pars2.replace('galaxy_mode = 1','galaxy_mode = '+str(int(def_gmode.get())))
	pars2 = pars2.replace('z1 = 10.0','z1 = '+str(def_Z1))
	pars2 = pars2.replace('h_V = 1.','h_V = '+str(def_hv))
	pars2 = pars2.replace('mode = 1#Adv','mode = '+str(def_mode)+'#Adv')
	pars2 = pars2.replace('adiabatic_losses = -1','adiabatic_losses = '+str((int(def_ablosses))))
	pars2 = pars2.replace('mu_diff = 0.5','mu_diff = '+str(def_mu))
	pars2 = pars2.replace('velocity_field = 1','velocity_field = '+str(def_velfield.get()))
	pars2 = pars2.replace('normalize_intensities = -1','normalize_intensities = '+str(def_norm.get()))
	pars2 = pars2.replace('epsilon = -1','epsilon = '+str(def_epsilon.get()))
	pars2 = pars2.replace('FWHM_effective_beam = 1.2','FWHM_effective_beam = '+str(FWHMfield.get()))
	parfile = open('parameters','w')
	parfile.write(pars2)
	parfile.flush()
	parfile.close()
	#run spinnaker
	if show_spin_out:
		os.system('./spectral.x')
	else:
		os.popen('./spectral.x')

	simdata = numpy.loadtxt('int.dat')
	sim_int1 = simdata[:,1]
	sim_int2 = simdata[:,2]
	sim_int3 = simdata[:,3]
	sim_int4 = simdata[:,4]
	sim_a12 = simdata[:,5]
	sim_a23 = simdata[:,6]
	sim_a34 = simdata[:,7]
	d2 = numpy.array(data)

	#get newarray
	sim_int1_temp = []	
	sim_int2_temp = []
	sim_int3_temp = []
	sim_int4_temp = []
	sim_a12_temp = []
	sim_a23_temp = []
	sim_a34_temp = []
	for i in range(0,len(sim_int1)+1,lim):
		sim_int1_temp.append(sim_int1[i])
		sim_int2_temp.append(sim_int2[i])
		sim_int3_temp.append(sim_int3[i])
		sim_int4_temp.append(sim_int4[i])
		sim_a12_temp.append(sim_a12[i])
		sim_a23_temp.append(sim_a23[i])
		sim_a34_temp.append(sim_a34[i])

	sim_int1 = sim_int1_temp
	sim_int2 = sim_int2_temp
	sim_int3 = sim_int3_temp
	sim_int4 = sim_int4_temp
	sim_a12 = sim_a12_temp
	sim_a23 = sim_a23_temp
	sim_a34 = sim_a34_temp


	chi_array = numpy.zeros(shape=d2.shape[0]+d2.shape[0]-1)

#	if haserror:
#		print 'has'

	for i in range(d2.shape[0]):
		handles, labels = plotfig.axes[i].get_legend_handles_labels()
		labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
		plotfig.axes[i].legend(handles,labels,fontsize=plotlabelsize,scatterpoints=1)

	if d2.shape[0]==1:
		simplot_i1.set_ydata(sim_int1)
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1,d2[0,:,2]/numpy.nanmax(d2[0,:,1]))
#		print round(chi_array[0],4)
		plotfig.axes[0].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[0],2)))


	if d2.shape[0]==2:
		simplot_a12.set_ydata(sim_a12)
		simplot_i1.set_ydata(sim_int1)
		simplot_i2.set_ydata(sim_int2)
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1,d2[0,:,2]/numpy.nanmax(d2[0,:,1]))
		chi_array[1] = chisq(d2[1,:,1]/numpy.nanmax(d2[1,:,1]),sim_int2,d2[1,:,2]/numpy.nanmax(d2[1,:,1]))
		chi_array[2] = chisq(spi1,sim_a12,spi1err)
		plotfig.axes[0].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[0],2)))
		plotfig.axes[1].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[1],2)))
		####### reverse legend for spectral index
		handles, labels = plotfig.axes[2].get_legend_handles_labels()
		handles[0],handles[1]=handles[1],handles[0]
		handles[1],handles[2]=handles[2],handles[1]
		labels[1],labels[2]=labels[2],labels[1]
		labels[2]='Sim'
		plotfig.axes[2].legend(handles,labels,fontsize=plotlabelsize,scatterpoints=1)
		plotfig.axes[2].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[2],2)))
		

	if d2.shape[0]==3:
		simplot_a12.set_ydata(sim_a12)
		simplot_a23.set_ydata(sim_a23)
		simplot_i1.set_ydata(sim_int1)
		simplot_i2.set_ydata(sim_int2)
		simplot_i3.set_ydata(sim_int3)
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1)
		chi_array[1] = chisq(d2[1,:,1]/numpy.nanmax(d2[1,:,1]),sim_int2)
		chi_array[2] = chisq(d2[2,:,1]/numpy.nanmax(d2[2,:,1]),sim_int3)
		chi_array[3] = chisq(spi1,sim_a12)
		chi_array[4] = chisq(spi2,sim_a23)

		for i in range(3,5):
			handles, labels = plotfig.axes[i].get_legend_handles_labels()
			handles[0],handles[2]=handles[2],handles[0]
			labels[0],labels[2]=labels[2],labels[0]
			plotfig.axes[i].legend(handles,labels,fontsize=plotlabelsize)

		plotfig.axes[0].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[0],2)))
		plotfig.axes[1].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[1],2)))
		plotfig.axes[2].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[2],2)))
		plotfig.axes[3].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[3],2)))
		plotfig.axes[4].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[4],2)))


	if d2.shape[0]==4:
		simplot_a12.set_ydata(sim_a12)
		simplot_a23.set_ydata(sim_a23)
		simplot_a34.set_ydata(sim_a34)
		simplot_i1.set_ydata(sim_int1)
		simplot_i2.set_ydata(sim_int2)
		simplot_i3.set_ydata(sim_int3)
		simplot_i4.set_ydata(sim_int4)
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1)
		chi_array[1] = chisq(d2[1,:,1]/numpy.nanmax(d2[1,:,1]),sim_int2)
		chi_array[2] = chisq(d2[2,:,1]/numpy.nanmax(d2[2,:,1]),sim_int3)
		chi_array[3] = chisq(d2[3,:,1]/numpy.nanmax(d2[3,:,1]),sim_int4)
		chi_array[4] = chisq(spi1,sim_a12)
		chi_array[5] = chisq(spi2,sim_a23)
		chi_array[6] = chisq(spi3,sim_a34)
		for i in range(4,7):
			handles, labels = plotfig.axes[i].get_legend_handles_labels()
			handles[0],handles[2]=handles[2],handles[0]
			labels[0],labels[2]=labels[2],labels[0]
			plotfig.axes[i].legend(handles,labels,fontsize=plotlabelsize)

		plotfig.axes[0].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[0],2)))
		plotfig.axes[1].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[1],2)))
		plotfig.axes[2].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[2],2)))
		plotfig.axes[3].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[3],2)))
		plotfig.axes[4].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[4],2)))
		plotfig.axes[5].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[5],2)))
		plotfig.axes[6].legend_.texts[0].set_text(r"$\chi$="+str(round(chi_array[6],2)))

	plotcanvas.draw()

def run_noplot(varstochange,vals):
	global filenames, def_mode, def_v0, def_fwhm, def_cri, def_b0, def_b1, def_hb1, def_hb2, def_D0, def_gmode, def_Z1, def_hv, def_mode, def_ablosses, def_mu, def_velfield, def_norm, def_epsilon, data, chi_array
	#copy to local variables
	global loc_def_v0
	global loc_def_cri
	global loc_def_b0
	global loc_def_b1
	global loc_def_hb1
	global loc_def_hb2
	global loc_def_D0
	global loc_def_Z1
	global loc_def_hv
	global loc_def_mu
	global loc_def_fwhm

	loc_def_v0 = def_v0
	loc_def_cri = def_cri
	loc_def_b0 = def_b0
	loc_def_b1 = def_b1
	loc_def_hb1 = def_hb1
	loc_def_hb2 = def_hb2
	loc_def_D0 = def_D0
	loc_def_Z1 = def_Z1
	loc_def_hv = def_hv
	loc_def_mu = def_mu
	loc_def_fwhm = def_fwhm

	#define parameters to run spinnaker
	parfile = open('parameters-template','r')
	parameters = parfile.read()
	parfile.close()
	pars2 = parameters
	d2 = numpy.array(data)
	halosize = numpy.max(d2[0,:,0])/1000.0
	gridsize = len(d2[0,:,0])-1

	###new gridsize
	# length of datapoints -1
#	data_length = len(d2[0,:,0])-1
#	lim = 200/data_length
#	gridsize = lim*data_length

	data_length = len(d2[0,:,0])-1
	lim = 200/data_length
	gridsize = lim*data_length
	steps = 1.0*max(d2[0,:,0])/gridsize

	try:
		nu1 = numpy.float(frequency1field.get())
		pars2 = pars2.replace('nu_1 = 1.5e9# [Hz]','nu_1 = '+str(nu1)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()
		pass

	try:
		nu2 = numpy.float(frequency2field.get())
		pars2 = pars2.replace('nu_2 = 6.0e9# [Hz]','nu_2 = '+str(nu2)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()
		pass

	try:
		nu3 = numpy.float(frequency3field.get())
		pars2 = pars2.replace('nu_3 = 6.20e9# [Hz]','nu_3 = '+str(nu3)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()
		pass

	try:
		nu4 = numpy.float(frequency4field.get())
		pars2 = pars2.replace('nu_4 = 8.40e9# [Hz]','nu_4 = '+str(nu4)+'# [Hz]')		
	except:
		if debug_out:
			traceback.print_exc()	
		pass
	
	#change variables
	for v in range(len(varstochange)):
		exec("loc_{0}={1}".format(varstochange[v],vals[v]),globals())


	pars2 = pars2.replace('z_halo = 8.','z_halo = '+str(halosize))	#enter here the parameter you wish to change
	pars2 = pars2.replace('grid_size = 200','grid_size = '+str(gridsize))
	pars2 = pars2.replace('V0 = 300.0e5','V0 = '+str(loc_def_v0)+'e5')
	pars2 = pars2.replace('gamma_in = 2.7','gamma_in = '+str(loc_def_cri))
	pars2 = pars2.replace('B0 = 13.5e-6','B0 = '+str(loc_def_b0)+'e-6')
	pars2 = pars2.replace('B1 = 8.6e-6','B1 = '+str(loc_def_b1)+'e-6')
	pars2 = pars2.replace('h_B1 = 0.6','h_B1 = '+str(loc_def_hb1))
	pars2 = pars2.replace('h_B2 = 5.0','h_B2 = '+str(loc_def_hb2))
	pars2 = pars2.replace('D0 = 2.0e28','D0 = '+str(loc_def_D0)+'e28')
	pars2 = pars2.replace('galaxy_mode = 1','galaxy_mode = '+str(int(def_gmode.get())))
	pars2 = pars2.replace('z1 = 10.0','z1 = '+str(loc_def_Z1))
	pars2 = pars2.replace('h_V = 1.','h_V = '+str(loc_def_hv))
	pars2 = pars2.replace('mode = 1#Adv','mode = '+str(def_mode)+'#Adv')
	pars2 = pars2.replace('adiabatic_losses = -1','adiabatic_losses = '+str((int(def_ablosses))))
	pars2 = pars2.replace('mu_diff = 0.5','mu_diff = '+str(loc_def_mu))
	pars2 = pars2.replace('velocity_field = 1','velocity_field = '+str(def_velfield.get()))
	pars2 = pars2.replace('normalize_intensities = -1','normalize_intensities = '+str(def_norm.get()))
	pars2 = pars2.replace('epsilon = -1','epsilon = '+str(def_epsilon.get()))
	pars2 = pars2.replace('FWHM_effective_beam = 1.2','FWHM_effective_beam = '+str(loc_def_fwhm))

	parfile = open('parameters','w')
	parfile.write(pars2)
	parfile.flush()
	parfile.close()
	#run spinnaker
	if show_spin_out:
#		print 'system'
		os.system('./spectral.x')
	else:
#		subprocess.Popen(['./spectral.x'],stdout=open(os.devnull,'w'))
		os.popen('./spectral.x')
#		print 'popen'
	try:
		simdata = numpy.loadtxt('int.dat')
		sim_int1 = simdata[:,1]
		sim_int2 = simdata[:,2]
		sim_int3 = simdata[:,3]
		sim_int4 = simdata[:,4]
		sim_a12 = simdata[:,5]
		sim_a23 = simdata[:,6]
		sim_a34 = simdata[:,7]
	except:
		return 999

	d2 = numpy.array(data)

	#get newarray
	sim_int1_temp = []	
	sim_int2_temp = []
	sim_int3_temp = []
	sim_int4_temp = []
	sim_a12_temp = []
	sim_a23_temp = []
	sim_a34_temp = []

#	print steps,lim
#	print d2[0,:,0]
	for i in range(0,len(sim_int1)+1,lim):
		sim_int1_temp.append(sim_int1[i])
		sim_int2_temp.append(sim_int2[i])
		sim_int3_temp.append(sim_int3[i])
		sim_int4_temp.append(sim_int4[i])
		sim_a12_temp.append(sim_a12[i])
		sim_a23_temp.append(sim_a23[i])
		sim_a34_temp.append(sim_a34[i])

	sim_int1 = sim_int1_temp
	sim_int2 = sim_int2_temp
	sim_int3 = sim_int3_temp
	sim_int4 = sim_int4_temp
	sim_a12 = sim_a12_temp
	sim_a23 = sim_a23_temp
	sim_a34 = sim_a34_temp

	chi_array = numpy.zeros(shape=d2.shape[0]+d2.shape[0]-1)

##	if haserror:
#		print 'has error'

	if d2.shape[0]==1:
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1)

	if d2.shape[0]==2:
#		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1)
#		chi_array[1] = chisq(d2[1,:,1]/numpy.nanmax(d2[1,:,1]),sim_int2)
#		chi_array[2] = chisq(spi1,sim_a12,spierror)
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1,d2[0,:,2]/numpy.nanmax(d2[0,:,1]))
		chi_array[1] = chisq(d2[1,:,1]/numpy.nanmax(d2[1,:,1]),sim_int2,d2[1,:,2]/numpy.nanmax(d2[1,:,1]))
		chi_array[2] = chisq(spi1,sim_a12,spi1err)

	if d2.shape[0]==3:
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1)
		chi_array[1] = chisq(d2[1,:,1]/numpy.nanmax(d2[1,:,1]),sim_int2)
		chi_array[2] = chisq(d2[2,:,1]/numpy.nanmax(d2[2,:,1]),sim_int3)
		chi_array[3] = chisq(spi1,sim_a12)
		chi_array[4] = chisq(spi2,sim_a23)

	if d2.shape[0]==4:
		chi_array[0] = chisq(d2[0,:,1]/numpy.nanmax(d2[0,:,1]),sim_int1)
		chi_array[1] = chisq(d2[1,:,1]/numpy.nanmax(d2[1,:,1]),sim_int2)
		chi_array[2] = chisq(d2[2,:,1]/numpy.nanmax(d2[2,:,1]),sim_int3)
		chi_array[3] = chisq(d2[3,:,1]/numpy.nanmax(d2[3,:,1]),sim_int4)
		chi_array[4] = chisq(spi1,sim_a12)
		chi_array[5] = chisq(spi2,sim_a23)
		chi_array[6] = chisq(spi3,sim_a34)

	res_chi = get_global_chi()
	if numpy.isnan(res_chi):
		print '999'
		return 999
	else:
		return get_global_chi()

############# setup config frame
configframe = Tk.Frame(master=root)	
configframe.place(relx=0.855,rely=0,relheight=1,relwidth=0.145)
configframe.config(bg=bgcolor)
#setup filename and frequency configarion button/labels
filename1 = Tk.StringVar()
filename1.set('None')
filename2 = Tk.StringVar()
filename2.set('None')
filename3 = Tk.StringVar()
filename3.set('None')
filename4 = Tk.StringVar()
filename4.set('None')

#array to count how many data sets are loaded
filenames = []
filenames.append(0)
filenames.append(0)
filenames.append(0)
filenames.append(0)

data = []


def openfile1():
	global data, filename1, openfile2_button, show_spin_out, debug_out
	fn1 = askopenfilename()
	try:
		d = numpy.loadtxt(fn1,delimiter=';')
		maxval = numpy.nanmax(d[:,1])
		if len(data)==0:
			data.append(d)
		else:
			data[0]=d
	except:
		tkMessageBox.showinfo("Error", "Fileformat not recognised!")
		if debug_out:
			traceback.print_exc()
		return 0

	filename1.set(fn1.split('/')[-1])
	filenames[0]=1
	initial_plot()
	openfile2_button.configure(state=Tk.NORMAL)
	savefig_button.configure(state=Tk.NORMAL)
	parsearch_button.configure(state=Tk.NORMAL)

def openfile2():
	global data, openfile3_button,show_spin_out, debug_out
	fn2 = askopenfilename()
	try:
		d = numpy.loadtxt(fn2,delimiter=';')
		maxval = numpy.nanmax(d[:,1])
		if len(data)==1:
			data.append(d)
		else:
			data[1]=d
	except:
		tkMessageBox.showinfo("Error", "Fileformat not recognised!")
		if debug_out:
			traceback.print_exc()
		return 0

	filename2.set(fn2.split('/')[-1])
	filenames[1]=1
	initial_plot()
	openfile3_button.configure(state=Tk.NORMAL)

def openfile3():
	global data, openfile4_button, show_spin_out, debug_out
	fn3 = askopenfilename()
	try:
		d = numpy.loadtxt(fn3,delimiter=';')
		maxval = numpy.nanmax(d[:,1])
		if len(data)==2:
			data.append(d)
		else:
			data[2]=d

	except:
		tkMessageBox.showinfo("Error", "Fileformat not recognised!")
		if debug_out:
			traceback.print_exc()
		return 0

	filename3.set(fn3.split('/')[-1])
	filenames[2]=1
	initial_plot()
	openfile4_button.configure(state=Tk.NORMAL)

def openfile4():
	global data, show_spin_out, debug_out
	fn4 = askopenfilename()
	try:
		d = numpy.loadtxt(fn4,delimiter=';')
		maxval = numpy.nanmax(d[:,1])
		if len(data)==3:
			data.append(d)
		else:
			data[3]=d

	except:
		tkMessageBox.showinfo("Error", "Fileformat not recognised!")
		if debug_out:
			traceback.print_exc()
		return 0
	filename4.set(fn4.split('/')[-1])
	filenames[3]=1
	initial_plot()

entrywidth=10
openfile1_button = Tk.Button(master=configframe,text='Data 1',command=openfile1)
openfile1_button.grid(row=0,column=0,sticky=Tk.W)
freq1label = Tk.Label(master=configframe,text='Nu 1 Hz:').grid(row=0,column=1,sticky=Tk.W)
openfile1_label = Tk.Label(master=configframe,textvariable=filename1,width=entrywidth).grid(row=1,column=0,columnspan=3,sticky=Tk.W+Tk.E)

openfile2_button = Tk.Button(master=configframe,text='Data 2',command=openfile2,state=Tk.DISABLED)
openfile2_button.grid(row=2,column=0,sticky=Tk.W)
freq2label = Tk.Label(master=configframe,text='Nu 2 Hz:').grid(row=2,column=1,sticky=Tk.W)
openfile2_label = Tk.Label(master=configframe,textvariable=filename2,width=entrywidth).grid(row=3,column=0,columnspan=3,sticky=Tk.W+Tk.E)

openfile3_button = Tk.Button(master=configframe,text='Data 3',command=openfile3,state=Tk.DISABLED)
openfile3_button.grid(row=4,column=0,sticky=Tk.W)
freq3label = Tk.Label(master=configframe,text='Nu 3 Hz:').grid(row=4,column=1,sticky=Tk.W)
openfile3_label = Tk.Label(master=configframe,textvariable=filename3,width=entrywidth).grid(row=5,column=0,columnspan=3,sticky=Tk.W+Tk.E)

openfile4_button = Tk.Button(master=configframe,text='Data 4',command=openfile4,state=Tk.DISABLED)
openfile4_button.grid(row=6,column=0,sticky=Tk.W)
freq4label = Tk.Label(master=configframe,text='Nu 4 Hz:').grid(row=6,column=1,sticky=Tk.W)
openfile4_label = Tk.Label(master=configframe,textvariable=filename4,width=entrywidth).grid(row=7,column=0,columnspan=3,sticky=Tk.W+Tk.E)

frequency1field = Tk.Entry(master=configframe,width=entrywidth)
frequency1field.grid(row=0,column=2,sticky=Tk.W)
frequency1field.insert(Tk.END,def_nu1)
frequency2field = Tk.Entry(master=configframe,width=entrywidth)
frequency2field.grid(row=2,column=2,sticky=Tk.W)
frequency2field.insert(Tk.END,def_nu2)
frequency3field = Tk.Entry(master=configframe,width=entrywidth)
frequency3field.grid(row=4,column=2,sticky=Tk.W)
frequency3field.insert(Tk.END,def_nu3)
frequency4field = Tk.Entry(master=configframe,width=entrywidth)
frequency4field.grid(row=6,column=2,sticky=Tk.W)
frequency4field.insert(Tk.END,def_nu4)

fwhm_var = Tk.DoubleVar()
fwhm_var.set(def_fwhm)
FWHMfield = Tk.Entry(master=configframe,width=entrywidth,textvariable=fwhm_var)
FWHMfield.grid(row=8,column=2,sticky=Tk.W)
Tk.Label(master=configframe,text='Effective Beam FWHM').grid(row=8,column=0,columnspan=2,sticky=Tk.W+Tk.E)

Tk.Frame(master=configframe,relief=Tk.SOLID,height=10,bg='red',width=entrywidth).grid(row=9,columnspan=3,sticky=Tk.W+Tk.E)
Tk.Label(master=configframe,text='Parameter Settings',relief=Tk.SUNKEN).grid(row=10,column=0,columnspan=3,sticky=Tk.W+Tk.E)
Tk.Label(master=configframe,text='Parameter').grid(row=11,column=0,columnspan=2,sticky=Tk.W+Tk.E)
Tk.Label(master=configframe,text='Delta').grid(row=11,column=2,columnspan=1,sticky=Tk.W+Tk.E)

##########################################################

def clicktab(event):
	clicked_tab = tabs.tk.call(tabs._w,"identify","tab",event.x,event.y)
	global def_modestr,def_mode
	if clicked_tab==0:
		def_modestr.set('Advection')
		def_mode = 1
	else:
		def_modestr.set('Diffusion')
		def_mode = 2
	run()
	plot_dummylegend()

############# setup parameter buttons
global tabs
tabs = ttk.Notebook(master=configframe)
tabs.bind("<Button-1>",clicktab)
tab_diff = ttk.Frame(tabs)
tab_adv = ttk.Frame(tabs)
tabs.add(tab_adv,text='Advection')
tabs.add(tab_diff,text='Diffusion')

if def_modestr.get()=='Advection':
	tabs.select(tab_adv)
else:
	tabs.select(tab_diff)

tabs.grid(column=0,row=21,rowspan=4,columnspan=3,sticky=Tk.W+Tk.E)


	
def parp(var,delta):
	exec('{0}+={1}'.format(var,delta),globals())
	run()
	plot_dummylegend()

def parm(var,delta):
	exec('{0}-={1}'.format(var,delta),globals())
	run()
	plot_dummylegend()	


def setdummval(val):
	#val no usage, just here so the function executes
	run()
	plot_dummylegend()

def_imgtype = Tk.StringVar()
def_imgtype.set('jpg')

def save_plot():
	global plotfig
	global def_imgtype
	number = 1
	while os.path.exists('plot_'+str(number)+'.'+def_imgtype.get().replace('\n','')):
		number+=1
	plt.savefig('plot_'+str(number)+'.'+def_imgtype.get())
	saveplotpars('plot_'+str(number)+'.'+def_imgtype.get()+'.pars')

def dummyrefresh(a):
	if sum(filenames)!=0:
		initial_plot()

#parameter definition for parameter space search
gamma_check_var = Tk.IntVar()
b0_check_var = Tk.IntVar()
b1_check_var = Tk.IntVar()
hb0_check_var = Tk.IntVar()
hb1_check_var = Tk.IntVar()
z1_check_var = Tk.IntVar()
#advection
v0_check_var = Tk.IntVar()
hv_check_var = Tk.IntVar()
#diffusion
d0_check_var = Tk.IntVar()
mu_check_var = Tk.IntVar()
searchres= Tk.StringVar()
searchres.set('coarse')

def search():
	searchpars = ['gamma_check_var','b0_check_var','b1_check_var','hb0_check_var','hb1_check_var','z1_check_var','v0_check_var','hv_check_var','d0_check_var','mu_check_var']
	pararray = ['def_cri', 'def_b0', 'def_b1', 'def_hb1', 'def_hb2', 'def_Z1', 'def_v0', 'def_hv', 'def_D0', 'def_mu']
	labelarray = ['def_cri', 'def_b0', 'def_b1', 'def_hb1', 'def_hb2', 'def_Z1', 'def_v0', 'def_hv', 'def_D0', 'def_mu']
	global tosearch
	tosearch = []
	for i in range(len(searchpars)):
		global b
		b=0
#		exec('print {0}.get()'.format(searchpars[i]))
		exec('if {0}.get()==1: b=1'.format(searchpars[i]),globals())
		if b==1:
			tosearch.append(pararray[i])
	#get current values
	currentvalues = []
	for i in tosearch:
		exec('a={0}'.format(i),globals())
		currentvalues.append(a)
		
	#determine search resolution for coarse mode
	searches = []

	searchdict = {}
	searchdict["def_cri"]=[0.1,5,0.5]
	searchdict["def_b0"]=[0.1,30,5]
	searchdict["def_b1"]=[0.1,30,5]
	searchdict["def_hb1"]=[0.1,15,2]
	searchdict["def_hb2"]=[0.1,15,2]
	searchdict["def_Z1"]=[0.1,20,2]
	searchdict["def_v0"]=[10,1000,100]
	searchdict["def_hv"]=[0.1,5,0.5]
	searchdict["def_D0"]=[1,5,0.5]
	searchdict["def_mu"]=[0.1,1,0.1]

	searchdict2 = {}
	searchdict2["def_cri"]=[0.1,5,1]
	searchdict2["def_b0"]=[0.1,30,10]
	searchdict2["def_b1"]=[0.1,30,10]
	searchdict2["def_hb1"]=[0.1,15,5]
	searchdict2["def_hb2"]=[0.1,15,5]
	searchdict2["def_Z1"]=[0.1,20,4]
	searchdict2["def_v0"]=[10,1000,300]
	searchdict2["def_hv"]=[0.1,5,1]
	searchdict2["def_D0"]=[1,5,1]
	searchdict2["def_mu"]=[0.1,1,0.1]

	for i in range(len(currentvalues)):
		if searchres.get()=='very coarse':
			searches.append((tosearch[i],searchdict2[tosearch[i]][0],searchdict2[tosearch[i]][1],searchdict2[tosearch[i]][2]))

		if searchres.get()=='coarse':
			searches.append((tosearch[i],searchdict[tosearch[i]][0],searchdict[tosearch[i]][1],searchdict[tosearch[i]][2]))

		if searchres.get()=='medium':
			delta = currentvalues[i]*0.50
			minval = currentvalues[i]-delta
			if minval<0:
				minval=delta
			maxval = currentvalues[i]+delta
			stepsize = 2*delta*0.1		#take 10% of total range as stepsize
			searches.append((tosearch[i],minval,maxval,stepsize))
		
		if searchres.get()=='fine':
			delta = currentvalues[i]*0.20
			minval = currentvalues[i]-delta
			if minval<0:
				minval=delta

			maxval = currentvalues[i]+delta
			stepsize = 2*delta*0.05		#take 5% of total range as stepsize
			searches.append((tosearch[i],minval,maxval,stepsize))

	varstochange = []
	mins = []
	maxs = []
#	steps = []
	vals = numpy
	for i in searches:
		varstochange.append(i[0])


	totaldict = {}
	for i in range(len(varstochange)):
		ar = []
		for y in numpy.arange(searches[i][1],searches[i][2]+searches[i][3],searches[i][3]):
			ar.append(y)
		totaldict["{0}".format(i)]=ar


	#sort dictionary
	ordereddict = collections.OrderedDict(sorted(totaldict.items()))
	combinations = list(itertools.product(*ordereddict.values()))
	resultdict = {}
	itercount = 0
	global progress
	global bestresult
	global bestresultchi
	#stuff for ETA calculation
	timearray = []
	Ncombinations = len(combinations)
	for i in combinations:
		t0 = time.time()
		result = run_noplot(varstochange,i)
	#	print i,result
		resultdict[str(i)]=result
		itercount+=1
		progress.set(str(round(100*itercount/Ncombinations,2))+'%')
		root.update_idletasks()
		t1 = time.time()
		timearray.append(t1-t0)
		avg_time = numpy.mean(timearray)
		eta_time = avg_time*(Ncombinations-itercount)
		h = int(eta_time/3600)
		m = int((eta_time-h*3600)/60)
		s = int(eta_time-h*3600-m*60)
		eta_label.set(str(int(h)).zfill(2)+':'+str(int(m)).zfill(2)+':'+str(s).zfill(2))

	#reset labels
	for v in pararray:
		exec('{0}_result_1_var.set("")'.format(v))
		exec('{0}_result_2_var.set("")'.format(v))
		exec('{0}_result_3_var.set("")'.format(v))

	orderedresult = sorted(resultdict.items(), key=operator.itemgetter(1))
	global best_results_array
	best_results_array = numpy.array(currentvalues)
	for v in range(len(tosearch)):
		exec('{0}={1}'.format('best_results_array['+str(v)+']',float(orderedresult[0][0][1:-1].split(',')[v]),globals()))
		exec('{0}_result_1_var.set({1})'.format(tosearch[v],float(orderedresult[0][0][1:-1].split(',')[v])))
		exec('{0}_result_2_var.set({1})'.format(tosearch[v],float(orderedresult[1][0][1:-1].split(',')[v])))
		exec('{0}_result_3_var.set({1})'.format(tosearch[v],float(orderedresult[2][0][1:-1].split(',')[v])))

	chi_result_1_label.set(round(orderedresult[0][1],5))
	chi_result_2_label.set(round(orderedresult[1][1],5))
	chi_result_3_label.set(round(orderedresult[2][1],5))
		
global progress
progress = Tk.StringVar()
global chi_result_1_label
chi_result_1_label=Tk.StringVar()
global chi_result_2_label
chi_result_2_label=Tk.StringVar()
global chi_result_3_label
chi_result_3_label=Tk.StringVar()

global def_cri_result_1_var
def_cri_result_1_var = Tk.StringVar()
global def_cri_result_2_var
def_cri_result_2_var = Tk.StringVar()
global def_cri_result_3_var
def_cri_result_3_var = Tk.StringVar()

global def_b0_result_1_var
def_b0_result_1_var = Tk.StringVar()
global def_b0_result_2_var
def_b0_result_2_var = Tk.StringVar()
global def_b0_result_3_var
def_b0_result_3_var = Tk.StringVar()

global def_b1_result_1_var
def_b1_result_1_var = Tk.StringVar()
global def_b1_result_2_var
def_b1_result_2_var = Tk.StringVar()
global def_b1_result_3_var
def_b1_result_3_var = Tk.StringVar()

global def_hb1_result_1_var
def_hb1_result_1_var = Tk.StringVar()
global def_hb1_result_2_var
def_hb1_result_2_var = Tk.StringVar()
global def_hb1_result_3_var
def_hb1_result_3_var = Tk.StringVar()

global def_hb2_result_1_var
def_hb2_result_1_var = Tk.StringVar()
global def_hb2_result_2_var
def_hb2_result_2_var = Tk.StringVar()
global def_hb2_result_3_var
def_hb2_result_3_var = Tk.StringVar()

global def_Z1_result_1_var
def_Z1_result_1_var = Tk.StringVar()
global def_Z1_result_2_var
def_Z1_result_2_var = Tk.StringVar()
global def_Z1_result_3_var
def_Z1_result_3_var = Tk.StringVar()

global def_v0_result_1_var
def_v0_result_1_var = Tk.StringVar()
global def_v0_result_2_var
def_v0_result_2_var = Tk.StringVar()
global def_v0_result_3_var
def_v0_result_3_var = Tk.StringVar()

global def_hv_result_1_var
def_hv_result_1_var = Tk.StringVar()
global def_hv_result_2_var
def_hv_result_2_var = Tk.StringVar()
global def_hv_result_3_var
def_hv_result_3_var = Tk.StringVar()

global def_D0_result_1_var
def_D0_result_1_var = Tk.StringVar()
global def_D0_result_2_var
def_D0_result_2_var = Tk.StringVar()
global def_D0_result_3_var
def_D0_result_3_var = Tk.StringVar()

global def_mu_result_1_var
def_mu_result_1_var = Tk.StringVar()
global def_mu_result_2_var
def_mu_result_2_var = Tk.StringVar()
global def_mu_result_3_var
def_mu_result_3_var = Tk.StringVar()

global eta_label
eta_label = Tk.StringVar()

global best_results_array
best_results_array = []
#global search_dif_button
#global search_adv_button

def setsearchmode_adv():
	global d0_check_var
	global mu_check_var
	d0_check_var.set(0)
	mu_check_var.set(0)
	global def_mode
	def_mode = 1
	global tabs
	tabs.select(tab_adv)
	initial_plot()
	global root
	root.update_idletasks()

def setsearchmode_dif():
	global v0_check_var
	global hv_check_var
	v0_check_var.set(0)
	hv_check_var.set(0)
	global def_mode
	def_mode = 2
	global tabs
	tabs.select(tab_diff)
	initial_plot()
	global root
	root.update_idletasks()

global applyvalues
applyvalues = 0

def apply_vals():
	searchpars = ['gamma_check_var','b0_check_var','b1_check_var','hb0_check_var','hb1_check_var','z1_check_var','v0_check_var','hv_check_var','d0_check_var','mu_check_var']
	pararray = ['def_cri', 'def_b0', 'def_b1', 'def_hb1', 'def_hb2', 'def_Z1', 'def_v0', 'def_hv', 'def_D0', 'def_mu']
	labelarray = ['def_cri', 'def_b0', 'def_b1', 'def_hb1', 'def_hb2', 'def_Z1', 'def_v0', 'def_hv', 'def_D0', 'def_mu']
	for i in range(len(tosearch)):
		exec('{0}={1}'.format(tosearch[i],best_results_array[i],globals()))
	global applyvalues
	applyvalues=1
	initial_plot()

def open_parsearch_window():
	global search_dif_button
	global search_adv_button

	parsearchframe = Tk.Toplevel(root)
	parsearchframe.wm_title('Parameter search')
	Tk.Label(master=parsearchframe,text='Choose parameters to iterate over').grid(row=0,column=0,columnspan=3)
	Tk.Label(master=parsearchframe,text='Parameter').grid(row=1,column=0)
	Tk.Label(master=parsearchframe,text='Best Value').grid(row=1,column=1)
	Tk.Label(master=parsearchframe,text='2nd Best').grid(row=1,column=2)
	Tk.Label(master=parsearchframe,text='3rd Best').grid(row=1,column=3)
	gamma_box = Tk.Checkbutton(parsearchframe,text=u"\N{GREEK SMALL LETTER GAMMA}", var=gamma_check_var).grid(row=2,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_cri_result_1_var).grid(row=2,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_cri_result_2_var).grid(row=2,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_cri_result_3_var).grid(row=2,column=3,stick=Tk.W)
	b0_box = Tk.Checkbutton(master=parsearchframe,text=u"B\u2080", var=b0_check_var).grid(row=3,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_b0_result_1_var).grid(row=3,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_b0_result_2_var).grid(row=3,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_b0_result_3_var).grid(row=3,column=3,stick=Tk.W)
	b1_box = Tk.Checkbutton(master=parsearchframe,text=u"B\u2081", var=b1_check_var).grid(row=4,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_b1_result_1_var).grid(row=4,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_b1_result_2_var).grid(row=4,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_b1_result_3_var).grid(row=4,column=3,stick=Tk.W)
	hb0_box = Tk.Checkbutton(master=parsearchframe,text=u"hB\u2081", var=hb0_check_var).grid(row=5,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hb1_result_1_var).grid(row=5,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hb1_result_2_var).grid(row=5,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hb1_result_3_var).grid(row=5,column=3,stick=Tk.W)
	hb1_box = Tk.Checkbutton(master=parsearchframe,text=u"hB\u2082", var=hb1_check_var).grid(row=6,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hb2_result_1_var).grid(row=6,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hb2_result_2_var).grid(row=6,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hb2_result_3_var).grid(row=6,column=3,stick=Tk.W)
	z1_box = Tk.Checkbutton(master=parsearchframe,text=u'Z\u2081', var=z1_check_var).grid(row=7,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_Z1_result_1_var).grid(row=7,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_Z1_result_2_var).grid(row=7,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_Z1_result_3_var).grid(row=7,column=3,stick=Tk.W)
	v0_box = Tk.Checkbutton(master=parsearchframe,text='v0', var=v0_check_var).grid(row=8,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_v0_result_1_var).grid(row=8,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_v0_result_2_var).grid(row=8,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_v0_result_3_var).grid(row=8,column=3,stick=Tk.W)
	hv_box = Tk.Checkbutton(master=parsearchframe,text='h_v', var=hv_check_var).grid(row=9,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hv_result_1_var).grid(row=9,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hv_result_2_var).grid(row=9,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_hv_result_3_var).grid(row=9,column=3,stick=Tk.W)
	d0_box = Tk.Checkbutton(master=parsearchframe,text='D0', var=d0_check_var).grid(row=10,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_D0_result_1_var).grid(row=10,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_D0_result_2_var).grid(row=10,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_D0_result_3_var).grid(row=10,column=3,stick=Tk.W)
	mu_box = Tk.Checkbutton(master=parsearchframe,text=u'\N{GREEK SMALL LETTER MU}', var=mu_check_var).grid(row=11,column=0,sticky=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_mu_result_1_var).grid(row=11,column=1,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_mu_result_2_var).grid(row=11,column=2,stick=Tk.W)
	Tk.Label(parsearchframe,textvariable=def_mu_result_3_var).grid(row=11,column=3,stick=Tk.W)

	Tk.Label(master=parsearchframe,text=u'\u03C7^2').grid(row=12,column=0,sticky=Tk.W)
	Tk.Label(master=parsearchframe,textvariable=chi_result_1_label).grid(row=12,column=1,sticky=Tk.W)
	Tk.Label(master=parsearchframe,textvariable=chi_result_2_label).grid(row=12,column=2,sticky=Tk.W)
	Tk.Label(master=parsearchframe,textvariable=chi_result_3_label).grid(row=12,column=3,sticky=Tk.W)
	Tk.Label(parsearchframe,text='Search resolution').grid(row=13,column=0,sticky=Tk.W)
	Tk.OptionMenu(parsearchframe, searchres, 'very coarse', 'coarse','medium','fine').grid(row=13,column=1)
	search_adv_button=Tk.Button(parsearchframe,text='Advection',command=setsearchmode_adv).grid(row=13,column=2)
	search_dif_button=Tk.Button(parsearchframe,text='Diffusion',command=setsearchmode_dif).grid(row=13,column=3)
	searchbutton = Tk.Button(parsearchframe,text='Search!',command=search).grid(row=14,column=0,columnspan=3)
	searchbutton = Tk.Button(parsearchframe,text='Apply results',command=apply_vals).grid(row=14,column=3,columnspan=1)
	Tk.Label(master=parsearchframe,text='Progress:').grid(row=15,column=0,sticky=Tk.W)
	Tk.Label(master=parsearchframe,textvariable=progress).grid(row=15,column=1,sticky=Tk.W)
	Tk.Label(master=parsearchframe,text='Time remaining:').grid(row=15,column=2,sticky=Tk.W)
	Tk.Label(master=parsearchframe,textvariable=eta_label).grid(row=15,column=3,sticky=Tk.W)




rownum=12
button_crip = Tk.Button(configframe, text=u"\N{GREEK SMALL LETTER GAMMA}+", command=lambda: parp('def_cri','dcri.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_crim = Tk.Button(configframe, text=u"\N{GREEK SMALL LETTER GAMMA}-", command=lambda: parm('def_cri','dcri.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)
dcri = Tk.DoubleVar()
dcri.set(0.1)
delta_cri = Tk.Entry(configframe,width=10,textvariable=dcri).grid(row=rownum,column=2)

rownum=13
button_b0p = Tk.Button(configframe, text=u"B\u2080+", command=lambda: parp('def_b0','db0.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_b0m = Tk.Button(configframe, text=u"B\u2080-", command=lambda: parm('def_b0','db0.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)
db0 = Tk.DoubleVar()
db0.set(0.5)
delta_b0 = Tk.Entry(configframe,width=10,textvariable=db0).grid(row=rownum,column=2)

rownum=14
button_b1p = Tk.Button(configframe, text=u"B\u2081+", command=lambda: parp('def_b1','db1.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_b1m = Tk.Button(configframe, text=u"B\u2081-", command=lambda: parm('def_b1','db1.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)
db1 = Tk.DoubleVar()
db1.set(0.5)
delta_b1 = Tk.Entry(configframe,width=10,textvariable=db1).grid(row=rownum,column=2)

rownum=15
button_hb0p = Tk.Button(configframe, text=u"hB\u2081"+"+", command=lambda: parp('def_hb1','dhb0.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_hb0m = Tk.Button(configframe, text=u"hB\u2081"+"-", command=lambda: parm('def_hb1','dhb0.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)
dhb0 = Tk.DoubleVar()
dhb0.set(0.5)
delta_b0 = Tk.Entry(configframe,width=10,textvariable=dhb0).grid(row=rownum,column=2)

rownum=16
button_hb1p = Tk.Button(configframe, text=u"hB\u2082+", command=lambda: parp('def_hb2','dhb1.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_hb1m = Tk.Button(configframe, text=u"hB\u2082-", command=lambda: parm('def_hb2','dhb1.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)
dhb1 = Tk.DoubleVar()
dhb1.set(0.5)
delta_b1 = Tk.Entry(configframe,width=10,textvariable=dhb1).grid(row=rownum,column=2)

rownum=17
button_z1p = Tk.Button(configframe, text=u'Z\u2081+', command=lambda: parp('def_Z1','dZ.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_z1m = Tk.Button(configframe, text=u'Z\u2081-', command=lambda: parm('def_Z1','dZ.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)
dZ = Tk.DoubleVar()
dZ.set(1)
delta_Z = Tk.Entry(configframe,width=10,textvariable=dZ).grid(row=rownum,column=2)

rownum=18
Tk.Label(configframe,text='Galaxy Mode').grid(row=rownum,column=0,columnspan=2,sticky=Tk.W+Tk.E)
option = Tk.OptionMenu(configframe, def_gmode, 1, -1,command=dummyrefresh).grid(row=rownum,column=2)

rownum=19
Tk.Label(configframe,text='Epsilon').grid(row=rownum,column=0,columnspan=2,sticky=Tk.W+Tk.E)
option = Tk.OptionMenu(configframe, def_epsilon, 1, -1,command=dummyrefresh).grid(row=rownum,column=2)

rownum=20
Tk.Label(configframe,text='Normalize').grid(row=rownum,column=0,columnspan=2,sticky=Tk.W+Tk.E)
option = Tk.OptionMenu(configframe, def_norm, 1, -1,command=dummyrefresh).grid(row=rownum,column=2)

def_scale = Tk.StringVar()
def_scale.set('Linear')
rownum=25
scalelabel = Tk.Label(configframe,text='Scaling: ').grid(row=rownum,column=0)
option = Tk.OptionMenu(configframe, def_scale, 'Linear', 'Log',command=dummyrefresh).grid(row=rownum,column=1,columnspan=2)

rownum=26
refresh_button = Tk.Button(configframe, text='Refresh Plot', command=initial_plot)
refresh_button.grid(row=rownum,column=0,columnspan=3,sticky=Tk.W+Tk.E)

rownum=27
savefig_button = Tk.Button(configframe, text='Save Plot', command=save_plot,state=Tk.DISABLED)
savefig_button.grid(row=rownum,column=0,columnspan=2,sticky=Tk.W+Tk.E)
option = Tk.OptionMenu(configframe, def_imgtype, 'jpg', 'png', 'eps','svg','pdf').grid(row=rownum,column=2)

rownum=28
parsearch_button = Tk.Button(configframe, text='Parameter Search', command=open_parsearch_window,state=Tk.DISABLED)
parsearch_button.grid(row=rownum,column=0,columnspan=2,sticky=Tk.W+Tk.E)

#label for x,y data
distlabelvar = Tk.StringVar()
distlabelvar.set('')
rownum=29
distlabel=Tk.Label(configframe,textvariable=distlabelvar)
distlabel.grid(row=rownum,column=0,columnspan=3,sticky=Tk.W+Tk.E)

#advection only options
rownum=0
dv0 = Tk.DoubleVar()
dv0.set(10)
delta_v0 = Tk.Entry(tab_adv,width=10,textvariable=dv0).grid(row=rownum,column=2)
button_v0p = Tk.Button(tab_adv, text='v0+', command=lambda: parp('def_v0','dv0.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_v0m = Tk.Button(tab_adv, text='v0-', command=lambda: parm('def_v0','dv0.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)

rownum=1
button_hvp = Tk.Button(tab_adv, text='h_v+', command=lambda: parp('def_hv','dhv.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_hvm = Tk.Button(tab_adv, text='h_v-', command=lambda: parm('def_hv','dhv.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)
dhv = Tk.DoubleVar()
dhv.set(1)
delta_hv = Tk.Entry(tab_adv,width=10,textvariable=dhv).grid(row=rownum,column=2)

rownum=2
Tk.Label(tab_adv,text='Adiabatic Losses').grid(row=rownum,column=0,columnspan=2,sticky=Tk.W+Tk.E)
option = Tk.OptionMenu(tab_adv, def_ablossesstr, "Yes", "No" ,command=setdummval).grid(row=rownum,column=2)

rownum=3
Tk.Label(tab_adv,text='Velocity Field').grid(row=rownum,column=0,columnspan=2,sticky=Tk.W+Tk.E)
option = Tk.OptionMenu(tab_adv, def_velfield, -1, 0, 1, 2,command=setdummval).grid(row=rownum,column=2)


#diffusion only options
rownum = 0
d0 = Tk.DoubleVar()
d0.set(0.1)
delta_d0 = Tk.Entry(tab_diff,width=10,textvariable=d0).grid(row=rownum,column=2)
button_d0p = Tk.Button(tab_diff, text='D0+', command=lambda: parp('def_D0','d0.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_d0m = Tk.Button(tab_diff, text='D0-', command=lambda: parm('def_D0','d0.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)

rownum = 1
mu = Tk.DoubleVar()
mu.set(0.1)
delta_mu = Tk.Entry(tab_diff,width=10,textvariable=mu).grid(row=rownum,column=2)
button_mup = Tk.Button(tab_diff, text=u'\N{GREEK SMALL LETTER MU}+', command=lambda: parp('def_mu','mu.get()')).grid(row=rownum,column=0,sticky=Tk.W+Tk.E)
button_mum = Tk.Button(tab_diff, text=u'\N{GREEK SMALL LETTER MU}-', command=lambda: parm('def_mu','mu.get()')).grid(row=rownum,column=1,sticky=Tk.W+Tk.E)

#run()
root.protocol("WM_DELETE_WINDOW", quitit)
Tk.mainloop()

