import os,sys,numpy
import Tkinter as Tk

class main():
	#define root window
	root = Tk.Tk()
	bgcolor = '#BFBFBF'
	Tk.Grid.columnconfigure(root,0,weight=1)
	Tk.Grid.rowconfigure(root,0,weight=1)
	w, h = root.winfo_screenwidth(), root.winfo_screenheight()
	root.geometry("%dx%d+0+0" % (w*0.4, h*0.75))
	root.wm_title("Spinteractive v2")
	root.config(bg=bgcolor)
	root.resizable(width=True,height=True)

	def quitit():
		pars = [self.def_cri,self.def_v0,self.def_hb1,self.def_hb2,self.def_b0,self.def_b1,self.def_Z1,self.def_hv,self.def_gmode.get(),self.def_D0,self.def_mu,self.def_modestr.get(),self.def_ablossesstr.get(),self.def_velfield.get()]
		numpy.savetxt('pars.last',pars,fmt='%s')
		sys.exit()


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
		print 'Parameters from last file loaded succesfully'
	except:
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
	##########################################################



	root.protocol("WM_DELETE_WINDOW", quitit)
	Tk.mainloop()

if __name__=="__main__":
	main()

