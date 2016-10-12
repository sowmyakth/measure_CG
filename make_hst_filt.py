import numpy as np

def main():
	filts = ['f606w', 'f814w']
	filts = ['f775w']
	for filt in filts:
		file_name = filt + '.UVIS1.tab'
		f = np.loadtxt(file_name).T
		w_A = f[1]
		t = f[2]
		header  = "{0}-band.\n\
		File taken from http://www.stsci.edu/hst/wfc3/ins_performance/throughputs/Throughput_Tables\n\
		Wavelength(nm)  Throughput(0-1)".format(filt)
		file_name = 'HST_'+ filt + '.dat'
		data = np.array([w_A/10.,t])
		np.savetxt(file_name, data.T, header=header) 



if __name__ == '__main__':
	main()