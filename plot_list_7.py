# Plots the dr7 spectrum and err from QSO in the file list.txt
# Uses Liam's code
# Thu 3 Nov
# Author: mjt91

from get_spectra import get_sdss_dr7_spec
#from get_spectra import get_boss_dr12_spec
import matplotlib.pyplot as plt

with open(r"list.txt", 'r') as infile:
    name = infile.read().strip()
    w7,dw7,f7,e7, mask = get_sdss_dr7_spec(name)
    print(name)

#de-redshift
z = float(input('Enter redshift z: '))
w7[:] = [x/(1+z) for x in w7]



#f = open('list.txt')
#name = (line.strip() for line in f)
#This might be needed later for multiple lines
#name = f

#plt.plot(w7,f7*mask, w7, e7)
plt.plot(w7,f7, w7, e7)
plt.xlabel('wavelength/angstroms')
plt.ylabel('flux, error')
plt.title(str(name)+', z = '+str(z))
plt.show()


