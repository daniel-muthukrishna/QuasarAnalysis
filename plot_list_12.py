# Plots the dr12 spectrum and err from QSO in the file list.txt
# Uses Liam's code
#  Thu 3 Nov 2016
# Author: mjt91

#from get_spectra import get_sdss_dr7_spec
from get_spectra import get_boss_dr12_spec
import matplotlib.pyplot as plt

with open(r"list.txt", 'r') as infile:
    name = infile.read().strip()
    w12,dw12,f12,e12, mask = get_boss_dr12_spec(name)
    print(name)

#De-redshift
z = float(input('Enter Redshift z: '))
w12[:] = [(x/(1+z)) for x in w12]


#f = open('list.txt')
#name = (line.strip() for line in f)
#This might be needed later for multiple lines
#name = f

plt.plot(w12,f12,w12,e12)
#plt.plot(w12,f12*mask, w12, e12)
plt.xlabel('wavelength / angstroms')
plt.ylabel('flux, error')
plt.title(str(name)+', z = '+str(z))
plt.show()

