import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import socket
import pandas as pd
import sys, os, glob, pickle

if (socket.gethostname() == "quirm.math.grinnell.edu"):
    prefix = '/home/christenc/Data/Sims/'
    outprefix = '/home/christenc/Figures/marvel/'
    dataprefix = '/home/christenc/Code/Datafiles/'
else:
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    dataprefix = '/home/christensen/Code/Datafiles/'

#Sandra
tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
objs_sandra = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_sandra.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_sandra = pd.DataFrame(objs_sandra)

#Plots stellar mass halo mass relation
plt.figure(1)
sandra_plt = plt.scatter(objs_pd_sandra['mvir'],objs_pd_sandra['mstar'])
plt.scatter(objs_pd_sandra['mvir'][objs_pd_sandra['mgas'] != 0],objs_pd_sandra['mstar'][objs_pd_sandra['mgas'] != 0]) #Plot only those halos with non-zero mass of gas
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Virial Mass [M$_\odot$]')
plt.ylabel(r'Stellar Mass [M$_\odot$]')
