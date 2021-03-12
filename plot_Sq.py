#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
fig = plt.figure(figsize=(10,10))

nqx=21
nqy=21

for i in range(21):
    ax=plt.subplot(5,5,i+1)
    cp="{:.1f}".format(0.2*i)
    qx,qy,Sq = np.loadtxt('CP='+cp+'/Sq.dat', unpack=True)
    S1 = Sq.reshape(nqx,nqy)
    S1 += S1.transpose()
    S1 += np.rot90(S1)
    cax = ax.imshow(S1, extent=(-np.pi,np.pi,-np.pi,np.pi), aspect = 'equal',origin="lower", cmap='PuRd', vmin=0)#,vmax=30)
    ax.set_xticks([-np.pi,0,np.pi])
    ax.set_yticks([-np.pi,0,np.pi])
    ax.xaxis.set_ticklabels([r'$-\pi$',0,r'$\pi$'])
    ax.yaxis.set_ticklabels([r'$-\pi$',0,r'$\pi$'])

plt.savefig('Sq.pdf', format='pdf', transparent=1, bbox_inches='tight', pad_inches = 0.02)