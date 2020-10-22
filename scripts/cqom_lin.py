import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss

# axis for mechanical mode
Omega_1     = np.linspace(0, 2, 301)
# axis for optical mode and laser
Omega_2     = np.linspace(1e8 - 1.25, 1e8 + 1.25, 501)
# mechanical mode
Omega_m     = ss.norm(1, 0.1).pdf(Omega_1)
Omega_m     = Omega_m/max(Omega_m)*0.5
# mechanical mode on red detuning
Omega_m_r   = ss.norm(1, 0.25).pdf(Omega_1)
Omega_m_r   = Omega_m_r/max(Omega_m_r)*0.1
# mechanical mode on blue detuning
Omega_m_b   = ss.norm(1, 0.025).pdf(Omega_1)
Omega_m_b   = Omega_m_b/max(Omega_m_b)
# laser mode for red detuning
Omega_2_r       = Omega_2[0:249]
Omega_l_r       = np.zeros(len(Omega_2_r)) 
Omega_l_r[50]   = 1
# laser mode for blue detuning
Omega_2_b       = Omega_2[250:501]
Omega_l_b       = np.zeros(len(Omega_2_b)) 
Omega_l_b[200]  = 1
# optical mode
Omega_o     = ss.norm(1e8, 0.25).pdf(Omega_2)
Omega_o     = Omega_o/max(Omega_o)

# colors for red and blue detunings
color_r     = '#d94732aa'
color_b     = '#4a77c2aa'

# figure layout
fig     = plt.figure(constrained_layout=True)
gs      = fig.add_gridspec(1, 6)
ax_1    = fig.add_subplot(gs[0, :2]) 
ax_2    = fig.add_subplot(gs[0, 2:], sharey=ax_1)

# first plot
ax_1.plot(Omega_1, Omega_m, c='k', alpha=0.5)
ax_1.plot(Omega_1, Omega_m_r, c=color_r)
ax_1.plot(Omega_1, Omega_m_b, c=color_b)
ax_1.tick_params(labelsize=16)
ax_1.set_xticks(np.arange(0, 3, 1))
# ax_1.axvspan(Omega_1[132], Omega_1[168], color='k', alpha=0.25)
# second plot
ax_2.plot()
ax_2.plot(Omega_2, Omega_o, c='k', alpha=0.5)
ax_2.plot(Omega_2_r, Omega_l_r, c=color_r)
ax_2.plot(Omega_2_b, Omega_l_b, c=color_b)
ax_2.axvspan(Omega_2[50], Omega_2[249], color=color_r, alpha=0.25)
ax_2.axvspan(Omega_2[251], Omega_2[450], color=color_b, alpha=0.25)
# ax_2.axvspan(Omega_2[192], Omega_2[308], color='k', alpha=0.25)
ax_2.tick_params(labelsize=16)
ax_2.set_xticks(np.arange(1e8 - 1, 1e8 + 2, 1))
ax_2.ticklabel_format(style='plain', useOffset=False, useMathText=True)

ax_1.spines['top'].set_visible(False)
ax_1.spines['right'].set_visible(False)
ax_1.yaxis.tick_left()
ax_2.spines['top'].set_visible(False)
ax_2.spines['left'].set_visible(False)
ax_2.yaxis.tick_right()
plt.show()