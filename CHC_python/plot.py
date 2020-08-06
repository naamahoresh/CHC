import numpy as np
import matplotlib.pyplot as plt


vGA_p1_11_inst= np.asarray([64,64,48,64,51,48,64,64,47,64,64])
vGA_p1_1_inst= np.asarray([64,64,64,64,64,64,64,64,64,64,64])

vGA_p2_11_inst= np.asarray([64,14,12,64,14,14,17,64,15,64,64])
vGA_p2_1_inst= np.asarray([64,64,64,64,64,64,64,64,64,64,64])


data_to_plot_p1 = [vGA_p1_1_inst, vGA_p1_11_inst]
data_to_plot_p2 = [vGA_p2_1_inst, vGA_p2_11_inst]

#
fig = plt.figure(1, figsize=(9, 6))
ax = fig.add_subplot(111)

bp = ax.boxplot(data_to_plot_p2,patch_artist=True,widths=(0.5, 0.5))
colors = ["#C8C2DE","#FFF3DD"]

plt.xticks([1, 2], ['(30,30) vGA: 1 instance', '(30,30) vGA: 11 instances'])
for patch, color_tmp in zip(bp['boxes'], colors):
    patch.set_facecolor(color_tmp)

for median in bp['medians']:
    median.set(color='#000000', linewidth=2)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)

ax.set_xlabel('Algorithm', fontsize=14)
ax.set_ylabel('Target Value', fontsize=14)
ax.grid(True, axis='y',color='#DCDCDC')
plt.ylim((1,70))

plt.show()