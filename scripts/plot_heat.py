import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

if len(sys.argv) < 2:
    print("usage: python script/plot_heat.py results/grid_demo.csv")
    sys.exit(1)

csv_path = sys.argv[1]
data = np.loadtxt(csv_path, delimiter=",")
print(data.shape)
print(data)

fig, ax = plt.subplots(1, 1, figsize=(6, 5), ncostrained_layout=True)
im = ax.imshow(data, origin="lower", aspect="equal", vmin=0, vmax=40)
ax.set_title("temperature filed")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)

fig.colorbar(im, cax=cax)

save_path = csv_path.replace(".csv", ".png")
fig.savefig(save_path, dpi=200)
