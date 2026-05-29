import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

heat_serial = np.loadtxt("results/heat2d/serial_heat.csv", delimiter=",")
heat_cuda = np.loadtxt("results/heat2d/cuda_heat.csv", delimiter=",")
heat_openmp = np.loadtxt("results/heat2d/openmp_heat.csv", delimiter=",")

diff_max_s2c = np.max(np.abs(heat_serial - heat_cuda))
diff_max_s2o = np.max(np.abs(heat_serial - heat_openmp))
diff_max_c2o = np.max(np.abs(heat_cuda - heat_openmp))

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
im = ax.imshow(heat_serial - heat_cuda, origin="lower", aspect="equal", cmap="coolwarm")
ax.set_title("Temperature difference: serial - CUDA")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)

cbar = fig.colorbar(im, cax=cax)
cbar.formatter = FormatStrFormatter("%.1e")
cbar.update_ticks()

save_path = "results/heat2d/dif_field.png"
fig.savefig(save_path, dpi=200)

print(f"max difference between serial and cuda is {diff_max_s2c}\n")
print(f"max difference between serial  and openmp is {diff_max_s2o}\n")
print(f"max difference between cuda  and openmp is {diff_max_c2o}\n")
print(f"max temperature : {np.max(np.abs(heat_cuda))}\n")
print(f"relative max error: {diff_max_c2o/np.max(np.abs(heat_cuda))*100}%")
