import matplotlib.pyplot as plt
import numpy as np

record = np.loadtxt("results/dif_record.csv", delimiter=",")
iteration = record[:, 0]
dif = record[:, 1]

fig, ax = plt.subplots(1, 2, figsize=(8, 5))
im1 = ax[0].semilogy(iteration, dif)
im2 = ax[1].plot(iteration, dif)

ax[0].set_title("semilogy_iteration_dif")
ax[1].set_title("plot_iteration_dif")

ax[0].set_xlabel("iteration")
ax[1].set_xlabel("iteration")
ax[0].set_ylabel("dif")
ax[1].set_ylabel("dif")
fig.savefig("results/iteration_dif.png", dpi=200)
