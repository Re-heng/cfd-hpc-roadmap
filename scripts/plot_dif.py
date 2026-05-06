import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory

record_jacobi = np.loadtxt("results/jacobi_convence.csv", delimiter=",")
iteration_jacobi = record_jacobi[:, 0]
dif_jacobi = record_jacobi[:, 1]

record_gauss = np.loadtxt("results/gauss_convence.csv", delimiter=",")
iteration_gauss = record_gauss[:, 0]
dif_gauss = record_gauss[:, 1]

record_SOR = np.loadtxt("results/SOR_convence.csv", delimiter=",")
iteration_SOR = record_SOR[:, 0]
dif_SOR = record_SOR[:, 1]

final_jacobi = iteration_jacobi[-1]
final_gauss = iteration_gauss[-1]
final_SOR = iteration_SOR[-1]

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
im1 = ax.semilogy(iteration_jacobi, dif_jacobi, label="jacobi", color="blue")
im1 = ax.semilogy(iteration_gauss, dif_gauss, label="gauss", color="orange")
im1 = ax.semilogy(iteration_SOR, dif_SOR, label="SOR(omega = 1.5)", color="green")

ax.axvline(x=final_jacobi, color="blue", linestyle="--", alpha=0.6)
ax.axvline(x=final_gauss, color="orange", linestyle="--", alpha=0.6)
ax.axvline(x=final_SOR, color="green", linestyle="--", alpha=0.6)
ax.axhline(y=1e-6, color="red", linestyle="--", alpha=0.6)

transformmix = blended_transform_factory(ax.transData, ax.transAxes)

ax.text(
    final_jacobi - 300,
    0.6,
    f"Iter = {final_jacobi}",
    transform=transformmix,
    color="blue",
    ha="right",
    va="center",
)
ax.text(
    final_gauss - 300,
    0.6,
    f"Iter = {final_gauss}",
    transform=transformmix,
    color="orange",
    ha="right",
    va="center",
)
ax.text(
    final_SOR - 300,
    0.6,
    f"Iter = {final_SOR}",
    transform=transformmix,
    color="green",
    ha="right",
    va="center",
)

ax.text(
    0,
    1e-6 * 1.1,
    f"goal_diff = 1e-6",
    transform=ax.transData,
    color="red",
    va="bottom",
    ha="left",
)
ax.set_title("three method convergence performance")
ax.set_xlabel("iteration")
ax.set_ylabel("diff")


ax.legend(
    loc="upper right",
    bbox_to_anchor=(final_jacobi - 500, 0.95),
    bbox_transform=transformmix,
)
ax.grid(True, which="both", linestyle="--", alpha=0.4)

fig.savefig("results/convence_dif.png", dpi=200)
