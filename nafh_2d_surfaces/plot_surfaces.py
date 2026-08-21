import numpy as np
import matplotlib.pyplot as plt

def load_grid(fname):
    # Because we wrote blank lines between rows, genfromtxt will skip them
    data = np.genfromtxt(fname)
    x = data[:, 0]  # R_NaH
    y = data[:, 1]  # R_HF
    z = data[:, 2]  # Energy

    # infer grid shape by counting repeats in x
    # assumes each row has nx points and y increases after each row
    unique_x = np.unique(x)
    unique_y = np.unique(y)
    nx = unique_x.size
    ny = unique_y.size

    X = unique_x.reshape(1, nx).repeat(ny, axis=0)
    Y = unique_y.reshape(ny, 1).repeat(nx, axis=1)
    Z = z.reshape(ny, nx)
    return X, Y, Z

for name in ["U11", "U12", "U22"]:
    X, Y, Z = load_grid(f"{name}.dat")

    plt.figure()
    plt.contourf(X, Y, Z, levels=40)
    plt.colorbar(label="Energy (Ha)")
    plt.xlabel("R_NaH (bohr)")
    plt.ylabel("R_HF (bohr)")
    plt.title(name)
    plt.tight_layout()
    plt.savefig(f"{name}.png", dpi=300)
    plt.close()

print("Saved U11.png, U12.png, U22.png")
