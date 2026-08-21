import numpy as np
import matplotlib.pyplot as plt

def load_surface(fname):
    """
    File format: columns = R_NaH, R_HF, Energy(Ha)
    Blank lines are allowed.
    """
    data = np.genfromtxt(fname)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"Bad data in {fname}. Expected 3 columns.")
    rNaH = data[:, 0]
    rHF  = data[:, 1]
    E    = data[:, 2]
    return rNaH, rHF, E

def slice_at_fixed_rhf(rNaH, rHF, E, rhf_slice, tol=None):
    """
    Return (x, y) for points with rHF ~ rhf_slice.
    If tol is None, auto-choose tol as half of the smallest rHF grid spacing.
    """
    # infer rHF grid spacing
    unique_rhf = np.unique(rHF)
    if unique_rhf.size < 2:
        raise ValueError("Not enough unique R_HF values to infer grid.")
    dr = np.min(np.diff(unique_rhf))

    if tol is None:
        tol = 0.5 * dr + 1e-12

    mask = np.abs(rHF - rhf_slice) <= tol
    x = rNaH[mask]
    y = E[mask]

    if x.size == 0:
        # pick the nearest available slice automatically
        nearest = unique_rhf[np.argmin(np.abs(unique_rhf - rhf_slice))]
        mask = np.abs(rHF - nearest) <= (0.5 * dr + 1e-12)
        x = rNaH[mask]
        y = E[mask]
        return x, y, nearest, tol

    return x, y, rhf_slice, tol

def sort_xy(x, y):
    idx = np.argsort(x)
    return x[idx], y[idx]

def main():
    # ---- USER SETTINGS ----
    RHF_slice = 2.0      # choose your fixed R_HF value (bohr)
    tol = None           # None = auto tolerance (recommended)
    out_png = f"slice_RHF_{RHF_slice:.3f}.png"

    files = {
        "U11": "U11.dat",
        "U12": "U12.dat",
        "U22": "U22.dat",
    }

    plt.figure(figsize=(8, 5))

    used_slice_value = None

    for label, fname in files.items():
        rNaH, rHF, E = load_surface(fname)
        x, y, used_rhf, used_tol = slice_at_fixed_rhf(rNaH, rHF, E, RHF_slice, tol=tol)
        x, y = sort_xy(x, y)

        if used_slice_value is None:
            used_slice_value = used_rhf

        plt.plot(x, y, linewidth=2, label=label)

    plt.xlabel("R_NaH (bohr)")
    plt.ylabel("Energy (Hartree)")
    plt.title(f"NaFH 1D slice at R_HF = {used_slice_value:.6f} bohr")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(out_png, dpi=300)
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
