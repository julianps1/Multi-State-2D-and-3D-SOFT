import numpy as np
import matplotlib.pyplot as plt

def load_surface(fname):
    """
    Reads columns: R_NaH, R_HF, Energy(Ha)
    Blank lines are OK.
    """
    data = np.genfromtxt(fname)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"{fname}: expected 3 columns (R_NaH, R_HF, E).")
    return data[:, 0], data[:, 1], data[:, 2]

def get_slice_fixed_rhf(rNaH, rHF, E, rhf_target, tol=None):
    """
    Extract points with R_HF ~= rhf_target.
    If tol=None -> auto tolerance from grid spacing.
    If none found -> snap to nearest grid R_HF.
    Returns sorted (r, E_slice, rhf_used).
    """
    uniq = np.unique(rHF)
    if uniq.size < 2:
        raise ValueError("Not enough unique R_HF to infer grid spacing.")
    dr = np.min(np.diff(uniq))

    if tol is None:
        tol = 0.5 * dr + 1e-12

    mask = np.abs(rHF - rhf_target) <= tol
    if not np.any(mask):
        rhf_used = float(uniq[np.argmin(np.abs(uniq - rhf_target))])
        mask = np.abs(rHF - rhf_used) <= (0.5 * dr + 1e-12)
    else:
        rhf_used = float(rhf_target)

    r = rNaH[mask]
    Es = E[mask]

    idx = np.argsort(r)
    return r[idx], Es[idx], rhf_used

def adiabatic_energies(V11, V22, V12):
    """
    Adiabatic eigenvalues of 2x2 symmetric Hamiltonian:
      [V11  V12]
      [V12  V22]
    Returns (E_low, E_high).
    """
    tr = 0.5 * (V11 + V22)
    rad = np.sqrt((0.5 * (V11 - V22))**2 + V12**2)
    return tr - rad, tr + rad

def main():
    # -------- USER SETTINGS --------
    RHF_SLICE = 2.0          # fixed R_HF slice (bohr) - change to what you want
    TOL = None               # None = auto
    OUTPNG = f"slice_diab_adiab_RHF_{RHF_SLICE:.3f}.png"

    # Axis window like your reference (optional):
    XLIM = (1.0, 4.0)        # set None to disable
    YLIM = (0.0, 0.2)        # set None to disable

    # --------------------------------

    r11, rhf11, V11_all = load_surface("U11.dat")
    r12, rhf12, V12_all = load_surface("U12.dat")
    r22, rhf22, V22_all = load_surface("U22.dat")

    r,  V11, rhf_used = get_slice_fixed_rhf(r11, rhf11, V11_all, RHF_SLICE, tol=TOL)
    r2, V12, _        = get_slice_fixed_rhf(r12, rhf12, V12_all, rhf_used, tol=TOL)
    r3, V22, _        = get_slice_fixed_rhf(r22, rhf22, V22_all, rhf_used, tol=TOL)

    # Ensure all curves are on same r-grid (interpolate if needed)
    if not np.array_equal(r, r2):
        V12 = np.interp(r, r2, V12)
    if not np.array_equal(r, r3):
        V22 = np.interp(r, r3, V22)

    # Adiabatic curves
    E_low, E_high = adiabatic_energies(V11, V22, V12)

    # -------- Plotting (style like reference) --------
    plt.figure(figsize=(5.0, 3.8))

    # Diabatic: V22 black solid, V11 red solid
    plt.plot(r, V22, "k-", lw=2.0, label=r"$V_{22}$")
    plt.plot(r, V11, "r-", lw=2.0, label=r"$V_{11}$")

    # Coupling: green dashed
    plt.plot(r, V12, "g--", lw=2.0, label=r"$V_{12}$")

    # Adiabatic: dashed black/red (no legend entries in your reference)
    plt.plot(r, E_high, "k--", lw=1.6, dashes=(6, 3))
    plt.plot(r, E_low,  "r--", lw=1.6, dashes=(6, 3))

    plt.xlabel(r"$r\ [a_0]$")
    plt.ylabel("Energy [hartree]")

    # (Optional title — your reference plot has none; keep it minimal)
    # plt.title(f"Slice at $R_{{HF}}$ = {rhf_used:.6f} bohr")

    plt.grid(True, alpha=0.25)
    plt.legend(frameon=False, loc="upper center")

    if XLIM is not None:
        plt.xlim(*XLIM)
    if YLIM is not None:
        plt.ylim(*YLIM)

    plt.tight_layout()
    plt.savefig(OUTPNG, dpi=300)
    print(f"Saved: {OUTPNG}")
    print(f"Used R_HF slice value: {rhf_used:.10f} bohr")

if __name__ == "__main__":
    main()
