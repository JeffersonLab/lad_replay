"""Proton energy deposition and scintillation light in the LAD bars (5.08 cm polyvinyltoluene).

Bethe stopping power (no shell or density corrections; fine for 5-400 MeV protons), integrated to a
CSDA range table, plus Birks' law for the light. Used by lad_edep_fit.py.

Conventions
  T        proton kinetic energy (MeV)
  t_bar    bar thickness along the track (g/cm^2)
  pre      material before the bar (g/cm^2, air + windows + GEMs; for back planes add the front bar)
"""
import numpy as np

M_P = 938.272  # MeV
M_E = 0.51099895
K = 0.307075  # MeV cm^2 / g
ZA_PVT = 0.54141  # polyvinyltoluene C9H10
I_PVT = 64.7e-6  # MeV
RHO_PVT = 1.032  # g/cm^3
THICK_CM = 5.08
T_BAR = THICK_CM * RHO_PVT  # g/cm^2 at normal incidence
ZA_AIR, I_AIR = 0.49919, 85.7e-6


def dedx(T, za=ZA_PVT, I=I_PVT):
    """Mass stopping power (MeV cm^2/g) for protons of kinetic energy T (MeV)."""
    T = np.asarray(T, dtype=float)
    g = 1 + T / M_P
    b2 = 1 - 1 / g**2
    tmax = 2 * M_E * b2 * g**2 / (1 + 2 * g * M_E / M_P + (M_E / M_P) ** 2)
    arg = 2 * M_E * b2 * g**2 * tmax / I**2
    return K * za / b2 * (0.5 * np.log(arg) - b2)


# CSDA range table in PVT; below 1 MeV use R ~ T^1.75 scaling (only matters for the last ~0.003 g/cm2)
_T = np.concatenate([np.linspace(0.05, 1, 40), np.logspace(0, np.log10(600), 600)[1:]])
_S = dedx(_T)
_R = np.concatenate([[0.0], np.cumsum(0.5 * (1 / _S[1:] + 1 / _S[:-1]) * np.diff(_T))])
_R += _T[0] / _S[0] / 1.75  # range of the first point


def rng(T):
    return np.interp(T, _T, _R, left=0.0)


def t_of_range(R):
    return np.interp(R, _R, _T, left=0.0)


def t_out(T, t):
    """Kinetic energy after t g/cm^2 of PVT (0 if the proton stops)."""
    R = rng(T) - t
    return np.where(R > 0, t_of_range(np.maximum(R, 0)), 0.0)


def edep(T, t=T_BAR):
    """Energy deposited in a bar of thickness t (g/cm^2) by a proton entering with T."""
    return T - t_out(T, t)


def punch_through_T(t=T_BAR):
    """Kinetic energy at which a proton just traverses thickness t."""
    return float(t_of_range(t))


def light(T, t=T_BAR, kB=0.0, nstep=200):
    """Birks light (MeV-electron-equivalent units) for a proton entering with T, thickness t.
    kB in g/(cm^2 MeV) (= kB[cm/MeV] * rho). kB=0 gives edep."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if kB == 0:
        return edep(T, t)
    out = np.empty_like(T)
    for i, Ti in enumerate(T):
        Tf = float(t_out(Ti, t))
        Es = np.linspace(Tf, Ti, nstep)
        S = dedx(np.maximum(Es, 0.05))
        # dL/dE = 1 / (1 + kB S)
        f = 1 / (1 + kB * S)
        out[i] = np.trapz(f, Es)
    return out


def t_from_invbeta(ib):
    ib = np.asarray(ib, dtype=float)
    b = 1 / np.maximum(ib, 1.000001)
    return M_P * (1 / np.sqrt(1 - b**2) - 1)


def invbeta_from_t(T):
    g = 1 + np.asarray(T, dtype=float) / M_P
    return 1 / np.sqrt(1 - 1 / g**2)


def air_loss_t(T, cm_air=600.0, extra=0.0):
    """Kinetic energy after cm_air of air plus extra g/cm^2 of PVT-equivalent material (approx)."""
    t_air = cm_air * 1.205e-3 * (ZA_AIR / ZA_PVT)  # crude PVT-equivalent
    return t_out(T, t_air + extra)


if __name__ == "__main__":
    for t in [T_BAR, 2 * T_BAR]:
        print(f"punch-through T for {t:.3f} g/cm2: {punch_through_T(t):.1f} MeV")
    for T in [10, 30, 50, 80, 100, 150, 200]:
        print(T, f"S={float(dedx(T)):.2f}", f"R={float(rng(T)):.3f}", f"edep={float(edep(T)):.1f}",
              f"light(kB=.013)={float(light(T, kB=0.0126*RHO_PVT)[0]):.1f}")
