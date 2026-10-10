# Isospectrality breaking in Reissner–Nordström black holes in dynamical Chern–Simons gravity

Code and data accompanying

> A. Rout, *Isospectrality breaking in the axial quasinormal modes of Reissner–Nordström black holes in dynamical Chern–Simons gravity* (2026).

The repository contains:

- the computer-algebra derivations of the axial couplings: the C-tensor, the linearized Pontryagin density, the reduction to master variables, and the consistency checks;
- the spectral solvers used for the quasinormal-mode (QNM) shifts;
- one script per numerical result in the paper, together with the data behind the tables and figures.

---

## Physics in one paragraph

A dynamical pseudoscalar coupled to the Pontryagin density, $`\tfrac{\alpha}{4}\vartheta\,{}^{\ast}\!R R`$, couples only to the axial (odd-parity) sector at linear order. It therefore shifts the axial QNMs of a Reissner–Nordström (RN) black hole while the polar ones stay unchanged, which breaks isospectrality.

The axial system closed here is

```math
\begin{aligned}
[\partial_{r_{\ast}}^2+\omega^2-V_{\rm RW}]\Psi &= C_Q\Psi_e + C_{\rm dCS}\Theta,\\
[\partial_{r_{\ast}}^2+\omega^2-V_{\rm EM}]\Psi_e &= C_Q\Psi + \chi\,C_{\rm dCS}\Theta,\\
[\partial_{r_{\ast}}^2+\omega^2-V_\Theta]\Theta &= C_{\Theta g}(\Psi+\chi\Psi_e),
\end{aligned}
```

with

- $C_{\rm dCS}=96\pi(Mr-Q^2)\Delta/r^8$,
- $C_{\Theta g}=6(Mr-Q^2)(l-1)l(l+1)(l+2)\Delta/(\beta r^8)$,
- $\chi=2Q/[r\sqrt{(l-1)(l+2)}]$,
- $\Psi=i\Psi_{\rm RW}/\omega$.

All shifts are quoted as $\beta\,\delta\omega$ in units $M=\alpha=1$ (with $\kappa=8\pi$). Since $\zeta=\alpha^2/(\beta\kappa M^4)=1/(8\pi\beta)$, they convert as $M\delta\omega = 8\pi(\beta\delta\omega)\,\zeta$.

---

## Repository layout

```
src/            core solvers (imported by the scripts)
  ecs.py          exterior-complex-scaling (ECS) Chebyshev solver for the coupled {Psi, Psi_e, Theta} system
  coef.py         linear coefficient beta*d(omega) via Richardson extrapolation in 1/beta
  branch.py       branch-reduced {Z_a, Theta} system (exact at O(zeta)); single-field RN branch/scalar spectra
  shoot.py        Frobenius (horizon) + complex-ray (infinity) shooting solver for near-extremal zero-damping modes
  zdm_tools.py    ZDM root finding on the imaginary axis, coupled ZDM matrix
symbolic/       SymPy derivations (Secs. V-VII, Apps. A-C)
scripts/        one script per numerical result (Secs. IX, XI, XIII-XIV, Apps. D-E)
figures/        scripts that regenerate Figs. 1 and 3 from data/
data/           JSON outputs used in the paper's tables and figures
```

> **Branch labels.** In `src/branch.py` and `src/zdm_tools.py`, `branch=1` / `'G'` is the gravitational-led branch (potential with $q_1$, called $Z_g=Z_2^{(-)}$ in the paper), and `branch=2` / `'EM'` is the electromagnetic-led branch ($q_2$, $Z_e=Z_1^{(-)}$).

---

## Installation

```bash
git clone https://github.com/AbhiRout-cosmo/RN-dCS.git
cd RN-dCS
pip install -r requirements.txt
```

Requires Python ≥ 3.9 with NumPy, SciPy, SymPy and Matplotlib. Every script runs from the repository root and finds `src/` and `data/` by itself.

---

## Reproducing the paper

### Symbolic derivations (`symbolic/`)

| Script | Paper | Result |
|---|---|---|
| `01_pontryagin_density.py` | Sec. VII A, App. C | $\delta({}^{\ast}RR)\propto L(Mr-Q^2)[r h_0'-2h_0+i\omega r h_1]Y_l/r^7$, exact at all $Q$ |
| `02_axial_einstein_maxwell.py` | Secs. IV, VI | linearized axial Einstein–Maxwell equations on RN |
| `03_cotton_covariant.py` | Sec. VI, App. B | covariant $C_{t\varphi}$, $C_{r\varphi}$; $C_{\theta\varphi}=0$ |
| `04_cotton_contravariant.py` | Sec. V, App. A | contravariant $C_{(1)}$, $C_{(2)}$ in the paper's conventions |
| `05_master_reduction.py` | Secs. VI–VII, App. B | $V_{\rm RW}$, $V_{\rm EM}$, $C_Q$, $C_{\rm dCS}$ and the Maxwell source $\chi C_{\rm dCS}$ |
| `06_tphi_consistency.py` | Apps. B–C | the unused $t\varphi$ equation vanishes on shell (exact rational check) |

Orientation: scripts 01, 03 and 05 use $\varepsilon^{0123}=+1/\sqrt{-g}$. The paper uses $\varepsilon_{tr\theta\varphi}=+\sqrt{-g}$, which flips the sign of every dCS coupling at once. Frequencies are unaffected.

### Numerical results (`scripts/`)

| Script | Paper | Output | Runtime |
|---|---|---|---|
| `01_gates_and_neutral_shift.py` | Secs. XIII–XIV, Eq. (82) | ECS gates; $\beta\delta\omega\vert_{Q=0}=-0.7733-0.3946i$, grid-independent | ~1 min |
| `02_charge_sweep.py` | Sec. XIV C, Table III, Fig. 3 | `data/sweep.json` ($l=2$–6, $Q/M\le0.98$, with/without $\chi$) | ~10 min |
| `03_convergence_checks.py` | App. E | convergence at high $Q$ and $l$ | ~5 min |
| `04_reduced_vs_full.py` | Sec. VIII A, App. E | reduced $\{Z_g,\Theta\}$ vs full system ($10^{-7}$) | ~3 min |
| `05_branch_sweeps.py` | Sec. IX, Fig. 1 | `data/c3_sweeps.json`: EM-led sign reversal, $Q_c(l)$, $l^3$ vs $l^2/Q$ | ~10 min |
| `06_em_crossing_spectra.py` | Sec. IX C | `data/spectra.json` (input to 07) | ~5 min |
| `07_avoided_crossing.py` | Sec. IX C | $\mathcal{O}(\sqrt\zeta)$ splitting near $Q_c$ | ~3 min |
| `08_zimmerman_mark_validation.py` | App. E | ZDM solver vs Zimmerman & Mark (2016), Fig. 10 | ~10 min |
| `09_zdm_shift.py` | Sec. XI C, Eq. (75) | `data/zdm_shift.json`; $C_0=481.7$ | ~5 min |
| `10_zdm_overtone.py` | Sec. XI C | `data/zdm_n1.json`; near-degenerate overtone | ~8 min |
| `11_local_wkb_proxy.py` | App. D, Table V | local-WKB proxy values | seconds |

Figures:
```bash
python figures/make_fig_branches.py   # Fig. 1 -> figures/fig_branches.pdf
python figures/make_fig_shift.py      # Fig. 3 -> figures/fig_shift_ecs.pdf
```

Runtimes are for a 2-core laptop.

---

## Numerical methods (summary of App. E)

**Exterior complex scaling.** In the $\Psi=i\Psi_{\rm RW}/\omega$ normalization all couplings are algebraic and independent of $\omega$, so the coupled problem is a linear eigenproblem in $\omega^2$.

- *Contour.* The tortoise coordinate is rotated, $r_*=e^{i\theta}s$ with $s\in[-S,S]$. QNM boundary conditions then become decay at both ends, imposed as Dirichlet conditions.
- *Radial map.* $r(s)$ is obtained by integrating $du/ds=e^{i\theta}(r-r_-)/r^2$ for $u=\ln(r-r_+)$, which avoids log branch cuts.
- *Discretization.* Chebyshev collocation followed by a dense eigensolve.
- *Shift extraction.* $\beta\delta\omega$ is Richardson-extrapolated from $\beta=2\times10^4$ and $4\times10^4$. Grid errors cancel because coupled and uncoupled eigenvalues come from the same grid.
- *Gates.*
  - Schwarzschild $l=2$ to $8\times10^{-9}$.
  - The Hatsuda–Kimura charge series to $8\times10^{-9}$ at $Q/M=0.1$.
  - Berti–Kokkotas at all charges.
- *Robustness.* Results are stable to six digits under changes of $(\theta,S,N)$.

**Zero-damping modes.** Near extremality $\omega/\kappa_+$ is close to an integer, so the horizon exponents nearly coincide and Chebyshev or Leaver methods fail. The shooting solver in `shoot.py` handles this:

- *Horizon side.* A matrix Frobenius series (order 60), whose Taylor coefficients come from an FFT on $|r-r_+|=(r_+-r_-)/2$.
- *Infinity side.* Outgoing data started on the ray $r=8M+t\,e^{3\pi i/4}$ and integrated inward.
- *QNM condition.* The determinant of horizon and outgoing solutions vanishes at $r_++2(r_+-r_-)$.
- *Validation.* It reproduces the ZDM residuals of Zimmerman & Mark.

---

## Key numbers

| Quantity | Value |
|---|---|
| Neutral shift, $l=2$ | $\beta\delta\omega=-0.7733-0.3946i$ |
| Large-$l$ law, $Q=0$ | $M\delta\omega\simeq-2.25\,\zeta\,l^3$ |
| Large-$l$ law, $Q>0$ | $\propto l^2/Q$ |
| EM-led crossing | $Q_c(l)\to2M/\sqrt{(l-1)(l+2)}$ ($0.661M$ at $l=2$) |
| Fundamental ZDM, $l=2$ | $\delta\omega/\omega\simeq6.05\times10^3\,\zeta(1-Q^2/M^2)$ |

---

## Citation

If you use this code, please cite the paper (see `CITATION.cff`).

## License

MIT. See `LICENSE`.
