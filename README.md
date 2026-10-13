# Isospectrality breaking in Reissner–Nordström black holes in dynamical Chern–Simons gravity

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23283164.svg)](https://doi.org/10.5281/zenodo.23283164)

Code and data accompanying

> A. Rout, *Isospectrality breaking in the axial quasinormal modes of Reissner–Nordström black holes in dynamical Chern–Simons gravity* (2026).


The repository contains:

- the computer-algebra derivations of the axial couplings: the C-tensor, the linearized Pontryagin density, the reduction to master variables, and the consistency checks;
- the spectral solvers used for the quasinormal-mode (QNM) shifts;
- one script per numerical result in the paper, together with the data behind the tables and figures;
- an independent first-principles re-derivation of the whole axial system (`symbolic/08_independent_audit/`).

---

## Physics in one paragraph

A dynamical pseudoscalar coupled to the Pontryagin density, $`\tfrac{\alpha}{4}\vartheta\,{}^{\ast}\!R R`$, couples only to the axial (odd-parity) sector at linear order. It therefore shifts the axial QNMs of a Reissner–Nordström (RN) black hole while the polar ones stay unchanged, which breaks isospectrality.

The axial system closed here is

```math
\begin{aligned}
\bigl[\partial_{r_{\ast}}^2+\omega^2-V_{\mathrm{RW}}\bigr]\Psi &= C_Q\Psi_e + C_{\mathrm{dCS}}\Theta, \\
\bigl[\partial_{r_{\ast}}^2+\omega^2-V_{\mathrm{EM}}\bigr]\Psi_e &= C_Q\Psi + \chi\,C_{\mathrm{dCS}}\Theta, \\
\bigl[\partial_{r_{\ast}}^2+\omega^2-V_\Theta\bigr]\Theta &= C_{\Theta g}(\Psi+\chi\Psi_e),
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
  ecs_pivot.py    ECS with movable rotation pivot, full 3-field operator (used by scripts 07 and 12-18)
symbolic/       SymPy derivations (Secs. V-IX, Apps. A-C) and an independent audit
scripts/        one script per numerical result (Secs. IX, XI, XIII-XIV, Apps. D-E)
figures/        scripts that regenerate the paper's figures from data/
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
| `07_exact_Qc.py` | Sec. IX C | at $Q_c=2M\sqrt{L+1}/(L+2)$ the scalar potential equals the polar EM-led potential, $V_\Theta\equiv V_e^{(+)}$, for generic $l$ ⇒ exact degeneracy |
| `08_independent_audit/` | Secs. II–VIII, App. E | first-principles linearization; master system substituted back (generic $l,Q,\omega$); $\delta V_\Theta$ at finite $Q$; identity $\nabla_\mu C^{\mu\nu}=-\tfrac18\nabla^\nu\vartheta\,{}^{\ast}RR$ (see its README) |

Orientation: scripts 01, 03 and 05 use $\varepsilon^{0123}=+1/\sqrt{-g}$. The paper uses $\varepsilon_{tr\theta\varphi}=+\sqrt{-g}$, which flips the sign of every dCS coupling at once. Frequencies are unaffected.

### Numerical results (`scripts/`)

| Script | Paper | Output | Runtime |
|---|---|---|---|
| `01_gates_and_neutral_shift.py` | Secs. XIII–XIV | ECS gates; $\beta\delta\omega\vert_{Q=0}=-0.7733-0.3946i$, grid-independent | ~1 min |
| `02_charge_sweep.py` | Sec. XIV C, Table III, shift figure | `data/sweep.json` ($l=2$–6, $Q/M\le0.98$, with/without $\chi$) | ~10 min |
| `03_convergence_checks.py` | App. E | convergence at high $Q$ and $l$ | ~5 min |
| `04_reduced_vs_full.py` | Sec. VIII A, App. E | reduced $\{Z_g,\Theta\}$ vs full system ($10^{-7}$) | ~3 min |
| `05_branch_sweeps.py` | Sec. IX, branch figure | `data/c3_sweeps.json`: EM-led sign reversal, $Q_c(l)$, $l^3$ vs $l^2/Q$ | ~10 min |
| `06_em_crossing_spectra.py` | Sec. IX C | `data/spectra.json` (decoupled spectra; nearest-eigenvalue tracking — use 15 for branch-resolved tracking) | ~5 min |
| `07_avoided_crossing.py` | Sec. IX C | `data/avoided.json`: splitting at the exact $Q_c$, $\vert\omega_+-\omega_-\vert\sqrt\beta=0.47995$ ($l=2$), $0.57617$ ($l=3$) for $\beta=10^3$–$10^{10}$ | ~3 min |
| `08_zimmerman_mark_validation.py` | App. E | ZDM solver vs Zimmerman & Mark (2016), Fig. 10 | ~10 min |
| `09_zdm_shift.py` | Sec. XI C | `data/zdm_shift.json`; $C_0=481.7$ | ~5 min |
| `10_zdm_overtone.py` | Sec. XI C | `data/zdm_n1.json`; near-degenerate overtone | ~8 min |
| `11_local_wkb_proxy.py` | App. D | local-WKB proxy values | seconds |
| `12_shift_convergence.py` | Sec. XIV A, App. E | `data/shift_convergence.json`: $\beta\delta\omega$ vs $(\theta,S,N)$ and rotation pivot; $\le4\times10^{-5}$ | ~40 min |
| `13_linearity_in_zeta.py` | App. E | `data/linearity_in_zeta.json`: 7-point fit in $1/\beta$ | ~15 min |
| `14_reduced_vs_full_table.py` | Sec. VIII, App. E | `data/reduced_vs_full_table.json`: $l=2,4,6$; difference $\propto\beta^{-2}$ | ~20 min |
| `15_avoided_crossing_branches.py` | Sec. IX C | `data/avoided_branches.json`: both modes through $Q_c$, character exchange | ~60 min |
| `16_shooting_vs_ecs.py` | App. E | `data/shooting_vs_ecs.json`: full-system shooting vs ECS, $1.4\times10^{-6}$ | ~4 min |
| `17_zdm_scaling.py` | Sec. XI C | `data/zdm_scaling.json`: free-exponent fit, $\delta\omega/\omega_0\propto b^{p}$, $p\to2$ | ~8 min |
| `18_eikonal_high_l.py` | Sec. IX A | `data/eikonal_high_l.json`: $l=12$–$30$ for the approach to the $l^3$, $l^2/Q$ laws | ~30 min |

Figures:
```bash
python figures/make_fig_branches.py   # branch-resolved shifts      -> figures/fig_branches.pdf
python figures/make_fig_eikonal.py    # approach to eikonal laws    -> figures/fig_R3_eikonal.pdf
python figures/make_fig_avoided.py    # avoided crossing at Q_c     -> figures/fig_R1_avoided.pdf
python figures/make_fig_zdm.py        # ZDM scaling                 -> figures/fig_R2_zdm.pdf
python figures/make_fig_shift.py      # charge dependence of shift  -> figures/fig_shift_ecs.pdf
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
- *Robustness.* The shift itself changes by $\le4\times10^{-5}$ under $\theta\in[0.5,0.8]$, $S\in[45,120]$, $N\in[220,460]$ and rotation pivot $r_*^c\in[-2,5]M$; full-system shooting agrees to $1.4\times10^{-6}$.

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
| EM-led/pseudoscalar degeneracy (exact) | $Q_c=2M\sqrt{L+1}/(L+2)$ ($0.66144M$ at $l=2$) |
| Splitting at $Q_c$ | $M\vert\Delta\omega\vert=2.41\sqrt\zeta$ ($l=2$), $2.89\sqrt\zeta$ ($l=3$) |
| Back-reaction potential | $\delta V_\Theta=576\pi L\alpha^2(Mr-Q^2)^2\Delta/(\beta r^{12})$ |
| Fundamental ZDM, $l=2$ | $\delta\omega/\omega\simeq6.05\times10^3\,\zeta(1-Q^2/M^2)$, $C_0=481.7\pm0.6$ |

---

## Version history

- **v1.1.0** (2026-10) — revision for the second referee report: exact $Q_c$ and splitting at $Q_c$ (script 07 rewritten; the v1.0.0 version used a mis-tracked minimum near $Q=0.655$), scripts 12–18, `symbolic/07_exact_Qc.py`, `symbolic/08_independent_audit/`, three new figures, `src/ecs_pivot.py`.
- **v1.0.0** (2026-10-10) — first release.

## Citation

If you use this code, please cite the paper (see `CITATION.cff`).

## License

MIT. See `LICENSE`.
