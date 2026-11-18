# qiskit-mini-hack

Team **Qiskiteers**' work for the BasQ Qiskit Fall Fest 2026 hackathon.

## Repository layout

| Folder | File | What it is |
|---|---|---|
| [`beginner/`](./beginner/) | `challenge_1_tutorial qiskiteers.ipynb` | The H₂ VQE tutorial: a single-geometry introduction. Start here. |
| | `heh_plus_dissociation.ipynb` | HeH⁺ dissociation curve |
| [`intermediate/`](./intermediate/) | `p1_qubits.ipynb` | LiH baseline and qubit-reduction study |
| | `p1_setup.py` | Shared LiH model, reduction levels, ansatz, resource counting and VQE helpers |
| | `p1_results.csv` | Resource and accuracy results (the table below) |
| | `p1_bond_scan.csv` | 6- and 4-qubit errors across LiH bond lengths |

## Beginner: HeH⁺ dissociation curve

After the H₂ tutorial, `heh_plus_dissociation.ipynb` scans the HeH⁺ bond length to estimate the equilibrium geometry and dissociation energy, and compares UCCSD/VQE with a particle-number-constrained exact reference.

It uses PySCF with the minimal STO-3G basis, Jordan–Wigner mapping, Hartree–Fock + UCCSD, and a noiseless statevector estimator. Its equilibrium and dissociation values are estimates within this model, not basis-set-limit or experimental values.

## Intermediate: LiH qubit reduction

We apply qubit-reduction techniques to LiH at 1.6 Å in the STO-3G basis, aiming to reduce quantum resources while staying within chemical accuracy (1.6 mHa). Every error is measured against the exact energy of the full molecule.

| Configuration | Qubits | Params | CNOTs | Depth | Error (mHa) | Chem. accurate |
|---|---:|---:|---:|---:|---:|:---:|
| Baseline (Jordan-Wigner) | 12 | 92 | 1572 | 3180 | 0.000 | Yes |
| Freeze core | 10 | 24 | 346 | 698 | 0.228 | Yes |
| Freeze core + parity mapping | 8 | 24 | 410 | 701 | 0.228 | Yes |
| Freeze core + parity + tapering | **6** | 10 | 118 | 242 | 0.228 | Yes |
| Naive lowest-three-orbital active space | 3 | 4 | 18 | 53 | 19.405 | No |
| Freeze core + parity + drop π pair | **4** | 8 | 72 | 162 | 1.252 | Yes |

The 6-qubit reduction keeps its accuracy with far fewer CNOTs, and the 4-qubit active-space model stays just inside chemical accuracy. The choice of orbitals matters: the naive lowest-three-orbital selection misses the accuracy target.

### Using the shared setup

Run notebooks with the repository root as the Jupyter working directory, and import the shared helper with:

```python
from intermediate.p1_setup import *

problem, mapper = reduced_lih(level="6q")  # "full" (12q), "6q" or "4q"
qubit_op = mapper.map(problem.hamiltonian.second_q_op())
ansatz = build_ansatz(problem, mapper)
```

## Requirements

Python 3.10+ with `qiskit`, `qiskit-nature`, `qiskit-algorithms`, `pyscf`, `numpy`, and `matplotlib`.

## Credits

The H₂ tutorial was written by Benjamin Tirado for the BasQ Qiskit Fall Fest 2026, organised by Basque Quantum (BasQ). Its text is licensed under CC BY 4.0 and its code under Apache 2.0.
