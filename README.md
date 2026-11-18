# qiskit-mini-hack

Team **Qiskiteers**' work for the BasQ Qiskit Fall Fest 2026 hackathon (Challenge 1: ground-state energy of a molecule with VQE).

## The challenge

The tutorial computes the ground-state energy of H₂ with the Variational Quantum Eigensolver (VQE), using a Jordan-Wigner mapping and a Hartree-Fock + UCCSD ansatz. For the intermediate challenge we apply the same approach to **lithium hydride (LiH)**. The goal is to cut the quantum resources (qubits, CNOTs, circuit depth) while staying within **chemical accuracy** (1.6 mHa) of the exact energy.

Settings: LiH at a bond length of 1.6 Å, `sto3g` basis.

## Results so far: reducing qubits

We applied reduction techniques one at a time and measured what each one saves. Every error is measured against the exact energy of the full molecule.

| Configuration | Qubits | Params | CNOTs | Depth | Error (mHa) | Chem. accurate |
|---|---:|---:|---:|---:|---:|:---:|
| Baseline (Jordan-Wigner) | 12 | 92 | 1572 | 3180 | 0.000 | ✅ |
| A: freeze core | 10 | 24 | 346 | 698 | 0.228 | ✅ |
| A+B: + parity mapping | 8 | 24 | 410 | 701 | 0.228 | ✅ |
| A+B+C: + symmetry tapering | **6** | 10 | 118 | 242 | 0.228 | ✅ |
| D (naive): lowest 3 orbitals | 3 | 4 | 18 | 53 | 19.405 | ❌ |
| A+B+C+D: drop the π orbital pair | **4** | 8 | 72 | 162 | 1.252 | ✅ |

What we found:

- **12 → 6 qubits** with almost no loss in accuracy (0.23 mHa error), and CNOTs drop from 1572 to 118.
- **12 → 4 qubits** stays just inside chemical accuracy (1.25 mHa) with only 72 CNOTs.
- **Which** orbitals you drop matters more than how many. Keeping the lowest-energy orbitals fails (19 mHa error). The π orbitals are perpendicular to the Li–H bond, so by symmetry they can't mix with it, and dropping them is the right choice.

For the 6- and 4-qubit versions, VQE reaches the exact ground state of the reduced problem.

## Repository contents

| File | What it is |
|---|---|
| `challenge_1_tutorial qiskiteers.ipynb` | The challenge tutorial (VQE for H₂), with our notes |
| `p1_qubits.ipynb` | The baseline on LiH and the qubit-reduction study |
| `p1_setup.py` | Shared code: molecule, reduction levels, ansatz, resource counting, VQE |
| `p1_results.csv` | The results table above, in CSV form |

### Using the shared setup

Teammates import `p1_setup.py` so everyone works on the same molecule, exact energy and resource counts:

```python
from p1_setup import *

problem, mapper = reduced_lih(level="6q")   # "full" (12q), "6q" or "4q"
qubit_op = mapper.map(problem.hamiltonian.second_q_op())
ansatz = build_ansatz(problem, mapper)

print(measure_resources(ansatz, qubit_op, label="6q"))
print(run_vqe(problem, mapper, ansatz))
```

## Requirements

Python 3.10+ with `qiskit`, `qiskit-nature`, `qiskit-algorithms` and `pyscf`.

## Credits

The tutorial notebook was written by Benjamin Tirado for the BasQ Qiskit Fall Fest 2026, organised by Basque Quantum (BasQ). Its text is licensed under CC BY 4.0 and its code under Apache 2.0.
