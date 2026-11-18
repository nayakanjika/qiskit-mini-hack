# qiskit-mini-hack

Team **Qiskiteers**' work for the BasQ Qiskit Fall Fest 2026 hackathon.

## Challenges

- [`beginner_challenge/`](./beginner_challenge/): H₂ VQE tutorial and a HeH⁺ dissociation-curve exercise.
- [`intermediate_challenge/`](./intermediate_challenge/): LiH qubit-reduction notebook, shared setup, and results.

The beginner exercise uses STO-3G and scans the HeH⁺ bond length to estimate the equilibrium geometry, dissociation energy, and UCCSD/VQE accuracy against a particle-number-filtered exact reference.

The intermediate study applies qubit-reduction techniques to LiH at 1.6 Å in the STO-3G basis, with the goal of reducing quantum resources while staying within chemical accuracy (1.6 mHa).

## Intermediate results: reducing qubits

Every error is measured against the exact energy of the full molecule.

| Configuration | Qubits | Params | CNOTs | Depth | Error (mHa) | Chem. accurate |
|---|---:|---:|---:|---:|---:|:---:|
| Baseline (Jordan-Wigner) | 12 | 92 | 1572 | 3180 | 0.000 | Yes |
| Freeze core | 10 | 24 | 346 | 698 | 0.228 | Yes |
| Freeze core + parity mapping | 8 | 24 | 410 | 701 | 0.228 | Yes |
| Freeze core + parity + tapering | **6** | 10 | 118 | 242 | 0.228 | Yes |
| Naive lowest-three-orbital active space | 3 | 4 | 18 | 53 | 19.405 | No |
| Freeze core + parity + drop π pair | **4** | 8 | 72 | 162 | 1.252 | Yes |

The 6-qubit reduction retains accuracy with far fewer CNOTs; the 4-qubit active-space model remains just inside chemical accuracy. The choice of orbitals matters: the naive lowest-three-orbital selection exceeds the accuracy target.

## Intermediate shared setup

Run notebooks with the repository root as the Jupyter working directory. Import the shared helper with:

```python
from intermediate_challenge.p1_setup import *

problem, mapper = reduced_lih(level="6q")  # "full" (12q), "6q" or "4q"
qubit_op = mapper.map(problem.hamiltonian.second_q_op())
ansatz = build_ansatz(problem, mapper)
```

## Requirements

Python 3.10+ with `qiskit`, `qiskit-nature`, `qiskit-algorithms`, `pyscf`, `numpy`, and `matplotlib`.

## Credits

The H₂ tutorial was written by Benjamin Tirado for the BasQ Qiskit Fall Fest 2026, organised by Basque Quantum (BasQ). Its text is licensed under CC BY 4.0 and its code under Apache 2.0.
