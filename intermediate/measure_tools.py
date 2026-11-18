"""Shot-based energy estimation for the Intermediate challenge.

Instead of an exact (infinite-shot) estimator, we simulate a real device: the Hamiltonian's Pauli
terms are split into qubit-wise commuting groups, each group is measured in one basis with a share
of the shot budget, and every shot is counted.

Speed: everything that doesn't depend on the circuit parameters (basis rotations and the ±1 sign
of every term for every possible outcome) is worked out once, so each energy evaluation is just
one statevector simulation plus a few small matrix products.
"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit import transpile
from qiskit.quantum_info import Statevector, Operator

SIM_GATES = ["cx", "rz", "ry", "rx", "sx", "x", "h", "s", "sdg"]


def measurement_groups(qubit_op):
    """Split the operator into qubit-wise commuting groups (each group = one measured circuit)."""
    return qubit_op.group_commuting(qubit_wise=True)


def _basis_change(group, n):
    """Circuit that rotates every qubit into the basis this group needs (X -> H, Y -> S† H)."""
    qc = QuantumCircuit(n)
    labels = group.paulis.to_labels()
    for q in range(n):
        idx = n - 1 - q                      # Qiskit labels are written q_{n-1} ... q_0
        ps = {l[idx] for l in labels} - {"I"}
        if ps == {"X"}:
            qc.h(q)
        elif ps == {"Y"}:
            qc.sdg(q)
            qc.h(q)
    return qc


def _sign_matrix(labels, n):
    """signs[t, o] = ±1: value of term t's Z-string on measurement outcome o (an integer bitstring)."""
    outcomes = np.arange(2 ** n)
    signs = np.empty((len(labels), 2 ** n))
    for t, l in enumerate(labels):
        mask = sum(1 << (n - 1 - i) for i, ch in enumerate(l) if ch != "I")
        parity = np.array([bin(o & mask).count("1") & 1 for o in outcomes])
        signs[t] = 1 - 2 * parity
    return signs


class ShotEstimator:
    """Estimate <H> with a fixed shot budget per energy evaluation, and count every shot.

    The budget is split equally between the measurement groups.
    """

    def __init__(self, ansatz, qubit_op, shots_per_eval, seed=None):
        self.ansatz = transpile(ansatz, basis_gates=SIM_GATES, optimization_level=1)  # simulate fast
        n = qubit_op.num_qubits
        self.identity = 0.0
        self.rotations, self.weights = [], []
        for g in measurement_groups(qubit_op):
            labels = g.paulis.to_labels()
            coeffs = g.coeffs.real
            keep = [i for i, l in enumerate(labels) if set(l) != {"I"}]
            self.identity += coeffs.sum() - coeffs[keep].sum()   # identity needs no measurement
            if not keep:
                continue
            self.rotations.append(Operator(_basis_change(g, n)).data)
            # coeffs @ signs: the energy contribution of each possible outcome, per shot
            self.weights.append(coeffs[keep] @ _sign_matrix([labels[i] for i in keep], n))
        self.shots = max(1, round(shots_per_eval / len(self.rotations)))
        self.rng = np.random.default_rng(seed)
        self.shots_used = 0
        self.evaluations = 0

    def energy(self, params):
        psi = Statevector(self.ansatz.assign_parameters(params)).data
        total = self.identity
        for U, w in zip(self.rotations, self.weights):
            probs = np.abs(U @ psi) ** 2
            counts = self.rng.multinomial(self.shots, probs / probs.sum())
            total += w @ counts / self.shots
        self.shots_used += self.shots * len(self.rotations)
        self.evaluations += 1
        return total


def exact_expectation(ansatz, qubit_op, params):
    """Noise-free <H> at given parameters (used to judge how good the optimised parameters are)."""
    circ = transpile(ansatz, basis_gates=SIM_GATES, optimization_level=1)
    return Statevector(circ.assign_parameters(params)).expectation_value(qubit_op).real
