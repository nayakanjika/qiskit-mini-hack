"""Shot-based energy estimation for the Intermediate challenge.

Instead of an exact (infinite-shot) estimator, we simulate a real device: the Hamiltonian's Pauli
terms are split into qubit-wise commuting groups, each group is measured in one basis with a fixed
number of shots, and every shot is counted towards the budget.
"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


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


class ShotEstimator:
    """Estimate <H> with a fixed shot budget per energy evaluation, and count every shot.

    allocation="uniform":  the budget is split equally between measurement groups.
    allocation="weighted": groups with bigger coefficients get more shots (shots ∝ sum |c| in the
                           group), which lowers the statistical error for the same budget.
    """

    def __init__(self, ansatz, qubit_op, shots_per_eval, allocation="uniform", seed=None):
        self.ansatz = ansatz
        self.n = qubit_op.num_qubits
        groups = measurement_groups(qubit_op)
        self.identity = sum(g.coeffs.real.sum() for g in groups
                            if all(set(l) == {"I"} for l in g.paulis.to_labels()))
        self.groups = [g for g in groups
                       if not all(set(l) == {"I"} for l in g.paulis.to_labels())]
        weights = (np.ones(len(self.groups)) if allocation == "uniform"
                   else np.array([np.abs(g.coeffs.real).sum() for g in self.groups]))
        self.shots = np.maximum(1, np.round(shots_per_eval * weights / weights.sum())).astype(int)
        self.rng = np.random.default_rng(seed)
        self.rotations = [_basis_change(g, self.n) for g in self.groups]
        self.shots_used = 0
        self.evaluations = 0

    def energy(self, params):
        state = Statevector(self.ansatz.assign_parameters(params))
        total = self.identity                                # identity term needs no measurement
        for group, rot, shots in zip(self.groups, self.rotations, self.shots):
            probs = state.evolve(rot).probabilities()        # index = bitstring as integer
            counts = self.rng.multinomial(shots, probs / probs.sum())
            self.shots_used += shots
            outcomes = np.nonzero(counts)[0]
            for l, c in zip(group.paulis.to_labels(), group.coeffs.real):
                mask = sum(1 << (self.n - 1 - i) for i, ch in enumerate(l) if ch != "I")
                signs = np.array([(-1) ** bin(o & mask).count("1") for o in outcomes])
                total += c * (signs @ counts[outcomes]) / shots
        self.evaluations += 1
        return total


def exact_expectation(ansatz, qubit_op, params):
    """Noise-free <H> at given parameters (used to judge how good the optimised parameters are)."""
    return Statevector(ansatz.assign_parameters(params)).expectation_value(qubit_op).real
