"""Shared setup for the Intermediate challenge (LiH), owned by Person 1.

Teammates can do:
    from intermediate_challenge.p1_setup import *
    problem, mapper = reduced_lih(level="6q")   # or "4q", "full"
    qubit_op = mapper.map(problem.hamiltonian.second_q_op())

so everyone works on exactly the same molecule, exact energy and resource counting.
"""
import numpy as np
from qiskit import transpile
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper
from qiskit_nature.second_q.transformers import FreezeCoreTransformer, ActiveSpaceTransformer
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_algorithms import VQE, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import SLSQP
from qiskit.primitives import StatevectorEstimator

# ---- Agreed settings (change here and everyone picks it up) ----
BOND_LENGTH = 1.6          # Angstrom, close to LiH equilibrium
BASIS = "sto3g"
CHEM_ACC = 1.6e-3          # Hartree
BASIS_GATES = ["cx", "rz", "sx", "x"]


def lih_problem(bond_length=BOND_LENGTH, basis=BASIS):
    """Full LiH electronic-structure problem (6 spatial orbitals, 12 spin-orbitals)."""
    driver = PySCFDriver(
        atom=f"Li 0 0 0; H 0 0 {bond_length}",
        basis=basis, charge=0, spin=0, unit=DistanceUnit.ANGSTROM,
    )
    return driver.run()


def exact_energy(problem, mapper=None):
    """Exact total ground-state energy of `problem` (classical diagonalisation)."""
    mapper = mapper or JordanWignerMapper()
    solver = GroundStateEigensolver(mapper, NumPyMinimumEigensolver())
    return solver.solve(problem).total_energies[0]


def reduced_lih(level="6q", bond_length=BOND_LENGTH):
    """Return (problem, mapper) for a chosen reduction level.

    "full" : Jordan-Wigner, no reduction             -> 12 qubits
    "6q"   : freeze core + parity + tapering         ->  6 qubits (very accurate)
    "4q"   : "6q" + drop the two pi orbitals         ->  4 qubits (just inside chem. accuracy)
    """
    full = lih_problem(bond_length)
    if level == "full":
        return full, JordanWignerMapper()
    fc = FreezeCoreTransformer().transform(full)
    if level == "6q":
        return fc, fc.get_tapered_mapper(ParityMapper(num_particles=fc.num_particles))
    if level == "4q":
        # After freezing, orbitals are [sigma, sigma*, pi_x, pi_y, high sigma].
        # The pi pair (2, 3) cannot mix with the sigma bond by symmetry, so we drop it.
        act = ActiveSpaceTransformer(num_electrons=2, num_spatial_orbitals=3,
                                     active_orbitals=[0, 1, 4]).transform(fc)
        return act, act.get_tapered_mapper(ParityMapper(num_particles=act.num_particles))
    raise ValueError("level must be 'full', '6q' or '4q'")


def build_ansatz(problem, mapper):
    """Hartree-Fock + UCCSD ansatz for a (possibly reduced) problem."""
    hf = HartreeFock(problem.num_spatial_orbitals, problem.num_particles, mapper)
    return UCCSD(problem.num_spatial_orbitals, problem.num_particles, mapper, initial_state=hf)


def measure_resources(ansatz, qubit_op, label=""):
    """The numbers that go in the shared results table."""
    tc = transpile(ansatz, basis_gates=BASIS_GATES, optimization_level=1)
    return {
        "config": label,
        "qubits": qubit_op.num_qubits,
        "params": ansatz.num_parameters,
        "CNOTs": tc.count_ops().get("cx", 0),
        "depth": tc.depth(),
        "pauli_terms": len(qubit_op),
        "meas_groups": len(qubit_op.group_commuting(qubit_wise=True)),
    }


def run_vqe(problem, mapper, ansatz, optimizer=None):
    """Noiseless VQE; returns the total energy (nuclear repulsion and frozen core included)."""
    vqe = VQE(StatevectorEstimator(), ansatz, optimizer or SLSQP(maxiter=300))
    vqe.initial_point = np.zeros(ansatz.num_parameters)
    return GroundStateEigensolver(mapper, vqe).solve(problem).total_energies[0]
