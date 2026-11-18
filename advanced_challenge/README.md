# Advanced challenge: quantum-selected configuration interaction

Open `h4_qsci.ipynb` for a reproducible QSCI demonstration on square H₄ at a stretched 2.0 Å side length. It samples computational-basis configurations from a shallow UCCSD state, selects configurations in the correct electron/spin sector, and classically diagonalizes the Hamiltonian in that selected subspace.

The notebook compares the unoptimized UCCSD trial state, the selected-subspace energy, and a particle-sector-filtered exact reference. This small STO-3G example is intended to explain the hybrid method, not demonstrate large-scale advantage: it uses a noiseless statevector for sampling and a dense Hamiltonian matrix. It explicitly reports the selected basis size against the full sector dimension.
