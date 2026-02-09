#!/usr/bin/env python3
# coding: utf-8

import itertools
import math
import sys
import numpy as np


def dagger(a: np.ndarray) -> np.ndarray:
    return a.conj().T


def kron_all(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def c2s(c: complex) -> str:
    if c == 0.0:
        return "0"
    if abs(c.imag) < 1e-15:
        return f"{c.real:g}"
    if abs(c.real) < 1e-15:
        return f"{c.imag:g}j"
    return f"{c.real:g}+{c.imag:g}j"


def decompose_pauli_string(H: np.ndarray, *, tol: float = 1e-12, max_qubits: int = 8) -> str:
    # --- sanity checks ---
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError(f"H must be a square matrix, got shape {H.shape}.")

    dim = H.shape[0]
    nbits = int(round(math.log2(dim)))
    if 2**nbits != dim:
        raise ValueError(f"Dimension {dim} is not a power of 2, so it can't be an n-qubit operator.")

    if nbits > max_qubits:
        print(f"Number of qubits: {nbits}")
        print(f"Refusing to decompose: >{max_qubits} qubits is exponentially hard (4^n Pauli terms).")
        sys.exit(1)

    if not np.allclose(H, dagger(H), atol=1e-10):
        print("Warning: H is not Hermitian (within tolerance). Decomposition still computed.")

    # --- Pauli basis ---
    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    id2 = np.array([[1, 0], [0, 1]], dtype=np.complex128)

    S = [id2, sx, sy, sz]
    labels = ["I", "X", "Y", "Z"]

    # --- Brute force decomposition (4^n terms) ---
    final_terms = []
    norm = 2**nbits

    for idxs in itertools.product(range(4), repeat=nbits):
        P = kron_all([S[i] for i in idxs])             # n-qubit Pauli string matrix
        a = (1.0 / norm) * np.trace(H @ P)             # coefficient

        if abs(a) > tol:
            pauli_str = "".join(labels[i] for i in idxs)
            final_terms.append(f"{c2s(a)}*{pauli_str}")

    return " + ".join(final_terms) if final_terms else "0"


if __name__ == "__main__":


    nbits = 6 
    if nbits > 8:
        print(f"Number of qubits: {nbits}")
        print("Refusing: >8 qubits is exponentially hard (4^n Pauli terms).")
        sys.exit(1)

    dim = 2**nbits
    rng = np.random.default_rng(0)

    # Random complex matrix -> make it Hermitian
    A = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    H = (A + A.conj().T) / 2.0

    Hps = decompose_pauli_string(H, tol=1e-12, max_qubits=8)
    print(Hps)
