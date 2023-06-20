#!/usr/bin/env python
# coding: utf-8


import itertools
import numpy as np
import math 
from qiskit.opflow import X, Y, Z, I
 
def decompose(H):
    
    def dagger(a):
        return np.transpose(a).conj()

    def kronecker_product(matrices):

        result = matrices[0]
        for i in range(1, len(matrices)):
            result = np.kron(result, matrices[i])
        return result

    def c2s(c):

        if c == 0.0:
            return "0"
        if c.imag == 0:
            return "%g" % c.real
        elif c.real == 0:
            return "%gj" % c.imag
        else:
            return "%g+%gj" % (c.real, c.imag)
    
    sx = np.array([[0, 1],  [ 1, 0]], dtype=np.complex128)
    sy = np.array([[0, -1j],[1j, 0]], dtype=np.complex128)
    sz = np.array([[1, 0],  [0, -1]], dtype=np.complex128)
    id = np.array([[1, 0],  [ 0, 1]], dtype=np.complex128)
    S = [id, sx, sy, sz]
    labels = ["I", "X", "Y", "Z"]
    final = ""

    if nbits  == 3:

        for i, j, k in itertools.product(range(4),range(4),range(4)):

            matrices = [S[i], S[j], S[k]]
            tmp = kronecker_product(matrices)
            if np.allclose(tmp, dagger(tmp)) == False:
                print ("ERROR")

            a = (1/2**nbits) * np.trace(H@tmp)

            if a != 0:
                term = c2s(a),'*',labels[i],'^',labels[j],'^',labels[k]
                t2 = str(term).replace(',', '')
                t2 = str(t2).replace("'", "")
                
                if final == "":
                    final = "".join((final, t2))
                else: 
                    final = "+".join((final, t2))
    
    return final


H = np.loadtxt("ham_HO.txt")
#, dtype='f', delimiter='\t'
nbits = int(math.log2(np.shape(H)[0]))
print("Number of qubits:", nbits) 
Hps = decompose(H)
print (Hps) 
