from __future__ import print_function, division

import sys,os
qspin_path = os.path.join(os.getcwd(),"../")
sys.path.insert(0,qspin_path)

from quspin.basis import spin_basis_general
from quspin.basis.transformations import square_lattice_trans
from quspin.operators import hamiltonian
import numpy as np




def test(S):
	Lx = 3
	Ly = 3

	N = Lx*Ly

	tr = square_lattice_trans(Lx,Ly)

	nmax = int(eval("("+S+")*2"))

	basis_dict = {}
	Nup_list=range(nmax*N+1)

	for Nup in Nup_list:
		basis_dict[(Nup,None,None)] = spin_basis_general(N,Nup=Nup)
		for kx in range(Lx):
			for ky in range(Ly):
				basis_dict[(Nup,kx,ky)] = spin_basis_general(N,Nup=Nup,kxblock=(tr.T_x,kx),kyblock=(tr.T_y,ky))
				# print basis_dict[(Nup,kx,ky)]._n

	J = [[1.0,i,tr.T_x[i]] for i in range(N)]
	J.extend([[1.0,i,tr.T_y[i]] for i in range(N)])

	static = [["zz",J],["xx",J],["yy",J]]

	E_symm = {}

	for key,basis in basis_dict.items():
		H = hamiltonian(static,[],basis=basis,dtype=np.complex128)
		E_symm[key] = H.eigvalsh()

	for Nup in Nup_list:
		E_full_block = E_symm[(Nup,None,None)]
		E_block = []
		for kx in range(Lx):
			for ky in range(Ly):
				E_block.append(E_symm[(Nup,kx,ky)])

		E_block = np.hstack(E_block)
		E_block.sort()
		np.testing.assert_allclose(E_full_block,E_block,atol=1e-13)
		print("passed Nup={} sector".format(Nup))


test("1/2")
test("1")
