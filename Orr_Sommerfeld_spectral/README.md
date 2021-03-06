Spectral code for finding the eigenvalue spectrum of the Orr Sommerfeld equation.

OS.m 		=> Orr Sommerfeld equation in the primitive variables
OS_2.m		=> velocity-vorticity formulation, as given in Schmid and Henningson, chapter 3

getD.m		=> Create Chebyechev differentiation matrices
dummyDo.m	=> Create Chebyechev differentiation matrices	
find_and_sort.m	=> takes matrices A and B, solves the generalized eigenvalue problem Ax=omega*Bx and then sorts the eigenvalues according to decreasing imaginary parts

References:
1. Schmid, P.J. and Henningson, D.S., 2012. Stability and transition in shear flows (Vol. 142). Springer Science & Business Media.
2. Orszag, S.A., 1971. Accurate solution of the Orr–Sommerfeld stability equation. Journal of Fluid Mechanics, 50(4), pp.689-703.
3. Srinivasan, S., Klika, M., Ludwig, M.H. and Ram, V.V., 1994. A beginner's guide to the use of the spectral collocation method for solving some eigenvalue problems in fluid mechanics.
