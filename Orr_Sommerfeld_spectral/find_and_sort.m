function [xs,es]= find_and_sort(A,B)
%takes matrices A and B, solves the generalized eigenvalue problem
%Ax=omega*Bx and then sorts the eigenvalues according to decreasing
%imaginary parts
[v,omega]=eig(A,B,'qz');

%sort the eigenvalues according to decreasing imaginary parts
%taken from Schmidt and Henningson's spectral code: function 'iord2.m'
    omega=diag(omega);
    [eimag,is]=sort(-imag(omega));
    xs=v(:,is); 
    es=omega(is);
