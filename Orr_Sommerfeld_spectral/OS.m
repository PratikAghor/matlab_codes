%*************************************************************************%
%Written by Pratik Aghor, Works only for even n
%Chebyeshev spectral collocation method to find eigenvalues of the Orr-Sommerfeld equation
%*************************************************************************%
Re=1e3;
ReInverse=1/Re;
zi=sqrt(-1);

n=100;
[D0,D1,D2,D4,y]=getD(n);
%*************************************************************************%
%*************************************************************************% 
%temporal problem, lambda1 and lambda3, wavenumbers in x and z direction
%are real and given. To find: unknown, complex omega

lambda1=1.0;%1.05886186;
lambda3=0.0;%0.283656;
k2=(lambda1*lambda1+lambda3*lambda3);
%*************************************************************************%
Z=zeros(n+1,n+1); %(1-y^2) on the diagonals
Y=zeros(n+1,n+1);%y on the diagonals
for k=1:n+1
    Z(k,k)=1-y(k)*y(k);
    Y(k,k)=y(k);
end

prompt = 'What is the base flow? Enter 0 for channel flow and 1 for Couette flow \n';
dummy = input(prompt);
if (dummy==0)
    U=Z;
    DU=-2*Y;
    D2U=-2*eye(n+1);
ax.imagmin=-1;
ax.imagmax=0.1;
ax.realmin=0;
ax.realmax=1;
end

if (dummy==1)
    U=Y;
    DU=eye(n+1);
    D2U=zeros(n+1,n+1);
ax.imagmin=-1;
ax.imagmax=0.1;
ax.realmin=-1;
ax.realmax=1;
end

%*************************************************************************%
A11=zeros(n+1,n+1); A12=zeros(n+1,n+1); A13=zeros(n+1,n+1);
A21=zeros(n+1,n+1); A22=zeros(n+1,n+1); A23=zeros(n+1,n+1);
A31=zeros(n+1,n+1); A32=zeros(n+1,n+1); A33=zeros(n+1,n+1);
B11=zeros(n+1,n+1); B12=zeros(n+1,n+1); B13=zeros(n+1,n+1);
B21=zeros(n+1,n+1); B22=zeros(n+1,n+1); B23=zeros(n+1,n+1);
B31=zeros(n+1,n+1); B32=zeros(n+1,n+1); B33=zeros(n+1,n+1);
%*************************************************************************%
% A11=ReInverse*(D4-2*k2*D2+k2*k2*eye(n+1))-zi*lambda1*(U*(D2-k2*eye(n+1))-D2*U);
A11=ReInverse*(D4-2*k2*D2+k2*k2*D0)+zi*lambda1*(-U*(D2-k2*D0)+D2U);

%include boundary conditions
A11(1,2:n+1)=0.0;
A11(n+1,1:n)=0.0;

A11(1,1)=1.0;
A11(n+1,n+1)=1.0;

A11(2,1:n+1)=D1(1,:);
A11(n,1:n+1)=D1(n+1,:);
%*************************************************************************%
A21=-zi*lambda3*DU;
% A21=zi*lambda3*(-2*Y);

%include boundary conditions
A21(1,1:n+1)=0.0; 
A21(n+1,1:n+1)=0.0;
%*************************************************************************%
% A22=zi*lambda3*(ReInverse*(D2-k2*eye(n+1))-zi*lambda1*U);
A22=zi*lambda3*(ReInverse*(D2-k2*D0)-zi*lambda1*U);

%include boundary conditions
A22(1,2:n+1)=0.0;
A22(n+1,1:n)=0.0;

A22(1,1)=1.0;
A22(n+1,n+1)=1.0;
%*************************************************************************%
% A23=-zi*lambda1*(ReInverse*(D2-k2*eye(n+1))-zi*lambda1*U);
A23=-zi*lambda1*(ReInverse*(D2-k2*D0)-zi*lambda1*U);

%include boundary conditions
A23(1,1:n+1)=0.0;
A23(n+1,1:n+1)=0.0;
%*************************************************************************%
A31=D1;
%include boundary conditions
A31(1,1:n+1)=0.0;
A31(n+1,1:n+1)=0.0;
%*************************************************************************%
% A32=zi*lambda1*eye(n+1);
A32=zi*lambda1*D0;

%include boundary conditions
A32(1,1:n+1)=0.0;
A32(n+1,1:n+1)=0.0;
%*************************************************************************%
% A33=zi*lambda3*eye(n+1);
A33=zi*lambda3*D0;

%include boundary conditions
A33(1,2:n+1)=0.0;
A33(n+1,1:n)=0.0;

A33(1,1)=1.0;
A33(n+1,n+1)=1.0;
%*************************************************************************%
% B11=-zi*(D2-k2*eye(n+1));
B11=-zi*(D2-k2*D0);

%include boundary conditions
B11(1,1:n+1)=0;
B11(n+1,1:n+1)=0;
B11(2,1:n+1)=0;
B11(n,1:n+1)=0;
%*************************************************************************%
% B22=-zi*lambda3*(zi*eye(n+1));
B22=-zi*(zi*lambda3*D0);

%include boundary conditions
B22(1,1:n+1)=0;
B22(n+1,1:n+1)=0;
%*************************************************************************%
% B23=zi*lambda1*(zi*eye(n+1));
B23=zi*(zi*lambda1*D0);
%include boundary conditions
B23(1,1:n+1)=0;
B23(n+1,1:n+1)=0;
%*************************************************************************%
%build matrices
A=[A11,A12,A13;A21,A22,A23;A31,A32,A33];
B=[B11,B12,B13;B21,B22,B23;B31,B32,B33];
%*************************************************************************%
[a_vecs,c_vals]=find_and_sort(A,B);       % a_vecs and c_vals defined pp 489 and 
%*************************************************************************%
%%%%% plotting
figure(1)
plot(real(c_vals),imag(c_vals),'o')
xlabel('\omega_r'); ylabel('\omega_i')
axis([ax.realmin ax.realmax ax.imagmin ax.imagmax])
grid on
title(['eigenmodes of the OSS-operator for lambda1= ',num2str(lambda1),', lambda3= '...
       ,num2str(lambda3),', Re = ',num2str(Re)])

%*************************************************************************%     
     
