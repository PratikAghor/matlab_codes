%*************************************************************************%
%velocity-vorticity formulation, as given in Schmid and Henningson, chapter 3
%Chebyeshev spectral collocation method to find eigenvalues of the Orr-Sommerfeld equation
%*************************************************************************%
Re=1e3;
ReInverse=1/Re;
zi=sqrt(-1);

n=100;
[D0,D1,D2,D4,y]=getD(n);
% [D0,D1,D2,D4,y]=Dmat(n);
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
    u=diag(Z);
    Du=(D1*u);
    DU=diag(Du);
    D2u=(D2*u);
    D2U=diag(D2u);

%     DU=-2*Y;
%     D2U=-2*eye(n+1);
ax.imagmin=-1;
ax.realmin=0;
ax.realmax=1;
end

if (dummy==1)
    U=Y;
    u=diag(Y);
    Du=(D1*u);
    DU=diag(Du);
    D2u=(D2*u);
    D2U=diag(D2u);
    
%     DU=eye(n+1);
%     D2U=zeros(n+1,n+1);
ax.imagmin=-1;
ax.realmin=-1;
ax.realmax=1;
end

%*************************************************************************%
A11=zeros(n+1,n+1); A12=zeros(n+1,n+1); 
A21=zeros(n+1,n+1); A22=zeros(n+1,n+1); 
B11=zeros(n+1,n+1); B12=zeros(n+1,n+1);
B21=zeros(n+1,n+1); B22=zeros(n+1,n+1); 
%*************************************************************************%
A11=ReInverse*(D4-2*k2*D2+k2*k2*D0)-zi*lambda1*(-D2U+U*(D2-k2*D0));
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
A22=ReInverse*(D2-k2*D0)-zi*lambda1*U;

%include boundary conditions
A22(1,2:n+1)=0.0;
A22(n+1,1:n)=0.0;

A22(1,1)=1.0;
A22(n+1,n+1)=1.0;
%*************************************************************************%
B11=-zi*(D2-k2*D0);

%include boundary conditions
B11(1,:)=0;
B11(n+1,:)=0;
B11(2,:)=0;
B11(n,:)=0;
%*************************************************************************%
B22=-zi*D0;
%include boundary conditions
B22(1,:)=0;
B22(n+1,:)=0;
%*************************************************************************%
%build matrices
A=[A11,A12;A21,A22];
B=[B11,B12;B21,B22];
%*************************************************************************%
[a_vecs,c_vals]=find_and_sort(A,B);       % a_vecs and c_vals defined pp 489 and 
%*************************************************************************%
%%%%% plotting
figure(1)
plot(real(c_vals),imag(c_vals),'o')
xlabel('\omega_r'); ylabel('\omega_i')
axis([ax.realmin ax.realmax ax.imagmin 0.1])
grid on
title(['eigenmodes of the OSS-operator for lambda1= ',num2str(lambda1),', lambda3= '...
       ,num2str(lambda3),', Re = ',num2str(Re)])

%*************************************************************************%     
