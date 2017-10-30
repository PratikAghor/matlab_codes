%DuFort-Frankel scheme for 1d heat conduction
clear;clc;
L=1;
k=1;
nx=20;
xi=0;xo=L;
x=linspace(xi,xo,nx);
dx=x(nx)-x(nx-1);
dt=1e-3;
t=0.0;
gamma=k*dt/(dx*dx);
epsilon=3.125e-2;
%*************************************************************************%
%initial condition
%*************************************************************************%
u=zeros(nx,1);
for j=1:nx
    u0(j)=1.0;
end
u0(1)=1-epsilon;
u0(nx)=1+epsilon;
%*************************************************************************%
%get the computational initial condition
%*************************************************************************%
for j=2:nx-1
    u1(j)=u0(j)+gamma*(u0(j-1)-2.0*u0(j)+u0(j+1));
end
u1(1)=1-epsilon;
u1(nx)=1+epsilon;

for n=1:70
    n
u2(1)=1-epsilon;
u2(nx)=1+epsilon;
   
    t=t+dt;
    
    %DuFort-Frankel
    for j=2:nx-1
        u2(j)=(1.0/(1+2.0*gamma))*(u0(j)+2.0*gamma*(u1(j-1)-u0(j)+u1(j+1)));
    end
    
%*************************************************************************%
%update
%*************************************************************************%
u0temp=u0; u1temp=u1;
u1=u2; u0=u1temp;

if (rem(n,5)==0)
    plot(u2(1,:),x(1,:)); hold on
end
end
%*************************************************************************%
