%Create Chebyeshev differentiation matrices
%*************************************************************************%
%Written by Pratik Aghor, Works only for even N
%*************************************************************************%
function [D0,D1,D2,D4,y]=getD(N) 
% input: no. of collocation points
% output: Chebyeshev differentiation matrices
% n=round(abs(N));
n=round(N);
D1=zeros(n+1,n+1);D2=zeros(n+1,n+1);D4=zeros(n+1,n+1);

if(rem(n,2)==0)
    D0=eye(N+1);

    c=ones(1,n+1);
    c(1)=2.0; c(n+1)=2.0; 

    y=zeros(1,n+1);
    for i=1:n+1
        y(i)=cos(pi*(i-1)/n);
    end    

    for i=1:n+1
        for j=1:n+1
            if(j~=i)
           D1(i,j)=(c(i)/c(j))*((-1)^(i+j))/(y(i)-y(j));
            end
        end
           D1(i,i)=-y(i)/(2*(1-y(i)*y(i)));
    end
    D1(1,1)=(2*n*n+1)/6; D1(n+1,n+1)=-(2*n*n+1)/6;

    D2=D1*D1;
    D4=D1*D1*D1*D1;
end
%*************************************************************************%
%taken from Treften's book: spectral methods in MATLAB, works for bothe
%even and odd N's
%*************************************************************************%
%   function [D0,D1,D2,D4,x] = getD(N)
%   if N==0, D0=0; D1=0; D2=0; D4=0; x=1; return, end
%   x = cos(pi*(0:N)/N)'; 
%   c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
%   X = repmat(x,1,N+1);
%   dX = X-X';           
%   
%     D0=eye(N+1);
% 
%   
%   D1  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
%   D1  = D1 - diag(sum(D1'));               % diagonal entries
%   D2=D1*D1;
%   D4=D1*D1*D1*D1;
%*************************************************************************%
