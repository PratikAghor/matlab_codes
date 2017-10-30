% create D0
%Create Chebyeshev differentiation matrices
%*************************************************************************%
%Written by Pratik Aghor, Works only for even N
%purpose: Directly create D1, and from there get D2,D3,D4,etc.  
%*************************************************************************%
n=100;

D1=zeros(n+1,n+1);D2=zeros(n+1,n+1);D4=zeros(n+1,n+1);

% D0= [] ;
% vec=(0:1:n)' ;
% for j=0:1:n
%     D0=[D0 cos(j*pi*vec/n)];
% end;

if(rem(n,2)==0)
    
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
           D1(i,i)=-y(i)/(2*(1-y(i)*y(i)));
        end
    end
    D1(1,1)=(2*n*n+1)/6; D1(n+1,n+1)=-(2*n*n+1)/6;

    D2=D1*D1;
    D4=D1*D1*D1*D1;
end

%*************************************************************************%
%   N=2;
%   if N==0, D1=0; D2=0; D4=0; x=1; return, end
%   x = cos(pi*(0:N)/N)'; 
%   c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
%   X = repmat(x,1,N+1);
%   dX = X-X';                  
%   D1  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
%   D1  = D1 - diag(sum(D1'));               % diagonal entries
%   D2=D1*D1;
%   D4=D2*D2;

