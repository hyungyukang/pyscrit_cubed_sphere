%-----------------------------------------------------------------------
%
%  Legendre-Gauss-Lobatto (LGL) points and weights 
%
%  JR.PARK
%  26 Sep. 2012
%
%  input : N = # of gridpoints
%  output: xg,wg
%          xg : gaussian points
%          wg : gaussian weights at each gaussian point 
%  
%-----------------------------------------------------------------------

function [xg, wg] = legendre_gauss_lobatto(a,b,N)
 
apb= (a+b)/2;
bma= (b-a)/2;

nh= floor(N+1)/2;
Np= N+1;

xg= zeros(1,N);
wg= zeros(1,N);
xx=0;

    xg(1)= -1;wg(1)=2/(N*(N-1));
    xg(N)= +1;wg(N)=2/(N*(N-1));
for i=2:N-1;
    %x= cos( (2*i-1)*pi/(2*N+1) );
    x=(1-(3*(N-2))/(8*(N-1)^3)) * cos(pi*(4*i-3)/(4*(N-1)+1));

    for j=1:100;
       P0=1; P1=0;
       dp=1;dpp=1;dppp=1;
    for n=1:N;
       P2=P1;
       P1=P0;
       dp1= dp; dpp1= dpp; dppp1= dppp;
      
       P0= (2*n-1)*x*P1/n - (n-1)*P2/n;
       dp= n*(P1-x*P0)/(1-x*x); % the same with pp= (n/(x*x-1))*(x*P0-P1)
      dpp= ( 2*x*dp - n*(n+1)*P0   )/(1-x*x); 
     dppp= ( 2*x*dpp-(n*(n+1)-2)*dp)/(1-x*x);
    end%n  
           xx= x-(2*dp1*dpp1)/(2*dpp1*dpp1 - dp1*dppp1);
    if(abs(xx-x)<10.E-15)
       break
    end           
       x= xx;
    end%j
    
       xg(Np-i)= x*bma+apb;
       xg(   i)= -xg(Np-i);
       wg(Np-i)= 2/(N*(N-1)*P1*P1);
       wg(   i)=  wg(Np-i);
end%i

    x=xg;
    w=wg;

