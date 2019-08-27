function[xmin, fval] = FSearch(FHandle, a, b, tol)
%  
%  FHandel: handle to function to be minimized
%  a      : lower limit
%  b      : upper limit
%  tol    : tolerance
%  
% Authors
% =======
% 
% Paul Rodriguez   prodrig@pucp.pe
% Brendt Wohlberg  brendt@lanl.gov
% 
%  
% Legal
% =====
% 
%  There is currently a patent pending that covers the incremental PCP method and
%  applications that is implemented in this source code.
%  
%  For non-commercial or academic use the source code of this program can be distributed 
%  and/or modified under the terms of the GNU Lesser General Public License (LGPL) version 3 
%  as published by the Free Software Foundation (http://opensource.org/licenses/lgpl-3.0.html).
%  
%  For other usage, please contact the authors.
%  

if(b < a)
  disp('upper limit < lower limit --> contradiction');
  xmin = nan;
  return;
end

fibo = FNumber(a, b, tol);

an = a;
bn = b;

N = length(fibo);


% --------------------------
% --------------------------

  delta = bn - an;
  
  rho = fibo(N-1)/fibo(N);
  
  x1 = an + (1-rho)*delta;
  x2 = an + rho*delta;

  fx1 = feval(FHandle, x1);
  fx2 = feval(FHandle, x2);

% --------------------------
% --------------------------

for n=N-1:-1:3
  
  rho = fibo(n-1)/fibo(n);
  
  % ========================
  if( fx1 > fx2 ) % --------
  
    an  = x1;
    x1  = x2;
    fx1 = fx2;
    
    delta = bn - an;
    
    x2 = an + rho*delta;
    fx2 = feval(FHandle, x2);
  
  else % -------------------
  
    bn  = x2;
    x2  = x1;
    fx2 = fx1;
    
    delta = bn - an;
    
    x1  = an + (1-rho)*delta;
    fx1 = feval(FHandle, x1);
  
  end % --------------------
  % ========================
  
  if( delta < tol )
    break;
  end
  
end  % _END_ FOR(n)


xmin = (an + bn)/2;

if( nargout == 2)
  fval = feval(FHandle, xmin);
end

return;

% ======================================

function[F, epsilon] = FNumber(a, b, tol, epsilon)

if nargin < 4
  epsilon = tol/20;
end

F(1) = 1;  % 1st fibonacci number
F(2) = 2;  % 2nd fibonacci number
n=2;

delta = b - a;

rel = tol/delta;
alpha = (1-2*epsilon);

while( alpha/F(n) > rel )
    
  n = n + 1;
  F(n) = F(n-1)+F(n-2);
  
end

return;