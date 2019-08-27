function[U S V] = rank1IncSVD(Uo, So, Vo, curFrame, flag)
%  
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


[Nrows Ncols] = size(Uo);

% -- U --

r = Uo'*curFrame;
z = curFrame - Uo*r;
rho = sqrt( sum( z.*z ) );

if( rho > 1e-8 )  % NOTE: if Uo is gpuArray, so is p
  p = z/rho;
else
  p = zeros(size(z));
end


St = [ So r; zeros(1,Ncols) rho ];


% this is cheap since Zp is really small
[Gu S1 Gv] = svd(St);     % flops Ncols^3

% adaptively increase rank
dS1 = diag(S1);
rel = 100*dS1(end)/sum(dS1(1:end-1));

if rel > 100    % change '100' to an actual threshold to use this
  flag = 0;
end

if flag == 0 % increment rank
  S  = S1;
  U  = [ Uo p ]*Gu;
  V  = [ Vo zeros(max(size(Vo)),1);  zeros(1,Ncols) 1 ]*Gv;
else
  S  = S1(1:Ncols, 1:Ncols);
  U  = Uo*Gu(1:Ncols,1:Ncols) + p*Gu(Ncols+1,1:Ncols);
  V  = [ Vo*Gv(1:Ncols,1:Ncols); Gv(Ncols+1,1:Ncols) ];
end


return

