function[U S V thresh] = rank1DwnSVD(Uo, So, Vo, k)
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
N = max(size(Vo));

% -- V --

r = Vo(k,:); r = r(:);
z = -Vo*r; z(k) = z(k)+1;

rho = sqrt( 1 - r'*r );

if( rho > 1e-8 )
  q = z/rho;
else
  q = zeros(size(z));
end


St = [ So-So*r*r'  -rho*So*r; zeros(1,Ncols+1) ];


% this is cheap since Zp is really small
[Gu S1 Gv] = svd(St);     % flops Ncols^3


if 1

  % rank is kept

  S  = S1(1:Ncols, 1:Ncols);
  U  = Uo*Gu(1:Ncols,1:Ncols);
  Vtmp  = Vo*Gv(1:Ncols,1:Ncols) + q*Gv(Ncols+1,1:Ncols); 
  thresh = sqrt( sum(Vtmp(k,:).*Vtmp(k,:)) );                
  if k<N
    if k > 1
      V  = [Vtmp(1:k-1,:); Vtmp(k+1:N,:)];
    else
      V  = Vtmp(2:N,:);
    end
  else
    V  = Vtmp(1:N-1,:);
  end

%%%%%%%%%%
else            %NOTE: Section 4.1, 4.2 of Brand, 'Fast low-rank mod' doesn't apply for partialSVD
%%%%%%%%%%


end

return

