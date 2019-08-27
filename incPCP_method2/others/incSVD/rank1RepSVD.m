function[U S V ] = rank1RepSVD(Uo, So, Vo, k, curFrame, newFrame)
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

rQ = Vo(k,:); rQ = rQ(:);
z = -Vo*rQ; z(k) = z(k)+1;

rhoQ = sqrt( 1 - rQ'*rQ );

if( rhoQ > 1e-8 )
  q = z/rhoQ;
else
  q = zeros(size(z));
end

% -- U --

rP = Uo'*newFrame;  % - So*rQ;
z = newFrame - Uo*rP;
rP = rP - So*rQ;    % OOO computation

rhoP = sqrt( sum( z.*z ) );

if( rhoP > 1e-8 )
  p = z/rhoP;
else
  p = zeros(size(z));
end


St = [ So+rP*rQ'  rhoQ*rP; rhoP*rQ' rhoP*rhoQ ];


% this is cheap since Zp is really small
[Gu S1 Gv] = svd(St);     % flops Ncols^3


if 1

  % rank is kept

  S  = S1(1:Ncols, 1:Ncols);
  U  = Uo*Gu(1:Ncols,1:Ncols) + p*Gu(Ncols+1,1:Ncols);
  V  = Vo*Gv(1:Ncols,1:Ncols) + q*Gv(Ncols+1,1:Ncols);


%%%%%%%%%%
else            %NOTE: Section 4.1, 4.2 of Brand, 'Fast low-rank mod' doesn't apply for partialSVD
%%%%%%%%%%


end

return

