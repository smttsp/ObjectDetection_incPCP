function[mask, tau] = unimodal(I, flag, cudaFlag)

if( nargin < 3)
  cudaFlag = 0;
  if( nargin < 2)
    flag = 0;
  end
end


vmaxI = max(I(:));
vminI = min(I(:));


if( abs(vmaxI - vminI) < 1e-2 )
  mask = 1.0;
  tau = 0.0;
  return;
end


if cudaFlag == 0
  [I_hist, pos] = hist(I(:),vminI:(vmaxI-vminI)/99:vmaxI);
else
  Ip = gather(I(:));
  [I_hist, pos] = hist(Ip(:),vminI:(vmaxI-vminI)/99:vmaxI);
end


[vmax pmax] = max(I_hist);      % no GPU

vend = I_hist(end);
%  pend = pos(end);
pend = length(pos);

m = (vend - vmax) / (pend - pmax);
alpha = pi/2 - atan(m);

k = pmax:pend;

dk = sqrt( (vmax-I_hist(pmax:pend)).^2 + (pmax - k).^2 );
mk = (I_hist(pmax:pend) - vmax) ./ (k-pmax);
alphak = pi/2 - atan(mk);

dpk = dk.*abs(sin(alpha-alphak));

%  figure; plot(dpk);

[d1 p1] = max(dpk);


if( (pmax + p1) <= length(pos) )
  tau = pos(pmax + p1);
else
  tau = pos(end);
end


if(flag > 0)

    tau = tau + flag*(pos(end)-tau);

end


if cudaFlag == 0
  mask = I >= tau;
else
  mask = I >= gpuArray(tau);
end
