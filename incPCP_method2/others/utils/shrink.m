function u = shrink(v, lambda)
    
  u = sign(v).*max(0, abs(v) - lambda);
  
return 

