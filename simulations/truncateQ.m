% function [num_trunc, flagO, flagS, decision] = truncateQ(num, delta, Ra, flagO, flagS)
function [num_trunc] = truncateQ(num, limit1, limit2, delta)
%     limit1 = abs(( (Ra - 1) * delta)) / 2 ; %for truncation
%     limit2 = -limit1;
if num > limit1
    num_trunc = limit1;
  elseif num < limit2
    num_trunc = limit2;
  else 
    num_trunc = delta * round( (num/delta) + 0.5);
  end
  
end
