function [num_trunc] = truncate(num, limit1, limit2)
  if num > limit1
    num_trunc = limit1;
%     count = count + 1;
  elseif num < limit2
    num_trunc = limit2;
%     count = count + 1;
  else 
    num_trunc = num;
  end
end
