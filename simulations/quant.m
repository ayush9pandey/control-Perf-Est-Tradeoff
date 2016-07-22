function num_q=quant(num,delta)
    % mid tread uniform quantizer
%       num_q = delta * round( (num/delta));
      num_q = delta * round( (num/delta) + 0.5);

end
