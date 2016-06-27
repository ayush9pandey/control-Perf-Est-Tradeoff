function weird_num = weird_map(num)
    LSB = sign(num);

    bits = dec2bin( abs(num) );
    SGN = bits(end);
    bits = [bits(1:end-1), dec2bin((-LSB + 1)/2)];
    
    weird_num = bin2dec(bits) * (-1)^(num2str(SGN)); 
end