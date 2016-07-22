sig = 1;
Deltas = linspace(0.1, 15, 100);
for i_delta = 1 : length(Deltas)
	delta = Deltas(i_delta);
	bound(i_delta) = delta^2 / 4 + sig^2 + (sig * delta)/sqrt(2 * pi) + 4 * sqrt(2/pi) * sig^ 3 / (3 * delta);
	summa = 0;
	for z = 1 : 10^5
	
		summa = summa + delta^2 * erfc( (delta / (2 * sqrt(2) * sig) ) * z)  * (z + 1)^2 / 4;
	
	end	
	val(i_delta) = summa;
		
end

hold on
plot(Deltas,val,'r');
plot(Deltas, bound);
legend('simulated', 'bound');
hold off

