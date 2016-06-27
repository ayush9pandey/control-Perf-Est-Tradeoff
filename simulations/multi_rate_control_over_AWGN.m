close all
clear all
clc
T=999; %total time instants (T-1)
s = randn(T+1,1);
SNR = 100; % 100 = 20 dB
n_1 = sqrt(1/SNR)*randn((T+1),1);
n_2 = sqrt(1/SNR)*randn((T+1),1);

deltas = linspace(.000001, 10,10000);
for i_delta = 1 : length(deltas)  %starting SNR to ensure stability according to condition given in the paper
   delta = deltas(i_delta);
   for t=1:1:(T+1)
       %encoder
       a_1(t) = quant(s(t),delta);
       a_2(t) = sqrt(12/(delta^2)) * (s(t) - a_1(t));
       %channel
       b_1(t) = a_1(t) + n_1(t);
       b_2(t) = a_2(t) + n_2(t);
      %decoder
       Q_b1(t) = quant(b_1(t),delta);
       s_hat(t) = Q_b1(t) + sqrt((delta^2)/12) * b_2(t);
       
       n_eff(t) = s(t) - s_hat(t);
   end
   n_eff_cov(i_delta) = cov(n_eff);
   
end
% 
% deltas_dB = linspace(-20, 20, 123);
% deltas    = 10.^(deltas_dB/20);
% 
% snr_dB = linspace(0, 20, 20);
% snr = 10.^(snr_dB/10);
% 
% n_samples = bitshift(1, 10);
% 
% s = randn(1, n_samples);
% [S, Deltas, SNR] = meshgrid(s, deltas, snr);
% 
% % Encoder
% A1 = Deltas .* round(S ./ Deltas);
% A2 = sqrt(12) ./ Deltas .* (S - A1);
% 
% % Channel
% N1 = 1./sqrt(SNR) .* randn(size(S));
% N2 = 1./sqrt(SNR) .* randn(size(S));
% B1 = A1 + N1;
% B2 = A2 + N2;
% 
% % Decoder
% Shat = Deltas .* round(B1 ./ Deltas) + Deltas / sqrt(12) .* B2;
% 
% % MSE calculation
% D = min( mean( (S - Shat).^2, 2), [], 1);
% D = reshape(D, 1, length(snr));
% sdr = 1 ./ D;
% sdr_dB = 10 * log10(sdr);
% 
% % OPTA
% sdr_OPTA = (1 + snr).^2 - 1;
% sdr_OPTA_dB = 10 * log10(sdr_OPTA);
% 
% % Plots
% figure
% hold on
% 	plot(snr_dB, sdr_OPTA_dB, 'bo');
% 	plot(snr_dB, sdr_dB, 'kp');
% 	xlabel('SNR [dB]')
% 	ylabel('SDR [dB]')
% 	legend('OPTA', 'Quantizer', 'location', 'best')
% hold off

%%% [Y,I] = min(n_eff_cov);
%%% delta_opt = deltas(I)
figure 
hold on
title('SDR vs Delta');
ylabel('SDR (in dB)');
xlabel('\Delta');
scatter(20*log10(deltas), -10 * log10(n_eff_cov),'p');
% scatter(10 * log10(SNRs), J_opt_sim);
% scatter(SNRs, J_opt,'+');
% scatter(SNRs, J_opt_sim);
% legend('Analytically computed','Simulated System');
hold off
