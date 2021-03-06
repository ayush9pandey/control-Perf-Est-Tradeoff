n = 1; % dimensions, also for rate matched case n= k
ip = 1; % number of inputs
T = 10^6 - 1; %time samples
Q = 1; B = 1; R = 1; V = 100; C = 1; F = 1;
S(T+1)=F;
% Cost = linspace(5,50,100);
SNR_dB = linspace(6,10,3);
SNRs = 10.^(SNR_dB/10);
% SNRs = linspace(8, 60,400);
% [w,b] = meshgrid(w_t, Cost);
A = 2; P(1) = 5; % system parameter and initial condition
%laplacian distribution for plant disturbance
W = 1; %power of disturbance signal
Wb = sqrt(W/2);
e = exp(1);
%% Gaussian Distribution
% 
w = sqrt(W) * randn((T+1),1); 
hw = log(sqrt(W) * sqrt(2*pi*e));

%% Laplacian Distribution
% w = laprnd((T+1), 1, 0, sqrt(W));
% hw = 0.5 * log(2*W*e^2);

%% Uniform Distribution
% lim1 = 0;
% lim2 = 1;
% w = lim1 + (lim2-lim1).*rand((T+1),1);
% hw = log(lim2-lim1);

Nw = (1/(2*pi*e)) * exp((2*hw)/n);

v = sqrt(V) * randn((T+1),1);


x(1) = sqrt(P(1))*randn(1,1);


%finding the theoretical bound to the minimum achievable cost



[Sd,~,~] = dare(A,B,Q,R,zeros(n,ip),eye(n));
M = Sd*B*((R+B'*Sd*B)\B')*Sd;
[Td,~,~] = dare(A',C',W,V,0,eye(n));
Kd = Td*C'/(C*Td*C' + V);
Pd = (eye(n) - Kd*C)*Td;
bminw = trace(W*Sd);
bminv = trace(Pd*A'*M*A);
bmin = bminw + bminv;

% nwv = (det(W*M))^1/n;
%computing converse bound
% 
% for i_cost = 1 : length(Cost)
%     b = Cost(i_cost);   
%     Rb = log(abs(det(A))) + (n/2) * log(1 + (Nw * (abs(det(M))^(1/n))) / ((b - bmin)/n));
%     snr(i_cost) = exp((2/n) * Rb) - 1;
% %     snr(i_cost) = double(snr(i_cost));
% end

%Simulation 


%Controller 
for t=(T+1):-1:2
    S(t-1)=((A^2)*R*S(t))/(S(t)+R) + Q; % controller ARE iterations
end
for t = (T+1):-1:2
    %controller gain
    L(t-1) = (A*S(t))/(S(t)+R);
end
for i_SNR = 1 : length(SNRs)  %starting SNR to ensure stability according to condition given in the paper
    SNR = SNRs(i_SNR);
% %     Rb = log(det(A)) + (n/2) * log(1 + (Nw * abs(det(M))^(1/n)) / ((b - bmin)/n));
    Rb = (n/2) * log(1+SNR);
    
    %% Fully observable
%     b(i_SNR) = bmin + (n * (Nw * (abs(det(M))^(1/n)))) / (exp( (2/n) * (Rb - log(abs(det(A))))) - 1);

    %% Partially Observable
    fact = (((abs(det(C * (Td - W) * C' + V))^(1/n)) + (abs(det(C'*C))^(1/n))*Nw)*(abs(det(Kd * M))^(1/n)))/(exp( (2/n) * (Rb - log(abs(det(A)))))-1);
    b(i_SNR) = bmin + n * fact;
%     
    %%
% %     
    for t=1:1:T
        P(t+1)=(A^2)*P(t)*(1-(P(t)/(P(t)+V))*(SNR/(1+SNR)))+W;
    end
    
% %%%%%%%%%simulating the system %%%%%%%%
 
    s_bar_tilde = sqrt(1/SNR)*randn((T+1),1);
    for t = 1:1:(T+1)
        %kalman gain
        K(t) = (SNR/(1+SNR))*(P(t)/sqrt(P(t)+V));
    end
    %simulate the flow
    y(1) = x(1)+v(1);
    x_hat_t(1) = 0;
    u(1) = 0;
    x(2) = A*x(1) + u(1) + w(1);
%     y(2)=x(2)+v(2);
    sumac = 0;
    for t = 2:1:(T+1)
        y(t) = x(t)+v(t);
        s(t) = y(t)-A*x_hat_t(t-1)-u(t-1);
        s_bar(t) = sqrt(1/(P(t)+V))*s(t);   
%         x_tilde=P(t)*randn(T,1);
%         y(t)=x(t)+v(t);
        %in another loop?
%         s(t)=x_tilde(t)+v(t);
%         s_bar(t)=(1/sqrt(P(t)+V))*s(t);
        y_tilde(t)=s_bar(t)-s_bar_tilde(t);
%         if t>=2
            %generate control signal
        x_hat(t)=A*x_hat_t(t-1)+u(t-1);
        x_hat_t(t)=x_hat(t)+K(t)*y_tilde(t);
        if t<(T+1)
            u(t)=-L(t)*x_hat_t(t);
            x(t+1)=A*x(t)+u(t)+w(t);
%             y(t+1)=x(t+1)+v(t+1);
            sumac = sumac + (Q*(x(t)^2)) + (R*(u(t))^2);
        end
    end
    J_opt_sim(i_SNR) = (  F*(x(T+1)^2) + sumac + (Q*(x(1)^2)) + (R*(u(1))^2) ) / (T+1);
%     J_opt_sim(SNR)=J_opt_sim(SNR)/(T-1);
%     J_opt_sim(SNR)=mean(J_opt_sim(SNR));
end

%% Plots

figure 
hold on
% title('Bound on best achievable minimum cost vs SNR  : Control over AWGN channel with Non Gaussian Disturbances');
ylabel('SNR (in dB)');
xlabel('Performance cost');
% plot(10*log10(snr), Cost,'kp');
plot(b,10*log10(SNRs),'r');
plot(J_opt_sim, 10*log10(SNRs),'--');
line([bmin bmin], [0 max(J_opt_sim)]);
legend('Computed lower bound','Simulated System');
hold off