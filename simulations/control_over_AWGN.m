close all
T=9999; %total time instants (T-1)
Q=1;R=1;F=1;V=1;W=1;
alpha=2; %system parameter
S(T+1)=F; %assumption from paper 
P(1)=5; %assuming Po, here denoted as P(1)
SNRs = linspace(5, 1000,500);
for i_SNR = 1 : length(SNRs)  %starting SNR to ensure stability according to condition given in the paper
    SNR = SNRs(i_SNR);
    for t=(T+1):-1:2
        S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q; % controller ARE iterations
    end
   
    for t=1:1:T
        P(t+1)=(alpha^2)*P(t)*(1-(P(t)/(P(t)+V))*(SNR/(1+SNR)))+W;
    end
    sig(1)=(P(1)*V)/(P(1)+V);
    for t=1:1:(T+1)
%        sig(t)=P(t)*(1-(P(t)/(P(t)+V))*(SNR/(1+SNR))) + W;
%        sig(t)=((P(t)*V)/(P(t)+V));%*(SNR/(1+SNR));
        sig(t)=(P(t)-W)/(alpha^2);
    end
    J_opt(i_SNR)=0;
    for t=2:1:(T+1)
        J_opt(i_SNR)=J_opt(i_SNR) + trace(Q*sig(t)) + trace(S(t)*((alpha^2)*sig(t-1) + W - sig(t)));
    end
    J_opt(i_SNR) = (J_opt(i_SNR) + trace(Q*sig(1)) + trace( S(1) * (P(1) - sig(1))))/ (T+1);
%     temp=0;
%     for t=2:1:(T+1)
%         temp=temp+((Q-S(t))*sig(t)+(alpha^2)*sig(t-1)*S(t) + S(t)*W);
% %         temp=temp+S(t)*W;
%     end
%     J_opt(i_SNR)=(((Q-S(1))*sig(1)+S(1)*P(1))+temp) / (T+1); %optimal performance cost value
%     J_opt(i_SNR)=(S(1)*P(1)+temp) / (T+1);
% %%%%%%%%%simulating the system %%%%%%%%

    w=sqrt(W)*randn((T+1),1); 
    v=sqrt(V)*randn((T+1),1);
    x(1)=sqrt(P(1))*randn(1,1);
    s_bar_tilde=sqrt(1/SNR)*randn((T+1),1);
    for t=(T+1):-1:2
        %controller gain
        L(t-1)=(alpha*S(t))/(S(t)+R);
    end
    for t=1:1:(T+1)
        %kalman gain
        K(t)=(SNR/(1+SNR))*(P(t)/sqrt(P(t)+V));
    end
    %simulate the flow
    y(1)=x(1)+v(1);
    x_hat_t(1)=0;
    u(1)=0;
    x(2)=alpha*x(1)+u(1)+w(1);
%     y(2)=x(2)+v(2);
    sumac=0;
    for t=2:1:(T+1)
        y(t)=x(t)+v(t);
        s(t)=y(t)-alpha*x_hat_t(t-1)-u(t-1);
        s_bar(t)=sqrt(1/(P(t)+V))*s(t);   
%         x_tilde=P(t)*randn(T,1);
%         y(t)=x(t)+v(t);
        %in another loop?
%         s(t)=x_tilde(t)+v(t);
%         s_bar(t)=(1/sqrt(P(t)+V))*s(t);
        y_tilde(t)=s_bar(t)-s_bar_tilde(t);
%         if t>=2
            %generate control signal
        x_hat(t)=alpha*x_hat_t(t-1)+u(t-1);
        x_hat_t(t)=x_hat(t)+K(t)*y_tilde(t);
        if t<(T+1)
            u(t)=-L(t)*x_hat_t(t);
            x(t+1)=alpha*x(t)+u(t)+w(t);
%             y(t+1)=x(t+1)+v(t+1);
            sumac = sumac + (Q*(x(t)^2)) + (R*(u(t))^2);
        end
    end
    J_opt_sim(i_SNR) = (  F*(x(T+1)^2) + sumac + (Q*(x(1)^2)) + (R*(u(1))^2) ) / (T+1);
%     J_opt_sim(SNR)=J_opt_sim(SNR)/(T-1);
%     J_opt_sim(SNR)=mean(J_opt_sim(SNR));
end
%calculating performance cost
figure 
hold on
title('Optimal cost vs SNR : Control over AWGN channel');
xlabel('SNR (in dB)');
ylabel('Performance cost');
scatter(10 * log10(SNRs), J_opt,'+');
scatter(10 * log10(SNRs), J_opt_sim);
% scatter(SNRs, J_opt,'+');
% scatter(SNRs, J_opt_sim);
legend('Analytically computed','Simulated System');
hold off

% figure
% 
% scatter(SNRs, J_opt-J_opt_sim,'.');