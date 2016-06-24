close all
T=100; %total time instants
Q=1;R=1;F=1;V=1;W=1;
alpha=2; %system parameter
S(T)=F; %assumption from paper 
P(1)=5; %assuming Po, here denoted as P(1)
SNRs = linspace(10, 1000, 1000);
for i_SNR = 1 : length(SNRs)  %starting SNR to ensure stability according to condition given in the paper
    SNR = SNRs(i_SNR);
    for t=T:-1:2
        S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q; % controller ARE iterations
    end
   
    for t=1:1:T-1
        P(t+1)=((alpha^2)*P(t)*(1-(P(t)*i_SNR)/((P(t)+V)*(1+i_SNR))))+W; % KF gain calculation
    end
    temp=0;
    for t=2:1:T
        temp=temp+((Q-S(1))*P(t)+(alpha^2)*P(t-1)*S(t)+W*S(t));
    end
    J_opt(i_SNR)=(Q*P(1)+temp) / T; %optimal performance cost value
% end
% i_SNR=1:1:50;
% %plotting
% scatter(SNR,J_opt)
% %

% 
% %%%%%%%%%simulating the system %%%%%%%%
% for SNR=4:1:50
    w=sqrt(W)*randn(T,1); 
    v=sqrt(V)*randn(T,1);
%     n=N*randn(T,1);
    x(1)=sqrt(P(1))*randn(1,1);
    s_bar_tilde=sqrt(1/SNR)*randn(T,1);
%     x_hat(1)=0;
    for t=T:-1:2
        %controller gain
        L(t-1)=(alpha*S(t))/(S(t)+R);
    end
    for t=1:1:T
        %kalman gain
        K(t)=(SNR/(1+SNR))*(P(t)/sqrt(P(t)+V));
    end
    %simulate the flow
    y(1)=x(1)+v(1);
    x_hat_t(1)=0;
    u(1)=0;
    x(2)=alpha*x(1)+u(1)+w(1);
    y(2)=x(2)+v(2);
    sumac=0;
    for t=2:1:T
        
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
        if t<T
            u(t)=-L(t)*x_hat_t(t);
            x(t+1)=alpha*x(t)+u(t)+w(t);
            y(t+1)=x(t+1)+v(t+1);
            sumac = sumac + (Q*(x(t)^2)) + (R*(u(t))^2);
        end
    end
    J_opt_sim(i_SNR) = (  F*(x(T)^2) + sumac  ) / T;
%     J_opt_sim(SNR)=J_opt_sim(SNR)/(T-1);
%     J_opt_sim(SNR)=mean(J_opt_sim(SNR));
end
%calculating performance cost

hold on
% scatter(10 * log10(SNRs), 10 * log10(J_opt),'+');
% scatter(10 * log10(SNRs), 10 * log10(J_opt_sim));
scatter(SNRs, J_opt,'+');
scatter(SNRs, J_opt_sim);
hold off

figure

scatter(SNRs, J_opt-J_opt_sim,'.');