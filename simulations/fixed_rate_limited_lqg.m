close all

T=10^5 - 1; %total time instants (T-1)
Q=1;R=1;F=1;V=1;W=1;
alpha=2; %system parameter
S(T+1)=F; %assumption from paper 
P(1)=5; %assuming Po, here denoted as P(1)

w=sqrt(W)*randn((T+1),1); 
v=sqrt(V)*randn((T+1),1);
x(1)=sqrt(P(1))*randn(1,1);
s_bar_tilde=zeros((T+1),1);

%% Initial setting of rate and delta
Ra = 3;
Deltas = linspace(0,10,100);


for t=(T+1):-1:2
    S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q; % controller ARE iterations
end

for t=1:1:T
    P(t+1)=(alpha^2)*P(t)*(1-(P(t)/(P(t)+V)))+W;  % KF ARE
end

% %%%%%%%%%simulating the system %%%%%%%%


for t=(T+1):-1:2
    %controller gain
    L(t-1)=(alpha*S(t))/(S(t)+R);
end
for t=1:1:(T+1)
    %kalman gain
    K(t)=(P(t)/sqrt(P(t)+V));
end
%simulate the flow
%% For delta 
for i_delta = 1 : length(Deltas)
    delta = Deltas(i_delta);
%     count2 = 0; count3 = 0;

%     delta = 2; flagO = 0; flagS = 0; decision = 0;
    limit1 = ( (Ra - 1) * delta) / 2 ; %for truncation
    limit2 = -limit1;
    y(1) = x(1) + v(1);
    x_hat_t(1) = 0;
    u(1) = 0;
    x(2) = alpha*x(1)+u(1)+w(1);
    sumac=0;
      for t=2:1:(T+1)
        y(t) = x(t) + v(t);
        s(t) = y(t) - alpha * x_hat_t(t-1) - u(t-1);
        s_bar(t) = sqrt(1/(P(t)+V)) * s(t);
        %% Quantizer here 
        s_barQ(t) = truncateQ(s_bar(t), limit1, limit2, delta);
%         [s_barQ(t), flagO, flagS, decision] = truncateQ(s_bar(t), delta, Ra, flagO, flagS);
%         if decision == 1
%             % stretch delta 
%             delta = delta ^ (flagO + 1);
% %             count2 = count2 + 1;
%         else
%             % squeeze delta
%             delta = delta / flagS;
% %             count3 = count3 + 1;
%         end
%                        
        y_tilde(t) = s_barQ(t) - s_bar_tilde(t);
        x_hat(t) = alpha * x_hat_t(t-1) + u(t-1);
        x_hat_t(t) = x_hat(t) + K(t)*y_tilde(t);
%         x_hat_tQ(t) = quant(x_hat_t(t), delta);
        if t<(T+1)
            u(t) = -L(t)*x_hat_t(t);
            x(t+1) = alpha*x(t)+u(t)+w(t);
            sumac = sumac + (Q*(x(t)^2)) + (R*(u(t))^2);
        end         
      end
      %% cost calculation
    J_opt_sim(i_delta) = (  F*(x(T+1)^2) + sumac + (Q*(x(1)^2)) + (R*(u(1))^2) ) / (T+1);
    
end
% J_opt_sim = (  F*(x(T+1)^2) + sumac + (Q*(x(1)^2)) + (R*(u(1))^2) ) / (T+1);

[minJ,ind] = min(J_opt_sim);
minJ
delta_stable = Deltas(ind)