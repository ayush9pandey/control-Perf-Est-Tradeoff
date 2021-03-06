close all

% T=10^5 - 1; %total time instants (T-1)
T = 1:10^5;
Q=1;R=1;F=1;V=0;W=1;
alpha=2; %system parameter
S(T+1)=F; %assumption from paper 

P = zeros(length(T)+1,1);
K = zeros(length(T),1);
s = zeros(length(T),1);
y = zeros(length(T),1);
s_bar = zeros(length(T),1);
s_barQ = zeros(length(T),1);
y_tilde = zeros(length(T),1);
x_hat = zeros(length(T),1);
x_hat_t = zeros(length(T),1);
u = zeros(length(T)-1,1);
sumac = zeros(length(T),1);
x = zeros(length(T),1);


P(1)=5; %assuming Po, here denoted as P(1)

w=sqrt(W)*randn(max(T),1); 
v=sqrt(V)*randn(max(T),1);
x(1)=sqrt(P(1))*randn(1,1);
s_bar_tilde=zeros(max(T),1);


%% Initial setting of rate and delta
Ra = 3; %choose a minimum rate
count2 = 0; count3 = 0;

delta = 1.2; flagO = 0; flagS = 0; decision = 0;
limit1 = ( (Ra - 1) * delta) / 2 ; %for truncation
limit2 = -limit1;

para1 = 2.2; % To stretch
para2 = 0.75; % To squeeze

for t=length(T):-1:2
    S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q; % controller ARE iterations
end

for t=1:length(T)
    P(t+1)=(alpha^2)*P(t)*(1-(P(t)/(P(t)+V)))+W;  % KF ARE
end

% %%%%%%%%%simulating the system %%%%%%%%


for t=length(T):-1:2
    %controller gain
    L(t-1)=(alpha*S(t))/(S(t)+R);
end
for t=1:1:length(T)
    %kalman gain
    K(t)=(P(t)/sqrt(P(t)+V));
end
%simulate the flow
y(1)=x(1)+v(1);
x_hat_t(1)=0;
u(1)=0;
x(2)=alpha*x(1)+u(1)+w(1);
sumac(1) = 0;
reverse_str = '';
for t=2:1:length(T)
    
    y(t) = x(t) + v(t);
    s(t) = y(t) - alpha * x_hat_t(t-1) - u(t-1);
    s_bar(t) = sqrt(1/(P(t)+V)) * s(t);
%% Quantizer here    
    if s_bar(t) > limit1
        s_barQ(t) = limit1;
        count2 = count2 + 1;
        delta = delta * para1;
%         para1 = para1 + 1/(1*max(T));
        limit1 = ( (Ra - 1) * delta) / 2 ; %for truncation
        limit2 = -limit1;
        flagS = 1;
    elseif s_bar(t) < limit2
        s_barQ(t) = limit2;
        count2 = count2 + 1;
        delta = delta * para1;
%         para1 = para1 + 1/(1*max(T));
        limit1 = ( (Ra - 1) * delta) / 2 ; %for truncation
        limit2 = -limit1;
        flagS = 1;
    else 
        s_barQ(t) = delta * round( (s_bar(t)/delta) );
        count3 = count3 + 1;
        delta = delta * para2;
%         para2 = para2 - 1e-6;
        limit1 = ( (Ra - 1) * delta) / 2 ; %for truncation
        limit2 = -limit1;
        flagS = 0;
    end                
    
    
    y_tilde(t) = s_barQ(t) - s_bar_tilde(t);
    x_hat(t) = alpha * x_hat_t(t-1) + u(t-1);
    if flagS == 1
        x_hat_t(t) = x_hat(t) + K(t)*y_tilde(t);
    else
        x_hat_t(t) = x_hat(t) + K(t)*y_tilde(t);
    end
    if t<max(T)
%         u(t) = -L(t)*truncate(x_hat_t(t), limit1, limit2);
%         if flagS == 1
%             u(t) = -L(t)*x_hat_t(t)*1.2;
%         else
%             u(t) = -L(t)*x_hat_t(t);
%         end
        u(t) = -L(t)*x_hat_t(t);
        x(t+1) = alpha*x(t)+u(t)+w(t);
%         sumac = sumac + (Q*(x(t)^2)) + (R*(u(t))^2);
        sumac(t+1) = (Q*(x(t)^2)) + (R*(u(t))^2);
       
    end

    % Progress bar
    if ~mod(t, 1000)
	msg = sprintf('%3.1f%%%% iterations done', 100 * t / length(T)); % Maybe 2 '%' would do
	fprintf([reverse_str, msg]);
	reverse_str = repmat(sprintf('\b'), 1, length(msg)-1);
    end
    
end
disp(' ')
%% Optimal cost

% J_opt_sim = (  F*(x(T+1)^2) + sumac + (Q*(x(1)^2)) + (R*(u(1))^2) ) / (T+1);
J_opt_sim = ( sum(sumac) +  F*(x(max(T))^2) + (Q*(x(1)^2)) + (R*(u(1))^2) ) / length(T);
cumsumac = cumsum(sumac) ./ T';
J_opt_sim_dB = 10*log10(J_opt_sim)


%% Calculating lower bound to cost
e = exp(1);
n = 1; B = 1; C = 1; ip = 1;
hw = log(sqrt(W) * sqrt(2*pi*e));
Nw = (1/(2*pi*e)) * exp((2*hw)/n);

%% classical LQG cost
[Sd,~,~] = dare(alpha,B,Q,R,zeros(n,ip),eye(n));
M = Sd*B*((R+B'*Sd*B)\B')*Sd;
[Td,~,~] = dare(alpha',C',W,V,0,eye(n));
Kd = Td*C'/(C*Td*C' + V);
Pd = (eye(n) - Kd*C)*Td;
bminw = trace(W*Sd);
bminv = trace(Pd*alpha'*M*alpha);
bmin = bminw + bminv;

 %Fully observable
%     Rb(i_delta) = H_q(i_delta) + log2(e) + log2(H_q(i_delta) + 1);
b = bmin + (n * (Nw * (abs(det(M))^(1/n)))) / (exp( (2/n) * (log(Ra) - log(abs(det(alpha))))) - 1);
b_dB = 10*log10(b)
%     
    % Partially Observable
%     fact = (((abs(det(C * (Td - W) * C' + V))^(1/n)) + (abs(det(C'*C))^(1/n))*Nw)*(abs(det(Kd * M))^(1/n)))/(exp( (2/n) * (H_q(i_delta) - log(abs(det(A)))))-1);
%     b(i_delta) = bmin + n * fact;


% %% Without Adaptive Quantizer
% 
% for t=2:1:length(T)
%     y(t) = x(t) + v(t);
%     s(t) = y(t) - alpha * x_hat_t(t-1) - u(t-1);
%     s_bar(t) = sqrt(1/(P(t)+V)) * s(t);
% % Quantizer here without adaptation
%     if s_bar(t) > limit1
%         s_barQ(t) = limit1;
%     elseif s_bar(t) < limit2
%         s_barQ(t) = limit2;
%     else 
%         s_barQ(t) = delta * round( (s_bar(t)/delta) );
%     end                
%     y_tilde(t) = s_barQ(t) - s_bar_tilde(t);
%     x_hat(t) = alpha * x_hat_t(t-1) + u(t-1);
%     x_hat_t(t) = x_hat(t) + K(t)*y_tilde(t);
%     if t<max(T)
%         u(t) = -L(t)*x_hat_t(t);
%         x(t+1) = alpha*x(t)+u(t)+w(t);
%         sumacQ(t+1) = (Q*(x(t)^2)) + (R*(u(t))^2);   
%     end
% end
% %% Optimal cost without adaptive quantizer
% 
% % J_opt_sim = (  F*(x(T+1)^2) + sumac + (Q*(x(1)^2)) + (R*(u(1))^2) ) / (T+1);
% J_opt_simQ = ( sum(sumacQ) +  F*(x(max(T))^2) + (Q*(x(1)^2)) + (R*(u(1))^2) ) / length(T) ;
% cumsumacQ = cumsum(sumacQ) ./ T;


%% Plots
figure 
hold on
% title('Entropy of Output of Quantizer vs Optimal cost bounds: Rate Limited LQG');
xlabel('Time');
ylabel('Cumulative Cost (in dB)');
plot(T(10^0:end),10*log10(cumsumac(10^0:end)));
% plot(T(10^0:10^2),10*log10(cumsumacQ(10^0:10^2)));
% line([0 max(10^2)],[10*log10(b) 10*log10(b)]);
axis tight;
% legend('With Adaptive Quantizer','Without Adaptive Quantizer');
hold off


