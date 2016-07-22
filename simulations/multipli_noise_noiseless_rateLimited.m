close all
T=10^5 - 1; %total time instants (T-1)
Q=1; R=0; F=1; V=0; W=0.5;
n = 1;
meanA = 2;
varA = 0.3^2;
e = exp(1);
A = meanA + sqrt(varA) * randn( (T + 1) , 1 ); % random system parameter modeling the multiplicative nosie
ha = log2(sqrt(varA) * sqrt(2*pi*e)); % For Gaussian A in bits
Na = (1/(2*pi*e)) * (2 ^ (2*ha/n));


S(T+1) = F; %assumption from paper 
P(1) = 1; %assuming Po, here denoted as P(1)

v = sqrt(V)*randn((T+1),1);
x(1) = sqrt(P(1))*randn(1,1);
s_bar_tilde = zeros((T+1),1);

% %Gaussian W
w = sqrt(W) * randn((T+1),1); 
hw = log2(sqrt(W) * sqrt(2*pi*e));

% Laplacian Distribution
% 
% w = laprnd((T+1), 1, 0, sqrt(W));
% hw = 0.5 * log2(2*W*e^2);

Nw = (1/(2*pi*e)) * 2^(2*hw/n);

%% Change Quantization Bin Size Here
Deltas = linspace(1e-2, 2, 100);

for i_delta = 1 : length(Deltas)
    delta = Deltas(i_delta);    
    for t=(T+1):-1:2
        S(t-1)=((A(t)^2)*R*S(t))/(S(t)+R) + Q; % controller ARE iterations
    end

    for t=1:1:T
        P(t+1)=(A(t)^2)*P(t)*(1-(P(t)/(P(t)+V)))+W;  % KF ARE
    end

    % %%%%%%%%%simulating the system %%%%%%%%


    for t=(T+1):-1:2
        %controller gain
        L(t-1)=(A(t)*S(t))/(S(t)+R);
    end
    for t=1:1:(T+1)
        %kalman gain
        K(t)=(P(t)/sqrt(P(t)+V));
    end
    %simulate the flow
    y(1)=x(1)+v(1);
    x_hat_t(1)=0;
    u(1)=0;
    x(2)=A(1)*x(1)+u(1)+w(1);
    sumac=0;
    for t=2:1:(T+1)
        y(t) = x(t) + v(t);
        s(t) = y(t) - A(t-1) * x_hat_t(t-1) - u(t-1);
        s_bar(t) = sqrt(1/(P(t)+V)) * s(t);
        s_barQ(t) = quant(s_bar(t), delta);
        y_tilde(t) = s_barQ(t) - s_bar_tilde(t);
        x_hat(t) = A(t-1) * x_hat_t(t-1) + u(t-1);
        x_hat_t(t) = x_hat(t) + K(t)*y_tilde(t);
%         x_hat_tQ(t) = quant(x_hat_t(t), delta);
        if t<(T+1)
            u(t) = -L(t)*x_hat_t(t);
            x(t+1) = A(t)*x(t)+u(t)+w(t);
            sumac = sumac + (x(t)^2);
        end         
    end
    ent_sum = 0;
    lenQ = (max(s_barQ));
    %% Improve probability estimation here
    for i = 1 : 100*lenQ
        xd1 = ((2*i)-1)*(delta/2);
        xd2 = ((2*i)+1)*(delta/2);
        phi1 = (1/2) * (1 + erf((xd1-mean(s_barQ))/sqrt(2*cov(s_barQ))));
        phi2 = (1/2) * (1 + erf((xd2-mean(s_barQ))/sqrt(2*cov(s_barQ))));
        pr(i) = phi2 - phi1;
        if pr(i) == 0
            addE = 0;
        else 
            addE = pr(i)*log2(1/pr(i));
        end
        
        ent_sum = ent_sum + addE;
        
    end
%     if isnan(ent_sum)
%         display(delta);
%     end
    phi1 = (1/2) * (1 + erf(((-delta/2)-mean(s_barQ))/sqrt(2*cov(s_barQ))));
    phi2 = (1/2) * (1 + erf(((delta/2)-mean(s_barQ))/sqrt(2*cov(s_barQ))));
    pr0 = phi2 - phi1;
    ent0 = pr0*log2(1/pr0);
    hq = log2(sqrt(cov(s_barQ)) * sqrt(2*pi*e));  %gaussian output of quanitzer
% %     hq = 0.5 * log2(2*cov(s_barQ)*e^2); %laplacian
    H_q(i_delta) = 2 * ent_sum + ent0;

        % Mean Square Cost only
    J_opt_sim(i_delta) = ( (x(T+1) ^ 2) + sumac + (x(1)^2) ) / (T+1);
    
end

%% Lower Bound Calculation 
Rate = linspace(1, 20, 100);
for i_rate = 1 : length(Rate)
    Rb = Rate(i_rate);
%     abc = (Nw*meanA*meanA*2^(-2*Rb) + W)/(1 - varA);
%     abcd = (Na*meanA*meanA*2^(-2*Rb))/(1 - varA);
%     abcde = 1 + abcd;
%     J_opt_lower(i_rate) = abc * abcd;
    J_opt_lower(i_rate) =( ( (Nw * (meanA ^ 2) * 2^(-2*Rb) ) + W ) / (1 - varA) ) * ( 1 + ( (Na * (meanA ^2) * 2^(-2*Rb)) / (1 - varA) ) );
end
% Cost = linspace(2, 20, 100);
% for i_cost = 1 : length(Cost)
%     J_lower = Cost(i_cost);
%     aa = Na*Nw*(meanA^4);
%     bb = (meanA^2) * ( Nw + W * Na - varA * Nw);
%     cc = J_lower*(1 - varA) - W + W * varA;
%     Rb(i_cost) = (1/2) * log2( (2*aa) / (-bb + sqrt( (bb^2) - 4*aa*cc ) ) );
% end

   
%     J_opt_lower(i_rate) =( ( (Nw * (meanA ^ 2) * 2^(-2*Rb) ) + W ) / (1 - varA) ) * ( 1 + ( (Na * (meanA ^2) * 2^(-2*Rb)) / (1 - varA) ) );
% Cost = linspace(min(J_opt_lower), max(J_opt_sim), 10);
% for i_cost= 1 : length(Cost)
%     Ratex(i_cost) = 0;
% end
% 
% Rate = cat(2,Rate,Ratex);
% J_opt_lower = cat(2,J_opt_lower,Cost);
%% Plots

figure 
hold on
title('Entropy of Output of Quantizer vs Optimal cost lower bound : Multiplicative Noise w/ Limited Rate');
xlabel('Cost');
ylabel('Quantizer Entropy, bits');
plot(J_opt_lower, Rate,'r');
plot(J_opt_sim, H_q,'--');
% scatter(J_opt_sim, H_q,'--');
% legend('Computed Lower Bound','Simulated System', 'Computer Upper Bound');
legend('Lower Bound','Simulated System');
hold off