close all
T=10^5 - 1; %total time instants (T-1)
Q=1;R=1;F=1;V=0;W=1; B=1; C = 1; n = 1; ip = 1;
alpha=2; %system parameter
S(T+1)=F; %assumption from paper 
P(1)=5; %assuming Po, here denoted as P(1)
e = exp(1);

v=sqrt(V)*randn((T+1),1);
x(1)=sqrt(P(1))*randn(1,1);
s_bar_tilde=zeros((T+1),1);


 %% Gaussian W
w = sqrt(W) * randn((T+1),1); 
hw = log(sqrt(W) * sqrt(2*pi*e));

c_1 = 1/W;
c_0 = 0;

%% Laplacian Distribution

% w = laprnd((T+1), 1, 0, sqrt(W));
% hw = 0.5 * log(2*W*e^2);
% 
% c_0 = sqrt(2)/sqrt(W);
% c_1 = 0;
% % 

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
%quantization bin size

Deltas = linspace(0.001, 5, 100);

for i_delta = 1 : length(Deltas)
    delta = Deltas(i_delta);    
   

    %% simulate the flow
    y(1)=x(1)+v(1);
    x_hat_t(1)=0;
    u(1)=0;
    x(2)=alpha*x(1)+u(1)+w(1);
    sumac=0;
    for t=2:1:(T+1)
        y(t) = x(t) + v(t);
        s(t) = y(t) - alpha * x_hat_t(t-1) - u(t-1);
        s_bar(t) = sqrt(1/(P(t)+V)) * s(t);
        s_barQ(t) = quant(s_bar(t), delta);
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
    ent_sum = 0;
    lenQ = (max(s_barQ));
    for i = 1 : 200*lenQ
        xd1 = ((2*i)-1)*(delta/2);
        xd2 = ((2*i)+1)*(delta/2);
        phi1 = (1/2) * (1 + erf((xd1-mean(s_barQ))/sqrt(2*cov(s_barQ))));
        phi2 = (1/2) * (1 + erf((xd2-mean(s_barQ))/sqrt(2*cov(s_barQ))));
        pr(i) = phi2 - phi1;
        if pr(i) == 0
            addE = 0;
        else 
            addE = pr(i)*log(1/pr(i));
        end
        
        ent_sum = ent_sum + addE;
        
    end
%     if isnan(ent_sum)
%         display(delta);
%     end
    phi1 = (1/2) * (1 + erf(((-delta/2)-mean(s_barQ))/sqrt(2*cov(s_barQ))));
    phi2 = (1/2) * (1 + erf(((delta/2)-mean(s_barQ))/sqrt(2*cov(s_barQ))));
    pr0 = phi2 - phi1;
    ent0 = pr0*log(1/pr0);
    hq = log(sqrt(cov(s_barQ)) * sqrt(2*pi*e));  %gaussian output of quanitzer
% %     hq = 0.5 * log(2*cov(s_barQ)*e^2); %laplacian
    H_q(i_delta) = 2 * ent_sum + ent0;

    J_opt_sim(i_delta) = (  F*(x(T+1)^2) + sumac + (Q*(x(1)^2)) + (R*(u(1))^2) ) / (T+1);
   
   

%     H_qc(i_delta) = hq - log(delta); % entropy of output of a uniform quantizer (approximate)
    
    %compute lower bound
%  %     Fully observable
%     Rb(i_delta) = H_q(i_delta) + log2(e) + log2(H_q(i_delta) + 1);
    b(i_delta) = bmin + (n * (Nw * (abs(det(M))^(1/n)))) / (exp( (2/n) * (H_q(i_delta) - log(abs(det(alpha))))) - 1);
%     
    % Partially Observable
%     fact = (((abs(det(C * (Td - W) * C' + V))^(1/n)) + (abs(det(C'*C))^(1/n))*Nw)*(abs(det(Kd * M))^(1/n)))/(exp( (2/n) * (H_q(i_delta) - log(abs(det(A)))))-1);
%     b(i_delta) = bmin + n * fact;


    

            
     % compute upper bound
%    
% %     alphaN = (n/2)*log( (2*e)/n ) + log( gamma((n/2) + 1));
% %     rho = ( (sqrt(pi) * ((n + 1)^(1/(2*n))) ) / ( (gamma((n/2) + 1)) ^ (1/n)) ) * sqrt( (n * (n + 2)) / (12 * (n + 1) ) );
%     c_0 = 0;
%     c_1 = 1/W;
%     sigM_max = max(svds(M));
%     sigM_min = min(svds(M));
%     sigMA_max = max(svds(M^(1/2) * alpha));
%     sigMA_min = min(svds(M^(1/2) * alpha));
% 
%     syms bU bminU c_1U sigMA_maxU sigM_minU WU c_0U eU bet RbU nU MU NwU piU alphaU
%     factor1 = (1/2) * log(( (2*piU*eU*(nU+2)) / (12*(nU+1)) ) * (nU + 1) ^(1/nU) );
%     factor0 = (nU/2) * log( (NwU * det(MU)^(1/nU)) / ((bU - bminU)/nU) );
%     bet = 2 * sqrt(bU - bminU) * log( eU*(c_1U*sigMA_maxU*( sqrt((bU - bminU)/sigM_minU) + sqrt(trace(WU)) ) + c_0U + sqrt(bU - bminU) ) ) + sqrt( (bU - bminU)/sigM_minU ) * log( eU * ( (c_1U/2)*sqrt(trace(WU)) + (c_1U/2) * sqrt(trace(WU) + ( (bU - bminU)/sigM_minU) ) + c_0U));
%     eqn = RbU == log(det(alphaU)) + factor1 + factor0 + bet;
%     rewrite(eqn,'log');
%     bU = solve(eqn,bU,'IgnoreAnalyticConstraints',1);
%    
    
%     syms bU
%     solx = solve( (H_qc(i_delta) - log( abs(det(alpha)) ) - n * log(rho * c_1)) == log(( ((n * Nw * (abs(det(M))^(1/n))) / (bU - bmin)) ^ (n/2) ) * ((e*(c_1 * sigMA_max * ( sqrt( (bU - bmin) / sigM_min) + sqrt(trace(W)) ) + c_0 + sqrt(bU - bmin))) ^ (2*sqrt(bU - bmin) )) * ( (e * ( (c_1/2)*sqrt(trace(W)) + (c_1/2)*sqrt(trace(W) + (bU - bmin)/sigM_min) + c_0) )^ sqrt((bU - bmin)/sigM_min) )), bU, 'Real', true);
%     bUpper(i_delta) = double(solx);    
%     
%     syms bU
%     solbU = solve( H_qc(i_delta) == log(det(alpha)) + (n/2) * log( (Nw * det(M)^(1/n)) / ((bU - bmin)/n) ) + (1/2)*log( (2*pi*e*(n+2)*((n+1)^(1/n))) / (12* (n+1)) ) + 2 * sqrt(bU - bmin) * log( e*(c_1*sigMA_max*( sqrt((bU - bmin)/sigM_min) + sqrt(trace(W)) ) + c_0 + sqrt(bU - bmin) ) ) + sqrt( (bU - bmin)/sigM_min ) * log( e * ( (c_1/2)*sqrt(trace(W)) + (c_1/2) * sqrt(trace(W) + ( (bU - bmin)/sigM_min) ) + c_0)), bU); 
%     bU = double(solbU);
end
sigM_max = max(svds(M));
sigM_min = min(svds(M));
sigMA_max = max(svds(M^(1/2) * alpha));
sigMA_min = min(svds(M^(1/2) * alpha));
maxb = max(J_opt_sim);
Cost = linspace((bmin + 0.01), (bmin + 20) , 100);
for i_Cost = length(Cost) : -1 : 1
    bU = Cost(i_Cost);
    %% Lower Bound
    
    %Fully Observable
%     Rb(i_Cost) = log(det(alpha)) + (n/2) * log(1 + (Nw * abs(det(M))^(1/n)) / ((bU - bmin)/n));
    
    %Partially Observable
    Rb(i_Cost) = log(det(alpha)) + (n/2) * log( 1 + ( n * ((abs(det( C * (Td - W) * C' + V))^(1/n)) + (abs(det(C'*C))^(1/n)) * Nw) * (abs(det(Kd*M))^(1/n))) / (bU - bmin) );
    
    %% Upper Bound
    
    factor1 = (1/2) * log(( (2*pi*e*(n+2)) / (12*(n+1)) ) * (n + 1) ^(1/n) );
%     factor1 = (1/2) * log( ( (2*pi*e*(n+2)) / (12*(n+1)) ) ) * (n + 1) ^(1/n) ;
    factor0 = (n/2) * log( (Nw * det(M)^(1/n)) / ((bU - bmin)/n) );
    bet = 2 * sqrt(bU - bmin) * log( e*(c_1*sigMA_max*( sqrt((bU - bmin)/sigM_min) + sqrt(trace(W)) ) + c_0 + c_1 * sqrt(bU - bmin) ) ) + sqrt( (bU - bmin)/sigM_min ) * log( e * ( (c_1/2)*sqrt(trace(W)) + (c_1/2) * sqrt(trace(W) + ( (bU - bmin)/sigM_min) ) + c_0));
%     bet = 2 * sqrt(bU - bmin) * log( (c_1*sigMA_max*( sqrt((bU - bmin)/sigM_min) + sqrt(trace(W)) ) + c_0 + c_1 * sqrt(bU - bmin) ) ) + sqrt( (bU - bmin)/sigM_min ) * log( ( (c_1/2)*sqrt(trace(W)) + (c_1/2) * sqrt(trace(W) + ( (bU - bmin)/sigM_min) ) + c_0));
    RbU(i_Cost) = log(det(alpha)) + factor1 + factor0 + bet;
%     RbU(i_Cost) = log(det(alpha)) + factor0 + factor1;

end

% 
 
% J_opt_sim = inpaint_nans(J_opt_sim);
% H_q = inpaint_nans(H_q,3);

% J_opt_simP = J_opt_sim(1:17);
% H_qP = H_q(1:17);
%% Plots
figure 
hold on
% title('Entropy of Output of Quantizer vs Optimal cost bounds: Rate Limited LQG');
xlabel('Cost');
ylabel('Quantizer Entropy, nats');
% scatter(10 * log10(SNRs), J_opt,'+');
% plot(10 * log10(SNRs), J_opt_sim);
plot(Cost, Rb,'r');
plot(J_opt_sim, H_q,'--');
% plot(J_opt_simP, H_qP,'--');
plot(Cost, RbU, 'm-.');
%Horizontal Line
line([0 max(Cost)], [log(abs(det(alpha))) log(abs(det(alpha)))]);
% line([0 max(b)], [1 1]);
%vertical Line
line([bmin bmin], [0 max(RbU)]);
% line([bmin bmin], [0 max(max(Rb))]);
legend('Computed Lower Bound','Simulated System', 'Computer Upper Bound');
% legend('Computed Lower Bound', 'Computer Upper Bound');
% legend('Lower (Converse) Bound','Simulated System');
hold off