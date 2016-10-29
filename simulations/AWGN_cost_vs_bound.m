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

