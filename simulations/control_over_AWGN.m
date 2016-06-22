T=1000; %total time instants
Q=1;R=1;F=1;V=0;W=1;
alpha=2;
S(T)=1;
for SNR=4:1:50
    for t=T:-1:2
        S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q;
    end
%     S(1)
    P(1)=5;
    for t=1:1:T
        P(t+1)=((alpha^2)*P(t)*(1-(P(t)*SNR)/((P(t)+V)*(1+SNR))))+W;
    end
    temp=0;
    for t=2:1:T
        temp=temp+((Q-S(1))*P(t)+(alpha^2)*P(t-1)*S(t)+W*S(t));
    end
    J_opt(SNR)=Q*P(1)+temp;
end
SNR=1:1:50;
scatter(SNR,J_opt)
%
% for t=T:-1:2
%         display('here')
%         S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q
% end