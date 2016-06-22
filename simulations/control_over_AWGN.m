T=1000; %total time instants
Q=1;R=1;F=1;V=0;W=1;
alpha=2; %system parameter
S(T)=F; %assumption from paper 
for SNR=4:1:50  %starting at SNR=4 to ensure stability according to condition given in the paper
    for t=T:-1:2
        S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q; % controller ARE iterations
    end
    P(1)=5; %assuming Po, here denoted as P(1)
    for t=1:1:T
        P(t+1)=((alpha^2)*P(t)*(1-(P(t)*SNR)/((P(t)+V)*(1+SNR))))+W; % KF gain calculation
    end
    temp=0;
    for t=2:1:T
        temp=temp+((Q-S(1))*P(t)+(alpha^2)*P(t-1)*S(t)+W*S(t));
    end
    J_opt(SNR)=Q*P(1)+temp; %optimal performance cost value
end
SNR=1:1:50;
%plotting
scatter(SNR,J_opt)
%
% for t=T:-1:2
%         display('here')
%         S(t-1)=((alpha^2)*R*S(t))/(S(t)+R) + Q
% end