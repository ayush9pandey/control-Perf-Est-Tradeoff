A = [0 1 0;0 0 1;1 0 0];    
B = [0.3 1;0 1;-0.3 0.9];
C = [1.9 1.3 1];  
D = [0.53 -0.61];
sys = ss(A,B,C,D);
%Define the noise covariance data and the weighting matrices 
nx = 3;    %Number of states
ny = 1;    %Number of outputs
Qn = [4 2 0; 2 1 0; 0 0 1]; %covariance matrix for w
for i=1:1:1
    lambda=i*1;
    Rn = 0.7*lambda; %covariance matrix for v
    R = [1 0;0 2]; %quadratic cost weight matrix R
    QXU = blkdiag(0.1*eye(nx),R);
    QWV = blkdiag(Qn,Rn);
    % QI = eye(ny);
    %Form the LQG regulator
    KLQG = lqg(sys,QXU,QWV)
    Req = TuningGoal.LQG(wname,zname,QWV,QXU) 
end
