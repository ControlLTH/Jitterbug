% Jitterbug example: simple.m
% ================================
% The very simplest Jitterbug example. An inverted pendulum is
% controlled by an LQG controller. The sampling is periodic, but
% there is a random input-output latency.
%
% Two timing nodes are used. In the first node, the sampler is
% executed. In the second node, the controller is executed.
%
% Play around with different sampling intervals (h), delay
% distributions (Ptau), delay compensation (tau), etc., and see how
% the cost changes.

s = tf('s');
P = 1/((s+1)*(s+10));   % The process (DC servo)

C = -50*(1+1/(1*s));

R1 = 1;            % Input noise
R2 = 0.0001;         % Measurement noise

Q = diag([1 0.0001]); % Cost J = E(y^2 + 0.001*u^2)

h = 0.1;         % Sampling period

C = c2d(C,h,'foh');

pvec = 0:0.01:1;
Jvec = [];

for p = pvec

tau = 0;           % Expected input-output latency

S = 1;                           % Sampler system


% Pd = c2d(P,h);                 
% margin(-C*Pd)                  % Open-loop Bode diagram
% impulse(feedback(-C*Pd,1))     % Closed-loop impulse response

dt = h;                       % Time-grain

%Ptau = 1;                         % zero delay
Ptau = [0 1-p p];
%Ptau = ones(1,round(h/dt));      % uniform random delay
%Ptau = Ptau/sum(Ptau);

N = initjitterbug(dt,h);         % Initialize Jitterbug
N = addtimingnode(N,1,Ptau,2);   % Add node 1 (the periodic node)
N = addtimingnode(N,2);          % Add node 2
N = addcontsys(N,1,P,3,Q,R1,R2); % Add sys 1 (P), input from sys 3
N = adddiscsys(N,2,S,1,1);       % Add sys 2 (S), input from 1, exec in 1
N = adddiscsys(N,3,C,2,2);       % Add sys 3 (C), input from 2, exec in 2
N = calcdynamics(N);             % Calculate the internal dynamics
J = calccost(N) % Calculate the cost

Jvec = [Jvec J];

end


plot(pvec,Jvec)



