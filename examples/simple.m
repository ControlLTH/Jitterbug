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
P = 1/(s^2-1);     % The process

R1 = 10;            % Input noise
R2 = 0.001;         % Measurement noise
Q = diag([10 0.001]); % Cost J = E(y^2 + 0.001*u^2)

h = 0.05;         % Sampling period
tau = 0;           % Expected input-output latency

S = 1;                           % Sampler system
CA = lqgdesign(P,Q,R1,R2,h);      % LQG controller

% Pd = c2d(P,h);                 
% margin(-C*Pd)                  % Open-loop Bode diagram
% impulse(feedback(-C*Pd,1))     % Closed-loop impulse response

dt = h/10;                       % Time-grain

%Ptau = 1;                         % zero delay
Ptau = [0 0 0 0 0 1.0 0 0 0 0 0];                 % constant delay
%Ptau = ones(1,round(h/dt));      % uniform random delay
%Ptau = Ptau/sum(Ptau);

N = initjitterbug(dt,h);         % Initialize Jitterbug
N = addtimingnode(N,1,Ptau,2);   % Add node 1 (the periodic node)
N = addtimingnode(N,2);          % Add node 2
N = addcontsys(N,1,P,3,Q,R1,R2); % Add sys 1 (P), input from sys 3
N = adddiscsys(N,2,S,1,1);       % Add sys 2 (S), input from 1, exec in 1
N = adddiscsys(N,3,CA,2,2);       % Add sys 3 (CA), input from 2, exec in 2
N = calcdynamics(N);             % Calculate the internal dynamics
J = calccost(N) % Calculate the cost
