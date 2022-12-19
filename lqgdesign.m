function [ctrl,L,obs,K,Kf,sysd,S] = lqgdesign(sys,Q,R1,R2,h,tau,nodir)
% [ctrl,L,obs,K,Kf,sysd] = lqgdesign(sys,Q,R1,R2,h,tau)
%
% Design a discrete-time LQG controller for a continuous-time plant
% with constant time delay, assuming a continuous-time cost function. 
% Note that the state/input noise is continuous-time, while the
% measurement noise is discrete-time.
%
% Arguments:
% sys      The continuous-time plant to be controlled. 
% Q        The continuous-time cost function is given by
%          [x(t);u(t)]'*Q*[x(t);u(t)] (state-space systems) or
%          [y(t);u(t)]'*Q*[y(t);u(t)] (transfer-function/zpk systems).
% R1       The continuous-time state/input noise covariance matrix.
% R2       The discrete-time measurement noise covariance matrix.
% h        The sampling period of the controller.
% 
% Optional arguments:
% tau      The assumed plant (or loop) time delay. May be larger than h.
% nodir    If this argument is non-zero, a controller without direct
%          term is produced
%
% Return values:
% ctrl     The LQG controller as an LTI system, assuming positive feedback.
% L        The state feedback gain vector.
% Obs      The observer as an LTI system.
% K, Kbar  The observer gains.
% sysd     The sampled time-delayed plant as an LTI system.

% Sanity checks
if nargin < 5
  error('To few arguments to function: lqgdesign(sys,Q,R1,R2,h[,tau])');
end
origclass = class(sys);
switch origclass
 case 'ss'
 case 'tf'
 case 'zpk'
  otherwise
  error(['System class ' origclass ' not supported. SYS must be' ...
		    ' either ss, tf, or zpk.']);
end

if ~isct(sys)
  error('System is not continuous time.');
end

if ~isproper(sys)
  error('System is not proper');
end

if (hasdelay(sys))
  disp('Waring: LTI system delay is ignored (use tau explicitly).');
end

sys = ss(sys);
A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;

if (max(max(abs(D))) > eps) 
  error('The continuous system has a direct term, which is not supported.');
else
  D = zeros(size(D));
end
  
totsize = size(A,2)+size(B,2);
outinsize = size(C,1)+size(B,2);

if isempty(Q)
  Q = zeros(totsize);
else
  switch origclass
   case 'ss'
    if size(Q,1) ~= totsize | size(Q,2) ~= totsize
      error(['For state-space systems, the cost Q should be a matrix' ...
	     ' punishing all states and inputs: [x;u]^T*Q*[x;u].'])
    end
   case {'tf','zpk'}
    if size(Q,1) ~= outinsize | size(Q,2) ~= outinsize
      error(['For transfer function systems, the cost Q should be a matrix' ...
	     ' punishing all outputs and inputs: [y;u]^T*Q*[y;u].'])
    else
      xutoyu = blkdiag(C,eye(size(B,2)));
      Q = xutoyu'*Q*xutoyu;
    end
  end
end

if isempty(R1)
  R1 = zeros(size(A,1));
else
  switch origclass
   case 'ss'
    if size(R1,1) ~= size(A,1) | size(R1,2) ~= size(A,2)
      error(['For state-space systems, the noise R1 should be an' ...
	     ' n*n matrix, where n is the number of states.'])
    end
   case {'tf','zpk'}
    if size(R1,1) ~= size(B,2) | size(R1,2) ~= size(B,2)
      error(['For transfer function systems, the noise R1 should' ...
	     ' be an r*r matrix, where r is the number of inputs.'])
    else
      R1 = B*R1*B';
    end
  end
end

if isempty(R2)
  R2 = zeros(size(C,1));
else
  if (size(R2,1) ~= size(C,1) | size(R2,2) ~= size(C,1))
    error(['The noise R2 should be an p*p matrix, where p is the' ...
	   ' number of outputs y.']);
  end
end

if nargin < 6 | isempty(tau)
  tau = 0;
end

if nargin < 7 | isempty(nodir)
  nodir = 0;
end

inttau = max(0,ceil(tau/h-1)); % nbr of whole samples extra delay
fractau = tau - inttau*h;

a1 = size(A,1);  % number of states
b2 = size(B,2);  % number of inputs
sysc = ss(A,B,C,D);
sysd = c2d(sysc,h);
[Phi,Gamma] = ssdata(sysd);

% Sample the plant, the cost, and the noise
[Phi2,r,Q2] = calcc2d([A B;zeros(b2,a1+b2)],blkdiag(R1,zeros(b2)),Q,fractau);
Gamma2 = Phi2(1:a1,a1+1:a1+b2);
[X,r,Q1] = calcc2d([A B;zeros(b2,a1+b2)],blkdiag(R1,zeros(b2)),Q,h-fractau);
Gamma2 = X(1:a1,1:a1)*Gamma2;
Gamma1 = X(1:a1,a1+1:a1+b2);
[Phi,R1,Q,Qconst] = calcc2d([A B;zeros(b2,a1+b2)],blkdiag(R1,zeros(b2)),Q,h);
R1 = R1(1:a1,1:a1);
Phi = Phi(1:a1,1:a1);

% Build extended state-space model
Phie = [Phi Gamma2; zeros(b2,a1+b2)];
Gammae = [Gamma1; eye(b2)];
Ce = [C zeros(size(C,1),b2)];
Ge = [eye(a1);zeros(b2,a1)];
Q1e = Q2+Phi2(1:a1,:)'*Q1(1:a1,1:a1)*Phi2(1:a1,:);
Q2e = Q1(a1+1:end,a1+1:end);
Q12e = [Phi2(1:a1,:)'*Q1(1:a1,a1+1:end)];

% Add additional integer delays
if inttau > 0
  Phie = blkdiag([Phie Gammae],eye((inttau-1)*b2));
  Gammae = zeros(size(Phie,1),b2);
  Phie = [Phie; zeros(b2,size(Phie,2))];
  Gammae = [Gammae; eye(b2)];
  Ce = [Ce zeros(size(Ce,1),inttau*b2)];
  Ge = [Ge; zeros(inttau*b2,size(Ge,2))];
  Q1e = [Q1e Q12e; Q12e' Q2e];
  Q1e = blkdiag(Q1e,zeros((inttau-1)*b2));
  Q2e = zeros(size(Q2e));
  Q12e = zeros(size(Q1e,1),size(Q2e,2));
end

[s,e,L] = dare(Phie,Gammae,Q1e,Q2e,Q12e);
sysk = ss(Phie,[Gammae Ge],Ce,0,h);
[kest,K,p,Kf] = kalman(sysk,R1,R2);
if nodir == 0
  ctrl = lqgreg(kest,L,'current');
else
  ctrl = lqgreg(kest,L);
end
obs = kest;
sysd = ss(Phie,Gammae,Ce,0,h);
S = s;
