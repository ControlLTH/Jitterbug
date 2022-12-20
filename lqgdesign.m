function [ctrl,L,obs,K,Kf,sysd] = lqgdesign(sys,Qc,R1c,R2,h,tau,nodir,impulse)
% [ctrl,L,obs,K,Kf,sysd] = lqgdesign(sys,Q,R1,R2,h,tau,nodir,impulse)
%
% Design a discrete-time LQG controller for a continuous-time plant with
% constant or random time delay, assuming a continuous-time cost function. Note
% that the state/input noise is continuous-time, while the measurement noise is
% discrete-time.
%
% Arguments:
% sys      The continuous-time plant to be controlled. 
% Qc       The continuous-time cost function is given by
%          [x(t);u(t)]'*Q*[x(t);u(t)] (state-space systems) or
%          [y(t);u(t)]'*Q*[y(t);u(t)] (transfer-function/zpk systems).
% R1c      The continuous-time state/input noise covariance matrix.
% R2       The discrete-time measurement noise covariance matrix.
% h        The sampling period of the controller.
% 
% Optional arguments:
% tau      For a fixed delay, tau is a scalar >= 0. For a random delay, tau
%          is a matrix, where each row specifies [delay probability], and the
%          delay variability may not be larger than h.
% nodir    If non-zero, a controller without direct term is produced,
%          u(k) = -L \hat x(k|k-1)
% impulse  If non-zero, a controller assuming impulse controls is produced
%
% Return values:
% ctrl     The LQG controller as an LTI system, assuming positive feedback.
% L        The state feedback gain vector.
% Obs      The observer as an LTI system.
% K, Kbar  The observer gains.
% sysd     The sampled time-delayed plant as an LTI system. For time-varying
%          delays, this is the average discrete-time evolution of the system.

% Sanity checks
if nargin < 5
  error('To few arguments to function: lqgdesign(sys,Qc,R1c,R2,h[,tau,nodir,impulse])');
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
  error('LTI system delay is not allowed (use tau explicitly).');
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

n = size(A,1);
r = size(B,2);
p = size(C,1);

if isempty(Qc)
  Qc = zeros(n+r,n+r);
else
  switch origclass
   case 'ss'
    if ~isequal(size(Qc),[n+r n+r])
      error(['For state-space systems, the cost Qc should be a (n+r)*(n+r)' ...
	     ' matrix punishing all states and inputs: [x;u]^T*Qc*[x;u].'])
    end
   case {'tf','zpk'}
    if ~isequal(size(Qc),[p+r p+r])
      error(['For transfer function systems, the cost Qc should be a (p+r)*(p+r)' ...
	     ' matrix punishing all outputs and inputs: [y;u]^T*Qc*[y;u].'])
    else
      xutoyu = blkdiag(C,eye(r));
      Qc = xutoyu'*Qc*xutoyu;
    end
  end
end

if isempty(R1c)
  R1c = zeros(n,n);
else
  switch origclass
   case 'ss'
    if ~isequal(size(R1c),[n n])
      error(['For state-space systems, the noise R1c should be an' ...
	     ' n*n matrix, where n is the number of states.'])
    end
   case {'tf','zpk'}
    if ~isequal(size(R1c),[r r])
      error(['For transfer function systems, the noise R1c should' ...
	     ' be an r*r matrix, where r is the number of inputs.'])
    else
      R1c = B*R1c*B';
    end
  end
end

if isempty(R2)
  R2 = zeros(p,p);
else
  if ~isequal(size(R2),[p p])
    error(['The noise R2 should be an p*p matrix, where p is the' ...
	   ' number of outputs.']);
  end
end

if nargin < 6 || isempty(tau)
  tau = 0;
end

if nargin < 7 || isempty(nodir)
  nodir = 0;
end

if nargin < 8 || isempty(impulse)
	impulse = 0;
end

if impulse > 0
	% impulse control: check the cost function and modify B
	if norm(Qc(n+1:end,n+1:end)) > 0 || norm(Qc(1:n,n+1:end)) > 0
		error('Control signal cost must be zero for impulse control')
    end
    % Modify B matrix to account for impulse control
	Bimpulse = B;  % save the B matrix
	B = 0*B;       % set the B matrix to zero when ZOH-sampling the process
end


if isequal(size(tau),[1 1])  %%% Constant delay design %%%
  
  inttau = max(0,ceil(tau/h-1)); % nbr of whole samples extra delay
  fractau = tau - inttau*h;
  
  % Sample the plant, the cost, and the noise
  [Phi,R1] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,h);
  [phi1,~,q1] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,fractau);
  [phi0,~,q0] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,h-fractau);

  Phi = Phi(1:n,1:n);
	R1 = R1(1:n,1:n);
  
	if impulse > 0
		% impulse control
		Gamma1 = 0*B;
		Gamma0 = phi0(1:n,1:n)*Bimpulse;
	else
		% zero-order hold control
		Gamma1 = phi0(1:n,1:n)*phi1(1:n,n+1:n+r);
		Gamma0 = phi0(1:n,n+1:n+r);
	end
    
  % Build extended state-space model
  Phie = [Phi Gamma1; zeros(r,n+r)];
  Gammae = [Gamma0; eye(r)];
  Ce = [C zeros(size(C,1),r)];
  Ge = [eye(n);zeros(r,n)];
  Q1e = q1+phi1(1:n,:)'*q0(1:n,1:n)*phi1(1:n,:);
  if impulse > 0
     % Modify cost matrices for impulse control
     Be = [Bimpulse; zeros(r)];
     Q12e = Q1e*Be;
     Q2e = Be'*Q1e*Be;
  else
    Q2e = q0(n+1:end,n+1:end);
    Q12e = [phi1(1:n,:)'*q0(1:n,n+1:end)];
  end
 
  % Add additional integer delays
  if inttau > 0
    Phie = blkdiag([Phie Gammae],eye((inttau-1)*r));
    Gammae = zeros(size(Phie,1),r);
    Phie = [Phie; zeros(r,size(Phie,2))];
    Gammae = [Gammae; eye(r)];
    Ce = [Ce zeros(size(Ce,1),inttau*r)];
    Ge = [Ge; zeros(inttau*r,size(Ge,2))];
    Q1e = [Q1e Q12e; Q12e' Q2e];
    Q1e = blkdiag(Q1e,zeros((inttau-1)*r));
    Q2e = zeros(size(Q2e));
    Q12e = zeros(size(Q1e,1),size(Q2e,2));
  end
  
  % Solve Riccati equation for state feedback gains 
  [S,~,L] = dare(Phie,Gammae,Q1e,Q2e,Q12e);
%  Lv = (Gammae'*S*Gammae)\(Gammae'*S)
  
  % Solve Riccati equation for Kalman filter gains
%  [P,~,Kt] = dare(Phie',Ce',R1,R2);
%  K = Kt'
%  Kf = P*Ce'/(Ce*P*Ce'+R2);
%  Kv = R12/(Ce*P*Ce'+R2);
    
  sysk = ss(Phie,[Gammae Ge],Ce,0,h);
  
  if nodir == 0
		[obs,K,~,Kf] = kalman(sysk,R1,R2);
  else
		[obs,K,~,Kf] = kalman(sysk,R1,R2,'delayed');
  end
  
  ctrl = lqgreg(obs,L);
  
  
  

  sysd = ss(Phie,Gammae,Ce,0,h);
  
else  %%% Random delay design %%%

  if size(tau,2) ~= 2
    error('For random delays, tau must be an mx2 matrix with [delay probability] pairs')
  end
  if (max(tau(:,1)) - min(tau(:,1))) > h + 10*eps
    error('Cannot handle delay variation larger than the sampling interval')
  end
  if (abs(sum(tau(:,2))-1) > 10*eps)
    warning('Delay probabilities have been corrected to sum up to one')
    tau(:,2) = tau(:,2)/sum(tau(:,2));
	end
  
  mintau = min(tau(:,1));
  inttau = max(0,ceil(mintau/h-1));  % nbr of whole samples extra delay
  remtau = h - mintau + inttau * h;
 
	[phih,rh,~,~] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,h);
	[~,rremtau,~,~] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,remtau);
	[phihremtau,~,~,~] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,h-remtau);

	 Phi = phih(1:n,1:n);
	 R12 = phihremtau(1:n,1:n)*rremtau(1:n,1:n);
	 Re = blkdiag([rh(1:n,1:n) R12; R12' rremtau(1:n,1:n)],0);
	
  for k = 1:length(tau(:,1))
		
    fractau = tau(k,1) - mintau;
    
    % Sample the plant, the cost, and the noise
    [phi1, ~,q1] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,fractau);
    [phi0,~,q0] = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,h-fractau);
	
    if impulse > 0
		% impulse control
		Gamma1 = 0*B;
		Gamma0 = phi0(1:n,1:n)*Bimpulse;
	else
		% zero-order hold control
		Gamma1 = phi0(1:n,1:n)*phi1(1:n,n+1:n+r);
		Gamma0 = phi0(1:n,n+1:n+r);
	end
    
    tauprim = min(fractau,remtau);
    phiprim = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,remtau);
    phi1prim = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,tauprim);
    phi0prim = calcc2d([A B;zeros(r,n+r)],blkdiag(R1c,zeros(r)),Qc,remtau-tauprim);

    Phiprim = phiprim(1:n,1:n);
    Gamma1prim = phi0prim(1:n,1:n)*phi1prim(1:n,n+1:n+r);
    Gamma0prim = phi0prim(1:n,n+1:n+r);
    
    Phie{k} = [Phi zeros(n) Gamma1; Phiprim zeros(n) Gamma1prim; zeros(r,2*n+r)];
    Gammae{k} = [Gamma0; Gamma0prim; eye(r)];
    Ce{k} = [zeros(p,n) C zeros(p,r)];
    Q1e{k} = blkdiag(q1+phi1(1:n,:)'*q0(1:n,1:n)*phi1(1:n,:),zeros(n));
    Q2e{k} = q0(n+1:end,n+1:end);
    Q12e{k} = [phi1(1:n,:)'*q0(1:n,n+1:end); zeros(n,r)];
		
	  % Add additional integer delays
		if inttau > 0
			Phie{k} = blkdiag([Phie{k} Gammae{k}],eye((inttau-1)*r));
			Gammae{k} = zeros(size(Phie{k},1),r);
			Phie{k} = [Phie{k}; zeros(r,size(Phie{k},2))];
			Gammae{k} = [Gammae{k}; eye(r)];
			Ce{k} = [Ce{k} zeros(size(Ce{k},1),inttau*r)];
			Q1e{k} = [Q1e{k} Q12e{k}; Q12e{k}' Q2e{k}];
			Q1e{k} = blkdiag(Q1e{k},zeros((inttau-1)*r));
			Q2e{k} = zeros(size(Q2e{k}));
			Q12e{k} = zeros(size(Q1e{k},1),size(Q2e{k},2));
		end
    
	end
	
	Re = blkdiag(Re,zeros(inttau*r));
		
  % Solve stochastic Riccati equation iteratively
  S = zeros(size(Phie{1},1));
  Snew = eye(size(S,1));
  while norm(S-Snew) > 1e-9
    S = Snew;
    X = zeros(size(S)+r);
    for k=1:size(tau,1)
      X = X + tau(k,2) * [Phie{k} Gammae{k}]'*S*[Phie{k} Gammae{k}] + [Q1e{k} Q12e{k}; Q12e{k}' Q2e{k}];
    end
    L = X(end-r+1:end,end-r+1:end) \ X(end-r+1:end,1:end-r);
    Snew = X(1:end-r,1:end-r) - L'*X(end-r+1:end,end-r+1:end)*L;
	 if norm(Snew) > 1e10
			error('S problem!')
		end
  end
  
  % Compute mean Phi and Gamma matrices, design Kalman filter
  Phibar = zeros(size(Phie{1}));
  Gambar = zeros(size(Gammae{1}));
  for k=1:size(tau,1)
    Phibar = Phibar + tau(k,2) * Phie{k};
    Gambar = Gambar + tau(k,2) * Gammae{k};
  end
	
  sysk = ss(Phibar,[Gambar eye(size(Gambar,1))],Ce{1},0,h);
  
  if nodir
      [obs,K,~,Kf] = kalman(sysk,Re,R2,'delayed');
  else
      [obs,K,~,Kf] = kalman(sysk,Re,R2,'current');
  end
 		
  ctrl = lqgreg(obs,L);
  
  sysd = ss(Phibar,Gambar,Ce{1},0,h);
  
end


