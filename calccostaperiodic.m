function [J,N,P,prob] = calccostaperiodic(N,options)
% J = calccostaperiodic(N)
% [J,N,P,prob] = calccostaperiodic(N,options)
%
% Calculate steady-state variance and cost for the state when system jumps
% between timing nodes defined by N.nodes. ITERATIVE method for APERIODIC
% systems only. (Note that there is no notion of "first node" or "delay since
% first node" in an aperiodic model.)
%
% Returns
% J -- the cost (Inf if unstable)
% N -- system including Markov states in N.states (for debugging)
% P -- covariance matrix between iterations (for debugging)
% prob - state probability distribution between iterations (for debugging)
%
% options: a struct containing any of the following fields:
%   maxiter - Maximum number of iterations (default = 1000)
%   tol - Relative tolerance (default = 1e-6)
%   printout - Print every 100 iterations if non-zero (default = 0)


%%% 1. Translate model into Markov chain with extra delay nodes

% Calculate total number of Markov states
nstates = 0;
nnodes = length(N.nodes);
for k = 1:nnodes
	if min(size(N.nodes{k}.nextprob)) >= 2
		error('Aperiodic solver cannot handle time-dependent random delays')
	end
	N.nodes{k}.startstate = nstates + 1;
	nstates = nstates + length(N.nodes{k}.Ptau);
end

% Translate timing nodes and random delays into Markov states
for k = 1:length(N.nodes)
	startstate = N.nodes{k}.startstate;
	
	% Create timing node state
	N.states{startstate} = [];
	N.states{startstate}.type = 1;         % Type 1 = execution state
	N.states{startstate}.E = N.nodes{k}.E;
	N.states{startstate}.R2 = N.nodes{k}.R2;
	exitstates = [];
	for l = 1:length(N.nodes{k}.next)
		exitstates = [exitstates N.nodes{N.nodes{k}.next(l)}.startstate];
	end
	if length(N.nodes{k}.Ptau) > 1
		nextstates = [startstate+1 exitstates];
		stayprob = sum(N.nodes{k}.Ptau(2:end));
		nextprobs = [stayprob (1-stayprob)*N.nodes{k}.nextprob'];
		exitprobs = N.nodes{k}.Ptau./circshift(1-cumsum(N.nodes{k}.Ptau),1); % NOTE: only valid for l>=2
	else
		nextstates = exitstates;
		nextprobs = N.nodes{k}.nextprob';		
	end
	keep = nextprobs > 0;
	nextstates = nextstates(keep);
	nextprobs = nextprobs(keep);
	N.states{startstate}.nextstates = nextstates;
	N.states{startstate}.nextprobs = nextprobs;

	% Create extra delay states between timing nodes
	for l = 2:length(N.nodes{k}.Ptau)
		N.states{startstate+l-1} = [];
		N.states{startstate+l-1}.type = 2;  % Type 2 = delay state
		N.states{startstate+l-1}.A = N.nodes{k}.A;
		N.states{startstate+l-1}.R1 = N.nodes{k}.R1;
		N.states{startstate+l-1}.Q = N.nodes{k}.Q;
		N.states{startstate+l-1}.Qconst = N.nodes{k}.Qconst;
		
		if l < length(N.nodes{k}.Ptau)
			nextstates = [startstate+l exitstates];
			nextprobs = [1-exitprobs(l) exitprobs(l)*N.nodes{k}.nextprob'];
		else
			nextstates = exitstates;
			nextprobs = N.nodes{k}.nextprob';
		end
		keep = nextprobs > 0;
		nextstates = nextstates(keep);
		nextprobs = nextprobs(keep);
		N.states{startstate+l-1}.nextstates = nextstates;
		N.states{startstate+l-1}.nextprobs = nextprobs;
	end
end

%%% 2. Iteratively calculate Markov state probabilities and total cost

% Default options
maxiter = 1000; 
tol = 1e-6;
printout = 0;

if nargin >= 2
	if isfield(options,'maxiter')
		maxiter = options.maxiter;
	end
	if isfield(options,'printout')
		printout = options.printout;
	end
	if isfield(options,'tol')
		tol = options.tol;
	end
end

ssdim = size(N.states{1}.E,1);              % Dimension of total state-space
P = zeros(ssdim,ssdim,nstates,maxiter+1);   % Covariance in each time step (saved for debugging)
prob = zeros(nstates,maxiter+1);            % Markov state probability distribution in each time step
prob(1,1) = 1;                              % Initial probability distribution

Jlast = 0; J = eps; time = 1;

while (time <= maxiter) && (abs(J-Jlast)/J > tol)
	
	Jlast = J;
	
	% (a) Update step. Execute timing nodes until probability of being in any timing node is zero
	P_update = P(:,:,:,time);
	prob_update = prob(:,time);

	done = 0;
	while done == 0
		done = 1;		
		P_new = zeros(ssdim,ssdim,nstates); % Accumulator for intermediate calculations
		prob_new = zeros(nstates,1);        % Accumulator for intermediate calculations
		
		for k = 1:nstates
			switch N.states{k}.type
				
				case 1 % Timing node: propagate prob and P forward
					
					if prob_update(k) > 0
						nextstates = N.states{k}.nextstates;
						prob_new(nextstates) = prob_new(nextstates) + prob_update(k)*N.states{k}.nextprobs';
						for i = 1:length(nextstates)
							P_new(:,:,nextstates(i)) = P_new(:,:,nextstates(i)) + ...
								N.states{k}.nextprobs(i) * (N.states{k}.E*P_update(:,:,k)*N.states{k}.E' + prob_update(k)*N.states{k}.R2);
						end
						prob_update(k) = 0;
						P_update(:,:,k) = 0*P_update(:,:,k);
						done = 0;
					end
					
				case 2 % Delay state: keep current value of prob and P
					
					prob_new(k) = prob_new(k) + prob_update(k);
					P_new(:,:,k) = P_new(:,:,k) + P_update(:,:,k);
					
			end
		
		end
		
		P_update = P_new;
		prob_update = prob_new;
		
	end
	
	% All probability mass is now in the delay states. Calculate total cost
	
	J = N.nodes{1}.Qconst;  % Qconst is the same for all nodes/states
	for k = 1:nstates
		if prob_update(k) > 0
			J = J + trace(N.states{k}.Q*P_update(:,:,k));
		end
	end
	J = J / N.dt;
	
	% (b) Delay step. Loop through all delay states and propagate probabilities
	% and covariances forward
	
	prob_next = zeros(nstates,1);
	P_next = zeros(ssdim,ssdim,nstates);
	for k = 1:nstates
		if prob_update(k) > 0
			nextstates = N.states{k}.nextstates;
			prob_next(nextstates) = prob_next(nextstates) + prob_update(k)*N.states{k}.nextprobs';
			for i = 1:length(nextstates)
				P_next(:,:,nextstates(i)) = P_next(:,:,nextstates(i)) + ...
				N.states{k}.nextprobs(i) * (N.states{k}.A*P_update(:,:,k)*N.states{k}.A' + prob_update(k)*N.states{k}.R1);
			end
		end
	end
	
	if (mod(time,100) == 0) && printout ~= 0
		fprintf('Time = %d, J = %d\n', time, J);
	end
	
	time = time + 1;
	prob(:,time) = prob_next;
	P(:,:,:,time) = P_next;
	
end

if time > maxiter
	warning('Maximum number of iterations reached')
end
