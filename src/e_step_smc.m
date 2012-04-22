function [ r h ] = e_step_smc( i, M, tau, delta, sigma, beta, lambda, b, w, data )
%e_step_smc Perform the sampling SMC E-step for one neuron
%   i - which neuron (scalar)
%   M - number of particles (scalar)
%   tau - decay time constant (scalar)
%   delta - discrete time interval size (scalar)
%   sigma - noise standard deviation (scalar)
%   beta - higher order interaction (N x N x S)
%   lambda - indirect influences (N x S)
%   b - baseline rates (1 x N)
%   w - connectivity matrix (N x N)
%   data - spike trains (N x T sparse)
% 
%   r - marginal sample weights (T x M x M)
%   h - samples (N x T x M)

[x_ S]  = size(lambda);
[N T] = size(data);

sd = sigma*sqrt(delta);

Nthr = M / 10; % wild guess

% Treat each neuron independently
% h(j,t,l)
h = zeros(N, T, M);
% At the first timestep, we have a mean-zero normal distribution.
% Draw M samples and give them uniform importance.
pf = ones(T, M) / M;
h(:,1,:) = normrnd(0, sd, N, M);

% Now draw samples using modified Equation (11) in [Mischenko11]
lastspike = -S;
for t = 2 : T
    % Sequential importance resampling    
    for m = 1 : M
        % Draw from the one-step-ahead proposal
        % Draw a new history term from the proposal (transition)
        % distribution        
        h_mean = (1 - delta)/tau * h(:,t-1,m) + data(:,t-1);        
        h(:,t,m) = normrnd(h_mean, sd);        

        % Compute J
        I = 0;
        if t - S > 0
            for s = 1 : S
                I = I + beta(i,:,s) * data(:,t-1) + lambda(i,s);
            end
        end        
        J = b(i) + I + w(i,:) * h(:,t,m);

        % Compute probabilities
        % q = normpdf(h(:,t,m), h_mean, sd)
        emission_param = 1 - exp(-exp(J) * delta);
        if data(i,t)
            emission_prob = emission_param;
        else
            emission_prob = 1 - emission_param;
        end        
        pf(t, m) = emission_prob * pf(t-1, m);
        %[norm(h(:,t,m)) emission_prob pf(t,m)]
    end
    pf(t,:) = pf(t,:) / sum(pf(t,:));

    % Stratified resampling explained in http://en.wikipedia.org/wiki/Particle_filter
    Neff = 1/(pf(t,:) * pf(t,:)');
    if Neff < Nthr
        % Resample
        % Draw P particles from the current particle set
        idxs = randsample(M, M, true, pf(t, :));
        % Replace the current particle set with a new one
        h(:,t,:) = h(:,t,idxs);
        % Reset weights to uniform
        pf(t, :) = ones(1, M) / M;
    end
    
    if data(i,t)
        lastspike = t;
    end
    if t <= lastspike + S
        E = pf(t, :) * squeeze(h(:,t,:))' / M;
        disp(['t = ' num2str(t) ' spike = ' num2str(data(i,t)) ' sample expectation norm ' num2str(norm(E))]);
    end

end

% Marginal smoothing
r = zeros(T, M, M);
% At t = T, the filtering and smoothing distributions are identical.
pb = pf(T,:);
pb_next = zeros(1, M);

for t = T-1:-1:2
    % Equations (12) and (13) in [Mischenko11]
    for m = 1 : M
        denom = 0;
        for mm = 1 : M
            denom = denom + mvnpdf(h(:,t,m), h(:,t-1,mm), sd*eye(N)) * pf(t-1,mm);
        end        
        for mm = 1 : M
            numer = mvnpdf(h(:,t,m), h(:,t-1,mm), sd*eye(N)) * pf(t-1, m);
            r(t,m,mm) = pb(m) * numer / denom;
        end
        pb_next(m) = sum(r(t,:,m));         
    end
    pb = pb_next;
    pb = pb / sum(pb);    
    r(t,:,:) = r(t,:,:) / sum(sum(r(t,:,:)));    
end

end

