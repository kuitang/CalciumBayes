function [ pb h ] = e_step_smc( i, M, tau, delta, sigma, params, data )
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
%   pb - backward sample weights (T x M)
%   h - samples (N x T x M)

beta = params.beta;
lambda = params.lambda;
w = params.w;
b = params.b;

[~, S]  = size(lambda);
[N, T] = size(data);

sd = sigma*sqrt(delta);

Nthr = M / 2; % From [Vogelstein09]

% Treat each neuron independently
% h(j,t,l)
h = zeros(N, T, M);
% At the first timestep, we have a mean-zero normal distribution.
% Draw M samples and give them uniform importance.
pf = ones(T, M) / M;
h(:,S+1,:) = normrnd(0, sd, N, M);

% Now draw samples using modified Equation (11) in [Mischenko11]
lastspike = -S;
% should start at S + 1
for t = S+1 : T %should we just start this at S+1??
    % Sequential importance resampling    
    for m = 1 : M
        % Draw from the one-step-ahead proposal
        % Draw a new history term from the proposal (transition)
        % distribution        
        h_mean = (1 - delta/tau) * h(:,t-1,m) + data(:,t-1);        
        h(:,t,m) = normrnd(h_mean, sd);        

        % Compute J
        I = 0;

        for s = 1 : S-1
            I = I + beta(i,:,s) * data(:,t-(s+1)) + lambda(i,s);
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
% We don't need to save the r's; we just keep one timestep
r = zeros(M, M);

% At t = T, the filtering and smoothing distributions are identical.
pb = zeros(T,M);
pb(T,:) = pf(T,:);

for t_index = 0:(T-S-2)
    
    t = T - t_index;
    
    % Equations (12) and (13) in [Mischenko11]
    for m = 1 : M
        denom = 0;
        for mm = 1 : M
            sigma_mat = sd^2*eye(N);
            %%%%%%%%%
            %%% THIS CALCULATION RETURNS VERY LOW PROB - 
            %%% LOOKING INTO IT - BS
            distance = h(:,t,m) - h(:,t-1,mm);       
            prob = det(6.28318530717959*sigma_mat)^(-.5) * exp(-.5 * (distance' * inv(sigma_mat) * distance));
            disp(prob);
            disp(mvnpdf(h(:,t,m), h(:,t-1,mm), sd^2*eye(N)) * pf(t-1,mm));
            denom = denom + mvnpdf(h(:,t,m), h(:,t-1,mm), sd^2*eye(N)) * pf(t-1,mm);
            
        end        
        for mm = 1 : M
            numer = mvnpdf(h(:,t,m), h(:,t-1,mm), sd^2*eye(N)) * pf(t-1, m);
            r(m, mm) = pb(t, m) * numer / denom;            
        end        
    end
    
    pb(t-1,:) = sum(r(:,:), 2);    
    
    pb(t-1,:) = pb(t-1,:) / sum(pb(t-1,:)) ;
end

end
