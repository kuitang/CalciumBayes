function [ Q ] = e_step_smc( M, tau, delta, sigma, beta, lambda, b, w, data )
%e_step_smc Perform the sampling SMC E-step
%   M - number of particles (scalar)
%   tau - decay time constant (scalar)
%   delta - discrete time interval size (scalar)
%   sigma - noise standard deviation (scalar)
%   beta - higher order interaction (N x N x S)
%   lambda - indirect influences (N x S)
%   b - baseline rates (1 x N)
%   w - connectivity matrix (N x N)
%   data - spike trains (N x T sparse)
% data(n,t) = 1 <=> neuron n fired at time t

[x_ S]  = size(lambda);
[N T] = size(data);

sd = sigma*sqrt(delta);

Nthr = M / 10; % wild guess

% Treat each neuron independently
% h(j,t,l)
for i = 1 : N
    h = zeros(N, T, M);
    % At the first timestep, we have a mean-zero normal distribution.
    % Draw M samples and give them uniform importance.
    pf = ones(T, M) / M;
    h(:,1,:) = normrnd(0, sd, N, M);
    % Now draw samples using modified Equation (11) in Mischenko
    for t = 2 : T
        % They say to use stratified resampling after doing this. Why
        % not just directly importance-sample?
        %
        % Read up on stratified resampling and understand why you do it
        % LATER.
        %
        
        % Draw from the one-step-ahead proposal
        for m = 1 : M
            % Draw a new history term from the proposal (transition)
            % distribution
    
            h_mean = (1 - delta)/tau * h(:,t-1,m) + data(:,t-1)';
            h(:,t,m) = normrnd(h_mean, sd);
            
            % Compute J
            I = 0;
            if t - S > 0
                for s = 1 : S
                    I = I + beta(i,:,s) * data(:,t-1) + lambda(i,s);
                end
            end
    
            J = b(i) + I + w(i,:) * h(:,t,m)'
            
            % Compute probabilities
            % q = normpdf(h(:,t,m), h_mean, sd)
            emission_param = 1 - exp(-exp(J) * delta);
            if data(i,t)
                emission_prob = emission_param;
            else
                emission_prob = 1 - emission_param;
            end
            pf(t, m) = emission_prob * pf(t-1, m);            
        end
        
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
    end
        
    % Now run the backwards pass
end

end

