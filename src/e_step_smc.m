function [ pb h ] = e_step_smc( i, M, tau, delta, sigma, beta, b, w, data )
%e_step_smc Perform the sampling SMC E-step for one neuron
%   i - which neuron (scalar)
%   M - number of particles (scalar)
%   tau - decay time constant (scalar)
%   delta - discrete time interval size (scalar)
%   sigma - noise standard deviation (scalar)
%   beta - higher order interaction for this neuron (N x S)
%   lambda - indirect influences (N x S)
%   b - baseline rate for this neuron (scalar)
%   w - connectivity vector (1 x N)
%   data - spike trains (N x T sparse)
% 
%   pb - backward sample weights (T x M)
%   h - samples (N x T x M)

[N, T] = size(data);
S = size(beta,2) + 1;

sd = sigma*sqrt(delta);

Nthr = M / 2; % From [Vogelstein09]

% Treat each neuron independently
% h(j,t,l)
h = zeros(N, T, M);
% At the first timestep, we have a mean-zero normal distribution.
% Draw M samples and give them uniform importance.
pf = ones(T, M) / M;
pf_new = pf;
h(:,S+1,:) = normrnd(0, sd, N, M);

% Now draw samples using modified Equation (11) in [Mischenko11]
%lastspike = -S;
% should start at S + 2
for t = S+2 : T 
    % Sequential importance resampling
    
    % I doesn't depend on samples, so compute here
    I_terms = beta .* data(:,(t-2):-1:(t-S));
    I = sum(I_terms(:));
    
    for m = 1 : M
        % Draw from the one-step-ahead proposal
        % Draw a new history term from the proposal (transition)
        % distribution        
        
        %actually we don't use one-step-ahead for now, we use "prior
        %sampler" - just as well though - apprently just less efficient but
        %I haven't figured out how to sample from conditional for one-ahead
        h_mean = (1 - delta/tau) * h(:,t-1,m) + data(:,t-1);        
        % normrnd wastes time error-checking        
        %h(:,t,m) = normrnd(h_mean, sd);        
        h(:,t,m) = randn(N, 1) * sd + h_mean;

        J = b + I + w * h(:,t,m);

        % Compute probabilities
        % q = normpdf(h(:,t,m), h_mean, sd)
        emission_param = 1 - exp(-exp(J) * delta);

        if (isnan(emission_param))
          disp('emission is nan');
          disp(J);
        end
        if data(i,t)
            emission_prob = emission_param;
        else
            emission_prob = 1 - emission_param;
        end        
        pf(t, m) = emission_prob * pf(t-1, m);
        if(isnan(pf(t,m)))
          disp(emission_prob);
        end
        %[norm(h(:,t,m)) emission_prob pf(t,m)]
    end
        if(isnan(sum(pf(t,:))))
            disp('pf before summing');
            disp(pf(t,:));
             return
             end
    
    pf_new(t,:) = pf(t,:) / sum(pf(t,:));

    if(isnan(sum(pf_new(t,:))))
      disp('pf before resampling');
      disp(pf(t,:));
      return
    end
    
    pf = pf_new;
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

    if(isnan(sum(pf(t,:))))
         disp('pf after  resampling');
         disp(pf(t,:));
         return;
     end

    
%     if data(i,t)
%         lastspike = t;
%     end
%     if t <= lastspike + S
%         E = pf(t, :) * squeeze(h(:,t,:))' / M;
%         disp(['t = ' num2str(t) ' spike = ' num2str(data(i,t)) ' sample expectation norm ' num2str(norm(E))]);
%     end

end

% Marginal smoothing
% We don't need to save the r's; we just keep one timestep
r = zeros(M, M);

% At t = T, the filtering and smoothing distributions are identical.
pb = zeros(T,M);
pb(T,:) = pf(T,:);

% Determinant of a diagonal matrix is the product of its diagonal elements.
% But here our diagonal is also constant.
%normpdf_normalizer = (2*pi)^(N/2) * (N*sd^2)^(-.5);

for t_index = 0:(T-S-2)
    
    t = T - t_index;
    
    % Equations (12) and (13) in [Mischenko11]
    for m = 1 : M
        %sigma_mat = sd^2*eye(N);
%         prob = zeros(M,1);
%        num_terms = zeros(M,1);
%         denom = 0;
%         denom1 = 0;
        
        % KT -- vectorized the inner loop and verified equivalent results
        % N x M
%        h_means   = bsxfun(@plus, (1 - delta/tau) * squeeze(h(:,t-1,:)), data(:,t-1));
%        distances = bsxfun(@minus, h(:,t,m), h_means);
        h_squeezed = reshape(h(:,t-1,:), N, M);
        distances   = bsxfun(@plus, (1 - delta/tau) * h_squeezed, data(:,t-1) - h(:,t,m));
        distances_sq = sum(distances .^ 2,1);
        num_terms = exp(-0.5 * sd^-2 * distances_sq) .* pf(t-1,:);
        %num_terms = normpdf(-1*d
        if(isnan(sum(num_terms)))
          disp('num_terms');
          disp(num_terms);
          disp('distances');
          disp(distances);
          disp('distances_sq');
          disp(distances_sq);
          disp('pf');
          disp(pf(t-1,:));
         end
%         for mm = 1 : M
%             % hottest line (83/180 secs)
%             h_mean = (1 - delta/tau) * h(:,t-1,mm) + data(:,t-1);
%             % third hottest line (25/180 secs)
%             distance = h(:,t,m) - h_mean;
%             % Don't bother with matrices; we have sigma * I form covariance
%             % Don't even need normalization, since we explicitly normalize
%             % at the end anyways.
%             distance'*distance
%             num_terms(mm) = exp(-0.5 * sd^-2 * (distance'*distance)) * pf(t-1,mm);
%             
%             % BS TOOK OUT MVNPDF TO SAVE TIME... CAN BE SLOW, TOO MUCH
%             % OVERHEAD
% %             if (m == mm)
% %                 disp(['explicit:' num2str(mvnp)]);
% %                 disp(mvnpdf(h(:,t,m), h(:,t-1,mm), sd^2*eye(N)));
% %             end
% %             disp(['first test ' num2str(mvnp == mvnpdf(h(:,t,m), h(:,t-1,mm), sd^2*eye(N)))]);
% %             denom1 = denom1 +  mvnp(mm) * pf(t-1,mm);
% %             denom = denom + mvnpdf(h(:,t,m), h(:,t-1,mm), sd^2*eye(N)) * pf(t-1,mm);
%             
%         end        
        r(m,:) = pb(t,m) * num_terms / sum(num_terms);      
    end
    
    pb(t-1,:) = sum(r, 1);%sum over 1 or 2 here? I THINK 1 - BS
    
    pb(t-1,:) = pb(t-1,:) / sum(pb(t-1,:));
<<<<<<< HEAD
    if(isnan(sum(sum(pb))))
=======
    if(isnan(sum(sum(b))))
>>>>>>> a049df48e7ea1a24f9a575f3988dd7cb06eead9d
      disp('PB WEIGHT IS NAN!!!');
    end
end

sample_expectation_mean = zeros(1, T);
for t = 1 : T
    % h is N x T x M => squeeze(h(:,t,:))' is M * N
    sample_expectation_mean(t) = mean(pb(t,:) * reshape(h(:,t,:), M, N));
end
% figure
% hold on
% plot(sample_expectation_mean)
% scatter(1:T, data(i,:))
% drawnow

end
