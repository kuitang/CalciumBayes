function q = q_single_neuron(b, w, beta, gamma, tau, sigma, h, n, i,delta)
%Q_SINGLE_NEURON 
% This evaluates the q function of the EM algorithm for our
% model. It takes in all of the parameters and sampled and observed values
% and returns optimzed parameter values b, w, beta, gamma, tau, and sigma.

% OPTIMIZED PARAMS
% b         the scalar parameter for this neuron's base firing rate
% w         n x 1 vector of this weights on this neuron from the other neurons
% beta      n x (S-1) matrix of the indirect weights on this neuron
% gamma     (S-1) x 1 vector of indirect weights from unobserved neurons
% tau       n x 1 vector of history decay constants for this neuron's history term
% sigma     std dev parameter for this neuron's history term

% CONSTANT PARAMS FOR THIS PART
% h         N X T X M matrix of latent history terms for this neuron for
% each particle
% n         N X T matrix of all observed spikes
% i         the index of this neuron


q_sum = 0;

for l = 1:M
    for t = S+1
        
        I = 0;
        if t - S > 0
            for s = 2 : S
                I = I + beta(i,:,s) * data(:,t-1) + lambda(i,s);
            end
        end 
        
        J = b + I + w * h(:,t,m);
        
        q_sum = q_sum + n(i,t)*log(1 - exp(-exp(J)*delta) + ... %if n(i,t) = 1
            (1 - n(i,t))*log(exp(-exp(J)*delta))l %if n(i,t) = 0
            
        for j = 1:N
           
            history_mean = (1 - delta/tau(j))*h(j,t,m) + n(j,t - delta);
            q_sum = q_sum + log(normpdf(history_mean, sigma(j)));
            
        end
    end
end

