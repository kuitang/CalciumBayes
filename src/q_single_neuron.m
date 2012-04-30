function q = q_single_neuron(theta_intrinsic, params, h, n, i, delta, tau, sigma, p_weights)

%Q_SINGLE_NEURON 
% This evaluates the negative of the q function of the EM algorithm for our
% model. It takes in all of the parameters and sampled and observed values
% and returns the value of q(theta,theta(old)).

% OPTIMIZED PARAMS
% b         the scalar parameter for this neuron's base firing rate
% w         n x 1 vector of this weights on this neuron from the other neurons
% beta      n x (S-1) matrix of the indirect weights on this neuron
% lambda     (S-1) x 1 vector of indirect weights from unobserved neurons
% tau       n x 1 vector of history decay constants for this neuron's history term
% sigma     std dev parameter for this neuron's history term

% CONSTANT PARAMS FOR THIS PART
% h         N X T X M matrix of latent history terms for this neuron for
% each particle
% n         N X T matrix of all observed spikes
% i         the index of this neuron



beta = params.beta;
lambda = params.lambda;
w = params.w;
[~, S]  = size(lambda);
[N, T] = size(n);
b_i = theta_intrinsic(1);
w_ii = theta_intrinsic(2);
beta_ii = theta_intrinsic(3:2+S);
lambda_i = theta_intrinsic(3+S:2*S+2);


reg_param = 1;
q_sum = 0;

for m = 1:M
    for t = S+1:T
        
        I = 0;
        if t - S > 0
            for s = 1 : S-1
               
                I = I + beta(i,:,s) * n(:,t-1) + lambda_i(s) - beta(i,i,s)*n(i,t-1) + ...
                    beta_ii(s)*n(i,t-1);
           
            end
        end 
        
        J = b_i + I + w(i,:) * h(i,:,t,m) - w(i,i) * h(i,i,t,m) + w_ii * h(i,i,t,m);
        
        q_sum = q_sum + n(i,t)*log(1 - exp(-exp(J)*delta)) + ... %if n(i,t) = 1
            (1 - n(i,t))*log(exp(-exp(J)*delta)); %if n(i,t) = 0
            
        for j = 1:N
           
            history_mean = (1 - delta/tau(j))*h(i,j,t,m) + n(j,t - delta);
            q_sum = q_sum + p_weights(t,l)*log(normpdf(history_mean, sigma(j)));
            
        end
    end
end

q = q_sum + reg_param * sum(w);

