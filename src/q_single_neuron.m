function [q g H] = q_single_neuron(theta_intrinsic, params, h, n, i, delta, tau, sigma, p_weights)

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

% Returns
% q         expectation at theta_intrinsic
% g         gradient at theta_intrinsic
% H         hessian at theta_intrinsic

beta = params.beta;
% lambda = params.lambda;
w = params.w;
S  = size(beta,3) + 1;
[N, T] = size(n);
M = size(h,3);
b_i = theta_intrinsic(1);
w(i,i) = theta_intrinsic(2);
beta(i,i,:) = reshape(theta_intrinsic(3:1+S),1,1,S-1);
disp(theta_intrinsic);
% lambda_i = theta_intrinsic(2+S:2*S+1);

reg_param1 = 1e1;
reg_param2 = 1e1;
q_sum = 0;

g = zeros(1, S + 1);
H = zeros(S + 1, S + 1);

disp('running objective function');

for t = S+1:T    
    % Partial derivatives of J
    % dJ(1) = dJ_i/db_i = 1
    % dJ(2) = DJ_i/dw_{ii} = h_{ii}(t)
    dJ = zeros(S + 1, 1);
    dJ(1) = 1;
    
    % Common gradient terms for this timestep (to multiply with dJ)
    % Both these techniques work due to distributive property
    dQ = 0;
    % Common Jacobian terms for this timestep (to multiply with dJ(n)*dJ(m))
    ddQ = 0;
    
    I = 0;
    for s = 2 : S
        dJ(1+s) = n(i,t-s);
        I = I + beta(i,:,s-1) * n(:,t-s);
    end
    
    for m = 1:M
        dJ(2) = dJ(2) + p_weights(t,m) * h(i,t,m);
        
            %squeeze or reshape here?
%         I = squeeze(beta(i,:,:),N,S-1) * reshape(n(:,t-S:t-2),N,S-1);
        
%         sum(sqeeze(beta(i,:,:),S,N) * ( 
        
        J = b_i + I + w(i,:) * h(:,t,m);
        
        eJd = exp(J)*delta;
        dQ1 = n(i,t) * 1/(1 - exp(-eJd)) * exp(-eJd) * eJd * delta;
        dQ2 = (1 - n(i,t)) * -eJd * delta;             
        dQ = dQ + p_weights(t,m) * (dQ1 + dQ2);
        
        ddQ1 = -eJd * delta * (1 - n(i,t));
        ddQ2 = exp(-2*eJd + 2*J) * delta^2 * n(i,t) / (1 - exp(-eJd))^2;
        ddQ3 = exp(-eJd + J) * delta * n(i,t) / (1 - exp(-eJd));
        ddQ4 = exp(-eJd + 2*J) * delta^2 * n(i,t) / (1 - exp(-eJd));
        ddQ = ddQ + p_weights(t,m) * (ddQ1 + ddQ2 + ddQ3 + ddQ4);        
                
%         J = b_i + I + w(i,i) * h(i,i,t,m);

        
%         sample_sum = sample_sum + n(i,t)*log(1 - exp(-exp(J)*delta)) + ... %if n(i,t) = 1
%             (1 - n(i,t))*log(exp(-exp(J)*delta)); %if n(i,t) = 0
%             
%         for j = 1:N
%            
%             history_mean = (1 - delta/tau)*h(i,j,t-1,m) + n(j,t - 1);
%             sample_sum = sample_sum + log(normpdf(h(i,j,t,m),history_mean, sigma));
%             
%         end
        
    %TOOK OUT normpdf BECAUSE NOT DEPENDENT ON PARAMS
%         history_mean = (1 - delta/tau).*h(i,:,t-1,m) + n(:,t - 1)';
%         disp(p_weights(t,m));
        sample_weighted = p_weights(t,m) * (n(i,t)*log(1 - exp(-exp(J)*delta)) + ... %if n(i,t) = 1
            (1 - n(i,t))*(-exp(J)*delta));% + ... %if n(i,t) = 0
%                 sum(log(normpdf(h(i,:,t,m),history_mean,sigma))));
      
        if(~isnan(sample_weighted))        
            q_sum = q_sum + sample_weighted;
        end
%         if (q_sum == -Inf)
%             break;
%         end

    end
    
    % Update the gradient and Hessian with information from this timeslice
    g = g + dQ * dJ;
    H = H + dQ * (dJ*dJ');
end

% Add regularization (to the variables only)
beta_vars = beta(i,i,:);
flat_beta_vars = beta_vars(:);

g = -g;
g(2) = g(2) + reg_param1 * sign(w);
g(3:end) = g(3:end) + reg_param2 * sign(flat_beta_subset);

H = -H;

% No regularization for H, since the L1 regularization terms have zero
% second derivative

q = -q_sum + reg_param1 * abs(w(i,i)) + reg_param2 * sum(abs(flat_beta_vars));

