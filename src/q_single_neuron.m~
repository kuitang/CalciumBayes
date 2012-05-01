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
% lambda = params.lambda;
w = params.w;
S  = size(beta,3) + 1;
[N, T] = size(n);
M = size(h,4);
b_i = theta_intrinsic(1);
w(i,i) = theta_intrinsic(2);
beta(i,i,:) = reshape(theta_intrinsic(3:1+S),1,1,S-1);
disp(theta_intrinsic);
% lambda_i = theta_intrinsic(2+S:2*S+1);


reg_param1 = 1e1;
reg_param2 = 1e1;
q_sum = 0;

disp('running objective function');

for m = 1:M
    
%     sample_sum = 0;
    
    for t = S+1:T
        
        I = 0;

        for s = 1 : S-1
               
            I = I + beta(i,:,s) * n(:,t-(s+1));
           
        end

            %squeeze or reshape here?
%         I = squeeze(beta(i,:,:),N,S-1) * reshape(n(:,t-S:t-2),N,S-1);
        
        
        
%         sum(sqeeze(beta(i,:,:),S,N) * ( 

        
        J = b_i + I + reshape(w(i,:),1,N) * reshape(h(i,:,t,m),N,1);
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
end

q = -q_sum + reg_param1 * sum(sum(abs(w))) + reg_param2 * sum(sum(sum(abs(beta))))

