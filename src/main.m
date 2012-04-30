% % %use the same random number sequence for debugging purposes
RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',2));

% load('../data/truncdata.mat')
% truncate the data even more!
% truncdata = truncdata(:,1:100);
% n = truncdata;

load('../data/testdata.mat')
n = spikes';


% Physical parameters (TODO: Set up and figure out scale!)
% Currently using unit (discrete) time and bullshit values
sigma = 0.5;
tau = 30;
delta = 1;



[N T] = size(n);
S = 20; % Indirect time window
M = 50; % size of particle sampler

% Parameter matrices
% TODO: Set up priors!
params.beta = zeros(N, N, S-1); % DO WE NEED BETA (N, N, 1)...
params.lambda = zeros(N, S);
params.b = zeros(1, N);
params.w = zeros(N, N);
params.w = .1*exprnd(.5,N);
for i=1:N/5
    params.w(i,:) = -exprnd(2.3,N,1);
end
params.w = params.w .* binornd(1,.1,N,N);%second arg is "sparesness"


h = zeros(1,N,T,M);

theta_intrinsic_thresh = .01; %???

beta = params.beta;
lambda = params.lambda;
w = params.w;
b = params.b;

%until the change in connectivity matrix w is below threshold change
% while(sum(sum(abs(w - w_prev))) > thresh_w)
    
    %parfor i = 1 : N
    for i = 1:N

        ll = -Inf;
        
        %theta_intrinsic = [b_i w_ii beta_ii(2:S)' lambda(1:S-1)]
        theta_intrinsic = [b(i) w(i,i) reshape(beta(i, i, :),1,S-1) lambda(i,:)];
        
        %e
        [p_weights h(i,:,:,:)] = e_step_smc(i, M, tau, delta, sigma, params, n);
        save('first_e_step_complete.mat');%save this if we complete an e_step so we can test m_step
        %m
        theta_intrinsic = m_step_smc(theta_intrinsic, params, h, n, i, delta, tau, sigma, p_weights);
        params.b(i) = theta_intrinsic(1);
        params.w(i,i) = theta_intrinsic(2);
        params.beta(i, i, :) = theta_intrinsic(3:2+S);
        params.lambda(i,:) = theta_intrinsic(3+S:2*S+2);
        
        %compute new log-likelihood
        nll = log_likelihood(params, h, n, delta);
        disp(['ll =' num2str(nll)]);
        
        while(ll + theta_intrinsic_thresh > nll)
            
            disp(['Neuron ' num2str(i) '/' num2str(N)]);
            
            %e
            [r h] = e_step_smc(i, M, tau, delta, sigma, params, n);
            %m
            m_step_smc(theta_intrinsic, params, h, n, i, delta, p_weights);
            params.b(i) = theta_intrinsic(1);
            params.w(i,i) = theta_intrinsic(2);
            params.beta(i, i, :) = theta_intrinsic(3:2+S);
            params.lambda(i,:) = theta_intrinsic(3+S:2*S+2);
            
        end
    end
    
    
%     m_step_joint(b, w, beta, lambda, tau, sigma, h, truncdata, i, delta, M);
    
% end



