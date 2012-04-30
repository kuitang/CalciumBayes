load('../data/truncdata.mat')

% Physical parameters (TODO: Set up and figure out scale!)
% Currently using unit (discrete) time and bullshit values
sigma = 0.1;
tau = 1;
delta = 0.1;

% truncate the data even more!
truncdata = truncdata(:,1:100);
n = truncdata;

[N T] = size(truncdata);
S = 20; % Indirect time window
M = 50; % size of particle sampler

% Parameter matrices
% TODO: Set up priors!
beta = zeros(N, N, S); % DO WE NEED BETA (N, N, 1)... no...
lambda = zeros(N, S);
b = zeros(1, N);
w = zeros(N, N);
w_prev = w;
h = zeros(N,N,T,M);


%until the change in connectivity matrix w is below threshold change
while(sum(sum(abs(w - w_prev))) > thresh_w)
    
    %parfor i = 1 : N
    for i = 1:N

        ll = -Inf;
        
        %theta_intrinsic = [b_i w_ii beta_ii(2:S)' lambda(1:S-1)]
        theta_intrinsic = [b(i) w(i,i) reshape(beta(i, i, :),1,S-1) lambda(i,:)];
        theta_intrinsic_prev = ones(size(theta_intrinsic)) * -Inf;
        
        %e
        [r h] = e_step_smc(i, M, tau, delta, sigma, beta, lambda, b, w, n);
        %m
        m_step_smc(theta_intrinsic, i, M, tau, delta, sigma, beta, lambda, b, w, n);
        
        %compute new log-likelihood
        nll = log_likelihood(b, w, beta, lambda, h, truncdata, delta);
        disp(['nll =' num2str(nll)]);
        
        % SHOULD WE COMPUTE LOG-LIKELIHOOD HERE FOR CONVERGENCE?
        while(ll + theta_intrinsic_thresh > nll)
            
            disp(['Neuron ' num2str(i) '/' num2str(N)]);
            
            %e
            [r h] = e_step_smc(i, M, tau, delta, sigma, beta, lambda, b, w, n);
            %m
            m_step_smc(theta_intrinsic, i, M, tau, delta, sigma, beta, lambda, b, w, n);
            
        end
    end
    
    
    m_step_joint(b, w, beta, lambda, tau, sigma, h, truncdata, i, delta, M);
    
end



