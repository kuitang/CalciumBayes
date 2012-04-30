load('../data/truncdata.mat')

% Physical parameters (TODO: Set up and figure out scale!)
% Currently using unit (discrete) time and bullshit values
sigma = 0.1;
tau = 1;
delta = 0.1;

% truncate the data even more!
truncdata = truncdata(:,1:100);

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

%until the change in connectivity matrix w is below threshold change
while(sum(sum(abs(w - w_prev))) > thresh_w)
    
    %parfor i = 1 : N
    for i = 1:N
        
        %theta_intrinsic = [b_i w_ii beta_ii(2:S)' lambda(1:S-1)]
        theta_intrinsic = [b(i) w(i,i) beta(i, i, :)
        
        while(theta_intrinsic - theta_intrinsic_prev > theta_intrinsic_thresh)
            
            disp(['Neuron ' num2str(i) '/' num2str(N)])
            
            [r h] = e_step_smc(i, M, tau, delta, sigma, beta, lambda, b, w, truncdata);
            m_step_smc(theta_intrinsic, 
            
        end
    end
    
    
       m_step_smc(b, w, beta, lambda, tau, sigma, h, truncdata, i, delta, M);
    
end