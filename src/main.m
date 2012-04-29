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
beta = zeros(N, N, S);
lambda = zeros(N, S);
b = zeros(1, N);
w = zeros(N, N);
w_prev = w;

%until the change in connectivity matrix w is below threshold change
while(sum(sum(abs(w - w_prev))) > thresh_w)
    
    %parfor i = 1 : N
    for i = 1:N
        while(intrinsic_params - intrinsic_params_prev > thresh_intrinsic_params)
            disp(['Neuron ' num2str(i) '/' num2str(N)])
            e_step_smc(i, M, tau, delta, sigma, beta, lambda, b, w, truncdata);
            m_step_smc(i, params);
        end
    end
end