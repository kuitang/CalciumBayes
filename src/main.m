clear
run_parallel = 1;

%% use the same random number sequence for debugging purposes
RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',2));

%% Load data
% load('../data/truncdata.mat')
% truncate the data even more!
% truncdata = truncdata(:,1:100);
% n = truncdata;

load('good_sim_data_01.mat')
n = sim.n(1:10,:);

%% Set optimization options
optim_options = optimset('LargeScale','on','Algorithm', ...
    'trust-region-reflective','GradObj','on','Hessian','user-supplied', 'MaxIter',100);

%% Set physical parameters
% Physical parameters (TODO: Set up and figure out scale!)
% Currently using unit (discrete) time and bullshit values
%sigma = 0.01;
sigma = 1e-14;
tau = .020; %set to match simulator
delta = .010;


[N T] = size(n);
S = 20; % Indirect time window
M = 50; % size of particle sampler

% Parameter matrices
% TODO: Set up priors!

%% Set codistributed arrays
if run_parallel
    spmd(N)
        codist = codistributor1d(1);
        beta = zeros(N, N, S-1, codist);
        % lambda = ones(N, S);
        b = ones(N, 1, codist);
        w = zeros(N, N, codist) .* .1*exprnd(.9,N);
        p_weights = zeros(N,T,M,codist);
        % for i=1:N/5
        %     rand_inh
        %     w(i,:) = -exprnd(2.3,N,1);
        % end
        h = zeros(N,N,T,M, codist);

        w = w .* binornd(1,.1,N,N);%second arg is "sparesness"


%        for i = drange(1:N)
%            w(i,i) = -abs(normrnd(.6,.2));
%        end
    end
else
    %% Set noncodistributed arrays
    beta = zeros(N, N, S-1);
    % lambda = ones(N, S);
    b = zeros(N, 1);
    w = ones(N, N) .* .1*exprnd(.9,N);
    p_weights = zeros(N,T,M);
    % for i=1:N/5
    %     rand_inh
    %     w(i,:) = -exprnd(2.3,N,1);
    % end

    w = w .* binornd(1,.1,N,N);%second arg is "sparesness"

    h = zeros(N,N,T,M);

    for i = 1:N
        w(i,i) = -abs(normrnd(.6,.2));
    end
end
% load('first_e_step_complete.mat');
% first = 1;

%% Main loop
%until the change in connectivity matrix w is below threshold change

w_prev = ones(size(w)) * 500;
thresh_w = .001;

ll = -Inf;

while(norm(w - w_prev) > thresh_w)    
    
    
    disp(w);    
    disp(['NORM: ' num2str(norm(w - w_prev))]); 
    disp('********BEGINNING OF LOOP********');    
    disp('*********************************');
    w_prev = w; 
    
    
    spmd(N)  
        for i = drange(1:N)
            disp(['Neuron ' num2str(i) '/' num2str(N)]);            

            %% Initialize the intrinsic parameters
            theta_intrinsic = [b(i) w(i,i) reshape(beta(i, i, :),1,S-1)];
            old_theta_intr = ones(size(theta_intrinsic)) * 500;

            %% Let the intrinsic parameters converge
            while(theta_intrinsic(2) - old_theta_intr(2) > .1 || norm(theta_intrinsic([1 3:end]) - old_theta_intr([1 3:end])) > 0.01)
                
                disp(['NEURON ' num2str(i) ' NORM: ' num2str(norm(theta_intrinsic([1 3:end]) - old_theta_intr([1 3:end])))]);
                
                %% E step (SMC) for one neuron                
                old_theta_intr = theta_intrinsic;
                beta_subset = reshape(beta(i,:,:), N, S - 1);
                [p_weights(i,:,:) h(i,:,:,:)] = e_step_smc(i, M, tau, delta, sigma, beta_subset, b(i), w(i,:), n);

                %% M step for the intrinsic parameters for one neuron
                theta_intrinsic = m_step_smc(theta_intrinsic, optim_options, beta_subset, w(i,:), squeeze(h(i,:,:,:)), n, i, delta, tau, sigma, squeeze(p_weights(i,:,:)));                
                b(i,1) = theta_intrinsic(1);
                w(i,i) = theta_intrinsic(2);
                beta(i, i, :) = reshape(theta_intrinsic(3:1+S),1,1,S-1);
                disp('new params:');
                disp(theta_intrinsic);
            end
           disp(['NEURON ' num2str(i) ' DONE!']);
       end
    end
    spmd(N)
        for i = drange(1:N)
            %% M step for all parameters
            theta = [b(i) w(i,:) reshape(beta(i, :, :),1,N*S-N)];
            beta_subset = reshape(beta(i,:,:), N, S - 1);
            theta = m_step_full(theta, optim_options, N, beta_subset, w(i,:), squeeze(h(i,:,:,:)), n, i, delta, tau, sigma, squeeze(p_weights(i,:,:)));

            b(i,1) = theta(1);
            w(i,:) = reshape(theta(2:N+1),1,N);
            beta(i,:,:) = reshape(theta(N+2:end), 1, N, (S - 1));
        end
    end % spmd
    

    
    %% Log likelihood for whole model
    %nll = log_likelihood(beta, b, w, h, n, delta, p_weights);
    %ll = [ll nll];
    %disp('nll =');
    %disp(nll); 

end

w_gather = gather(w);
b_gather = gather(b);
beta_gather = gather(beta);
h_gather = gather(h);

save('finished_run.mat');


