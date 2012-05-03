clear

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
optim_options = optimset('LargeScale','on','Dispay','iter','Algorithm', ...
    'trust-region-reflective','GradObj','on','Hessian','user-supplied'); 

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

spmd(N)
    codist = codistributor1d(1);
    beta = zeros(N, N, S-1, codist);
    % lambda = ones(N, S);
    b = zeros(N, 1, codist);
    w = ones(N, N, codist) .* .1*exprnd(.9,N);
    % for i=1:N/5
    %     rand_inh
    %     w(i,:) = -exprnd(2.3,N,1);
    % end

    w = w .* binornd(1,.1,N,N);%second arg is "sparesness"

    h = zeros(N,N,T,M, codist);

    for i = drange(1:N)
        w(i,i) = -abs(normrnd(.6,.2));
    end
end

% load('first_e_step_complete.mat');
% first = 1;

theta_intrinsic = cell(N,1);

for i = 1:N
    
    theta_intrinsic{i} = [b(i) w(i,i) reshape(beta(i, i, :),1,S-1)];
    
end

%% Main loop
%until the change in connectivity matrix w is below threshold change
% while(norm(w - w_prev) > thresh_w)
    

    %parfor i = 1 : N
    ll = -Inf;
    spmd(N)        
        for i = drange(1:N)

            disp(['Neuron ' num2str(i) '/' num2str(N)]);            

            
            iter = 1;

            %theta_intrinsic = [b_i w_ii beta_ii(2:S)' lambda(1:S-1)] % NO
            %LAMBDA
            %% Initialize the intrinsic parameters

            theta_intrinsic = [b(i) w(i,i) reshape(beta(i, i, :),1,S-1)];
            old_theta_intr = ones(size(theta_intrinsic)) * 500;

            %% Let the intrinsic parameters converge
            while( norm(theta_intrinsic - old_theta_intr) > 1e-4)
                %% E step (SMC) for one neuron
                iter = iter + 1;
                old_theta_intr = theta_intrinsic;
                beta_subset = reshape(beta(i,:,:), N, S - 1);
                [p_weights h(i,:,:,:)] = e_step_smc(i, M, tau, delta, sigma, beta_subset, b(i), w(i,:), n);

                %% M step for the intrinsic parameters for one neuron
                theta_intrinsic = m_step_smc(theta_intrinsic, optim_options, beta_subset, w(i,:), squeeze(h(i,:,:,:)), n, i, delta, tau, sigma, p_weights);                
                b(i,1) = theta_intrinsic(1);
                w(i,i) = theta_intrinsic(2);
                beta(i, i, :) = reshape(theta_intrinsic(3:1+S),1,1,S-1);
                disp('new params:');
                disp(theta_intrinsic);
                disp('abs diff of params');
                disp(abs(theta_intrinsic - old_theta_intr));
            end
            
            %% M step for all parameters
            
            %% Compute the new log-likelihood
    %         save('new_m_step_complete.mat');
        end
    end
    
    nll = log_likelihood(beta, b, w, h, n, delta, ones(T,M));
    %nll = log_likelihood(gather(beta), gather(b), gather(w), gather(h), n, delta, ones(T,M));
    ll = [ll nll];
    disp('ll =');
    disp(ll);
%         disp('******');
%         disp(theta_intrinsic);
%         disp(old_theta_intr);
%         disp(theta_intrinsic - old_theta_intr);
%         disp(abs(theta_intrinsic - old_theta_intr));
%         disp(sum(abs(theta_intrinsic - old_theta_intr)));

%         while(sum(abs(theta_intrinsic - old_theta_intr)) > .01) %we shoudl make this bigger...
        %... maybe this: abs(theta_intrinsic - old_theta_intr) > .01
        %... that is, don't sum, just test each one



    
    %% blah blah blah
    save('finished_run.mat');
    
    
%     m_step_joint(b, w, beta, lambda, tau, sigma, h, truncdata, i, delta, M);
    



