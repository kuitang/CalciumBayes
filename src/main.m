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
params.beta = zeros(N, N, S-1);
% params.lambda = ones(N, S);
params.b = zeros(1, N);
params.w = zeros(N, N);
params.w = .1*exprnd(.9,N);
% for i=1:N/5
%     rand_inh
%     params.w(i,:) = -exprnd(2.3,N,1);
% end

params.w = params.w .* binornd(1,.1,N,N);%second arg is "sparesness"
for i = 1:N   
    params.w(i,i) = -abs(normrnd(.6,.2));
end

theta_intrinsic_thresh = .01; %???

beta = params.beta;
% lambda = params.lambda;
w = params.w;
b = params.b;

h = zeros(N,N,T,M);

% load('first_e_step_complete.mat');
% first = 1;

%% Main loop
%until the change in connectivity matrix w is below threshold change
% while(sum(sum(abs(w - w_prev))) > thresh_w)
    
    %parfor i = 1 : N
    for i = 1:N
                
        disp(['Neuron ' num2str(i) '/' num2str(N)]);

        ll = -Inf;
        iter = 1;
        
        %theta_intrinsic = [b_i w_ii beta_ii(2:S)' lambda(1:S-1)] % NO
        %LAMBDA
        %% Initialize the intrinsic parameters
        theta_intrinsic = [b(i) w(i,i) reshape(beta(i, i, :),1,S-1)];
        old_theta_intr = ones(size(theta_intrinsic)) * 500;
        
        %% Let the intrinsic parameters converge
        while( any(abs(theta_intrinsic - old_theta_intr) > .1))
            %% E step (SMC) for one neuron
            iter = iter + 1;
            old_theta_intr = theta_intrinsic;
            [p_weights h(i,:,:,:)] = e_step_smc(i, M, tau, delta, sigma, params, n);

    %         if(~first)
    %             [p_weights h(i,:,:,:)] = e_step_smc(i, M, tau, delta, sigma, params, n);
    %             first = 0;
    %         end
    %         save('new_e_step_complete.mat');%save this if we complete an
    %         e_step so we can test m_step
    %         load('first_e_step_complete.mat');                            
                                    
            %m
            
            %% M step for the intrinsic parameters for one neuron
            theta_intrinsic = m_step_smc(theta_intrinsic, params, squeeze(h(i,:,:,:)), n, i, delta, tau, sigma, p_weights);
            params.b(i) = theta_intrinsic(1);
            params.w(i,i) = theta_intrinsic(2);
            params.beta(i, i, :) = reshape(theta_intrinsic(3:1+S),1,1,S-1);
            disp('new params:');
            disp(theta_intrinsic);
            disp('abs diff of params');
            disp(abs(theta_intrinsic - old_theta_intr));
%             params.lambda(i,:) = theta_intrinsic(3+S:2*S+2);

%             nll = log_likelihood(params, h, n, delta, p_weights);
%             ll = [ll nll];
%             disp('ll =');
%             disp(ll);
            
        end
                
        %e

        
        %% Compute the new log-likelihood        
%         save('new_m_step_complete.mat');
                
        nll = log_likelihood(params, h, n, delta, p_weights);
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
                  

    end
    
    %% blah blah blah
    save('finished_run.mat');
    
    
%     m_step_joint(b, w, beta, lambda, tau, sigma, h, truncdata, i, delta, M);
    
% end



