function [ output_args ] = m_step_full(theta_intrinsic, options, N, varargin)
%M_STEP_FULL Summary of this function goes here
%   Detailed explanation goes here


    lb = ones(size(theta_intrinsic)) * -2;
    ub = ones(size(theta_intrinsic)) * 2;
    lb(1) = -5; lb(2:N+1) = -10;
    ub(1) = 5; ub(2:N+1) = 10;
    [theta_intrinsic,fval,exitflag,output] = fmincon('q_single_neuron_full',theta_intrinsic, [], [], [], [], lb, ub,[],options, varargin{:});

end

