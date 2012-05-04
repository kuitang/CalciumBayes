function [ theta ] = m_step_full(theta, options, N, varargin)
%M_STEP_FULL Summary of this function goes here
%   Detailed explanation goes here


    lb = ones(size(theta)) * -5;
    ub = ones(size(theta)) * 5;
    lb(1) = 0; lb(2:N+1) = -5;
    ub(1) = 5; ub(2:N+1) = 5;
    [theta,fval,exitflag,output] = fmincon('q_single_neuron_full',theta, [], [], [], [], lb, ub,[],options, varargin{:});

end

