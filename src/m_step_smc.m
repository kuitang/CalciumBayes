function [theta_intrinsic] = m_step_smc(theta_intrinsic, options, varargin)

    lb = ones(size(theta_intrinsic)) * -2;
    ub = ones(size(theta_intrinsic)) * 2;
    lb(1) = -5; lb(2) = -10;
    ub(1) = 5; ub(2) = 10;
    [theta_intrinsic,fval,exitflag,output] = fmincon('q_single_neuron',theta_intrinsic, [], [], [], [], lb, ub,[],options, varargin{:});

end

