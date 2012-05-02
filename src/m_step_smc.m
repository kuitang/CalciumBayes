function [theta_intrinsic] = m_step_smc(theta_intrinsic, params, h, n, i, delta, tau, sigma, p_weights)


    options = optimset('LargeScale','off');
    options = optimset(options,'Display','iter');
    options = optimset(options,'Algorithm','trust-region-reflective');
    options = optimset(options,'Hessian','user-supplied');
    options = optimset(options,'TolFun',.5);
    options = optimset(options,'MaxIter',3);

    lb = ones(size(theta_intrinsic)) * -.3;
    ub = ones(size(theta_intrinsic)) * .1;
    lb(1) = 0; lb(2) = -5;
    ub(1) = 1; ub(2) = 5;
    [theta_intrinsic,fval,exitflag,output] = fmincon('q_single_neuron',theta_intrinsic, [], [], [], [], lb, ub,[],options, params, h, n, i, delta, tau, sigma, p_weights)

end

