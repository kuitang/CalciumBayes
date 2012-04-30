function [theta_intrinsic] = m_step_smc(theta_intrinsic, params, h, n, i, delta, tau, sigma, p_weights)


    options = optimset('LargeScale','on');
    options = optimset(options,'Display','iter');
    options = optimset(options,'Algorithm','interior-point');
    lb = ones(size(theta_intrinsic)) * -4;
    ub = ones(size(theta_intrinsic)) * 4;
    [theta_intrinsic,fval,exitflag,output] = fmincon('q_single_neuron',theta_intrinsic, [], [], [], [], lb, ub,[],options, params, h, n, i, delta, tau, sigma, p_weights);

end

