function [theta_intrinsic] = m_step_smc(theta_intrinsic, params, h, n, i, delta, tau, sigma, p_weights)


    options = optimset('LargeScale','off');
    [theta_intrinsic,fval,exitflag,output] = fminunc('q_signle_neuron',theta_intrinsic,options, params, h, n, i, delta, tau, sigma, p_weights);

end

