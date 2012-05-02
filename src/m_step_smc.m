function [theta_intrinsic] = m_step_smc(theta_intrinsic, params, h, n, i, delta, tau, sigma, p_weights)


    options = optimset('LargeScale','on'); 
    options = optimset(options,'Display','iter');
    options = optimset(options,'Algorithm','trust-region-reflective');
    options = optimset(options,'GradObj','on');
    options = optimset(options,'Hessian','user-supplied');    
    options = optimset(options,'TolFun',.001);
%     options = optimset(options,'MaxIter',100);
    %DRASTIC SLOWDOWN -- run this on cluster machine!
    options = optimset(options,'DerivativeCheck', 'on');

    lb = ones(size(theta_intrinsic)) * -2;
    ub = ones(size(theta_intrinsic)) * 2;
    lb(1) = -5; lb(2) = -10;
    ub(1) = 5; ub(2) = 10;
    [theta_intrinsic,fval,exitflag,output] = fmincon('q_single_neuron',theta_intrinsic, [], [], [], [], lb, ub,[],options, params, h, n, i, delta, tau, sigma, p_weights)

end

