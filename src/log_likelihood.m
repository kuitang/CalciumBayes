function ll = log_likelihood(params, h, n, delta)

%log(P(data | params))

beta = params.beta;
lambda = params.lambda;
w = params.w;
b = params.b;

ll = 0;

for i = 1:N
    for t = S+1:T
        
        I = 0;
        for s = 2 : S
            I = I + beta(i,:,s) * n(i,:,t-1) + lambda(s);
        end
        
        J = b + I + w * h(i,:,t,m);
        
        ll = ll + n(i,t)*log(1 - exp(-exp(J)*delta)) + ... %if n(i,t) = 1
            (1 - n(i,t))*log(exp(-exp(J)*delta)); %if n(i,t) = 0
            
    end
end

end