function ll = log_likelihood(params, h, n, delta, p_weights)

%log(P(data | params))

beta = params.beta;
% lambda = params.lambda;
w = params.w;
b = params.b;
S  = size(beta,3) + 1;
[N, T] = size(n);
M = size(h,4);

ll = 0;

for i = 1:N
    
  for t = S+1:T
      
    n_t_prob = 0;
    
    for m = 1:M
        
        I = 0;
        
        for s = 1 : S-1
               
            I = I + beta(i,:,s) * n(:,t-(s+1));
           
        end 
        J = b(i) + I + reshape(w(i,:),1,N) * reshape(h(i,:,t,m),N,1);
        
%         disp(m)
%         disp(p_weights(t,m));
%         disp('calc=')
%         disp(((1 - exp(-exp(J)*delta))^(n(i,t)) * ... 
%             (1 - (1 - exp(-exp(J)*delta)))^(1 - n(i,t))));
        n_t_prob = n_t_prob + p_weights(t,m)*((1 - exp(-exp(J)*delta))^(n(i,t)) * ... %if n(i,t) = 1
            (1 - (1 - exp(-exp(J)*delta)))^(1 - n(i,t))); %if n(i,t) = 0
            
    end
    
    ll = ll + log(n_t_prob);
    
  end
  
end

end