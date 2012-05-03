function ll = log_likelihood(beta, b, w, h, n, delta, p_weights)

% error('TODO: Fix the log_likelihood interface to only use the history
% terms of one neuron at a time')
% isn't this already the case??

%log(P(data | params))

S  = size(beta,3) + 1;
[N, T] = size(n);
M = size(h,4);

spmd
    llv = codistributed.zeros(1, N);
    for i = drange(1:N)
        beta_i = reshape(beta(i,:,:), N, S - 1);
        w_i    = reshape(w(i,:), 1, N);        

        for t = S+1:T
                        
            I_terms = beta_i .* n(:,(t-2):-1:(t-S));
            I = sum(I_terms(:));
            
            for m = 1:M                
                J = b(i) + I + w_i * h(i,:,t,m)';
                eJd = exp(J)*delta;                
                if n(i,t)
                    eeJd = exp(-eJd);
                    Qm = log(1 - eeJd);                    
                else
                    Qm = -eJd;                    
                end
                llv(i) = llv(i) + p_weights(t,m) * Qm;
            end
        end                
    end
end
ll = sum(llv);

end
