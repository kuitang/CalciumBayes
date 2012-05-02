function Sim = simulator(N, T)

Sim.sparesness = .1;
Sim.N = N;
Sim.T = T;
Sim.delta = .001;
Sim.tau = .020;
Sim.sigma = .05;
Sim.h=zeros(N,N,T);
Sim.f_J = zeros(N,T);
Sim.J=zeros(N,T);

Sim.lambda = exprnd(.85,[N,1]);
Sim.b = normrnd(5,.5,[N,1]);


Sim.f_J(:,1)= Sim.lambda;
% Sim.f_J(:,1) = zeros(N,1); % if you want no input

for i = 1:N
    if (Sim.f_J(i,1) > 1)
        Sim.f_J(i,1) = .95;
    end
end



Sim.w = exprnd(.5,N);

for i = 1:5:N    
    Sim.w(:,i) = -exprnd(2.3,[N,1]);
end

Sim.w = Sim.w .* (binornd(1,Sim.sparesness,N,N));
for i = 1:N   
    Sim.w(i,i) = -abs(normrnd(.6,.2));
end

Sim.n = zeros(N,T);
for i = 1:N
    Sim.n(i,1) = binornd(1,Sim.f_J(i,1));
end

Sim.h(:,:,1) = normrnd(0,Sim.sigma);

gam = 1-(Sim.delta/Sim.tau);


for t = 2:T
    for i = 1:N    
        
        Sim.J(i,t) = Sim.b(i);
        
        for j = 1:N
            
            Sim.h(i,j,t) = gam * Sim.h(i,j,t-1) + Sim.n(j,t-1) + Sim.sigma * sqrt(Sim.delta) * normrnd(0,1);
            Sim.J(i,t) = Sim.J(i,t) + Sim.w(i,j) * Sim.h(i,j,t);

        end
        
    Sim.f_J(i,t) = 1 - exp(-exp(Sim.J(i,t))*Sim.delta);
    Sim.n(i,t) = binornd(1, Sim.f_J(i,t));
    
    end
end

end