function [sim] = simulator(N, T)

Sim.N = N;
Sim.T = T;
Sim.delta = 1;
Sim.tau = 20;
Sim.H=zeros(T,N);
f_J=zeros(T,N);
J=zeros(T,N);

f_J(1,:)=;

Sim.w = exprnd(.5,V.Ncells);

        Sim.W=(1-Sim.P.WS)*(rand(V.Ncells)-1);
        for i=1:5:V.Ncells
            Sim.W(i,:)=rand(V.Ncells,1);
        end
    case 5
        Sim.W=(1-Sim.P.WS)*(exprnd(.5,V.Ncells));
        for i=1:V.Ncells/5
            Sim.W(i,:)=-exprnd(2.3,V.Ncells,1);
        end
        Sim.W=Sim.W.*(binornd(1,sparseness,V.Ncells,V.Ncells));
    case 6 %random plus along daisy-chain
        dgnl=(rand(V.Ncells-1,1));
        Sim.W=diag(dgnl,1);

spikes=zeros(V.T,V.Ncells);
spikes(1,:)=binornd(1,FJ(1,:),[1,V.Ncells]);

for i=2:V.T
    H(i,:)=gam*H(i-1,:)+spikes(i-1,:);
    for j=1:V.Ncells
        J(i,j)=beta(j)+Sim.W(:,j)'*H(i,:)';%spikes(i-1,:)')+J(i-1,j)');
        FJ(i,j)=1-exp((-exp(J(i,j))*V.dtI));
    end
    spikes(i,:)=binornd(1,FJ(i,:),[1,V.Ncells]);
end


V.T=input('number of time steps: ');
V.Ncells=input('number of cells: ');
tau=.02; %PSP decay time
% tauI=.02; % inhibitory PSP decay time
% tauE=.01; %excitatory PSP decay time
gam=1-(V.dtI/tau); %gam=1-(V.dtI/tauI); gamE=1-(V.dtI/tauE);
sr=2;%input('mean SR: ');
lam=exprnd((sr*V.dtI),[1,V.Ncells]);
% alpha=1;
beta=log(-log(1-lam)/V.dtI);
% snr=input('SNR: ');
% sig=1/snr;




end