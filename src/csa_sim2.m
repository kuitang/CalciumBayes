%generate spike trains with variable firing rate from a poisson
%distribution

clear

V.dtI=.010;%input('dt: ');
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


% Vsave=('V.mat');
% Psave=('P.mat');
% save(Vsave,'V');
% save(Psave,'tau','gam','sr', 'lam', 'alpha', 'beta', 'snr', 'sig', 'gamI');
%%
Sim.P.S=2;%input(['Network Type:' '\n' '0= promiscuous' '\n' '1= daisy chained' '\n' '2= sparse' '\n']);
if isempty(Sim.P.S)
    Sim.P.S=2;
end
if Sim.P.S == 0
    Sim.P.W=input(['Weight Type:' '\n' '0=all or nothing network (plus/minus)' '\n' '1=random network' '\n' '2=neighbors sequentially connected' '\n' '3=25% inhibitory interneurons' '\n' '4=20% "excitatory" interneurons' '\n']);
    if isempty(Sim.P.W)
        Sim.P.W=1;
    end
elseif Sim.P.S == 1
    Sim.P.W=input(['Weight Type:' '\n' '5=random positive along daisy-chain' '\n' '6=plus/minus along daisy chain' '\n']);
    if isempty(Sim.P.W)
        Sim.P.W=6;
    end
elseif Sim.P.S == 2
    sparseness=input(['How connected are they? [0-1] (default .1)' '\n']);
    if isempty(sparseness)
        sparseness=.1;
    end
    Sim.P.W=5;
end
if Sim.P.W ~= 1
    Sim.P.WS=1;%input('Additional Strength of Selected Weights [0-1]: ');
    if isempty(Sim.P.WS)
        Sim.P.WS=.9;
    end
end

%weights
Sim.W=zeros(V.Ncells);
switch Sim.P.W
    
    case 0 %diagnostic all or nothing weights (plus/minus)
        Sim.W=round(2*(rand(V.Ncells)-.5));
    case 1 %random
        Sim.W=2*(rand(V.Ncells)-.5);
    case 2 %neighbors sequentially connected (plus)
        Sim.W=(1-Sim.P.WS)*(rand(V.Ncells)-.5);
        for i=1:V.Ncells-1
            Sim.W(i,i+1)=1;
        end
    case 3 % 20% inhibitory interneurons
        Sim.W=(1-Sim.P.WS)*(rand(V.Ncells));
        for i=1:5:V.Ncells
            Sim.W(i,:)=-1*rand(V.Ncells,1);
        end
    case 4 % 20% excitatory interneurons
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
    case 7 %plus or minus along daisy-chain
        dgnl=round((rand(V.Ncells-1,1)+1));
        for i=1:length(dgnl)
            if dgnl(i) ==2
                dgnl(i)=-1;
            end
        end
        Sim.W=diag(dgnl,1);
        
end

% for i=1:V.Ncells
%     Sim.W(i,i)=0;           %each cell has a weight of zero for itself
% end

% save('WeightMatrix.mat','Sim');
%%

H=zeros(V.T,V.Ncells);
FJ=zeros(V.T,V.Ncells);
J=zeros(V.T,V.Ncells);

FJ(1,:)=lam;

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


V.n=spikes;
%
% T=V.T*V.dtI;
% V.dt=.02;
% V.Tnew=T/V.dt;
% conv=V.dt/V.dtI;
% eps=normrnd(0,sig,[V.Tnew,V.Ncells]);
%
% traces=zeros(V.Tnew,V.Ncells);
%
% for i=1:V.Tnew-1
%     for j=1:V.Ncells
%     traces(i+1,j)=gam*traces(i,j)+sum(spikes((i-1)*conv+1:i*conv,j));
%     end
% end
% %
% F=((alpha*traces)+beta+eps)';

% F=(downsample(F,50))';
%
% [n p_best v]=fast_oopsi(F,V)