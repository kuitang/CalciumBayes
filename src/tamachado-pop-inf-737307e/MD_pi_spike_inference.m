function [spikes struct] = MD_pi_spike_inference(signals,V,spikes)
% inference = pi_spike_inference(inputs):
%
% inputs (only first argument is required):
%  signal      - fluorescence traces in a matrix (nNeurons x nFrames)
%  V.epath     - path to associated .paq file containing ephys data (path)
%  V.alpha     - threshold for extracting spikes from voltage data (0-1)
%  V.flip      - boolean indicating whether signal should be flipped (boolean)
%  V.channel   - headstage channel containing voltage data (1-4)
%  V.figure    - figure handle to plot stuff in
%  V.indices   - a logical vector showing selected neurons
%  V.fast      - use fast or smc
%  V.inf_plot  - generate spike inference plot
%  spikes      - pass in spike inference here instead of recomputing
%
% output:
%  the chosen spike inference algorithm is run on the input
%  signal. the ephys data is aligned with the fluorescence data.
%  the results are plotted on screen to a new figure.
%
%  inference - vector of length(signal) containing spike inference
%
%  tamachado (5/10)



%% MD TESTING COMBO of Mor, Tanya and Tims codes....

%% RUN FOOPSI
fr = input('Frame Rate (Hz): ');  % V has all sorts of unnecessary stuff by now....
T.dt = 1/fr;
P.lam=input('estimated firing rate: ');
tau=3.5;
P.gam = (1-T.dt/tau)';
rawF = signals;

spikes = cell(size(rawF,1),1); % spike inference output vector
struct = cell(size(rawF,1),1); % spike inference output struct
indices = ones(size(rawF,1),1);

Ncells = size(rawF,1);
n = zeros(size(rawF,2),1);

for i = 1:size(rawF, 1)
    %  i=17;
    fprintf('\nneuron %d\n',i);
    Fcell = rawF(i,:);
    Fcell = preprocessdata(Fcell);
%     Fcell=mapab(Fcell,eps,1);
    P.b = median(Fcell);
    fprintf('running spike inference...\n');
    [n test.P test.V]= fast_oopsi(Fcell,T,P);
    struct{i}.P = test.P;
    struct{i}.n = n;
    spikes{i}=struct{i}.n;
    spikes{i}=spikes{i}/max(spikes{i});
    V.F{i}=Fcell(1,:);
end

% make a new figure if desired
if ~isfield(V,'figure'),  handle = figure; else handle=figure(V.figure); end

%% draw the gui
if V.inf_plot
    if ~exist('times','var'), times = 1:size(rawF,2); end;
    figure(handle);
    k=1; kmin=1; kmax=Ncells;
    plot_callback;
    set(gcf,'Color','w','Toolbar','figure');
    guidata(handle,indices);
    if Ncells > 1
        hb = uicontrol(...
            'Style', 'togglebutton',...
            'String', 'Exclude',...
            'Units','normalized',...
            'Position', [0 0 .1 .04],...
            'Callback',@clicked_callback);
        ha = uicontrol(gcf,...
            'Style','slider',...
            'Min' ,kmin,'Max',kmax,...
            'Units','normalized',...
            'Position',[.1 0 .9 .04],...
            'Value', k,...
            'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
            'Callback',@plot_callback);
        % wait until the window is closed before exiting
        uiwait;
    end
end

%% make the interactive plotting window
    function k=plot_callback(handle, eventdata, handles) %#ok<INUSD>
        
        % move the scroll bar
        if exist('ha','var')
            indices = guidata(gcbo);
            k = round(get(ha,'Value'));
        else
            k = 1;
        end
        
        % truncate fluorescence if there are too few camera triggers
        if length(times) < length(V.F{k}) && ~isempty(times)
            V.F{k} = V.F{k}(1:length(times));
            spikes{k} = spikes{k}(1:length(times));
        end
        
        % plot spike times and fluorescence
        ax(1) = subplot(2,1,1);
        cla; plot(V.F{k}./max(V.F{k}),'k');
        mstats = 0;
        
        % indicate whether the current neuron is excluded
        if ~exist('indices','var') || length(indices) < k || indices(k) == 1
            tt = sprintf('neuron %d',k);
            eSnr = struct{k}.P.a/struct{k}.P.sig;
            if mstats == 1
                title(sprintf('%s, AUC %.2f, eSNR %.2f',tt,roc.AUC,eSnr),'FontSize',14);
            else
                title(sprintf('%s, eSNR %.2f',tt,eSnr),'FontSize',14);
            end
        else
            title(sprintf('neuron %d: excluded',k),'FontSize',14);
        end
        
        % set x axis appropriately
        nTicks = 40;
        xt = [1:round(length(V.F{k})/nTicks):length(V.F{k})];
        if numel(times) > 0
            xl = round(times(xt) - times(1));
            set(gca,'XTick',xt,'YTick',[],'XTickLabel',xl,'YTickLabel',[]);
            xlabel('time (s)');
        end
        
        % plot inference output
        ax(2) = subplot(2,1,2);
        cla; bar(spikes{k},'k');
        
        % link the axes
        linkaxes(ax,'x');
        
    end

    function c=clicked_callback(handle, eventdata, handles) %#ok<INUSD>
        % flip the bit at the appropriate index
        c = get(hb,'Value');
        k = get(ha,'Value');
        k = round(k);
        % update the data
        indices = guidata(gcbo);
        % save out our modifications
        if indices(k), indices(k) = 0; else indices(k) = 1; end;
        guidata(gcbo,indices);
        plot_callback;
        assignin('base','indices', indices);
    end
end