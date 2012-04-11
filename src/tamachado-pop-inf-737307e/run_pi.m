% run_pi
%
% this script runs population inference on raw fluorescence traces
% and displays the results
%
% in addition to the estimated connectivity matrix, the spike rate
% nonlinearity for each neuron is plotted. this function should be
% increasing (ie the predicted firing rate given the model is correlated
% with the actual firing rate of that neuron). run the code on the 
% simulated data to get an example for what this should look like.
%
% tamachado 5/10

% clear the workspace
% clear all;

% add the subfolders to the path
% addpath([pwd '/oopsi'])
% addpath([pwd '/glmnet/'])

% set up parameters (or pi_infer_connectivity will prompt for params)

        % MD - Data
        V.flip = false;
        [f p]=uigetfile('D:\PiperData\Movies','MultiSelect','on');

        % If only a single trace file is selected, load traces.
        if ischar(f)
            V.path = [p f];
            temp = load(V.path);
            name = fieldnames(temp);
            name = name{1};
            V.traces = temp.(name);
        else
            % If multiple trace files are selected combine them and load
            % them.  BE AWARE - this combines them in sequential order!
            for i=1:length(f)
                path = strcat(p,f(i));
                temp = load(path{1});
                name = fieldnames(temp);
                name = name{1};
                numcells = size(temp.(name),1);
                numframes = length(temp.(name));
                if i ==1
                    traces = temp.(name); 
                else
                    traces(1:numcells,length(traces)+1:length(traces)+numframes) = temp.(name);
                end

            end
            V.traces = traces;
        end
            
            V.inf_plot = 1;
            numcells=size(V.traces,1);
%             V.posspiketimes = 1000:1000:size(V.traces,2);
%             V.posspiketimes = 1:200:size(V.traces,2);
%             V.dt = 0.0222;
        
% run connectivity inference
% O = pi_infer_connectivity(V);
O = MD_pi_infer_connectivity(V);

% plot spike inference summary plot
if exist('indices','var'), O.indices = indices; disp(indices'), end
pi_plot_inference(O.F,O.N,O.indices);

% plot connectivity matrix for lam_max
O.Phat.omega=O.Phat.omega';
figure; imagesc(O.Phat.omega); colormap gray;
title('connectivity matrix'); xlabel('postsynaptic'); ylabel('presynaptic');

% plot nonlinearity to show predictive power of model fit
% pi_plot_nonlinearity(O);
% clear;

[maxpre maxpost] = find(O.Phat.omega==max(O.Phat.omega(:)))
[minpre minpost] = find(O.Phat.omega==min(O.Phat.omega(:)))