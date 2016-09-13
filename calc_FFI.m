function win_FR = calc_FFI(exp_paths)

% exp_paths = {'H:\Tlx3project\T9\T9_driftgrating_trains_160125_123706',...
%             'H:\Tlx3project\T14\T14_driftgratingtrains_160329_211516',...
%             'H:\Tlx3project\T15\3-31-16\T15_driftgratingtrains_160324_214136',...
%             'H:\Tlx3project\T19\5-2-16\T19_trains_160502_222135',...
%             'H:\Tlx3project\T18\5-3-16\T18_trains_160503_192729'};

% make rasters for multiple experiments
for i = 1:length(exp_paths)
    exp_path = exp_paths{i};
    load(sprintf('%s\\su_spikedata.mat',exp_path))
    for n=1:length(spiketimes)          % for each unit in that experiment
        spikes = spiketimes{n};
        spike_raster{i}(n,:,:) = make_raster(spikes,1000,2.5);
    end
end

% for each experiment, calculate FRs in first and last pulse of long train
count = 1;
for i = 1:length(exp_paths)
    load(sprintf('%s\\light_params.mat',exp_paths{i}))
    light_conds = unique(all_light);
    light_trials = find(all_light ==light_conds(2));         % TEMP
    time = 0:1/1000:2.5;
    lightstart = find(time==round(av_light_start));        % SAMPLE at which light started
    win = 20;                  % number of samples to use for FR windows
    pre_win = [lightstart-(win):lightstart-1];
    onset_win = [lightstart:lightstart+win-1];
    post_win = [lightstart+win:lightstart+2*win-1];
    for n =1:size(spike_raster{i},1)
        win_FR(count,1) = sum(sum(spike_raster{i}(n,light_trials,pre_win)))/(length(light_trials)*(win/1000));
        win_FR(count,2) = sum(sum(spike_raster{i}(n,light_trials,onset_win)))/(length(light_trials)*(win/1000));
        win_FR(count,3) = sum(sum(spike_raster{i}(n,light_trials,post_win)))/(length(light_trials)*(win/1000));
        count = count+1;
    end
end
        
end