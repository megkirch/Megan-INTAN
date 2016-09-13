function kilo_spike_analysis(exp_path,unittype,exp_type)

% modified from phy_spike_analysis_v2 to make more efficient and to 
    % separate out graphing - MAK 3/14/16

% Takes in data from Intan, analyzer file, and phy .kwik file to look at
% spike times and spike frequency from isolated single and multi unit
% clusters
% unittype = 'su' or 'mu'    (single or multi units)
% pulse_dur = duration of light (in sec)

% 6/7/16 - removed prelight and postlight parameters, replaced with
% av_light_start from get_lightstim.m
% 7/17 - changed to get_lightstim_v2.m

cd(exp_path)
if exist(sprintf('%s/data.mat',exp_path),'file')
    load(sprintf('%s/data.mat',exp_path))      % data from intanphy2matlab.m
else
    intanphy2matlab_v2(exp_path);      % data from intanphy2matlab.m
    load(sprintf('%s/data.mat',exp_path))      % data from intanphy2matlab.m
end


% trials - num_trials x 2 matrix of trial start and end times (in seconds)
% field_trials - trial num x 2 matrix of trial start and end samples (in 
    % 1000 Hz sampling rate) - starts from 1!
% time_index - vector of time stamps of each sample (in seconds), with 1000
    % samples per second - starts from 0!
% amp_sr - original amplifier sampling rate
% epoc - vector of 1s and 0s indicating whether each sample was part of a 
    % trial (1) or not (0). Still in original sampling rate (amp_sr)
% mmvt -  vector of 1s and 0s indicating whether mouse was running (1) or  
    % not (0) at time of sample. Still in original sampling rate (amp_sr)
% encd - vector of analog output from movement encoder. Still in original 
    % sampling rate (amp_sr)
% LED - is 0 if not an optogenetics experiment, but otherwise, vectors of
    % whether LED was on (1) or off (0)
% re - vector for the timestamps (rising phase of photodiode signal). you 
    % need to remove the first 2 and the last 2 timestamps. Timestamps 
    % signify the end of a four second black to gray period.
% Photo - vector of analog output from photodiode. In 1000Hz sampling rate

    
% load kilosort data
spike_times = readNPY('spike_times.npy');
spike_times = floor(spike_times/(amp_sr/1000));   % change sample #s to account for downsampling
clusters = readNPY('spike_clusters.npy');

% load analyzer
s = dir; 
for i=1:length(s)
    if strfind(s(i).name,'.analyzer') 
        analyze_file = s(i).name;
    end
end
load(sprintf('%s/%s',exp_path,analyze_file),'-mat')     % load analyzer file with stimulus info
prestim = Analyzer.P.param{1}{3};
poststim = Analyzer.P.param{2}{3};
if poststim == 0
    extra = .25*1000;   % add extra samples for pseudo-poststim period
else
    extra = 0;
end
stimtime = Analyzer.P.param{3}{3};
totaltime = prestim+poststim+stimtime;
onset = .2;     % onset time (in seconds) - this will be subtracted from the beginning of the stim period to determine evoked firing rate


% optogenetics?
if sum(LED) > 0
    light_exp = 1;          % if using optogenetics
else
    light_exp = 0;
end

% downsample
LN              = length(epoc);
div             = amp_sr/1000;
zx              = 1:div:LN;
izx             = floor(zx);
epoch = epoc(izx);
% move = mmvt(izx);
% encode = encd(izx);
if light_exp
    light = LED(izx);
end

% extract trial types
num_trials = size(trials,1);
for v = 1:length(Analyzer.loops.conds{1}.val)   % for each variable
    conds{v} = cellfun(@(x) x.val{v}, Analyzer.loops.conds, 'UniformOutput',false);   % all possible conditions
    if strcmp(Analyzer.loops.conds{1}.symbol{v},'ori')
        blanks = cellfun(@(x) find(isempty(x)), conds{v}, 'UniformOutput', false);
        orivar = v;
    end
end
num_conds = length(conds{1});    % number of different conditions 
oris = unique(cell2mat(conds{orivar}));   % actual orientations (e.g. 0, 30, 60, etc.)
num_oris = length(oris);
trial_type = zeros(1,num_trials);
for c = 1:num_conds
    ori_trials{c} = cellfun(@(x) x.trialno, Analyzer.loops.conds{c}.repeats);     % store trials from different conditions in separate columns
    trial_type(ori_trials{c}) = c;        % use when using ALL trials
%     trials2use = ori_trials(find(ori_trials(:,c)>num_trials),c);      % uncomment if you want to use only second half of experiment
%     trial_type(trials2use-300) = c;
end

% in event that analyzer files had to be manually combined (i.e. something
% went weird in experiment)
if exist(sprintf('%s/trial_type.mat',exp_path),'file')
    clear trial_type
    load(sprintf('%s/trial_type.mat',exp_path))
end

blank_conds = [];
for bl = 1:length(blanks)
    if blanks{bl}
        blank_conds = [blank_conds bl];
    end
end
% isolate blank trials
if ~isempty(blank_conds)    % if there were blank trials
    blank_trials = find(ismember(trial_type,blank_conds));
    visual_trials = find(~ismember(trial_type,blank_conds));
else
    visual_trials = 1:length(trial_type);       
end
% find light trials
if light_exp
    lightvar = find(strcmp(Analyzer.loops.conds{1}.symbol,'light_bit'));
    diff_lightconds = unique(cellfun(@(x) x.val{lightvar}, Analyzer.loops.conds));
    num_lightconds = length(find(diff_lightconds));
    all_light_trials = zeros(1,num_trials);
    for cc = 1:num_lightconds        % for each condition WITH light stimulation
        all_light_conds(cc,:) = find(cellfun(@(x) x.val{lightvar}, Analyzer.loops.conds)==diff_lightconds(cc+1)); 
        all_light_trials(find(ismember(trial_type,all_light_conds(cc,:))))=cc;
    end
    nolight_conds = find(cellfun(@(x) x.val{lightvar}, Analyzer.loops.conds)==0); 
    light_conds = sort(all_light_conds(:));
    light_trials = find(ismember(trial_type,light_conds));
    nolight_trials = find(~ismember(trial_type,light_conds));
    
    % The following also works for identifying light trials IF the LED
    % signal was digital
%     for t = 1:num_trials  % for each trial
%         if sum(light(field_trials(t,1):field_trials(t,2)))  % if LED turned on during trial
%             light_on(t) = 1;
%         else
%             light_on(t) = 0;
%         end
%     end
%     light_trials = find(light_on);

    % get info about light parameters
    [all_light_hz, pulse_dur, lighttime, av_light_start] = get_lightstim_v2(exp_path,exp_type);
    lightconds = unique(all_light_hz);
    [type,idx] = sort(all_light_hz);    % sort trials by light conditions
end

% find trials when mouse was running
if exist('encdA','var')
    move_trials = Intan_digital_movement(field_trials,encdA,encdB,0);
elseif exist('mouse','var')    % using optical mouse for old recordings - might want to change to encoder?
    move_trials = zeros(1,size(field_trials,1));
    moveVec = move_trials;
    for i = 1:size(field_trials,1)
        moveVec(i) = sum(mouse(field_trials(i,1)*20:field_trials(i,2)*20));         % *20 because mouse is still in original amp_sr
    end
    move_trials(find(moveVec)) = 1;
end

prestim_samps = prestim*1000;
poststim_samps = poststim*1000;
stimtime_samps = stimtime*1000;
lighttime_samps = lighttime*1000;
spike_times = spike_times + 1; % has to be +1 because spike_times starts at 0, but the min possible field_trials value could be 1    

fid = fopen('cluster_groups.csv');
% read column headers
C_text = textscan(fid, '%s', 1, 'delimiter', ',');
% read group info
grp_dat = textscan(fid, '%f %s', 'delimiter', ',');
fclose(fid);
units = grp_dat{1};             % cluster numbers
cluster_group = grp_dat{2};     % 'good', 'noise', etc.
num_units = length(units); 

% determine which group each cluster belongs to (good, MUA, etc.)
single_units = units(find(strcmp(cluster_group,'good')));
multi_units = units(find(strcmp(cluster_group,'mu')));          % CHECK

%% extract spike times for units of interest (su or mu)
if strcmp(unittype,'su')
    good_units = single_units;
elseif strcmp(unittype,'mu')
    good_units = multi_units;
else
    error('Didnt specify singleunits (su) or multiunits (mu)')
end

% if units to analyze are indicated in text file
if exist(sprintf('%s/good_units.txt',exp_path),'file')
    good_units = load(sprintf('%s/good_units.txt',exp_path));
end

% get each cluster's spikes by trial, and count
spiketimes = cell(1,length(good_units));
spikes_prestim = nan(num_trials,length(good_units));
spikes_ev = spikes_prestim;
spikes_onset = spikes_prestim;
spikes_all = spikes_prestim;
for n = 1:length(good_units)
    unit = find(units == good_units(n));
    unit_times = spike_times(find(clusters==good_units(n))); % timestamps of unit's spikes
    spiketimes{n} = spikes_by_trial(unit_times,field_trials,extra,time_index);
    
    % count spikes during prestim, stim, poststim, and onset periods
    for t = 1:length(spiketimes{n})
        spikes_all(t,n) = length(spiketimes{n}{t});
        spikes_prestim(t,n) = length(find(spiketimes{n}{t}<prestim));
        spikes_ev(t,n) = length(find((spiketimes{n}{t}>=av_light_start)&(spiketimes{n}{t}<=av_light_start+lighttime)));
        spikes_onset(t,n) = length(find((spiketimes{n}{t}>=prestim)&(spiketimes{n}{t}<prestim+onset)));
        spikes_prestim_onsettime(t,n) = length(find((spiketimes{n}{t}>=prestim-onset)&(spiketimes{n}{t}>prestim-onset)<prestim));       % for kruskal-wallis test of visual modulation
    end
    
end

%% Calculate firing rates
for n = 1:length(good_units)
    for i = 1:length(lightconds)
        [spikerate_prestim(n,i), spikerateSE_prestim(n,i)] = calc_firing_rates(spikes_prestim(:,n),intersect(find(all_light_hz == lightconds(i)),visual_trials),prestim); 
        [spikerate_ev(n,i), spikerateSE_ev(n,i)] = calc_firing_rates(spikes_ev(:,n),intersect(find(all_light_hz == lightconds(i)),visual_trials),lighttime);
        [spikerate_onset(n,i), spikerateSE_onset(n,i)] = calc_firing_rates(spikes_onset(:,n),intersect(find(all_light_hz == lightconds(i)),visual_trials),onset);

    %     % calculate firing rates by condition
    %     for c = 1:num_conds
    %         which_trials = find(trial_type == c);
    %         [spikerate_ev_bycond(c,n), spikerateSE_ev_bycond(c,n)] = calc_firing_rates(spikes_ev(:,n),which_trials,lighttime);
    %     end
        
        if ~isempty(blank_conds)    % if there were blank trials
            % calculate firing rates for blanks
            [spikerate_blank(n,i) spikerateSE_blank(n,i)] = calc_firing_rates(spikes_ev(:,n),intersect(find(all_light_hz == lightconds(i)),blank_trials),lighttime);
        else
            spikerate_blank = [];
            spikerateSE_blank = [];
        end
        
        if sum(move_trials)         % if mouse ran at all
            [spikerate_run(n,i) spikerateSE_run(n,i)] = calc_firing_rates(spikes_ev(:,n),intersect(find(all_light_hz == lightconds(i)),find(move_trials)),lighttime);
            [spikerate_stat(n,i) spikerateSE_stat(n,i)] = calc_firing_rates(spikes_ev(:,n),intersect(find(all_light_hz == lightconds(i)),find(~move_trials)),lighttime);
        else
            spikerate_run(n,i) = 0;
            spikerateSE_run(n,i) = 0;
            [spikerate_stat(n,i) spikerateSE_stat(n,i)] = calc_firing_rates(spikes_ev(:,n),intersect(find(all_light_hz == lightconds(i)),find(~move_trials)),lighttime);
        end
    end
end

%% save data
save(sprintf('%s_spikedata.mat',unittype), 'good_units','spiketimes','spikes_all', ...
    'spikes_prestim', 'spikes_ev', 'spikes_onset', 'spikes_prestim_onsettime',...
    'spikerate_prestim', 'spikerateSE_prestim', 'spikerate_ev',...
    'spikerateSE_ev', 'spikerate_onset', 'spikerateSE_onset', ...
    'spikerate_blank', 'spikerateSE_blank','spikerate_run','spikerateSE_run',...
    'spikerate_stat','spikerateSE_stat')
    %     'spikerate_ev_bycond', 'spikerateSE_ev_bycond', ...

return



function [spikerate, spikerate_SE] = calc_firing_rates(spikes,...
    which_trials,period)

% spikes = column vector of number of spikes by trial for unit of interest
% which_trials = vector of which trials you want to include
% period = length of period (in seconds) that number of spikes was counted 
    % from (e.g. evoked: 1.8s)

spikerate = mean(spikes(which_trials)/period);
spikerate_SE = std(spikes(which_trials)/period)/sqrt(length(which_trials));

return
