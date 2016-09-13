function [all_light, pulse_dur, light_dur, av_light_start] = get_lightstim_v2(exp_path,exp_type)

% get_lightstim_v2.m - use LED input (analog or digital) to determine the parameters of the
% optogenetic stimulation
% updated by MAK 7/14/16
% Inputs:
% exp_path - e.g. 'H:\Tlx3project\T19\5-2-16\T19_diffintensities_160502_210213'
% exp_type - 'ramp', 'trains', 'intensities' or 'size'
% Outputs:
% all_light - 1xnum_trials vector of values 1:num_conds (0 = no light, 1+ =
% light conditions)
% pulse_dur - duration of light pulse in sec (single pulse, not total
% duration of light!)
% light_dur = total duration of light (from beginning of first pulse to end
% of last pulse)
% av_light_start - average time from visual stimulus onset that light
% turned on

cd(exp_path)
if exist(sprintf('%s/data.mat',exp_path),'file')
    load(sprintf('%s/data.mat',exp_path))      % data from intanphy2matlab.m
else
    [~,field_trials,time_index,amp_sr,...
    epoc,~,~,LED,~,~] = intanphy2matlab_v2(exp_path);      % data from intanphy2matlab.m
end

% field_trials - trial num x 2 matrix of trial start and end samples (in 
    % 1000 Hz sampling rate) - starts from 1!
% time_index - vector of time stamps of each sample (in seconds), with 1000
    % samples per second - starts from 0!
% amp_sr - original amplifier sampling rate
% encd - vector of analog output from movement encoder. Still in original 
    % sampling rate (amp_sr)
% LED - is 0 if not an optogenetics 
    % experiment, but otherwise, vectors of whether LED was on (1) or off (0)



LN              = length(epoc);
div             = amp_sr/1000;
zx              = 1:div:LN;
izx             = floor(zx);
light = LED(izx);       % put into 1000hz sampling rate

%% first, find light trials
for t = 1:size(field_trials,1)      % for each trial
    samps_per_t = max(diff(field_trials,[],2));                     % do this in case each trial doens't have exact same number of samples
    light_out(t,:) = light(field_trials(t,1):field_trials(t,1)+samps_per_t);      % get actual values of LED input to Intan
end

[~,ex_trial] = min(mean(light_out,2));        % get example of a trial in which there was definitely no light
light_thresh = 5*std(light_out(ex_trial,:))+max(light_out(ex_trial,:));       % considered light trial if max value  exceeds 5 standard deviations of ex_trial above the max value of ex_trial
light_trials = find(max(light_out,[],2)>light_thresh)';       % which trials the light turned on in
light_on = zeros(1,size(field_trials,1));
light_on(light_trials) = 1;                         % matrix of all trials: 0 if no light, 1 if light

% check against analyzer file
s = dir; 
for i=1:length(s)
    if strfind(s(i).name,'.analyzer') 
        analyze_file = s(i).name;
    end
end
load(sprintf('%s/%s',exp_path,analyze_file),'-mat')     % load analyzer file with stimulus info

% get light conds
lightcond = find(strcmp(Analyzer.loops.conds{1}.symbol,'light_bit'));
conds = cellfun(@(x) x.val{lightcond}, Analyzer.loops.conds);
diff_conds = unique(conds);
num_conds = length(diff_conds);
% get trial types
trial_type = zeros(1,size(field_trials,1));
for c = 1:length(conds)
    ori_trials{:,c} = cellfun(@(x) x.trialno, Analyzer.loops.conds{c}.repeats);     % store trials from different conditions in separate columns
    trial_type(ori_trials{:,c}) = c;        
end
trial_type = trial_type(1:length(light_on));
light_trials_check = find(ismember(trial_type,find(conds)));
min_length = min(length(light_trials),length(light_trials_check));      % in case number of light trials are different between analyzer and LED input to Intan (suggests there's probably a problem...)
if sum(light_trials(1:min_length)~=light_trials_check(1:min_length))    % if numbers of light trials don't agree
    err_count = 1;
    fprintf('Analyzer file and LED inputs do not match\n')
    off_trial = find(light_trials(1:min_length)~=light_trials_check(1:min_length),1,'first');       % check first conflicting trial
    h = figure; plot(light_out(off_trial,:));
    check = input('Was this a light trial? Y or N: ','s');
    if strcmp(check,'N')
        light_trials(off_trial) = [];   
    end
    close(h)
    if length(light_trials) > length(light_trials_check)        % if more light trials were detected by light input to Intan than in analyzer file
        offsetA = 1;
        offsetB = 0;
    elseif length(light_trials) < length(light_trials_check)        % vice versa
        offsetA = 0;
        offsetB = 1;
    else                    % correct number of trials, just disagreement about which ones were light trials
        offsetA = 0;
        offsetB = 0;
    end
    if  sum(offsetA+offsetB) && sum(light_trials(off_trial+offsetA:end)~=light_trials_check(off_trial+offsetB:end))     % check the next conflicting trial
%         fprintf('Analyzer file and LED inputs are not same length\n')
        off_trial2 = find(light_trials(off_trial+1:end)~=light_trials_check(off_trial:end),1,'first');
        h = figure; plot(light(field_trials(off_trial+off_trial2,1):field_trials(off_trial+off_trial2,2)))
        check = input('Was this a light trial? Y or N: ','s');
        if strcmp(check,'N')
            light_trials(off_trial+off_trial2) = [];
        end
        close(h)
    else
        fprintf('LED input conflicts with analyzer file. Going with LED input\n')
    end
else
    err_count = 0;      % no conlict between LED input and analyzer file
end



%% determine light parameters given the type of light pulse

if strcmp(exp_type,'ramp')
    for t = 1:length(light_trials)
        tt = light_trials(t); 
        light_end(t) = find(light_out(tt,:)>light_thresh,1,'last');         % find last sample when light was on
        light_dur = 1;      % assumes ramp lasts 1sec
        light_start(t) = light_end(t)-light_dur*1000+1;                     % go backwards to find the start
    end
    light_hz = 1;
    pulse_dur = round(mean(light_dur));
    all_light = zeros(1,size(field_trials,1));       % includes trials w/o light stim (as 0)
    all_light(light_trials) = light_hz; 
    
elseif strcmp(exp_type,'trains')
    for t = 1:length(light_trials)
        tt = light_trials(t); 
        light_start(t) = find(light_out(tt,:)>light_thresh,1,'first');
        light_pulsedur(t) = find(light_out(tt,light_start(t):end)<=light_thresh,1,'first')-1;       % finds duration of first light pulse in the trial
        light_end(t) = find(light_out(tt,:)>light_thresh,1,'last');
        light_dur(t) = length(light_start(t):light_end(t))/1000;
        light_sum(t) = sum(light_out(tt,:)>light_thresh);     % sum of light pulses - used to figure out frequency
    end
    pulse_dur = round(mean(light_pulsedur));       % in sec
    light_dur = round(max(light_dur));              % in case of trials with only single pulse, still consider light duration of trials with multiple pulses
    light_hz = round(floor(light_sum./(pulse_dur)));     % get the light frequencies for each trial (use floor in case ex. sum of 40Hz trial is 81 instead of 80)
    pulse_dur = pulse_dur/1000;         % change from ms to sec
    all_light = zeros(1,size(field_trials,1));       % includes trials w/o light stim (as 0)
    all_light(light_trials) = light_hz; 
    
else                % simple step function
    for t = 1:length(light_trials)
        tt = light_trials(t); 
        light_start(t) = find(light_out(tt,:)>light_thresh,1,'first');
        light_end(t) = find(light_out(tt,:)>light_thresh,1,'last');
        light_dur(t) = length(light_start(t):light_end(t))/1000;
        light_inten(t) = roundn(mean(light_out(tt,light_start(t):light_end(t))),-1);          % round intensity to nearest hundreth
%         light_inten(t) = roundn(mean(light_out(tt,light_start:light_end)),-1);          % round intensity to nearest hundreth
    end
    intensities = unique(light_inten);
    light_hz = 1;
    pulse_dur = round(mean(light_dur));
    light_dur = pulse_dur;
    all_light = zeros(1,size(field_trials,1));       % includes trials w/o light stim (as 0)
    all_light(light_trials) = light_hz; 

%     % TEMP b/c TiC2
%     diff_durs = unique(light_inten);
%     for i=1:length(diff_durs)
%         which_dur = ismember(light_inten,diff_durs(i));            % find trials that were at a particular light intensity
%         all_light(light_trials(which_dur)) = diff_durs(i);
%     end
    
    if length(intensities) > 1              % if this experiment involved light steps at different intensities - use LED input to determine intensities
        small_diffs = find(diff(intensities) < .5);         % to account for minor variability amongst intensities
        if ~isempty(small_diffs)
    %         find((light_inten == intensities(small_diffs)) | (light_inten == intensities(small_diffs+1)))
            orig_intensities = intensities;
            intensities(small_diffs) = [];
            count = 1;
            for i = 1:length(intensities)
                diff_intensities{i} = orig_intensities(count);
                count = count+1;
                if i == small_diffs
                    diff_intensities{i} = orig_intensities(small_diffs:small_diffs+1);
                    count = count+1;
                end
            end
        else
            for i = 1:length(intensities)
                diff_intensities{i} = intensities(i);
            end
        end

        all_light = zeros(1,size(field_trials,1));       % includes trials w/o light stim (as 0)
        for i = 1:length(intensities)
            which_inten = ismember(light_inten,diff_intensities{i});            % find trials that were at a particular light intensity
            all_light(light_trials(which_inten)) = intensities(i);
        end
    else                                % if single intensity experiment
        all_light(light_trials) = 1;     
    end
end

%% Get final parameters to save 
av_light_start = time_index(round(mean(light_start)));      % average time at which light started (1 less than sample # b/c starting with 0!)
% pulse_dur and all_light already made

% save parameters
save('light_params.mat','all_light', 'pulse_dur', 'av_light_start');

return
