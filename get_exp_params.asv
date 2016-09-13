function [prestim,poststim,stimtime,trial_type,IVs] = get_exp_params(exp_path,exp_type)

% from analyzer file, determine the different independent variables (IVs) 
% in the experiment and create a matrix 'trial_type' of num_trials x 
% num_IVs in which values reflect the level of each variable
% also output the prestim, poststim, and stim times
% Inputs:
%   exp_path - e.g. 'H:\Tlx3project\T19\5-2-16\T19_diffintensities_160502_210213'
%   exp_type - 'ramp', 'trains', 'intensities' or 'size'
% created 8/11/16 by MAK

cd(exp_path)

% load analyzer
s = dir; 
for i=1:length(s)
    if strfind(s(i).name,'.analyzer') 
        analyze_file = s(i).name;
    end
end
load(sprintf('%s/%s',exp_path,analyze_file),'-mat')     % load analyzer file with stimulus info
prestim = Analyzer.P.param{1}{3};       % in seconds
poststim = Analyzer.P.param{2}{3};
if poststim == 0
    extra = .25*1000;   % add extra samples for pseudo-poststim period
else
    extra = 0;
end
stimtime = Analyzer.P.param{3}{3};
totaltime = prestim+poststim+stimtime;
onset = .2;     % onset time (in seconds) - this will be subtracted from the beginning of the stim period to determine evoked firing rate
num_trials = Analyzer.M.NumTrials;          % here, num_trials is defined from the analyzer file rather than from the actual experiment

% convert to ms
prestim_samps = prestim*1000;
poststim_samps = poststim*1000;
stimtime_samps = stimtime*1000;

% extract trial types
% num_trials = size(trials,1);
IVs = Analyzer.loops.conds{1}.symbol;       % independent variables
trial_type = zeros(num_trials,length(IVs));     % matrix of trials by IVs
for v = 1:length(IVs)   % for each variable
    conds{v} = cellfun(@(x) x.val{v}, Analyzer.loops.conds, 'UniformOutput',false);   % all possible conditions
    levs{v} = unique(cell2mat(conds{v}));
    num_levs(v) = length(levs{v});   % number of levels of IV
    if sum(cellfun('isempty',conds{v}));        % if '[]' was a condition (e.g. for blank trials)
        num_levs(v) = num_levs(v)+1;
        levs{v} = [levs{v} 999];        % blank trials indicated by 999 instead of empty cell (tried with NaN but was too annoying)
    end
    for i = 1:num_levs(v)
        if levs{v}(i)==999
            which_conds = find(cellfun('isempty',conds{v}));    % indices of which conds were blank conds
        else
            which_conds = find([conds{v}{:}] == levs{v}(i));    % indices of which non-blank conds have the same level of variable v
        end
        cond_trials = [];
        for ii = 1:length(which_conds)
            cond_trials = [cond_trials cellfun(@(x) x.trialno, Analyzer.loops.conds{which_conds(ii)}.repeats)];     % store trials from different conditions in separate columns
        end
        cond_trials = sort(cond_trials);
        trial_type(cond_trials,v) = levs{v}(i);        % use when using ALL trials
    end
end

% import data from intan2matlab.m and make sure trial_type reflects the
% number of trials that actually happened (may be less than num_trials
% because analyzer file is created before the experiment actually starts)
load(sprintf('%s/data.mat',exp_path))
real_num_trials = size(field_trials,1);
trial_type = trial_type(1:real_num_trials,:);

% check analyzer-determined light trials against the LED input to Intan
if exist('LED','var')             % if it was an optogenetics experiment
    lightvar = find(strcmp(IVs,'light_bit'));
    [all_light,~,~,~] = get_lightstim_v2(exp_path,exp_type);
    trial_type(:,lightvar) = all_light';
end

% add column to indicate whether trials were visual (1) or blank (0) 
orivar = find(strcmp(IVs,'ori'));
blank_trials = trial_type(:,orivar)<=360;
trial_type = [blank_trials trial_type];         % add column to FRONT of matrix
IVs = {'visual' IVs{:}};                   % keep track of variable names

% determine run vs. stationary trials
if exist('encdA','var')
    move_trials = Intan_digital_movement(field_trials,encdA,encdB,0)';
elseif exist('mouse','var')    % using optical mouse for old recordings - might want to change to encoder?
    move_trials = zeros(1,size(field_trials,1));
    moveVec = move_trials;
    for i = 1:size(field_trials,1)
        moveVec(i) = sum(mouse(field_trials(i,1)*20:field_trials(i,2)*20));         % *20 because mouse is still in original amp_sr
    end
    move_trials(find(moveVec)) = 1;
else
    move_trials = zeros(1,size(field_trials,1));
end
trial_type = [trial_type move_trials'];          % add column to trial_type indicating stationary (0) or run (1) trials
IVs{length(IVs)+1} = 'running';

end