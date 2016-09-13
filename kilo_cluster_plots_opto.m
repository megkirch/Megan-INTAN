function kilo_cluster_plots_opto(exp_path,unittype,normalizing,exp_type)

% modified from phy_spike_analysis_v2.m - takes in data from
% phy_spike_analysis_v3.m
% MAK - 3/14/16
% assumes this is an optogenetics experiment
% unittype = 'su' or 'mu'    (single or multi unit)
% normalizing =  'blank' or 'prestim' to designate how to correct firing
    % rates
% exp_type = none, or 'ramp', 'trains', 'intensities' or 'size'. default is none.

% 6/7/16 - removed prelight and postlight parameters, replaced with
% av_light_start from get_lightstim.m
% 6/13/16 - added exp_type input and plots for running
% 7/17/16 - changed get_lightstim to get_lightstim_v2

color_mat = [0 0 0; 0 0.2 .7; 0 0.3 .5; 0 0.4 .4; 0 0 1]; % for graphing purposes (first is black, last is blue)

% set up experiment directory
exp_dir = 'H:\Tlx3project\ChR2';     
% animal_name = regexp(exp_path,'T\d+','match');
% animal_name = animal_name{1};
% run_name = regexp(exp_path,'run\d','match');
% if isempty(run_name)
%     exp_name = animal_name;
% else
%     exp_name = strcat(animal_name,run_name{1});
% end
exp_name = 'LP';              % TEMP

% load data from phy_spike_analysis_v3
cd(exp_path)
% if exist(sprintf('%s/%s_spikedata.mat',exp_path,unittype),'file')
%     load(sprintf('%s/%s_spikedata.mat',exp_path,unittype))      
% else
    kilo_spike_analysis(exp_path,unittype,exp_type);
    load(sprintf('%s/%s_spikedata.mat',exp_path,unittype)) 
% end

% get amplifier sampling rate (from intanphy2matlab.m)
load(sprintf('%s/data.mat',exp_path),'amp_sr')      % data from intanphy2matlab.m
div = amp_sr/1000;  % samples per millisecond

% make directory to store figures (if one doesn't already exist)
fig_dir = sprintf('%s/Figures',exp_path);
if ~exist(fig_dir)
    mkdir(fig_dir);
end
all_fig_dir = sprintf('%s/All_clusters',fig_dir);
if ~exist(all_fig_dir)
    mkdir(all_fig_dir);
% else
%     rmdir(all_fig_dir,'s')          % clear figures folder if it already exists
%     mkdir(all_fig_dir);
%     delete(sprintf('%s/Figures/*.png',all_fig_dir))

end

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
    poststim = .25;
end
stimtime = Analyzer.P.param{3}{3};
totaltime = prestim+poststim+stimtime;
onset = .2;     % onset time (in seconds) - this will be subtracted from the beginning of the stim period to determine evoked firing rate

% Determine which probe was used
tmp = dir(exp_path);
possible_probes = {'A2x16_2.prb','A1x32.prb','NeuroNexus_8FCS.prb','A1x32_A48.prb'};
for i = 1:length(tmp)
    isprobename = strcmp(tmp(i).name,possible_probes);
    if sum(isprobename)
        probetype = find(isprobename);
    end
end

% get necessary trial types
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
for c = 1:num_conds
    ori_trials{c} = cellfun(@(x) x.trialno, Analyzer.loops.conds{c}.repeats);     % store trials from different conditions in separate columns
    trial_type(ori_trials{c}) = c;        % use when using ALL trials
%     trials2use = ori_trials(find(ori_trials(:,c)>num_trials),c);      % uncomment if you want to use only second half of experiment
%     trial_type(trials2use-300) = c;
end
blank_conds = [];
for bl = 1:length(blanks)
    if blanks{bl}
        blank_conds = [blank_conds bl];
    end
end

% in event that analyzer files had to be manually combined (i.e. something
% went weird in experiment)
if exist(sprintf('%s\\trial_type.mat',exp_path),'file')
    clear trial_type
    load(sprintf('%s\\trial_type.mat',exp_path))
end

if ~isempty(blank_conds)    % if there were blank trials
    blank_trials = ismember(trial_type,blank_conds);
    visual_trials = ~ismember(trial_type,blank_conds);
else
    visual_trials = ones(1,length(trial_type));       
end
save('blank_trials.mat','blank_trials')

% get info about light parameters
% if exist(sprintf('%s\\light_params.mat',exp_path),'file')
%     load(sprintf('%s\\light_params.mat',exp_path))   
% else
    [all_light_hz, pulse_dur, lighttime, av_light_start] = get_lightstim_v2(exp_path,exp_type);
% end
lightconds = unique(all_light_hz);
[type,idx] = sort(all_light_hz);    % sort trials by light conditions
for n = 1:length(lightconds)
    trials_per_cond(n) = length(find(all_light_hz==lightconds(n)));
end
% make appropriate figure legends
if strcmp(exp_type,'trains')
    for i = 2:length(lightconds)        % for graphing purposes
        legend_labels{i-1} = sprintf('%dHz',lightconds(i));
    end
elseif strcmp(exp_type,'intensities')
    legend_labels = {'Low light','Medium light','High light'};
    full_legend_labels = {'Light OFF','Low light','Medium light','High light'};
else
    legend_labels = {'Light ON'};
end

%% Run CSD (if you haven't already)
% if exist(sprintf('%s\\layers2.mat',exp_path),'file')
%     layers = importdata(sprintf('%s\\layers.mat',exp_path));
% else
%     IntanLFP_CSD(exp_path)
% end

%% find running trials
load(sprintf('%s/data.mat',exp_path),'field_trials','encdA','encdB','mouse');
if exist('encdA','var')
    move_trials = Intan_digital_movement(field_trials,encdA,encdB,0);
elseif exist('mouse','var')    % using optical mouse for old recordings - might want to change to encoder?
    move_trials = zeros(1,size(field_trials,1));
    moveVec = move_trials;
    for i = 1:size(field_trials,1)
        moveVec(i) = sum(mouse(field_trials(i,1)*20:field_trials(i,2)*20));     % *20 because mouse is still in original amp_sr
    end
    move_trials(find(moveVec)) = 1;
end

%% First, make array of important plots for all units
% load kilosort data
spike_times = readNPY('spike_times.npy');
spike_times_ds = floor(spike_times/(amp_sr/1000));   % change sample #s to account for downsampling
clusters = readNPY('spike_clusters.npy');

% first, extract waveforms (do it only once because it's time consuming)
for n = 1:length(good_units)
    load_waveforms(n) = exist(sprintf('Cluster_%s_waveforms.mat',num2str(good_units(n))),'file');
end

if length(find(load_waveforms)) < length(good_units)     % if any waveforms are missing, reload them all
    for n = 1:length(good_units)
        fprintf(sprintf('Extracting waveforms for cluster %s \n',num2str(good_units(n))))
        [waveforms_microV, max_ch] = readWaveformsFromDat(sprintf('%s/amplifier.dat',exp_path),32,spike_times(find(clusters==good_units(n))),[-16 16],[],4);
        save(sprintf('Cluster_%s_waveforms.mat',num2str(good_units(n))),'waveforms_microV','max_ch');
    end
end

for n = 1:length(good_units)
  
    fig_title = ['Cluster: ' num2str(good_units(n))];
    clust_fig = figure('name', fig_title); %title(sprintf('Cluster %s',clusters2analyze(n))); hold on
    unit = good_units(n);
    unit_times = spiketimes{n}; % timestamps of unit's spikes
    num_trials = length(unit_times);
    
    % 1) Raster plot
    subplot(241)
    ylim([0 num_trials])
    xlim([-prestim stimtime+poststim])
    % first, draw light pulses (Raster will go on top)
    x1 = av_light_start-prestim;       % when the light starts
    for c = 1:length(lightconds)
        start_patch = find(type==lightconds(c),1,'first')-1;
        end_patch = find(type==lightconds(c),1,'last')-1;
        for p = 1:(all_light_hz(idx(end_patch)))
            space = 1/(all_light_hz(idx(end_patch)));
            patch([x1+((p-1)*space) x1+((p-1)*space) x1+pulse_dur/1000+((p-1)*space) x1+pulse_dur/1000+((p-1)*space) x1+((p-1)*space)],[start_patch end_patch end_patch start_patch start_patch], [0.9 0.9 0.9], 'LineStyle', 'none')
        end
    end
    for t = 1:num_trials        % for each trial
        tt = idx(t);
        for i=1:length(unit_times{tt})
            if type(t)  % if it's a light condition, make it blue
                line([unit_times{tt}(i)-prestim unit_times{tt}(i)-prestim], [t-1 t]','Color','b')                  
            else
                line([unit_times{tt}(i)-prestim unit_times{tt}(i)-prestim], [t-1 t]','Color','k')
            end
        end
    end
    xlabel('Time (sec)')
    ylabel('Trial #')
    line([0 0], [0 num_trials]','Color','r','LineStyle','--')
%     line([onset onset], [0 num_trials]','Color','r','LineStyle','--')
    line([stimtime stimtime], [0 num_trials]','Color','r','LineStyle','--')
    for ii = 1:length(unique(type))
        line([-prestim stimtime+poststim], [trials_per_cond(ii)*(ii-1) trials_per_cond(ii)*(ii-1)], 'Color','b','LineStyle','--')
    end
    title('Raster plot')

    % 2) PSTH
    binsize = .025;      % in seconds   
    edges = [0:binsize:totaltime-binsize];
    psth = make_psth(binsize,edges,visual_trials(1:num_trials),spiketimes{n},all_light_hz);
    if ~isempty(blank_conds)
        psth_blanks = make_psth(binsize,edges,blank_trials(1:num_trials),spiketimes{n},all_light_hz);
    else
        psth_blanks = [];
    end

    subplot(242)
%     edges_stim = [-prestim:step:(stimtime+poststim-step)]'; % x signifies the timepoint of the START of the bin
    edges_stim = [-prestim:binsize:(stimtime+poststim-binsize)]'; % x signifies the timepoint of the START of the bin
    for c = 1:size(psth,2)
        plot(edges_stim,psth(:,c),'color',color_mat(c,:),'linewidth',2)
        hold on
    end
%     plot(edges_stim,psth(:,1),'k','linewidth',2)
%     hold on
%     plot(edges_stim,psth(:,2),'b','linewidth',2)

%     xlim([-prestim stimtime+poststim-binsize])  % because points mark the START of the bin
    xlim([-prestim stimtime+poststim-binsize])  % because points mark the START of the bin
    set(gca,'XMinorTick','on')
    yax = get(gca,'YLim');
    line([0 0], [0 yax(2)]','Color','r','LineStyle','--')
    line([stimtime stimtime], [0 yax(2)]','Color','r','LineStyle','--')
    xlabel('Time (sec)')
    ylabel('spikes/sec','Fontsize',12)
%     leg = legend(full_legend_labels);
    title('PSTH')
    xx = [x1 x1 x1+pulse_dur x1+pulse_dur x1];
    yy = [0 yax(2) yax(2) 0 0];
    patch(xx, yy, -1 * ones(size(xx)), [0.9 0.9 0.9], 'LineStyle', 'none')
    
    % 3) PSTH - Blanks
    subplot(243)
        
    if ~isempty(blank_conds)
        plot(edges_stim,psth_blanks(:,1),'k','linewidth',2)
        hold on
        plot(edges_stim,psth_blanks(:,2),'b','linewidth',2)
        xlim([-prestim stimtime+poststim-binsize])  % because points mark the START of the bin
        set(gca,'XMinorTick','on')
        yax = get(gca,'YLim');
        line([0 0], [0 yax(2)]','Color','r','LineStyle','--')
        line([stimtime stimtime], [0 yax(2)]','Color','r','LineStyle','--')
        xlabel('Time (sec)')
        ylabel('spikes/sec','Fontsize',12)
%         leg = legend('Light OFF', sprintf('Light (%s)',legend_labels{1}));
        title('PSTH - Blanks')
        xx = [x1 x1 x1+pulse_dur x1+pulse_dur x1];
        yy = [0 yax(2) yax(2) 0 0];
        patch(xx, yy, -1 * ones(size(xx)), [0.9 0.9 0.9], 'LineStyle', 'none')
    end
    
    % either zoomed PSTH, or size tuning plot
    subplot(244)
    
    if strcmp(exp_type,'size')
        elseif  find(strcmp(Analyzer.loops.conds{1}.symbol,'mask_radius'))    % TEMP - if no blank trials, plot size tuning instead
        sizevar = find(strcmp(Analyzer.loops.conds{1}.symbol,'mask_radius'));
        sizes = unique(cell2mat(conds{sizevar}));   % actual orientations (e.g. 0, 30, 60, etc.)
        num_sizes = length(sizes);
       for i = 1:length(lightconds)
            for o = 1:length(sizes)
                which_sizetrials = find(ismember(trial_type,find(cell2mat(conds{sizevar})==sizes(o))));   % indices of trials at a particular orientation
                which_light_sizetrials = intersect(which_sizetrials,find(all_light_hz == lightconds(i))); % and in diff light conditions
                [spikerate_ev_bysize(i,o), spikerateSE_ev_bysize(i,o)] = calc_firing_rates(spikes_ev(:,n),which_light_sizetrials,lighttime);
            end
        end

        % normalize firing rates
%         if strcmp(normalizing,'blank')
%             norm_spikerate_size = spikerate_ev_bysize - repmat(spikerate_blank(n,:)',1,length(sizes));
%         elseif strcmp(normalizing,'prestim')
%             norm_spikerate_size = spikerate_ev_bysize - repmat(spikerate_prestim(n,:)',1,length(sizes));
%         else
%             norm_spikerate_size = spikerate_ev_bysize;
%         end
        norm_spikerate_size = spikerate_ev_bysize./repmat(max(spikerate_ev_bysize,[],2),1,length(sizes));   % normalize to the peak response in each condition
        SE_transform_size = spikerate_ev_bysize./spikerateSE_ev_bysize;                                          % get relationship b/w FR and SE to adjust SEs for normalized FRs
        norm_spikerateSE_size = spikerateSE_ev_bysize./SE_transform_size;
        
        shadedErrorBar(sizes,norm_spikerate_size(1,:),norm_spikerateSE_size(1,:),{'Color',color_mat(1,:),'linewidth',2},1);
        hold on
        shadedErrorBar(sizes,norm_spikerate_size(end,:),norm_spikerateSE_size(end,:),{'Color',color_mat(end,:),'linewidth',2},1);
        ylabel('Normalized Firing rate (Hz)')
        xlabel('size (degrees)')
        xlim([0 max(sizes)])
        ylim([0 1])
        set(gca,'XTick',sizes(1:2:end));
        yax = get(gca,'YLim');
        title('Size tuning during sustained stimulus')
%         legend('Light OFF', sprintf('Light (%s)',legend_labels{end}),'Location','SouthEast');

    else
        
        % first, draw light pulses (Raster will go on top)
        for c = 1:length(lightconds)
            start_patch = find(type==lightconds(c),1,'first')-1;
            end_patch = find(type==lightconds(c),1,'last')-1;
            for p = 1:(all_light_hz(idx(end_patch)))
                space = 1/(all_light_hz(idx(end_patch)));
                patch([x1+((p-1)*space) x1+((p-1)*space) x1+pulse_dur/1000+((p-1)*space) x1+pulse_dur/1000+((p-1)*space) x1+((p-1)*space)],[start_patch end_patch end_patch start_patch start_patch], [0.9 0.9 0.9], 'LineStyle', 'none')
            end
        end
        for t = 1:num_trials        % for each trial
            tt = idx(t);
            for i=1:length(unit_times{tt})
                if type(t)  % if it's a light condition, make it blue
                    line([unit_times{tt}(i)-prestim unit_times{tt}(i)-prestim], [t-1 t]','Color','b')                  
                else
                    line([unit_times{tt}(i)-prestim unit_times{tt}(i)-prestim], [t-1 t]','Color','k')
                end
            end
        end
        ylim([0 num_trials])
        xlim([.4 .8])
        set(gca,'XMinorTick','on')
        xlabel('Time (sec)')
        ylabel('Trial #')
        for ii = 1:length(unique(type))
            line([.4 1], [trials_per_cond(ii)*(ii-1) trials_per_cond(ii)*(ii-1)], 'Color','b','LineStyle','--')
        end
        title('Raster plot (zoomed in)')
        
    

    end
    
    % 5) Running vs. stationary plot
    subplot(245)
    runbar = bargraph([spikerate_run(n,:); spikerate_stat(n,:)],...
        [spikerateSE_run(n,:); spikerateSE_stat(n,:)]);
    set(get(gca,'YLabel'),'String','Mean FR (spikes/s)','Fontsize',12)
    set(gca,'XTicklabel','Running| Stationary')
    for i = 1:length(lightconds)
        if i == length(lightconds)
            set(runbar(i),'FaceColor',color_mat(end,:),'EdgeColor',color_mat(end,:)); % make sure last lightcond is bright blue
        else
            set(runbar(i),'FaceColor',color_mat(i,:),'EdgeColor',color_mat(i,:));
        end
    end
    legend('Light OFF',legend_labels)
    title('Firing rate - running vs. stationary')
    legend('off')
    
    % 4) Firing rates - Bar plot
    subplot(246)
    FR = [spikerate_prestim(n,:); spikerate_onset(n,:); spikerate_ev(n,:)];
    SE = [spikerateSE_prestim(n,:); spikerateSE_onset(n,:); spikerateSE_ev(n,:)];
    if ~isempty(blank_conds)
        FR = [FR; spikerate_blank(n,:)];
        SE = [SE; spikerateSE_blank(n,:)];
        xcondslabel = 'Prestim| Onset| Evoked| Blanks';
    else
        xcondslabel = 'Prestim| Onset| Evoked';
    end
    evokedbar = bargraph(FR,SE);
    set(get(gca,'YLabel'),'String','Mean FR (spikes/sec)','Fontsize',12)
    set(gca,'XTicklabel',xcondslabel)
    for i = 1:length(lightconds)
        if i == length(lightconds)
            set(evokedbar(i),'FaceColor',color_mat(end,:),'EdgeColor',color_mat(end,:)); % make sure last lightcond is bright blue
        else
            set(evokedbar(i),'FaceColor',color_mat(i,:),'EdgeColor',color_mat(i,:));
        end
    end
    legend('Light OFF',legend_labels)
    title('Firing rate')
    legend('off')

    % 5) Tuning Curve (evoked FR)
    % calculate firing rates by light and orientation conditions
    trialsperori = sum(visual_trials)/length(oris);
    trialsperoricond = trialsperori/length(lightconds); 
    trialsbyori = nan(length(lightconds),trialsperoricond,length(oris));
    
    for i = 1:length(lightconds)
        for o = 1:length(oris)
            which_oritrials = find(ismember(trial_type,find(cell2mat(conds{orivar})==oris(o))));   % indices of trials at a particular orientation
            which_light_oritrials = intersect(which_oritrials,find(all_light_hz == lightconds(i)));
            which_oritrials_run = intersect(which_oritrials,find(move_trials));
            which_oritrials_norun = intersect(which_oritrials,find(~move_trials));
            which_light_oritrials_run = intersect(which_oritrials_run,find(all_light_hz == lightconds(i))); % and in diff light conditions
            which_light_oritrials_norun = intersect(which_oritrials_norun,find(all_light_hz == lightconds(i))); % and in diff light conditions
            if length(which_light_oritrials) < trialsperoricond             % if there are unequal numbers of trials per condition
                missing_trials = trialsperoricond - length(which_light_oritrials);      % add the mean FR of the trials that there are as pseudo-trials
                add_trials = repmat(mean(spikes_ev(which_light_oritrials,n)./lighttime),missing_trials,1);  
                fprintf(sprintf('Lightcond%d ori%d missing %d trial(s) - adding fake trials for Hotellings test\n',i,o,missing_trials))
            elseif length(which_light_oritrials) > trialsperoricond
                toomany_trials = length(which_light_oritrials) - trialsperoricond;
                which_light_oritrials = which_light_oritrials(1:trialsperoricond);      % if too many trials, drop the last one(s)
                add_trials = [];
                fprintf(sprintf('Lightcond%d ori%d has %d too many trial(s) - dropping the last one(s) for Hotellings test\n',i,o,toomany_trials))
            else
                add_trials = [];
            end
            trialsbyori(i,:,o) = [spikes_ev(which_light_oritrials,n)./lighttime; add_trials];       % for T2Hot1 - using all trials regardless of running
            [spikerate_ev_bycond_run(i,o), spikerateSE_ev_bycond_run(i,o)] = calc_firing_rates(spikes_ev(:,n),which_light_oritrials_run,lighttime);
            [spikerate_ev_bycond_norun(i,o), spikerateSE_ev_bycond_norun(i,o)] = calc_firing_rates(spikes_ev(:,n),which_light_oritrials_norun,lighttime);
%             spikes = count_spikes(spiketimes{n},prestim+prelight,prestim+prelight+.25);
%             [spikerate_ev_bycond(i,o), spikerateSE_ev_bycond(i,o)] = calc_firing_rates(spikes,which_light_oritrials,.25);
        end
    end
    
    for i = 1:length(lightconds)
        for o = 1:length(oris)/2
            tuning_curve_nodir(:,o)   = mean([trialsbyori(i,:,o)' trialsbyori(i,:,o+length(oris)/2)'],2);         % average FRs of same orientations but different directions in NO LIGHT condition
        end
        if T2Hot1(tuning_curve_nodir,0.05) < .05
            tuned_cell(n,i) = 1;
        else
            tuned_cell(n,i) = 0;
        end
    end

    
    % normalize firing rates
    if strcmp(normalizing,'blank')
        norm_spikerate_ev = spikerate_ev_bycond_norun - repmat(spikerate_blank(n,:)',1,length(oris));   %  stationary trials 
        norm_spikerateSE_ev = spikerateSE_ev_bycond_norun;
    elseif strcmp(normalizing,'prestim')
        norm_spikerate_ev = spikerate_ev_bycond_norun - repmat(spikerate_prestim(n,:)',1,length(oris));
        norm_spikerateSE_ev = spikerateSE_ev_bycond_norun;
    else
        norm_spikerate_ev = spikerate_ev_bycond_norun;
        norm_spikerateSE_ev = spikerateSE_ev_bycond_norun;
    end
%     norm_spikerate_ev = spikerate_ev_bycond_norun./repmat(max(spikerate_ev_bycond_norun,[],2),1,length(oris));   % normalize to the peak response in each condition
%     norm_spikerate_ev(isnan(norm_spikerate_ev)) = 0;
%     SE_transform = spikerate_ev_bycond_norun./spikerateSE_ev_bycond_norun;                                          % get relationship b/w FR and SE to adjust SEs for normalized FRs
%     norm_spikerateSE_ev = norm_spikerate_ev./SE_transform;
%     norm_spikerateSE_ev(isnan(norm_spikerateSE_ev)) = 0;
%     
%     norm_spikerate_ev_run = spikerate_ev_bycond_run./repmat(max(spikerate_ev_bycond_run,[],2),1,length(oris));   % normalize to the peak response in each condition
%     norm_spikerate_ev_run(isnan(norm_spikerate_ev_run)) = 0;
%     SE_transform_run = spikerate_ev_bycond_run./spikerateSE_ev_bycond_run;                                          % get relationship b/w FR and SE to adjust SEs for normalized FRs
%     norm_spikerateSE_ev_run = norm_spikerate_ev_run./SE_transform_run;
        
    
    
    subplot(247)
    shadedErrorBar(oris,norm_spikerate_ev(1,:),norm_spikerateSE_ev(1,:),{'Color',color_mat(1,:),'linewidth',2},1);  % plot ori tuning in NO RUN trials
    hold on
    shadedErrorBar(oris,norm_spikerate_ev(2,:),norm_spikerateSE_ev(2,:),{'Color',color_mat(2,:),'linewidth',2},1);

    ylabel('Normalized Firing rate (Hz)')
    xlabel('Orientation (degrees)')
    xlim([0 max(oris)])
%     ylim([0 1.1])
    set(gca,'XTick',oris(1:2:end));
    yax = get(gca,'YLim');
    title('Orientation tuning during sustained stimulus')
%     legend('Light OFF', sprintf('Light (%s)',legend_labels{end}),'Location','SouthEast');

    
    % 6) Plot waveforms
    unitname = good_units(n);
    subplot(248)
%     [waveforms_microV{n},t2p(n),max_ch(n)] = phy_waveforms(unitname,probetype,amp_sr);
%     plot([0:size(waveforms_microV{n},2)-1]/div,waveforms_microV{n}','LineWidth',2);
    load(sprintf('Cluster_%s_waveforms.mat',num2str(good_units(n))));
    t = linspace(0,(size(waveforms_microV,1)-1)/20,size(waveforms_microV,1));    
    plot(t,waveforms_microV,'LineWidth',2);
    xlim([0 max(t)])
    hold on 
    title(sprintf('Average waveform of spikes (Max ch: %d)',max_ch))
    ylabel('Amplitude (uV)')
    xlabel('Time (ms)')
%     layers = importdata(sprintf('%s\\layers.mat',exp_path));
%     unit_layer(n) = layers(max_ch(n));
    
    % get unit layer
%     layers(26:30) = 5.5;
%     unit_layer(n) = layers(max_ch);

  % save figs
    cd(all_fig_dir)
    xSize = 24; ySize = 11;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[0 0 xSize*50 ySize*50])
%     if unit_layer(n) == 5
%         layer_name = 'L5A';
%     elseif unit_layer(n) == 5.5
%         layer_name = 'L5B';
%     elseif unit_layer(n) == 2.5
%         layer_name = 'L23';
%     else
%         layer_name = sprintf('L%s',num2str(unit_layer(n)));
%     end
    layer_name = 'LP';
    save_clust_name= sprintf('%s_Cluster_%d_%s',upper(unittype),good_units(n),layer_name);
    print(clust_fig,'-dpng',save_clust_name)
    print2eps(save_clust_name,clust_fig)
%     export_fig (save_clust_name, '-png','-r600','-opengl')
    
%     % test for significant visual modulation
%      onset_spikes_lightoff = spikes_onset(intersect(find(all_light_hz == lightconds(1)),find(visual_trials)),n);
%     onset_spikes_lighton = spikes_onset(intersect(find(all_light_hz == lightconds(2)),find(visual_trials)),n);
%     ev_spikes_lightoff = spikes_ev(intersect(find(all_light_hz == lightconds(1)),find(visual_trials)),n);
%     ev_spikes_lighton = spikes_onset(intersect(find(all_light_hz == lightconds(2)),find(visual_trials)),n);

% kruskal-wallis test to test for significant light- and visual-modulation
    if kruskalwallis([spikes_onset(find(visual_trials(1:num_trials)),n) spikes_prestim_onsettime(find(visual_trials(1:num_trials)),n)],[],'off') < .05
        vis_cell(n) = 1;
    else
        vis_cell(n) = 0;
    end
    ev_spikes_lightoff = spikes_ev(find(~all_light_hz),n);
    ev_spikes_lighton = spikes_ev(find(all_light_hz==lightconds(2)),n);     % currently using second light condition!
    ev_spikes = [ev_spikes_lightoff' ev_spikes_lighton'];
    ev_spikes_group = [ones(1,length(ev_spikes_lightoff)) 2*ones(1,length(ev_spikes_lighton))];
%     if length(ev_spikes_lightoff) ~= length(ev_spikes_lighton)
%         min_len = min(length(ev_spikes_lightoff),length(ev_spikes_lighton));
%         ev_spikes_lightoff = ev_spikes_lightoff(1:min_len);
%         ev_spikes_lighton = ev_spikes_lighton(1:min_len);
%     end
%     if kruskalwallis([ev_spikes_lightoff ev_spikes_lighton],[],'off') < .05    % currently testing against 2nd light condition
    if kruskalwallis(ev_spikes,ev_spikes_group,'off') < .05
        light_cell(n) = 1;
    else
        light_cell(n) = 0;
    end
    
    % COME BACK TO THIS
%    %determine statistical significance of orientation tuning 
%    if length(deg) == 8
%         raw_new = ([ mean([raw_lightoff_sep{1}(1:length(c)); raw_lightoff_sep{5}]); mean([raw_lightoff_sep{2}; raw_lightoff_sep{6}]); mean([raw_lightoff_sep{3}; raw_lightoff_sep{7}]); mean([raw_lightoff_sep{4}; raw_lightoff_sep{8}]) ]) - spontfr_mean;  
%    else
%        raw_new = ([ mean([raw_lightoff_sep{1}(1:length(c)); raw_lightoff_sep{7}]); mean([raw_lightoff_sep{2}; raw_lightoff_sep{8}]); mean([raw_lightoff_sep{3}; raw_lightoff_sep{9}]); mean([raw_lightoff_sep{4}; raw_lightoff_sep{10}]); mean([raw_lightoff_sep{5}; raw_lightoff_sep{11}]); mean([raw_lightoff_sep{6}; raw_lightoff_sep{12}]) ]) - spontfr_mean;
%    end
%             
%    P = T2Hot1(raw_new', 0.05);
    
    % save data
    good_units_results(n).clustnum = unit;
%     good_units_results(n).layer = unit_layer(n);
    good_units_results(n).peakchannel = max_ch;
    good_units_results(n).visuallymodulated = vis_cell(n);
    good_units_results(n).lightmodulated = light_cell(n);
    good_units_results(n).tuned = tuned_cell(n,:);
    good_units_results(n).PSTH_visual = psth;
    if ~isempty(blank_conds)
        good_units_results(n).PSTH_blank = psth_blanks;
    end
    good_units_results(n).FR.label = {'Light OFF','Light ON'};
    good_units_results(n).FR.prestim = FR(1,:);    % NO LIGHT ALWAYS BEFORE LIGHT!
    good_units_results(n).FR.onset = FR(2,:);
    good_units_results(n).FR.evoked = FR(3,:);
    if ~isempty(blank_conds)
        good_units_results(n).FR.blanks = FR(4,:);
    end
    good_units_results(n).orituning.oris = oris;
    good_units_results(n).orituning.FR = spikerate_ev_bycond_norun;
    good_units_results(n).orituning.SE = spikerateSE_ev_bycond_norun;
    good_units_results(n).orituning.FRcorrected = norm_spikerate_ev;
    good_units_results(n).orituning.SEcorrected = norm_spikerateSE_ev;
    if find(strcmp(Analyzer.loops.conds{1}.symbol,'mask_radius'))
        good_units_results(n).sizetuning.sizes = sizes;
        good_units_results(n).sizetuning.FR = spikerate_ev_bysize;
        good_units_results(n).sizetuning.SE = spikerateSE_ev_bysize;
        good_units_results(n).sizetuning.FRcorrected = norm_spikerate_size;
        good_units_results(n).sizetuning.SEcorrected = norm_spikerateSE_size;
    end
%     good_units_results(n).trough2peak = t2p(n);

%     good_units_results(n).slope = poly(n);
    good_units_results(n).waveforms = waveforms_microV;
    
    close all
    
    cd(exp_path)
end

save(strcat(exp_name,'_results.mat'),'good_units_results')
cd(exp_dir)             % also save to cumulative experiment folder
save(strcat(exp_name,'_results.mat'),'good_units_results')
    
  

