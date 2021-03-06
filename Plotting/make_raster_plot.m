function make_raster_plot(spikes,Fs,prestim,totaltime,light_trialtypes,light_start,pulse_dur)

% spikes = cell array of a single unit's spike times with length of number of trials.
% Fs = sampling rate (in sec)
% prestim = time before visual stimulus onset (in sec)
% totaltime = total time of trial(in sec)
% light_trialtypes = 1xnumtrials vector defining each trial's light
    % condition
% light_start = time when the light started (in sec)
% pulse_dur = duration of light pulse (in sec)
    % pulse_dur in trains experiments)

spike_raster = make_raster(spikes,Fs,totaltime);
num_trials = length(spikes);
    color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2]; % for graphing purposes (first is black, last is green)

time_vec = linspace(-prestim,totaltime-prestim,size(spike_raster,2));      % for x-axis 
cond_colors = ones(1,num_trials);      % just in case you have multiple conditions

% if light experiment, add blue shading on raster
if length(unique(light_trialtypes)) > 1         % it was an opto experiment if there's more than one "light" condition (one condition would be nolight)
    [cond,idx] = sort(light_trialtypes);    % sort trials by light conditions
    diff_conds = unique(cond);
    new_spike_rast = spike_raster(idx,:);             % reordered according to light condition
    new_spike_rast(find(~new_spike_rast)) = nan;        % so that you don't plot absence of spikes as y=0
    for c = 1:length(diff_conds)
        cond_colors(find(cond==diff_conds(c)))=c;
    end
    start_patch = find(cond>0,1,'first');     % first trial with light
    end_patch = num_trials;         % assumes all trials where cond > 0 had light
    x1 = light_start-prestim;       % when the light starts

    for c = 1:length(diff_conds)
        start_patch = find(cond==diff_conds(c),1,'first')-1;
        end_patch = find(cond==diff_conds(c),1,'last')-1;
        if sum(diff_conds>10)         % if it was a trains experiment
            for p = 1:(light_trialtypes(idx(end_patch)))
                space = 1/(light_trialtypes(idx(end_patch)));
                patch([x1+((p-1)*space) x1+((p-1)*space) x1+pulse_dur+((p-1)*space) x1+pulse_dur+((p-1)*space) x1+((p-1)*space)],[start_patch end_patch end_patch start_patch start_patch], [0.9 0.9 0.9], 'LineStyle', 'none')
            end
        elseif diff_conds(c)         % any other light condition
            patch([x1 x1 x1+pulse_dur x1+pulse_dur x1],[start_patch end_patch end_patch start_patch start_patch], [0.9 0.9 0.9], 'LineStyle', 'none')
        end
        line([-prestim totaltime-prestim], [end_patch end_patch]','Color','r','LineStyle','--')
    end
    

%         patch([x1 x1 border border x1],[start_patch end_patch end_patch start_patch start_patch], [0.9 0.9 0.9], 'LineStyle', 'none')
    hold on
end

if exist('xrange','var')        
    rasterrange = [find(time_vec==.4):find(time_vec==.8)];
else
    rasterrange = 1:size(new_spike_rast,2);
end

% draw raster
for t = 1:num_trials     % for each trial
    new_spike_rast(t,find(new_spike_rast(t,:))) = new_spike_rast(t,find(new_spike_rast(t,:)))+(t-1);        % plot 'ones' in vector so that y = trial number
    plot(time_vec(rasterrange),new_spike_rast(t,rasterrange),'.','Color',color_mat(cond_colors(t),:),'MarkerSize',6)    
    hold on
end
ylim([0 num_trials])
xlim([-prestim totaltime-prestim])

xlabel('Time from visual stimulus onset (sec)','fontsize',14)
ylabel('Trial (by condition)','fontsize',14)

return