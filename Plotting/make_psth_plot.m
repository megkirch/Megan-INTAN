function make_psth_plot(binsize,spikes,which_trials,prestim,stimtime,totaltime,light_trialtypes,light_start,light_dur)

% binsize = size of bins in which to count spikes (in sec)
% spikes = cell array of a single unit's spike times with length of number of trials.
% which_trials = 1xnum_trials vector of 1s and 0s, 1s indicating trials to look
% prestim = time before visual stimulus onset (in sec)
% stimtime = duration of visual stimulus (in sec)
% totaltime = total time of trial(in sec)
% light_trialtypes = 1xnumtrials vector defining each trial's light
    % condition
% light_start = time when the light started (in sec) (if not an opto
    % experiment, should be [])
% light_dur = duration of light pulse (in sec) (if not an opto
    % experiment, should be [])
    
    color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2]; % for graphing purposes (first is black, last is green)


edges = [0:binsize:totaltime-binsize];
psth = make_psth(binsize,edges,which_trials,spikes,light_trialtypes);

edges_stim = [-prestim:binsize:(totaltime-prestim-binsize)]'; % x signifies the timepoint of the START of the bin
for c = 1:size(psth,2)
    plot(edges_stim,psth(:,c),'color',color_mat(c,:),'linewidth',2)
    hold on
end

xlim([-prestim totaltime-prestim-binsize])  % because points mark the START of the bin
set(gca,'XMinorTick','on')
yax = get(gca,'YLim');
line([0 0], [0 yax(2)]','Color','r','LineStyle','--')
line([stimtime stimtime], [0 yax(2)]','Color','r','LineStyle','--')
xlabel('Time (sec)','fontsize',14)
ylabel('spikes/sec','Fontsize',14)

% if it's not an opto experiment, draw patch
if ~isempty(light_start)     
    x1 = light_start - prestim;     % when the light starts
    xx = [x1 x1 x1+light_dur x1+light_dur x1];
    yy = [0 yax(2) yax(2) 0 0];
    patch(xx, yy, -1 * ones(size(xx)), [0.9 0.9 0.9], 'LineStyle', 'none')
end

end
    