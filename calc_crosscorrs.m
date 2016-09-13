function calc_crosscorrs

exp_paths = {'H:\Tlx3project\T18\5-3-16\T18_diffintensities_160503_171447',...
'H:\Tlx3project\T18\5-5-16\T18_run2_diffintensities_160505_092232',...
'H:\Tlx3project\T19\5-2-16\T19_diffintensities_160502_210213',...
'H:\Tlx3project\T19\5-3-16\T19_run2_diffintensities_160503_221201'};

% make rasters for multiple experiments
for i = 1:length(exp_paths)
    exp_path = exp_paths{i};
    load(sprintf('%s\\su_spikedata.mat',exp_path))
    
    % get experiment name
    animal_name = regexp(exp_paths{i},'T\d+','match');
    animal_name = animal_name{1};
    run_name = regexp(exp_paths{i},'run\d','match');
    if isempty(run_name)
        exp_name = animal_name;
    else
        exp_name = strcat(animal_name,run_name{1});
    end
    
    results = importdata(sprintf('%s\\%s_results.mat',exp_paths{i},exp_name));
    for n=1:length(spiketimes)          % for each unit in that experiment
        spikes = spiketimes{n};
        spike_raster{i}(n,:,:) = make_raster(spikes,1000,2.5);
        unit{i}(n) = results(n).clustnum;
        layer{i}(n) = results(n).layer;
    end
end

cd 'H:\Tlx3project\CumulativeData\T18&T19\CrossCorrs'

% for each experiment, calculate FRs in first and last pulse of long train
count = 1;
for i = 1:length(exp_paths)
    load(sprintf('%s\\light_params.mat',exp_paths{i}))
    light_conds = unique(all_light);
    light_trials = find(all_light ==light_conds(end));         % TEMP
    nolight_trials = find(all_light==light_conds(1));
    time = 0:1/1000:2.5;
    lightstart = find(time==round(av_light_start));        % SAMPLE at which light started
    [layer_out,inds] = sort(layer{i});
    light_raster(:,:) = mean(spike_raster{i}(:,light_trials,lightstart:lightstart+999),2);
    new_lightrast = light_raster(inds,:);
    [C,lags] = xcorr(new_lightrast','coeff');
    nolight_raster(:,:) = mean(spike_raster{i}(:,nolight_trials,lightstart:lightstart+999),2);
    new_nolightrast = nolight_raster(inds,:);
    [D,lags] = xcorr(new_nolightrast','coeff');
    test = reshape(C(1000,:),size(spike_raster{i},1),size(spike_raster{i},1));
    test2 = reshape(D(1000,:),size(spike_raster{i},1),size(spike_raster{i},1));
    fig1=figure;imagesc(test);colorbar
    fig2=figure;imagesc(test2);colorbar
    fig3=figure;imagesc(test-test2);colorbar
    print(fig1,'-dpng',sprintf('Exp%d_Lighton',i))
    print(fig2,'-dpng',sprintf('Exp%d_Lightoff',i))
    print(fig3,'-dpng',sprintf('Exp%d_Lighton-Lightoff',i))
    close all
    clear light_raster nolight_raster
end
        
end