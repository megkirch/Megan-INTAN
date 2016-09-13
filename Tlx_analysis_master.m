function Tlx_analysis_master(exp_type)

main_dir = 'H:\Tlx3project\Augustresults';
cd(main_dir)
fig_dir =  sprintf('H:\\Tlx3project\\Augustresults\\%s_figs',exp_type);

if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

if strcmp(exp_type,'intensities')

    exp_paths = {'H:\Tlx3project\T18\5-3-16\T18_diffintensities_160503_171447',...
    'H:\Tlx3project\T18\5-5-16\T18_run2_diffintensities_160505_092232',...
    'H:\Tlx3project\T19\5-2-16\T19_diffintensities_160502_210213',...
    'H:\Tlx3project\T19\5-3-16\T19_run2_diffintensities_160503_221201',...
    'H:\Tlx3project\T22\7-26-16\T22_diffintensities_160726_171349'};

elseif strcmp(exp_type,'ramp')
    exp_paths = {'H:\Tlx3project\T14\T14_driftgratingwithramp_160329_175645',...
        'H:\Tlx3project\T15\4-1-16\T15_driftgratingramp_160325_205545',...
        'H:\Tlx3project\T18\5-3-16\T18_ramp_160503_182443',...
        'H:\Tlx3project\T18\5-5-16\T18_run2_ramp_160505_103140',...
        'H:\Tlx3project\T19\5-3-16\T19_run2_ramp_160503_211551',...
        'H:\Tlx3project\T22\7-26-16\T22_ramp_160726_161333',...
        };

elseif strcmp(exp_type,'trains')

exp_paths = {'H:\Tlx3project\T9\T9_driftgrating_trains_160125_123706',...
        'H:\Tlx3project\T14\T14_driftgratingtrains_160329_211516',...
        'H:\Tlx3project\T15\3-31-16\T15_driftgratingtrains_160324_214136',...
        'H:\Tlx3project\T19\5-2-16\T19_trains_160502_222135',...
        'H:\Tlx3project\T18\5-3-16\T18_trains_160503_192729',...
        'H:\Tlx3project\T22\7-26-16\T22_trains_160726_183114'};
    
elseif strcmp(exp_type,'inhibit')
    exp_paths = {'H:\Tlx3project\TiC1\TiC1_step_160901_141812',...
        'H:\Tlx3project\TiC2\8-31-16\TiC2_step_160831_121353'};
    
end

count = 1;
for i = 1:length(exp_paths)
    exp_path = exp_paths{i};
    % get experiment name
    out = regexp(exp_path,'\\','split');
    inds = regexp(out{end},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
    exp_name = out{end}(1:inds(1)-1);
    exp_dir = strcat(main_dir,'\',exp_name);
    fprintf(sprintf('Processing experiment %s\n',exp_name))
%     if i == 1
%         intan_analysis_master(exp_path,exp_type);
%     end
    
    % load results
    cd(exp_dir)
    s = dir; 
    for ii=1:length(s)
        if strfind(s(ii).name,'_results') 
            results_file = s(ii).name;
        end
    end
    S{i} = importdata(results_file,'-mat');     % load results mat
    load(sprintf('%s\\light_params.mat',exp_path))
    lightconds = unique(all_light);
end

count = 0;
for s = 1:length(S)
    for n = 1:S{s}.numunits
        if length(S{s}.lightconds)>4        % TEMP - HACK - skip 10hz condition
            which_ones = [1 3 4];
        else
            which_ones = [1:length(S{s}.lightconds)-1];
        end
        % which light condition to use for comparisons
%         if strcmp(exp_type,'trains')
%             which_lightcond = length(S{s}.lightconds);    % use highest frequency   % change to make best for particular unit!!
%         else
%             which_lightcond(count) = 2;     % use first light condition (e.g. lowest intensity)
%         end

%         if S{s}.FRs(n).ev(1) > .25          % pretty liberal cut-off for which cells to include
            count = count+1;
            name{count} = S{s}.unitinfo(n).name;
            layer(count) = S{s}.unitinfo(n).layer;
            vissig(count) = S{s}.unitinfo(n).vissig;
            vissig2(count) = S{s}.unitinfo(n).vissig2;
            lightsig_vis(count,:) = S{s}.unitinfo(n).lightsig_vis(which_ones);
            lightsig_blnk(count,:) = S{s}.unitinfo(n).lightsig_blnk(which_ones);
            if length(S{s}.lightconds)>4        % TEMP - HACK
                FRev(count,:) = S{s}.FRs(n).ev([1 which_ones+1]);     % need to fix - currently contains both FRs and SEs
            else
                FRev(count,:) = S{s}.FRs(n).ev(1:length(S{s}.FRs(n).ev)/2);     % need to fix - currently contains both FRs and SEs
            end

            for ii = 2:length(FRev(count,:))
                evdiff(ii-1) = diff(FRev(count,[1 ii]),[],2);       % get which light condition resulted in largest change in FR for each unit
            end
            [~,which_lightcond(count)] = max(abs(evdiff),[],2);
            which_lightcond(count) = which_lightcond(count)+1;

            baseline(count) = S{s}.FRs(n).baseline;
            if ~isempty(S{s}.FRs(n).SALT10ms)
                SALT10ms(count,:) = S{s}.FRs(n).SALT10ms(which_ones);
                SALT5ms(count,:) = S{s}.FRs(n).SALT5ms(which_ones);
                SALT2ms(count,:) = S{s}.FRs(n).SALT2ms(which_ones);
            end
            OSI(count,:) = S{s}.tuning(n).OSI(([1 which_ones+1]));
            OSI_CV(count,:) = S{s}.tuning(n).OSI_CV([1 which_ones+1]);
            DSI(count,:) = S{s}.tuning(n).DSI([1 which_ones+1]);
            DSI_CV(count,:) = S{s}.tuning(n).DSI_CV([1 which_ones+1]);
            sig_Hot(count,:) = S{s}.tuning(n).sig_Hot([1 which_ones+1]);
            sig_Hot2(count,:) = S{s}.tuning(n).sig_Hot2([1 which_ones+1]);
            tuning_curve_collapse = zeros(length(S{s}.lightconds),size(S{s}.tuning(n).curve,2)/2);      % reset
            for o = 1:size(S{s}.tuning(n).curve,2)/2
                tuning_curve_collapse(:,o)   = mean([S{s}.tuning(n).curve(:,o) S{s}.tuning(n).curve(:,o+size(S{s}.tuning(n).curve,2)/2)],2);         % average evoked FR of same orientations but different directions in NO LIGHT condition
            end
            [~,pref_deg] = max(tuning_curve_collapse(1,:));
            if length(pref_deg) > 1
                [~,prefpref] = max(tuning_curve_collapse(2,pref_deg));
                pref_deg = pref_deg(prefpref);
            end
            prefFR_delta(count) =  abs(diff(tuning_curve_collapse([1 which_lightcond(count)],pref_deg)));      % when more than two light conditions, using only first light condition vs no light
            if pref_deg <= size(S{s}.tuning(n).curve,2)/4
                orth_deg = pref_deg+size(S{s}.tuning(n).curve,2)/4;
            else
                orth_deg = pref_deg-size(S{s}.tuning(n).curve,2)/4;
            end
            % check if cells should be excluded from tuning analysis (no FR
                % at any orientation <-.5 or >-.5)
            if abs(tuning_curve_collapse(1,pref_deg))<.5 || abs(tuning_curve_collapse(2,pref_deg)) < .5
                include(count) = 0;
            else
                include(count) = 1;
            end
            
            orthFR_delta(count) =  abs(diff(tuning_curve_collapse([1 which_lightcond(count)],orth_deg)));
    %         sig_Hot_polar(count,:) = S{s}.tuning(n).sig_Hot_polar([1 which_ones+1]);
            t2p_t(count) = S{s}.waveforms(n).t2p_t;
            t2p_r(count) = S{s}.waveforms(n).t2p_r;
            fwhm(count) = S{s}.waveforms(n).fwhm;

            lightmod(count) = diff(FRev(count,[1 which_lightcond(count)]))./sum(FRev(count,[1 which_lightcond(count)]));
%         end
    end
end
        

% get different cell types
visual_cells = find((vissig < .05)|(vissig2 < .05));
light_cells= find(min(lightsig_vis,[],2)<.05);   % find cells with significant effect in any light condition
tuned_cells = find(sig_Hot(:,1) < .05);
tuned_cells2 = find(sig_Hot2(:,1) < .05);
%  tuned_cells(find(FRev(tuned_cells,1)<.5)) = [];        % exclude cells with <.5 baseline FR from being considered "tuned"
tuned_cells = tuned_cells(ismember(tuned_cells,find(include)));
if strcmp(exp_type,'trains')
    onset_cells = find(mean(SALT10ms,2)<.05);
    quikonset_cells = find(mean(SALT5ms,2)<.05);
    tlx_cells = find(mean(SALT2ms,2)<.05);
elseif strcmp(exp_type,'intensities')
    onset_cells = find(min(SALT10ms,[],2)<.05);
    quikonset_cells = find(min(SALT5ms,[],2)<.05);
    tlx_cells = find(min(SALT2ms,[],2)<.05);
end
FS_cells = find(t2p_t < .65);
reg_cells = find(t2p_t >= .65);

chs_23 = find(layer == 2.5);
chs_4 = find(layer == 4);
chs_5 = find(layer == 5);
chs_55 = find(layer == 5.5);
chs_6 = find(layer == 6);
infgran = [chs_5 chs_55 chs_6];

% get means and medians for overall (1st row), supra (2nd) gran (3rd) and infragran (4th) layers for
%reg-spiking (first column) and fast-spiking (second column)
median_sum = [nanmedian(OSI_CV(reg_cells,1)) nanmedian(OSI_CV(FS_cells,1));
                 nanmedian(OSI_CV(intersect(chs_23,reg_cells),1)) nanmedian(OSI_CV(intersect(chs_23,FS_cells),1));
                 nanmedian(OSI_CV(intersect(chs_4,reg_cells),1)) nanmedian(OSI_CV(intersect(chs_4,FS_cells),1));
                 nanmedian(OSI_CV(intersect(infgran,reg_cells),1)) nanmedian(OSI_CV(intersect(infgran,FS_cells),1))];
mean_sum = [nanmean(OSI_CV(reg_cells,1)) nanmean(OSI_CV(FS_cells,1));
                 nanmean(OSI_CV(intersect(chs_23,reg_cells),1)) nanmean(OSI_CV(intersect(chs_23,FS_cells),1));
                 nanmean(OSI_CV(intersect(chs_4,reg_cells),1)) nanmean(OSI_CV(intersect(chs_4,FS_cells),1));
                 nanmean(OSI_CV(intersect(infgran,reg_cells),1)) nanmean(OSI_CV(intersect(infgran,FS_cells),1))];
groupn = [length(~isnan(OSI_CV(reg_cells,1))) length(~isnan(OSI_CV(FS_cells,1)));
 length(isnan(OSI_CV(intersect(chs_23,reg_cells),1))) length(~isnan(OSI_CV(intersect(chs_23,FS_cells),1)));
 length(isnan(OSI_CV(intersect(chs_4,reg_cells),1))) length(~isnan(OSI_CV(intersect(chs_4,FS_cells),1)));
 length(isnan(OSI_CV(intersect(infgran,reg_cells),1))) length(~isnan(OSI_CV(intersect(infgran,FS_cells),1)))];
             
cd(fig_dir)
barplot_by_layer(lightmod,layer,ones(1,length(layer)),'Light modulation index','Lightmodulation (evoked)')
barplot_by_layer(OSI_CV',layer,ones(1,length(layer)),'Orientation selectivity index (CV)','OSI (CV)')
barplot_by_layer(DSI_CV',layer,ones(1,length(layer)),'Direction selectivity index (CV)','DSI (CV)')
barplot_by_layer(OSI',layer,ones(1,length(layer)),'Orientation selectivity index','OSI')
barplot_by_layer(DSI',layer,ones(1,length(layer)),'Direction selectivity index','DSI')
barplot_by_layer(diff(OSI_CV(tuned_cells,:),[],2)',layer(tuned_cells),ones(1,length(layer(tuned_cells))),'OSI CV change (light-no light)','OSI change')
barplot_by_layer(diff(DSI_CV(tuned_cells,:),[],2)',layer(tuned_cells),ones(1,length(layer(tuned_cells))),'DSI CV change (light-no light)','DSI change')
if ~strcmp(exp_type,'ramp')
    barplot_by_layer(lightmod(onset_cells),layer(onset_cells),ones(1,length(onset_cells)),'Light modulation index','Lightmodulation (onset cells)')
end

% FS cells
FSFR_fig = figure('name','FS firing rate: No Light vs. Light');
plot(FRev(intersect(chs_23,FS_cells),1),FRev(intersect(chs_23,FS_cells),mode(which_lightcond)),'r.','MarkerSize',24)
hold on
plot(FRev(intersect(chs_4,FS_cells),1),FRev(intersect(chs_4,FS_cells),mode(which_lightcond)),'b.','MarkerSize',24)
plot(FRev(intersect(chs_5,FS_cells),1),FRev(intersect(chs_5,FS_cells),mode(which_lightcond)),'g.','MarkerSize',24)
plot(FRev(intersect(chs_55,FS_cells),1),FRev(intersect(chs_55,FS_cells),mode(which_lightcond)),'k.','MarkerSize',24)
plot(FRev(intersect(chs_6,FS_cells),1),FRev(intersect(chs_6,FS_cells),mode(which_lightcond)),'c.','MarkerSize',24)
xlabel('Firing rate (spks/s) - No Light','Fontsize',24)
ylabel('Firing rate (spks/s) - Light','Fontsize',24)
xax = get(gca,'XLim');
yax = get(gca,'YLim');
x = [min(xax(1),yax(1)):1:max(xax(2),yax(2))];
y=x;
plot(x,y,'k');
xlim([min(xax(1),yax(1)), max(xax(2),yax(2))])
ylim([min(xax(1),yax(1)), max(xax(2),yax(2))])
% xlim([0 25])
% ylim([0 25])
set(gca,'fontsize',18)
l=legend('L2/3','L4','L5A','L5B','L6','Location','SouthEast');
set(l,'fontsize',18)
print(FSFR_fig, '-dpng','FiringRate_FS')

% RS cells
RSFR_fig = figure('name','RS firing rate: No Light vs. Light');
plot(FRev(intersect(chs_23,reg_cells),1),FRev(intersect(chs_23,reg_cells),mode(which_lightcond)),'r.','MarkerSize',24)
hold on
plot(FRev(intersect(chs_4,reg_cells),1),FRev(intersect(chs_4,reg_cells),mode(which_lightcond)),'b.','MarkerSize',24)
plot(FRev(intersect(chs_5,reg_cells),1),FRev(intersect(chs_5,reg_cells),mode(which_lightcond)),'g.','MarkerSize',24)
plot(FRev(intersect(chs_55,reg_cells),1),FRev(intersect(chs_55,reg_cells),mode(which_lightcond)),'k.','MarkerSize',24)
plot(FRev(intersect(chs_6,reg_cells),1),FRev(intersect(chs_6,reg_cells),mode(which_lightcond)),'c.','MarkerSize',24)
xlabel('Firing rate (spks/s) - No Light','Fontsize',24)
ylabel('Firing rate (spks/s) - Light','Fontsize',24)
xax = get(gca,'XLim');
yax = get(gca,'YLim');
x = [min(xax(1),yax(1)):1:max(xax(2),yax(2))];
y=x;
plot(x,y,'k');
xlim([min(xax(1),yax(1)), max(xax(2),yax(2))])
ylim([min(xax(1),yax(1)), max(xax(2),yax(2))])
% xlim([0 25])
% ylim([0 25])
set(gca,'fontsize',18)
l=legend('L2/3','L4','L5A','L5B','L6','Location','SouthEast');
set(l,'fontsize',18)
print(RSFR_fig, '-dpng','FiringRate_RS')

osi_fig = figure('name','OSI Change (light-nolight) vs. Light modulation');
plot(diff(OSI_CV(intersect(tuned_cells,chs_23),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_23)),'r.','MarkerSize',24)
hold on
plot(diff(OSI_CV(intersect(tuned_cells,chs_4),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_4)),'b.','MarkerSize',24)
plot(diff(OSI_CV(intersect(tuned_cells,chs_5),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_5)),'g.','MarkerSize',24)
plot(diff(OSI_CV(intersect(tuned_cells,chs_55),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_55)),'k.','MarkerSize',24)
plot(diff(OSI_CV(intersect(tuned_cells,chs_6),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_6)),'c.','MarkerSize',24)
xlabel('change in OSI CV (light - no light)','Fontsize',24)
ylabel('Light modulation index','Fontsize',24)
line([0 0],[-1 1],'Color','k')
line([-1 1],[0 0],'Color','k')
set(gca,'fontsize',18)
l=legend('L2/3','L4','L5A','L5B','L6','Location','SouthEast');
set(l,'fontsize',18)
print(osi_fig, '-dpng','OSIvslight')
% print2eps('OSIvslight', osi_fig)

osi_cellfig = figure('name','OSI Change (light-nolight) vs. Light modulation');
plot(diff(OSI_CV(intersect(reg_cells,tuned_cells),[1 mode(which_lightcond)]),[],2),lightmod(intersect(reg_cells,tuned_cells)),'r.','MarkerSize',24)
hold on
plot(diff(OSI_CV(intersect(FS_cells,tuned_cells),[1 mode(which_lightcond)]),[],2),lightmod(intersect(FS_cells,tuned_cells)),'b.','MarkerSize',24)
xlabel('change in OSI CV (light - no light)','Fontsize',24)
ylabel('Light modulation index','Fontsize',24)
line([0 0],[-1 1],'Color','k')
line([-1 1],[0 0],'Color','k')
set(gca,'fontsize',18)
l=legend('Regular spiking','Fast spiking');
set(l,'fontsize',18)
print(osi_cellfig, '-dpng','OSIvslight_bycelltype')

dsi_cellfig = figure('name','DSI Change (light-nolight) vs. Light modulation');
plot(diff(DSI_CV(intersect(reg_cells,tuned_cells),[1 mode(which_lightcond)]),[],2),lightmod(intersect(reg_cells,tuned_cells)),'r.','MarkerSize',24)
hold on
plot(diff(DSI_CV(intersect(FS_cells,tuned_cells),[1 mode(which_lightcond)]),[],2),lightmod(intersect(FS_cells,tuned_cells)),'b.','MarkerSize',24)
xlabel('change in DSI CV (light - no light)','Fontsize',24)
ylabel('Light modulation index','Fontsize',24)
line([0 0],[-1 1],'Color','k')
line([-1 1],[0 0],'Color','k')
set(gca,'fontsize',18)
l=legend('Regular spiking','Fast spiking');
set(l,'fontsize',18)
print(dsi_cellfig, '-dpng','DSIvslight_bycelltype')

dsi_fig = figure('name','OSI Change (light-nolight) vs. Light modulation');
plot(diff(DSI_CV(intersect(tuned_cells,chs_23),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_23)),'r.','MarkerSize',24)
hold on
plot(diff(DSI_CV(intersect(tuned_cells,chs_4),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_4)),'b.','MarkerSize',24)
plot(diff(DSI_CV(intersect(tuned_cells,chs_5),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_5)),'g.','MarkerSize',24)
plot(diff(DSI_CV(intersect(tuned_cells,chs_55),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_55)),'k.','MarkerSize',24)
plot(diff(DSI_CV(intersect(tuned_cells,chs_6),[1 mode(which_lightcond)]),[],2),lightmod(intersect(tuned_cells,chs_6)),'c.','MarkerSize',24)
xlabel('change in DSI CV (light - no light)','Fontsize',24)
ylabel('Light modulation index','Fontsize',24)
line([0 0],[-1 1],'Color','k')
line([-1 1],[0 0],'Color','k')
set(gca,'fontsize',18)
l=legend('L2/3','L4','L5A','L5B','L6','Location','SouthEast');
set(l,'fontsize',18)
print(dsi_fig, '-dpng','DSIvslight')
% print2eps('OSIvslight', dsi_fig)

plot_OSI(layer,OSI_CV,'OSI (CV)')      
plot_OSI(layer,DSI_CV,'DSI (CV)') 

% additive vs multiplicative changes
change_by_celltype = figure('name','Change in orthogonal vs preferred FR');
plot(orthFR_delta(reg_cells),prefFR_delta(reg_cells),'r.','MarkerSize',24)
hold on
plot(orthFR_delta(FS_cells),prefFR_delta(FS_cells),'b.','MarkerSize',24)
xax = get(gca,'XLim');
yax = get(gca,'YLim');
x = [min(xax(1),yax(1)):1:max(xax(2),yax(2))];
y=x;
plot(x,y,'k');
coeffs(1,:) = polyfit(orthFR_delta(reg_cells), prefFR_delta(reg_cells), 1);
coeffs(2,:) = polyfit(orthFR_delta(FS_cells), prefFR_delta(FS_cells), 1);
% Get fitted values
fittedX = linspace(min(orthFR_delta), max(orthFR_delta), 200);
for i = 1:size(coeffs,1)
    fittedY(i,:) = polyval(coeffs(i,:), fittedX);
end
plot(fittedX,fittedY(1,:),'r')
plot(fittedX,fittedY(2,:),'b')
xlabel('Change in orthogonal FR','fontsize',24)
ylabel('Change in preferred FR','fontsize',24)
set(gca,'fontsize',18)
l=legend('Regular spiking','Fast spiking');
set(l,'fontsize',18)
print(change_by_celltype, '-dpng','FRchange_bycelltype')

% diff layers
change_by_layer = figure('name','Change in orthogonal vs preferred FR');
plot(orthFR_delta(chs_23),prefFR_delta(chs_23),'r.','MarkerSize',24)
hold on
plot(orthFR_delta(chs_4),prefFR_delta(chs_4),'b.','MarkerSize',24)
plot(orthFR_delta(chs_5),prefFR_delta(chs_5),'g.','MarkerSize',24)
plot(orthFR_delta(chs_55),prefFR_delta(chs_55),'k.','MarkerSize',24)
plot(orthFR_delta(chs_6),prefFR_delta(chs_6),'c.','MarkerSize',24)
xax = get(gca,'XLim');
yax = get(gca,'YLim');
x = [min(xax(1),yax(1)):1:max(xax(2),yax(2))];
y=x;
plot(x,y,'k--');
coeffs(1,:) = polyfit(orthFR_delta(chs_23), prefFR_delta(chs_23), 1);
coeffs(2,:) = polyfit(orthFR_delta(chs_4), prefFR_delta(chs_4), 1);
coeffs(3,:) = polyfit(orthFR_delta(chs_5), prefFR_delta(chs_5), 1);
coeffs(4,:) = polyfit(orthFR_delta(chs_55), prefFR_delta(chs_55), 1);
coeffs(5,:) = polyfit(orthFR_delta(chs_6), prefFR_delta(chs_6), 1);
% Get fitted values
fittedX = linspace(min(orthFR_delta), max(orthFR_delta), 200);
for i = 1:size(coeffs,1)
    fittedY(i,:) = polyval(coeffs(i,:), fittedX);
end
plot(fittedX,fittedY(1,:),'r')
plot(fittedX,fittedY(2,:),'b')
plot(fittedX,fittedY(3,:),'g')
plot(fittedX,fittedY(4,:),'k')
plot(fittedX,fittedY(5,:),'c')
xlabel('Change in orthogonal FR','fontsize',24)
ylabel('Change in preferred FR','fontsize',24)
set(gca,'fontsize',18)
l=legend('L2/3','L4','L5A','L5B','L6','Location','SouthEast');
set(l,'fontsize',18)
print(change_by_layer, '-dpng','FRchange_bylayer')

end
