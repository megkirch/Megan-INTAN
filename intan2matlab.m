function intan2matlab(exp_path)

% Takes in .dat files in 'one file per signal type' format from Intan
% (see http://www.intantech.com/files/Intan_RHD2000_data_file_formats.pdf)
%
% MAK - 7/9/16
% modified from intanphy2matlab_v2 so that it automatically recognizes 
% which inputs are LED, epoch, etc. 
%
% exp_path - establish directory that contains .kwik file as well as files 
% from intan
% ex. 'E:\MK\9-4-15\09042015_189_150904_145916'; 

cd(exp_path)
kwik_file = sprintf('%s\\amplifier.kwik',exp_path); 

% read header file
[amplifier_channels,aux_input_channels,board_adc_channels,...
    board_dig_in_channels,spike_triggers,supply_voltage_channels,...
    frequency_parameters] = read_Intan_RHD2000_file;    % script from Intan
amp_sr = frequency_parameters.amplifier_sample_rate;    % amplifier sampling rate

% read time.dat file from intan
fileinfo = dir('time.dat');
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
fid = fopen('time.dat', 'r');
time = fread(fid, num_samples, 'int32');     
fclose(fid);
time = time / amp_sr; % time index of sample (in seconds)

% % read auxillary.dat from intan
% num_channels = length(aux_input_channels); % aux input channel info from header file
% fileinfo = dir('auxiliary.dat');
% num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
% fid = fopen('auxiliary.dat', 'r');
% aux_v = fread(fid, [num_channels, num_samples], 'uint16');
% fclose(fid);
% aux_v = aux_v * 0.0000374; % convert to volts

% % read supply.dat from intan
% fileinfo = dir('supply.dat');
% num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
% fid = fopen('supply.dat', 'r');
% supply_v = fread(fid, [num_channels, num_samples], 'uint16');
% fclose(fid);
% supply_v = supply_v * 0.0000748; % convert to volts

% read analogin.dat from intan
num_channels = length(board_adc_channels); % ADC input info from header file
fileinfo = dir('analogin.dat');
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen('analogin.dat', 'r');
adc = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
adc = adc * 0.000050354; % convert to volts

% read digitalin.dat from intan
num_board_dig_in_channels = length(board_dig_in_channels); % ADC input info from header file
fileinfo = dir('digitalin.dat');
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen('digitalin.dat', 'r');
dig_in_raw = fread(fid, num_samples, 'uint16');
dig_in = zeros(num_board_dig_in_channels, num_samples);
fclose(fid);
% Extract digital input channels to separate variables.
for i=1:num_board_dig_in_channels
   mask = 2^(board_dig_in_channels(i).native_order) * ones(size(dig_in_raw));
   dig_in(i, :) = (bitand(dig_in_raw, mask) > 0);
end

% assign inputs to appropriate vectors
possible_inputs = {'photodiode','led','LED','laser','epoch','opt_mouse',...
    'encoder','encdA','encdB'};     % for string comparison
var_names = {'photo','LED','LED','LED',...
    'epoc','mouse','encdA','encdA','encdB'};
diff_vars = unique(var_names);
used_vars = zeros(1,length(var_names));   % keep track of which var_names were used and whether they're digital (1) or analog (2)
disp('Assign inputs to appropriate variables')
% analog inputs
for i = 1:size(adc,1)
    which_input = strfind(possible_inputs,board_adc_channels(i).custom_channel_name) ;  % find which possible input corresponds to current input
    ind = find(cellfun(@(x) ~isempty(x),which_input,'UniformOutput',1));
    used_vars(ind) = 2;
    if ~isempty(ind)   % if this input matches any of the predetermined possible inputs
        fprintf(1, 'Analog input %d: %s', i, board_adc_channels(i).custom_channel_name);
        fprintf(1, '\n');
        eval([sprintf('%s_input',var_names{ind}),sprintf('=%d',i)]);
    else            % if names don't match
        fprintf(1, 'Analog input %d: %s', i, board_adc_channels(i).custom_channel_name);
        fprintf(1, '\n');
        disp('Doesnt match possible variable names. Choose from the following:')
        fprintf(1, '\n');
        for ii = 1:length(diff_vars)
            fprintf(1, '%s: %d', diff_vars{ii},ii);
            fprintf(1, '\n');
        end
        chk = input('Which variable should this input be assigned to?: ','s');
        eval([sprintf('%s_input',diff_vars{chk}),sprintf('=%d',i)]);   
        used_vars(ind) = 2;
    end
end
% digital inputs
for i = 1:size(dig_in,1)
    which_input = strfind(possible_inputs,board_dig_in_channels(i).custom_channel_name) ;  % find which possible input corresponds to current input
    ind = find(cellfun(@(x) ~isempty(x),which_input,'UniformOutput',1));
    used_vars(ind) = 1;
    if ~isempty(ind)   % if this input matches any of the predetermined possible inputs
        fprintf(1, 'Digital input %d: %s', i, board_dig_in_channels(i).custom_channel_name);
        fprintf(1, '\n');
        eval([sprintf('%s_input',var_names{ind}),sprintf('=%d',i)]);
    else            % if names don't match
        fprintf(1, 'Digital input %d: %s', i, board_dig_in_channels(i).custom_channel_name);
        fprintf(1, '\n');
        fprintf(1, 'Doesnt match possible variable names. Choose from the following:');
        fprintf(1, '\n');
        for ii = 1:length(diff_vars)
            fprintf(1, '%s: %d', diff_vars{ii},ii);
            fprintf(1, '\n');
        end
        chk = input('Which variable should this input be assigned to?: ','s');
        eval([sprintf('%s_input',diff_vars{chk}),sprintf('=%d',i)]);   
        used_vars(ind) = 1;
    end
end


% extract appropriate inputs
epoc = dig_in(epoc_input,:); 
photo = adc(photo_input,:);   
if exist('LED_input','var')
    if sum(used_vars(2:4)) == 1 % because indices 2-4 are all for LED_input
        LED = dig_in(LED_input,:);
    else
        LED = adc(LED_input,:);
    end
end
if exist('encdA_input','var')     % assumes that if encdA exists, so does encdB
    if sum(used_vars(7:8)) == 1
        encdA = dig_in(encdA_input,:); % indices 7-8 are for encdA_input
    else
        encdA = adc(encdA_input,:);
    end
end
if exist('encdB_input','var')
    encdB = dig_in(encdB_input,:); 
end
if exist('mouse_input','var')
    if used_var(6) == 1
        mouse = dig_in(mouse_input,:);
    else
        mouse = adc(mouse_input,:);
    end
end

clear dig_in adc dig_in_raw 

% Extract event information (taken from intan2matlab_ver8.m)
disp('find epoc times')
disp('search var and find trans points')
st=find(epoc(1,:),1,'first');
en=find(epoc(1,:),1,'last');
counter=0;      
for i=st:1:en;  
    if counter==0
        if epoc(1,i)==1 && epoc(1,i-1)==0
            counter=counter+1;           
            fprintf('The number of tranisitons = %d\n',counter);    
            trans_timestamps2(counter,1)=i;
            trans_timestamps2(counter,2)=time(i);
            x=epoc(i:end);
            trans_timestamps2(counter,3)=find(x==0,1,'first')+st;
            trans_timestamps2(counter,4)=time(trans_timestamps2(counter,3));
        end
    elseif counter>0
        if epoc(1,i)==1 && epoc(1,i-1)==0
            counter=counter+1;           
            fprintf('The number of tranisitons = %d\n',counter);    
            trans_timestamps2(counter,1)=i;
            trans_timestamps2(counter,2)=time(i);
            x=epoc(i:end);
            if length(find(x==1))==length(x)
                disp('no tranisition')
                trans_timestamps2(counter,:)=[];
            else
            trans_timestamps2(counter,3)=find(x==0,1,'first')+i;
            trans_timestamps2(counter,4)=time(trans_timestamps2(counter,3));
            end
        end
    end
end
clear en st x i counter trials ans

trials(:,1)=trans_timestamps2(:,2);
trials(:,2)=trans_timestamps2(:,4);  

LN              = length(time);
div             = amp_sr/1000;
zx              = 1:div:LN;
izx             = floor(zx);
time_index      = time(izx);    % downsample from 20000 hz to 1000 hz

for i=1:size(trials,1)
    field_trials (i,1) = find(time_index>=trials(i,1),1,'first');   % in 1000 hz samples
    field_trials (i,2) = find(time_index>=trials(i,2),1,'first');
end 

[re]= getSyncTimesRevCorr_AC(photo',(1/amp_sr)); %flip vector      


save('data.mat', 'trials','field_trials','time_index','amp_sr',...
    var_names{find(used_vars>0)},'re')

return
