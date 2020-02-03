 LFP=TETFP09(:,1);
time = TETFP09(:,2);

succ = SUCC(EVT04,EVT01);
[succrate] = SUCCSESSION(EVT04,EVT01);

for i = 1:1:length(Start)
  
    
    num=[];
    num = find(EVT04>Start(i)&EVT04<Stop(i));
    
    if succrate(i) >= 0.6
        
    Trialslab(num) = ones(length(num),1) + 1;
    
    else
        
          Trialslab(num) = ones(length(num),1);
    
    end
end


for i = 1:1:length(find(succ==1))
    tnum = find(succ==1);
    temp = EVT04(tnum(i));
    tempLFP = find(time > temp - 1.5 & time < temp+1.5);
     succLFP(1:2999,i) = LFP(tempLFP(1:2999));
end

for i = 1:1:length(find(succ==0))
     tnum = find(succ==0);
     temp = EVT04(tnum(i));
    tempLFP = find(time > temp - 1.5 & time < temp+1.5);
     falseLFP(1:2999,i) = LFP(tempLFP(1:2999));
     
end

tLFP{1} = succLFP;
tLFP{2} = falseLFP;

clear falseLFP
clear succLFP
% Create sine wave

duration = 3.0;        % total epoch length in seconds
srate    = 1000;       % Sampling rate in Hz
amp      = 1;          % peak amplitude of sine wave in uVolts      
freq     = 40;         % frequency of sine wave in Hz

nTrials = 81;       % number of trials
% perform time-frequency analysis
winsize_sec = 0.256;     % window size in seconds

winsize  = winsize_sec * srate;     % window size in points
padratio = 2;

[ersp,itc,powbase,times,freqs,erspboot,itcboot] = modified_newtimef(...
                tLFP,...                        % The EEG data for this channel
                2999, ...                       % Number of samples per epoch 
                [0 2999], ...    % Epoch start & end time in miliseconds 
                srate,...                       % Samples per second               
                0,...                           % varwin. Value ignored when 'cycles' parameter is provided
                'cycles',[3 0.5],...            % [0] = FFT; [3 0.5] = wavelets with number of cycles increasing with frequency.
                'winsize',winsize,...           % window size in samples. 
                'padratio',padratio, ...        % zero padding. 1=no padding, 2=double total length by zero padding
                'baseline',duration/2*1000,...  % endpoint (in ms) of baseline interval. NaN = no baseline;duration/2*1000
                'scale','log', ...              % results in uV^2/Hz rather than dB 
                'freqs',[0 120], ...            % range of frequencies to report
                'plotphase','on',...           % Do not plot ITC phase
                'plotersp','on',...            % Do not plot ERSP
                'plotitc','on');               % Do not plot ITC 
              


% Integrate powbase over frequency 
% NOTES: 
%   ASSUMES linear frequency spacing.
%   ASSUMES sine wave frequency is well within the range of frequencies
%   extracted.
% 
% FrqInterval = freqs(2) - freqs(1);              % Hz per frequency bin
% 
% TotalPower = sum(powbase * FrqInterval);
% 
% fprintf('\nTotal Power in Baseline Interval averaged across trials = %f\n',TotalPower);
% 
% 
% % Dislay the power spectral density of the original signal
% 
% figure();
% 
% plot(freqs,powbase);
% title('Power Spectral Density');
% xlabel('Frequency (Hz)');
% ylabel('Power (uV^2/Hz)');
% ylim([0 0.1]);
