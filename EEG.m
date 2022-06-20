
%
% EEG.m
%
% Front-end for sleep analysis programs
%

%[~, hostname] = system('hostname');
%if hostname(length(hostname)) == 10
    % Note: last character of hostname may be ASCII #10, i.e. a newline.
    % Truncate this so strcmps will work.
%    hostname = hostname(1:length(hostname) - 1);
%end
%if strcmp(hostname, 'TomJhou-Office')
%    DropBoxDir = 'C:\Users\TomJhou\Dropbox\';
%elseif strcmp(hostname, 'NeuroSysLLC')
%    DropBoxDir = 'C:\Users\TomJhou\Dropbox\';
%end
% cd(DropBoxDir);
% if strcmp(hostname, 'TJ64-PC')
%     cd(DropBoxDir);% 'D:\Dropbox\Jhou lab shared items'; % folder at MUSC
% %    DropBoxDir = 'C:\Users\TJ64\Documents\My Dropbox\';
% elseif strcmp(hostname, 'TomJhou-Office')
%     cd 'E:\_projects\_LHb raw data\'; % folder at MUSC
% %    DropBoxDir = 'E:\My Dropbox\';
% end

global thresholdEEG;
global thresholdEMG;
global SleepState;
global SleepStateHourly;
global EEG_histogram_bins;
global EEG_histogram_counts;


fprintf('\n\nWhat would you like to do?\n');
fprintf('(1): Read EEG data from PLX file into memory, then score.\n');
fprintf('(2): Load EEG data from MAT file into memory, then score. (*.mat).\n');
fprintf('(3): Score EEG data in memory (must run options #1 or #2 first).\n');
fprintf('(4): Save EEG data from memory to *.mat file on disk. Scoring will be saved also.\n');
fprintf('(5): Clear previous scoring.\n');
fprintf('Simply press enter to exit\n');
answer = input('Select option: ');

clear skipprompt;       % This ensures that scripts will prompt user for folder

scoring = 0;
saving = 0;

if isempty(answer)

    return;

elseif answer == 1 || answer == 3

    if answer == 1
        [EEG_data, EMG_data, adfreq, DateTime, OpenedFileName] = EEG_read;
        centralEpoch = 1;
        CreateGlobalSleepVars;
    end
    if exist('EEG_data', 'var')
        if ~isempty(EEG_data)
            scoring = 1;
        end
    end

    if ~scoring
        fprintf('\nNo data found.\n\n');
    end
    
elseif answer == 2

    clear EEG_data EMG_data adfreq;     % Clear old variables so we will know if uiload actually loaded anything.
    centralEpoch = 1;
    OpenedFileName = '';
    CreateGlobalSleepVars;

    uiload;

    if exist('EEG_data', 'var')
        if ~isempty(EEG_data)
            scoring = 1;
        end
    end

elseif answer == 4

    saving = 1;
    
elseif answer == 5

    % Clear old scoring
    clear global

else
    
    fprintf('Invalid option\n');
    return;
    
end

if scoring
    
    if ~exist('centralEpoch', 'var')
        centralEpoch = 1;
    end
    
    [SleepState, SleepStateHourly, thresholdEEG, thresholdEMG, EEG_histogram_bins, EEG_histogram_counts, centralEpoch] = EEG_score(OpenedFileName, EEG_data, EMG_data, adfreq, SleepState, thresholdEEG, thresholdEMG, EEG_histogram_bins, EEG_histogram_counts, centralEpoch);
    fprintf('\n\nDo you want to save scored data (this will take about 10-15 seconds for a 24-hour recording)?\n');
    answer = input('Type "y" or "n": ', 's'); % Second argument must be 's', to indicate string argument that is not evaluated as an expression
    if answer == 'y' || answer == 'Y'
        saving = 1;
    end

end

if saving
    defaultFileName = strrep(OpenedFileName, '.plx', '.mat');
    if ~exist('centralEpoch', 'var')
        centralEpoch = 1;
    end
    uisave({'EEG_data', 'EMG_data', 'adfreq', 'DateTime', 'EEG_histogram_bins', 'EEG_histogram_counts', 'OpenedFileName', 'SleepState', 'SleepStateHourly', 'thresholdEEG', 'thresholdEMG', 'centralEpoch'}, defaultFileName);
%    fprintf('Finished saving.\n');
end

