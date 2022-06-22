
%
% Returns two column vectors
% 
%

function [EEG_data, EMG_data, adfreq1, DateTime, OpenedFileName] = EEG_read()

    EEG_CHANNEL_LEFT = 2;  % This is ignored.
    EEG_CHANNEL_RIGHT = 4; % This is the highest channel to start looking. Will search descending until data is found.
    EMG_CHANNEL = 6;       % Will also accept chan 5

    % OpenedFileName = '2016-1203-1608_Box1_Ephys.plx';


    [FileName2,PathName2,~] = uigetfile('*.plx', 'MultiSelect','on');
    OpenedFileName = [PathName2, FileName2];

    EEG_data = [];
    EMG_data = [];
    adfreq1 = 0;
    DateTime = [];

    if FileName2 == 0
        % User did not select a file
        fprintf('No file selected\n');
        return
    end
    
    fprintf('Selected file %s.\n', OpenedFileName);
    
    [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(OpenedFileName);


    % Normally IDs are sequential (1, 2, 3 ...)
    % but there could be non-consecutive chans, or map could start at 0 instead of 1.
    % ID number is used to retrieve raw data.
    [n, channel_ID] = plx_ad_chanmap(OpenedFileName);

    if channel_ID(1) == 0
        % Channel map starts with channel 0 (not used), so we have to add one to
        % channel ID to get map position.
    else
        offset = 0;
    end
    
    % get counts of # of events for spike, waveform, event, and continuous data channels
    [~, ~, ~, continuous_counts] = plx_info(OpenedFileName,1);

    for x = EEG_CHANNEL_RIGHT:-1:1
        % Read all data from EEG and EMG channels
        numPoints = continuous_counts(x + offset);
        if numPoints == 0
            continue
        end
        [adfreq1, n1, EEG_data] = plx_ad_span_v(OpenedFileName, x, 1, numPoints);
        fprintf('Detected EEG signal in channel %d\n', x);
        break
    end
    
    % Convert from mV to uV
    EEG_data = EEG_data * 1000;
    haveEMG = false;

    for x = (EMG_CHANNEL-1):EMG_CHANNEL
        % Read all data from EEG and EMG channels
        n2 = continuous_counts(x + offset);
        if n2 == 0
            continue
        end
        [adfreq2, n2, EMG_data] = plx_ad_span_v(OpenedFileName, x, 1, n2);
        haveEMG = true;
        break
    end

    if haveEMG
        % Find the shorter of the two waveforms (generally they should match
        % exactly, so we're being ultra cautious here)
        numPoints = min(numPoints, n2);

        % Convert from mV to uV
        EMG_data = EMG_data * 1000;
    end
    
    % Close all open files
    plx_close('');
    
    err = 0;
    
    if haveEMG
        if adfreq1 ~= adfreq2
            fprintf('Error: different sample frequencies %d and %d for EEG and EMG.\n', adfreq1, adfreq2);
            err = 1;
        end

        if n1 ~= n2 || n1 ~= numPoints
            fprintf('Error: different sample counts %d and %d for EEG and EMG.\n', n1, n2);
            err = 1;
        end
    else
        fprintf('Warning: no EMG data found.\n');
    end    
  
    fprintf('\nSuccessfully read EEG/EMG from file \"%s\"\n', OpenedFileName);
    fprintf('Data has %d samples, sample rate is %d Hz.\n', n1, adfreq1);

    if adfreq1 >= 250
        factor = 2;
        if adfreq1 >= 500
            factor = 4;
            if (adfreq1 >= 1000)
                factor = 8;
            end
        end
        
        fprintf('Sample rate %d is higher than necessary. Downsample by %dx? ', adfreq1, factor);
        answer = input('Enter y or n: ', 's');
        
        if answer == 'y' || answer == 'Y'
            if rem(n1, factor) > 0
                % Make sure array sizes are integer multiples of
                % downsampling factor
                n1 = n1 - rem(n1, factor);
                EEG_data = EEG_data(1:n1);
                if haveEMG
                    n2 = n2 - rem(n2, factor);
                    EMG_data = EMG_data(1:n2);
                end
            end
            
            % Downsample
            EEG_data = mean(reshape(EEG_data, factor, n1/factor))';
            if haveEMG
                EMG_data = mean(reshape(EMG_data, factor, n1/factor))';
            end
        end
        
        adfreq1 = adfreq1 / factor;
        
        fprintf('Downsampled EEG/EMG from %dHz to %dHz\n', adfreq1 * factor, adfreq1);
    end
    
end


function  [tscounts, wfcounts, evcounts, contcounts] = plx_info(filename, fullread)
% plx_info(filename, fullread) -- read and display .plx or .pl2 file info
%
% [tscounts, wfcounts, evcounts, contcounts] = plx_info(filename, fullread)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   fullread - if 0, reads only the file header
%              if 1, reads the entire file
%               for .pl2 files, this parameter is ignored
%
% OUTPUT:
%   tscounts - 2-dimensional array of timestamp counts for each unit
%      tscounts(i, j) is the number of timestamps for channel j-1, unit i
%                                (see comment below)
%   wfcounts - 2-dimensional array of waveform counts for each unit
%     wfcounts(i, j) is the number of waveforms for channel j-1, unit i
%                                (see comment below)
%   evcounts - 1-dimensional array of external event counts
%     evcounts(i) is the number of events for event channel i
%
%   contcounts - 1-dimensional array of sample counts for continuous channels
%     contcounts(i) is the number of continuous for slow channel i-1
%
% Note that for tscounts, wfcounts, the channel indices i,j are off by one.
% That is, for channels, the count for channel n is at index n+1.
% 
% Unit index 1 is unsorted, 2 = unit a, 3 = unit b, etc
% The dimensions of the tscounts and wfcounts arrays are
%   (NChan+1) x (MaxUnits+1)
% where NChan is the number of spike channel headers in the plx file, and
% MaxUnits is 4 if fullread is 0, or 26 if fullread is 1. This is because
% the header of a .plx file can only accomodate 4 units, but doing a
% fullread on the file may show that there are actually up to 26 units
% present in the file. Likewise, NChan will have a maximum of 128 channels
% if fullread is 0.
% The dimension of the evcounts and contcounts arrays is the number of event
% and continuous (slow) channels. 
% The counts for slow channel 0 is at contcounts(1)

tscounts = [];
wfcounts = [];
evcounts = [];
contcounts = [];

if nargin ~= 2
    error 'expected 2 input arguments';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    numSpikeChannels = numel(pl2.SpikeChannels);
    % pl2 files support up to 256 units, but we limit to 
    % 26 sorted plus 1 unsorted to be compatible with plx_info
    tscounts = zeros(27,numSpikeChannels+1);
    for i=1:numSpikeChannels
        tscounts(:,i+1) = pl2.SpikeChannels{i}.UnitCounts(1:27);
    end
    wfcounts = tscounts;

    numAnalogChannels = numel(pl2.AnalogChannels);
    contcounts = zeros(1,numAnalogChannels);
    for i=1:numAnalogChannels
        contcounts(1,i) = pl2.AnalogChannels{i}.NumValues;
    end

    evcounts = pl2.EventCounts;
    return
end

[tscounts, wfcounts, evcounts, contcounts] = mexPlex(4, filename, fullread);

end


function  [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(filename)
% plx_information(filename) -- read extended header infromation from a .plx file
%
% [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
% OpenedFileName    - returns the filename (useful if empty string is passed as filename)
% Version -  version code of the plx file format
% Freq -  timestamp frequency for waveform digitization
% Comment - user-entered comment
% Trodalness - 0,1 = single electrode, 2 = stereotrode, 4 = tetrode
% Number of Points Per Wave - number of samples in a spike waveform
% Pre Threshold Points - the sample where the threshold was crossed
% SpikePeakV - peak voltage in mV of the final spike A/D converter
% SpikeADResBits - resolution of the spike A/D converter (usually 12 bits)
% SlowPeakV - peak voltage of mV of the final analog A/D converter
% SlowADResBits - resolution of the analog A/D converter (usually 12 bits)
% Duration - the duration of the file in seconds
% DateTime - date and time string for the file

[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = mexPlex(13,filename);

end


function [adfreq, n, ad] = plx_ad_span_v(filename, channel, startCount, endCount)
% plx_ad_span_v(filename, channel): Read a span of a/d data from a .plx file
%
% [adfreq, n, ad] = plx_ad_span_v(filename, channel, startCount, endCount)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   startCount - index of first sample to fetch
%   endCount - index of last sample to fetch
%   channel - 0 - based channel number
%
% OUTPUT:
%   adfreq - digitization frequency for this channel
%   n - total number of data points 
%   ad - array of a/d values converted to mV

if nargin < 4
    error 'Expected 4 input arguments';
end
if (isempty(filename))
   [fname, pathname] = uigetfile('*.plx', 'Select a Plexon .plx file');
   if isequal(fname,0)
     error 'No file was selected'
   end
   filename = fullfile(pathname, fname);
end

[adfreq, n, ad] = mexPlex(18, filename, channel, startCount, endCount);

end


function  [n,adchans] = plx_ad_chanmap(filename)
% plx_ad_chanmap(filename) -- return map of raw continuous channel numbers for each channel
%
% [n,adchans] = plx_ad_chanmap(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   n - number of continuous channels
%   adchans - 1 x n array of continuous channel numbers
%
% Normally, there is one channel entry in the .plx for for each raw continuous channel,
% so the mapping is trivial adchans[i] = i-1 (because continuous channels start at 0).
% However, for certain .plx files saved in some ways from OFS (notably after
% loading data files from other vendors), the mapping can be more complex.
% E.g. there may be only 2 non-empty channels in a .plx file, but those channels
% correspond to raw channel numbers 7 and 34. So in this case NChans = 2, 
% and adchans[1] = 7, adchans[2] = 34.
% The plx_ routines that return arrays always return arrays of size NChans. However,
% routines that take channels numbers as arguments always expect the raw  
% channel number.  So in the above example, to get the data from  
% the second channel, use
%   [adfreq, n, ts, fn, ad] = plx_ad(filename, adchans[2])

if nargin < 1
    error 'Expected 1 input argument';
end
if (isempty(filename))
   [fname, pathname] = uigetfile('*.plx', 'Select a Plexon .plx file');
   if isequal(fname,0)
     error 'No file was selected'
   end
   filename = fullfile(pathname, fname);
end

[n,adchans] = mexPlex(27,filename);

end


function [n] = plx_close(filename)
% plx_close(filename): close the .plx file
%
% [n] = plx_close(filename)
%
% INPUT:
%   filename - if empty string, will close any open files
%
% OUTPUT:
%   n - always 0

[n] = mexPlex(22, filename);

end

