
%
% Reads PLX file, and returns EEG_data variable of size nPoints x nChans,
% where nChans is the number of EEG channels detected, and nPoints is the
% length of the channels (all must be of same length).
%

function [EEG_data, adfreq, DateTime, OpenedFileName] = EEG_read_multichannel()

    [FileName2,PathName2,FilterIndex2] = uigetfile('*.plx', 'MultiSelect','on');
    OpenedFileName = [PathName2, FileName2];

    EEG_data = [];
    adfreq = -1;
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
        % Channel map starts with channel 0 (not used), so we have to add
        % one to channel ID to get position in map. Conversely, we have to
        % subtract one from map position to get channel ID.
        offset = 1;
    else
        offset = 0;
    end
    
    % get counts of # of events for spike, waveform, event, and continuous data channels
    [~, ~, ~, continuous_counts] = plx_info(OpenedFileName,1);

    numPoints = -1;

    for x = 1:n % Loop variable is map position, with is usually = chan ID but not always!
        numPoints1 = continuous_counts(x);
        if numPoints1 == 0
            continue
        end

        [adfreq1, n1, EEG_data1] = plx_ad_span_v(OpenedFileName, x - offset, 1, numPoints1);

        if adfreq < 0
            adfreq = adfreq1;
        end

        if numPoints < 0
            numPoints = numPoints1;
        end

        if adfreq1 ~= adfreq
            fprintf('Warning: different channels has unexpected sample rate (%d instead of %d).\n', adfreq1, adfreq);
        end

        if n1 ~= numPoints
            fprintf('Warning: channel %d has unexpected length (%d instead of %d), will skip.\n', x, numPoints1, numPoints);
            continue;
        end

        EEG_data = [EEG_data, EEG_data1];
    end
    
    % Close all open files
    plx_close('');
        
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

