
function [SleepState, SleepStateHourly, thresholdEEG, thresholdEMG, EEG_histogram_bins, EEG_histogram_counts, centralEpoch] = EEG_score(OpenedFileName, EEG_data, EMG_data, adfreq, SleepState, thresholdEEG, thresholdEMG, EEG_histogram_bins, EEG_histogram_counts, centralEpoch)

% This function takes EEG and EMG time series, attempts to decide
% wake/NREM/REM for every 12-second epoch, and displays EEG/EMG and
% sleep status (W/N/R) on screen in a user-friendly graphical user
% interface. User can manually change auto scoring for each epoch.
%
% Input arguments:
%   isFirstScore    True if this data file is freshly loaded from PLX, and
%                   has never been scored before. Causes outlier removal
%                   and 
%   EEG_data        Vector of EEG voltages at specified sampling rate
%   EMG_data        Vector of EMG voltages at specified sampling rate
%   adfreq          Sampling rate in Hz. Must be the same for both EEG and EMG
%
% Output arguments:
%   SleepStateMatrix      matrix with hourly wake, NREM, and REM percentages
%   EEG_histogram_bins    delta power histogram bins
%   EEG_histogram_counts  delta power histogram counts
%

    % change this to 1 if you want to autoscore.
    autoscore = 0;
    
    BIN_SECONDS = 12;
    DELTA_SCALE_FACTOR = 20;    % Delta powers are about an order of magnitude lower than voltage amplitudes

    STATE_WAKE = 0;
    STATE_NREM = 1;
    STATE_REM = 2;
    STATE_A = 3;
    STATE_B = 4;
    STATE_C = 5;
    
    
    global    strings;

    strings = { 'wake', 'NREM', 'REM', 'A', 'B', 'C' };

    fprintf('\n\n*************************************\n');
    fprintf('SLEEP SCORING INSTRUCTIONS:\n\n');
    fprintf('Data is shown in 12 second epochs. Central epoch (outlined in green) is the "active" epoch,\n     surrounding epochs are shown for context.\n');
    fprintf('Left window shows FFT for active epoch. Orange lines are top and bottom quartiles of amplitudes.\n');
    fprintf('Use left/right arrow keys to scroll through epochs.\n');
    fprintf('Type W, N, and R to change state to wake, non-REM, and REM.\n');
    fprintf('Use space bar to zoom in/out of active epoch.\n');
    fprintf('Type "t" to adjust thresholds for automatic scoring algorithm.\n\n');
    fprintf('Type "x" to exit.\n\n');
    fprintf('IMPORTANT: Do NOT try to exit by killing the scoring window - If you do, then typing "x" will no longer work,\n    and you will have to exit with Control-C, and will be unable to save your work!!!!\n\n');
    fprintf('*************************************\n');
    
    
    % How many epochs before and after to plot in scrolling window
    EPOCH_PADDING = 5;

    haveEMG = length(EMG_data) > 0;
    
    numPoints = length(EEG_data);
    if haveEMG && numPoints ~= length(EMG_data)
        fprintf('Error: EEG and EMG data are different lengths: %d vs %d\n', numPoints, length(EMG_data));
        return;
    end

    [EEG_data, EMG_data] = RemoveOutliers(EEG_data, EMG_data, adfreq);
    
    % Replace NaN with zero, or else all math will barf.
    EEG_data(isnan(EEG_data)) = 0;
    EMG_data(isnan(EMG_data)) = 0;

    [EEG_data, EMG_data] = FilterEEG(EEG_data, EMG_data, adfreq);

    % Number of time points in each 12-second bin
    if ceil(adfreq) == floor(adfreq)
        numSamplesPerEpoch = adfreq * BIN_SECONDS;
    else
        % Non-integer sample frequency, requires slight adjustment to BIN
        oldBin = BIN_SECONDS;
        numSamplesPerEpoch = round(adfreq * BIN_SECONDS);
        BIN_SECONDS = numSamplesPerEpoch / adfreq;
        fprintf('Warning: Sample frequency is not an integer, so Epoch size has been adjusted slightly from %f to %.5f seconds\n', oldBin, BIN_SECONDS);
    end
    
    % How many 12-second bins are there?
    numEpochs = floor(length(EEG_data) / numSamplesPerEpoch);

    % Slightly reduce numPoints so it is exact multiple of binsize
    numPoints = numEpochs * numSamplesPerEpoch;

    % Truncate ends of array so it is exact multiple of binsize
    EEG_data = EEG_data(1:numPoints);
    if haveEMG
        EMG_data = EMG_data(1:numPoints);
    end
    
    % Organize vector into matrix.
    % # of rows is # of samples in each 12-second bin
    % # of columns is # of 12-second bins in recording
    EEG_data = reshape(EEG_data, [numSamplesPerEpoch numEpochs]);
    if haveEMG
        EMG_data = reshape(EMG_data, [numSamplesPerEpoch numEpochs]);
    end
    
    % Subtract mean of each column from EVERY element.
    % This will remove DC components from each bin individually.
%    EEG_data = bsxfun(@minus,EEG_data,mean(EEG_data));
%    EMG_data = bsxfun(@minus,EMG_data,mean(EMG_data));

    [delta_power, theta_power, EMG_amplitude, EEG_spectrum, FFT_frequencies] = calcDeltaTheta(EEG_data, EMG_data, numSamplesPerEpoch, adfreq);
    
    if isempty(thresholdEEG)
        thresholdEEG = median(delta_power); % findThreshold(delta_power, rangeEEG);
    end
    
    median_delta = median(delta_power);
    [EEG_histogram_bins, EEG_histogram_counts] = plotHistAmplitude(delta_power, median_delta, median_delta * 5);
    
    title('Wake/NREM threshold in EEG delta power histogram');
    xlabel('EEG delta power')

    if isempty(thresholdEMG)
        if haveEMG
    %        rangeEMG = 5 * median(EMG_amplitude);
            thresholdEMG = median(EMG_amplitude); % findThreshold(EMG_amplitude, rangeEMG);
            plotHistAmplitude(EMG_amplitude, thresholdEMG, thresholdEMG * 5);
            title('EMG threshold');
            xlabel('EMG amplitude')
        else
            thresholdEMG = 0;
        end
    end
    
    MAX_EEG_AMPLITUDE = max(thresholdEEG) * 5 * DELTA_SCALE_FACTOR;
    if haveEMG
        MAX_EMG_AMPLITUDE = max(thresholdEMG) * 10;
    end
    
    avg_spectrum = mean(EEG_spectrum, 2);
    spectrum_bottom_quartile = quantile(EEG_spectrum, 0.25, 2);
    spectrum_top_quartile = quantile(EEG_spectrum, 0.75, 2);

    % centralEpoch will increment or decrement as user scrolls back and forth. It
    % identifies which epoch will appear in center of screen.
%    centralEpoch = 1;
    
    % Create main figure, and put filename in titlebar
    h = figure('Name', OpenedFileName);

    h.WindowState = 'maximized';
    
    % Create main figure, and put filename in titlebar, then make fullscreen
%    figure('Name', OpenedFileName, 'units','normalized','outerposition',[0 0 1 1])
    
    % Use fancy new method of reading keyboard.
    set(h,'KeyPressFcn',@KeyPressCb);

    global val;
    global s;   
    
    % Determines width of scoring windows
    PLOT_W = 7;
    
    max_spectrum = max(max(EEG_spectrum));
    
    % Accelerator when you hold down page up/down key
    consec_page_keys = 0;
    
    lastPageKeyTime = 0;
    % Start timer that detects interval between page up/down keys
    tic;
    
    % Convert thresholds to vectors so that thresholds can be different for each epoch
    if length(thresholdEEG) == 1
        thresholdEEG = ones(1, numEpochs) * thresholdEEG;
        thresholdEMG = ones(1, numEpochs) * thresholdEMG;
    end
    
    if isempty(SleepState)
        if autoscore
            SleepState = DecideState(delta_power, EMG_amplitude, thresholdEEG, thresholdEMG);
        else
            SleepState = ones(length(delta_power), 1) * -1;
        end
    else
        if size(SleepState, 2) > 1
            SleepState = SleepState(:, 1);
        end
        
        % Create scatterplot figure
        h2 = figure('Name', 'Delta/Theta scatterplot');
        hold on;
        
        WakeEpochs = find(SleepState == STATE_WAKE);
        NREM_epochs = find(SleepState == STATE_NREM);
        REM_epochs = find(SleepState == STATE_REM);
        
        if length(WakeEpochs) > 0
            scatter3(EMG_amplitude(WakeEpochs), delta_power(WakeEpochs), theta_power(WakeEpochs), '.', 'k');
        end
        if length(NREM_epochs) > 0
            scatter3(EMG_amplitude(NREM_epochs), delta_power(NREM_epochs), theta_power(NREM_epochs), '.', 'b');
        end
        if length(REM_epochs) > 0
            scatter3(EMG_amplitude(REM_epochs), delta_power(REM_epochs), theta_power(REM_epochs), '.', 'r');
        end

        xlabel('Amplitude')
        ylabel('Delta power')
        zlabel('Theta power')
        title('Delta/Theta correlations for wake, NREM, REM (black, blue, red)')

        % Revert back to original figure
        figure(h);
    end


    % Write instructions
    subplot(4, PLOT_W, 1);
    title("Instructions: IMPORTANT, MUST READ!!!")
    axis off
    text(0,0.9, "Type 'x' to exit (do NOT close window)");
    text(0,0.8, "Type 'n', 'w', 'r' to score NREM, wake, or REM");
    text(0,0.7, "Home/End  jumps to beginning/end");
    text(0,0.6, "Left/Right moves forward/back one epoch");
    text(0,0.5, "Up/Down   changes vertical scaling");


    % Now generate main window, and loop until user quits.
    while 1

        % x values and tick labels have to be calculated every time, in
        % case we have zoomed in or out.
        
        % We are showing a range of epochs around the central epoch. Build
        % the data from each item in the range by "stacking" the data
        % vertically on top of each other in the following "stacked"
        % arrays:
        EEG_stacked = [];
        EMG_stacked = [];
        delta_stacked = [];
        theta_stacked = [];
        EMG_amplitude_stacked = [];
        threshold_EEG_stacked = [];
        threshold_EMG_stacked = [];
        
        % "Stack" several epochs together for display purposes
        for offset = -EPOCH_PADDING:EPOCH_PADDING
            epoch = centralEpoch + offset;
            if epoch < 1 || epoch > numEpochs
                % Insert blank data if we are past the very end or before the
                % very beginning of file
                blankData = zeros(numSamplesPerEpoch,1);
                EEG_stacked = [EEG_stacked; blankData];
                threshold_EEG_stacked = [threshold_EEG_stacked; blankData];
                if haveEMG
                    EMG_stacked = [EMG_stacked; blankData];
                    threshold_EMG_stacked = [threshold_EMG_stacked; blankData];
                    EMG_amplitude_stacked = [EMG_amplitude_stacked; blankData];
                end
                delta_stacked = [delta_stacked; blankData];
                theta_stacked = [theta_stacked; blankData];
            else
                EEG_stacked = [EEG_stacked; EEG_data(:, epoch)];
                threshold_EEG_stacked = [threshold_EEG_stacked; ones(numSamplesPerEpoch, 1) * thresholdEEG(epoch)];
                if haveEMG
                    EMG_stacked = [EMG_stacked; EMG_data(:, epoch)];
                    threshold_EMG_stacked = [threshold_EMG_stacked; ones(numSamplesPerEpoch, 1) * thresholdEMG(epoch)];
                    EMG_amplitude_stacked = [EMG_amplitude_stacked; ones(numSamplesPerEpoch, 1) * EMG_amplitude(epoch)];
                end
                delta_stacked = [delta_stacked; ones(numSamplesPerEpoch, 1) * delta_power(epoch)];
                theta_stacked = [theta_stacked; ones(numSamplesPerEpoch, 1) * theta_power(epoch)];
            end
        end
        
        % Plot FFT
        subplot(4,PLOT_W,1 + PLOT_W * 2);
        hold off;
        loglog(FFT_frequencies, smooth_copy(EEG_spectrum(:,centralEpoch)));
        hold on;
        loglog(FFT_frequencies, smooth_copy(spectrum_bottom_quartile), 'color', [1 .4 0]);
        loglog(FFT_frequencies, smooth_copy(spectrum_top_quartile), 'color', [1 .4 0]);

        % Limit the X-axis to 1- 0, and y-axis to 4 log units.
        axis([1 100 max_spectrum/10000 max_spectrum]);
        title('FFT');
        xlabel('Frequency (Hz)');

        % Convert sample # to seconds
        min_x = -numSamplesPerEpoch * EPOCH_PADDING;
        max_x = numSamplesPerEpoch * (EPOCH_PADDING + 1);
        x_values_seconds = (min_x:max_x - 1)/adfreq;

        % First row, EEG
        subplot(4, PLOT_W, [2:PLOT_W]);
        hold off;   % So we can erase previous plot
        EEG_PLOT_PANEL(BIN_SECONDS, adfreq, centralEpoch, numSamplesPerEpoch, numEpochs, ...
            x_values_seconds, SleepState, ...
            EEG_stacked, threshold_EEG_stacked, delta_stacked, theta_stacked, ...
            DELTA_SCALE_FACTOR, EPOCH_PADDING, MAX_EEG_AMPLITUDE);
        xlabel('Seconds');

        % Second row, EMG
        if haveEMG
            subplot(4, PLOT_W, [PLOT_W+2:PLOT_W*2]);
            hold off;
            EMG_PLOT_PANEL(BIN_SECONDS, adfreq, centralEpoch, numSamplesPerEpoch, numEpochs, ...
                x_values_seconds, SleepState, ...
                EMG_stacked, threshold_EMG_stacked, EMG_amplitude_stacked, EPOCH_PADDING, MAX_EMG_AMPLITUDE);

            xlabel('Seconds');
        end
        
        % Third row, EEG
        subplot(4, PLOT_W, [PLOT_W*2 + 2:PLOT_W*3]);
        hold off;   % So we can erase previous plot
        firstSample = EPOCH_PADDING * numSamplesPerEpoch;
        lastSample = firstSample + numSamplesPerEpoch - 1;

        x_values_seconds = (0:numSamplesPerEpoch - 1)/adfreq;

        EEG_PLOT_PANEL(BIN_SECONDS, adfreq, centralEpoch, numSamplesPerEpoch, numEpochs, ...
            x_values_seconds, SleepState, ...
            EEG_stacked(firstSample:lastSample), threshold_EEG_stacked(firstSample:lastSample), delta_stacked(firstSample:lastSample), ...
            theta_stacked(firstSample:lastSample), DELTA_SCALE_FACTOR, 0, MAX_EEG_AMPLITUDE);

        xlabel('Seconds');

        % Final (fourth) row, EMG
        if haveEMG
            subplot(4, PLOT_W, [PLOT_W*3 + 2:PLOT_W*4]);

            hold off;

            EMG_PLOT_PANEL(BIN_SECONDS, adfreq, centralEpoch, numSamplesPerEpoch, numEpochs, ...
                x_values_seconds, SleepState, ...
                EMG_stacked(firstSample:lastSample), threshold_EMG_stacked(firstSample:lastSample), EMG_amplitude_stacked(firstSample:lastSample), ...
                0, MAX_EMG_AMPLITUDE);

            xlabel('Seconds');
        end
        
        % Wait for uiresume, which is in the keyboard callback function
        % below this one. It sets the values of 's' and 'val', which 
        % contain key code and ASCII value, respectively.
        uiwait;

        advanceOne = 0;

        if strcmp(s, 'leftarrow')
            % Left arrow
            centralEpoch = centralEpoch - 1;
            if centralEpoch < 1
                % Avoid scrolling before beginning
                centralEpoch = 1;
            end
        elseif strcmp(s, 'rightarrow')
            % Right arrow
            advanceOne = 1;
            
        elseif strcmp(s, 'uparrow')
            % Up arrow
            MAX_EEG_AMPLITUDE = MAX_EEG_AMPLITUDE / sqrt(2);
        elseif strcmp(s, 'downarrow')
            % Down arrow
            MAX_EEG_AMPLITUDE = MAX_EEG_AMPLITUDE * sqrt(2);
        elseif strcmp(s, '+') || strcmp(s, '=') % Need strcmp in case s is a string
            MAX_EMG_AMPLITUDE = MAX_EMG_AMPLITUDE / sqrt(2);
        elseif s == '-'
            MAX_EMG_AMPLITUDE = MAX_EMG_AMPLITUDE * sqrt(2);
        elseif s == 'x'
            % letter x = 24th letter in alphabet
            close all;
            break;
        elseif s == 's'
            msgbox('Please use "N" instead of "S" to score epoch as non-REM sleep');
        elseif s == 'n'
            SleepState(centralEpoch) = 1;
            advanceOne = 1;
        elseif s == 'w'
            SleepState(centralEpoch) = 0;
            advanceOne = 1;
        elseif s == 'r'
            SleepState(centralEpoch) = 2;
            advanceOne = 1;
        elseif s == 'a'
            SleepState(centralEpoch) = 3;
            advanceOne = 1;
        elseif s == 'b'
            SleepState(centralEpoch) = 4;
            advanceOne = 1;
        elseif s == 'c'
            SleepState(centralEpoch) = 5;
            advanceOne = 1;
        elseif strcmp(s, 'N') || strcmp(s, 'W') || strcmp(s, 'R') || strcmp(s, 'T') || strcmp(s, 'X')
            msgbox('Please don''t use upper case letters. Do you have caps lock on?');
        elseif strcmp(s, 'pagedown')==1

            if toc < 0.5
                consec_page_keys = consec_page_keys+1;
            else
                % If more than 0.5 seconds have passed
                % since last page up/down key, then we reset
                % counter
                consec_page_keys = 0;
            end
            tic;

            if consec_page_keys > 50
                mult = 30;
            elseif consec_page_keys > 25
                mult = 10;
            else
                mult = 1;
            end
            
            centralEpoch = centralEpoch + (EPOCH_PADDING*2 + 1) * mult;
            if centralEpoch > numEpochs
                centralEpoch = numEpochs;
            end
            % Need a continue so we don't reset the consecutive page key
            % counter
            continue;
        elseif strcmp(s, 'pageup')==1
            
            if toc < 0.5
                consec_page_keys = consec_page_keys+1;
            else
                % If more than 0.5 seconds have passed
                % since last page up/down key, then we reset
                % counter
                consec_page_keys = 0;
            end
            tic;

            if consec_page_keys > 50
                mult = 30;
            elseif consec_page_keys > 25
                mult = 10;
            else
                mult = 1;
            end
            
            centralEpoch = centralEpoch - (EPOCH_PADDING*2 + 1) * mult;
            if centralEpoch < 1
                centralEpoch = 1;
            end
            % Need a continue so we don't reset the consecutive page key
            % counter
            continue;
        
        elseif strcmp(s, 'home')

            centralEpoch = 1;
            
        elseif strcmp(s, 'end')

            centralEpoch = numEpochs;
        
        elseif s == 't'

            [v,ok] = listdlg('PromptString','What threshold to change?',...
                'SelectionMode','single', 'ListString',{'EEG threshold for entire file','EMG threshold for entire file','EEG threshold starting with current epoch','EMG threshold starting with current epoch'},...
                'Name','Select option','ListSize',[500 130]);
            
            if ok
                if v == 1 || v == 3
                    prompt={sprintf('EEG threshold for current epoch is %0.2f, enter new value:',thresholdEEG(centralEpoch) * DELTA_SCALE_FACTOR)};
                    threshold_title='Input EEG threshold'; 
                else
                    prompt={sprintf('EMG threshold for current epoch is %0.2f, enter new value:',thresholdEMG(centralEpoch))};
                    threshold_title='Input EMG threshold'; 
                end
                
                num_lines = 1;
    %             defaultans={num2str(thresholdEEG),num2str(thresholdEMG)};
                answer=inputdlg(prompt,threshold_title, num_lines);
            
                if isempty(answer)
                    continue;
                end

                if v == 1 || v == 3
                    thresholdEEGtmp = (str2double(answer)) / DELTA_SCALE_FACTOR;
                    if v == 1
                        thresholdEEG(1:end) = thresholdEEGtmp;
                    else
                        thresholdEEG(centralEpoch:end) = thresholdEEGtmp;
                    end
                elseif v == 2 || v == 4
                    thresholdEMGtmp = str2double(answer);
                    if v == 2
                        thresholdEMG(1:end) = thresholdEMGtmp;
                    else
                        thresholdEMG(centralEpoch:end) = thresholdEMGtmp;
                    end
                end
            
                SleepStateTmp = DecideState(delta_power, EMG_amplitude, thresholdEEG, thresholdEMG);

                % Update auto-scoring only for epochs including and after the
                % current one. Previous epochs will remain the same.
                if v == 1 || v == 2
                    SleepState = SleepStateTmp;
                else
                    SleepState(centralEpoch:end) = SleepStateTmp(centralEpoch:end);
                end
            end
            
        end % closes large if elseif statement
        
        if advanceOne
            % Advance one epoch. Make sure we don't scroll past end of
            % file.
            centralEpoch = centralEpoch + 1;
            if centralEpoch > numEpochs
                % Avoid scrolling past end
                centralEpoch = numEpochs;
            end
            
        end

        % Reset the consecutive page up/down key counter so we won't
        % accelerate
        consec_page_keys = 0;

    end %closes the while loop
    
     %put percentage calc here
    nhrs = floor(length(SleepState)/(5*60));
    npts = 300;
    SleepStatePerHour = reshape(SleepState(1:(nhrs*npts)),[npts,nhrs])';
    
    for p = 1:nhrs % size(SleepStatePerHour,1)
        percSleep(p) = length(find(SleepStatePerHour(p,:)==1))/npts;
        percWake(p) = length(find(SleepStatePerHour(p,:)==0))/npts;
        percREM(p) = length(find(SleepStatePerHour(p,:)==2))/npts;
        percA(p) = length(find(SleepStatePerHour(p,:)==3))/npts;
        percB(p) = length(find(SleepStatePerHour(p,:)==4))/npts;
        percC(p) = length(find(SleepStatePerHour(p,:)==5))/npts;
%        percSleepandREM(p) = length(find(SleepStatePerHour(p,:)==1 | SleepStatePerHour(p,:)==2))/npts;
    end
    
    if nhrs > 0
        SleepStateHourly = [percSleep' percWake' percREM' percA' percB' percC'].*100;
    else
        fprintf('Need at least one hour of sleep to calculate SleepStateMatrix.\n');
        SleepStateHourly = [];
    end
    
    %Figure, plot(SleepStateMatrix)
    %legend('Sleep', 'Wake', 'REM', 'Sleep+REM') 
end % closes EEG_score function


% Decide sleep state based on delta power, EMG amplitude, and thresholds
function SleepState = DecideState(delta_power, EMG_amplitude, thresholdEEG, thresholdEMG)

    SleepState = double(delta_power > thresholdEEG);
    PutativeWake = EMG_amplitude > thresholdEMG;

    SleepState(find(PutativeWake)) = 0;
    
    numEpochs = length(delta_power);

    for x = 1:numEpochs

        if length(EMG_amplitude) >= x
            if EMG_amplitude(x) > thresholdEMG(x)
                % EMG is high, we are awake
                SleepState(x) = 0;
                continue;
            end
        end

        if delta_power(x) > thresholdEEG(x)
            % EMG low, delta high, we are asleep
            SleepState(x) = 1;
        else
            % EMG was low, and so is delta power. Is this REM or
            % wake?
            if x == 1
                SleepState(x) = 0;
                continue;
            end

            if SleepState(x - 1) == 0
                % Previously awake, so we stay awake
                SleepState(x) = 0;
            else
                % Previously sleep or REM, so go to REM
                SleepState(x) = 2;
            end
        end
    end
end

function y = KeyPressCb(~,evnt)
%    fprintf('key pressed: %s\n',evnt.Key);
    global s;
    global val;
    global anyKey;
    
    if strcmp(evnt.Key,'rightarrow') || strcmp(evnt.Key, 'leftarrow') || strcmp(evnt.Key, 'uparrow') || strcmp(evnt.Key, 'downarrow')
        s = evnt.Key;
    elseif strcmp(evnt.Key, 'pageup') || strcmp(evnt.Key, 'pagedown')
        s = evnt.Key;
    elseif strcmp(evnt.Key,'space')
        s = evnt.Key;
    elseif strcmp(evnt.Key, 'home') || strcmp(evnt.Key, 'end')
        s = evnt.Key;
    else
        % Letter, symbol or other character
        val = -1;
        s = evnt.Character;  % Use evnt.Character instead of evnt.Key, as the latter does not account for Shift
        if strcmp(evnt.Key, 'shift') && isempty(evnt.Character)
            % Shift with no char pressed. Just skip.
        else
            uiresume
        end
        return;
    end
    
    uiresume
    
end

% Calculate delta and theta power
function [delta_power, theta_power, EMG_amplitude, EEG_spectrum, frequencies] = calcDeltaTheta(EEG_data, EMG_data, samplesPerEpoch, adfreq)

    % EEG and EMG data must be in matrix form.
    % # of rows is # of samples in each 12-second epoch
    % # of columns is # of 12-second epoch in recording
    
    DELTA_START = 1.0;
    DELTA_CUTOFF = 4.0;
    
    THETA_START = 6.0;
    THETA_CUTOFF = 10.0;

    % Calculate spectrum for each column (i.e. each 12-second epoch), remove phase information using ABS()
    EEG_spectrum = abs(fft(EEG_data)) / samplesPerEpoch;

    % Remove DC component
    EEG_spectrum(1,:) = 0;

    % Calculcate which bin corresponds to 4 Hz cutoff
    delta_freq_bin1 = floor(DELTA_START / adfreq * samplesPerEpoch);
    delta_freq_bin2 = floor(DELTA_CUTOFF / adfreq * samplesPerEpoch);

    % Calculcate which bin corresponds to 4 Hz to 8 Hz cutoff
    theta_freq_bin1 = floor(THETA_START / adfreq * samplesPerEpoch);
    theta_freq_bin2 = floor(THETA_CUTOFF / adfreq * samplesPerEpoch);

    % Calculate delta power
    delta_power = mean(EEG_spectrum(delta_freq_bin1:delta_freq_bin2,:));

    % Calculate theta power
    theta_power = mean(EEG_spectrum(theta_freq_bin1:theta_freq_bin2,:));

    % EMG amplitude is simply average rectified voltage.
    EMG_amplitude = mean(abs(EMG_data));

    % Identify outliers (extremely high delta/theta power or EMG)
    EEG_outliers_delta = delta_power > 100;
    EEG_outliers_theta = theta_power > 100; 
    EMG_outliers = EMG_amplitude > 1000;

    % Remove outliers (turn them into zeros)
    delta_power = delta_power .* (1 - EEG_outliers_delta);
    theta_power = theta_power .* (1 - EEG_outliers_theta); 
    EMG_amplitude = EMG_amplitude .* (1 - EMG_outliers);

    % Eliminate frequencies above Nyquist
    maxBin = round(samplesPerEpoch/2);
    EEG_spectrum(maxBin+1:end,:) = [];
    % Calculate real frequency (in Hertz) for each bin produced by FFT
    frequencies = (0:maxBin-1) / samplesPerEpoch * adfreq;

end

function plotDeltaTheta(delta_power, theta_power, EMG_amplitude)

    numBins = length(delta_power);
    
    % Convert bin # to hours
    time_in_minutes = (1:numBins) * (BIN_SECONDS / 60);
    time_in_hours = time_in_minutes / 60;

    % Plot versus hours
    plotyy(time_in_hours, delta_power, time_in_hours, EMG_amplitude);
    plotyy(time_in_hours, delta_power, time_in_hours, theta_power); 

    legend('delta', 'EMG');
    legend('delta', 'theta');

end

% Remove outlier data points
function [EEG_data, EMG_data] = RemoveOutliers(EEG_data, EMG_data, adfreq)

    % Remove values near +/-32768 (usually related to signal loss)
    EEG_outliers = abs(EEG_data) > 30000;
    EEG_data = EEG_data .* (1 - EEG_outliers);

    % Calculate moving average so that movement artifact outliers can be identified
    kern_pad = 25;
    kern_pad2 = kern_pad * 2 + 1;
    kern1 = ones(kern_pad2,1)/(kern_pad2 - 1);
    kern1(kern_pad + 1) = 0;
    EEG_moving_average = conv(EEG_data, kern1);
    EEG_moving_average = EEG_moving_average(kern_pad+1:(length(EEG_moving_average)-kern_pad));

    % Calculate 10 quantiles
    Y = 10 * quantile(EEG_data, 10);
        
    % Remove outliners that are more than 10x the 10% and 90%
    % percentiles (this assumes that data is zero-centered)
    EEG_outliers = EEG_data < Y(1) | EEG_data > Y(10);

    numPoints = length(EEG_data);
    
    if sum(EEG_outliers) > 0
        fprintf('Warning: removing %d outlying points from EEG (%.2f%% out of %d total points)\n', sum(EEG_outliers), sum(EEG_outliers)/numPoints*100, numPoints);
    end
    
    % Remove outliers (turn them into zeros)
    EEG_data(EEG_outliers) = NaN;

    if length(EMG_data) > 0
        
        % Calculate 10 quantiles
        Y = 10 * quantile(EMG_data, 10);
        
        % Remove outliners that are more than 4x the 10% and 90%
        % percentiles (this assumes that data is zero-centered)
        EMG_outliers = EMG_data < Y(1) | EMG_data > Y(10);
        
        numPointsRemoved = sum(EMG_outliers);
        EMG_data = EMG_data .* (1 - EMG_outliers);

        if sum(EMG_outliers) > 0
            fprintf('Warning: removing %d outlying points from EMG (%.2f%% out of %d total points)\n', numPointsRemoved, numPointsRemoved/numPoints*100, numPoints);
        end

        EMG_data(EMG_outliers) = NaN;

    end
end

function [EEG_data, EMG_data] = FilterEEG(EEG_data, EMG_data, adfreq)
    
    ad2 = adfreq/2;

    % Filter to remove extremely low frequencies (below 0.1Hz) from EEG.
    % This helps remove slow baseline drifts
    [b, a] = butter(2, 0.1/(ad2), 'high');
    EEG_data = filter(b, a, EEG_data);
    
    % Remove <100Hz frequencies from EMG
    if (adfreq > 250)
        [b, a] = butter(2, 100/(ad2), 'high');
        EMG_data = filter(b, a, EMG_data);
        
        % Comb filter to remove 60Hz and harmonics from EMG
        answer = input('Use comb filter to remove 60Hz? Be careful - this can severely exacerbate movement artifact noise: ', 's');
        if strcmp(answer, 'y') || strcmp(answer, 'Y')
            [b, a] = iirnotch(60/(ad2), 5/ad2);
            EMG_data = filter(b, a, EMG_data);
        end
    end
    

    % Filter to gradually attenuate higher frequencies (above 1 Hz) from EEG
    % This low cutoff is needed to compensate for hardware filtering
    % in pre-amplifier, which attenuated the lower frequencies
    fprintf('Was data recorded with hardware high-pass filter? I.e. box 1, prior to end of November 2017?\n');
    answer = input('Enter y or n: ', 's');
    if strcmp(answer, 'y') || strcmp(answer, 'Y')
        [b, a] = butter(1, 1.0/(ad2));
        % Add back 1/10 of original signal so that high frequencies will
        % not disappear entirely.
        EEG_data = filter(b, a, EEG_data) + EEG_data / 10;
    end
end

%
% Takes a list of delta powers, calculates histogram,
% and tries to find point in histogram where density is about 30%
% of peak. This is often "break" between wake and sleep.
function [s] = findThreshold(e, range)

    e = reshape(e, numel(e), 1);
    [h,x] = hist(e, 0:(range/250):range);

    % Last bin might have a lot of datapoints, so zero that out so that we 
    % don't erroneously detect it as being peak.
    h(length(h)) = 0;
    
    % Find histogram peak, which will usually be waking delta power
    maxh = max(h);
    
    if maxh == 0
        % If all bins are zero, then something went wrong. Most likely
        % the range is too low, and all delta powers are all outside of
        % this range.
        s = 0.3;
        return
    end
    
    peakpointIndex = find(h == maxh);
    
    % If histogram has more than one peak, take the first one.
    % This will tend to happen more often if there are fewer epochs
    peakpointIndex = peakpointIndex(1);
    
    if peakpointIndex == 1
        % If peak is at the first bin, it is usually DC
        % offset. Get the second highest.
        h = h(2:end);
        x = x(2:end);
        peakpointIndex = find(h == max(h));
        peakpointIndex = peakpointIndex(1);
    end

    % Normalize entire histogram to peak value
    h = h/max(h);
    
    % Find all indexes where histogram is below peak
    % by specified amount
    index = find(h < .5);

    % Remove indexes that are lower than midpoint
    index(index < peakpointIndex) = NaN;

    % Calculate voltage corresponding to this threshold
    % where histogram first dips below .01
    s = x(min(index));
    
    if isempty(s) || s == 0
        % If histogram is all zero, then just default to 0.35
        s = 0.35;
    end
   
end

function [binBoundaries, binCounts] = plotHistAmplitude(data, threshold, range)

    data = reshape(data, numel(data), 1);
    [binCounts,binBoundaries] = hist(data, 0:(range/250):range);

    % Zero out the last bin, which might have an inordinate number of items
    binCounts(end:end) = 0;
    
    figure
    plot(binBoundaries,binCounts)
    hold on;
    plot([threshold threshold], [0 max(binCounts)]);
%    axis([-100 100 .00001 1]);

end

function plotHistVoltage(data, newFig)

    data = reshape(data, numel(data), 1);
    [h,x] = hist(data, 5000);

    if newFig
        figure
    else
        hold on;
    end
    
    plot(x,h)
%    hold on;
%    plot([threshold threshold], [0 max(h)]);
    axis([-100 100 .00001 * max(h) max(h)]);

end

%time_in_hours=[1:size(SleepStateMatrix,1)]'
%DATA = cat(2,time_in_hours, SleepStateMatrix)
%col_header={'Time in hours', 'Sleep','Awake','REM','Sleep+REM'};     %Row cell array (for column labels)
%filename = 'SleepState.xlsx'    %Join cell arrays
%xlswrite(filename,DATA,'Sheet1','A2');     %Write data
%xlswrite(filename,col_header,'Sheet1','A1');     %Write column header
%xlswrite(filename,row_header,'Sheet1','A2');      %Write row heade