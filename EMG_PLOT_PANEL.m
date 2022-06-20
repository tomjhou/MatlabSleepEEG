
function EMG_PLOT_PANEL(BIN_SECONDS, adfreq, centralEpoch, numSamplesPerEpoch, numEpochs, ...
    x_values_seconds, SleepState, ...
    EMG_stacked, threshold_EMG_stacked, EMG_amplitude_stacked, ...
    EPOCH_PADDING, MAX_EMG_AMPLITUDE)

global strings;

% Plot EMG
plot(gca, x_values_seconds, EMG_stacked);
% Get highest value plus one
max_x = x_values_seconds(end) + 1/adfreq;
axis([x_values_seconds(1) max_x -MAX_EMG_AMPLITUDE MAX_EMG_AMPLITUDE]);

title('EMG');

% Turn on hold so that new traces will overlay old
hold on;

if EPOCH_PADDING > 0
    lineX1 = 0;
    lineX2 = numSamplesPerEpoch/adfreq;
    lineY1 = -MAX_EMG_AMPLITUDE;
    lineY2 =  MAX_EMG_AMPLITUDE;

    % Plot vertical lines
    plot(gca, [lineX1 lineX1], [lineY1 lineY2], 'g', 'LineWidth', 3);
    plot(gca, [lineX2 lineX2], [lineY1 lineY2], 'g', 'LineWidth', 3);
end

% plot EMG threshold line
plot(gca, x_values_seconds, threshold_EMG_stacked, '--k');
plot(gca, x_values_seconds,  EMG_amplitude_stacked, 'm');
plot(gca, x_values_seconds, -EMG_amplitude_stacked, 'm');
% plot axis (y=0 line)
plot(gca, x_values_seconds, zeros(length(EMG_amplitude_stacked), 1), 'k');

% Place tick marks every 12 seconds on graph (default is 10 seconds,
% but we want tick to match epoch boundaries)
r = round(BIN_SECONDS);
if EPOCH_PADDING == 0
    xTicks = [0:r]';
else
    xTicks = [(-r * EPOCH_PADDING):r:(r*(EPOCH_PADDING+1))]';
end

set(gca,'XTick',xTicks);

for offset = -EPOCH_PADDING:EPOCH_PADDING
    epoch = centralEpoch + offset;
    if epoch < 1 || epoch > numEpochs
        continue;
    end
    
    if SleepState(epoch) >= 0
        sleep_state_string = strings(1 + SleepState(epoch));
    else
        sleep_state_string = '-';
    end
    
    text(offset * BIN_SECONDS, MAX_EMG_AMPLITUDE * .9, sleep_state_string);
end