
function EEG_PLOT_PANEL(BIN_SECONDS, adfreq, centralEpoch, numSamplesPerEpoch, numEpochs, ...
    x_values_seconds, SleepState, ...
    EEG_stacked, threshold_EEG_stacked, delta_stacked, theta_stacked, DELTA_SCALE_FACTOR, ...
    EPOCH_PADDING, MAX_EEG_AMPLITUDE)

global strings;

% Plot EEG (erasing old plot)
plot(gca, x_values_seconds, EEG_stacked);

% Get highest value plus one
max_x = x_values_seconds(end) + 1 / adfreq;

% Control y-axis scaling
axis([x_values_seconds(1) max_x -MAX_EEG_AMPLITUDE MAX_EEG_AMPLITUDE]);

% Title is # of seconds into file
startEpoch = centralEpoch - EPOCH_PADDING - 1;
endEpoch = centralEpoch + EPOCH_PADDING;
if startEpoch < 0
    startEpoch = 0;
end
if endEpoch >= numEpochs
    endEpoch = numEpochs - 1;
end

title(['EEG: ' num2str(startEpoch * BIN_SECONDS / 60) ' - ' num2str(endEpoch * BIN_SECONDS / 60) ' minutes']);

% Turn on hold so that overlays will work
hold on;

% Calculate positions of vertical lines
if EPOCH_PADDING > 0
    lineX1 = 0;
    lineX2 = numSamplesPerEpoch/adfreq;
    lineY1 = -MAX_EEG_AMPLITUDE;
    lineY2 =  MAX_EEG_AMPLITUDE;
    % Plot vertical lines
    plot(gca, [lineX1 lineX1], [lineY1 lineY2], 'g', 'LineWidth', 3);
    plot(gca, [lineX2 lineX2], [lineY1 lineY2], 'g', 'LineWidth', 3);
end

% plot EEG threshold line (black dashed line)
plot(gca, x_values_seconds, threshold_EEG_stacked * DELTA_SCALE_FACTOR, '--k');
% plot delta power
plot(gca, x_values_seconds, delta_stacked * DELTA_SCALE_FACTOR, 'm');
% plot theta power
plot(gca, x_values_seconds, theta_stacked * DELTA_SCALE_FACTOR, 'r');
% plot axis (y=0 line)
plot(gca, x_values_seconds, zeros(length(delta_stacked), 1), 'k');

% Place tick marks every 12 seconds on graph (default is 10 seconds,
% but we want tick to match epoch boundaries)
r = round(BIN_SECONDS);
if EPOCH_PADDING == 0
    xTicks = [0:r]';
else
    xTicks = [(-r * EPOCH_PADDING):r:(r*(EPOCH_PADDING+1))]';
end

xLabels = cellstr(num2str(xTicks));

% Default ticks are every 10 units. Change to every 12 seconds.
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

    text(offset * BIN_SECONDS + 1, MAX_EEG_AMPLITUDE * .9, sleep_state_string);
    text(offset * BIN_SECONDS, MAX_EEG_AMPLITUDE * -.8, num2str(epoch));
    
end

end

