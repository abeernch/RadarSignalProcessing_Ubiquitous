%% POI Analytical Modeling and simulated verification for a scanning receiver
%% POI Modeling
% The analysis invloves modeling two independent scanning systems, a
% transmitter and a receiver.
% 1. Receiver: 
% Wide band receiver that scans 1-18GHz with an IBW of 1GHz having a defind
% scan tiem for each band and a switching time during which the frequency
% bins are changed.

% 2. Transmitter:
% THe transmitter is modelled using the time specific parameters such as
% the PW, PRI, Scan time/rate and beamwidth, so as to establish a 
% probabalistic relationship between the two asynchronous systems.

%% Modeling:
% The simulation determines how many pulses from the radar are intercepted 
% by the receiver while:
% - The radar beam is actually pointing toward the receiver (beam-on)
% - The receiver is tuned to the correct frequency band at that exact 
% moment (frequency overlap)
% - The radar pulse in time overlaps the receiver dwell window (temporal 
% overlap)

% The analytical estimates for the POI are broken down into 3 sequential
% probabilities which are multiplied to see the full effect.
% 1. Probability on beam:
% This accounts for the fact that a scanning radar is not always
% illuminating the receiver. It is modeled as a fraction of the dwell time
% in the receivers direction to the total scan time, representing the
% likelihood that the receiver is/is not illuminated.
% 2. Pulse-Dwell Overlap:
% This models the probability the receiver is in the band in which a
% transmitter is transmitting. It is the ratio of the time of a single bin
% scan and the total time to scan the 17 bins.
% 3. Pulse overlap
% Pulses may partially overlap with the active frequency bin. This time
% overlap is modeled as well together with 2.

%% Simulation validation
%% Phase Aware Simulation
% A time scheduling simulation with randomized phases for both the receiver
% and the transmitter is implemented by modeling scan window timing, PWs
% and phase offsets. This counts the total number of pulses that are
% intercepted given the receiver and transmitter properties also while
% randomizing the initial phases to capture the randomness of
% synchronization between the two systems.

% The initial phase randomization can be turned off. In that event all or
% none of the pulses may overlap given the PRI,rxDwell ratios.

% The two outputs can be compared to counter-verify the accuracy of the 
% analytical and simulated implementations for a given set of receiver and
% transmitter parameters.
% Author: Abeer Chaudhry

clc;clear;close all

%% Receiver parameters
rxIBW = 1e9;                  % Hz 
rxMinFreq = 1e9; rxMaxFreq = 18e9;    % scanned band (Hz)
B = (rxMaxFreq - rxMinFreq)/rxIBW;      % bins

rxDwell = 30000e-6;                      % band scan time sec
rxSwitch = 100e-6;                       % receiver switching time between bins (sec)
rxCycle = B*(rxDwell + rxSwitch);    % sec

%% Radar transmitter parameters 

simTime = 60;             % total simulation time (sec)
txBeamwidth = 3;         % deg (beamwidth)
txScanRate = 6*6;          % RPM->deg/sec
txScanTime = 360/(txScanRate); % sec
txDwell = txBeamwidth/(6*txScanRate);
txPRI = 1.4e-3;           % PRI (sec)
txPW = 1e-6;                % pulse width (s)
pulses_per_dwell = txDwell / txPRI;
N_total = simTime/txPRI; % Total pulses transmitted per simulation

% Simulation control 
omni = 0; % Enable if the beam on is irrelevant 
phaseAvg = 1;    % if true, average over several random phase offsets
n_phase_samples = 50;       % number of random phase samples to average when random phase

%% Analytical Estimate 

% Probability of beam-on
P_beam = txDwell/txScanTime; % probability emulating the radar beam on with the receiver given its configuration params

% Actual no. of pulses received at the receiver
if omni 
    N_emit = N_total;
else
    N_emit = N_total*P_beam; 
end

% probability of dwell-pulse overlap
POI_overlap = min( (rxDwell + txPW) / rxCycle, 1 );

N_overlap = N_emit * POI_overlap;

fprintf('Analytical Estimates \n');
fprintf('B = %d bins, T_cycle = %.3e s\n', B, rxCycle);
fprintf('Probability intercept p = %.5f\nIntercepted Pulses = %.2f pulses of %g pulses\n\n', POI_overlap, N_overlap,N_emit);

%% Time Scheduling simulation 
% Build receive windows for the target band
% Assume the target radar is in a single known IBW bin index 
% For example, if radar occupies 1-2 GHz and the first bin is 1-2 GHz, target_bin_idx = 1

target_bin_idx = 12;

% receiver cycle windows for bin i inside cycle:
% offset to target bin within each cycle:
bin_offset = (target_bin_idx - 1)*(rxDwell + rxSwitch);  % seconds from cycle start to target bin start

% For deterministic counting, construct a list of pulse start times 
% Assume the N_emit pulses are all directed at receiver and occur at times t_k = t0 + k*T_PRI
% Choose simulation length to include all pulses:
T_sim_from_emit = (N_emit-1)*txPRI + txPW;


% run deterministic counting either single-phase or average random phases
if phaseAvg
    counts = zeros(n_phase_samples,1);
    rng(0); % reproducible
    for s=1:n_phase_samples
        rx_phase = rand()*rxCycle;        % random start time of RX scanning cycle
        radar_phase = rand()*txPRI;      % random offset for radar pulse train
        counts(s) = count_for_phases(rx_phase, radar_phase, T_sim_from_emit, ...
                        rxCycle, bin_offset, rxDwell, txPW, N_emit, txPRI);
    end
    mean_count = mean(counts);
    std_count = std(counts);
    fprintf('Phase-aware (random-phase avg over %d samples)\n', n_phase_samples);
    fprintf('Mean counted intercepted pulses = %.2f\n', mean_count);
else
    rx_phase = 0;
    radar_phase = 0;
    Ncount = count_for_phases(rx_phase, radar_phase, T_sim_from_emit, ...
                        rxCycle, bin_offset, rxDwell, txPW, N_emit, txPRI);
    fprintf('--- Deterministic / Phase-aware (single run) ---\n');
    fprintf('Counted intercepted pulses = %d\n', Ncount);
end

% Fn to count pulses for given phase offsets (rx_phase, radar_phase)
function Ncount = count_for_phases(rx_phase, radar_phase, T_sim_from_emit, T_cycle, bin_offset, T_dwell_rx, T_PW, N_emit, T_PRI)
    Ncount = 0;
    % Precompute receive windows up to T_sim_from_emit
    n_cycles = ceil((T_sim_from_emit + rx_phase)/T_cycle) + 1;
    win_start = ((0:(n_cycles-1))*T_cycle) + rx_phase + bin_offset;
    win_end   = win_start + T_dwell_rx;
    % now for each pulse:
    for k = 0:(N_emit-1)
        t_p = radar_phase + k*T_PRI;          % pulse start
        t_p_end = t_p + T_PW;                   % pulse end
        % if pulse beyond simulation horizon, break
        if t_p > (T_sim_from_emit + max(rx_phase, radar_phase)), break; end
        % check overlap: does [t_p, t_p_end] overlap any [win_start, win_end]?
        % find candidate windows whose start <= t_p_end and end >= t_p
        % do vectorized check in MATLAB: (faster in compiled MATLAB) - keep simple loop here
        % but we can quick check if any window satisfies:
        idx = find( (win_start <= t_p_end) & (win_end >= t_p), 1 );
        if ~isempty(idx)
            Ncount = Ncount + 1;
        end
    end
end