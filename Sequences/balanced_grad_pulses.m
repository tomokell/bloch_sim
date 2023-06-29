% Return balanced gradient pulses for e.g. slice selection or readout.
%
% Tom Okell, June 2023
%
% [G, N] = balanced_grad_pulses(dt, flat_top_time, G_s, G_amp, axis)
%
% This function returns balanced (zeroth moment = 0) gradient pulses (for
% e.g. slice selection or readout), sampled at time intervals dt (s). The
% constant gradient section, with amplitude G_amp (mT/m) has duration =
% flat_top_time (s).  The gradient slew rate is set by G_s (mT/m/ms) and
% pre/re-phasing gradients are played out as quickly as possible before and
% after the flat top gradient, subject to rounding to dt.  The gradient is
% played along a direction defined by a 3-element vector, axis. The
% function returns the 3xN gradient waveforms and N, the number of time
% points generated.

function [G, N, flat_top_start_idx] = balanced_grad_pulses(dt, flat_top_time, G_s, G_amp, axis)

if nargin < 5; axis = [0 0 1]'; end

% Normalise axis and convert to a column vector
axis = axis(:) / sqrt(sum(axis(:).^2));

% Convert slew rate into mT/m/s
G_s = G_s * 1e3;

% Round flat_top_time to the nearest dt
flat_top_time = ceil(flat_top_time/dt)*dt;

% Determine the ramp up time, rounding up to the nearest dt
t_ramp = ceil(G_amp/G_s/dt)*dt;

% Adjust slew rate for rounding
G_s_ramp = G_amp / t_ramp;

% Calculate the zeroth moment of the positive lobe
Zeroth_mom = G_amp * (t_ramp + flat_top_time);

% Calculate the pre-phasing/refocussing time required
t_refoc = sqrt(2*Zeroth_mom / G_s);

% Round up to the nearest 2*dt, so the peak aligns with dt also
t_refoc = ceil(t_refoc / (2*dt))*2*dt;

% Calculate the maximum (abs) negative gradient required to match the
% zeroth moment and the slew rate required for this (needs to be adjusted
% due to rounding above)
G_neg = Zeroth_mom / t_refoc;
G_s_neg = G_neg / (t_refoc/2);

% Determine the total time required and number of data points:
T = 2*t_refoc + 2*t_ramp + flat_top_time;
N = round(T/dt);

% Calculate an array of time relative to the start time
t = (0:(N-1))*dt;

% Initialise the array
G = zeros(1,N);

% Loop through time, calculating the appropriate gradient size at
% this point
for ii=1:N
    if t(ii) <= t_refoc/2  % Pre-phaser ramp up
        G(ii) = -t(ii)*G_s_neg;

    elseif t(ii) <= t_refoc % Pre-phaser ramp down
        G(ii) = (t(ii) - t_refoc)*G_s_neg;

    elseif t(ii) <= t_refoc + t_ramp  % Flat top ramp up
        G(ii) = (t(ii) - t_refoc)*G_s_ramp;

    elseif t(ii) <= t_refoc + t_ramp + flat_top_time  % Flat top
        G(ii) = G_amp;

    elseif t(ii) <= t_refoc + t_ramp + flat_top_time + t_ramp % Flat top ramp down
        G(ii) = -(t(ii) - (t_refoc + t_ramp + flat_top_time + t_ramp))*G_s_ramp;

    elseif t(ii) <= t_refoc + t_ramp + flat_top_time + t_ramp + t_refoc/2  % Re-phaser ramp up
        G(ii) = -(t(ii) - (t_refoc + t_ramp + flat_top_time + t_ramp))*G_s_neg;

    else % Re-phaser ramp down
        G(ii) = (t(ii) - T)*G_s_neg;

    end
end

% Apply the gradient to the correct axis, making it a 3xN matrix
G = G .* axis;

% Return the index of the start of the flat top gradient
flat_top_start_idx = find(t>=t_refoc + t_ramp,1,'first');
