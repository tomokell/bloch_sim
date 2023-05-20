% This function returns a Gaussian pulse with pulse amplitude, RF_amp,
% temporal standard deviation, sig, phase phi (radians) and duration T
% seconds, sampled at time intervals dt.  If HannWindow is true, a Hann
% window is also applied. The number of data points, N, and central data
% point, mid, are also returned.
%
% Tom Okell, May 2023
%
% [RF, N, mid] = gaussian_RF(RF_amp, sig, phi, T, dt, HannWindow)

function [RF, N, mid] = gaussian_RF(RF_amp, sig, phi, T, dt, HannWindow)
  
  if nargin < 6; HannWindow = false; end

  % Determine the required number of data points:
  N = round(T/dt);
  
  % Round up to an odd number to make the sampling symmetric
  if mod(N, 2) == 0  % N is even
    N = N + 1;
  end
  
  % Determine the central point of the array
  mid = round(N/2); % NB. rounds up odd numbers/2 i.e.round(5/2)=3
  
  % Initialise the array
  RF = zeros(2,N); % NB. Two rows for amplitude and phase
  RF(2,:) = phi;  % Define the phase as constant throughout the pulse
  
  % Give a warning if the requested duration is less than 4 * sig and Hann
  % windowing is not being applied or dt is large compared to sig
  if T < (4 * sig) && ~HannWindow
    disp(['!!! Warning: Requested RF duration is less than 4 standard ' ...
          'deviations of the Gaussian profile - clipping will occur !!!'])
  end
  
  if sig < (4*dt)
    disp('!!! Warning, dt too large to sample RF properly !!!');
  end
  
  % Calculate the Gaussian function
  % NB. at this point the pulse is not normalised
  t = ((1:N) - mid) * dt;  % time relative to pulse centre
  RF(1,:) = exp(-0.5 * ( t/sig ).^2 );  % Calculate the unnormalised Gaussian function
    
  % Normalise to give the correct RF amplitude
  RF(1,:) = RF(1,:) * RF_amp;
  
  % Apply a Hann window if requested
  if HannWindow
    RF(1,:) = RF(1,:) .* Hann_window(N);
  end