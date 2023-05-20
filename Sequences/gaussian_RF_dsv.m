% This function returns a Gaussian pulse taken from a DSV file from IDEA
% SDE simulation output, of duration 600 us, normalised to RF_amp.  The
% pulse will be squashed/stretched to the RF_dur requested.
%
% Tom Okell, May 2023
%
% [RF, N, mid] = gaussian_RF_dsv(dt, RF_dur, RF_amp)

function [RF, N, mid] = gaussian_RF_dsv(dt, RF_dur, RF_amp)

  % Load the DSV data
  load ~/Documents/Bloch_Simulation/Data/Gaussian_RF_Pulse_DSV.mat;
  
  % Set up the sampling time array and a time array corresponding to the
  % DSV data if the given pulse length had been set.
  t = 0:dt:RF_dur;  % Sampling time array
  DSVdt = RF_dur / (size(Gaussian_RF_Pulse_DSV.t,1)-1); % separation of DSV sampling times
  DSVt = 0:DSVdt:RF_dur;  % DSV sampling times
  
  % Interpolate the DSV data
  RF = interp1(DSVt, Gaussian_RF_Pulse_DSV.RF, t);
    
  % Set the phase to zero since this is set later in other functions
  RF(2,:) = 0;  
  
  % Normalise to give the correct RF amplitude
  RF(1,:) = RF(1,:) * RF_amp;
  
  % Return the size of the array and the midpoint
  N = size(RF,2); mid = round(N/2);
  