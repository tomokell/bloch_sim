% This function returns an RF since pulse with flip angle alph
% (degrees), phase phi (radians) and duration T seconds, sampled
% at time intervals dt.  The shape and size of the pulse is determined by the
% gradient strength, G_amp, gamma for protons, and the slab width, delta.  The number
% of data points, N, and central data point, mid, are also returned.
%
% Tom Okell, May 2023
%
% [RF, N, mid] = sinc_RF(alph, phi, T, dt, G_amp, gamma, delta)

function [RF, N, mid] = sinc_RF(alph, phi, T, dt, G_amp, gamma, delta)
  
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
  
  % Loop through the array, calculating the sinc function at that point
  % NB. at this point the pulse is not normalised
  RFsum = 0;
  for ii = 1:N
    t = (ii - mid) * dt;  % current time relative to pulse centre
    
    if t ~= 0  % The following division will fail at t = 0 so check here
      sinc_arg = 0.5 * gamma * G_amp * delta * t;
      RF(1,ii) = sin( sinc_arg ) / sinc_arg;
      
    else  % At t = 0, define sin(0)/0 = 1
      RF(1,ii) = 1;
      
    end
    
    % Keep a running total of the sum of RF amplitudes for normalisation
    % to the RF flip angle alpha.
    RFsum = RFsum + RF(1,ii);
  end
  
  % Normalise to give the correct flip angle
  RF(1,:) = RF(1,:) * todeg2rad(alph) / ( gamma * RFsum * dt );
  