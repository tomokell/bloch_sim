% This function returns an RF sinc pulse with flip angle alph (degrees),
% phase phi (radians) and duration T (seconds), sampled at time intervals
% dt (s). The shape and size of the pulse is determined by the gradient
% strength, G_amp (mT/m), gamma for protons (rad/s/mT), and the slab width,
% delta (m).  The number of data points, N, and central data point, mid,
% are also returned. If HannWindow is set to true, a Hann window is also
% applied to the sinc pulse to reduce ringing in the slice profile, at a
% cost of a broader transition band.
%
% Tom Okell, June 2023
%
% [RF, N, mid] = sinc_RF(alph, phi, T, dt, G_amp, gamma, delta, HannWindow)

function [RF, N, mid] = sinc_RF(alph, phi, T, dt, G_amp, gamma, delta, HannWindow)
  
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
    
  end

  % Apply a Hann Window if requested
  if HannWindow
      RF(1,:) = RF(1,:) .* Hann_window(N);
  end
  
  % Keep a running total of the sum of RF amplitudes for normalisation
  % to the RF flip angle alpha
  RFsum = sum(RF(1,:));
  
  % Normalise to give the correct flip angle
  RF(1,:) = RF(1,:) * todeg2rad(alph) / ( gamma * RFsum * dt );

  