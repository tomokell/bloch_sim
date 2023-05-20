% This function returns a gradient pulse with duration T seconds, sampled
% at time intervals dt.  The shape and size of the z pulse is determined by
% the gyromagnetic ratio, g (in rad/s/mT), the RF duration, RF_dur (in s),
% the gradient strength, G_amp, and mean gradient, meanGz, in mT/m.  The
% shape and size of the xy pulse is determined by the gradient direction,
% Gxy_dir, and vessel separation, D (in m). The number of data points, N is
% also returned.
%
% Tom Okell, May 2023
%
% [G, N] = VEPCASL_grad_pulse_meanGz_Gamp(g, T, dt, RF_dur, meanGz, G_amp, Gxy_dir, D)

function [G, N] = VEPCASL_grad_pulse_meanGz_Gamp(g, T, dt, RF_dur, meanGz, ...
                                             G_amp, Gxy_dir, D)
  
  % Determine the required number of data points:
  N = round(T/dt);
  
  % Initialise the array
  G = zeros(3,N);

  % Calculate ramp time
  t_r = (T - RF_dur)/2;
  
  % Calculate an array of time relative to the start time
  t = (0:(N-1))*dt;

  % Calculate the required slew rate
  G_s = (G_amp - meanGz) * T / (t_r^2);
  
  % Loop through the ramp, calculating the appropriate z gradient size at
  % this point
  for ii=1:N
    if t(ii) < RF_dur  % Flat top
      G(3,ii) = G_amp;
      
    elseif t(ii) < RF_dur + t_r  % Ramp down
      G(3,ii) = G_amp + (t(ii) - RF_dur)*(-G_s);
      
    else  % Ramp up
      G(3,ii) = G_amp + (t(ii) - T)*G_s;
      
    end
    
  end
  
  
  % Now calculate the required area and thus amplitude of the Gxy pulse
  Gxy_area = pi / (g * D);
  Gxy_amp = 2 * Gxy_area / (T - RF_dur);
  
  Gxy = zeros(1,N);
  
  % Calculate Gxy at each time point
  for ii=1:N
    if t(ii) < RF_dur  % No gradient
      Gxy(ii) = 0;
      
    elseif t(ii) < RF_dur + t_r  % Ramp up
      Gxy(ii) = Gxy_amp * (t(ii) - RF_dur) / t_r;
      
    else  % Ramp down
      Gxy(ii) = Gxy_amp * (T - t(ii)) / t_r;
      
    end
    
  end
  
  % Split Gxy into x and y components
  Gxy_dir = Gxy_dir / vecmod( Gxy_dir );  % Normalise
  G(1,:) = Gxy_dir(1) * Gxy;
  G(2,:) = Gxy_dir(2) * Gxy;