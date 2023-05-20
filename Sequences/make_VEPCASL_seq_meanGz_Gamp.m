% This function creates RF and gradients for a (vessel-encoded) PCASL type
% sequence based on an RF pulse shape, RF_shape, with shape parameters,
% RF_shape_params, RF peak amplitude in mT, RF_amp, pulse duration, RF_dur
% (in s), pulse separation, RF_sep (in s), the time step, dt (in s), total
% time for RF and gradient pulses to occur, T (in s), mean z gradient,
% meanGz (in mT/m), the 2D position vectors of vessels A and B within the
% labelling plane, Pa and Pb, the z-offset of the labelling plane,
% z_offset, and a cycle number between 1 and 4: 1 = all inverted (standard
% PCASL label), 2 = none inverted (standard PCASL control), 3 = vessel A
% inverted whilst vessel B is controlled, 4 = vessel B inverted whilst
% vessel A is controlled. For standard PCASL simulations, only cycles 1 and
% 2 are required and Pa and Pb can be set arbitrarily.
%
% Tom Okell, May 2023
%
% [G, RF] = make_VEPCASL_seq_meanGz_Gamp(RF_shape, RF_shape_params, RF_amp, RF_dur, RF_sep, dt, T, meanGz, G_amp, Pa, Pb, z_offset, cycle)

function [G, RF] = make_VEPCASL_seq_meanGz_Gamp(RF_shape, RF_shape_params, RF_amp, ...
                                     RF_dur, RF_sep, dt, T, meanGz, G_amp, ...
                                     Pa, Pb, z_offset, cycle)

  % Declare gamma
    g = GetGamma;  % rad/s/mT

  % Check that the pulse duration is shorter than the required separation
    if RF_dur > RF_sep
      error('!!! Requested RF duration is more than the RF separation...!!!')
    end
    
  % Generate the time array, rounding down to the nearest pulse
    N_pulses = floor( T/RF_sep );
    disp([ num2str(N_pulses) ' pulses are being generated for this sequence'])
    t = 0:dt:T;
    
  % Generate a single RF pulse of the required shape
  pulse_samples = round(RF_dur/dt);
  
  switch lower(RF_shape)
   case 'hard'
    RF_pulse = ones(1,pulse_samples) * RF_amp;
    RF_pulse(2,:) = 0;
    
   case 'gaussian'
    % gaussian_RF(pulse_amp, temporal SD, phi, samp_time, dt, HannWindow)
    RF_pulse = gaussian_RF(   RF_amp, RF_shape_params(1),   0,    RF_dur, dt, false);
    
   case 'gaussianhann'
    % gaussian_RF(pulse_amp, temporal SD, phi, samp_time, dt, HannWindow)
    RF_pulse = gaussian_RF(   RF_amp, RF_shape_params(1),   0,    RF_dur, dt, true);

   case 'sinc'
    % Create using an arbitrary flip angle before rescaling
    RF_pulse = sinc_RF(10, 0, RF_dur, dt, G_amp, g, RF_shape_params(1) );
    RF_pulse(1,:) = RF_pulse(1,:) / max(RF_pulse(1,:)) * RF_amp;
    
   case 'hanning'
    RF_pulse = Hann_window(pulse_samples)*RF_amp;
    RF_pulse(2,:) = 0;
    
   case 'hamming'
    RF_pulse = Hamming_window(pulse_samples)*RF_amp;
    RF_pulse(2,:) = 0;
    
   case 'gaussian_dsv' % From DSV data, stretched to fit RF_dur
	
	RF_pulse = gaussian_RF_dsv(dt, RF_dur, RF_amp);
	
   otherwise
    disp('>>> !!! Error - RF pulse shape not recognised... !!! <<<')
    disp('Returning zeros as RF pulse...')
    RF_pulse = zeros(2,pulse_samples);
    
  end
  
  % Check the duration of the RF pulse (some functions add or remove one data
  % point to make the pulse symmetric about its centre
  pulse_samples = size(RF_pulse,2);

  % Calculate the vessel separation and normalised gradient direction
  D = vecmod(Pb - Pa);
  Gxy = (Pb - Pa) / D;
  Da = sum(Pa .* Gxy);
  Db = sum(Pb .* Gxy);  % Distance of each vessel from isocentre along Gxy
  
  % Generate a gradient pulse
  [G_pulse, N] = VEPCASL_grad_pulse_meanGz_Gamp(g, RF_sep, dt, RF_dur, meanGz, ...
                                     G_amp, Gxy, D);
  
  % Find the mean Gz for later RF phase calculations
  mean_Gz = mean(G_pulse(3,:) );
  
  % Create N_pulses copies of this RF shape and the associated gradient
  % waveform
  RF = zeros(2,size(t,2));    G = zeros(3,size(t,2));
  
  for ii = 1:N_pulses
    pulse_start = 1 + round( (ii-1)*RF_sep/dt );
    pulse_end = round( ii*RF_sep/dt );
    pulse_length = min(pulse_end - pulse_start + 1, N);
    
    RF(:, pulse_start:(pulse_start + pulse_samples - 1) ) = RF_pulse;
  
    % Calculate the RF phase at this time
    phi_z = ii * g * mean_Gz * t(ii) * z_offset;  % Accounts for offset of
                                                % labelling plane from
                                                % isocentre along z
    % Note that the phase terms below use i-1 rather than i since matlab
    % counts from the first pulse rather than the zeroth pulse, and
    % phi_xyA and phi_xyB are only required after the first pulse
    phi_xyA = mod(ii-1,2) * pi * Da / D;  % Keeps in phase with vessel A
    phi_xyB = mod(ii-1,2) * pi * Db / D;  % Keeps in phase with vessel B
    
    % Check which gradient and RF phase cycling scheme is used
    switch cycle
     case 1
      Gxy_mod = 0;
      RF_ph = phi_z;
      
     case 2
      Gxy_mod = 0;
      RF_ph = phi_z + mod(ii-1,2)*pi;
      
     case 3
      Gxy_mod = ( 2*mod(ii,2) - 1 );  % Alternately positive and negative
      RF_ph = phi_z + phi_xyA;
      
     case 4
      Gxy_mod = ( 2*mod(ii,2) - 1 );  % Alternately positive and negative
      RF_ph = phi_z + phi_xyB;
      
     otherwise
      error('WARNING! Cycle number not set between 1-4!!!'); 
      
    end
    
    % Insert the pulse into the array
    G(:,pulse_start:(pulse_start + pulse_length - 1) ) = G_pulse(:,1:pulse_length);
    
    % Modulate the Gxy pulse
    G(1:2, pulse_start:(pulse_start + pulse_length - 1) ) = ... 
        G(1:2, pulse_start:(pulse_start + pulse_length - 1) ) * Gxy_mod;
    
    % Insert the RF phase
    RF(2,pulse_start:(pulse_start + pulse_samples - 1) ) = RF_ph;
    
  end
  