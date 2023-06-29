% Calculate the ASL signal contrast between tag and control conditions in
% the presence of an SSFP/SPGR readout using Bloch simulation.
% 
% Tom Okell, June 2023
%
% [C, tC, Mtag, Mcntl, P, G, RF, t, TEIdx, Ccpx, Cmag, ADC] = test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses, LHR, VFAType, G_s, delta, ADC_dur)
%
% This function runs two simulations: one with the magnetisation starting
% after a perfect inversion [0, 0, -1] and the other after a perfect
% control [0, 0, 1]. It then calculates the absolute value of the complex
% difference signal, C, at each echo time, tC. Also returned are the full
% tag and control magnetisation simulations, Mtag and Mcntl, the position
% array of the spins, P, gradient and RF waveforms, G and RF, a time array
% for the full simulation, t, an index of echo times within t, TEIdx. In
% addition, the complex magnetisation difference, Ccpx, can also be
% returned, along with a magnitude subtraction, Cmag (i.e. abs(Mcntl) -
% abs(Mtag)), and a logical array showing when the ADC is on, ADC.
%
% Inputs: many sequence parameters (FlipAng, RF_dur, TR, G_amp, SeqType,
% NoLinIncPulses, VFAType, G_s, delta, ADC_dur) are defined in
% make_SSFP_seq.m. Other parameters:
%   v               =   Spin velocity (m/s) in the z direction
%   T               =   Total sequence duration (s)
%   Ps              =   Initial spin position [x y z] (m)
%   dt              =   Time increment for waveform sampling (s)
%   T1              =   Longitudinal relaxation time (s)
%   T2              =   Transverse relaxation time (s)
%   ORFreq          =   Off-resonance frequency (Hz)
%   LHR             =   Use left-hand rule rotations for Bloch simulation

function [C, tC, Mtag, Mcntl, P, G, RF, t, TEIdx, Ccpx, Cmag, ADC] = ...
            test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, ...
                                   Ps, dt, T1, T2, ORFreq, SeqType, ...
                                   NoLinIncPulses, LHR, VFAType, ...
                                   G_s, delta, ADC_dur)

  if nargin < 12; SeqType = 'HalfAngle';  end
  if nargin < 13; NoLinIncPulses = 10;    end
  if nargin < 14; LHR = true;             end
  if nargin < 15; VFAType = 'CFA';        end
  if nargin < 16; G_s = 200;              end
  if nargin < 17; delta = 0.2;            end
  if nargin < 18; ADC_dur = TR/3;         end
  
  % Calculate the magnetisation over time for the tag condition
  disp('Simulating tag condition...')
  Ms = [0 0 -1]';  % Start with a perfect inversion
  [Mtag, t, P, G, RF, TEIdx, ADC] = test_SSFP_seq(v, Ms, G_amp, FlipAng, T, RF_dur, TR, ...
	                               Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses, LHR, VFAType, ...
                                   G_s, delta, ADC_dur);
							   
  % Calculate the magnetisation over time for the control condition
  disp('Simulating control condition...')
  Ms = [0 0 1]'; % Start at equilibrium
  [Mcntl, ~, ~, ~, ~, ~, ~] = test_SSFP_seq(v, Ms, G_amp, FlipAng, T, RF_dur, TR, ...
	                               Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses, LHR, VFAType, ...
                                   G_s, delta, ADC_dur);	
							   

  % Calculate the sampling time array
  disp('Calculating contrast...')
  tC = t(TEIdx); % Times to calculate the ASL contrast at
  
  % Calculate the contrast (magnitude of the complex difference in
  % transverse magnetisation)
  C = abs( Mcntl(1,TEIdx) + 1i*Mcntl(2,TEIdx) - Mtag(1,TEIdx) - 1i*Mtag(2,TEIdx) );

  % Also calculate the complex and magnitude only contrast
  Ccpx = Mcntl(1,TEIdx) + 1i*Mcntl(2,TEIdx) - Mtag(1,TEIdx) - 1i*Mtag(2,TEIdx);
  Cmag = abs(Mcntl(1,TEIdx) + 1i*Mcntl(2,TEIdx)) - abs(Mtag(1,TEIdx) - 1i*Mtag(2,TEIdx));