% This function creates RF and gradients for a CASL type sequence based
% on a time array, t, a B1 field strength in mT, and gradient field
% strength, G_amp in mT/m.
%
% Tom Okell, May 2023
%
% [G, RF] = make_CASL_seq(t, RF_amp, G_amp)

function [G, RF] = make_CASL_seq(t, RF_amp, G_amp)
  
  % Generate a continuous RF waveform of size B1
  RF = zeros(2,size(t,2));
  RF(1,:) = RF_amp;
  
  % Generate a continuous gradient in the z direction of amplitude G_amp
  G = zeros(3,size(t,2));
  G(3,:) = G_amp;