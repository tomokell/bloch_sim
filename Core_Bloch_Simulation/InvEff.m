% Calculate inversion efficiency of an ASL simulation.
%
% Tom Okell, May 2023
%
% IE = InvEff(Final_Mz, T_final, T_inv, T1)
%
% This function calculates the inversion efficiency of a sequence given the
% final z magnetisation, Final_Mz, the time at which Final_Mz was measured,
% T_final, the nominal inversion time (i.e. the time at which the spin
% passes through the centre of the labelling plane), T_inv, and the T1.

function IE = InvEff(Final_Mz, T_final, T_inv, T1)
  
  IE = 0.5 * (1 - Final_Mz) .* exp( (T_final - T_inv) / T1 );