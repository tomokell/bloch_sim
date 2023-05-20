% This function calculates the modulus of the supplied vector, v
%
% Tom Okell, May 2023
%
% M = vecmod(v)

function M = vecmod(v)
  
  M = sqrt( sum( v.^2 ) );