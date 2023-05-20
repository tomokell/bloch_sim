% This function calculates the flip angle (in degrees) due to a general RF
% pulse given in the RF array, with sampling time dt.
%
% Tom Okell, May 2023
%
% FlipAng = FlipAngle(RF_pulse, dt)

function FlipAng = FlipAngle(RF_pulse, dt)

% Define gamma
g = GetGamma;  % rad/s/mT            

% Theta = integral{ gamma * B1 * dt }
FlipAng = torad2deg(g * dt * sum(RF_pulse(1,:)) );