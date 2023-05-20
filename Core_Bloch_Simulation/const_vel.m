% Generate a constant velocity position vector array
%
% Tom Okell, May 2023
%
% P = const_vel(start_pos, v, t)
%
% Generates a position waveform for the case of constant velocity for use
% in Bloch simulation code.  Input parameters are the starting position
% vector at t = 0 (start_pos), the velocity vector (v), and the array of
% time points to evaluate (t).  Units are m and s.

function P = const_vel(start_pos, v, t)
  
   % Initialise the P array to the correct size
   P = zeros(3, size(t,2));
   
   % For each time point, fill the P array with the appropriate position
   for ii = 1:size(t,2)
     
     P(:,ii) = start_pos + v * t(ii);
     
   end
  