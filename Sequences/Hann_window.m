% This function returns a Hanning window shaped pulse with N samples and
% unit amplitude
%
% Tom Okell, May 2023
%
% RF_pulse = Hann_window(N)

function RF_pulse = Hann_window(N)
  
% Set up an index
  i = 0:(N-1);
    
  % Calculate the Hann window function
  % Ref: https://en.wikipedia.org/wiki/Hann_function
  RF_pulse = 0.5 * ( 1 - cos ( 2*pi*i / (N - 1) ) );