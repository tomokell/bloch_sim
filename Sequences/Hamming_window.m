% This function returns a Hamming window shaped pulse with N samples and
% unit amplitude
% 
% Tom Okell, May 2023
%
% RF_pulse = Hamming_window(N)

function RF_pulse = Hamming_window(N)
  
  % Set up an index
  i = 0:(N-1);
    
  % Calculate the Hamming window function
  % Ref: http://en.wikipedia.org/wiki/Hamming_function#Hamming_window  19/5/08
  RF_pulse = 0.53836 - 0.46164 * cos ( 2*pi*i / (N - 1) );