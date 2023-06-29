% VFA scheme for N pulses that ensures any attenuation from the previous
% pulse is cancelled out by the increase in FA for the next pulse, giving
% constant signal in the absence of T1 decay, similar to Wang MRM 1991, but
% in backwards recursive form so you can define the final flip angle,
% FAmax (degrees).
%
% Tom Okell, June 2022
%
% FAs = MaintainVFA(N,FAmax)

function FAs = MaintainVFA(N,FAmax)

% Initialise
FAs = zeros(1,N);

% Recursive calculation
FAs(N) = FAmax;
for ii = (N-1):-1:1
    FAs(ii) = atan(sin(FAs(ii+1)*pi/180))*180/pi;
end