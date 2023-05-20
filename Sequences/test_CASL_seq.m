% Tests the adiabatic inversion caused by the CASL sequence for a spin of
% velocity v (cm/s), maximum z position zmax (m), gradient amplitude G_amp
% (mT/m), RF amplitude RF_amp (mT), off-resonance frequency ORFreq (Hz),
% time interval for simulation dt (s), longitudinal relaxation time T1 (s),
% transverse relaxation time T2 (s), boolean to set negative z gradient
% NegGz, and boolean to set left-hand rule rotations LHR.
%
% Tom Okell, May 2023
% 
% [M, t, P, G, RF] = test_CASL_seq(v, zmax, G_amp, RF_amp, ORFreq, dt, T1,T2,NegGz,LHR)

function [M, t, P, G, RF] = test_CASL_seq(v, zmax, G_amp, RF_amp, ORFreq, dt, T1,T2,NegGz,LHR)

if nargin < 9; NegGz = false; end
if nargin < 10; LHR = false; end

% Calculate the required time for the simulation
T = 2*zmax / v;

% Round up T to the nearest dt and calculate the number of time steps
T = T - mod(T,dt) + dt;
N = T/dt;
disp(['Number of time steps required is: ' num2str(N) ] );

% Generate the time array
disp('Generating time array...')
t = 0:dt:T;

% Set up the position array of the spin: starting at -zmax and ending at +zmax
disp('Generating position array...')
P = generate_position('const_vel', [0 0 -zmax]', [0 0 v]', t);

% Generate the gradient and RF arrays
disp('Generating gradient and RF arrays...')
[G, RF] = make_CASL_seq(t, RF_amp, G_amp);

% Negate Gz if required
if NegGz; G(3,:) = -G(3,:); end

% Generate an Off-resonance array
OffRes = ones(1,size(t,2))*ORFreq;

% Run the simulation
disp('Starting simulation...')
Ms = [0 0 1]'; % Starting magnetisation
M = bloch_sim(P, G, RF, OffRes, t, T1, T2, Ms, LHR);
