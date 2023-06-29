% This script sets default values for an SSFP type readout

%% Declare gamma
g = GetGamma; % Gyromagnetic ratio of hydrogen

%% Define T1 and T2
T1 = 1.65; T2 = 0.150; % As per oxford_asl

%% Spin parameters
v = 0.0; % Velocity (m/s): assume zero here to avoid the spin 
         % moving outside the imaging region, but could be modified 
         % to investigate flow sensitivity
Ms = [0 0 1]'; % Initial magnetisation state (relative to M0)
Ps = [0 0 0]'; % Initial position
ORFreq = 0; % Off-resonance frequency (Hz)

%% Sequence parameters
G_amp = [2 10]; % Gradient amplitude during the [RF readout] (mT/m)
FlipAng = 40; % Excitation flip angle (degs)
RF_dur = 1e-3; % RF pulse duration (s)
TR = 18e-3; % Repetition time (s)
Spoil = false; % Use spoiling (false for bSSFP)
G_s = 200; % Gradient slew rate (T/m/s = mT/m/ms)
delta = 0.2; % Excitation slab thickness (m)
ADC_dur = TR/3; % ADC duration (s)

%% Simulation parameters
T = 2; % Total simulation time (s)
dt = 0.02e-3; % Time step for simulation (s)


