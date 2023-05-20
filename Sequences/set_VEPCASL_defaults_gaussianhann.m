%% Sets default values for the parameters required to simulate the (VE)PCASL sequence
% 
% Tom Okell, May 2023

%% Declare gamma
g = GetGamma;  % rad/s/mT                    

%% Define the velocity and positions
v = 20 / 100; %  Spin velocity, cm/s to m/s conversion
zmax = 3 / 100; % Maximum z location away from the labelling plane to simulate, cm to m conversion
Pa = [-0.05 0]'; % Position of vessel "A" within the labelling plane (see Wong MRM 2007)
Pb = [0 0]'; % Position of vessel "B" within the labelling plane
Ps = Pa; % Position of the spin to simulate within the labelling plane
z_offset = 0; % z offset of the labellign plane (assume isocentre - not well tested off-isocentre)

%% Define the off resonance frequency
ORFreq = 0;

%% Define VEPCASL cycle number
cycle = 3; % Label at vessel A and control at vessel B (see Wong MRM 2007)
% NB. 1 and 2 can be used for conventional PCASL tag and control,
% respectively

%% Define the gradient parameters
% Mean labelling gradient
meanGz = 0.8; % mT/m
                                   
% Gradient slew rate
G_s = 150;  % mT/m/ms

% Gradient amplitude during RF pulses                               
G_amp = 6; % mT/m

%% Define T1 and T2 of blood
T1 = 1.932; T2 = 0.275;  % @ 3T, ref: Stanisz, MRM 2005

%% Define RF parameters for the PCASL sequence
RF_amp = 0.04 * 1e-4 * 1000; % Gauss to mT conversion
                                  % Wong = 0.04 G
RF_shape = 'gaussianhann';  % Hann windowed Gaussian pulse
RF_dur = 600e-6;  % RF pulse duration, in s
RF_sep = 960e-6;  % RF pulse repetition time (interval), in s
RF_shape_params = RF_dur/1.5; % Gives result very similar to standard Okell MRM 2010 implementation

%% Define dt to sample RF adequately
dt = 10e-6;

%% Rescale the RF amplitude to get 20 degree flip angle pulses
Desired_FlipAngle_degs = 20;
Act_FlipAngle_degs = FlipAngle(gaussian_RF(RF_amp,RF_shape_params(1),0,RF_dur,dt,true),dt);
RF_amp = RF_amp * Desired_FlipAngle_degs / Act_FlipAngle_degs;

%% Define bipolar vs. unipolar VEPCASL
Unipolar = false; % Bipolar

%% Define the use of negative Gz waveform
NegGz = false;

%% Define the angle of the spin's motion relative to the normal vector of the labelling plane
ang = 0;
xory = 'x';