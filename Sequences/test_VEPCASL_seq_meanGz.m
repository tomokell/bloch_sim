% Tests the adiabatic inversion caused by the (VE)PCASL sequence.
% 
% Tom Okell, May 2023
% 
% function [M, t, P, G, RF] = test_VEPCASL_seq_meanGz(v, zmax, z_offset, meanGz, G_amp, RF_shape, RF_shape_params, RF_amp, RF_dur, RF_sep, Pa, Pb, cycle, Ps, dt, T1, T2, ORFreq)
% 
% Runs a Bloch simulation using a (VE)PCASL sequence for a single spin.
% Required inputs:
%   v               =   spin speed (m/s)
%   zmax            =   maximum z location to simulate (m)
%   z_offset        =   z offset of the labelling plane from isocentre (m):
%                       WARNING: Non-zero values not well tested...
%   meanGz          =   Mean labelling gradient (mT/m)
%   G_amp           =   Gradient amplitude during RF pulses (mT/m)
%   RF_shape        =   'hard': Rectangular (hard) pulse - no parameters
%                       'gaussian': Gaussian pulse - 1 parameter - temporal
%                                   SD (s)
%                       'gaussianhann': As for Gaussian, but multiplied by a Hann window
%                       'sinc': Sinc pulse - 1 parameter - slice width (m)
%                       'hanning': Hanning pulse - no parameters
%                       'hamming': Hamming pulse - no parameters
%                       'gaussian_dsv': Gaussian pulse from a DSV file - no
%                                       parameters
%   RF_shape_params =   As needed for each RF pulse (see above)
%   RF_amp          =   RF peak amplitude (mT)
%   RF_dur          =   RF pulse duration (s)
%   RF_sep          =   RF pulse interval (s)
%   Pa              =   2D position of the "vessel A" within the labelling
%                       plane (m). Only relevant for vessel-encoding. See Wong MRM 2007.
%   Pb              =   2D position of the "vessel B" within the labelling
%                       plane (m). Only relevant for vessel-encoding. See Wong MRM 2007.
%   cycle           =   (VE)PCASL cycle number:
%                           1: Label all vessels (standard PCASL label)
%                           2: Control all vessels (standard PCASL control)
%                           3: Label at location A, control at location B
%                           4: Label at location B, control at location A
%   Ps              =   2D location within the labelling plane of the
%                       simulated spin (m)
%   dt              =   Time step (s)
%   T1              =   Longitudinal relaxation time (s)
%   T2              =   Transverse relaxation time (s)
%   ORFreq          =   Off-resonance frequency (Hz)
%
% Optional inputs:
%   Unipolar        =   true: Unipolar VEPCASL; false: bipolar VEPCASL
%   NegGz           =   Negate the z gradient
%   ang             =   Angle of the spin relative to the normal vector of
%                       the labelling plane
%   xory            =   Angle the spin so that it moves in the xz plane
%                       ('x') or yz plane ('y')
%
% Outputs: 
%   M               =   Magnetisation vector at each time point (3xNt),
%                       normalised to M0 = 1
%   t               =   Simulated time points (1xNt) (s)
%   P               =   Position vector at each time point (3xNt) (m)
%   G               =   Gradient at each time point (3xNt) (mT/m)
%   RF              =   RF waveform at each time point ([amplitude in mT/phase in rads],Nt)
%   T               =   Total simulation time (s)

function [M, t, P, G, RF, T] = test_VEPCASL_seq_meanGz(v, zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                           RF_shape_params, RF_amp, RF_dur, ...
                                           RF_sep, Pa, Pb, cycle, Ps, dt, T1, ...
                                           T2, ORFreq, Unipolar, NegGz, ang, xory)
                                       
if nargin < 19; Unipolar = false; end
if nargin < 20; NegGz = false; end
if nargin < 21; ang = 0; end
if nargin < 22; xory = 'x'; end

% Calculate the velocity vector of the spin
if xory == 'x'  % Orient velocity vector in the xz plane
  vel = v * [sin(ang * pi/180) 0 cos(ang*pi/180)]';
else  % Orient the velocity vector in the yz plane
  vel = v * [0 sin(ang * pi/180) cos(ang*pi/180)]';
end

% Calculate the required time for the simulation
T = 2*zmax / vel(3);

% Round up T to the nearest dt and calculate the number of time steps
T = T - mod(T,dt) + dt;
N = T/dt;
disp(['Number of time steps required is: ' num2str(N) ] );

% Generate the time array
t = 0:dt:T;

% Generate the off resonance frequency array
OffRes = ones(1,size(t,2))*ORFreq;

% Set up the position array of the spin: starting at -zmax and ending at
% +zmax, going through [Ps(1) Ps(2) 0] at t = T/2
P = generate_position('const_vel', [(Ps(1) - vel(1)*T/2) (Ps(2)- vel(2)*T/2) -zmax]', vel, t);

% Generate the gradient and RF arrays
[G, RF] = make_VEPCASL_seq_meanGz_Gamp(RF_shape, RF_shape_params, RF_amp, ...
                                   RF_dur, RF_sep, dt, T, meanGz, G_amp, ...
                                   Pa, Pb, z_offset, cycle);


% Modify to make unipolar if requested
if Unipolar
    G(1:2,:) = abs(G(1:2,:));
end

% Modify to make Gz negative if requested
if NegGz
    G(3,:) = -G(3,:);
end

% Run the simulation
M = bloch_sim(P, G, RF, OffRes, t, T1, T2);
