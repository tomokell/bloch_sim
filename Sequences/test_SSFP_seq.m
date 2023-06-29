% Simulates the effect of an SSFP/SPGR sequence on spins with or without
% spoiling.
%
% Tom Okell, June 2023
%
% [M, t, P, G, RF, TEIdx, ADC] = test_SSFP_seq(v, Ms, G_amp, FlipAng, T, RF_dur, TR, Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses, LHR, VFAType, G_s, delta, ADC_dur)
%
% This function runs a Bloch simulation of the magnetisation during an
% SSFP/SPGR pulse sequence, returning the magnetisation, M, at time points
% t, position array, P, gradient and RF waveforms, G and RF, an index of
% times when the TE occurs, TEIdx, and a logical array showing when the ADC
% is turned on, ADC.
%
% Inputs: many sequence parameters (FlipAng, RF_dur, TR, G_amp, SeqType,
% NoLinIncPulses, VFAType, G_s, delta, ADC_dur) are defined in
% make_SSFP_seq.m. Other parameters:
%   v               =   Spin velocity (m/s) in the z direction
%   Ms              =   Starting magnetisation vector (relative to M0)
%   T               =   Total sequence duration (s)
%   Ps              =   Initial spin position [x y z] (m)
%   dt              =   Time increment for waveform sampling (s)
%   T1              =   Longitudinal relaxation time (s)
%   T2              =   Transverse relaxation time (s)
%   ORFreq          =   Off-resonance frequency (Hz)
%   LHR             =   Use left-hand rule rotations for Bloch simulation

function [M, t, P, G, RF, TEIdx, ADC] = test_SSFP_seq(v, Ms, G_amp, FlipAng, T, RF_dur, TR, ...
                                           Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses, ...
                                           LHR, VFAType, G_s, delta, ADC_dur)

if nargin < 13; SeqType = 'HalfAngle';  end
if nargin < 14; NoLinIncPulses = 10;    end
if nargin < 15; LHR = true;             end
if nargin < 16; VFAType = 'CFA';        end
if nargin < 17; G_s = 200;              end
if nargin < 18; delta = 0.2;            end
if nargin < 19; ADC_dur = TR/3;         end

% Round up T and TR to the nearest dt and calculate the number of time steps
if mod(T,dt) ~=0;  T  = T  - mod(T,dt)  + dt; end
if mod(TR,dt) ~=0; TR = TR - mod(TR,dt) + dt; end
N = T/dt + 1;

disp(['Number of time steps required is: ' num2str(N) ] );

% Generate the time array
t = 0:dt:T;

% Generate the off resonance frequency array
if numel(ORFreq) == 1 % Same off-resonance for the whole simulation
    OffRes = ones(1,size(t,2))*ORFreq;
else  % Different off-resonance over time
    OffRes = ORFreq(1:length(t));
end

% Set up the position array of the spin
P = generate_position('const_vel', Ps(:), [0 0 v]', t);

% Generate the gradient and RF arrays 
[G, RF, Spoil, TEIdx, ~, ADC] = make_SSFP_seq(FlipAng, RF_dur, TR, dt, T, G_amp, ...
                                              SeqType, NoLinIncPulses, VFAType, ...
                                              G_s, delta, ADC_dur);

% Run the simulation
M = bloch_sim(P, G, RF, OffRes, t, T1, T2, Ms, LHR, Spoil);
