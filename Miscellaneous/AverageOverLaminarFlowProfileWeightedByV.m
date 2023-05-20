% Average over a laminar flow profile, weighting by velocity to give an
% average over the volume of fluid passing through a plane
% 
% Tom Okell, May 2023
%
% X_av = AverageOverLaminarFlowProfileWeightedByV(v,X,v_av)
% 
% This function averages quantity X, defined at the velocities found in v,
% over a laminar flow profile with average velocity v_av.  Note that max(v)
% should be greater than 2*v_av and v should go as close to zero as
% possible, otherwise values close to zero and 2*v_av will need to be
% extrapolated. This version also weights by velocity so the final result
% should give the average across a volume of fluid passing through the
% cross section of the vessel in a given period of time. Here we use a
% direct method to integrate across the laminar profile.

function X_av = AverageOverLaminarFlowProfileWeightedByV(v,X,v_av)

% Set up an array describing normalised radius
dr = 0.001;
r_norm = dr/2:dr:(1-dr/2);

% Calculate normalised velocities (i.e. max(v) = 1) at each radius
v_norm = 1 - r_norm.^2;

% Set up the velocity array for averaging purposes
% NB. The maximum velocity is twice the average velocity
v_int = 2 * v_av * v_norm;

% Calculate inversion efficiencies at each velocity, corresponding to each radius
% NB. Extrapolate using linear interpolation since can't always simulate
% down to zero velocity
X_int = interp1(v,X,v_int,'linear','extrap');

% Calculate the average over the whole area of the artery, weighting by the
% area of each circular area of artery and the velocity.
X_av = sum( X_int .* r_norm * dr .* v_int) / sum( r_norm * dr .* v_int);

% We use a simple sum here but almost identical results come from using
% e.g. trapezoidal integration.

% NB. This gives almost identical results to the Maccotta et al. NMR Biomed
% 1997 method (see AverageOverLaminarFlowProfileWeightedByVMaccottaMethod.m)

