% Average a quantity over a laminar flow profile.
%
% Tom Okell, May 2023
%
% X_av = AverageOverLaminarFlowProfile(v,X,v_av)
%
% This function averages quantity X, defined at the velocities found in v,
% over a laminar flow profile with average velocity v_av.  Note that max(v)
% should be greater than 2*v_av and v should go as close to zero as
% possible, otherwise values close to zero and 2*v_av will need to be
% extrapolated. 
 
function X_av = AverageOverLaminarFlowProfile(v,X,v_av)

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
% area of each circular area of artery.
X_av = sum( X_int .* r_norm * dr ) / sum( r_norm * dr );



