% General Bloch simulation function.
%
% Tom Okell, May 2023
%
% M = bloch_sim(P, G, RF, OffRes, t, T1, T2, Ms, LHR, Spoil)
% 
% Runs the Bloch Simulation, returning the magnetisation vector (M) at each
% time point ([Mx,My,Mz],t) corresponding to the input data: spin position,
% P([x,y,z],t), gradients G([Gx,Gy,Gz],t), and RF([amplitude/mT,
% phase/rads],t), resonance offset, OffRes/Hz (t), time points, t,
% relaxation times (T1 and T2), starting magnetisation state, Ms (default
% is equilibrium = [0 0 1]'), whether to use the left-hand rule for
% rotations (default = true). If the vector Spoil is specified, the
% transverse magnetisation is perfectly spoiled (set to zero) whenever
% Spoil is true.

function M = bloch_sim(P, G, RF, OffRes, t, T1, T2, Ms, LHR, Spoil)

if nargin < 8; Ms = [0 0 1]'; end
if nargin < 9; LHR = true; end
if nargin < 10;Spoil = false(size(t)); end

if size(Ms,1) == 1; Ms = Ms'; end

% Declare gamma
g = GetGamma;  % rad/s/mT

% Initialise the magnetisation vector
M = zeros(3,size(t,2));
M(:,1) = Ms;

% Determine the step size in time (assumed constant)
dt = t(2) - t(1);

% Determine T1 and T2 decay factors for this time interval
T1_decay_factor = exp( -dt/T1 );
T2_decay_factor = exp( -dt/T2 );

% Run through time steps, determining rotation matrix each time.  Apply and
% save as output, but keep running combined matrix and compare final
% results
for ii = 1:(size(t,2)-1)
  
  % Calculate the effective B field at this time point
  Beff = Calc_Beff( P(:,ii), G(:,ii), RF(:,ii), OffRes(ii), g );
  
  % Apply Beff for the time period dt
  M(:,ii+1) = Apply_Beff( M(:,ii), Beff, dt, g, LHR);
  
  % Apply T1 and T2 decay
  M(1,ii+1) = M(1,ii+1) * T2_decay_factor;
  M(2,ii+1) = M(2,ii+1) * T2_decay_factor;
  M(3,ii+1) = 1 + ( M(3,ii+1) - 1 ) * T1_decay_factor;

  % Apply spoiling (set transverse magnetisation to zero)
  if Spoil(ii) == true
    M(1:2,ii+1) = 0;
  end  
end
