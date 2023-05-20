% Apply an effective magnetic field to rotate the magnetisation
%
% Tom Okell, May 2023
%
% [M_new, R] = Apply_Beff(M, Beff, dt, g, LHR)
%
% This function takes a magnetisation vector (M), effective magnetic field
% (Beff), time step (dt) and gyromagnetic ratio (g), and rotates the
% magnetisation vector by the required amount, returning both the rotation
% matrix (R) and the new magnetisation vector (M_new). If LHR is set to
% "true" (the default), then left-hand rule rotations are applied.

function [M_new, R] = Apply_Beff(M, Beff, dt, g, LHR)

if nargin < 5; LHR = true; end

% Determine the angle through which M should be rotated about Beff
% omega/rad per s = gamma * B so theta = g * B * dt
th = g * sqrt(sum(Beff.^2)) * dt;

% By default, the rotation matrix below calculates right-handed rotations,
% so negate the rotation angle for left-hand rule rotations
if LHR
    th = -th;
end

% Calculate the rotation matrix
R = RotMat3D(Beff, th);

% Apply the rotation matrix to M
M_new = R * M;