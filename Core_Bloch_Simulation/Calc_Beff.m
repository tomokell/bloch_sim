% Calculate the effective magnetic field
%
% Tom Okell, May 2023
%
% Beff = Calc_Beff(Pt, Gt, RFt, OffRes, g)
%
% This function calculates the effective magnetic field, Beff, in the
% rotating frame, given the spin position (Pt), gradients (Gt), RF (RFt) and off
% resonance frequency/Hz (OffRes) at a single time point.  Note that a negative
% OffRes value correctly translates to a negative Beff contribution. 
% The gyromagnetic ratio, Gamma (g), is also required.

function Beff = Calc_Beff(Pt, Gt, RFt, OffRes, g)

% The effective B field is a combination of precession about the z axis due 
% to the gradients, the RF and the off-resonance frequency.  In the frame rotating at the Larmor
% frequency, W0 = gamma*B0, so we have Beff = (0,0,G.P) + B1 + (0,0,2*pi*OffRes/gamma)

% B1 must be calculated from the amplitude, RFt(1), and phase, RFt(2), of
% the RF at this point in time, with zero phase corresponding to the x'
% axis.  Note that we have used the convention that increasing RF phase
% rotates in a right hand sense about z, which would need to be modified
% if using the standard left hand sense convention for precession of the
% magnetisation about z.
B1 = [ RFt(1)*cos( RFt(2) ); RFt(1)*sin( RFt(2) ); 0];

% Now return Beff
Beff = [0; 0; Gt'*Pt] + B1 + [0; 0; 2*pi*OffRes/g];