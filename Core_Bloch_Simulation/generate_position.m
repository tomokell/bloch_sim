% Generates a position waveform for Bloch simulation code.
%
% Tom Okell, May 2023
%
% P = generate_position(move_type, start_pos, params, t)
%
% Input parameters are the movement type (move_type [string] = 'static' or
% 'const_vel'), starting position (start_pos [3x1]), parameters (params,
% e.g. velocity vector for the constant speed option) and an array of time
% points, t (in s)

function P = generate_position(move_type, start_pos, params, t)

% Check which movement type is required and call the appropriate function
switch move_type
    case 'static'
        P = repmat(start_pos, 1, size(t,2));

    case 'const_vel'
        P = const_vel(start_pos, params, t);

    otherwise
        error('Movement type not recognised!')        

end
