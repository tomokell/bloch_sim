%  This function returns a 3D rotation matrix about the specified axis for
%  the rotation angle given in radians
%
%  Tom Okell, May 2023
%
%  R = RotMat3D(axis, angle)

function R = RotMat3D(axis, angle)

% Check that axis is non-zero
if sum(axis.^2) == 0  % Axis is zero, so return the identity matrix
    R = eye(3);
else

    % Check that axis is in the right format and transpose if necessary
    if size(axis) == [1 3]  % Transpose axis to turn into a column vector
        axis = axis';

    elseif size(axis) == [3 1]  % Correct format so do nothing

    else  % Incorrect format - reject
        error(['Incorrect vector format of size ' ns(size(axis)) ' given'])
    end

    % If axis is not a unit vector, convert it to a unit vector
    mag = sqrt( sum(axis.^2) );

    if mag ~= 1
        axis = axis/mag;
    end

    % Use symbols for the axis for simplicity
    x = axis(1);  y = axis(2); z = axis(3);

    % Return the general rotation matrix
    R = [0 -z y ; z 0 -x ; -y x 0 ] * sin(angle) + (eye(3) - axis * axis') * cos(angle) + axis * axis';

end