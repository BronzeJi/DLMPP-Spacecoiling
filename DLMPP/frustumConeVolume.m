function volume = frustumConeVolume(radius1, radius2, height)
%this function calculate the volume of a frustumcone
    % Calculate areas of the bases
    A1 = pi * radius1^2;
    A2 = pi * radius2^2;

    % Calculate frustum volume
    volume = (height / 3) * (A1 + A2 + sqrt(A1 * A2));
end

