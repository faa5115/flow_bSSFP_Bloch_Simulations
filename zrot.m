function Rz=zrot(phi)
%precession about z axis.

%phi is in radians.

Rz = [cos(phi) -sin(phi) 0 ; sin(phi) cos(phi) 0; 0 0 1];

