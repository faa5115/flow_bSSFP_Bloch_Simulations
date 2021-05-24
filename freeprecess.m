function [Afp,Bfp]=freeprecess(T,T1,T2,df)

phi = 2 * pi*df * T;

flip = pi/3; %radians
FlipY = yrot(flip);

Afp = [exp(-T/T2) 0 0; 0 exp(-T/T2) 0; 0 0 exp(-T/T1)] * zrot(phi);
Bfp = [0 0 (1 - exp(-T/T1))]';

