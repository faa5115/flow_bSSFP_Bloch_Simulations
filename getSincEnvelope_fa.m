function B1e = getSincEnvelope_fa(gamma, RFDur, NumOfTP, TBW, flipAngle)

%gamma - gyromagnetic ratio
%RFDur - the duration of the RF pulse you want.
%NumOfTP - number of time samples within the RF pulse. 
%TBW - time-bandwidth product
%flipAngle - the flip angle .
t = linspace(-RFDur/2,RFDur/2,NumOfTP); %msec.
ExcBW = TBW / RFDur; %Hz.   ----> FWHM of the slice profile.


t0 = 1/ExcBW; %s.  half width of the central pulse.
%to great approximation: bandwidth ~1/t0.

dt = RFDur/NumOfTP;
tempSincFunc = t0 * sin(pi*t/t0)./(pi*t); %temporary sinc function to determine the amplitude.
sumOfTempSinc = sum(tempSincFunc);
A = flipAngle/gamma/(sumOfTempSinc * dt); %T.  B1 field's amplitude.  
B1e = A * t0 * sin(pi*t/t0)./(pi*t); %envelope function.



