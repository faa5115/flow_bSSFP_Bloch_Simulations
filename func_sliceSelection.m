function M = func_sliceSelection (B1e, sliceDir, RFBW, sliceThickness, RFDur, M0, gamma, sliceOffCenter, RFPhase0, Meq, T1, T2, dfArr )

%inputes: 
%
%           B1e:             RF's envelope function.  Tesla.
%           RFPhase:         Phase of RF pulse.  Radians.
%           sliceDir:        a 1D array that spans space along the ss dir.  mm.
%           RFBW:              excitation bandwidth.  Hz.
%           sliceThickness:  the slice thickness.  mm.
%           RFDur:           RF pulse's duration.  s.
%           gamma:           Gyromagnetic ratio. rad/s/T.
%           M0:              num of isochromats along slice x 3 (mx, my,
%                            mz)-->stores the initial M for each isochromat 
%                            along the ss direction. 
%           Meq:             length(sliceDir) x 1 vector telling you the
%                            equilibrium magnetization values. 
%           dfArr:           length(sliceDir) x 1 .an array of off resonance frequencies along
%                            the slice select direction. Hz.

NumOfTP = length(B1e);
tRF = linspace(-RFDur/2,RFDur/2,NumOfTP); %s.
dt = RFDur / length(B1e);

%gamma = 2*pi*42.577*10^6; %rad/T. Gyromagnetic ratio.

Gslice = RFBW/(gamma*sliceThickness/(2*pi)); %T/mm

%dB --> off resonance felt along each location.
%       its 1x(Number of points along slice direction)
dB = Gslice * sliceDir + dfArr/ (gamma/ (2*pi)); %T




B1ea = B1e;
for n = 1 : size(B1e,2)
    B1ea(1,n) = B1e(1,n) * (cos(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n))) + 1i*(sin(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n)))));
end

B1e = B1ea ;
clear B1ea;

%our magnetization vector.
%       [slice position, logical axes, time]
M = zeros(length(sliceDir), 3, length(tRF));
for k = 1 :length(sliceDir)
    M(k, :, 1) = squeeze(M0(k, :));
end

%both of these are 1x(Number of time points of RF duration).
%they tell us the component of the B1 field in x and y respectively.
%later on...we can add concamitant field stuff...but not now.
vecToRotateAboutX = real(B1e) ;
vecToRotateAboutY = imag(B1e) ;



%iterate through each time point

%take into account the decay/recovery.
%we don't need to include off resonance here..because it's taken care of in
%the excitation process. 
[A, B] = freeprecess(dt, T1, T2, 0);
for iter_time = 2: length(tRF)
    %we update this for each timepoint.  
    unitVecToRotateAbout = zeros(length(sliceDir),3);
    amountToRotate = ones( length(sliceDir), 1);
    
    for sliceIter = 1: length(sliceDir)
        amountToRotate(sliceIter,1) = gamma * sqrt((vecToRotateAboutX(1, iter_time))^2 + (vecToRotateAboutY(1, iter_time))^2 + (dB(sliceIter,1))^2) * dt*180/pi;
        unitVecToRotateAbout(sliceIter, 1) = squeeze(vecToRotateAboutX(1,iter_time))/(squeeze(sqrt( (vecToRotateAboutX(1, iter_time))^2 + (vecToRotateAboutY(1, iter_time))^2 + (squeeze(dB(sliceIter,1)))^2)));
        unitVecToRotateAbout(sliceIter, 2) = squeeze(vecToRotateAboutY(1,iter_time))/(squeeze(sqrt( (vecToRotateAboutX(1, iter_time))^2 + (vecToRotateAboutY(1, iter_time))^2 + (squeeze(dB(sliceIter,1)))^2)));
        unitVecToRotateAbout(sliceIter, 3) = squeeze(dB(sliceIter,1))/(squeeze(sqrt( (vecToRotateAboutX(1, iter_time))^2 + (vecToRotateAboutY(1, iter_time))^2 + (squeeze(dB(sliceIter,1)))^2)));
        
        M(sliceIter,:,iter_time) = rotVecAroundArbAxis(squeeze(M(sliceIter,:,iter_time-1)),squeeze(unitVecToRotateAbout(sliceIter, :) ),squeeze(amountToRotate(sliceIter,1)));
        M(sliceIter,:,iter_time) = A * squeeze(M(sliceIter, :, iter_time)') + B * squeeze(Meq(sliceIter));
    end
    
    %M(:,:,iter_time) = rotVecAroundArbAxis(squeeze(M(:,:,iter_time-1)),unitVecToRotateAbout,amountToRotate);
    
    %now take decay/recovery into account.
    %{
    for sliceIter = 1 : length(sliceDir)
        %the df part of freeprecess should always be zero here...because we
        %took care of the off resonance above. 
        [A, B] = freeprecess((iter_time - 1) * dt, T1, T2, 0);
        M(sliceIter, :, iter_time) = A * squeeze(M(sliceIter, :, iter_time)') + B * squeeze(Meq(sliceIter));
    end
    %}
end