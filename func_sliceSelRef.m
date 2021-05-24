function MRef = func_sliceSelRef (RFDur, RFBW, sliceThickness,  B1e, MLast, sliceDir, Meq, T1, T2, dfArr)


gamma = 2*pi*42.577*10^6; %rad/T. Gyromagnetic ratio.
Gslice = RFBW/(gamma*sliceThickness/(2*pi)); %T/mm
GsliceRef = -Gslice;
dBRef = GsliceRef * sliceDir + dfArr/(gamma/(2*pi)); %T
slRefDur = (RFDur/2);
NumOfTPRef = round(length(B1e)/2) ;
tRef = linspace(-slRefDur/2,slRefDur/2,NumOfTPRef); %s.
dtRef = slRefDur/NumOfTPRef; %s

MRef = zeros(length(sliceDir), 3, length(tRef));


%{
for k = 1 :length(sliceDir)
    MRef(k, :, 1) = MLast(k, :);
end
%}

%both of these are 1x(Number of time points of RF duration).
%Now we have nothing applied along the  x and y respectively.
%later on...we can add concamitant field stuff...but not now.
vecToRotateAboutX = zeros(1,length(tRef));
vecToRotateAboutY = zeros(1,length(tRef));

%take into account the decay/recovery.
%we don't need to include off resonance here..because it's taken care of in
%the excitation process. 
[A, B] = freeprecess(dtRef, T1, T2, 0);


%iterate through each time point
for iter_time = 1: length(tRef)
    %we update this for each timepoint.  
    unitVecToRotateAbout = zeros(length(sliceDir),3);
    amountToRotate = ones( length(sliceDir), 1);
    
    for sliceIter = 1: length(sliceDir)
        amountToRotate(sliceIter,1) = gamma * dBRef(sliceIter,1) * dtRef*180/pi;
        unitVecToRotateAbout(sliceIter, 1) = squeeze(vecToRotateAboutX(1,iter_time))/(squeeze(dBRef(sliceIter,1) ));
        unitVecToRotateAbout(sliceIter, 2) = squeeze(vecToRotateAboutY(1,iter_time))/(squeeze(dBRef(sliceIter,1) ));
        unitVecToRotateAbout(sliceIter, 3) = squeeze(dBRef(sliceIter,1))/(squeeze(dBRef(sliceIter,1)  ));
        if (amountToRotate(sliceIter,1) == 0)
            %doing this because there is no rotation ---> but that makes
            %this unit vector be zero...which will result in an error
            %below.  
            unitVecToRotateAbout(sliceIter, 1) = 0;
            unitVecToRotateAbout(sliceIter, 2) = 0;
            unitVecToRotateAbout(sliceIter, 3) = 1; 
        end
        
            if iter_time == 1
                MRef(sliceIter,:,iter_time) = rotVecAroundArbAxis(squeeze(MLast(sliceIter,:)),squeeze(unitVecToRotateAbout(sliceIter,:)),squeeze(amountToRotate(sliceIter, :)));
            else
                MRef(sliceIter,:,iter_time) = rotVecAroundArbAxis(squeeze(MRef(sliceIter,:,iter_time-1)),squeeze(unitVecToRotateAbout(sliceIter, :)),squeeze(amountToRotate(sliceIter, :)));
            end 
            
            MRef(sliceIter,:,iter_time) = A * squeeze(MRef(sliceIter, :, iter_time)') + B * squeeze(Meq(sliceIter));
    end
    %{
    if iter_time == 1
        MRef(:,:,iter_time) = rotVecAroundArbAxis(squeeze(MLast(:,:)),unitVecToRotateAbout,amountToRotate);
    else
        MRef(:,:,iter_time) = rotVecAroundArbAxis(squeeze(MRef(:,:,iter_time-1)),unitVecToRotateAbout,amountToRotate);
    end 
    %}
end