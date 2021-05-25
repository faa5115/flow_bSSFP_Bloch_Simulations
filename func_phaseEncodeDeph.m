function MPE = func_phaseEncodeDeph (FOV, M, lengthPE, stepNum, dfThetaArr)
%dfThetaArr should have the same length as lengthPE.  It's in radians. 

%--------------------------------------------------------------------------
%PEDur is in seconds.  we really could just calculate the intravoxel dephasing
%appropriately without it, but i want to include what effects off-resonance
%has...and for that we'll need to the amount of time. 

%nevermind i removed it.
%--------------------------------------------------------------------------

dK  = 1/FOV; %1/mm. 

%so.... we will be taking an "image" of our excitation profile.  
%let's say we want our image to have say...NPixelz points...




%just imagination.  we sample this "excited" data.  but we don't know
%what's within the rang of our actual FOV or not. below is if we sampled
%appropriately with the above sampling rate.  we'd need a FOV of 
% size(sampledImageData, 1) * dz.  
%sampledImageData = MRef( 1:sampleRate:size(MRef,1) ,1,length(tRef));
%sampledImageData = MRef( 1:sampleRate:size(M,1) ,:,length(tRef));


%sampledImageLocation = lengthPE(1:sampleRate:size(lengthPE,1));

rotAngle = zeros(length(lengthPE), 1);
MPE = zeros(length(lengthPE),3);

for n = 1: length(lengthPE)
    rotAngle(n, 1) = squeeze(lengthPE(n, 1)) * stepNum * dK * 2 * pi + squeeze(dfThetaArr(n, 1));
    
    %rotMatrix = zrot(-squeeze( rotAngle(n, 1)) );
    rotMatrix = zrot(squeeze( rotAngle(n, 1)) );
    
    MPE(n, : ) = rotMatrix * [squeeze( M(n, 1)); squeeze( M(n, 2)); squeeze( M(n, 3))]  ;
end










