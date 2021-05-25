%This file gives Bloch simulation results for an array for spatially
%arranged set of spins that shift down by a fractional spin replacement
%rate per TR, dS.  The simulation is done for varying dS and off-resonances
%(for rotations from -2 * pi rad/TR to 2 * pi rad/TR) after 300 excitations.   
clear;

%you want to control the duration, the thickness, excitation profile...
RFDur = 800e-6; %s %original
gamma = 2*pi*42.577*10^6; %rad/s/T. Gyromagnetic ratio.
TBW = 3;%original
%TBW = 4; 
NumOfTPRF = 128; %for the RF profile. 
tRF = linspace(-RFDur/2,RFDur/2,NumOfTPRF); %s.
sliceThickness = 6.00; %mm
flipAngle = 60 * pi/180; %radians. 
%RFPhase = pi/2; %radians.
RFPhase = 0; %radians.
%sliceFOV = 30; %mm.
%IntElementsPerMM = 10;
%let's say we want Ns subslices...


Ns = 20;  %Ns must be an integer multiple of sliceThickness.
%so now...our IntElementsPerMM is: 
IntElementsPerMM = Ns/sliceThickness;
subSliceThickness = sliceThickness/Ns; %mm/sub slice

TR = 4e-3; %s
TE = TR/2; %s

T2 = 150e-3;  %s
T1 = 1000e-3;  %s

%T1 = 500e-3; %s
%T2 = 500e-3; %s

%T2Fac = 4; %unitless.  just saying after how many T2s will we  begin ignoring signal. 
T2Fac = 3;

%dS is going to be a parameter that defines how quickly the spins flow
%between TRs.  dSMin is always going to be one element shift / TR
%therefore it'll always be 1/Ns.  
%otherwise ... dS will be zero (stationary spins).  
%dS * Ns is the number of elements an isochromat shifted between TRs.  
dSMin = 1/Ns; %smallest nonzero spin exchange percentage.  always going to be 1/Ns.

%NPart = 6;
NPart = 1; %Don't worry about it. 
NEffecgtiveFOVElements = Ns * NPart;
%dS = 0.3; %fraction. 
%dS = 0.2; %example
%dS = 0;
%dS = 1200 * TR/sliceThickness;
%dS = 1.0;
dS = 0.16;
NShift = round(dS * Ns);

NOs =  round(T2Fac * dS * T2 * Ns /TR);
%sliceFOV = (Ns + NOs)*subSliceThickness; %mm
sliceFOV = (NEffecgtiveFOVElements + 2* NOs)*subSliceThickness; %mm


sliceDir = (linspace(-sliceFOV/2, sliceFOV/2, IntElementsPerMM * sliceFOV))'; %mm

%intElementsOnSlice = 100; %number of elements along one axis of a slice.  


RFBW = TBW/RFDur; %Hz %this is the RF's bandwidth. ...
Gslice = RFBW/(gamma*sliceThickness/(2*pi)); %T/mm

B1e = getSincEnvelope_fa(gamma, RFDur, NumOfTPRF, TBW, flipAngle);
%B1efirst = getSincEnvelope_fa(gamma, RFDur, NumOfTPRF, TBW, flipAngle/2);

sliceOffCenter = 0;
%RFPhase0 = pi/2;
RFPhase0 = pi;
%RFPhase0 = 0;
%RFPhase0 = pi;
NumOfExc = 300;
M0     = zeros(size(sliceDir,1), 3);
MbSSFP = zeros(size(sliceDir,1), 3, NumOfExc );

Meq = zeros(size(sliceDir,1), 1);
for n = 1 : size(sliceDir, 1)
    M0(n, :) = [0 0 1]';
    Meq(n)   = norm(squeeze(M0(n,:) ));
end
%{
B1ea = B1e;
B1eb = B1e;

for n = 1 : size(B1e,2)
    B1ea(1,n) = B1e(1,n) * (cos(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n))) + 1i*(sin(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n)))));
end

B1e = B1ea ; 
clear B1ea;
%}

phaseInc = pi;

NumOfOffRes = 75;
PhiPerTR = linspace(-2*pi, 2*pi, NumOfOffRes);
OffRes   = PhiPerTR/(2 * pi * TR);

offresonanceones = ones(size(sliceDir));
%%


MAllExcitationsEachOffRes = zeros(size(M0,1), size(M0,2), NumOfExc, NumOfOffRes);
MPrev = zeros(size(M0));

%NShift = 0;
%OffResIter = round(3*length(PhiPerTR)/4);
for OffResIter = 1 : NumOfOffRes

    %now prepare the decay matrices...
    df = squeeze(OffRes(1,OffResIter));
    [Atr, Btr] = freeprecess(TR,T1,T2,df);

    MSliceProfile0 = func_sliceSelection (B1e, sliceDir, RFBW, sliceThickness, RFDur, M0, gamma, sliceOffCenter, RFPhase0, Meq, T1, T2, df *  offresonanceones );

    %debug
    %{
    figure, 
    subplot(2, 1, 1)
    hold on
    plot(sliceDir,abs( squeeze(MSliceProfile0(:, 1, NumOfTPRF)) + ...
        1i* squeeze(MSliceProfile0(:, 2, NumOfTPRF))), 'Linewidth', 5.0);
    plot(sliceDir,angle( squeeze(MSliceProfile0(:, 1, NumOfTPRF)) + ...
        1i* squeeze(MSliceProfile0(:, 2, NumOfTPRF)))/ pi,'Linewidth', 5.0);
    
legend('Magnitude','Phase (Fraction of \pi)')
    xlabel('Length (mm)')
    title('M Magnitude and Phase Distribution Immediately After RF Pulse')
    hold off
    
    
     subplot(2, 1, 2)
    hold on
    plot(sliceDir,( squeeze(MSliceProfile0(:, 1, NumOfTPRF))), 'Linewidth', 5.0);
    plot(sliceDir,( squeeze(MSliceProfile0(:, 2, NumOfTPRF))), 'Linewidth', 5.0);
    legend('Real Component','Imaginary Component')
    xlabel('Length (mm)')
    title('Transverse Magnetization Immediately After RF Pulse')
    hold off
    %}
    
    MSliceProfile=func_sliceSelRef (RFDur, RFBW, sliceThickness,  B1e, squeeze(MSliceProfile0(:, :, size(MSliceProfile0,3))), sliceDir, Meq, T1, T2, df *  offresonanceones);
  
    %{
    figure, 
    subplot(2, 1, 1)
    hold on
    plot(sliceDir,abs( squeeze( MSliceProfile( :, 1, round(NumOfTPRF/2) ) ) + ...
        1i* squeeze(MSliceProfile( :, 2, round(NumOfTPRF/2) ))), 'Linewidth', 5.0);
    plot(sliceDir,angle( squeeze(MSliceProfile(:, 1, round(NumOfTPRF/2))) + ...
        1i* squeeze(MSliceProfile(:, 2, round(NumOfTPRF/2))))/ pi,'Linewidth', 5.0);
    
legend('Magnitude','Phase (Fraction of \pi)')
    xlabel('Length (mm)')
    title('M Magnitude and Phase Distribution ')
    hold off
    
    
     subplot(2, 1, 2)
    hold on
    plot(sliceDir,( squeeze(MSliceProfile(:, 1, round(NumOfTPRF/2)))), 'Linewidth', 5.0);
    plot(sliceDir,( squeeze(MSliceProfile(:, 2, round(NumOfTPRF/2)))), 'Linewidth', 5.0);
    legend('Real Component','Imaginary Component')
    xlabel('Length (mm)')
    title('Transverse Magnetization ')
    hold off
    %}
    
    %{
    figure, plot( sliceDir,angle (squeeze(MSliceProfile(:,1,size(MSliceProfile,3))) +1i* squeeze(MSliceProfile(:,2,size(MSliceProfile,3)))))

    figure, plot( sliceDir,abs ( squeeze(MSliceProfile(:,3,size(MSliceProfile,3)))))


    figure, cplot(  (squeeze(MSliceProfile(:,1,size(MSliceProfile,3))) +1i* squeeze(MSliceProfile(:,2,size(MSliceProfile,3))))  )
    %}

    % calculate the flip angle and the phase of each spin-ensemble.   --------------------------------------------

    localNetRFPhase = zeros(length(sliceDir), 1);
    localNetRFFA    = zeros(length(sliceDir), 1);

    %in the loop below ... you changed the signs and swapped the im/re
    %components because rf phase is 90 deg out of phase of the magnetization's
    %phase.  
    for n = 1: length(sliceDir)
        localNetRFPhase(n, 1) = angle(  squeeze(MSliceProfile(n,1,size(MSliceProfile,3))) + 1i *squeeze(MSliceProfile(n,2,size(MSliceProfile,3))));
        localNetRFFA   (n, 1) = angle( 1i*abs(squeeze(MSliceProfile(n,1,size(MSliceProfile,3))) +1i* squeeze(MSliceProfile(n,2,size(MSliceProfile,3)))) +  (squeeze(MSliceProfile(n,3,size(MSliceProfile,3)))  )     );
    end
    
    %{
    figure, subplot(3, 1, 1), plot(sliceDir, lphase1 * 180 /pi, 'Linewidth',5.0), xlabel('Length (mm)'), ylabel('Transverse Magnetization Phase (Degrees)'), title('Odd Excitations'),  
    subplot(3, 1, 2), plot(sliceDir, localNetRFPhase * 180 /pi, 'Linewidth',5.0), xlabel('Length (mm)'), ylabel('Transverse Magnetization Phase (Degrees)'), title('Even Excitations'), 
    subplot(3, 1, 3), plot(sliceDir, localNetRFFA * 180 /pi, 'Linewidth',2.0), xlabel('Length (mm)'), ylabel('Flip Angle (Degrees)'), title('Local Flip Angle')
    %}
    
    % calculate the proper rotation matrix for each isochromat.  
    estExcMat = zeros(3, 3, length(sliceDir)); %estimated excitation matrix. 
    for n = 1 : length(sliceDir)
        estExcMat(:, :, n) = throt(squeeze ( localNetRFFA(n, 1) ), squeeze( localNetRFPhase(n, 1) ) );
    end

    localNetRFPhase0 = localNetRFPhase;
    localNetRFFA0 = localNetRFFA;
    estExcMat0 = estExcMat;
    % loop over several off resonances and isochromats.  

    % end calculate the flip angle and the phase of each isochromat.   -----------------------------------------
    %-------------------------- -------------------------- -------------------------- --------------------------  

    %the slice select function outputs a Num along slice x 3 x num of rf time
    %points dimension matrix, detailing the behavior of each isochromat along
    %the slice, for each time point that the rf pulse is active. 

    %but remember...this is *balanced* ssfp.  so you're going to want a slice
    %prephaser.  so... you can use your func_sliceSelRef function for
    %prephasing and refocusing. 

    %------------------------CHANGE FOR II VERSION OF THIS CODE----------------
    %    instead of using the entire gradient waveform calculation for the
    %    excitation profile of bssfp, we will use the rotation matrices above.
    %    we will add pi radieanse to each phase too.  
    %--------------------------------------------------------------------------

    %okay...let's iterate through...

    MAllExcitations = zeros(size(M0,1), size(M0,2), NumOfExc);



        localNetRFPhase = localNetRFPhase0;
        localNetRFFA = localNetRFFA0;
        estExcMat = estExcMat0;
        MNextCycle = M0;




        for excInc = 1 : NumOfExc

            if (excInc == 1)
                MPrev = M0;
            else
                MPrev = squeeze(MAllExcitationsEachOffRes(:, :, excInc-1, OffResIter));
            end

            %right here is a good spot to "move" the isochromats...
            %MPrev2 = MPrev;
            MPrev2 = shiftFromFlow(MPrev, NShift);


            for n = 1 : length(sliceDir)

                %prepare the rotation excitation matrics...

                estExcMat(:, :, n) = throt(squeeze ( localNetRFFA(n, 1) ), squeeze( localNetRFPhase(n, 1) ) );

                %increment the phase for the following excitation...
                localNetRFPhase(n, 1) = localNetRFPhase(n, 1) + phaseInc;



                MAllExcitations(n, :, excInc) = Atr * squeeze(estExcMat(:, :, n)) * squeeze(MPrev2(n,:))' + abs(squeeze(Meq(n,1))) * Btr;
                MAllExcitationsEachOffRes(n, :, excInc, OffResIter) = MAllExcitations(n, :, excInc);

                %MAllExcitations(all locations, :, some time point) is already
                %decayed. 
            end

            %disp(excInc)
        end

    disp(OffResIter)
end

%the below commented out was just for looking at trends.  
%{
%MAllExcitations has dimensions length(sliceDir) x 3 (mxmymyz) x NumOfExc.

MAllExcitationsxy = squeeze(MAllExcitations(:, 1, :)) + 1i*  squeeze(MAllExcitations(:, 2, :));

%let's sum across all excitations at each time point.  
MAllExcitationsxySUM = zeros(NumOfExc,1);
for n = 1 : NumOfExc
    for m = 1 : size(sliceDir,1)
        MAllExcitationsxySUM(n, 1) = MAllExcitationsxySUM(n, 1) + MAllExcitationsxy(m, n );
    end
end
%}

%-------------------------- -------------------------- -------------------------- -------------------------- -------------------------- -------------------------- 
%% Sum up the complex signal for each (excInc, OffResIter) pair. 

MCompSumExcOffRes = zeros(NumOfExc, NumOfOffRes); 

for excInc = 1 : NumOfExc
    for OffResIter = 1 : NumOfOffRes
        for n = 1 : length(sliceDir)
            MCompSumExcOffRes(excInc, OffResIter) = MCompSumExcOffRes(excInc, OffResIter) + MAllExcitationsEachOffRes(n, 1, excInc, OffResIter) + 1i * MAllExcitationsEachOffRes(n, 2, excInc, OffResIter);
        end
    end
end
%% plot the excitation profile for a particular off resonance, for a particular excitation.
indexOffRes = round( length(PhiPerTR)/2);
%indexOffRes = round( round(length(PhiPerTR)/2) +round(length(PhiPerTR)/4) );

indexExc = NumOfExc;
%indexExc = 1;
figure, 
plot(sliceDir, abs(MAllExcitationsEachOffRes(:, 1, indexExc, indexOffRes  ) + 1i* MAllExcitationsEachOffRes(:, 2, indexExc, indexOffRes  ))  )
xlabel('Length Along Slice Select Direction (mm) ')
ylabel('Transverse Magnetization ')
title(strcat('Excitation Profile After:  ', num2str(indexExc),' Repetitions' ))

%hold on
%for indexOffRes = 1 : length(PhiPerTR)
%plot(sliceDir, abs(MAllExcitationsEachOffRes(:, 1, indexExc, indexOffRes  ) + 1i* MAllExcitationsEachOffRes(:, 2, indexExc, indexOffRes  ))  )
%end
%hold off
%%
figure, imagesc(PhiPerTR, 1 : NumOfExc , abs(MCompSumExcOffRes)), colormap('gray')
ylabel('number of excitations'), xlabel('isochromat off resonant phases')
%% look at the bulk signal distribution over each repetition for a SPECIFIED OFF-res.
figure, plot(abs(MCompSumExcOffRes(1:200, 38))), title('zero')
figure, plot(abs(MCompSumExcOffRes(1:200, 57))), title('zero')
%% look at the bulk signal distribution over steady state over all dephasing.
figure, plot(PhiPerTR,abs(MCompSumExcOffRes(NumOfExc, : ) ) )
%% look at the bulk phase distribution over steady state over all dephasing.
figure, plot(PhiPerTR,angle(MCompSumExcOffRes(NumOfExc, : ) ) )
%% look at the bulk signal's approach to steady state for a particular isochromat
figure, plot(abs(MCompSumExcOffRes(:, round(length(PhiPerTR)/4)  ) ) )
%% look at the signal distribution in steady state at the center isochromat of the slice.
figure, plot(PhiPerTR, abs(squeeze(MAllExcitationsEachOffRes(length(round(sliceDir/2)), 1, NumOfExc, : ) ) + 1i*squeeze(MAllExcitationsEachOffRes(length(round(sliceDir/2)), 2, NumOfExc, : ) )  ))
%% look at the signal variation over each repition for the center isochromat for no off res. 
figure, plot( abs(squeeze(MAllExcitationsEachOffRes(length(round(sliceDir/2)), 1, :, round(NumOfOffRes/2) ) ) + 1i*squeeze(MAllExcitationsEachOffRes(length(round(sliceDir/2)), 2, :, round(NumOfOffRes/2) ) )  ))
%% Phase Encode the Signal.

indexOffRes = round( round(length(PhiPerTR)/2) +round(length(PhiPerTR)/4) );
indexExc = NumOfExc;

M = zeros(length(sliceDir), 3);
for n = 1 : length(sliceDir)
   M(n, :) =  MAllExcitationsEachOffRes(n, :, indexExc,indexOffRes);
end

%*************************
NActualPhaseEncodes = 12;
%*************************
NPixel = NActualPhaseEncodes *  Ns; %for NActualPhaseEncodes phase encodes. 
SamplingFOV = NPixel / IntElementsPerMM;

%old phaseEncoding function
%FourierDataNoEffPE1 = func_phaseEncode_B (SamplingFOV, NPixel, M, sliceDir);
%fftFourierDataNoEffPE1 = ifftnc( FourierDataNoEffPE1 );


dfThetaArr = zeros(length(sliceDir), 1) * PhiPerTR(indexOffRes );

FourierDataNoEffPE = func_phaseEncode_B (SamplingFOV, NPixel, M, sliceDir, dfThetaArr );
fftFourierDataNoEffPE = ifftnc( FourierDataNoEffPE );
%% plot the results
figure, 
%plot(linspace(-SamplingFOV/2, SamplingFOV/2,length(fftFourierDataNoEffPE) ), abs(fftFourierDataNoEffPE),  'LineWidth',2.0)
cplot2(linspace(-SamplingFOV/2, SamplingFOV/2,length(fftFourierDataNoEffPE) ), (fftFourierDataNoEffPE1))
title(strcat(num2str(NActualPhaseEncodes),' Phase Encoding Steps  '))
xlabel('Encoded Length Along the Slice Select Direction (mm) ', 'FontSize',12,'FontWeight','bold')

%% Write some code to have only the "in-slice" contribution.  
%for each off-resonance accumulated/TR, you're going to want to look at the
%initial excitation profile. 
%the line below commented out shows the profile (z mag) for the
%PhiPerTR(34) off res/TR.  for in slice weight it by 1-Mz value.
%figure, plot(abs( MAllExcitationsEachOffRes(:, 3,1, 34   ) ) )

wSliceProfileOffRes = zeros(length(sliceDir), NumOfOffRes);
wSliceMAllExcitationsEachOffRes = zeros(size(MAllExcitationsEachOffRes,1), size(MAllExcitationsEachOffRes,3), size(MAllExcitationsEachOffRes,4));

for n = 1 : NumOfOffRes
    wSliceProfileOffRes(:,n) = 1 - abs( MAllExcitationsEachOffRes(:, 3,1, n   ) );
end

