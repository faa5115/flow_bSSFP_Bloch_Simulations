%We will be doing this to test our phase-encoding function.
%we will be generating the data from scratch for a particular dS and
%offRes.

%FOR SOME REASON THIS CODE PROPERLY WORKS WITH AN EVEN NUMBER OF ENCODING
%STEPS/ EVEN VPS.  

clear
%you want to control the duration, the thickness, excitation profile...

sliceOffCenter = 0;
%RFPhase0 = pi/2;
RFPhase0 = pi;
%RFPhase0 = pi; 


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
T2Decay= 4;
TR = 4e-3;   %s
TE = TR/2; %s

T2 = 150e-3;  %s
T1 = 1000e-3;  %s

%T1 = 500e-3; %s
%T2 = 500e-3; %s

T2Fac = 3; %unitless.  just saying after how many T2s will we  begin ignoring signal. 

NumOfOffRes = 75;
PhiPerTR = linspace(-2*pi, 2*pi, NumOfOffRes);
OffRes   = PhiPerTR/(2 * pi * TR); %Hz. 


%dS is going to be a parameter that defines how quickly the spins flow
%between TRs.  dSMin is always going to be one element shift / TR
%therefore it'll always be 1/Ns.  
%otherwise ... dS will be zero (stationary spins).  
%dS * Ns is the number of elements an isochromat shifted between TRs.  
dSMin = 1/Ns; %smallest nonzero spin exchange percentage.  always going to be 1/Ns.

%NPart = 6;
NPart = 4;
NEffecgtiveFOVElements = Ns * NPart;

TotCyclicTimeAcqOG = 1; %seconds.



%-----------NON-CYCLIC DATA------------------------------------------------
%{
intendedVPS = 25;
NActualPhaseEncodes = 25;

numOfTP_og = 100;
subTractFromTP_og = mod(numOfTP_og, intendedVPS);
numOfTP = numOfTP_og - subTractFromTP_og;

%dS = 0.3; %fraction. 
%dS = 0.2; %example

%{
dS = 0;
OffResIter = round(length(PhiPerTR)/2);
%}


dS = 1.0;
OffResIter = round(length(PhiPerTR)/2) + round(length(PhiPerTR)/4);


NShift = round(dS * Ns);


dSsim_pulse = dS * ones(1, numOfTP);
dSsim_total = dS * ones(1, numOfTP);
NShift = round(dS * Ns);

NumOfPrepExc = length(dSsim_pulse) * 4;
PeakVelExc = round(length(dSsim_pulse)/2); %doesn't matter what value it is because the flow is constant...
%}
%--------------------------------------------------------------------------

%----CYCLIC DATA-----------------------------------------------------------


%indexOffRes = round(length(PhiPerTR)/2) + round(length(PhiPerTR)/4); %pi


intendedVPS = 1;
NActualPhaseEncodes = 1;

numOfTP_og = 250;
subTractFromTP_og = mod(numOfTP_og, intendedVPS);
numOfTP = numOfTP_og - subTractFromTP_og;

NumOfExcPerBeat = numOfTP;
t_pulse = 0 : TR : (NumOfExcPerBeat - 1) * TR; %ms
%vMax = 1600; %mm/s
vMax = 1500; %mm/s
dSMax = TR * vMax/sliceThickness; 
dS = dSMax;
%PeakVelExc = round(NumOfExcPerBeat / 3); 
PeakVelExc = round( numOfTP * 0.18 );
stdDevExc = 10;  
dSsim_pulse = dSMax * exp(- ((t_pulse - PeakVelExc * TR).^2)/(2*((stdDevExc*TR)^2)));
NumOfHeartBeats = 1;

dSsim_total = zeros(1, NumOfHeartBeats * length(dSsim_pulse) );

for n = 1 : NumOfHeartBeats
    dSsim_total(1, (n - 1) * NumOfExcPerBeat + 1 : n * NumOfExcPerBeat  ) = dSsim_pulse;
end

%OffResIter = round(length(PhiPerTR)/2);
OffResIter = round(length(PhiPerTR)/2) + round(length(PhiPerTR)/4); %pi
indexExc = PeakVelExc; %peak
numOfCompletedCycles = length(dSsim_total) / length(dSsim_pulse); 

%NumOfExc = 1000;
%NumOfPrepExc = NumOfExc + 1;
NumOfPrepExc = length(dSsim_pulse) * 4;

%--------------------------------------------------------------------------



NOs =  round(T2Fac * dS * T2 * Ns /TR);
%sliceFOV = (Ns + NOs)*subSliceThickness; %mm
sliceFOV = (NEffecgtiveFOVElements + 2* NOs)*subSliceThickness; %mm


sliceDir = (linspace(-sliceFOV/2, sliceFOV/2, IntElementsPerMM * sliceFOV))'; %mm

%intElementsOnSlice = 100; %number of elements along one axis of a slice.  


RFBW = TBW/RFDur; %Hz %this is the RF's bandwidth. ...
Gslice = RFBW/(gamma*sliceThickness/(2*pi)); %T/mm

B1e = getSincEnvelope_fa(gamma, RFDur, NumOfTPRF, TBW, flipAngle);
%B1efirst = getSincEnvelope_fa(gamma, RFDur, NumOfTPRF, TBW, flipAngle/2);


M0     = zeros(size(sliceDir,1), 3);
%MbSSFP = zeros(size(sliceDir,1), 3, NumOfExc );

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
offresonanceones = ones(size(sliceDir));


%---------COIL SENSITIVITY----------------------------
sdcs  = 20;
g = exp( -0.5 * ( sliceDir / sdcs ).^2 );
%-----------------------------------------------------


%% Prep pulses. 

%just an array of the rad/TR throughout each spatial location. 
dfThetaArr = squeeze(PhiPerTR(OffResIter)) ...
    * ones(length(sliceDir), 1);

%now prepare the decay matrices...
df = squeeze(OffRes(1,OffResIter));
[Atr, Btr] = freeprecess(TR,T1,T2,df);

MSliceProfile0 = func_sliceSelection (B1e, sliceDir, RFBW, sliceThickness, RFDur, M0, gamma, sliceOffCenter, RFPhase0, Meq, T1, T2, df *  offresonanceones );

MSliceProfile=func_sliceSelRef (RFDur, RFBW, sliceThickness,  B1e, squeeze(MSliceProfile0(:, :, size(MSliceProfile0,3))), sliceDir, Meq, T1, T2, df *  offresonanceones);

% calculate the flip angle and the phase of each isochromat.   --------------------------------------------

localNetRFPhase = zeros(length(sliceDir), 1);
localNetRFFA    = zeros(length(sliceDir), 1);

for n = 1: length(sliceDir)
    localNetRFPhase(n, 1) = angle(  squeeze(MSliceProfile(n,1,size(MSliceProfile,3))) + ...
        1i *squeeze(MSliceProfile(n,2,size(MSliceProfile,3))));
    
    localNetRFFA   (n, 1) = angle( 1i*abs(squeeze(MSliceProfile(n,1,size(MSliceProfile,3))) + ...
        1i* squeeze(MSliceProfile(n,2,size(MSliceProfile,3)))) + ...
        (squeeze(MSliceProfile(n,3,size(MSliceProfile,3)))  )     );
end
    
        % calculate the proper rotation matrix for each isochromat.  
    estExcMat = zeros(3, 3, length(sliceDir)); %estimated excitation matrix. 
    for n = 1 : length(sliceDir)
        estExcMat(:, :, n) = throt(squeeze ( localNetRFFA(n, 1) ), squeeze( localNetRFPhase(n, 1) ) );
    end
        localNetRFPhase0 = localNetRFPhase;
    localNetRFFA0 = localNetRFFA;
    
    estExcMat0 = estExcMat;


MAllExcitations = zeros(size(M0,1), size(M0,2), NumOfPrepExc);



localNetRFPhase = localNetRFPhase0;
localNetRFFA = localNetRFFA0;
estExcMat = estExcMat0;
MNextCycle = M0;

excInc = 1;


MPrev = M0;




for n = 1 : length(sliceDir)

    %prepare the rotation excitation matrics...

    estExcMat(:, :, n) = throt(squeeze ( localNetRFFA(n, 1) ), squeeze( localNetRFPhase(n, 1) ) );

    %increment the phase for the following excitation...
    localNetRFPhase(n, 1) = localNetRFPhase(n, 1) + phaseInc;



   %  MAllExcitations(n, :, excInc) = Atr * squeeze(estExcMat(:, :, n)) * squeeze(MPrev2(n,:))' + abs(squeeze(Meq(n,1))) * Btr;

end

prevExcitations = 0; %number of excitations we've already done from the 
prevExcitationsPrep = 0;
%func_FLOW_PartEncode_B_II_sp function.
spKSpace = NumOfPrepExc / 2;
%spKSpace * deltaK is our k-space coordiante that we start from when 
%func_FLOW_PartEncode_B_II_sp is called. 
sp = NumOfPrepExc; %irrelevant...idk why I still have it. 

firstTarget = 1;
%firstTarget indicates which point within the cycle we start. 

%dummy
NPixel = NumOfPrepExc;
%dummy
SamplingFOV = NPixel * sliceThickness; %mm
%SamplingFOV = NPixel * subSliceThickness; %mm

%dummy
VPS = NumOfPrepExc;
%dummy
linesInSegment = VPS;

MStartPrep = MPrev;

%pre-kspace ---------------------------------------------------------------

%RFPhaseArrPrep is the central spin-ensemble's phase over time. 
        [FlowFourierDataPrep, MTREndPrep, MTRBegPrep, nextFirstTarget,RFPhaseArrPrep , prevExcitationsPrep] = ...
                                func_FLOW_PartEncode_B_II_sp (prevExcitationsPrep, spKSpace, sp,  ...
                                                        firstTarget, localNetRFFA,...
                                                 localNetRFPhase, phaseInc, ...
                                                 T1, T2, TR, TE, dSsim_pulse, ...
                                                  Ns, SamplingFOV, linesInSegment, ... 
                                                 MStartPrep, sliceDir, dfThetaArr, Meq);
%pre-kspace  end ----------------------------------------------------------


%save(strcat('CyclicMTREndPrep_',num2str(PhiPerTR(OffResIter)),'radPerTR',num2str(VPS),'.mat'),'MTREndPrep','-v7.3');

nextFirstTargetPrep = nextFirstTarget;
%% Debug the pre-kspace steps...
%MTREnd's size:  # elements along slice x #pixels (TRs) x Mxyz.  
MTREndBegCompSum = zeros( size(MTREndPrep, 2), 1);
for n = 1 : NumOfPrepExc
    for m = 1 : length(sliceDir)
        MTREndBegCompSum(n, 1) = MTREndBegCompSum(n, 1) + ...
            squeeze(MTREndPrep(m, n, 1 )) + 1i*squeeze(MTREndPrep(m, n, 2 ));
    end
end
            %%  Plot the abs(complex sum) so far...
            %{
            %figure, cplot(MTREndBegCompSum)
            figure, plot(abs(squeeze(MTREndBegCompSum(:, 1))), 'LineWidth', 5.0)
            figure, plot((1 : length(dSsim_total)) , dSsim_total,'LineWidth',5.0 )
            title('Simulated Pulse')
            ylabel('\Delta s' )
            xlabel('Time (#TRs)')
            %}
            
            %{ 
            figure, imshow( )
            %}
            %% Plot the excitation's signal and phase profile.
            
            %take the last cycle's contribution ("equilibrium").
            MTREndPrepLastCycle = squeeze(MTREndPrep(:,NumOfPrepExc - length(dSsim_pulse) + 1 : NumOfPrepExc , : ));
            
            %indexExc = NumOfPrepExc;
            indexExc = PeakVelExc;
            %indexExc = PeakVelExc + 80;
            %indexExc = 2;
            %indexExc = 250;
           %indexExc = 20;
           % indexExc = 220;
           %indexExc = 110;
            %{
            indexExc = 43;
            indexExc = 55;
            indexExc = 54;
            indexExc = 2;
            %}
            

            
            figure, 
            hold on

            magImage = abs(squeeze(MTREndPrepLastCycle(:, indexExc, 1)) + 1i * squeeze(MTREndPrepLastCycle(:, indexExc, 2)));
            phsImage = angle(squeeze(MTREndPrepLastCycle(:, indexExc, 1)) + 1i * squeeze(MTREndPrepLastCycle(:, indexExc, 2)));
            % FADILALI -------- KEEP THIS IFFTSHIFT-IFFTNC ARRANGMENT.  IT WILL BE THE
            % OPPOSITE FLOW.  
            plot(sliceDir, ... 
                magImage / max(magImage(:)), 'LineWidth',5.0)


            xline(3, 'LineWidth',2.0), xline(-3,'LineWidth',5.0)
            xlim([-10 10])


            plot(sliceDir, ...
                phsImage / max(phsImage(:)), 'LineWidth',5.0)

            title(strcat('Analog Signal Time Point:  ', 32, num2str(indexExc)))

            %% Plot the RF's phase. 
            %{
            figure, plot(RFPhaseArrPrep)
            %}
%% Now for the actual k-space acquisition.
%you're going to want to change a lot of the previous variables' values. 

%--------------------------------------------------------------------------
%                           FADIL ALI
%Fix this section when you return in the morning.
%you must loop over segments and frames...
%and appropriately set the prevFirstTargetARR...
%--------------------------------------------------------------------------

%set the proper localNetRFPhase:
%   in the func_FLOW_PartEncode_B_II_sp called above, RFPhaseArrPrep is an
%   outputted array that stores the center spin-ensemble's RF-Phase for 
%   each excitation.  
%   therefore for the start of the k-space sampling, we must increment 
%   all of localNetRFPhase0's elements by (length(RFPhaseArrPrep) + 1) *
%   phaseInc.
for n = 1 : length(sliceDir)
    localNetRFPhase(n, 1) = squeeze(localNetRFPhase0(n, 1)) + (length(RFPhaseArrPrep) + 1) * phaseInc;
end

TotCycleTime = TR * length(dSsim_pulse);

%{
NActualPhaseEncodes  = 100; 
%NActualPhaseEncodes  = 100;
NPixel = NActualPhaseEncodes;
VPS = 20; 
%VPS = 100;
%}

%NActualPhaseEncodes = 1;
%NActualPhaseEncodes = 6;
%NActualPhaseEncodes  = 12; 
%NActualPhaseEncodes  = 20;
%NActualPhaseEncodes = 30;
%NActualPhaseEncodes = 50;
%NActualPhaseEncodes  = 100;
%NActualPhaseEncodes  = 240;
NPixel = NActualPhaseEncodes;

%set  VPS
%VPS = 240; 
%VPS = 100;
%VPS = 12;
%VPS = 2;
%VPS = 1;
VPS = intendedVPS;
%SamplingFOV = NPixel * subSliceThickness; %mm
SamplingFOV = NPixel * sliceThickness; %mm



if (VPS > NPixel) 
     error(strcat('VPS cannot be higher than NPixel') );
end
tempRes = VPS * TR; 
if (tempRes > TotCycleTime) 
     error(strcat('tempRes', 32, num2str(tempRes), 32,', cannot be higher than TotCycleTime,',32, ...
         num2str(TotCycleTime),'.') );
     
end

%FADILALI you have an Acquisition Window bug...go to the flowPE function
%below and fix it...
AcqWindow = TotCycleTime; 
%AcqWindow = 350 * TR; 
NumOfFrames = floor(AcqWindow / tempRes);  

%our outer loop is going to be the number of segments. 
%then we nest a loop over frames.
%finally we loop over views within that frame's segment. 

if mod(NPixel, VPS) == 0
    numOfSegments = NPixel / VPS; 
    linesInFinalSegment = VPS; 
   
else
    numOfSegments = floor( NPixel / VPS ) + 1; 
    linesInFinalSegment = NPixel - VPS * floor ( NPixel / VPS ); 
    %or ... you could just have:  mod(NPixel, VPS);
end

rawData = zeros(NPixel, NumOfFrames);
FormattedSpeedArr = zeros(NPixel, NumOfFrames);
%               slice's length    , VPS,  Num of Frames,  numOfSegments, mxxyz
MTREndAll = zeros(length(sliceDir), VPS ,NumOfFrames , numOfSegments, 3 );

MTRBegAll = zeros(length(sliceDir), VPS ,NumOfFrames , numOfSegments, 3 );

FlowFourierDataAll = zeros(VPS, NumOfFrames, numOfSegments );

%just for debugging purposes...store the dS values. 
StoreSpeedArr = zeros(VPS, NumOfFrames, numOfSegments);

incPECountDebug = 1; 
%asOrDes = 'd';
sp = NPixel;
spKSpace = NPixel/2;

firstTarget = nextFirstTarget;

MStart = squeeze( MTREndPrep(:, NumOfPrepExc, : ) );

prevFirstTargetArr = zeros(numOfSegments,1 );
prevFirstTargetArr(1, 1) = firstTarget - 1;

incSegFrames = 1; 
RFPhaseArrTot = zeros (NPixel * NumOfFrames , 1);
disp(strcat(num2str(numOfSegments),' .  Total Segments'))

for incSegments = 1  : numOfSegments
    
    disp( strcat('Segment ',32, num2str(incSegments) ) )
    
    if incSegments == numOfSegments
        linesInSegment = linesInFinalSegment;
    else
        linesInSegment = VPS; 
    end
    
    for incFrames = 1 : NumOfFrames
        disp( strcat('Frame ',32, num2str(incFrames) ) )
        %{
        for incLines = 1 : linesInSegment
            %should end up being 1 + NPixel * NumOfFrames.
             incPECountDebug = incPECountDebug + 1;
        end
        %}
        
        %this is what we originally had...and it worked.. did not violate
        %ss. but we need to generalize this.  because we may not have bSSFP
        %in future simulations. 
        
        %******************************************************************
        %{
        RFPhaseFactor = NumOfExc + firstTarget + prevFirstTargetArr(incSegments, 1 ) ; 
        firstTarget = prevFirstTargetArr(incSegments, 1 )  + 1;
        for n = 1 : length(sliceDir)
            localNetRFPhase(n, 1) = localNetRFPhase0(n, 1) + RFPhaseFactor * phaseInc;
        end
        %}
  
        %FADILALI MAJOR CHANGE...WE ARE COMMENTING THE LINE BELOW OUT!
        %firstTarget = prevFirstTargetArr(incSegments, 1 )  + 1;
        
        %we are going to do this for each segment/frame's iteratoin to safely guarantee the
        %RF's proper phase.
        %this ramps up to the RF phase we want when going into the
        %func_FLOW_PartEncode_B_II_sp function.  if you want to see the
        %total RF Phase's plot, you'll have to plot the final iteration's 
        %RFPhaseArrTot from the loop below along plus the RFPhaseArr output
        %from the function below. 
        %{
        RFPhaseFactor = 0;
        for incPrevExc = 1 : NumOfExc + prevExcitations + 1
            RFPhaseFactor = RFPhaseFactor + 1;
            for n = 1 : length(sliceDir)
                localNetRFPhase(n, 1) = localNetRFPhase0(n, 1) + (RFPhaseFactor ) * phaseInc;
                if n == round(length(sliceDir) / 2) 
                    RFPhaseArrTot(incPrevExc, 1) = localNetRFPhase(n, 1);
                    disp(strcat('RFPhaseArrTot ', 32, num2str( squeeze( RFPhaseArrTot(incPrevExc, 1) ) ) ))
                end
            end
        end
        %}
        
        %debug
        %{
        figure, plot(RFPhaseArrTot)
        figure, cplot(localNetRFPhase0)
        figure, cplot(localNetRFPhase)
        %}
        %******************************************************************
        
        
            
        %MSpatial is the total spatial information for each phase encoding
        %line.  
        [FlowFourierData, MTREnd, MTRBeg, nextFirstTarget,RFPhaseArr , prevExcitations, SpeedArr] = ...
                                func_FLOW_PartEncode_B_II_sp_CoilSense (prevExcitations, spKSpace, sp,  ...
                                                        firstTarget, localNetRFFA,...
                                                 localNetRFPhase, phaseInc, ...
                                                 T1, T2, TR, TE, dSsim_pulse, ...
                                                  Ns, SamplingFOV, linesInSegment, ... 
                                                 MStart, sliceDir, dfThetaArr, Meq, g);
                                             
        StoreSpeedArr(:, incFrames, incSegments) = SpeedArr;
                                             
        %rawData( sp - linesInSegment + 1 : sp, incFrames) =  FlowFourierData;
        FlowFourierDataAll(:, incFrames, incSegments ) = FlowFourierData;
        
        %commenting the below two lines out...does not match the dimensions
        %we reserved for them above...
        %{
        MTREndAll(:, sp - linesInSegment + 1 : sp, incSegFrames, :) = MTREnd; 
        MTRBegAll(:, sp - linesInSegment + 1 : sp,  incFrames, incSegments, :) = MTRBeg; 
        %}
        
        %{
            figure, 
            hold on
            plot(1 : length(FlowFourierData), abs( ( (FlowFourierData))), 'LineWidth',2.0)
            plot(1 : length(FlowFourierData), angle( ( (FlowFourierData))),'LineWidth',2.0 )
            legend('Magnitude', 'Phase')
            title(strcat('Frame',32, num2str(incFrames),32,'Segment',32,num2str(incSegments) ) )
            hold off
            %}
        
        %I think this is better.
        %       slice's length, VPS,  Num of Frames,  numOfSegments, mxxyz
        MTREndAll(:, :,  incFrames, incSegments, :) = MTREnd; 
        MTRBegAll(:, :,  incFrames, incSegments, :) = MTRBeg;        
        
        %MTREndAll(:, :, incSegFrames, :) = MTREnd;
        %MTRBegAll(:, :, incFrames, incSegments, :) = MTRBeg;
        MStart = squeeze(MTREnd(:, size(MTREnd,2), :));
        %this is how we keep cycling through the dS-values array.
        %idk why you have "firstTarget - 1" ...
        %prevFirstTargetArr(incSegments, 1 ) = firstTarget - 1 ;
        prevFirstTargetArr(incSegments, 1 ) = firstTarget;
        firstTarget = nextFirstTarget;
        incSegFrames = incSegFrames + 1;
        
        RFPhaseArrTot( (prevExcitations - VPS  + 1) : prevExcitations  ...
            , 1) = RFPhaseArr ;
        
        for n = 1 : length(sliceDir)
            localNetRFPhase(n, 1) = squeeze(localNetRFPhase(n, 1)) + linesInSegment * phaseInc;
        end
        

        
        %disp(strcat('incFrames', num2str(incFrames)))
    end
    %once you finish all frames, you increment to the next starting 
    %PE line. 
    spKSpace = spKSpace - linesInSegment;
    sp = sp  - linesInSegment ;
    %disp(strcat('sp', num2str(sp)))
    
end
%% Let's try to rotate the FlowFourierDataAll to the ADC's phase. 
incTR = 1;
%we will update FlowFourierDataAll to rotate to the ADC's phase.  will
%store the previous values as FlowFourierDataAllPreADC.
FlowFourierDataAllPreADC = FlowFourierDataAll;
for incSegments = 1 : numOfSegments 
    for incFrames = 1 : NumOfFrames
        for incPixels = 1 : VPS
            RotToADCPhase = zrot( (NumOfPrepExc + incTR - 1) * phaseInc );
            vecToRot = [real(squeeze(FlowFourierDataAllPreADC(incPixels, incFrames, incSegments))); ...
                        imag(squeeze(FlowFourierDataAllPreADC(incPixels, incFrames, incSegments))); ...
                        0];
            newVec = RotToADCPhase * vecToRot;
            FlowFourierDataAll(incPixels, incFrames, incSegments) = squeeze(newVec(1,1)) + 1i* squeeze(newVec(2,1));
          
            incTR = incTR + 1;
        end
        
        %{
        figure, 
        hold on
        plot(1 : length(FlowFourierData), abs( squeeze(FlowFourierDataAll(:, incFrames, incSegments ) ) ), 'LineWidth',2.0)
        plot(1 : length(FlowFourierData), angle( squeeze(FlowFourierDataAll(:, incFrames, incSegments ) ) ),'LineWidth',2.0 )
        legend('Magnitude', 'Phase')
        title(strcat('Frame',32, num2str(incFrames),32,'Segment',32,num2str(incSegments) ) )
        hold off
        %}
    end
end

%% Properly Arrange FlowFourierDataAll into K-Space Format. and StoreSpeedArr

%{
%just to make sure FlowFourierDataAll has the values that I expect it to.
for incSegments = 1 : numOfSegments
   for  incFrames = 1 : NumOfFrames
                   figure, 
            hold on
            plot(1 : length(FlowFourierData), abs( ( (FlowFourierDataAll(: ,incFrames, incSegments )))), 'LineWidth',2.0)
            plot(1 : length(FlowFourierData), angle( ( (FlowFourierDataAll(: ,incFrames, incSegments )))),'LineWidth',2.0 )
            legend('Magnitude', 'Phase')
            title(strcat('Frame',32, num2str(incFrames),32,'Segment',32,num2str(incSegments) ) )
            hold off
   end
end
%}

for incFrames = 1 : NumOfFrames
    sp = NPixel;
    for incSegments = 1 : numOfSegments
        rawData( sp - VPS + 1 : sp, incFrames) = squeeze(FlowFourierDataAll(:, incFrames, incSegments));
        FormattedSpeedArr(sp - VPS + 1 : sp, incFrames) = squeeze(StoreSpeedArr(:, incFrames, incSegments));
        sp = sp - VPS;
    end
end

%{
            save(strcat('CyclicFormattedSpeedArr_',num2str(NPixel),'PhaseEncodes_',num2str(VPS),'VPS_StartsAt', ...
                num2str(nextFirstTargetPrep),'_',num2str(NumOfFrames),'Frames_', ...
                num2str(numOfSegments),'Segments_','.mat'),'FormattedSpeedArr','-v7.3');
%}

% figure, plot(abs(squeeze(rawData(:, 1))))
            %% Plot:  debug the FormattedSpeedArr.
            indexFrame = round(PeakVelExc / intendedVPS);
         % indexFrame = 29;
            %indexFrame = 114;
            figure, 
            plot(flip(squeeze(FormattedSpeedArr(:, indexFrame ))))
            title(strcat('speed array for each frame (reversed) time point ',32,num2str(indexFrame)))
            %% Plot:  debug the FlowFourierData within the loop

            %{
            plot(linspace(-SamplingFOV/2, SamplingFOV/2,length(FlowFourierData) ) ...
                    , abs((fftnc(FlowFourierData))), 'LineWidth',2.0)
            %}
            %    , abs( ( (FlowFourierData))), 'LineWidth',2.0)

            %{
            plot(linspace(-SamplingFOV/2, SamplingFOV/2,length(FlowFourierData) ) ...
                ,angle( ( (FlowFourierData))),'LineWidth',2.0 )
            %}
            %   ,angle((fftnc(FlowFourierData))),'LineWidth',2.0 )


            %FADILALI we incorporated this in the section above...
            %           it's a loop's output.
            %{
            figure, 
            hold on
            plot(1 : length(FlowFourierData), abs( ( (FlowFourierData))), 'LineWidth',2.0)
            plot(1 : length(FlowFourierData), angle( ( (FlowFourierData))),'LineWidth',2.0 )
            legend('Magnitude', 'Phase')
            title(strcat('Frame',32, num2str(NumOfFrames),32,'Segment',32,num2str(numOfSegments) ) )
            hold off
            %}

            %%
            %1 vps test
            indexFrame = round(PeakVelExc/intendedVPS) ;
          %indexFrame = 110;
           %indexFrame = 114;

            %10 vps test
            %indexFrame = 16;

            %OKAY HERE IS YOUR PROBLEM...EACH LINE IS HAVING THE SAME MAGNITUDE...FIX
            %THAT!
            figure, 
            hold on

            % FADILALI -------- KEEP THIS IFFTSHIFT-IFFTNC ARRANGMENT.  IT WILL BE THE
            % OPPOSITE FLOW. 
            magImage = abs( flip((ifftnc( squeeze( rawData(:, indexFrame) ) ))) );
            phsImage = angle( flip((ifftnc( squeeze( rawData(:, indexFrame) ) ))) );
            plot(linspace(-round(NPixel/2), round(NPixel/2),size(rawData,1) ) ...
            ,magImage / max(magImage(:)), 'LineWidth',2.0)


           


            plot(linspace(-round(NPixel/2), round(NPixel/2),size(rawData,1) ) ...
                ,phsImage / max(phsImage(:)), 'LineWidth',2.0)

            title(strcat('Total Sum In Each Partition:',32,32,'Frame#', 32, num2str(indexFrame)))
            hold off
            legend('Magnitude', 'Phase')
%% Check the steady-state. 
MTREndAllCompSum = zeros(NPixel * NumOfFrames, 1);

incTR = 1; 

%the commented out section below is irrelevant now because you changed 
%MTREndAll's structure.
%               slice's length    , VPS,  Num of Frames,  numOfSegments, mxxyz
%{
for frameSegInc = 1 : size(MTREndAll,3)
    for excVPS = 1 : size(MTREndAll,2)
        for sliceIter = 1 : length(sliceDir)
            MTREndAllCompSum( (frameSegInc-1) * size(MTREndAll,2) + excVPS,1) =  ...
                squeeze(MTREndAllCompSum((frameSegInc-1) + excVPS, 1)) ...
                + squeeze(MTREndAll(sliceIter, excVPS, frameSegInc, 1)) + ...
                1i * squeeze (MTREndAll(sliceIter, excVPS, frameSegInc, 2 ));
        end
    end
end
%}

%{
for incFrame = 1 : size(MTREndAll,3)
    for incSeg = 1 : size(MTREndAll, 4)
        for excVPS = 1 : size(MTREndAll,2)
            for sliceIter = 1 : length(sliceDir)
                MTREndAllCompSum( (incFrame * incSeg-1) * size(MTREndAll,2) + excVPS,1) =  ...
                    squeeze(MTREndAllCompSum((incFrame * incSeg-1) * size(MTREndAll,2) + excVPS, 1)) ...
                    + squeeze(MTREndAll(sliceIter, excVPS, incFrame, incSeg, 1)) + ...
                    1i * squeeze (MTREndAll(sliceIter, excVPS, incFrame, incSeg, 2 ));
            end
        end
    end  
end
%}

%MTREndAll's structure.
%               slice's length    , VPS,  Num of Frames,  numOfSegments, mxxyz

incTR = 1;
for incSeg = 1 : size(MTREndAll, 4)
    for incFrame = 1 : size(MTREndAll,3)
        for excVPS = 1 : size(MTREndAll,2)
            for sliceIter = 1 : length(sliceDir)
                MTREndAllCompSum( incTR,1) = squeeze(MTREndAllCompSum(incTR, 1)) + ...
                     squeeze(MTREndAll(sliceIter, excVPS, incFrame, incSeg, 1)) + ...
                    1i * squeeze (MTREndAll(sliceIter, excVPS, incFrame, incSeg, 2 ));
            end
            incTR = incTR + 1;
        end
    end  
end
%by the end incTR should = NPixel * NumOfFrames + 1. 
%{
save(strcat('CyclicMTREndAll_',num2str(NPixel),'PhaseEncodes_',num2str(VPS),'VPS_StartsAt', ...
                num2str(nextFirstTargetPrep),'_',num2str(NumOfFrames),'Frames_', ...
                num2str(numOfSegments),'Segments_','.mat'),'MTREndAll','-v7.3');
%}

%% Plot the total abs(Complex sum) 

figure, plot(  abs(squeeze(MTREndBegCompSum(:, 1))), 'LineWidth', 5.0)
figure, plot( [ abs(squeeze(MTREndBegCompSum(:, 1)));  abs(squeeze(MTREndAllCompSum(:, 1))) ], 'LineWidth', 5.0)
figure, plot(abs(squeeze(MTREndAllCompSum(:, 1))), 'LineWidth', 5.0)

%% Plot the RF's phase. 
%{
figure, plot(RFPhaseArrPrep, 'LineWidth', 5.0), title('prep only')
figure, plot([RFPhaseArrPrep; RFPhaseArrTot], 'LineWidth', 5.0), title('prep + kspace')
figure, plot(RFPhaseArrTot, 'LineWidth', 5.0), title('kspace only')
%}
                %% PLOT to check squeeze(MTRBegAll(:, incRepetitions, incFrames, incSegments, : ));
                %{
            %IMPORTANT DEBUG:  used this to check MTRBegAll is proper.
                for incFrames = 1 : NumOfFrames
                    for incSegments = 1 : numOfSegments
                        for incRepetitions = 1 : VPS
                            figure, 
                            hold on

                            % FADILALI -------- KEEP THIS IFFTSHIFT-IFFTNC ARRANGMENT.  IT WILL BE THE
                            % OPPOSITE FLOW.  
                            plot(sliceDir, abs(  squeeze( MTRBegAll(:, incRepetitions, incFrames, incSegments, 1) )  + ...
                                1i*squeeze( MTRBegAll(:, incRepetitions, incFrames, incSegments, 2)  )), 'LineWidth',2.0)





                            plot(sliceDir,angle(  squeeze( MTRBegAll(:, incRepetitions, incFrames, incSegments, 1) )  + ...
                                1i*squeeze( MTRBegAll(:, incRepetitions, incFrames, incSegments, 2)  )), 'LineWidth',2.0)

                            xline(3, 'LineWidth',2.0), xline(-3,'LineWidth',2.0)
                            title(strcat('Frame#', 32, num2str(indexFrame), 32, 'repetition# ', 32, num2str(incRepetitions), 32, 'Segment# ', 32, num2str(incSegments)))
                            hold off
                            legend('Magnitude', 'Phase')
                        end
                    end
                end
                %}
%% Now phase-encode within the phase encoded lines to resolve the subslices. 

NSubSlicePixels = Ns;
SubSliceSamplingFOV = (NSubSlicePixels * NPixel ) * subSliceThickness;

SubFlowFourierDataAll = zeros(NSubSlicePixels,VPS, NumOfFrames, numOfSegments );

%just the off-resonant frequency throughout space. 
OffResArr   = dfThetaArr/(2 * pi * TR);

for incFrames = 1 : NumOfFrames
    spKSpace = NSubSlicePixels * NPixel / 2;
    for incSegments = 1 : numOfSegments
        
        %WARNING: FADIL ALI CHECK TO SEE IF BELOW SHOULD BE NPIXEL OR VPS!!
        for incRepetitions = VPS : -1 : 1
            
            
            %MInput = squeeze(MTRBegAll(:, incRepetitions, incFrames, incSegments, : ));
            %...let's test ... let's focus on one repetition's data.
            MInput = squeeze(MTRBegAll(:, incRepetitions, incFrames, incSegments, : ));
            
            MInputTE = zeros(size(MInput));
            for incSlice = 1 : length(sliceDir)
                [Ate, Bte] = freeprecess(TE ,T1,T2,squeeze(dfThetaArr(n)));
                %you don't have to include the Off-resonance in our phase
                %encoding call because we are already including it in the TE
                %decay. 
                RotToADCPhase = zrot( (incRepetitions - 1) * phaseInc );
                MInputTE(incSlice, :) = Ate * RotToADCPhase * squeeze(MInput(incSlice, :))' + Bte * squeeze( Meq(incSlice, 1) );
            end
            SubFlowFourierData = func_phaseEncode_B_sp_CoilSense(spKSpace, SubSliceSamplingFOV, NSubSlicePixels, ...
                MInputTE, sliceDir, zeros(length(dfThetaArr), 1), g);
            
            SubFlowFourierDataAll(:, incRepetitions, incFrames, incSegments) = SubFlowFourierData;
            
            %IMPORTANT DEBUG:  used this to check to see if you organized
            %the plots properly in the section below
            %{
            figure,
            hold on
            plot(abs(squeeze(SubFlowFourierData(:, 1))))
            plot(angle(squeeze(SubFlowFourierData(:, 1))))
            hold off
            title(strcat('TR ', 32, num2str(incRepetitions), 32, ' Segment ', 32, num2str(incSegments), 32, ' Frame ', 32, num2str(incFrames)))
            %}
            spKSpace = spKSpace - NSubSlicePixels;
        end
    end
    
end
%% SubSpace RawData
SubRawData = zeros(NSubSlicePixels * NPixel, NumOfFrames);

%FADILALI THE LOOPING WORKS FOR SINGLE SHOT...CHECK TO SEE IF IT WORKS FOR
%MULTIPSHOT/SEGMENTS!!!! YOU MAY HAVE TO CHANGE WERE sp = NSubSlicePixels * NPixel; IS.
for incFrames = 1 : NumOfFrames
    sp = NSubSlicePixels * NPixel;
    for incSegments = 1 : numOfSegments
        for incRepetitions = VPS : -1 : 1
            SubRawData( sp - NSubSlicePixels  + 1 : sp, incFrames) = squeeze(SubFlowFourierDataAll(:, incRepetitions,  incFrames, incSegments));
            sp = sp - NSubSlicePixels;
        end
        
    end
end
%% PLOT subrawdata

            indexFrame = round(PeakVelExc/intendedVPS);
           % indexFrame = 20;
          % indexFrame = 2;
           % indexFrame = 43;
           % indexFrame = 1;

            %10 vps test
            %indexFrame = 16;

            %OKAY HERE IS YOUR PROBLEM...EACH LINE IS HAVING THE SAME MAGNITUDE...FIX
            %THAT!
            figure, 
            hold on
            magImage = abs( flip( squeeze( SubRawData(:, indexFrame) ) ));
            phsImage = angle(angle( flip(squeeze( SubRawData(:, indexFrame) ) )));
            % FADILALI -------- KEEP THIS IFFTSHIFT-IFFTNC ARRANGMENT.  IT WILL BE THE
            % OPPOSITE FLOW.  
            plot(linspace(-SubSliceSamplingFOV/2, SubSliceSamplingFOV/2,size(SubRawData,1) ) ...
            , magImage / max(magImage(:)), 'LineWidth',2.0)


            


            plot(linspace(-SubSliceSamplingFOV/2, SubSliceSamplingFOV/2,size(SubRawData,1) ) ...
                ,phsImage / max(phsImage(:)), 'LineWidth',2.0)

            xline(3, 'LineWidth',2.0), xline(-3,'LineWidth',2.0)
            title(strcat('K Space Sub Partitions  Frame#', 32, num2str(indexFrame)))
            hold off
            legend('Magnitude', 'Phase')
            

%% PLOT THE IFFT

            indexFrame = round(PeakVelExc/intendedVPS)  ;
           % indexFrame = 1;
           % indexFrame = 20;
           %indexFrame = 2;
          % indexFrame = 250;
           % indexFrame = PeakVelExc + 80;
           % indexFrame = 220;
           %indexFrame = 55; 
          % indexFrame = round(PeakVelExc/2);
           
            %10 vps test
            %indexFrame = 16;

            %OKAY HERE IS YOUR PROBLEM...EACH LINE IS HAVING THE SAME MAGNITUDE...FIX
            %THAT!
            figure, 
            hold on
            
            maxImage = abs( flip((ifftnc( squeeze( SubRawData(:, indexFrame) ) ))) );
            phsImage = angle( flip((ifftnc( squeeze( SubRawData(:, indexFrame) ) ))) );
            % FADILALI -------- KEEP THIS IFFTSHIFT-IFFTNC ARRANGMENT.  IT WILL BE THE
            % OPPOSITE FLOW.  
            plot(linspace(-SubSliceSamplingFOV/2, SubSliceSamplingFOV/2,size(SubRawData,1) ) ...
            , ((maxImage / max(maxImage(:)))), 'LineWidth',5.0)


            


            plot(linspace(-SubSliceSamplingFOV/2, SubSliceSamplingFOV/2,size(SubRawData,1) ) ...
                ,((phsImage / max(phsImage(:)))), 'LineWidth',2.0)

            xline(3, 'LineWidth',2.0), xline(-3,'LineWidth',2.0)
            xlim([-10 10])
            title(strcat('Image Space Sub Partitions  Frame#', 32, num2str(indexFrame)))
            hold off
            legend('Magnitude', 'Phase')
            %% Save everything...
            %rawData, subRawData, and MTREndAll
            %CyclicRawData2 has a cycle which peaks at around .18 ms.

                
            save(strcat('CyclicRawData2FLOWCOMP_',num2str(squeeze(PhiPerTR(OffResIter))),'radPerTR',num2str(NPixel),'PhaseEncodes_',num2str(VPS),'VPS_StartsAt', ...
                num2str(nextFirstTargetPrep),'_',num2str(NumOfFrames),'Frames_', ...
                num2str(numOfSegments),'Segments_CoilSense2','.mat'),'rawData','-v7.3');
            
            save(strcat('CyclicSubRawData2FLOWCOMP_',num2str(squeeze(PhiPerTR(OffResIter))),'radPerTR',num2str(NPixel),'PhaseEncodes_',num2str(VPS),'VPS_StartsAt', ...
                num2str(nextFirstTargetPrep),'_',num2str(NumOfFrames),'Frames_', ...
                num2str(numOfSegments),'Segments_CoilSense2','.mat'),'SubRawData','-v7.3');
            
            save(strcat('CyclicMTREndAll2FLOWCOMP_',num2str(squeeze(PhiPerTR(OffResIter))),'radPerTR',num2str(NPixel),'PhaseEncodes_',num2str(VPS),'VPS_StartsAt', ...
                num2str(nextFirstTargetPrep),'_',num2str(NumOfFrames),'Frames_', ...
                num2str(numOfSegments),'Segments_CoilSense2','.mat'),'MTREndAll','-v7.3');
            
            %FAdil ali change this to the formatted speed arr.  
            save(strcat('CyclicFormattedSpeedArr2FLOWCOMP_',num2str(squeeze(PhiPerTR(OffResIter))),'radPerTR',num2str(NPixel),'PhaseEncodes_',num2str(VPS),'VPS_StartsAt', ...
                num2str(nextFirstTargetPrep),'_',num2str(NumOfFrames),'Frames_', ...
                num2str(numOfSegments),'Segments_CoilSense2','.mat'),'FormattedSpeedArr','-v7.3');
            
            disp('Finished saving...')
            
