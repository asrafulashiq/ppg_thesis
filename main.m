clc;
clear;
close all;



fSampling = 125 ; % sampling frequency of the data
multiplier = round(fSampling/128);


for fileNo = 1:13
    
    [sig, bpm0 ] =  input_file(fileNo);
    
    ecgSignal  = sig(1,:); % original ecg signal
    ppgSignal1 = sig(2,:); % ppg signal 1
    ppgSignal2 = sig(3,:); % ppg signal 2
    
    % acceleration data of x,y,z
    accDataX = sig(4,:);
    accDataY = sig(5,:);
    accDataZ = sig(6,:);
    
    ppgSignalAverage = (ppgSignal1 + ppgSignal2) / 2;
    
    
    % rls filtering
    lParameterOfRls = 40; % rls filter parameter
    [~,ex] = filter(adaptfilt.rls(lParameterOfRls),accDataX,...
        ppgSignalAverage);
    [~,exy] = filter(adaptfilt.rls(lParameterOfRls),accDataY,ex);
    [~,exyz] = filter(adaptfilt.rls(lParameterOfRls),accDataZ,exy);
    
    rRaw = exyz;  % exyz can be regarded as a denoised signal rRaw(n)
    % which is assumed to have no correlation with the
    % acceleration.
    
    
    % filtering all data to frequency range for human HR
    rN       = myBandPass(rRaw,fSampling);
    accDataX = myBandPass(accDataX,fSampling);
    accDataY = myBandPass(accDataY,fSampling);
    accDataZ = myBandPass(accDataZ,fSampling);
    
    fPrev = initialize( rN(1:1000), sig(4:6,1:1000), fSampling ); % intial value of bpm
    
    
    % bandpassing all the data of signal
    filterObj = fdesign.bandpass( 70/(fSampling*60), 80/(fSampling*60),...
                400/(fSampling*60), 410/(fSampling*60), 80, 0.01, 80  );
    D = design(filterObj,'iir');
    for i=2:6
        sig(i,:)=filter(D,sig(i,:));
    end
    
    
    % now doing the emd portion 
    
    iStart = 1;
    iStep  = 250 * multiplier ;
    iStop  = length(rN);
    
    delta_count = 0 ; % need in EEMD
    
    for iSegment = iStart : iStep : iStop
       
        currentSegment = iSegment : ( iSegment + 1000 * multiplier - 1 );
        
        [freqEstimates,peaks] = doEEMD(sig(:,currentSegment),fPrev,delta_count,fSampling);
        
        % we construct Simf \3 Sa, 0.5, and from this set, we 
        % take the peak location nearest to fprev
        [minimum,loc] = min(abs(peaks - fPrev));
        
        if freqEstimates ~= -1
           % tracking from AC 
           delta_count = 0;
           fprintf('tracking from AC');
           
        elseif minimum <= 7 % If its distance from fprev is within 7 BPM
                            % tracking from emd 2
            freqEstimates = peaks(loc);
            if delta_count>0
               delta_count = delta_count - 0.5; 
            end
            fprintf('tracking from emd ');
        else % track from rls
            delta_count = delta_count + 1;
            
            y_cropped = rN(currentSegment);
            
            % we put its dominant peaks (80% of the maximum) in 
            % Srls
            S_rls = maxFindFromThreshold(y_cropped,0.8,fSampling);
            
            % We also construct Sa, 0.6 by taking the dominant peaks 
            %(60% of the maximum, a moderate threshold for tracking 
            % purpose) from a? (n)
            
             
            
            
        end
        
    end
    
    
    
end