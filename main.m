clear;
close all;

total_file_no = 13;

fSampling = 125 ; % sampling frequency of the data
multiplier = round(fSampling/125);

bpm_estimeates = zeros(total_file_no,1);
avg_error = zeros(total_file_no,1);

% write the results in a new file
% create new file
fileToSaveResult = 'result.txt';
fileID = fopen(fileToSaveResult,'w');
fprintf(fileID,'## RESULT : \n----------------\n');
fclose(fileID);

% open file for appendind
fileID = fopen(fileToSaveResult,'a');


for fileNo = 2:total_file_no
    
    fprintf(fileID,'file no : %d\n',fileNo);
    
    ii = 0;
    avg = 0;
    
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
    iStop  = length(rN) - 1000 * multiplier ;
    
    delta_count = 0 ; % need in EEMD
    
    
    iCounter = 1;
    for iSegment = iStart : iStep : iStop
        
        currentSegment = iSegment : ( iSegment + 1000 * multiplier - 1 );
        
        [freqEstimates,peaks] = doEEMD(sig(:,currentSegment),fPrev,delta_count,fSampling);
        
        % we construct Simf \3 Sa, 0.5, and from this set, we
        % take the peak location nearest to fprev
        [minimum,loc] = min(abs(peaks - fPrev));
        
        freq_td = frequency_estimate( sig(2,currentSegment),...
            sig(3,currentSegment),sig(4,currentSegment),...
            sig(5,currentSegment),sig(6,currentSegment),...
            fPrev,multiplier);
        
        if freqEstimates ~= -1
            % tracking from AC
            delta_count = 0;
            fprintf('tracking from AC : ');
            
        elseif minimum <= 7 % If its distance from fprev is within 7 BPM
            % tracking from emd 2
            freqEstimates = peaks(loc);
            if delta_count>0
                delta_count = delta_count - 0.5;
            end
            fprintf('tracking from emd : ');
        else % track from rls
            delta_count = delta_count + 1;
            
            y_cropped = rN(currentSegment);
            
            % we put its dominant peaks (80% of the maximum) in
            % Srls
            S_rls = maxFindFromThreshold(y_cropped,0.8,fSampling);
            
            % We also construct Sa, 0.6 by taking the dominant peaks
            %(60% of the maximum, a moderate threshold for tracking
            % purpose) from a? (n)
            
            S_a = [];
            for iAcc = 4:6
                
                dominantPeaks = maxFindFromThreshold(sig(iAcc,currentSegment),...
                    0.6,fSampling);
                if length(dominantPeaks)>2
                    dominantPeaks = dominantPeaks(1:2);
                end
                S_a = [S_a , dominantPeaks];
                
            end
            
            % set Srls\3 Sa,0.6
            f_rls_set = [];
            
            for iRls = S_rls
                
                if min( abs( S_a - iRls ) ) > 3
                    f_rls_set = [f_rls_set , iRls];
                end
                
            end
            
            f_rls_set = f_rls_set( f_rls_set>40 & f_rls_set<200 );
            
            if length(f_rls_set)==1 && abs(f_rls_set-fPrev)<25
                freqEstimates = f_rls_set(1);
                %fprintf('abs cause\n');
                
            elseif freq_td ~= -1 && abs(freq_td - fPrev ) < 12
                % from td
                fprintf('tracking from td : ');
                freqEstimates = freq_td;
                
            else
                % strongest peak in Srls is looked for such that it lies
                % close to fprev within a range 7 - 12 BPM
                
                f_ = maxFind(y_cropped,fSampling);
                if abs(f_ - fPrev) <= 9
                    freqEstimates = f_;
                end
                
                if freqEstimates == -1
                    
                    % If the above steps fail to provide with the crude estimate
                    % f , then we consider all the peak locations attainable
                    % from the periodograms of yi (n) and array them together
                    % in a set Sorg
                    
                    S_org = findSignalPeaks(sig(2:3,currentSegment),fPrev ,5,fSampling);
                    S_a0  = findSignalPeaks(sig(4:6,currentSegment),fPrev,10,fSampling);
                    
                    if abs(S_org - S_a0)>3
                        freqEstimates = S_org;
                    else
                        freqEstimates = fPrev;
                    end
                    
                    
                end
                
                
            end
        end
        
        [freqEstimates, ii]=call_nlms2(sig(2,currentSegment),sig(3,currentSegment),...
            multiplier,freqEstimates,4,ii);
        
        
        fPrev = freqEstimates;
        
        
        bpm_estimeates(fileNo,iCounter) = freqEstimates;
        err = freqEstimates - bpm0(iCounter);
        
        avg=(avg*(iCounter-1)+abs(err))/iCounter;
        
        fprintf('\naverage : %.2f\n',avg);
        
        fprintf(fileID,' #%d : error : %.2f , average : %.2f\n'...
            ,iCounter,err,avg);
        
        iCounter = iCounter + 1;
        
    end
    
    avg_error(fileNo) = avg;
    
    fprintf(fileID,'Average Error : %.2f \n',avg);
    fprintf(fileID,'--------------------');
    fprintf(fileID,'--------------------');
    
    
end

fclose(fileID);
