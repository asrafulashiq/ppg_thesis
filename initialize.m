function fPrev = initialize(sig, accData, fSampling)
% initialize estimate initial fequency
    % parameters 
    % ----------
    %   sig     : ppg signal data
    %   accData : accelaration data (x,y,z)
    
    % output
    % --------
    %   fPrev   : returns initial value of frequency
    
w=linspace(50,150,4000);
ww = 2 * pi * w / (fSampling * 60); 
fSig = abs(freqz(sig, 1, ww));

[~,locs]=findpeaks(fSig, 'MINPEAKHEIGHT', 0.8 * max(fSig),...
                    'SORTSTR', 'descend');
                
if length(locs)==1
    fPrev = w(locs);
    return;
end

% if more than one frequencies have value above 80% 
% then we look for their second harmonics

peakFrequencies = w(locs);
secondHarmonicMax = -1;

for i = 1:length(locs)
    for j = 1:length(locs)
        if i ~= j && abs( 2 * peakFrequencies(i) - peakFrequencies(j)) < 5 
            if peakFrequencies(i) > secondHarmonicMax
               secondHarmonicMax =  peakFrequencies(i);
            end         
        end
    end
end

if secondHarmonicMax ~= -1  % |2x?y| < 5 BPM found
   fPrev = secondHarmonicMax;
   return;
end

%  taking all the dominant peaks (80% of the maximum peak) 
%  from the acceleration signals 
dominantPeaksOfAcc = [];
for iAcc = 1:3
    fAccSig = abs( freqz( accData(iAcc,:), 1, ww) );
    [~,locsAcc]=findpeaks(fAccSig, 'MINPEAKHEIGHT', 0.8*max(fAccSig),...
                        'SORTSTR', 'descend');
    dominantPeaksOfAcc = [ dominantPeaksOfAcc, w(locsAcc)];
end

accCloseFreqMax = -1 ; 

for iFreqCount = 1 : length(locs)
    for jAccFreqCount = 1 : length(dominantPeaksOfAcc)      
        if abs(peakFrequencies(iFreqCount) - dominantPeaksOfAcc(jAccFreqCount)) > 5
            if peakFrequencies(iFreqCount) > accCloseFreqMax
               accCloseFreqMax = peakFrequencies(iFreqCount); 
            end
        end
    end
end

if accCloseFreqMax ~= -1
   fPrev = accCloseFreqMax;
   return;
end

% else return strongest peak?s location
[~ , loc] = max(fSig);
fPrev = w(loc);

end