function [freqEstimates,peaks] = doEEMD(sig,fPrev,delta_count,fSampling)
% do EEMD returns frequency estimate using EEMD algorithm


% ensemble of NE signals is created from the given signal
% by adding white Gaussian noise of 30db
NE = 8;
SNR = 30;

imfs1 = {} ; % imfs of the first channel
imfs2 = {} ; % imfs of the second channel

freqEstimates = -1; % invalid frequency estimate return

% determine which imf to chose based on sampling frequency
imfToChose = [];
if fSampling == 25
    imfToChose = 1;
elseif fSampling == 125
    imfToChose = 2;
elseif fSampling == 250
    imfToChose = 3;
elseif fSampling == 500
    imfToChose = 4;
else
    error('Sampling frequency must be 25,125,250 or 500 Hz');
end

for i = 1:NE
    
    % imf for first channel
    iSig = 2;
    sigWithNoise = awgn( sig(iSig,:), SNR, 'measured' );
    tmp_imfs = nwem( sigWithNoise );
    if length(tmp_imfs) >= 2
        imfs1{i} = tmp_imfs{imfToChose};
    end
    
    % imf for second channel
    iSig = 3;
    sigWithNoise = awgn( sig(iSig,:), SNR, 'measured' );
    tmp_imfs = nwem( sigWithNoise );
    if length(tmp_imfs) >= 2
        imfs2{i} = tmp_imfs{imfToChose};
    end
    
end

% determine ensemble average of imfs of each channel
imfsAverage = {};

% for channel 1
sumOfImfs = zeros(1,length(sig(2,:)));
for i = 1:length(imfs1)
    sumOfImfs = sumOfImfs + imfs1{i};
end
imfsAverage{1} = sumOfImfs / length(imfs1);

% for channel 2
sumOfImfs = zeros(1,length(sig(3,:)));
for i = 1:length(imfs2)
    sumOfImfs = sumOfImfs + imfs2{i};
end
imfsAverage{2} = sumOfImfs / length(imfs2);

% two selected IMFs obtained from both channels
% are inspected in Fourier domain  and their maximum
% peaks? locations are put in a set S_imf

S_imf = zeros(1,2);
for iChannel = 1:2
    pks = maxFind(imfsAverage{iChannel},fSampling);
    S_imf(iChannel) = pks;
end

% we construct another set Sa, 0.5 by taking
% all the dominant peaks (50% of the maximum peak,
% a low threshold to ensure capturing all MA peaks)
% from the acceleration signals a?(n)
threshold = 0.5; % 50%

S_a = [];

for iAcc = 4:6
    accData = sig(iAcc,:);
    pks_vector = maxFindFromThreshold(accData,threshold,fSampling);
    if length(pks_vector) > 3
        pks_vector = pks_vector(1:3);
    end
    S_a = [S_a , pks_vector];
end


% if G is such a group, then max(G) ? min(G) ? 2 BPM
% and |G \2 Sa,0.5| ?= 0. Groupwise averages are taken
% and from these g averages, the one closest to
% fprev, say fAC is considered. If |fAC ? fprev| < ?AC with
% ?AC being a suitable range

S_imf_a_3 = [] ; % Simf \3 Sa,0.5
G_2 = [] ;

for j = 1:length(S_imf)
    
    if min( abs( S_imf(j) - S_a ) ) > 3
        S_imf_a_3 = [S_imf_a_3 S_imf(j)];
    end
    
    if min( abs( S_imf(j) - S_a ) ) > 2
        G_2  = [G_2 S_imf(j)];
    end
    
end

peaks = S_imf_a_3;

if isempty(G_2)
    %fprintf('G_2 is empty');
    return;
end

f_AC = [];
% group wise average
if max(G_2)-min(G_2) <= 2
    f_AC = mean(G_2);
else
    [~,loc] = min(abs( G_2 - fPrev ));
    f_AC = G_2(loc);
end

% if f is not assigned any value from this step m times
% in a row, then we will have for the next time window
% ?AC =?0 +m?d

delta_0 = 5;
delta_d = 1;

delta_AC = delta_0 + delta_count * delta_d;

if abs( f_AC - fPrev ) < delta_AC
    freqEstimates = f_AC;
end


