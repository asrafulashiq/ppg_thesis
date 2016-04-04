function [freqEstimates,peaks] = doEEMD(sig,fPrev)
% do EEMD returns frequency estimate using EEMD algorithm


% ensemble of NE signals is created from the given signal
% by adding white Gaussian noise of 30db
NE = 5;
SNR = 30;

imfs1 = {} ; % imfs of the first channel
imfs2 = {} ; % imfs of the second channel


for i = 1:5
    
    % imf for first channel
    iSig = 2;
    sigWithNoise = awgn( sig(iSig,:), SNR, 'measured' );
    tmp_imfs = nwem( sigWithNoise );
    if length(tmp_imfs) >= 2
        imfs1{i} = tmp_imfs{2};
    end
    
    % imf for second channel
    iSig = 3;
    sigWithNoise = awgn( sig(iSig,:), SNR, 'measured' );
    tmp_imfs = nwem( sigWithNoise );
    if length(tmp_imfs) >= 2
        imfs2{i} = tmp_imfs{2};
    end
      
end

% determine ensemble average of imfs of each channel

% for channel 1
sumOfImfs = zeros(1,length(sig(2,:)));
for i = 1:length(imfs1)
    sumOfImfs = sumOfImfs + imfs1{i};
end
imfs1Average = sumOfImfs / length(imfs1);

% for channel 2
sumOfImfs = zeros(1,length(sig(3,:)));
for i = 1:length(imfs2)
    sumOfImfs = sumOfImfs + imfs2{i};
end
imfs2Average = sumOfImfs / length(imfs2);



