
total_file_no = 13;

for fileNo =  1:total_file_no
    
    fprintf('%d',fileNo);
    
    [sig, bpm0 ] =  input_file(fileNo);
    
    NE = 8;
    SNR = 30;
    
    imfs1 = {} ; % imfs of the first channel
    imfs2 = {} ; % imfs of the second channel
    
    
    % determine which imf to chose based on sampling frequency
    imfToChose = 2;
    
    
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
    
    filename = sprinf('%d.mat',fileNo);
    
    save(filename,imfsAverage);
    
    
end

