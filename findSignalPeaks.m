function fPeaks = findSignalPeaks(Y,fPrev,bound,fSampling)

multiplier = fSampling/125;
w=2*pi*linspace(fPrev-bound,fPrev+bound,100)/(fSampling*60);

l=size(Y,1);
locs=[];
for i=1:l
    fy=abs(freqz(Y(i,:),1,w)).^2;
    [~,locs1]=findpeaks(fy);
    locs=[locs,locs1];
end
cand=w(locs)*(125*60*multiplier)/(2*pi);
if(~isempty(locs))
    [~,pos]=min(abs(cand-fPrev));
    fPeaks=cand(pos);
else
    fPeaks=fPrev;
end

end