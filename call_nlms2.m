function [z,ii]=call_nlms2(p,q,mul,frq_est,boun,iii)
    w=2*pi*linspace(frq_est-boun,frq_est+boun,4000)/(125*60*mul);
    fy1=abs(freqz(p,1,w)).^2;
    fy2=abs(freqz(q,1,w)).^2;
    [~,locs1]=findpeaks(fy1);
    [~,locs2]=findpeaks(fy2);
    locs=[locs1,locs2];
    cand=w(locs)*(125*60*mul)/(2*pi);
    if(~isempty(locs))
        [~,pos]=min(abs(cand-frq_est));
        z=cand(pos);
        ii=iii+1;
    else
        z=frq_est;
        ii=iii;
    end
end