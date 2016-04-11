function frq_est=frequency_estimate(p,q,ax,ay,az,prev,mul)
    frq_est1=-1;
    frq_est2=-1;
    ch1=time_domain(p,ax,ay,az,mul);
    ch2=time_domain(q,ax,ay,az,mul);

    if(length(ch1)~=1)
        w=linspace(0,200,1000);
        w=(w/(125*60*mul))*2*pi;

        yy1=ch1;
        yy1=abs(freqz(yy1,1,w)).^2;
        [faul1,loc]=findpeaks(yy1,'Sortstr','descend');

        w=w(loc);
        w=(w*125*60*mul)/(2*pi);
        [faul,loc]=min(abs(w-prev));
        frq_est1=w(loc);
    end
    if(length(ch2)~=1)

        w=linspace(0,200,1000);
        w=(w/(125*60*mul))*2*pi;
        yy2=ch2;
        yy2=abs(freqz(yy2,1,w)).^2;
        [faul2,loc]=findpeaks(yy2,'Sortstr','descend');

        w=w(loc);
        w=(w*125*60*mul)/(2*pi);
        [faul,loc]=min(abs(w-prev));
        frq_est2=w(loc)
    end

    
    frq_est=-1;


    bd=12;

    if (abs(frq_est1-prev)<=bd && abs(frq_est2-prev)<=bd)
        if( abs(frq_est1-prev) < abs(frq_est2-prev))
            frq_est=frq_est1;
        else
            frq_est=frq_est2;
        end

    elseif (abs(frq_est1-prev)<=bd)
        frq_est=frq_est1;

    elseif (abs(frq_est2-prev)<=bd)
        frq_est=frq_est2;

    end

end