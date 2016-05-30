function y=time_domain(x,bx,by,bz,mul)

    en=bx.^2+by.^2+bz.^2;

    tt=.1*max(en);

    bin=en>tt;

    cl=0;
    mxl=0;
    mxst=0;
    st=1;
     
    for j=1:200*5*mul
         if ~bin(j)
             cl=cl+1;
             if(cl>=mxl)
                 mxl=cl;
                 mxst=st;
             end
         else
             cl=0;
             st=j+1;
         end
     end
     if(mxst+mxl-1<=200*5*mul)
         I=mxst:mxst+mxl-1;
        if (length(I)>80*5*mul) 
            kkk=I(end);
            I=((kkk-400*mul+1):kkk); 
        end

        if    (  ( length(I)>27*5*mul && I(end)>150*5*mul ) || length(I)>60*5*mul  )
            y=x(I);
        else
            y=-1;
        end
     else
         y=-1;
     end
     
end




     
     
   
