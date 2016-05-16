fSampling  = 125;

f = [0, 70/(fSampling*60) , 80/(fSampling*60),...
                           400/(fSampling*60), 410/(fSampling*60) ,1];
                       
a = [0 0 1 1 0 0];

n = 150;

h = firls(n,f,a);

fvtool(h);

save 'filter_coeff/mypass1.dat' h -ascii




% fSampling = 125;
% fcuts = [ 30/(60) , 40/(60),...
%                             400/(60), 410/(60) ];
% mags = [0 1 0];
% devs = [0.001 0.05 0.001];
% 
% [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fSampling);
% hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
% 
% freqz(hh)