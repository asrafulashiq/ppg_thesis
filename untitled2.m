% 
% b = [1,2,3,1,4];
% w = 0:4;
% n = length(b);
% 
% h = [];
% 
% e = exp(1j * w);
% 
% h = abs(polyval(b,e)./exp(1j * w * (n-1)))

% rng(0);
% SNR = 30;
% a = 10^(SNR/10);
% 
% t = 1:10;
% s = sin(t);
% 
% p = sum(s.^2)/length(s)
% 
% sa = s + sqrt(1/a * p) * randn(1,length(t))
