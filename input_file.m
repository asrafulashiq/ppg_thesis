function [sig,BPM0] = input_file(fileNo)
% input_file returns the signal and bpm of file
    % [sig,bpm] = input_file(fileNo) , where fileNo
    % is the number of the file from 1 to 13
    
string1='Training_data/DATA_0';
string2='Training_data/DATA_';
string3='_TYPE02.mat';
string4='_TYPE02_BPMtrace.mat';
sp1='Training_data/DATA_01_TYPE01.mat';
sp2='Training_data/DATA_01_TYPE01_BPMtrace.mat';
if fileNo==1
    datafile=sp1;
    ansfile=sp2;
elseif fileNo>=10
    datafile=[string2 num2str(fileNo) string3];
    ansfile=[string2 num2str(fileNo) string4];
else
    datafile=[string1 num2str(fileNo) string3];
    ansfile=[string1 num2str(fileNo) string4];
end
if exist(datafile,'file')==0 || exist(ansfile,'file')==0
   error('datafile or ansfile of %d does not exist',fileNo);
end

load(datafile); %The main data file, returns sig ;
load(ansfile); % calculated BPM Data file, returns BPM0 ;

% load('Original_JOSS.mat');
% BPM=Original(fileNo,find(Original(fileNo,:)~=0));
% shift=length(BPM0)-length(BPM);
% sig=sig(:,(1+shift*250):length(sig));
% BPM0=BPM;


end