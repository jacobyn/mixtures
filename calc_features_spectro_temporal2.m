
% compute features
close all;clc;clear all;
NFB=50; %number of subbands
Ms=[1,2]; %num of spks in mixture Ms are vector of all m=#spk to compute
ITER=1000;

datadir='data-inst\';
outdir='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\';
cd (outdir);

addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\');
%bank range
% LF=300;HF=2500;
LF=40;HF=8000; 

FS0=16000; %fs (TIMIT)
WIND=50/1000; %window time in s
% WIND=5/1000; %window time in s
% WIND=.5/1000; %window time in s
JMP=10; % take every JMP samples (reduce correlation, dimensionality)
% we span a window of WIND sec every JMP samples 

NTP=length(1:JMP:(WIND*FS0)); %temporal window size
NTPFACT=5; %delute the number of samples from each wav
BEG=round(200/1000*FS0); %remove the begining
% BEG=1;


is_show_snap=false; %for debuging showing snapshots

MN=length(Ms);
FEATURES=cell(MN,1);
INFO=cell(MN,1);


tsold=-1;
for m=1:MN,
    M=Ms(m);
    INFO{m}=cell(2,1);
    cnt=0;
    
    fnames=dir(sprintf('%s*.M.%d.*',datadir,m));
    for KK=1:ITER,
        mfname=sprintf('%s%s',datadir,fnames(KK).name);
        display(sprintf('%s',mfname));
        
        [ts,fs]=audioread(mfname);
        assert(fs==FS0);
        ts=ts/sqrt(mean(ts.^2))*0.03; %normalize wav with rms
        
        if size(ts,1)~=size(tsold,1)
            [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF); %make filter bank
            tsold=ts;
            xlgnd=1000*(1:JMP:(WIND*FS0))./FS0;
            ylgnd=SubbandFrequencies;
        end
        
        Cgrm= generate_subbands(ts,FilterBank)';
        tt1=1:size(Cgrm,2);        
        Cgrm_jmp=Cgrm(:,BEG:JMP:end);
        tt2=tt1(BEG:JMP:end);
        
        acum=[];
        
        for II=1:NTP*NTPFACT:(size(Cgrm_jmp,2)-NTP),
            snap=Cgrm_jmp(:,(II):(II+NTP-1));
            tt3=tt2((II):(II+NTP-1));
            if is_show_snap
                         if mod(II,100)==1
                             imagesc(snap);title(mfname)
                             pause
                         end
            end
            snapr=reshape(snap,1,NTP*NFB);
            snapr=10*log10((snapr.*snapr)+eps);
%             snapr=snapr/(max(abs(snapr)));
            acum=[acum;snapr];
            cnt=cnt+1;
            INFO{m}{cnt,1}.fname=mfname;
            INFO{m}{cnt,1}.range=[min(tt3),max(tt3)];
      
        end
            FEATURES{m}=[FEATURES{m};acum];  
        
    end
end
cd (outdir);
% save('FEATURES-N-v4.mat');
save('FEATURES-NI-v1.mat');

figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
% figure(1);imagesc([log(abs(FEATURES{2}));log(abs(FEATURES{1}))]);axis xy