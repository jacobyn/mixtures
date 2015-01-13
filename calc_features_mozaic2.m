
% compute features
close all;clc;clear all;
NFB=50; %number of subbands

TRAIN='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\timit-train';
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\');

cd (TRAIN);
ITER=1000;

files=dir('*122*.wav');
ITER=min(ITER,length(files));

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

Ms=[1];
MN=length(Ms);
FEATURES=cell(MN,1);
INFO=cell(MN,1);
ORIG=cell(MN,1);

FEATURES2=cell(MN,1);
INFO2=cell(MN,2);
ORIG2=cell(MN,2);


tsold=-1;
for m=1:MN,
    M=Ms(m);
    INFO{m}=cell(2,1);
    cnt=0;
    for KK=1:ITER,
        mfname=files(KK).name;
        display(sprintf('%s [%d]',mfname,length(FEATURES{m})));
        
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
        acum_orig=[];
        
        for II=1:NTP*NTPFACT:(size(Cgrm_jmp,2)-NTP),
            snap=Cgrm_jmp(:,(II):(II+NTP-1));
            tt3=tt2((II):(II+NTP-1));
            if is_show_snap
                         if mod(II,100)==1
                             imagesc(snap);
                             title(mfname)
                             pause
                         end
            end
            snapr1=reshape(snap,1,NTP*NFB);
            snapr=10*log10((snapr1.*snapr1)+eps);
%             snapr=snapr/(max(abs(snapr)));
            acum=[acum;snap];
            acum_orig=[acum_orig;snapr1];
            
            cnt=cnt+1;
            INFO{m}{cnt,1}.fname=mfname;
            INFO{m}{cnt,1}.range=[min(tt3),max(tt3)];
      
        end
            FEATURES{m}=[FEATURES{m};acum];  
            ORIG{m}=[ORIG{m};acum];  
        
    end
    
end
%%
% m=1;
% FEATURES2{m}=nan(size(FEATURES{m},1).^2,size(FEATURES{m},2));
% ORIG2{m,1}=nan(size(FEATURES{m},1).^2,size(FEATURES{m},2));
% ORIG2{m,2}=nan(size(FEATURES{m},1).^2,size(FEATURES{m},2));
% INFO2{m,1}=nan(size(FEATURES{m},1).^2,size(FEATURES{m},2));
% INFO2{m,2}=nan(size(FEATURES{m},1).^2,size(FEATURES{m},2));
% 
% for I=1:size(FEATURES{m},1),
%         display(sprintf('%d %d \t %3.3g',I,size(FEATURES{m}),I*100/size(FEATURES{m})));
%     for J=1:size(FEATURES{m},1),
%         snap1=ORIG{m}(I,:);
%         snap2=ORIG{m}(J,:);
%         snap=snap1+snap2;
%         snapr=10*log10((snap.*snap)+eps);
%         FEATURES2{m}((I-1)*size(FEATURES{m},1)+J+1,:)=snapr;
%         ORIG2{m,1}((I-1)*size(FEATURES{m},1)+J+1,:)=ORIG{m}(I,:);
%         ORIG2{m,2}((I-1)*size(FEATURES{m},1)+J+1,:)=ORIG{m}(J,:);
%         INFO2{m,1}((I-1)*size(FEATURES{m},1)+J+1,:)=INFO{m}(I,:);
%         INFO2{m,2}((I-1)*size(FEATURES{m},1)+J+1,:)=INFO{m}(J,:);
%     end
% end


save('FEATURES-N-v4.mat');
%%

% figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
% figure(1);imagesc([log(abs(FEATURES{2}));log(abs(FEATURES{1}))]);axis xy