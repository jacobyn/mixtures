% compute features
close all;clc;clear all;
NFB=50; %number of subbands

TRAIN='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\timit-train';
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\');

cd (TRAIN);
ITER=1000;

files=dir('*11*.wav');
ITER=min(ITER,length(files));

%bank range
LF=40;HF=8000; 

FS0=16000; %fs (TIMIT)
FS1=1600;
WIND=50/1000; %window time in s
NW=round(WIND*FS0);
NTP=round(WIND*FS1); %temporal window size
NTPFACT=5; %delute the number of samples from each wav
BEG=round(10/1000*FS1); %remove the begining

is_show_snap=false;
%for debuging showing snapshots

Ms=[1];
MN=length(Ms);
FEATURES=cell(MN,1);
INFO=cell(MN,1);
ORIG=cell(MN,1);


tsold=-1;
for m=1:MN,
    M=Ms(m);
    INFO{m}=cell(2,1); 
    cnt=0;
    for KK=1:ITER,
        mfname=files(KK).name;
        display(sprintf('%s [%d]',mfname,size(FEATURES{m},1)));
        
        [ts,fs]=audioread(mfname);
        assert(fs==FS0);
        tsr=ts/sqrt(mean(ts.^2))*0.03; %normalize wav with rms
        
        if size(tsr,1)~=size(tsold,1)
            [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(tsr),fs,NFB-2,LF,HF); %make filter bank
            tsold=tsr;
            xlgnd=1000*(1:NTP)./FS1;
            ylgnd=SubbandFrequencies;
        end
        
        Cgrm= generate_subbands(tsr,FilterBank)';
        

        acum=[];
        acum_orig=[];
        
        for II=1:NW:(size(Cgrm,2)-NW-1),
            snap=Cgrm(:,(II):(II+NW-1));
            tt3=(II):(II+ NW-1);
            if is_show_snap
                         if mod(II,100)==1
                             imagesc(snap);
                             title(mfname)
                             
                             p=audioplayer(tsr(tt3),fs);p.play;
                             pause 
                         end
            end
            snapr=nan(size(snap(:,1),1),NTP);
            
            func=@(x)resample(x,FS1,FS0);
            snapr = blkproc(snap,[1 size(snap,2)] ,func);
%             
%             for JJ=1:size(snap(:,1),1)
%                 snapr(JJ,:)=resample(snap(JJ,:),FS1,FS0);
%             end
            assert(size(snapr,2)==NTP);
            
            snapr1=reshape(snapr,1,NTP*NFB);
            snapr2=10*log10((snapr1.*snapr1)+eps);
            acum=[acum;snapr2];
            
            cnt=cnt+1;
            INFO{m}{cnt,1}.fname=mfname;
            INFO{m}{cnt,1}.range=[min(tt3),max(tt3)];
      
        end
            FEATURES{m}=[FEATURES{m};acum];  
            assert(size(FEATURES{m},1)==length(INFO{m}))
        
    end
    
end
%%


save('FEATURES-N-v4.mat');
%%

% figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
% figure(1);imagesc([log(abs(FEATURES{2}));log(abs(FEATURES{1}))]);axis xy