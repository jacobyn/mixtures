clear all;close all;clc
%load('~/data/mixture-res/FEATURES-timit-mask-mix-SPEC-WHATTYPE-2-100.mat');
load('~/data/mixture-res/FEATURES-timit-mask-mix-SPEC-WHATTYPE-1-100.mat');

cd '/Users/jacoby/Dropbox (PPCA)/Research MIT/mixtures'
WF=load('WFR-200-MIX-WHATTYPE-A2.mat');
NS=size(INFO{M},1);
scores=nan(2,NS);
features={nan(NS,640),nan(NS,640)}
for M=1:2
    for I=1:NS
        fprintf('I=%d of %d  M=%d\n',I,NS,M);
        y=INFO{M}{I}.audio;
        fs=INFO{M}{I}.fs;
        S = nori_measure_texture_stats(y, fs); %measure stats
        S=S.S;
        score=sum(sum(S.mod_power.*WF.wfr));
        scores(I,M)=score;
        
        S = nori_measure_texture_stats(y, fs); %measure stats
        
        
        S=S.S;
        
        feature=reshape(S.mod_power,[1 (length(S.Hz_mod_cfreqs)*length(S.audio_cutoffs_Hz))]); % rehshape stats
        features{M}(I,:)=feature;
    end
    
    
end
%%
figure(1);clf;
plot(scores(:,1),'bo');hold on;
plot(scores(:,2),'rx');
NAMES={'single','mixtures'};
ANAMES=cell(2,2);
set(gca,'FontSize',14);
legend(NAMES,'FontSize',14);
HM=20;
all_sound=cell(2,2);
featuresE=cell(2,2);
for M=1:2,
    
    
    
    fprintf('%s that are most singles\n',NAMES{M});
    ANAMES{M,1}=sprintf('%s that are most singles\n',NAMES{M});
    
    [vals,idxs]=sort(scores(:,M));
    
    
    for II=1:HM,
        I=idxs(II);
        
        
        y=INFO{M}{I}.audio;
        fs=INFO{M}{I}.fs;
        %nori_doplay(y,fs);
        all_sound{M,1}=[all_sound{M,1};y;zeros(round(fs/0.5),1)];
        featuresE{M,1}(II,:)=features{M}(I,:);
    end
    %pause
    
    fprintf('%s that are most mixtures\n',NAMES{M});
    ANAMES{M,2}=sprintf('%s that are most mixtures\n',NAMES{M});
    
    for II=1:HM,
        I=idxs(end-II);
        
        y=INFO{M}{I}.audio;
        fs=INFO{M}{I}.fs;
        %nori_doplay(y,fs);
        all_sound{M,2}=[all_sound{M,2};y;zeros(round(fs/0.5),1)];
        featuresE{M,2}(II,:)=features{M}(I,:);
    end
end

mmm=cell(2,2);
sss=cell(2,2);

for t=1:2,
 for M=1:2,
     mmm{M,t}=mean(featuresE{M,t});
     sss{M,t}=std(featuresE{M,t});
     
 end
end
   

mm=cell(2,1);
ss=cell(2,1);

for t=1:2,
 for M=1:2,
     mm{M}=mean(features{M});
     ss{M}=std(features{M});
     
 end
end
   
%%
for M=1:2,
    for t=1:2,
        %nori_doplay(all_sound{M,t},fs);
        fname=sprintf('MAR15c-grp-%d-max-%d.wav',M,t)
        audiowrite(fname,all_sound{M,t},fs)
    end
end
%%
xlgnd=S.Hz_mod_cfreqs;
ylgnd=S.audio_cutoffs_Hz;

M=1;
figure(10);clf;
for M=1:2,
    subplot(2,2,M);
    imagesc(reshape(mean(features{M}),[length(ylgnd),length(xlgnd)]));
end

subplot(2,2,3);
imagesc(reshape((mean(features{2})-mean(features{1}))./(std(features{1})+std(features{2})),[length(ylgnd),length(xlgnd)]));


%%
figure(11);clf
cnt=0;
for M=1:2,
    for t=1:2,
        cnt=cnt+1;
        subplot(4,2,cnt);
        
        %imagesc(reshape(mean(featuresE{M,t}),[length(ylgnd),length(xlgnd)]));axis xy;
        imagesc(reshape(mmm{M,t}./sss{M,t},[length(ylgnd),length(xlgnd)]));axis xy;
        title(ANAMES{M,t});
        %imagesc(reshape(mean(featuresE{M,t})./std(featuresE{M,t}),[length(ylgnd),length(xlgnd)]));axis xy;
        
        
        
        %nori_log_imagesc(xlgnd,ylgnd,reshape((m2-m1)./(s1+s2),[length(ylgnd),length(xlgnd)]),[],[]);
        %nori_log_imagesc(xlgnd,ylgnd,reshape((m2-m1)./(s1+s2),[length(ylgnd),length(xlgnd)]),[],[]);
    end
end

cnt=cnt+1;
subplot(4,2,cnt);title('filter');
nori_log_imagesc(xlgnd,ylgnd,reshape(WF.wfr,[length(ylgnd),length(xlgnd)]),[],[]);

cnt=cnt+1;
subplot(4,2,cnt);title('extreme mixtures classified as mixtures');
nori_log_imagesc(xlgnd,ylgnd,reshape((mmm{2,2}./sss{2,2}-mmm{1,1}./sss{1,1}),[length(ylgnd),length(xlgnd)]),[],[]);
%%
figure(11);clf
cnt=0;
for M=1:2,
        cnt=cnt+1;
        subplot(4,2,cnt);
         imagesc(reshape(mm{M},[length(ylgnd),length(xlgnd)]));axis xy;
       title(sprintf('mean of all samples grp %d',M));
end


for M=1:2,
        cnt=cnt+1;
        subplot(4,2,cnt);
         imagesc(reshape(mm{M}./ss{M},[length(ylgnd),length(xlgnd)]));axis xy;
        title(sprintf('mean/std of all samples grp %d',M));
end
cnt=cnt+1;
subplot(4,2,cnt);
nori_log_imagesc(xlgnd,ylgnd,reshape(WF.wfr,[length(ylgnd),length(xlgnd)]),[],[]);title('filter');

cnt=cnt+1;
subplot(4,2,cnt);
nori_log_imagesc(xlgnd,ylgnd,reshape((mm{2}-mm{1})./(ss{2}+ss{1}),[length(ylgnd),length(xlgnd)]),[],[]);title('d prime');

cnt=cnt+1;
subplot(4,2,cnt);
nori_log_imagesc(xlgnd,ylgnd,reshape((mm{2}./ss{2}-mm{1}./ss{1}),[length(ylgnd),length(xlgnd)]),[],[]);title('mean/std {2}- mean/std{1}');


%%

    cnt=cnt+1;
    subplot(4,2,cnt);
    nori_log_imagesc(xlgnd,ylgnd,reshape((mmm{1,2}-mmm{1,1})./(sss{1,2}+sss{1,1}),[length(ylgnd),length(xlgnd)]),[],[]);
 
    cnt=cnt+1;
    subplot(4,2,cnt);
    nori_log_imagesc(xlgnd,ylgnd,reshape((mmm{2,1}-mmm{1,1})./(sss{2,1}+sss{1,1}),[length(ylgnd),length(xlgnd)]),[],[]);
    
    cnt=cnt+1;
    subplot(4,2,cnt);
    nori_log_imagesc(xlgnd,ylgnd,reshape((mmm{2,1}-mmm{1,1})./(sss{2,1}+sss{1,1})-(mmm{2,2}-mmm{1,2})./(sss{2,2}+sss{1,2}),[length(ylgnd),length(xlgnd)]),[],[]);

cnt=cnt+1;
subplot(4,2,cnt);
nori_log_imagesc(xlgnd,ylgnd,reshape(WF.wfr,[length(ylgnd),length(xlgnd)]),[],[]);
%%
%wfr=reshape(wf,[length(ylgnd),length(xlgnd)]);
%     imagesc(xlgnd,ylgnd,wfr);axis xy;


%data=nori_cell_array_unvectorize(wf,vecformat);is_colormap_negative=true;
clim=[-max(abs(wf)),max(abs(wf))];
%clim=[-1,1];
%clim=[];
%clim=[-5*sqrt(mean(wf.*wf)),5*sqrt(mean(wf.*wf))];
nori_log_imagesc(xlgnd,ylgnd,wfr,clim,10)
