clear all;close all;clc
load('~/data/mixture-res/FEATURES-timit-mask-mix-SPEC-WHATTYPE-2-100.mat');
%load('~/data/mixture-res/FEATURES-timit-mask-mix-SPEC-WHATTYPE-1-100.mat');

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

xlgnd=S.Hz_mod_cfreqs;
ylgnd=S.audio_cutoffs_Hz;
%%

mm=cell(2,1);
ss=cell(2,1);
ma=mean([features{1};features{2}]);
sa=std([features{1};features{2}]);
for t=1:2,
 for M=1:2,
     mm{M}=mean((features{M}-repmat(ma,size(features{M},1),1))./repmat(sa,size(features{M},1),1));
     ss{M}=std((features{M}-repmat(ma,size(features{M},1),1))./repmat(sa,size(features{M},1),1));
     
 end
end

figure(11);clf
cnt=0;
for M=1:2,
        cnt=cnt+1;
        subplot(4,2,cnt);
         nori_log_imagesc(xlgnd,ylgnd,reshape(mm{M},[length(ylgnd),length(xlgnd)]),[],[]);axis xy;
       title(sprintf('mean of all samples grp %d',M));
end


for M=1:2,
        cnt=cnt+1;
        subplot(4,2,cnt);
         nori_log_imagesc(xlgnd,ylgnd,reshape(mm{M}./ss{M},[length(ylgnd),length(xlgnd)]),[],[]);axis xy;
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

