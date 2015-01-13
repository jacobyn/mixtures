clear all;close all;clc;
  load('FEATURES-N-v3.mat'); %high resolution


Nab=2;
FEATURESab=cell(Nab,1);
FEATURESab{1}=FEATURES{1};
FEATURESab{2}=FEATURES{2};
FEATURESall=[FEATURESab{1};FEATURESab{2}];
scale_min=min(FEATURESall);
scale_max=max(FEATURESall);

%%%% scale features to range 0..1 according to the recomendation of libsvm
for I=1:Nab
    FEATURESab{I}=(FEATURESab{I}-repmat(scale_min,size(FEATURESab{I},1),1))./(repmat(scale_max,size(FEATURESab{I},1),1)-repmat(scale_min,size(FEATURESab{I},1),1));
end

addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\libsvm');
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\libplsda');

NTRAIN=round(size(FEATURESab{1},1)*9/10);

vals=[];labels=[];
for I=1:Nab
    vals=[vals;FEATURESab{I}(1:NTRAIN,:)];
    labels=[labels;ones(size(FEATURESab{I}(1:NTRAIN,:),1),1)*I];
end

valsT=[];labelsT=[];
for I=1:Nab
    valsT=[valsT;FEATURESab{I}((NTRAIN+1):end,:)];
    labelsT=[labelsT;ones(size(FEATURESab{I}((NTRAIN+1):end,:),1),1)*I];
end

Xp=vals;
yp=labels*2-3;
Ap=15;
XpT=valsT;
ypT=labelsT*2-3;
%%

LDA=plslda(Xp,yp,Ap);
figure(17);clf;figure(18);clf;
for I=1:min(16,size(LDA.W,2)),
    figure(17);subplot(4,4,I);
    wf=LDA.Xloadings(:,I);
    xlgnd=1000*(1:JMP:(WIND*FS0))./FS0;
    ylgnd=SubbandFrequencies;
    wfr=reshape(wf,[NFB,NTP]);
%     figure(11);clf;
%      subplot(2,2,1);
    imagesc(xlgnd,ylgnd,wfr);hold on;axis xy;
    ylabel('Fq (Hz)');
    xlabel('time(ms)');
    title(sprintf('%d',I));
    
    figure(18);subplot(4,4,I);
    plot(ylgnd,sum(reshape(LDA.Xloadings(:,I),[NFB,NTP]),2));
    xlabel('Fq (Hz)');
%     ylabel('weight');
    title(sprintf('marginal of comp %d',I));
    
end

figure(17);subplot(4,4,16);
plot(1-LDA.R2X,'rx--','LineWidth',2);hold on;
plot(1-LDA.R2Y,'bo-','LineWidth',2);hold on;
xlabel('coef');
ylabel('explained var');
legend('explain X','explain Y','Location','SouthEast');
figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy

figure(17);

