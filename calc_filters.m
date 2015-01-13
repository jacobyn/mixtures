clear all;close all;clc;
  load('FEATURESv2.mat'); %high resolution
%   is_flat=false;
%  load('FEATURESv3.mat')
% load('FEATURESv5.mat') % low resolution
% load('FEATURES-v2.mat');
  is_flat=false;

%  load('FEATURES-flat-v3.mat');
%  is_flat=true;

LF=300;HF=2500;
[ts,fs]=wavread('data/mix-v1.M.1.32.wav');
[FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF);

NTRAIN=900;
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


% logcs=-10:1:10;
 logcs=0;
log2g=inf;bestg=inf
mygrid=nan(length(logcs),1);
bestcv = 0;
for I=1:length(logcs),
    display(sprintf('I=%d of %d',I,length(logcs)));
    log2c=logcs(I);
    
    %     cmd = ['-v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
    cmd = ['-v 5 -c ', num2str(2^log2c), ' -t 0 '];
    cv = svmtrain(labels, vals, cmd);
    if (cv >= bestcv),
        bestcv = cv; bestc = 2^log2c;
    end
    mygrid(I)=cv;
    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
end
figure(10);plot(logcs,mygrid)




bestc=2^5;

cmd = ['-c ', num2str(bestc), ' -t 0 '];
model = svmtrain(labels, vals, cmd);

% get w and b
w = model.SVs' * model.sv_coef;
b = -model.rho;

if model.Label(1) == -1
    w = -w;
    b = -b;
end

[predict_label, accuracy, dec_values] = svmpredict(labelsT, valsT, model); % test the training data
%%
figure(7);clf;
plot(1:length(dec_values(labelsT==2)),0*(1:length(dec_values(labelsT==2))),'g--','LineWidth',2);hold on;
plot(1:length(dec_values(labelsT==1)),dec_values(labelsT==1),'bo','MarkerFaceColor','b');hold on;
plot(1:length(dec_values(labelsT==2)),dec_values(labelsT==2),'rs','MarkerFaceColor','r');

title(sprintf('preediction accuracy=%g',accuracy(1)));
xlabel('item #'); ylabel('decision score');

%%
if ~is_flat
    wf=reshape(w,NFB,NFB);
    % bf=reshape(b,NFB,NFB);
    
    
    figure(11);clf
    
    subplot(2,2,1);
    wf_plus=abs(wf).*(wf>0);
    wf_minus=abs(wf).*(wf<0);
    
    wf_pm=wf_plus/max(max(wf_plus))-wf_minus/max(max(wf_minus));
    cmp=zeros(64,3)
    for I=1:32,
        cmp(I,:)=[0,0,1-I/32];
        cmp(32+I,:)=[I/32,0,0];
        
    end
    
    figure(12);
    imagesc(SubbandFrequencies,SubbandFrequencies,wf_pm);title('wf');axis xy
    colormap(cmp)
    
    figure(11);
    
    imagesc(SubbandFrequencies,SubbandFrequencies,wf_pm);title('wf');axis xy
    subplot(2,2,3);
    imagesc(SubbandFrequencies,SubbandFrequencies,wf_plus);title('wf plus');axis xy
    colormap('default')
    subplot(2,2,4);
    imagesc(SubbandFrequencies,SubbandFrequencies,wf_minus);title('wf minus');axis xy
    colormap('default')
    % subplot(2,2,2);
    % imagesc(SubbandFrequencies,SubbandFrequencies,abs(wf).*(wf<0));title('wf minus');
    
    % figure(5);imagesc(logcs,loggs,mygrid);axis xy;
else
    wf=w+b;
    figure(11);clf;
     subplot(2,2,1);
    plot(SubbandFrequencies,wf,'g--','LineWidth',2);hold on;
    xlabel('Fq (Hz)');
    ylabel('Filter Weight');
%     plot(SubbandFrequencies,ones(size(wf))*(-b),'g--','LineWidth',2);
    
    
     subplot(2,2,3);
%      errorbar(SubbandFrequencies,mean(valsT(labelsT==1,:)),std(valsT(labelsT==1,:)),'c-','LineWidth',2);hold on;
%     errorbar(SubbandFrequencies,mean(valsT(labelsT==2,:)),std(valsT(labelsT==2,:)),'m-','LineWidth',2);
    
    plot(SubbandFrequencies,mean(valsT(labelsT==1,:)),'b-','LineWidth',2);hold on;
    plot(SubbandFrequencies,mean(valsT(labelsT==2,:)),'r-','LineWidth',2);
    
%     subplot(2,2,3);
%     plot(SubbandFrequencies,(valsT(labelsT==1,:)),'b-','LineWidth',.5);hold on;
%     plot(SubbandFrequencies,(valsT(labelsT==2,:)),'r-','LineWidth',.5);
end
%%
figure(90);clf;    
plot(mean(valsT(labelsT==1,:)),'b-','LineWidth',2);hold on;
plot(mean(valsT(labelsT==2,:)),'r-','LineWidth',2);
