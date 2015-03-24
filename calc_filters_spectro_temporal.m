clear all;close all;clc;
  load('FEATURES-N-v1.mat'); %high resolution



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
%%

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

%%
%  logcs=-1:1:3;
  logcs=0;
log2g=inf;bestg=inf;
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



%%
 bestc=8;

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

title(sprintf('prediction accuracy=%g testsize=%d trainsize=%d features=%d',accuracy(1),length(labelsT),length(labels),size(vals,2)));
xlabel('item #'); ylabel('decision score');

%%
%     wf=w+b;
    wf=w;
    wfr=reshape(wf,[NFB,NTP]);
    figure(11);clf;
     subplot(2,2,1);
    imagesc(wfr);hold on;axis xy;
    ylabel('Fq (Hz)');
    xlabel('time(s)');


    
%%
figure(90);clf;    
subplot(3,1,1);
plot(mean(valsT(labelsT==1,:)),'b-','LineWidth',2);hold on;
plot(mean(valsT(labelsT==2,:)),'r-','LineWidth',2);
subplot(3,1,2);
plot(mean(valsT(labelsT==1,:))-mean(valsT(labelsT==2,:)),'g-','LineWidth',2);

subplot(3,1,3);
plot(wfr,'g-','LineWidth',2)
% linkaxes(ax);