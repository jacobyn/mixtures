% assert(1==0);
%clear all;close all;clc;
%close all
%rng(8000,'twister') % For reproducibility

DO_EXMP=true;
fprintf('reading and formatting data...\n');

tic
%load('~/data/mixture-res/FEATURES-timit-modpower-25.mat');
load('~/data/mixture-res/FEATURES-timit-envC-25.mat');
toc
fprintf('formatting...\n');
tic
Nab=2;
FEATURESab=cell(Nab,1);
A=1;B=2;
FEATURESab{1}=FEATURES{A};
FEATURESab{2}=FEATURES{B};
FEATURESab{1}(isnan(FEATURESab{1}))=0;
FEATURESab{2}(isnan(FEATURESab{2}))=0;

INFOab=cell(Nab,1);
INFOab{1}=INFO{A};
INFOab{2}=INFO{B};



% 1 Env Mean 	 start 1 end 33
% 2 Env StdDev/Mean 	 start 33 end 65
% 3 Env Skewness 	 start 65 end 97
% 4 Mod. C2 	 start 97 end 481
% 5 Orig. Mod. Power 	 start 481 end 1121
% 6 Env C 	 start 1121 end 2145
% 7 Mod C1 (3.1 Hz) 	 start 2145 end 3169
% 8 Mod C1 (6.2 Hz) 	 start 3169 end 4193
% 9 Mod C1 (12.5 Hz) 	 start 4193 end 5217
% 10 Mod C1 (25.0 Hz) 	 start 5217 end 6241
% 11 Mod C1 (50.0 Hz) 	 start 6241 end 7265
% 12 Mod C1 (100.0 Hz) 	 start 7265 end 8289


if ispc()
    % addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\libsvm');
    addpath(genpath('C:\Users\user\Dropbox (PPCA)\Research MIT\toolboxes'));
elseif ismac()
    addpath(genpath('~/ResearchMIT/toolboxes/'))
    
end

NTRAIN=round(size(FEATURESab{1},1)*9/10);

vals=[];labels=[];
for I=1:Nab
    vals=[vals;FEATURESab{I}(1:NTRAIN,:)];
    labels=[labels;ones(size(FEATURESab{I}(1:NTRAIN,:),1),1)*I];
    
end

info=cell(length(labels),1);cnt=0;
for I=1:Nab
    for J=1:NTRAIN
        cnt=cnt+1;
        info{cnt,1}.fname=INFOab{I}{J}.fname;
        info{cnt,1}.range=INFOab{I}{J}.range;
        
        if isfield(INFOab{I}{J,1},'fs')
            fs=INFOab{I}{J,1}.fs;
        else
            fs=16000;
        end
        
        info{cnt,1}.audio=INFOab{I}{J}.audio;
        info{cnt,1}.fs=INFOab{I}{J}.fs;
        
        
        if isfield(INFOab{I}{J},'xlgnd')
            xlgnd=INFOab{I}{J}.xlgnd;
            xlgnd_name=INFOab{I}{J}.xlgnd_name;
            
            
        end
        
        if isfield(INFOab{I}{J},'ylgnd')
            ylgnd=INFOab{I}{J}.ylgnd;
            ylgnd_name=INFOab{I}{J}.ylgnd_name;
            
            
        end
        
        
    end
end
assert(length(labels)==length(info));
toc
%%
valsN=vals;

[vals,MU,SIGMA] = zscore(vals);


valsT=[];labelsT=[];
for I=1:Nab
    valsT=[valsT;FEATURESab{I}((NTRAIN+1):end,:)];
    labelsT=[labelsT;ones(size(FEATURESab{I}((NTRAIN+1):end,:),1),1)*I];
end
valsNT=valsT;

valsT=(valsT-repmat(MU,size(valsT,1),1))./repmat(SIGMA,size(valsT,1),1);
valsT(isnan(valsT))=0;

infoT=cell(length(labelsT),1);cnt=0;
for I=1:Nab
    for J=(NTRAIN+1):(size(FEATURESab{I},1))
        cnt=cnt+1;
        infoT{cnt,1}.fname=INFOab{I}{J}.fname;
        infoT{cnt,1}.range=INFOab{I}{J}.range;
        
        
    end
end
assert(length(labelsT)==length(infoT));
%%
Xp=vals;
yp=labels*2-3;
Ap=1; %maximal number of dimensions in projections
XpT=valsT;
ypT=labelsT*2-3;

% disp(sprintf('doing cross-validation...'));
% tic
% CV=plscv(Xp,yp,Ap);
% toc

% obj = fitcdiscr(vals,labels);

% disp(sprintf('doing pls pn all data...'));
% tic
% PLS=pls(Xp,yp,Ap);
% acc=sum(sign(CV.Ypred(:,:))==repmat(yp,1,size(CV.Ypred,2)))/length(yp)
% toc
% fprintf('doing lin-disc\n');
% tic
% obj = fitcdiscr(vals,labels);
Mdl = fitcdiscr(vals,labels,'SaveMemory','on','FillCoeffs','off');

[err,gamma,delta,numpred] = cvshrink(Mdl,'NumGamma',29,'NumDelta',29,'Verbose',2);
minerr = min(min(err));
%%
NUM_COEF=200;
low_limit = min(min(err(numpred <= NUM_COEF)));
[p,q] = find((err == low_limit) & (numpred <= NUM_COEF),1);
minerr=err(p,q);
mygamma= gamma(p);
mydelta=delta(p,q);
mynumpred=numpred(p,q);
obj = fitcdiscr(vals,labels,'SaveMemory','on','FillCoeffs','on','Gamma',mygamma,'delta',mydelta);
fprintf('best gamma %3.3g best delta %3.3g err %3.3g numpred %d\n',mygamma,mydelta,minerr,mynumpred);

fprintf('predict on train set...\n');
tic
[prdc,score] = predict(obj,vals);
linCoeffs=obj.Coeffs(2,1).Linear;
acc=sum(prdc==labels)/length(labels);
toc;
fprintf('predict on test set...\n');
tic
[prdcT,scoreT] = predict(obj,valsT);
accT=sum(prdcT==labelsT)/length(labelsT);
toc
fprintf('acc on train %g acc on test %g\n',acc,accT);
%%
figure (30);clf;
indx = repmat(1:size(delta,2),size(delta,1),1);
subplot(2,2,1)
imagesc(err);hold on;
plot(p,q,'ro','MarkerFaceColor','k');
colorbar;
colormap('jet')
title 'Classification error';
xlabel 'Delta index';
ylabel 'Gamma index';

subplot(2,2,2)
imagesc(numpred);hold on;
plot(p,q,'ro','MarkerFaceColor','k');
colorbar;
title 'Number of predictors in the model';
xlabel 'Delta index' ;
ylabel 'Gamma index' ;

subplot(2,2,3);
plot(err,numpred,'k.');hold on;
plot(err(p,q),numpred(p,q),'ro');

xlabel('Error rate');
ylabel('Number of predictors');
subplot(2,2,4);
b=bar([acc,accT]*100);b(1).FaceColor='blue';
ylim([50 ,100]);




%%

fprintf('showing filters...\n');
cmp=zeros(64,3);
for I=1:32,
    cmp(32+I,:)=[0,0,I/32];
    cmp(I,:)=[1-I/32,0,0];
    
end
for I=1:min(15,Ap),
    figure(190+I);
    %     subplot(4,4,I);
    wf=linCoeffs;
    %wf=nan(size(wf));wf(1:1000)=A(1:1000);wf(isnan(wf))=0;
    
    %     xlgnd=audio_cutoffs_Hz;
    %     ylgnd=Hz_mod_cfreqs;
    wfr=reshape(wf,[length(ylgnd),length(xlgnd)]);
    
    %     imagesc(xlgnd,ylgnd,wfr);axis xy;
    
    
    %data=nori_cell_array_unvectorize(wf,vecformat);is_colormap_negative=true;
    clim=[-max(abs(wf)),max(abs(wf))];
    %clim=[-1,1];
    %clim=[];
    %clim=[-5*sqrt(mean(wf.*wf)),5*sqrt(mean(wf.*wf))];
    nori_log_imagesc(xlgnd,ylgnd,wfr,clim)
    
    
    colormap(cmp)
    xlabel(xlgnd_name);
    ylabel(ylgnd_name);
    %title(sprintf('%d',I));
    
end


%%
%title(sprintf('%d',I));



figure(18);clf;
subplot(2,2,1);
ylabel('Accuracy (%%)');
h=legend('cv-train','all-test','Location','NorthEastOutside');
set(h,'FontSize',8);


plot(acc*100,'kx-','LineWidth',2);hold on;
plot([1 Ap],[accT,accT]*100,'gd--','LineWidth',2);hold on;


% plot([15 15],[0,100],'r--','LineWidth',1);hold on;
axis([0 2 50 100]);
title(sprintf('Labeling accuracy %g %g',acc*100,accT*100));
xlabel('Number of coefficents');


ylabel('Accuracy (%%)');
h=legend('cv-train','all-test','Location','SouthEast');
set(h,'FontSize',8);


subplot(2,2,2);
scores=Xp*wf;
scoresT=XpT*wf;



plot([scores;scoresT]);hold on;
plot(1:length([scores;scoresT]),0*([scores;scoresT]),'g-');
pos=find((labels==1));
plot(pos,mean(scores(pos))*ones(size(pos)),'m-');
pos=find((labels~=1));
plot(pos,mean(scores(pos))*ones(size(pos)),'m-');
pos=find((labelsT==1));
plot(length(scores)+pos,mean(scoresT(pos))*ones(size(pos)),'m-');
pos=find((labelsT~=1));
plot(length(scores)+pos,mean(scoresT(pos))*ones(size(pos)),'m-');%%
fprintf('displaying raw data as a matrix...\n');

figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
%%
figure(20);

if DO_EXMP
    WHO=[1,2];
    for J=1:length(WHO),
        
        MAXSOUNDS=9;
        fprintf('displaying examples of stimuli\n');
        for I=1:MAXSOUNDS,
            %figure(J*10000+I);
            subplot(floor(sqrt(MAXSOUNDS)+0.5),floor(sqrt(MAXSOUNDS)+0.5),I)
            idx=find(labels==WHO(J));
            
            nori_log_imagesc(xlgnd,ylgnd, (reshape(vals(idx(I),:),[length(ylgnd),length(xlgnd)]))',[]);axis xy;
            %         data=nori_cell_array_unvectorize(info{idx(idx(I))}.jstats.data,vecformat);is_colormap_negative=false;
            %nori_figure_stat_summary_as_cell_array(info{idx(I)}.jstats.data,info{idx(I)}.jstats.xleg,info{idx(I)}.jstats.yleg,info{idx(I)}.jstats.xlabels,info{idx(I)}.jstats.ylabels,info{idx(I)}.jstats.titles,[],[]);
            
            %         set(gca,'Xtick',xlgnd(1:4:end));
            %         set(gca,'XtickLabel',round(xlgnd(1:4:end)));
            %         set(gca,'Ytick',ylgnd(1:5:end));
            %         set(gca,'YtickLabel',round(ylgnd(1:4:end)));
        end
    end
end
%%




%%
if DO_EXMP
    fprintf('concatenating typical audio examples...\n');
    
    figure(2);clf;
    
    hold on;
    % WHO=[4]
    MAXSOUNDS=16; 
    Nplay=16;
    NSHOW=16;
    I1=1;
    
    [val,idx]=sort(scores);
    Nplay=400;Njump=1;
    RANGES={[1:Njump:Nplay],[length(idx):-Njump:(length(idx)-Nplay)]};
    %pos=(length(idx)/2):(Nplay+length(idx)/2); RANGES={pos(labels(idx(length(idx)/2:(Nplay+length(idx)/2)))==1),pos(labels(idx(length(idx)/2:(Nplay+length(idx)/2)))==2)};
    
    
    for k=1:length(RANGES),
        
        playme=zeros(fs,1);
        
        cnt=0;
        hashy=inf;
        MYHASHES=[];
        for J=RANGES{k}
            if cnt>=MAXSOUNDS;
                break
            end
            hashy=sum(double(info{idx(I)}.fname));
            if isempty(find(MYHASHES==hashy, 1));
                cnt=cnt+1
                NYHASHES=[MYHASHES,hashy];
                %hashy= sum(vals(idx(J),:));
                %nori_figure_stat_summary_as_cell_array(info{idx(idx(J))}.jstats.data,xleg,yleg,xlabels,ylabels,titles,[],[]);
                if cnt<=NSHOW
                    subplot(floor(sqrt(MAXSOUNDS)+0.5),floor(sqrt(MAXSOUNDS)+0.5),min(cnt,MAXSOUNDS));
                    nori_log_imagesc(xlgnd,ylgnd,reshape(vals(idx(J),:),[length(ylgnd),length(xlgnd)]),[]);axis xy;
                    axis xy;
                    %                 axis off;
                    
                end
                
                
                myrng=info{idx(J)}.range;
                
                
                Y=info{idx(J)}.audio;
                fs=info{idx(J)}.fs;
                
                
                Yh=Y;
                playmeold=playme;
                playme=[playme;Yh];
                
            end
        end
        
        
        drawnow;
        
        ofname=sprintf('NLIN-10ktimi15k-onematirx-%d.direc-%d.wav',I1,k);
        
        audiowrite(ofname,playme,fs)
        pause(0.8*length(playme)/fs);
         pause(2);
    end
end

toc
%fprintf('saving\n');
%clear myFEATURES;clear FEATURES; clear FEATURESab; clear myINFO; clear info; clear INFO; clear INFOab;

fprintf('DONE!\n');
%%
