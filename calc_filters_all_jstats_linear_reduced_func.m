function calc_filters_all_jstats_linear_reduced_func(LN_SELECT,moutfname)
% assert(1==0);
%clear all;close all;clc;
%close all
%rng(8000,'twister') % For reproducibility

FEATNUM=[10,20,30,50,100,200,500,1000];
DO_PLOT=false;
DO_EXMP=false;
fprintf('reading and formatting data...\n');
% load('FEATURES-N-v4.mat');
%  load('FEATURES-NP-v2.mat'); %timit mixture speech
%  load('FEATURES-cello-v1.mat');
%  load('FEATURES-AUDIO-BOOKS-ALL')
%     load('FEATURES-AUDIO-BOOKS-ENGLISH.mat');
%    load('FEATURES-AUDIO-BOOKS-MAND.mat');
%   load('FEATURES-AUDIO-BOOKS-ITAL.mat');
%  load('FEATURES-AUDIO-BOOKS-TOG.mat'); % 4 books together
% load('FEATURES-AUDIO-BOOKS-HUN.mat');
%  load('FEATURES-CELLO2.mat');
%  load('FEATURES-FLUTE.mat');
%  load('FEATURES-PIANO.mat');
% load('FEATURES-3INST.mat'); % 3 instrument togehter.
%  load('FEATURES-C1-1.mat');
% load('FEATURES-C1-2.mat');
% load('FEATURES-C1-3.mat');
% load('FEATURES-C1-4.mat');
%  load('FEATURES-C1-5.mat');
%   load('FEATURES-C1-6.mat');
%     load('FEATURES-C.mat');
%    load('FEATURES-CELLO-C.mat');
%     load('FEATURES-C2.mat');
%  load('FEATURES-CELLO-C2.mat');
%  load('FEATURES-PIANO-C2.mat');
%load('FEATURES-test.mat');
tic
%load('FEATURES-ALLJSTATS-ENGLISH-10k.mat');
%load('FEATURES-timit-10kb.mat');
%load('FEATURES-timit-norm-10k.mat');

%load('FEATURES-cello-10k.mat');
%%%%%load('~/data/mixture-res/FEATURES-timit-jstat-RMS15-10k.mat');
%load('~/data/mixture-res/FEATURES-timit-jstat-fail80-10kb.mat');
fname='~/data/mixture-res/FEATURES-timit-jstat-25db.mat';

load(fname);
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

%LN_SELECT=[1:12];

%LNA=1;LNB=2144;
%LNA=1;LNB=32;


%LN_RANGE=LNA:LNB;
%LN_RANGE=1:8288;
%LN_RNAGE=[1:481,1121:2145];
%LN_SELECT=[1 2 3 4 6];
%LN_SELECT=[1:12];
%LN_RNAGE=[1:33,1121:2145];




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
LNA=1;LNB=8288;
LN_ALL=zeros(8288,1);
LN_ALL(1:32)=1;
LN_ALL(33:64)=2;
LN_ALL(65:96)=3;
LN_ALL(97:480)=4;
LN_ALL(481:1120)=5;
LN_ALL(1121:2144)=6;
LN_ALL(2145:3168)=7;
LN_ALL(3169:4192)=8;
LN_ALL(4193:5216)=9;
LN_ALL(5217:6240)=10;
LN_ALL(6241:7264)=11;
LN_ALL(7265:8288)=12;
LN_RANGE=zeros(8288,1);
for K=1:length(LN_SELECT),
    LN_RANGE=LN_RANGE+ (LN_ALL==LN_SELECT(K));
end
LN_RANGE=find(LN_RANGE>0);



if ispc()
    % addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\libsvm');
    addpath(genpath('C:\Users\user\Dropbox (PPCA)\Research MIT\toolboxes'));
elseif ismac()
    addpath(genpath('~/ResearchMIT/toolboxes/'))
    
end

NTRAIN=round(size(FEATURESab{1},1)*9/10);

vals=[];labels=[];
for I=1:Nab
    vals=[vals;FEATURESab{I}(1:NTRAIN,LN_RANGE)];
    labels=[labels;ones(size(FEATURESab{I}(1:NTRAIN,LN_RANGE),1),1)*I];
    
end

info=cell(length(labels),1);cnt=0;
for I=1:Nab
    for J=1:NTRAIN
        cnt=cnt+1;
        info{cnt,1}.fname=INFOab{I}{J}.fname;
        info{cnt,1}.range=INFOab{I}{J}.range;
        %         info{cnt,1}.Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        %         info{cnt,1}.audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        if isfield(INFOab{I}{J,1},'fs')
            fs=INFOab{I}{J,1}.fs;
        else
            fs=16000;
        end
        
        info{cnt,1}.audio=INFOab{I}{J}.audio;
        info{cnt,1}.fs=INFOab{I}{J}.fs;
        
        if isfield(INFOab{I}{J},'jstats')
            info{cnt,1}.jstats=INFOab{I}{J}.jstats;
            vecformat=INFOab{I}{J}.jstats.vecformat;
            
            
            
            xleg=INFOab{I}{J}.jstats.xleg;
            yleg=INFOab{I}{J}.jstats.yleg;
            xlabels=INFOab{I}{J}.jstats.xlabels;
            ylabels=INFOab{I}{J}.jstats.ylabels;
            titles=INFOab{I}{J}.jstats.titles;
            
            ovecformat=vecformat;
            oxlabels=xlabels;
            oylabels=ylabels;
            otitles=titles;
            vecformat=cell(length(LN_SELECT),1);
            xlabels=cell(length(LN_SELECT),1);
            ylabels=cell(length(LN_SELECT),1);
            for K=1:length(LN_SELECT),
                vecformat{K}=ovecformat{LN_SELECT(K)};
                xlabels{K}=oxlabels{LN_SELECT(K)};
                ylabels{K}=oylabels{LN_SELECT(K)};
                titles{K}=otitles{LN_SELECT(K)};
                
            end
        end
        
        
        if isfield(INFOab{I}{J},'xlgnd')
            xlgnd=INFOab{I}{J}.xlgnd;
            xlgnd_name=INFOab{I}{J}.xlgnd_name;
            
        else
            xlgnd=Hz_mod_cfreqs;
            xlgnd_name='mod fq(Hz)';
        end
        
        if isfield(INFOab{I}{J},'ylgnd')
            ylgnd=INFOab{I}{J}.ylgnd;
            ylgnd_name=INFOab{I}{J}.ylgnd_name;
            
        else
            ylgnd=audio_cutoffs_Hz;
            ylgnd_name='fq(Hz)';
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
    valsT=[valsT;FEATURESab{I}((NTRAIN+1):end,LN_RANGE)];
    labelsT=[labelsT;ones(size(FEATURESab{I}((NTRAIN+1):end,LN_RANGE),1),1)*I];
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
        %         infoT{cnt,1}.Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        %         infoT{cnt,1}.audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        
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
%%
minerr = min(min(err));
NN=size(vals,2);
FEATNUM=sort([FEATNUM,NN*0.1,NN*0.2,NN*0.3,NN*0.4,NN*0.5,NN*0.7,NN*0.9,NN*1.1,inf]);
ps=nan(size(FEATNUM));
qs=nan(size(FEATNUM));

for K=1:length(FEATNUM),
    NUM_COEF=FEATNUM(K);
    low_limit = min(min(err(numpred <= NUM_COEF)));
    [mp,mq] = find(err == low_limit,1);
    ps(K)=mp;
    qs(K)=mq;
end

temp=unique([ps;qs]','rows');temp=temp(~isnan(temp(:,1)),:);
ps=temp(:,1);
qs=temp(:,2);

%%
accTs=nan(size(length(ps)));
accs=nan(size(length(ps)));

for K=1:length(ps),
    p=ps(K);
    q=qs(K);
    %[p,q] = find(err == minerr,1);
    mynum=numpred(p,q);
    myperr=err(p,q);
    mygamma= gamma(p);
    mydelta=delta(p,q);
    obj = fitcdiscr(vals,labels,'SaveMemory','on','FillCoeffs','on','Gamma',mygamma,'delta',mydelta);
    %fprintf('best gamma %3.3g best delta %3.3g err %3.3g\n',mygamma,mydelta,minerr);
    
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
    accs(K)=acc;
    accTs(K)=accT;
    toc
    fprintf('RESULTS acc on train %g acc on test %g my_num %g myperr %g mygamma %g mydelta %g| %s |%s \n',acc,accT,mynum,myperr,mygamma,mydelta,sprintf('%d ',LN_SELECT),fname);
end

%%
if DO_PLOT
    figure (30);clf
    indx = repmat(1:size(delta,2),size(delta,1),1);
    subplot(2,2,1)
    imagesc(err);hold on;
    for K=1:length(ps),
        p=ps(K);
        q=qs(K);
        plot(p,q,'ro','MarkerFaceColor','k');
    end
    
    
    colorbar;
    colormap('jet')
    title 'Classification error';
    xlabel 'Delta index';
    ylabel 'Gamma index';
    
    subplot(2,2,2)
    imagesc(log2(numpred));hold on;
    for K=1:length(ps),
        p=ps(K);
        q=qs(K);
        plot(p,q,'ro','MarkerFaceColor','k');
    end
    
    colorbar;
    title 'log2 of Number of predictors in the model';
    xlabel 'Delta index' ;
    ylabel 'Gamma index' ;
    
    subplot(2,2,3);
    plot(err,numpred,'k.');hold on;
    for K=1:length(ps),
        p=ps(K);
        q=qs(K);
        
        plot(err(p,q),numpred(p,q),'ro');
    end
    
    xlabel('Error rate');
    ylabel('Number of predictors');
    subplot(2,2,4);
    b=bar([accs;accTs]'*100);b(1).FaceColor='blue';
    ylim([50 ,100]);
    legend('train','test');
    set(gca,'XtickLabel',diag(numpred(ps,qs)));
end


% % %%
% % toc
% %
% % fprintf('apply predict on train set\n');
% % tic
% % [prdc,score] = predict(obj,vals);
% % linCoeffs=obj.Coeffs(2,1).Linear;
% % acc=sum(prdc==labels)/length(labels);
% % toc
% %
% % disp(sprintf('predict on test set...'));
% % tic
% % [prdcT,scoreT] = predict(obj,valsT);
% % accT=sum(prdcT==labelsT)/length(labelsT);
% %
% % % [ypredT]=plsval(PLS,XpT,ypT,Ap);
% % % accT=sum(sign(ypredT)==ypT)/length(ypT);
% % toc
%%
% disp(sprintf('displaying data projected to sub spaces...'));
% figure(50);clf;
% NTS=5;
% for I=1:NTS,
%     for J=(I):NTS,
%         subplot(NTS,NTS,(I-1)*NTS+J);
%         if I==J
%             plot(PLS.X_scores(yp>0,J)+randn(size(PLS.X_scores(yp>0,J))),PLS.X_scores(yp>0,I)+randn(size(PLS.X_scores(yp>0,I))),'.b');hold on;
%             plot(PLS.X_scores(yp<0,J)+randn(size(PLS.X_scores(yp<0,J))),PLS.X_scores(yp<0,I)+randn(size(PLS.X_scores(yp<0,I))),'.r');
%         else
%             plot(PLS.X_scores(yp>0,J),PLS.X_scores(yp>0,I),'.b');hold on;
%             plot(PLS.X_scores(yp<0,J),PLS.X_scores(yp<0,I),'.r');
%         end
%         xlabel(sprintf('%d',J));
%         ylabel(sprintf('%d',I));
%     end
% end
% figure(51);clf;
%
% hold on;I1=1;I2=2;
% pos=(labels~=1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.r');
% pos=(labels==1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.b');
%
% [val,idx]=sort(PLS.X_scores(:,I1));

%%

fprintf('showing filters...\n');

if DO_PLOT
    cmp=zeros(64,3);
    for I=1:32,
        cmp(32+I,:)=[0,0,I/32];
        cmp(I,:)=[1-I/32,0,0];
        
    end
    for I=1:min(15,Ap),
        figure(170+I);
        %     subplot(4,4,I);
        wf=linCoeffs;
        %wf=nan(size(wf));wf(1:1000)=A(1:1000);wf(isnan(wf))=0;
        
        %     xlgnd=audio_cutoffs_Hz;
        %     ylgnd=Hz_mod_cfreqs;
        %     wfr=reshape(wf,[length(ylgnd),length(xlgnd)]);
        
        %     imagesc(xlgnd,ylgnd,wfr);axis xy;
        
        %     nori_log_imagesc(xlgnd,ylgnd,wfr,[])
        data=nori_cell_array_unvectorize(wf,vecformat);is_colormap_negative=true;
        clim=[-max(abs(wf)),max(abs(wf))];
        %clim=[-1,1];
        %clim=[];
        %clim=[-5*sqrt(mean(wf.*wf)),5*sqrt(mean(wf.*wf))];
        nori_figure_stat_summary_as_cell_array(data,xleg,yleg,xlabels,ylabels,titles,is_colormap_negative,clim)
        
        
        
        colormap(cmp)
        xlabel(xlgnd_name);
        ylabel(ylgnd_name);
        %title(sprintf('%d',I));
        
    end
end
%%
if DO_PLOT
    figure(172);
    dp=sqrt(2)*(mean(vals(labels==1,:))-mean(vals(labels~=1,:)))./(std(mean(vals(labels==1,:))) + std(mean(vals(labels~=1,:))));
    dataDP=nori_cell_array_unvectorize(dp,vecformat);is_colormap_negative=true;
    %     clim=[-max(abs(wf)),max(abs(wf))];
    %     clim=[-3,3];
    % clim=[];
    clim=[-3*sqrt(mean(dp.*dp)),3*sqrt(mean(dp.*dp))];
    nori_figure_stat_summary_as_cell_array(dataDP,xleg,yleg,xlabels,ylabels,titles,is_colormap_negative,clim)
    colormap(cmp)
    %xlabel(xlgnd_name);
    %ylabel(ylgnd_name);
end

%%
if DO_PLOT
    dp1= (mean(vals(labels==1,:)));
    dp2= (mean(vals(labels==2,:)));
    
    %     clim=[-max(abs(wf)),max(abs(wf))];
    %     clim=[-3,3];
    %clim=[];
    clim=[-3*sqrt(mean(dp1.*dp1)),3*sqrt(mean(dp1.*dp1))];
    figure(173);dataDP1=nori_cell_array_unvectorize(dp1,vecformat);is_colormap_negative=true;
    
    nori_figure_stat_summary_as_cell_array(dataDP1,xleg,yleg,xlabels,ylabels,titles,is_colormap_negative,clim);
    
    figure(174);dataDP2=nori_cell_array_unvectorize(dp2,vecformat);is_colormap_negative=true;
    nori_figure_stat_summary_as_cell_array(dataDP2,xleg,yleg,xlabels,ylabels,titles,is_colormap_negative,clim);
    
    
    figure(175);dataDP3=nori_cell_array_unvectorize(dp1-dp2,vecformat);is_colormap_negative=true;
    nori_figure_stat_summary_as_cell_array(dataDP3,xleg,yleg,xlabels,ylabels,titles,is_colormap_negative,clim);
    
    %colormap(cmp)
    %xlabel(xlgnd_name);
    %ylabel(ylgnd_name);
end

%%
%title(sprintf('%d',I));


% figure(17);subplot(4,4,16);
%
% plot(acc*100,'kx-','LineWidth',2);hold on;
% plot([1 Ap],[accT,accT]*100,'g--','LineWidth',2);hold on;
% % plot([15 15],[0,100],'r--','LineWidth',1);hold on;
% axis([1 Ap 50 100]);
% title('Labeling accuracy');
% xlabel('Number of coefficents');
if DO_PLOT
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
end
%%

if DO_EXMP
    WHO=[1,2];
    for J=1:length(WHO),
        
        MAXSOUNDS=5;
        fprintf('displaying examples of stimuli\n');
        for I=1:MAXSOUNDS,
            figure(J*10000+I);
            %         subplot(floor(sqrt(MAXSOUNDS)+0.5),floor(sqrt(MAXSOUNDS)+0.5),I)
            idx=find(labels==WHO(J));
            
            %         nori_log_imagesc(xlgnd,ylgnd, (reshape(vals(idx(I),:),[length(ylgnd),length(xlgnd)]))',[]);axis xy;
            %         data=nori_cell_array_unvectorize(info{idx(idx(I))}.jstats.data,vecformat);is_colormap_negative=false;
            nori_figure_stat_summary_as_cell_array(info{idx(I)}.jstats.data,info{idx(I)}.jstats.xleg,info{idx(I)}.jstats.yleg,info{idx(I)}.jstats.xlabels,info{idx(I)}.jstats.ylabels,info{idx(I)}.jstats.titles,[],[]);
            
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
    WHO=[1:Ap]
    % WHO=[4]
    MAXSOUNDS=15;
    Nplay=10;
    NSHOW=3;
    for I1=WHO,
        % I1=6;
        % I2=1;
        % pos=(labels~=1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.r');
        % pos=(labels==1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.b');
        
        [val,idx]=sort(scores);
        Nplay=400;Njump=1;
        RANGES={[1:Njump:Nplay],[length(idx):-Njump:(length(idx)-Nplay)]};
        
        for k=1:length(RANGES),
            %     if k==1
            %         playme=cos(2*pi*(1:8000)/fs*300)';
            %     else
            %         playme=0.3*cos(2*pi*(1:8000)/fs*300)'+0.3*cos(2*pi*(1:8000)/fs*450)';
            %     end
            playme=zeros(fs,1);
            
            cnt=1;
            hashy=inf;
            MYHASHES=[];
            for J=RANGES{k},
                if cnt>MAXSOUNDS;
                    break
                end
                hashy=sum(double(info{idx(I)}.fname));
                if isempty(find(MYHASHES==hashy, 1));
                    %                 subplot(floor(sqrt(MAXSOUNDS)+0.5),floor(sqrt(MAXSOUNDS)+0.5),min(cnt,MAXSOUNDS));
                    cnt=cnt+1;NYHASHES=[MYHASHES,hashy];
                    %hashy= sum(vals(idx(J),:));
                    %nori_figure_stat_summary_as_cell_array(info{idx(idx(J))}.jstats.data,xleg,yleg,xlabels,ylabels,titles,[],[]);
                    if cnt<NSHOW
                        figure(1000+cnt+100*k);clf
                        %nori_figure_stat_summary_as_cell_array(info{idx(J)}.jstats.data,xleg,yleg,xlabels,ylabels,titles,[],[]);
                        nori_figure_stat_summary_as_cell_array(info{idx(I)}.jstats.data,info{idx(I)}.jstats.xleg,info{idx(I)}.jstats.yleg,info{idx(I)}.jstats.xlabels,info{idx(I)}.jstats.ylabels,info{idx(I)}.jstats.titles,[],[]);
                        
                    end
                    
                    %                 nori_log_imagesc(xlgnd,ylgnd,reshape(vals(idx(J),:),[length(ylgnd),length(xlgnd)]),[]);axis xy;
                    %                 axis xy;
                    %                 axis off;
                    myrng=info{idx(J)}.range;
                    %                 myrng=[round(max(info{idx(J)}.range)*1/3),round(max(info{idx(J)}.range)*2/3)];
                    
                    % do this if these are not saved
                    %                 [Y, fs]=audioread(info{idx(J)}.fname,myrng );
                    
                    Y=info{idx(J)}.audio;
                    fs=info{idx(J)}.fs;
                    
                    %                 HASH=sum(Y)+sum(round(Y*1000))/1000;
                    %                 if sum(find(MYHASHES==HASH))==0
                    %                     MYHASHES=[MYHASHES,HASH];
                    
                    % % %                 Yh=Y.*hanning(length(Y));
                    % % %                 n3=round(length(Yh)/3);
                    % % %                 playmeold=playme;
                    % % %                  playme=zeros(length(playmeold)+3*n3,1);
                    % % %                 playme(1:length(playmeold))=playmeold;
                    % % %                 playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))=playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))+Yh;
                    % % %
                    Yh=Y;
                    playmeold=playme;
                    playme=[playme;Yh];
                    %                 playme=zeros(length(playmeold)+3*n3,1);
                    %                 playme(1:length(playmeold))=playmeold;
                    %                 playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))=playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))+Yh;
                    %                 else
                    %                     fprintf('ignoring %g\n',HASH);
                    %
                    %                 end
                end
            end
            
            
            drawnow;
            p=audioplayer(playme,fs);p.play;
            ofname=sprintf('NLIN-10ktimi15k-concat-d-%d.direc-%d.wav',I1,k);
            
            audiowrite(ofname,playme,fs)
            pause(0.8*length(playme)/fs);
            pause
        end
    end
end
toc
fprintf('saving\n');
clear myFEATURES;clear FEATURES; clear FEATURESab; clear myINFO; clear info; clear INFO; clear INFOab;
save(moutfname);
fprintf('DONE!\n');
%%
