addpath(genpath('~/ResearchMIT/toolboxes/'))
%load('~/NJO/data/LIN-job-1-2-5.mat');
tic
fprintf('loading data...\n');
load('~/NJO/data/LIN-job-1-2-3-9.mat');
%load('~/NJO/data/LIN-job-1-2-3-4-5-6-7-8-9-10-11-12.mat');
toc

if (length(xleg)==12)
otitles=titles;
oxleg=xleg;
oyleg=yleg;
titles=cell(size(LN_SELECT));
xleg=cell(size(LN_SELECT));
yleg=cell(size(LN_SELECT));
for K=1:length(LN_SELECT),
    
    titles{K}=otitles{LN_SELECT(K)};
    xleg{K}=oxleg{LN_SELECT(K)};
    yleg{K}=oyleg{LN_SELECT(K)};
    
    
end
end
%%
NUM_COEF=4200;
low_limit = min(min(err(numpred <= NUM_COEF)));
[p,q] = find(err == low_limit,1);
mynum=numpred(p,q);
myperr=err(p,q);
mygamma= gamma(p);
mydelta=delta(p,q);
obj = fitcdiscr(vals,labels,'SaveMemory','on','FillCoeffs','on','Gamma',mygamma,'delta',mydelta);

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
%%
figure(171);clf;
wf=linCoeffs;
data=nori_cell_array_unvectorize(wf,vecformat);is_colormap_negative=true;
clim=[-max(abs(wf)),max(abs(wf))];
%clim=[-1,1];
clim=[];

%clim=[-5*sqrt(mean(wf.*wf)),5*sqrt(mean(wf.*wf))];
is_colormap_negative=true;

nori_figure_stat_summary_as_cell_array(data,xleg,yleg,xlabels,ylabels,titles,is_colormap_negative,clim)
cmp=nori_create_cmp();
colormap(cmp)
xlabel(xlgnd_name);
ylabel(ylgnd_name);


