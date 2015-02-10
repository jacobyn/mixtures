Mdl = fitcdiscr(vals,labels,'SaveMemory','on','FillCoeffs','off');
rng(8000,'twister') % For reproducibility
[err,gamma,delta,numpred] = cvshrink(Mdl,'NumGamma',29,'NumDelta',29,'Verbose',1);
minerr = min(min(err));
[p,q] = find(err == minerr);
mygamma= gamma(p);
mydelta=delta(p,q);
fprintf('best gamma %3.3g best delta=%3.3g err=%g',mygamma,mydelta,minerr);

disp(sprintf('predict on train set...'));
tic
[prdc,score] = predict(obj,valsN);
linCoeffs=obj.Coeffs(2,1).Linear;
acc=sum(prdc==labels)/length(labels);
toc
obj = fitcdiscr(vals,labels,'SaveMemory','on','FillCoeffs','on','Gamma',mygamma,'delta',mydelta);
disp(sprintf('predict on test set...'));
tic
[prdcT,scoreT] = predict(obj,valsNT);
accT=sum(prdcT==labelsT)/length(labelsT);
toc

%%
figure;
plot(err,numpred,'k.')
xlabel('Error rate');
ylabel('Number of predictors');
%%

low250 = min(min(err(numpred <= 250)));
lownum = min(min(numpred(err == low250)));
[low250 lownum]

% 
low1000 = min(min(err(numpred <= 1000)));
lownum = min(min(numpred(err == low1000)));
[low1000 lownum]
[r,s] = find(err == low1000,1)

Mdl.Gamma = gamma(r);
Mdl.Delta = delta(r,s);

%
obj = fitcdiscr(vals,labels,'SaveMemory','on','FillCoeffs','on','Gamma',mygamma,'delta',mydelta);
%%

[prdc,score] = predict(obj,valsN);
linCoeffs=obj.Coeffs(2,1).Linear;
acc=sum(prdc==labels)/length(labels);
toc

disp(sprintf('predict on test set...'));
tic
[prdcT,scoreT] = predict(obj,valsNT);
accT=sum(prdcT==labelsT)/length(labelsT);

%%
% Create the Delta index matrix
indx = repmat(1:size(delta,2),size(delta,1),1);
figure
subplot(1,2,1)
imagesc(err);
colorbar;
colormap('jet')
title 'Classification error';
xlabel 'Delta index';
ylabel 'Gamma index';

subplot(1,2,2)
imagesc(numpred);
colorbar;
title 'Number of predictors in the model';
xlabel 'Delta index' ;
ylabel 'Gamma index' ;