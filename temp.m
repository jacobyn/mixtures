%SVMModel = fitcsvm(vals(1:2:end,:),labels(1:2:end),'Standardize',true,'KernelFunction','linear','BoxConstraint',1e1);
BoxL=[-2:9];
SigL=[-2:9];
%BoxL=[-5:9];
%SigL=[-1];

BoxS=10.^BoxL;
SigS=10.^SigL;
err=zeros(length(BoxS),length(SigS));



for I=1:length(BoxS),
       for J=1:length(SigS),
           sig=SigS(J);
           box=BoxS(I);
           SVMModel = fitcsvm(vals(1:2:end,:),labels(1:2:end),'Standardize',true,'KernelFunction','rbf','BoxConstraint',box,'KernelScale',sig);
           %SVMModel = fitcsvm(vals(1:2:end,:),labels(1:2:end),'Standardize',true,'KernelFunction','linear','BoxConstraint',box);
           
           [cl,sc]=predict(SVMModel,vals(2:2:end,:));
           er=sum(cl==labels(2:2:end,:))/length(cl)
           err(I,J)=er;
           figure(1);imagesc(BoxL,SigL,err',[0 1]);axis xy; colormap jet;colorbar;
           ylabel('sig');xlabel('box');
           
           drawnow;
       end
end
%%

SVMModel = fitcsvm(vals(1:2:end,:),labels(1:2:end),'Standardize',true,'KernelFunction','linear','BoxConstraint',1e1);
           [cl,sc]=predict(SVMModel,vals(2:2:end,:));
           er=sum(cl==labels(2:2:end,:))/length(cl)

%%
vals=vals(1:2:end,:);
labels=labels(1:2:end);
%%
SVMModel = fitcsvm(vals,labels,'Standardize',true,'KernelFunction','linear','BoxConstraint',1e1);

 valsM= (vals-repmat(SVMModel.Mu,size(vals,1),1))./repmat(SVMModel.Sigma,size(vals,1),1);
 valsTM=(valsT-repmat(SVMModel.Mu,size(valsT,1),1))./repmat(SVMModel.Sigma,size(valsT,1),1);
 y= (labels==1)*(1)+(labels==2)*(-1);
 idx=find(SVMModel.IsSupportVector);
 score=((y(idx).*SVMModel.Alpha)')*(valsM(idx,:)*(valsTM'))-SVMModel.Bias;
[lb,sc] = predict(SVMModel,vals);
ws=((y(idx).*SVMModel.Alpha)')*(valsM(idx,:)*(eye(size(vals,2))));
wfs=reshape(ws,[length(ylgnd),length(xlgnd)]);

 clim=[-max(max(abs(wfs))),max(max(abs(wfs)))];
   
    cmp=zeros(64,3);
for I=1:32,
    cmp(32+I,:)=[0,0,I/32];
    cmp(I,:)=[1-I/32,0,0];
    
end
figure(1);clf;nori_log_imagesc(xlgnd,ylgnd,wfs,clim,10)
colormap(cmp)
    xlabel(xlgnd_name);
    ylabel(ylgnd_name);
nori_log_imagesc(xlgnd,ylgnd,reshape(wfs,[length(ylgnd),length(xlgnd)]),[-max(abs(wfs)),max(abs(wfs))],[]);axis xy;
%%

% CVSVMModel = crossval(SVMModel);
% kfoldLoss(CVSVMModel)

sum(cl==labelsT)/length(labelsT)
% ws=((y(idx).*SVMModel.Alpha)')*(valsM(idx,:)*(eye(size(vals,2))));
% wfs=reshape(ws,[length(ylgnd),length(xlgnd)]);
% 
%  clim=[-max(max(abs(wfs))),max(max(abs(wfs)))];
%    
%     cmp=zeros(64,3);
% for I=1:32,
%     cmp(32+I,:)=[0,0,I/32];
%     cmp(I,:)=[1-I/32,0,0];
%     
% end
% figure(1);clf;nori_log_imagesc(xlgnd,ylgnd,wfs,clim,10)
% colormap(cmp)
%     xlabel(xlgnd_name);
%     ylabel(ylgnd_name);
%nori_log_imagesc(xlgnd,ylgnd,reshape(wfs,[length(ylgnd),length(xlgnd)]),[-max(abs(wfs)),max(abs(wfs))],[]);axis xy;

%%
NV=2000;
c = cvpartition(NV,'KFold',10);
minfn = @(z)kfoldLoss(fitcsvm(vals(1:NV,:),labels(1:NV),'CVPartition',c,...
    'KernelFunction','rbf','BoxConstraint',exp(z(2)),...
    'KernelScale',exp(z(1))));
opts = optimset('TolX',5e-4,'TolFun',5e-4);
[searchmin fval] = fminsearch(minfn,randn(2,1),opts);
z = exp(searchmin);
%%
SVMModel = fitcsvm(vals(1:NV,:),labels(1:NV),'KernelFunction','rbf',...
    'KernelScale',z(1),'BoxConstraint',z(2));
CVSVMModel = crossval(SVMModel);
kfoldLoss(CVSVMModel);
[lbT,scT] = predict(SVMModel,vals,'Standardize',true)