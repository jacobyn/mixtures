close all ;
W=100/1000*fs;
X=log(myCgrm{1}.^2)';
figure(1);
nori_log_imagesc(1:length(X),SubbandFrequencies,X,[]);
set(gca,'xticklabel',[]);

ttt=(1:W)/fs -50/1000;
figure;plot(ttt);


G=ttt.*exp((-(ttt).^2)/(2*(40/1000)^2));
ttt=ttt/sum(ttt);
ttt=ttt-mean(ttt);
plot(ttt,G);
figure(1);hold on;
for I=1:size(X,1),
    x=X(I,:);
    %plot(tt,x/max(x),'r');hold on;
    sx=conv(x,G,'same');
    sx=conv(sx,ones(1,100),'same');
    %subplot(3,1,2);
    %plot(tt,sx/max(sx),'xk');
    %hold on;
    [ymax2,imax2,ymin2,imin2] = extrema(sx/max(sx));
    pos1=ymax2>0.1;
    pos2=ymin2<-0.1;
    vec1=[(imax2(pos1))];
    vec2=[(imin2(pos2))];
    %plot(vec1,I,'r*');
    %plot(vec2,I,'g*');
    
end

%%


%%
figure;
tt=(1:size(X,2))/fs;
subplot(3,1,1);
x=X(50,:)
plot(tt,x/max(x),'r');hold on;
sx=conv(X(50,:),G,'same');
sx=conv(sx,ones(1,100),'same');
subplot(3,1,2);
plot(tt,sx/max(sx),'xk');
hold on;
[ymax2,imax2,ymin2,imin2] = extrema(sx/max(sx));
pos1=ymax2>0.1;
pos2=ymin2<-0.1;

plot(tt(imax2(pos1)),ymax2(pos1),'r*',tt(imin2(pos2)),ymin2(pos2),'g*')

%%

subplot(3,1,3);

%%
;