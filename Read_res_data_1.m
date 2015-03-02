clear all;close all;clc
F=fopen('~/NJO/mixtures/RES3.RES','rt');
resnames={};
resdata={};
reslen=[];
while ~(feof(F))
    
    l=fgetl(F);
    mls=strsplit(strrep(l,'|',' '),' ');
    acc=str2double(mls{5});
    accT=str2double(mls{9});
    numP=str2double(mls{11});
    str=mls{18};
    mygamma=str2double(mls{15});
    mydelta=str2double(mls{17});
    fprintf('reading: %s\n',l);
    %[acc,accT,numP,mygamma,mydelta]
    
    isfound=false;
    K=-1;
    for K=1:size(resnames,1),
        if strcmp(resnames{K,1},str)
            isfound=true;
            break;
        end
    end
    if ~isfound
        pos=size(resnames,1)+1;
        resnames{pos,1}=str;
        resdata{pos,1}=[acc,accT,numP,mygamma,mydelta];
        reslen(pos)=length(str)+0.01*str2num(str(1));
        
    else
        pos=K;
        assert(strcmp(resnames{pos,1},str));
        resdata{pos}=[resdata{pos};acc,accT,numP,mygamma,mydelta];
        
    end
end
fclose (F);




figure(1);clf;
for I=1:length(resnames);
    for J=1:(size(resdata{I},1));
       placex=100*resdata{I}(J,2);
       placey=I+0+0.25*rand(1,1);
       plot(placex,placey,'.');hold all;
       text(placex+0.05,placey+0.05,sprintf('%d',resdata{I}(J,3)),'Fontsize',8);hold all;
    end
end
ylim([0.5 ,size(resnames,1)+0.5]);
set(gca,'Ytick',1:size(resnames,1));
set(gca,'YtickLabel',resnames);
