function MYTHRESH=compute_rms(files,VOLUMEITER,percent_remove,DUR)
if isempty(percent_remove)
    percent_remove=.18;
end
NF=length(files);
rmss=nan(VOLUMEITER,1);

for KK=1:VOLUMEITER,
    fprintf('running volumeiter %d of %d\n',KK,VOLUMEITER);
    tsmpls=0;tfs=0;
    while (tsmpls-tfs*DUR)<=0
        tiD=randi(NF,1,1);
        
        tfname=files(tiD).name;
        tinfo=audioinfo(tfname);
        tsmpls=tinfo.TotalSamples;
        tfs=tinfo.SampleRate;
        if (tsmpls-tfs*DUR)<0
            fprintf('filename %s too short (%g sec)\n', tfname,tsmpls/tfs)
        end
        
    end
    
    tmypos=randi(tsmpls-tfs*DUR,1,1);
    tmyrange=[tmypos, tmypos+tfs*DUR];
    
    [tY, ~]=audioread(tfname, tmyrange);
    if size(tY,2)==2
        tY=sum(tY,2);
    end
    rmss(KK)=sqrt(mean(tY.^2));
    
end
srms=sort(rmss);
MYTHRESH=srms(round(VOLUMEITER*percent_remove));