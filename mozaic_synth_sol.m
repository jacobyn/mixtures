function playme=mozaic_synth_sol(sol,info)


playme=[];


for III=1:length(sol),
    IDX=sol(III);
    [Y, fs]=audioread(info{IDX}.fname, info{IDX}.range);
    if isempty(playme)
        playme=Y;
    else
    Yh=Y.*hanning(length(Y));
    n33=round(length(Yh)/3);
    playmeold=playme;
    playme=zeros(length(playmeold)+n33,1);
    playme(1:length(playmeold))=playmeold;
    playme((length(playmeold)-(2*n33)+1):(length(playmeold)-(2*n33)+length(Yh)))=playme((length(playmeold)-(2*n33)+1):(length(playmeold)-(2*n33)+length(Yh)))+Yh;
    end
    
end
