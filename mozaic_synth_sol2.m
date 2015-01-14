function playme=mozaic_synth_sol2(sol,info)


playme=[];


for III=1:size(sol,1),
    IDX1=sol(III,1);
    IDX2=sol(III,2);
    [Y1, ~]=audioread(info{IDX1}.fname, info{IDX1}.range);
    [Y2, ~]=audioread(info{IDX2}.fname, info{IDX2}.range);
    
    Y=Y1+Y2;
    
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
