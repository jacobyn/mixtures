clear all;close all;clc;
isWIN=false;
DO_CORPUS=true;
ITER=1000; % number of random
Ms=[1:9]; %all possible number of speakers

if isWIN
    base_corpus='M:\shared\Sounds\Speech\timit\train';
    output_data='N:\mixtures\data';
    fname_corpus='timit_corpus_WIN';
else
    base_corpus='~/timit/train';
    output_data='~/mixtures/data';
    fname_corpus='timit_corpus_UNIX';
    
end

length_min=1;
FS0=16000;

if DO_CORPUS
    disp('reading corpus...');
    cd (base_corpus);
    dialects=dir ('dr*');
    ND=length(dialects);
    DNAMES=cell(ND,1);
    NG=2;
    GNAMES={'f','m'};
    for I=1:ND,
        DNAMES{I}=dialects(I).name;
    end
    NS=14;
    SPEAKERS=cell(ND,NG,NS);
    NSPEAKERS=nan(ND,NG);
    
    NW=10;
    WAVS=cell(ND,NG,NS,NW);
    NWAVS=nan(ND,NG,NS);
    FILES=cell(ND,NG,NS,NW);
    for ID=1:ND,
        for IG=1:NG,
            if isWIN
                mdir=sprintf('%s\\%s',base_corpus,DNAMES{I});
            else
                mdir=sprintf('%s/%s',base_corpus,DNAMES{I});
            end
            
            cd (mdir);
            spks=dir(sprintf('%s*0',GNAMES{IG}));
            NSPEAKERS(ID,IG)=length(spks);
            display(sprintf('reading dialect %2d of %2d (%s=%d speakers)...',ID,ND,GNAMES{IG},NSPEAKERS(ID,IG)));
            for IS=1:NSPEAKERS(ID,IG)
                SPEAKERS{ID,IG,IS}=spks(IS).name;
                if isWIN
                    mdir=sprintf('%s\\%s\\%s',base_corpus,DNAMES{I},SPEAKERS{ID,IG,IS});
                else
                    mdir=sprintf('%s/%s/%s',base_corpus,DNAMES{I},SPEAKERS{ID,IG,IS});
                end
                
                
                cd (mdir);
                mwavs=dir('*.wav');
                NWAVS(ID,IG,IS)=length(mwavs);
                
                display(sprintf('\t speaker %s %d wavs',SPEAKERS{ID,IG,IS},NWAVS(ID,IG,IS)));
                for IW=1:NWAVS(ID,IG,IS),
                    WAVS{ID,IG,IS,IW}=mwavs(IW).name;
                    %             display(sprintf('%s',WAVS{ID,IG,IS}));
                    if isWIN
                        FILES{ID,IG,IS,IW}=sprintf('%s\\%s\\%s\\%s',base_corpus,DNAMES{I},SPEAKERS{ID,IG,IS},WAVS{ID,IG,IS,IW});
                    else
                        FILES{ID,IG,IS,IW}=sprintf('%s/%s/%s/%s',base_corpus,DNAMES{I},SPEAKERS{ID,IG,IS},WAVS{ID,IG,IS,IW});
                    end
                    
                    display(sprintf('\t\t%s',FILES{ID,IG,IS,IW}));
                    
                end
            end
        end
    end
    TIMIT.ND=ND;
    TIMIT.NG=NG;
    TIMIT.NS=NS;
    TIMIT.NW=NW;
    TIMIT.DNAMES=DNAMES;
    TIMIT.GNAMES=GNAMES;
    TIMIT.SPEAKERS=SPEAKERS;
    TIMIT.NSPEAKERS=NSPEAKERS;
    TIMIT.WAVS=WAVS;
    TIMIT.NWAVS=NWAVS;
    TIMIT.FILES=FILES;
    cd (output_data);
    save(fname_corpus,'TIMIT');
    
else
    cd (output_data);
    load(fname_corpus);
    
end

%     rng('shuffle');


for mm=1:length(Ms),
    M=Ms(mm);
    
    cd (output_data);
    for KK=1:ITER,
        display(sprintf('creating file %d of %d\t M=%d',KK,ITER,M));
        mixme=cell(M,1);
        mynames=cell(M,1);
        
        mlength=length_min;
        for I=1:M,
            iD=randi(TIMIT.ND,1,1);
            iG=randi(TIMIT.NG,1,1);
            iS=randi(TIMIT.NSPEAKERS(iD,iG),1,1);
            iW=randi(TIMIT.NWAVS(iD,iG,iS),1,1);
            if isWIN
                [R,FS]=audioread(TIMIT.FILES{iD,iG,iS,iW});
            else
                [R,FS]=wavread(TIMIT.FILES{iD,iG,iS,iW});
            end
            assert(FS==FS0);
            %     AP=audioplayer(R1,FS);
            %     AP.play
            lengthR=size(R,1)/FS;
            
            
            mixme{I}=0.03*(R/sqrt(mean((R.^2))));
            mlength=min(mlength,lengthR);
            mynames{I}=TIMIT.WAVS{iD,iG,iS,iW};
        end
        
        soundM=zeros(round(FS*mlength),1);
        for I=1:M,
            soundM(1:round(FS*mlength))=soundM(1:round(FS*mlength))+ mixme{I}(1:round(FS*mlength));
        end
        soundM=0.03*(soundM/sqrt(mean((soundM.^2))));
        mfname=sprintf('mix-v1.M.%d.%d.wav',M,KK);
        if isWIN
            audiowrite(mfname,soundM,FS0)
        else
            wavwrite(soundM,FS0,mfname);
        end
    end
end
% AP=audioplayer(soundM,FS);
% AP.play






