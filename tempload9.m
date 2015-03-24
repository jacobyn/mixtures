%[y,fs]=wavread('norm_SE2-67 Insects In A Swamp.wav');
[y,fs]=wavread('norm_SSE-41 Lapping Waves.wav');

 S1 = nori_measure_texture_stats(y, fs); %measure stats
%S = measure_texture_stats(sample_sound, P, varargin)

        S1=S1.S;
        figure(20);clf
        nori_log_imagesc(xlgnd,ylgnd,reshape(S1.mod_power,[length(ylgnd),length(xlgnd)]),[],[]);title('waves');colormap jet
%%

        [y,fs]=wavread('pink.wav');

 S2 = nori_measure_texture_stats(y, fs); %measure stats
%S = measure_texture_stats(sample_sound, P, varargin)

        S2=S2.S;
        figure(21);clf
        nori_log_imagesc(xlgnd,ylgnd,reshape(S2.mod_power,[length(ylgnd),length(xlgnd)]),[],[]);title('pink');colormap jet
        %%
        
        figure(23);
        nori_log_imagesc(xlgnd,ylgnd,reshape(10*log10(S1.mod_power)-10*log10(S2.mod_power),[length(ylgnd),length(xlgnd)]),[],[]);title('pink');colormap jet
        