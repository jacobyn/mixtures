clear all;close all;clc
addpath('/Users/jacoby/Dropbox (PPCA)/Research MIT/toolboxes/Sound_Texture_Synthesis_Toolbox');
%[y,fs]=audioread('Wind.wav');
[y,fs]=audioread('Insects.wav');
%[y,fs]=audioread('Bwater.wav');

y=0.03*y/sum(y.^2);

nori_doplay(y,fs);
S0 = nori_generate_modsubands(y, fs);



y1=randn(size(y));
y1=0.03*y1/sum(y1.^2);
ITER=30;
for k=1:ITER,
    S1 = nori_generate_modsubands(y1, fs);
    
    if mod(k,3)~=1
        A=eye(size(S1.all_mod_subbands_envs,1));
    else
        A=orth(randn(size(S1.all_mod_subbands_envs,1)));
    end
    A=eye(size(S1.all_mod_subbands_envs,1));%NOTENTOE
    
    new_all_mod_subbands_envs=nan(size(S1.all_mod_subbands_envs));
    n1=size(S0.all_mod_subbands_envs,1);
    n2=size(S0.all_mod_subbands_envs,2);
    n3=size(S0.all_mod_subbands_envs,3);
    
    S0_all_mod_subbands_envs=S0.all_mod_subbands_envs;
    S1_all_mod_subbands_envs=S1.all_mod_subbands_envs;
    
    S0_all_mod_subbands_envs=reshape(A*reshape(S0_all_mod_subbands_envs,[n1,n2*n3]),[n1,n2,n3]);
    S1_all_mod_subbands_envs=reshape(A*reshape(S1_all_mod_subbands_envs,[n1,n2*n3]),[n1,n2,n3]);
    
    if mod(k,3)~=2
        B=eye(size(S1.all_mod_subbands_envs,3));
    else
        B=orth(randn(size(S1.all_mod_subbands_envs,3)));
    end
    B=eye(size(S1.all_mod_subbands_envs,3)); %NOTENTOE
    
    S0_all_mod_subbands_envs=reshape(reshape(S0_all_mod_subbands_envs,[n1*n2,n3])*B,[n1,n2,n3]);
    S1_all_mod_subbands_envs=reshape(reshape(S1_all_mod_subbands_envs,[n1*n2,n3])*B,[n1,n2,n3]);
    %
    
    
    for I=1:size(S1_all_mod_subbands_envs,1),
        for J=1:size(S1_all_mod_subbands_envs,3),
            v1=S1_all_mod_subbands_envs(I,:,J);
            v2=S0_all_mod_subbands_envs(I,:,J);
            v2n=transform_vals_according_to_hist(v1,v2);
            new_all_mod_subbands_envs(I,:,J)=v2n;
        end
    end
    new_all_mod_subbands_envs=reshape(reshape(new_all_mod_subbands_envs,[n1*n2,n3])*(B'),[n1*n2,n3]);
    new_all_mod_subbands_envs=reshape((A')*reshape(new_all_mod_subbands_envs,[n1,n2*n3]),[n1,n2,n3]);
    figure(4);clf;nhist({v1,v2,v2n});title(sprintf('iteration %d',k));drawnow
    
    
    %new_mod_subbands=(S1.all_mod_subbands./S1.all_mod_subbands_envs).*S0.all_mod_subbands_envs;
    
    new_mod_subbands=(S1.all_mod_subbands./S1.all_mod_subbands_envs).* new_all_mod_subbands_envs;
    
    
    [collapse_subband_envs_n_up,collapse_subband_envs_n]=nori_collapse_modsubands_to_subbands(new_mod_subbands,S1.mod_filts,S1.P.audio_sr,S1.P.env_sr);
    
    %     for J=1:S1.all_mod_subbands_envs
    
    orig_subband_envs_1_up=resample(S1.subband_envs,S1.P.audio_sr,S1.P.env_sr);
    collapse_subbands_1=(S1.subbands./(orig_subband_envs_1_up+eps)).*collapse_subband_envs_n_up;
    collapse_sound_1=collapse_subbands(collapse_subbands_1,S1.audio_filts);
    z=collapse_sound_1;
    y1=z;
    
    
    
    
    S3 = nori_generate_modsubands(z, fs);
    S0s = nori_measure_texture_stats_simplify(y, fs);
    S3s = nori_measure_texture_stats_simplify(z, fs);
    
    K=3;J=12;
    figure(1);clf;
    plot((S0.all_mod_subbands_envs(K,:,J)));hold all;
    plot(S1.all_mod_subbands_envs(K,:,J));hold all;
    plot(S3.all_mod_subbands_envs(K,:,J));hold all;
    title(sprintf('iteration %d',k));
    
    figure(2);
    subplot(2,2,1);
    imagesc(S0s.mod_power');axis xy;
    title(sprintf('iteration %d',k));
    subplot(2,2,2);
    imagesc(S3s.mod_power');axis xy;
    title(sprintf('iteration %d',k));
    subplot(2,2,3);
    imagesc(S0s.subband_envs');axis xy;
    title(sprintf('iteration %d',k));
    subplot(2,2,4);
    imagesc(S3s.subband_envs');axis xy;
    title(sprintf('iteration %d',k));
    
    figure(3);clf;
    plot(S0s.env_mean);hold all;
    plot(S3s.env_mean);
    title(sprintf('iteration %d',k));
    drawnow;
    mysound=z;
    p2 = audioplayer(mysound/max(mysound), fs);p2.play
    %pause(1.1*length(mysound)/fs +0.5);
end
%%


mysound=z;
p2 = audioplayer(mysound/max(mysound), fs);p2.play
pause(1.1*length(mysound)/fs +0.1);

mysound=y;
p2 = audioplayer(mysound/max(mysound), fs);p2.play
pause(1.1*length(mysound)/fs +0.1);

mysound=z;
p2 = audioplayer(mysound/max(mysound), fs);p2.play



