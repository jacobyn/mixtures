clear all;close all;clc
addpath('/Users/jacoby/Dropbox (PPCA)/Research MIT/toolboxes/Sound_Texture_Synthesis_Toolbox');
%[y,fs]=audioread('Wind.wav');
[y,fs]=audioread('Insects.wav');
y=0.03*y/sum(y.^2);

%nori_doplay(y,fs);
S0 = nori_generate_modsubands(y, fs);

y1=randn(size(y));
y1=0.03*y1/sum(y1.^2);
ITER=1;
S1 = nori_generate_modsubands(y1, fs);


figure(100);clf
cnt=1;

collapse_subband_envs=zeros(size(S0.all_mod_subbands,2),size(S0.all_mod_subbands,1));
for j=1:size(S0.all_mod_subbands,1) %go through subbands
    my_mod_subbands=reshape(S0.all_mod_subbands(j,:,:),size(S0.all_mod_subbands,2),size(S0.all_mod_subbands,3));
    my_colapse_envs=collapse_subbands(my_mod_subbands, S0.mod_filts);
    collapse_subband_envs(:,j)=my_colapse_envs;
end




for I=1:size(collapse_subband_envs,2)
    subplot(6,6,cnt);
    plot(S0.subband_envs(:,I),'b-','LineWidth',2);hold on;
    plot(collapse_subband_envs(:,I),'r-');
    cnt=cnt+1;
end


%%

collapse_subband_envs_up=resample(collapse_subband_envs,S0.P.audio_sr,S0.P.env_sr);


for I=1:size(collapse_subband_envs_n_up,2)
    subplot(6,6,cnt);
    plot(S0.subband_envs(:,I),'b-');hold on;
    plot(collapse_subband_envs_n(:,I),'r-');
    cnt=cnt+1;
end
%%
figure(100);clf;cnt=1;
for I=1:size(collapse_subband_envs_n_up,2)
    subplot(6,6,cnt);
    plot(S0.subband_envs_up(:,I).^(1/S0.P.comp_exponent),'b-','LineWidth',2);hold on;
    plot(S0.subbands(:,I),'r-');
    cnt=cnt+1;
end

