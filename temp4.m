%function S = nori_measure_texture_stats_simplify(orig_sound, fs)

% Crutial parameters
P.N_audio_channels = 30; %this is the number excluding lowpass and highpass filters on ends of spectrum
P.low_audio_f = 20; %Hz
P.hi_audio_f = 10000; %Hz
P.use_more_audio_filters = 0; % should be 1 if 2x overcomplete, 2 if 4x overcomplete
%P.lin_or_log_filters = 1; %1--> log acoustic & mod; 2--> log acoust, lin mod; 3--> lin acoust, log mod; 4--> lin acoust, lin mod

P.env_sr = 400;
P.N_mod_channels = 20; %These next four parameters control the modulation filterbank from which modulation power is measured
% P.low_mod_f = 0.5; %Hz
P.low_mod_f = 1; %Hz %%%%%%%%%%NOTE NOTE NOTE!!! I HAD TO CHANGE IT BECAUSE I WORK WITH SHORT SAMPLES
P.hi_mod_f = 200; %Hz


% this make sure we heve even number of samples (change that!)
my_length=floor(length(orig_sound)/2)*2; % Library equires a muliple of 2!! (correct that inside the library, this should not be like that)
orig_sound=orig_sound(1:my_length);


P.audio_sr = fs;


P.measurement_windowing = 2; %1 --> unwindowed (circular boundary handling), 2 --> global window
P.imposition_windowing = 1; %1 --> unwindowed (circular boundary handling), 2 --> global window
P.win_steepness = .5; % must be between 0 and 1; smaller means steeper window edges

P.use_more_mod_filters=0; % should be 1 if 2x overcomplete,
P.mod_filt_Q_value = 2;

P.use_zp = 0;% 0 means circular convolution; 1 means zeropadding (for modulation filtering)
P.low_mod_f_C12=1; %Hz - this is the lowest frequency in the octave-spaced modulation filterbank used for the C1 and C2 correlations

% always do for power compression
P.comp_exponent = .3;
P.log_constant = 10^-12; %this is the constant added prior to taking the log

P.n_hist_bins = 128;

P.env_ac_intervals_smp = [1 2 3 4 5 6 7 9 11 14 18 22 28 36 45 57 73 92 116 148 187 237 301]; %in samples

% % %for subband autocorrelation measurement
P.sub_ac_undo_win = 1; % divide by ac of window
P.sub_ac_win_choice = 2; % type of window
P.num_sub_ac_period = 5; %num periods of subband cf over which to match sub ac


%make filters
[audio_filts, audio_cutoffs_Hz] = make_erb_cos_filters(length(orig_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);%%

%generate subbands and subband envelopes from which statistics are measured
subbands = generate_subbands(orig_sound, audio_filts);
%analytic_subbands = hilbert(subbands);
subband_envs = abs(hilbert(subbands));
%power compression
subband_envs = subband_envs.^P.comp_exponent;%
%downsampling
ds_factor=P.audio_sr/P.env_sr;
subband_envs = resample(subband_envs,1,ds_factor);
subband_envs(subband_envs<0)=0;
mod_filt_length = length(subband_envs);%%
measurement_win = ones(length(subband_envs),1);%


%for comparison: [audio_filts, audio_cutoffs_Hz] = make_erb_cos_filters(length(sample_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
% [mod_filts,Hz_mod_cfreqs,mod_freqs] = make_constQ_cos_filters(mod_filt_length, P.env_sr, P.N_mod_channels, P.low_mod_f, P.hi_mod_f, P.mod_filt_Q_value);%%
% make_lin_cos_filters(length(sample_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
[mod_filts,Hz_mod_cfreqs,mod_freqs] = make_lin_cos_filters(mod_filt_length, P.env_sr, P.N_mod_channels, P.low_mod_f, P.hi_mod_f);
all_mod_subbands=zeros(size(audio_filts,2),size(subband_envs,1),size(mod_filts,2));
all_mod_subbands_envs=zeros(size(audio_filts,2),size(subband_envs,1),size(mod_filts,2));


for j=1:size(audio_filts,2) %go through subbands
   mys=subband_envs(:,j);
    %if P.use_zp    mod_subbands = generate_subbands_zp(s, mod_filts);
   my_mod_subbands = generate_subbands(mys, mod_filts);
   all_mod_subbands(j,:,:)=my_mod_subbands;
   all_mod_subbands_envs(j,:,:)=abs(hilbert(my_mod_subbands));
end

S.all_mod_subbands_envs=all_mod_subbands_envs;
S.all_mod_subbands=all_mod_subbands;
S.mod_filts=mod_filts;
S.Hz_mod_cfreqs=Hz_mod_cfreqs;
S.mod_freqs=mod_freqs;
S.subbands=subbands;
S.subband_envs=subband_envs;
S.P=P;

% 
% %%%%
% %do collaps
% %%%
% 
% collapse_subband_envs=zeros(size(subband_envs));
% for j=1:size(audio_filts,2) %go through subbands
%     my_mod_subbands=reshape(all_mod_subbands(j,:,:),size(all_mod_subbands,2),size(all_mod_subbands,3));
%     my_colapse_envs=collapse_subbands(my_mod_subbands, mod_filts);
%     collapse_subband_envs(:,j)=my_colapse_envs;
% end



