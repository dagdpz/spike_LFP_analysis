% Finding an optimal number of observation (pairs of spike-field phases)
% when using ppc1 for spike-field coherence analysis

clear all
% fr_vec = 0.25:0.25:10;
% fr_vec = 0.1:0.125:5;
% for frate = 1:numel(fr_vec),
% 
% clear data trialinfo keys morlet

N_cycles = 6; 

cfg = [];
cfg.numtrl	= 30;
cfg.fsample     = 1000; % Hz
cfg.trllen      = 4; % s

cfg.s1.freq     = 10; % Hz
cfg.s1.phase    = 'random'; % 'random' or 0
cfg.s1.ampl     = 1;

cfg.s2.freq = 25; % Hz
cfg.s2.ampl = 0.5;
cfg.s2.phase = 0;

cfg.s3.freq = 60; % Hz
cfg.s3.ampl = 0.2;
cfg.s3.phase = 0;

cfg.noise.ampl = 0.2;

cfg.method	= 'superimposed';
cfg.output	= 'all'; % mixed or all

data = ft_freqsimulation(cfg);
   
lfp_freq = [cfg.s1.freq cfg.s2.freq cfg.s3.freq];
lfp_amp =  [cfg.s1.ampl cfg.s2.ampl cfg.s3.ampl];

% spikes
% spikeRate = fr_vec(frate); % Hz, overall spike rate
spikeRate = [5]; % Hz, overall spike rate
spikePhaseFreq = 1; % 1 or 2 or 3, s1 or s2 or s3 components
spikePhaseMean = [-pi]; % align to trough
spikePhaseStd  = [pi/60]; % rad, spread around spikePhaseMean
spikeLockProb = 0.2; % probability of spike to be locked, the phase-locked firing rate would be max(spikeRate*spikeLockProb,lfp_freq(spikePhaseFreq))

lockedSpikes = zeros(cfg.numtrl,cfg.trllen*cfg.fsample);
Spikes = MS_test_FT_simulate_spike_train(spikeRate,cfg.trllen,cfg.numtrl); % non-locked spikes
for t = 1:cfg.numtrl,
	 data.trial{t} = [data.trial{t}; Spikes(t,:)]; 
end

disp(sprintf('spikeRate %d Hz, spikeLockProb %.2f, mean rate %.2f, mean locked rate %.2f %d spikes',spikeRate,spikeLockProb,mean(sum(Spikes,2)/cfg.trllen),mean(sum(lockedSpikes,2)/cfg.trllen),sum(sum(Spikes))));
data.label = {'lfp1', 's1', 's2', 's3', 'noise', 'spk1'};

%% Here i compute the number of spike pairs even before the convolution

spikecount = 0; 
for trr = 1:numel(data.trial),
    trialinfo{1}(spikecount+1:(spikecount + sum(data.trial{trr}(6, :))), 1) = trr;
    spikecount = spikecount + sum(data.trial{trr}(6, :)); 
end

n_p=0;
for temp_t=unique(trialinfo{1})'
    n_p=n_p+sum(trialinfo{1}~=temp_t)*sum(trialinfo{1}==temp_t)/2;
end
n_p

%% MORLET convolution settings
keys.spike_field               = [];
keys.spike_field.method        = 'ppc1'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
keys.spike_field.spikechannel  = 'spk1';
keys.spike_field.channel       = 'lfp1'; % selected LFP channels
keys.spike_field.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
keys.spike_field.timwin        = 'all'; % compute over all available spikes in the window
keys.spike_field.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
keys.spike_field.N_cycles      = 6;


keys.LFP.frequencies        =10.^(0.6:0.05:2);
keys.LFP.binsize            =1.1*keys.spike_field.N_cycles./keys.LFP.frequencies;
keys.LFP.Morlet_bandwidth   =2*(keys.spike_field.N_cycles./(2*pi*keys.LFP.frequencies)).^2;
keys.LFP.SR                 =1017.252604166667;

%morlet_borders=ceil(keys.LFP.binsize/2*keys.LFP.SR)/keys.LFP.SR;
for f=1:numel(keys.LFP.frequencies)
    fr=keys.LFP.frequencies(f);
    morlet_borders=ceil(keys.LFP.binsize(f)/2*keys.LFP.SR)/keys.LFP.SR;
    keys.LFP.Morlet_Matrix{f}=cmorwavf(-morlet_borders,morlet_borders,2*morlet_borders*keys.LFP.SR,keys.LFP.Morlet_bandwidth(f),fr);
end

% Morlet wavelets convlution
nspikesaccum=0;

for t=1:numel(data.trial)
    clear spectrum
    for f=1:numel(keys.LFP.Morlet_Matrix)
        filtered(f,:)=conv(data.trial{t}(1,:),keys.LFP.Morlet_Matrix{f},'same');
    end
    spectrum(:,:,1)=filtered(:,data.trial{t}(end,:)==1);
    nspikesend= nspikesaccum+size(spectrum,2);
    spectrum_perm(nspikesaccum+1:nspikesend,1,:)=permute(spectrum,[2,3,1]);
    timevec(nspikesaccum+1:nspikesend,1) = find(data.trial{t}(end,:)==1)/1000;
    trialvec(nspikesaccum+1:nspikesend,1) = t;
    trialtimevec(t, :) = [min(data.time{t}), max(data.time{t})] ;
    nspikesaccum=nspikesend;
end

morlet.lfplabel = {keys.spike_field.channel};
morlet.label = {keys.spike_field.spikechannel};
morlet.freq = keys.LFP.frequencies;
morlet.fourierspctrm = {spectrum_perm}; 
morlet.time = {timevec};
morlet.trial = {trialvec};
morlet.dimord = '{chan}_spike_lfpchan_freq';
morlet.trialtime = trialtimevec; 

%% Here allocate the data structures for ppc, FR, number of spike pairs
iter_count = 1;
spike_count = numel(morlet.trial{1});
step = 5; %remove 5 spikes at each iteration

while spike_count >= 15,
   tmp_ind = randi(spike_count,1,5);
   spike_count = spike_count - 5;
   
   morlet.time{1}(tmp_ind) = [];
   morlet.trial{1}(tmp_ind) = [];
   morlet.fourierspctrm{1}(tmp_ind, :, :) = [];

%% So then the number of independent pairs will be n*(n-1)/2 and the number of pairs for computing ppc considering the removal 
%  of same-trial spikes will be  n*(n-1)/2 - sum( nt*(nt-1)/2) where n is the number of spikes across all the trials, 
%  nt - nuber of spikes in each trial
    n_pairs=0;
    for temp_t=unique(morlet.trial{1})'
        n_pairs=n_pairs+sum(morlet.trial{1}~=temp_t)*sum(morlet.trial{1}==temp_t)/2;
    end

% % Well, here I doublecheck the number of spike pairs becase I am stupid and
% % cannot read the Lukas` code properly
% indep_across = numel(morlet.trial{1})*( numel(morlet.trial{1}) - 1) /2; %2087946
% excl_seld = 0; 
% for temp_t=unique(morlet.trial{1})'
%     excl_seld = excl_seld + (sum(morlet.trial{1}==temp_t) * (sum(morlet.trial{1}==temp_t)-1) /2 ) ;
% end
% out_num_sp = indep_across - excl_seld;


%% compute ppc
    cfg               = [];
    cfg.method        = keys.spike_field.method; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
    cfg.spikechannel  = morlet.label{1};
    cfg.channel       = morlet.lfplabel(1); % selected LFP channels
    cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
    cfg.timwin        = 'all'; % compute over all available spikes in the window
    cfg.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
    statSts           = ft_spiketriggeredspectrum_stat(cfg,morlet);

    data_iter.ppc(iter_count, :) = statSts.ppc1;
    data_iter.pairs(iter_count, 1) = n_pairs;
    data_iter.spikenum(iter_count, 1) = spike_count;
    data_iter.rate(iter_count, 1) = spike_count / 120; % 30 trials * 4 sec trial length
    
    iter_count = iter_count + 1;
end

% overall title - trial duration, 

% subplot(5,8,frate);
% %% plot the results
% % figure('Name','Spike-triggered spectrum');
% % hold on
% plot(statSts.freq,statSts.ppc1','r');
% % xlabel('frequency');
% % ylabel(cfg.method);
% xlim([4 100])
% % vline(statSts.freq(signif))
% set(gca, 'XScale', 'log')
% title(sprintf('FR %.2f Hz, %d spike pairs', spikeRate, n_p));
% % legend({'FFT','morlet convolution'});
% ylim([-0.03 0.03])
% text(5, 0.02, sprintf('%f\n%f', mean(statSts.ppc1), std(statSts.ppc1)))
% end

keys.path_to_save = pwd;

FigName_short2= sprintf('spike_pairs_plot7'); %name of the file
FigName_Long2=sprintf('Same data, remove 5 spikes at a step, trial duration 4 sec'); % overall title of the figure
set(0,'DefaultTextInterpreter','none');
wanted_papersize= [20 12];
h1=figure('outerposition',[0 0 1900 1200],'name', FigName_Long2);
figure(h1)

clrs = colormap(jet(numel(data_iter.rate)));
for ll = 1:numel(data_iter.rate),
    plot(morlet.freq, data_iter.ppc(ll, :), 'Color', clrs(ll, :)); hold on; 
end    
% legend(data_iter.pairs);
set(gca, 'XScale', 'log')
xlim([4 100])
ylim([min(min(data_iter.ppc)) max(max(data_iter.ppc))])
%val = round(linspace(1, numel(data_iter.pairs), 13));
val = linspace(10, 110, 11);
c = colorbar('YTickLabel', {data_iter.pairs(val)} );
title(sprintf('Same data, remove 5 spikes at a step, trial duration 4 sec'));
ylabel('ppc1')
xlabel('frequency, log')
set(gca,'XScale','log', 'XTick', [10 25 60]);

set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
mtit(figure(h1),  [FigName_Long2 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 14,'Interpreter', 'none');
export_fig(figure(h1), [keys.path_to_save filesep FigName_short2], '-pdf','-transparent') % pdf by run