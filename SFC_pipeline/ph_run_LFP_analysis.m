function ph_run_LFP_analysis(varargin)
%ph_run_population_analysis('PPC_pulv_eye_hand',{'LIP_dPul_inj_working'},{'pop','ons','sct','ccs'})
% run from C:\Users\mslashcheva\Dropbox\DAG\DAG_toolbox\spike_analysis
% run ph_run_LFP_analysis('PPC_pulv_eye_hand',{'MIP_dPul_inj_working_post_sfn'})

population_analysis_to_perform={'pop','ons','sct','ccs'}; %,'pop'};

if nargin>2
    population_analysis_to_perform=varargin{3};
end
keys=struct;
project=varargin{1};
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);

if nargin>1
    keys.project_versions=varargin{2};
end

keys.spike_field               = [];
keys.spike_field.method        = 'ppc1'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
keys.spike_field.spikechannel  = 'SPK';
keys.spike_field.channel       = 'LFP'; % selected LFP channels
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


for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !
    if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
    end
    version_folder=keys.project_version;
    keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    %    keys.monkeys=keys.batching.monkeys;
    %     if keys.batching.combine_monkeys
    %         keys.batching.monkeys={''};
    %     end
    %
    %     for m=1:numel(keys.batching.monkeys)
    %         keys.monkey=keys.batching.monkeys{m};
    keys.anova_table_file=[keys.drive '\Projects\' project '\ephys\' keys.project_version '\tuning_table_combined_CI.mat'];
    %         if any(ismember(population_analysis_to_perform,{'ons','pop','gaz'}))
    %             population=ph_load_population([keys.drive filesep keys.basepath_to_save filesep keys.project_version],['population_' keys.monkey]);
    %             population=ph_assign_perturbation_group(keys,population);
    %             population=ph_epochs(population,keys);
    %         else
    %             population=[];
    %         end
    %         for t=1:numel(keys.batching.targets)
    %             target=keys.batching.targets{t};
    %
    %             for subregion=1:keys.batching.n_Subregions
    %                 %population_selection={};
    %                 if keys.batching.Subregions_separately
    %                     population_selection            = {'target',target;'Subregion', subregion};
    %                 else
    %                     population_selection            = {'target',target};
    %                 end
    %
    for a=1:numel(keys.position_and_plotting_arrangements)
        
        keys.arrangement=keys.position_and_plotting_arrangements{a};
        
        spike_field=compute_per_session(keys);
        save([keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'spike_field_' keys.arrangement '_' keys.spike_field.method],'spike_field','keys');
        
    end
    %             end
    %         end
    %end
end
end

function spike_field=compute_per_session(keys)
session_folder=[keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep];
dir_session_folder=dir([session_folder '*.mat']);
names={dir_session_folder.name};
population_names=names(cellfun(@(x) any(strfind(x,['population_'])),names));

tic
n_site=1;
for sess=1:numel(population_names)
    clear population sites
    load([session_folder population_names{sess}]);
    load([session_folder 'sites' population_names{sess}(11:end)]);
    
    %% get epochs & perturbation group
    population=ph_assign_perturbation_group(keys,population);
    population=ph_epochs(population,keys);
    
    %All_site_IDs={sites.site_ID};
    for s=1:numel(sites);
        clear units_with_site
        keep_site=0;
        
        site=sites(s);
        [site.trial.accepted]=deal(true); %% needs to be defined somewhere before!
        
        % takes around 1 min
        complete_LFP=[site.trial.LFP];
        for f=1:numel(keys.LFP.Morlet_Matrix)
            filtered_LFP=conv(complete_LFP,keys.LFP.Morlet_Matrix{f},'same');
            n_sample_start=1;
            for t=1:numel(site.trial)
                n_sample_end=n_sample_start+numel(site.trial(t).LFP)-1;
                site.trial(t).filtered(:,f)=filtered_LFP(n_sample_start:n_sample_end);
                n_sample_start=n_sample_end+1;
            end
        end
        clear complete_LFP filtered_LFP
        
        %         for t=1:numel(site.trial)
        %                 site.trial(t).filtered(:,f)=conv(site.trial(t).LFP,keys.LFP.Morlet_Matrix(f,:),'same');
        %        end
        n_unit=0;
        for u=1:numel(population)
            unit=population(u);
            keep_unit=0;
            s_tr_mat=[site.trial.block; site.trial.run; site.trial.n];
            u_tr_mat=[unit.trial.block; unit.trial.run; unit.trial.n];
            s_idx=ismember(s_tr_mat',u_tr_mat','rows');
            u_idx=ismember(u_tr_mat',s_tr_mat','rows');
            
            if ~any(s_idx) || ~any(u_idx)
                continue
            end
            
            unitC.trial=unit.trial(u_idx);
            siteC.trial=site.trial(s_idx);
            
            % for trials that are excluded either from spikes or LFP
            accepted=[unitC.trial.accepted] & [siteC.trial.accepted] & [unitC.trial.success] & ~isnan([unitC.trial.perturbation]);
            if ~any(accepted)
                continue
            else
                n_unit=n_unit+1;
                keep_site=1;
            end
            unitC=ph_arrange_positions_and_plots(unitC.trial(accepted),keys);
            %unitC.trial=unitC.trial(accepted);
            siteC.trial=siteC.trial(accepted);
            
            nspikesaccum=0;
            clear stsConvol states_per_spike
            stsConvol.lfplabel={keys.spike_field.channel};
            stsConvol.label={keys.spike_field.spikechannel};
            stsConvol.freq=keys.LFP.frequencies;
            stsConvol.dimord='{chan}_spike_lfpchan_freq';
            sz=size(vertcat(unitC.trial.arrival_times));
            stsConvol.fourierspctrm{1}  = NaN([sz numel(keys.LFP.Morlet_Matrix)]); % preallocation... saves TONS of time!!
            %             stsConvol.trial{1}          = NaN(sz);
            %             stsConvol.time{1}           = NaN(sz);
            %             stsConvol.trialtime         = NaN(numel(siteC.trial),2);
            for t=1:numel(siteC.trial)
                timeline=siteC.trial(t).TDT_LFPx_tStart:(1/siteC.trial(t).TDT_LFPx_SR):(numel(siteC.trial(t).LFP)/siteC.trial(t).TDT_LFPx_SR+siteC.trial(t).TDT_LFPx_tStart);
                AT=unitC.trial(t).arrival_times;
                AT(AT>=timeline(end-1))=[];
                AT(AT<timeline(1))=[];
                timestamps=round((AT-siteC.trial(t).TDT_LFPx_tStart)*siteC.trial(t).TDT_LFPx_SR)+1;
                nspikesend= nspikesaccum+size(AT,1);
                % input for FT
                stsConvol.fourierspctrm{1}(nspikesaccum+1:nspikesend,1,:)=siteC.trial(t).filtered(timestamps,:);
                stsConvol.trial{1}(nspikesaccum+1:nspikesend,1)=t;
                stsConvol.time{1}(nspikesaccum+1:nspikesend,1)=AT;
                [states_per_spike(nspikesaccum+1:nspikesend).states]=deal(unitC.trial(t).states);%,1,nspikesend-nspikesaccum);
                [states_per_spike(nspikesaccum+1:nspikesend).states_onset]=deal(unitC.trial(t).states_onset);%,1,nspikesend-nspikesaccum);
                stsConvol.trialtime(t,:)=[timeline(1) timeline(end)];
                nspikesaccum=nspikesend;
            end
            stsConvol.fourierspctrm{1}=stsConvol.fourierspctrm{1}(1:nspikesend,1,:);
            %             stsConvol.trial{1}        =stsConvol.trial{1}(1:nspikesend,1);
            %             stsConvol.time{1}         =stsConvol.time{1}(1:nspikesend,1);
            %
            % compute the statistics on the phases
            % statSts           = ft_spiketriggeredspectrum_stat(keys.spike_field,stsConvol);
            
            %% conditions and epochs
            clear unique_conditions trial_matrix
            split_conditions={'type','effector','reach_hand','choice','hemifield','perturbation'}; % hemifield missing, wait for arrange_positions...
            for par=1:numel(split_conditions)
                fn=split_conditions{par};
                unique_conditions{par}=unique([unitC.trial.(fn)]);
                unique_conditions{par}(isnan(unique_conditions{par}))=[]; % perturbation not defined everywhere??
                trial_matrix(:,par)=[unitC.trial.(fn)]';
            end
            
            condition_matrix    = combvec(unique_conditions{:})';
            condition_matrix    = condition_matrix(ismember(condition_matrix,trial_matrix,'rows'),:);
            % loop through conditions
            sts_con=stsConvol;
            sts_epo=stsConvol;
            clear per_condition
            for con=1:size(condition_matrix,1)
                clear per_epoch
                for par=1:numel(split_conditions)
                    per_condition(con).(split_conditions{par})=condition_matrix(con,par);
                end
                tr_ind=find(ismember(trial_matrix,condition_matrix(con,:),'rows'));
                sp_ind=ismember(stsConvol.trial{1},tr_ind);
                sts_con.trial{1}=stsConvol.trial{1}(sp_ind); %% re-order trial numbers...
                sts_con.time{1}=stsConvol.time{1}(sp_ind);
                sts_con.fourierspctrm{1}=stsConvol.fourierspctrm{1}(sp_ind,1,:);
                %sts_con.trialtime=stsConvol.trialtime(tr_ind,:);
                
                
                event_onsets=vertcat(states_per_spike(sp_ind).states_onset);
                events=vertcat(states_per_spike(sp_ind).states);
                
                LFP_per_condition=siteC.trial(tr_ind); %% only needed for LFP
                EPOCHS=keys.EPOCHS_PER_TYPE{condition_matrix(con,1)};
                for e=1:size(EPOCHS,1)
                    event=EPOCHS{e,2};
                    event_onset=event_onsets(events==event);
                    ep_ind=sts_con.time{1}>=event_onset+EPOCHS{e,3} & sts_con.time{1}<=event_onset+EPOCHS{e,4};
                    [states_per_spike(nspikesaccum+1:nspikesend).states_onset]=deal(unitC.trial(t).states_onset);
                    sts_epo.trial{1}=sts_con.trial{1}(ep_ind);
                    sts_epo.time{1}=sts_con.time{1}(ep_ind);
                    sts_epo.fourierspctrm{1}=sts_con.fourierspctrm{1}(ep_ind,1,:);
                    % minumum pairs criterion
                    n_pairs=0;
                    for temp_t=unique(sts_epo.trial{1})'
                        n_pairs=n_pairs+sum(sts_epo.trial{1}~=temp_t)*sum(sts_epo.trial{1}==temp_t)/2;
                    end
                    
                    statStsEp           = ft_spiketriggeredspectrum_stat(keys.spike_field,sts_epo);
                    
                    current_epoch.(keys.spike_field.method)=statStsEp.(keys.spike_field.method);
                    current_epoch.epoch=EPOCHS{e,1};
                    current_epoch.n_pairs=n_pairs;
                    if numel(current_epoch.(keys.spike_field.method))==1 %|| n_pairs<100 %% no spike
                        current_epoch.(keys.spike_field.method)=NaN(1,numel(keys.LFP.Morlet_Matrix));
                    elseif any(isinf(current_epoch.(keys.spike_field.method))) %% one spike?
                        disp(['INF: ' num2str(numel(isinf(current_epoch.(keys.spike_field.method))))]);
                        current_epoch.(keys.spike_field.method)(isinf(current_epoch.(keys.spike_field.method)))=NaN;
                    end
                    
                    %% LFP part (to take out later, because it does not belong here at all!!
                    clear power_spect
                    valid_t=true(numel(LFP_per_condition),1);
                    for temp_t=1:numel(LFP_per_condition)
                        state_on=LFP_per_condition(temp_t).states_onset(LFP_per_condition(temp_t).states==event)-LFP_per_condition(temp_t).TDT_LFPx_tStart;
                        if isnan(state_on)
                           valid_t(temp_t)=false; 
                           continue
                        end
                        temp_indexes=round((state_on+EPOCHS{e,3}:1/LFP_per_condition(temp_t).TDT_LFPx_SR:state_on+EPOCHS{e,4})*LFP_per_condition(temp_t).TDT_LFPx_SR);
                        power_spect(temp_t,:)=nanmean(abs(LFP_per_condition(temp_t).filtered(temp_indexes,:)),1);
                    end
                    current_epoch.power_sepctrum=nanmean(power_spect(valid_t,:),1);
                    
                    per_epoch(e)=current_epoch;
                end
                per_condition(con).per_epoch=per_epoch;
            end
            % if keep_unit
                unit=rmfield(unit,'trial');
                %unit.(keys.spike_field.method)=statSts.(keys.spike_field.method);
                unit.per_condition=per_condition;
                units_with_site(n_unit)=unit;
%             else
%                 n_unit=n_unit-1;
%             end
        end
        
        site=rmfield(site,'trial');
        
        if ~keep_site
            continue
        else
            site.unit=units_with_site;
            spike_field(n_site)=site;
            n_site=n_site+1;
        end
    end
    
end
toc

end
%
% function  [amp,phase]=spike_freq_coherence(arrival_times,LFP,SR,s1,s2,kernel)
%
% samples =0:(1/SR):(size(LFP,2)-s1)/SR;
% SD= conv(hist(arrival_times,samples),normpdf(-5*kernel:1/SR:5*kernel,0,kernel),'same');
% filtered_cut=LFP(:,s1:s2);
% amp=abs(filtered_cut)*SD'/numel(arrival_times);
% phase=angle(exp(1i*angle(filtered_cut))*SD');
% end
