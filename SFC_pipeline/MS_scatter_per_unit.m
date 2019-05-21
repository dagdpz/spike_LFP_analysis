function MS_scatter_per_unit( spike_field_location, PPC_method, hemi )
%	Inputs: SFfile_location, PPC_method, hemi

%   PPC_method is 'ppc0' or 'ppc1'
%   hemi is MIP_L, MIP_R, MIP_R_MIP_L, MIP_L_MIP_R

% Output: Scatter plots of ppc values before and after inactivaion for each
% condition and epoch separately PER UNIT (averaged across LFP channels)

load(spike_field_location);

keys.bands = {...
    'theta', 4, 8 ;...
    'alpha', 8, 12;...
    'beta', 12, 30;...
    'gamma', 30, 100;...
    };

spike_pairs = 2500; 
bands_v = [4, 8, 12, 30, 100];
SF_combinations = [spike_field.unit];

siteIDs_tmp = arrayfun(@(x) repmat({x.site_ID}, numel(x.unit),1), spike_field, 'uniformoutput',false);
LFP_siteIDs = vertcat(siteIDs_tmp{:});

hemi_tmp = arrayfun(@(x) repmat({x.target}, numel(x.unit),1), spike_field, 'uniformoutput',false);
LFP_hemiIDs = (vertcat(hemi_tmp{:}));

% getting info about the same site unit and LFP, 1-same, 0-different
same_site_ind = strcmp(LFP_siteIDs, {SF_combinations.site_ID}');

% getting info about hemisphere - ACROSS - LFP and spikes are in different
% hemispheres
for ind=1:numel(LFP_hemiIDs),
    if strcmp(LFP_hemiIDs(ind), SF_combinations(ind).target),
        same_hemi_ind_t(ind) = LFP_hemiIDs(ind);
    else
        same_hemi_ind_t(ind)= {[SF_combinations(ind).target, '_', LFP_hemiIDs{ind}]}; %{'ACROSS'};
    end
end
same_hemi_ind = same_hemi_ind_t';

% choose hemisphere AND remove same-site pairs
subset_ind = strcmp(same_hemi_ind, {hemi}) & ~same_site_ind;

% subset of all the pairs from indicated hemi and different sites
SF_difchan_hemi = SF_combinations(subset_ind);
LFP_siteID_difchan_hemi = LFP_siteIDs(subset_ind); %LFP channels
unit_ids_all = [{SF_difchan_hemi.unit_ID}]';

% [b,m,n] = unique([{SF_difchan_hemi.unit_ID}]');
% single_rate = [SF_difchan_hemi(m).Single_rating];
% sum(single_rate==1) %SU
% numel(single_rate) - sum(single_rate==1) %MU

%% Conditions to compare keys
keys.conditions_to_compare{1}(1).reach_hand=1;
keys.conditions_to_compare{1}(1).hemifield=-1;
keys.conditions_to_compare{1}(1).choice=0;
keys.conditions_to_compare{1}(1).perturbation=0;
keys.conditions_to_compare{1}(1).color='b';
keys.conditions_to_compare{1}(1).title='LHLS';

keys.conditions_to_compare{1}(2).reach_hand=1;
keys.conditions_to_compare{1}(2).hemifield=-1;
keys.conditions_to_compare{1}(2).choice=0;
keys.conditions_to_compare{1}(2).perturbation=1;
keys.conditions_to_compare{1}(2).color='r';
keys.conditions_to_compare{1}(2).title='LHLS';


keys.conditions_to_compare{2}(1).reach_hand=1;
keys.conditions_to_compare{2}(1).hemifield=1;
keys.conditions_to_compare{2}(1).choice=0;
keys.conditions_to_compare{2}(1).perturbation=0;
keys.conditions_to_compare{2}(1).color='b';
keys.conditions_to_compare{2}(1).title='LHRS';

keys.conditions_to_compare{2}(2).reach_hand=1;
keys.conditions_to_compare{2}(2).hemifield=1;
keys.conditions_to_compare{2}(2).choice=0;
keys.conditions_to_compare{2}(2).perturbation=1;
keys.conditions_to_compare{2}(2).color='r';
keys.conditions_to_compare{2}(2).title='LHRS';


keys.conditions_to_compare{3}(1).reach_hand=2;
keys.conditions_to_compare{3}(1).hemifield=-1;
keys.conditions_to_compare{3}(1).choice=0;
keys.conditions_to_compare{3}(1).perturbation=0;
keys.conditions_to_compare{3}(1).color='b';
keys.conditions_to_compare{3}(1).title='RHLS';

keys.conditions_to_compare{3}(2).reach_hand=2;
keys.conditions_to_compare{3}(2).hemifield=-1;
keys.conditions_to_compare{3}(2).choice=0;
keys.conditions_to_compare{3}(2).perturbation=1;
keys.conditions_to_compare{3}(2).color='r';
keys.conditions_to_compare{3}(2).title='RHLS';


keys.conditions_to_compare{4}(1).reach_hand=2;
keys.conditions_to_compare{4}(1).hemifield=1;
keys.conditions_to_compare{4}(1).choice=0;
keys.conditions_to_compare{4}(1).perturbation=0;
keys.conditions_to_compare{4}(1).color='b';
keys.conditions_to_compare{4}(1).title='RHRS';

keys.conditions_to_compare{4}(2).reach_hand=2;
keys.conditions_to_compare{4}(2).hemifield=1;
keys.conditions_to_compare{4}(2).choice=0;
keys.conditions_to_compare{4}(2).perturbation=1;
keys.conditions_to_compare{4}(2).color='r';
keys.conditions_to_compare{4}(2).title='RHRS';

epochs = [keys.EPOCHS_PER_TYPE{4}(:,1)'];

clear all_cond
Parameters={'type','effector','reach_hand','choice','perturbation','hemifield'};
for par=1:numel(Parameters),
    SFidx.(Parameters{par})=arrayfun(@(x) isfield(x.per_condition,Parameters{par}) && numel(unique([x.per_condition.(Parameters{par})]))>1,SF_difchan_hemi);
end

limit_pairs_by_condition={'perturbation'};
valid_SF_combinations=true(size(SF_difchan_hemi));
for L=1:numel(limit_pairs_by_condition)
    valid_SF_combinations = valid_SF_combinations &   SFidx.(limit_pairs_by_condition{L});
end

n_unique_units=numel(unique({SF_difchan_hemi(valid_SF_combinations).unit_ID}));
n_unique_sites=numel(unique({SF_difchan_hemi(valid_SF_combinations).site_ID}));
disp(['unique_sites:' num2str(n_unique_sites) ' ,unique_units:' num2str(n_unique_units)]);

LFP_siteID_difchan_hemi = LFP_siteID_difchan_hemi(valid_SF_combinations);
unit_ids_all = unit_ids_all(valid_SF_combinations);
uniq_units = unique(unit_ids_all);

%% Finding the same set of SF pairs within each epoch
SF_difchan_hemi_tmp = SF_difchan_hemi(valid_SF_combinations);
SF_same_set = zeros(sum(valid_SF_combinations), numel(epochs));
for un=1:numel(SF_difchan_hemi_tmp),
    un_tmp = [SF_difchan_hemi_tmp(un).per_condition];
    un_tmp = un_tmp(logical(~[un_tmp.choice]));
    num_p = zeros(numel(un_tmp), numel(epochs));
    for cnd=1:numel(un_tmp),
        ep_tmp = [un_tmp(cnd).per_epoch];
        num_p(cnd, :) = [ep_tmp.n_pairs];
    end
    SF_same_set(un, :) = all(num_p > spike_pairs, 1);
end
% sum(SF_same_set)

%
keys.path_to_save = sprintf('%s\\%s\\SFC_scatter_%s', [keys.drive, keys.basepath_to_save], keys.project_version, PPC_method);

if exist(keys.path_to_save) == 0;
    mkdir(keys.path_to_save);
end

%%% loop for bands goes here
for band = 1:numel(keys.bands(:,1)),
    band_name = keys.bands(band,1);
    disp(sprintf('Plotting scatter for %s band', band_name{:}));
    
    FigName_short2= sprintf('%s_%s_scatter_%s', hemi,PPC_method, band_name{:}); %name of the file
    FigName_Long2=sprintf('%s_%s, %d SF pairs, %s band (%d-%d), %d units, %d sites', hemi, PPC_method, ...
            sum(valid_SF_combinations), band_name{:}, keys.bands{band, 2}, keys.bands{band, 3}, n_unique_units, n_unique_sites); % overall title of the figure
    set(0,'DefaultTextInterpreter','none');
    wanted_papersize= [40 25];
    h1=figure('outerposition',[0 0 1900 1200],'name', FigName_Long2);
    figure(h1)

    clear H

%% PLOT

    for ep=1:numel(epochs)

        unique_un = numel(unique({SF_difchan_hemi_tmp(logical(SF_same_set(:, ep))).unit_ID}));
        SF_set_epoch = SF_difchan_hemi_tmp(logical(SF_same_set(:, ep)));
        SF_comb_per_cond = [SF_set_epoch.per_condition];

        LFP_siteID_ep = LFP_siteID_difchan_hemi(logical(SF_same_set(:, ep)));
        unit_ids_ep = unit_ids_all(logical(SF_same_set(:, ep)));
        uniq_units = unique(unit_ids_ep);

        rand_col = colormap(hsv(unique_un));
        %%%

        clear all_cond
        for par=1:numel(Parameters),
            all_cond(:, par) = [SF_comb_per_cond.(Parameters{par})];
            paridx.(Parameters{par}) = par;
        end

        %%% prepare structure here
        % or do this per unit here
        clear one_epoch_scatter
        min_ppc = 0;
        max_ppc = 0;
        for row=1:numel(keys.conditions_to_compare)
            condition_fieldnames=fieldnames(keys.conditions_to_compare{row});
            condition_fieldnames=condition_fieldnames(ismember(condition_fieldnames,Parameters));

            ppc_per_plot = zeros(numel(uniq_units),2);
            for con=1:numel(keys.conditions_to_compare{row}) % this is ugly, but the same loop again
                current_index=true(size(all_cond,1),1);
                for FN=condition_fieldnames'
                    current_index=current_index & all_cond(:, paridx.(FN{:}))==keys.conditions_to_compare{row}(con).(FN{:});
                end
                SF_comb_per_cond_per_epoch = vertcat(SF_comb_per_cond(current_index).per_epoch);
                SF_comb_per_cond_per_epoch = SF_comb_per_cond_per_epoch(:, ep); % not very smart
                % 116 pairs - here i can select 

                %temptemp=arrayfun(@(x) numel(x.(PPC_method)),SF_comb_per_cond_per_epoch(:,ep))>1;
                PPC=vertcat(SF_comb_per_cond_per_epoch.(PPC_method));
                PPC(isinf(PPC))=NaN;
                
                % choose specific band now
                band_ind = keys.LFP.frequencies > keys.bands{band, 2} & keys.LFP.frequencies <= keys.bands{band, 3};                    
                PPC_band = PPC(:, band_ind); % select bands here
                
                ppc_unit = zeros(numel(uniq_units),  numel(PPC_band(1, :)));
                % now let`s average units
                for uu=1:numel(uniq_units), %unique unit
                    ind_st = find(strcmp(unit_ids_ep, uniq_units(uu)));
                    ppc_unit(uu, :) = mean(PPC_band(ind_st, :));
                end     
                ppc_per_plot(:, con) = mean(ppc_unit, 2);   
            end

            if min(min(ppc_per_plot)) < min_ppc,
                min_ppc = min(min(ppc_per_plot));
            end
            if max(max(ppc_per_plot)) > max_ppc,
                max_ppc = max(max(ppc_per_plot));
            end

            one_epoch_scatter{row} = ppc_per_plot;
        end

        for rrow = 1:numel(keys.conditions_to_compare),
            ax(rrow) = subplot(numel(keys.conditions_to_compare), numel(epochs), numel(epochs)*(rrow-1)+ep); hold on;
            per_plot = one_epoch_scatter{rrow};
            for uu=1:numel(uniq_units), %unique unit
                scatter(per_plot(uu, 1), per_plot(uu, 2), [], [rand_col(uu, :)]);  hold on; 
            end
            plot(min_ppc:0.0001:max_ppc, min_ppc:0.0001:max_ppc)
            ylim([min_ppc max_ppc])
            xlim([min_ppc max_ppc])
            if ep == 1;
                ylabel([keys.conditions_to_compare{rrow}(con).title]);
            end

            if rrow == 1;
                %title(sprintf('%s, %d ms', epochs{ep}, round((keys.EPOCHS_PER_TYPE{4}{ep, 4} - keys.EPOCHS_PER_TYPE{4}{ep, 3})*1000) )); % ???????
                title(sprintf('%s, %d ms \nSF pairs: %d; Units: %d \n', epochs{ep}, round((keys.EPOCHS_PER_TYPE{4}{ep, 4} - keys.EPOCHS_PER_TYPE{4}{ep, 3})*1000),...
                        numel(PPC(:, 1)), unique_un));
            end
        end

    end

    set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
    mtit(figure(h1),  [FigName_Long2 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
    export_fig(figure(h1), [keys.path_to_save filesep FigName_short2], '-pdf','-transparent') % pdf by run

    close all
end

end
