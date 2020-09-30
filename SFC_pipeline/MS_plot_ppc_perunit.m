function MS_plot_ppc_perunit( spike_field_location, PPC_method, hemi )
%	Inputs: SFfile_location, PPC_method, hemi

%   PPC_method is 'ppc0' or 'ppc1'
%   hemi is MIP_L, MIP_R, MIP_R_MIP_L, MIP_L_MIP_R

% Output: plots of ppc values before and after inactivaion for each
% condition and epoch separately PER UNIT (averaged across LFP channels)

% load(spike_field_location);
[spike_field keys]= ph_load_population_ppc(spike_field_location, 'spike_field_');

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

%% Conditions to compare keys
keys.conditions_to_compare{1}(1).reach_hand=1;
keys.conditions_to_compare{1}(1).hemifield=-1;
keys.conditions_to_compare{1}(1).choice=0;
keys.conditions_to_compare{1}(1).perturbation=0;
keys.conditions_to_compare{1}(1).color='b';
keys.conditions_to_compare{1}(1).title='LHLS';

keys.conditions_to_compare{1}(2).reach_hand=1;
keys.conditions_to_compare{1}(2).hemifield=-1;
keys.conditions_to_compare{1}(2).choice=1;
keys.conditions_to_compare{1}(2).perturbation=0;
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
keys.conditions_to_compare{2}(2).choice=1;
keys.conditions_to_compare{2}(2).perturbation=0;
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
keys.conditions_to_compare{3}(2).choice=1;
keys.conditions_to_compare{3}(2).perturbation=0;
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
keys.conditions_to_compare{4}(2).choice=1;
keys.conditions_to_compare{4}(2).perturbation=0;
keys.conditions_to_compare{4}(2).color='r';
keys.conditions_to_compare{4}(2).title='RHRS';

epochs = [keys.EPOCHS_PER_TYPE{4}(:,1)'];


clear all_cond
Parameters={'type','effector','reach_hand','choice','perturbation','hemifield'};
for par=1:numel(Parameters),
    SFidx.(Parameters{par})=arrayfun(@(x) isfield(x.per_condition,Parameters{par}) && numel(unique([x.per_condition.(Parameters{par})]))>1,SF_difchan_hemi);
end

limit_pairs_by_condition={'hemifield'};
valid_SF_combinations=true(size(SF_difchan_hemi));
for L=1:numel(limit_pairs_by_condition)
    valid_SF_combinations = valid_SF_combinations &   SFidx.(limit_pairs_by_condition{L});
end

% Plotting population ppc
n_unique_units=numel(unique({SF_difchan_hemi(valid_SF_combinations).unit_ID}));
n_unique_sites=numel(unique({SF_difchan_hemi(valid_SF_combinations).site_ID}));
disp(['unique_sites:' num2str(n_unique_sites) ' ,unique_units:' num2str(n_unique_units)]);
% SF_comb_per_cond = [SF_difchan_hemi(valid_SF_combinations).per_condition];
% 
% for par=1:numel(Parameters),
%     all_cond(:, par) = [SF_comb_per_cond.(Parameters{par})];
%     paridx.(Parameters{par}) = par;
% end

SF_difchan_hemi_valid = SF_difchan_hemi(valid_SF_combinations);
LFP_siteID_difchan_hemi = LFP_siteID_difchan_hemi(valid_SF_combinations);
unit_ids_all = unit_ids_all(valid_SF_combinations);
uniq_units = unique(unit_ids_all);


for uniq_un=1:numel(uniq_units),
    un_ind = find(strcmp(unit_ids_all, uniq_units(uniq_un)));
    SF_difchan_hemi_tmp = SF_difchan_hemi_valid(un_ind);
    
    %% Finding the same set of SF pairs within each epoch

    SF_same_set = zeros(numel(un_ind), numel(epochs));
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
    if all(sum(SF_same_set) == 0),
        continue
    end
%
    keys.path_to_save = sprintf('%s%s\\SFC_perunit_%s\\Mean_per_unit', [keys.basepath_to_save], keys.project_version, hemi);
    %keys.path_to_save = sprintf('%s\\%s\\SFC_perunit_%s\\Per_site', [keys.drive, keys.basepath_to_save], keys.project_version, hemi);

    if exist(keys.path_to_save) == 0;
        mkdir(keys.path_to_save);
    end

    FigName_short2= sprintf('%d_%s', uniq_un, uniq_units{uniq_un}); %name of the file
    FigName_Long2=sprintf('%d. %s, %s, %d pairs, unit %s, rating %d/%d/%d', ...
        uniq_un, hemi, PPC_method, numel(un_ind), uniq_units{uniq_un}, ...
        SF_difchan_hemi_tmp(1).SNR_rating, SF_difchan_hemi_tmp(1).Single_rating, SF_difchan_hemi_tmp(1).stability_rating); % overall title of the figure
    set(0,'DefaultTextInterpreter','none');
    wanted_papersize= [40 25];
    h1=figure('outerposition',[0 0 1900 1200],'name', FigName_Long2);
    figure(h1)

    clear H

    %% PLOT

    for ep=1:numel(epochs)
        if sum(SF_same_set(:, ep)) == 0,
            continue
        end
        y_max = 0;
        y_min = 0;
        y_bar = 0;
        max_cond = 1;
        unique_un = numel(unique({SF_difchan_hemi_tmp(logical(SF_same_set(:, ep))).unit_ID}));

        SF_comb_per_cond = [SF_difchan_hemi_tmp(logical(SF_same_set(:, ep))).per_condition];

        clear all_cond
        for par=1:numel(Parameters),
            all_cond(:, par) = [SF_comb_per_cond.(Parameters{par})];
            paridx.(Parameters{par}) = par;
        end

        for row=1:numel(keys.conditions_to_compare)
            condition_fieldnames=fieldnames(keys.conditions_to_compare{row});
            condition_fieldnames=condition_fieldnames(ismember(condition_fieldnames,Parameters));

            for con=1:numel(keys.conditions_to_compare{row}) % this is ugly, but the same loop again
                current_index=true(size(all_cond,1),1);
                for FN=condition_fieldnames'
                    current_index=current_index & all_cond(:, paridx.(FN{:}))==keys.conditions_to_compare{row}(con).(FN{:});
                end
                SF_comb_per_cond_per_epoch = vertcat(SF_comb_per_cond(current_index).per_epoch);
                SF_comb_per_cond_per_epoch = SF_comb_per_cond_per_epoch(:, ep); % not very smart

                ax(row) = subplot(numel(keys.conditions_to_compare), numel(epochs), numel(epochs)*(row-1)+ep); hold on;

                %temptemp=arrayfun(@(x) numel(x.(PPC_method)),SF_comb_per_cond_per_epoch(:,ep))>1;
                PPC=vertcat(SF_comb_per_cond_per_epoch.(PPC_method));
                PPC(isinf(PPC))=NaN;

                if max(mean(PPC, 1)) > y_max, 
                    y_max = max(mean(PPC, 1) + 3*max(sem(PPC,1))); 
                    y_bar = max(mean(PPC, 1) + 1*max(sem(PPC,1))); 
                    max_cond = row; 
                end
                if min(mean(PPC, 1)) < y_min, y_min = min(mean(PPC, 1) - max(sem(PPC,1))); end
                
                %plot(keys.LFP.frequencies, PPC, keys.conditions_to_compare{row}(con).color);
                H(ep) = shadedErrorBar(keys.LFP.frequencies, mean(PPC, 1), sem(PPC,1), keys.conditions_to_compare{row}(con).color,1);
                y_lim=get(gca,'ylim');
                set(gca,'XScale','log') 

                % text(min(keys.LFP.frequencies),y_lim(2)-con*diff(y_lim)/10,['Excluded:' num2str(sum(all(isnan(PPC),2)) + sum(~temptemp))],'color', keys.conditions_to_compare{row}(con).color);
                if row==1 && con==2
                    title(sprintf('%s, %d ms \nSF pairs: %d; Units: %d \n', epochs{ep}, round((keys.EPOCHS_PER_TYPE{4}{ep, 4} - keys.EPOCHS_PER_TYPE{4}{ep, 3})*1000),...
                        numel(PPC(:, 1)), unique_un));
    %                 text(min(keys.LFP.frequencies),y_lim(2)-con*diff(y_lim)/10,...
    %                     ['SF_pairs:' num2str(numel(PPC(:, 1))) newline 'Units:' num2str(unique_un)], 'Color', 'r');
                end

                if con==1 && ep==1
                    ylabel([keys.conditions_to_compare{row}(con).title]);
                end

                if con == 2,
                    set(gca,'XTick',[bands_v]); %round(10.^(0.6:0.2:2))
    %                 set(gca,'YGrid', 'off');
                end
                xlim([min(keys.LFP.frequencies), max(keys.LFP.frequencies)]);
            end
        end
        set(ax(1:4), 'ylim', [y_min y_max]);
        for sp = 1:row,
            subplot(ax(sp));
            z = vline([bands_v]); %, 'r:'
            set(z, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'LineStyle', ':');
        end
    end

    set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
    mtit(figure(h1),  [FigName_Long2 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
    export_fig(figure(h1), [keys.path_to_save filesep FigName_short2], '-pdf','-transparent') % pdf by run

    close all
    
end
end
