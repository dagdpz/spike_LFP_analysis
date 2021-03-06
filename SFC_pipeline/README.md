# Spike-Field Coherence analysis (e.g. of the macaque parietal cortex activity within and across hemispheres)  

## Main scripts  
**ph_run_LFP_analysis.m** - the function for ppc calculation (part of 'Lukas' DAG ephys pipeline) 

Data in specific folder on the server, e.g. Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn_SF (here SF stands for spike-field)

It contains, for each session, separate mat files for spikes and LFPs (which are the output of ...LUKAS), output is spike_field_hands_inactivation_ppc1.mat

## Plotting functions
**MS_plot_PPC.m** - plot ppc before and after inactivation for each epoch and condition  
**MS_plot_ppc_perunit.m** - plot ppc before and after inactivation for each epoch and condition PER UNIT (averaged across LFP channels)  
**MS_scatter_per_unit.m** - plot *scatter plots* (one dot - one unit) before and after inactivation to see the overall change of tuning  

## Additional functions (called from MS_plot_PPC)  
**MS_perm_test.m** - Cluster-adjusted permutation t-test (as described in Maris & Oostenveld, 2007)  
**MS_find_clusters.m** - accessory function for finding clusters in permutation test  
plus
**MS_FDR_test.m** - False Discovery Rate (FDR) correction for multiple comparisons


### Example usage:  

1. run the spike analysis pipeline (with desired settings) from github. e.g Run ph_initiation('PPC_pulv_eye_hand',{'MIP_dPul_inj_working_post_sfn_Sarath'})



2. Run ph_run_LFP_analysis('PPC_pulv_eye_hand',{'MIP_dPul_inj_working_post_sfn_Sarath'})

3. For first time, set calculate_ttest = 1 inside MS_plot_PPC.m and run the following code snippet: 


4.
```
Loc = 'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn_Sarath\'; 
PPC_method='ppc1';  
hemis = {'MIP_L', 'MIP_R', 'MIP_L_MIP_R', 'MIP_R_MIP_L'};  

for hemi = hemis,  
    disp(hemi)  
    hemi = hemi{:};  
    MS_plot_PPC(Loc, PPC_method, hemi)  
    MS_plot_ppc_perunit(Loc, PPC_method, hemi)  
    MS_scatter_per_unit(Loc, PPC_method, hemi)  
end
```
