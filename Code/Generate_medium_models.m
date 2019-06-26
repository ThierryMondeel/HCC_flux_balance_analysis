%% Settings and setup
format compact

initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

%% Load the patched Recon3D
Recon3DModel = readCbModel('../Models/Recon3DModel_301_patched.mat');

%% Save medium specific models
media = {'1','2','3','4','5','6'};
for i = 1:6
    medium = media{i};
    
    R3 = Recon3DModel; % reset R3
    
    % Update specific bounds defined by our experimental media
    allowed_IDs = {'EX_gly[e]'; 'EX_arg_L[e]'; 'EX_cys_L[e]'; ...
                'EX_his_L[e]'; 'EX_ile_L[e]'; 'EX_leu_L[e]'; 'EX_lys_L[e]'; 'EX_met_L[e]'; ...
                'EX_phe_L[e]'; 'EX_ser_L[e]'; 'EX_thr_L[e]'; 'EX_trp_L[e]'; 'EX_tyr_L[e]'; ...
                'EX_val_L[e]'; 'EX_chol[e]'; 'EX_pnto_R[e]'; 'EX_fol[e]'; 'EX_ncam[e]'; ...
                'EX_pydxn[e]'; 'EX_ribflv[e]'; 'EX_thm[e]'; 'EX_inost[e]'; 'EX_ca2[e]'; ...
                'EX_cl[e]'; 'EX_fe3[e]'; 'EX_so4[e]'; 'EX_k[e]'; 'EX_na1[e]'; 'EX_hco3[e]'; ...
                'EX_pi[e]'; 'EX_h2o[e]'; 'EX_o2[e]';'EX_co2[e]'}; 
    allowed_lbs = -1*[0.4; 0.398; 0.201; ...
        0.2; 0.801; 0.801; 0.798; 0.201; ... 
        0.4; 0.4; 0.798; 0.0784; 0.398; ...
        0.803; 0.028; 0.0083; 0.0091; 0.0328; ...
        0.0194; 0.0011; 0.011; 0.04; 1.801;
        117.4; 0.00024; 0.814; 5.33; 155.0; 44.04; ...
        0.9; 55560; 1000; 1000];

    R3 = changeRxnBounds(R3, allowed_IDs, allowed_lbs, 'l');
    
    switch medium
        case '1' % high glucose, glutamine
            R3 = changeRxnBounds(R3, 'EX_glc_D[e]', -25, 'l');
            R3 = changeRxnBounds(R3, 'EX_gln_L[e]', -4, 'l');
        case '2' % high glucose, no glutamine
            R3 = changeRxnBounds(R3, 'EX_glc_D[e]', -25, 'l');
            R3 = changeRxnBounds(R3, 'EX_gln_L[e]', 0, 'l');
        case '3' % low glucose, glutamine
            R3 = changeRxnBounds(R3, 'EX_glc_D[e]', -5.55, 'l');
            R3 = changeRxnBounds(R3, 'EX_gln_L[e]', -4, 'l');
        case '4' % low glucose, no glutamine
            R3 = changeRxnBounds(R3, 'EX_glc_D[e]', -5.55, 'l');
            R3 = changeRxnBounds(R3, 'EX_gln_L[e]', 0, 'l');
        case '5' % no glucose, glutamine
            R3 = changeRxnBounds(R3, 'EX_glc_D[e]', 0, 'l');
            R3 = changeRxnBounds(R3, 'EX_gln_L[e]', -4, 'l');
        case '6' % no glucose, no glutamine
            R3 = changeRxnBounds(R3, 'EX_glc_D[e]', 0, 'l');
            R3 = changeRxnBounds(R3, 'EX_gln_L[e]', 0, 'l');
        
        otherwise
            fprintf('Undefined medium.\n')
    end
        
    % find all exchange and uptake reactions (including sinks and demand reactions)
    [selExc, selUpt] = findExcRxns(R3);
    
    % block ALL additional uptake except through sink_ reactions
    selUpt_indx = find(selUpt); 
    selUpt_indx_additional = selUpt_indx(~ismember(R3.rxns(selUpt_indx), allowed_IDs)); % filter out our medium components
    selUpt_indx_additional = selUpt_indx_additional(~ismember(R3.rxns(selUpt_indx_additional), {'EX_glc_D[e]','EX_gln_L[e]'})); % filter out glucose and glutamine
    selUpt_indx_additional_non_sinks = selUpt_indx_additional(~contains(R3.rxns(selUpt_indx_additional),'sink_')); % filter out sinks
    R3 = changeRxnBounds(R3, R3.rxns(selUpt_indx_additional_non_sinks), 0, 'l');  
        
    % block (all) amino acid sinks
    DMs = find(contains(R3.rxns,'DM_'));
    sinks = find(contains(R3.rxns,'sink_'));
    aa_sinks = sinks(contains(R3.rxns(sinks),'_L['));
    R3 = changeRxnBounds(R3, R3.rxns(aa_sinks), 0, 'b');
    % this leaves out glyceine
    gly_indx = find(contains(R3.rxns,'sink_gly[c]'));
    R3 = changeRxnBounds(R3, R3.rxns(gly_indx), 0, 'b');
        
    % block inward flux for all sink reactions
    R3 = changeRxnBounds(R3, R3.rxns(sinks), 0, 'l');

    % block sinks and DM_s in all directions
    R3 = changeRxnBounds(R3, R3.rxns(DMs), 0, 'b');
    R3 = changeRxnBounds(R3, R3.rxns(sinks), 0, 'b');
    
    % test growth
    sol = optimizeCbModel(R3,'max');
    fprintf('Growth rate after setting medium bounds for medium #%s: %f.\n',medium, sol.f)
    
    writeCbModel(R3, 'mat',['../Models/medium_only/Recon3DModel_301_patched_M' medium]);
end

%% Write an Excel file 
R3 = readCbModel('../Models/medium_only/Recon3DModel_301_patched_M1.mat');

writeCbModel(R3,'format','xlsx','fileName','../Tables/Recon3DModel_301_patched_medium_1.xlsx')
