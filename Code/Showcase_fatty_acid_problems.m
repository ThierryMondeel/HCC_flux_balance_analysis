%% Settings and setup
format compact

initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

%% Load Recon3D
Recon3DModel = readCbModel('../Models/Recon3DModel_301.mat');

% make sure the set full biomass as objective
Recon3DModel.c = zeros(size(Recon3DModel.c));
Recon3DModel.c(strcmp(Recon3DModel.rxns,'biomass_reaction')) = 1;


%% Biomass flux in default Recon3D version
% This is ~750
sol = optimizeCbModel(Recon3DModel,'max');
fprintf('Growth rate on the vanilla R3 with biomass as the objective: %f.\n',sol.f);

%% Impose our medium and block all other exchanges inward
% Now biomass flux becomes ~368 which is way higher than one would expect.
% This is due to the presence of sources and sinks.

R3 = Recon3DModel; % shorthand

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

% Set glucose and glutamine equivalent to medium 6 in the manuscript
% This is the most limited medium and therefore growth here means it is
% also possible in the other media
R3 = changeRxnBounds(R3, 'EX_glc_D[e]', 0, 'l');
R3 = changeRxnBounds(R3, 'EX_gln_L[e]', 0, 'l');

% find all exchange and uptake reactions (including sinks and demand reactions)
[selExc, selUpt] = findExcRxns(R3);

% block ALL additional uptake except through sink_ reactions
selUpt_indx = find(selUpt); 
selUpt_indx_additional = selUpt_indx(~ismember(R3.rxns(selUpt_indx), allowed_IDs)); % filter out our medium components
selUpt_indx_additional = selUpt_indx_additional(~ismember(R3.rxns(selUpt_indx_additional), {'EX_glc_D[e]','EX_gln_L[e]'})); % filter out glucose and glutamine
selUpt_indx_additional_non_sinks = selUpt_indx_additional(~contains(R3.rxns(selUpt_indx_additional),'sink_')); % filter out sinks
R3 = changeRxnBounds(R3, R3.rxns(selUpt_indx_additional_non_sinks), 0, 'l');  

% test growth
sol = optimizeCbModel(R3,'max','one');
fprintf('Growth rate after setting medium bounds and blocking other exchanges inward: %f.\n', sol.f)

% find active sink and demand reactions and show their flux
% sol.v contains the fluxes
DMs = contains(R3.rxns,'DM_'); % boolean vector 
sinks = contains(R3.rxns,'sink_');
DMs_and_sinks = DMs + sinks; % boolean vector
DMs_and_sinks_fluxes_inw = find(sol.v < 0 & DMs_and_sinks);

table(R3.rxns(DMs_and_sinks_fluxes_inw), sol.v(DMs_and_sinks_fluxes_inw), R3.rxnNames(DMs_and_sinks_fluxes_inw),'VariableNames',{'ID', 'Flux', 'Name'})

%% Impose blocks on demand and amino acid sink reactions

% block (all) amino acid sinks
DMs = find(contains(R3.rxns,'DM_'));
sinks = find(contains(R3.rxns,'sink_'));
aa_sinks = sinks(contains(R3.rxns(sinks),'_L['));
R3 = changeRxnBounds(R3, R3.rxns(aa_sinks), 0, 'b');
% this leaves out glyceine
gly_indx = find(contains(R3.rxns,'sink_gly[c]'));
R3 = changeRxnBounds(R3, R3.rxns(gly_indx), 0, 'b');

% test growth
sol = optimizeCbModel(R3,'max','one');
fprintf('Growth rate after blocking amino acid sinks: %f.\n', sol.f)

% find active sink and demand reactions and show their flux
% sol.v contains the fluxes
DMs = contains(R3.rxns,'DM_'); % boolean vector 
sinks = contains(R3.rxns,'sink_');
DMs_and_sinks = DMs + sinks; % boolean vector
DMs_and_sinks_fluxes_inw = find(sol.v < 0 & DMs_and_sinks);

table(R3.rxns(DMs_and_sinks_fluxes_inw), sol.v(DMs_and_sinks_fluxes_inw), R3.rxnNames(DMs_and_sinks_fluxes_inw),'VariableNames',{'ID', 'Flux', 'Name'})

%% Block all inward sinks

% block inward flux for all sink reactions
R3 = changeRxnBounds(R3, R3.rxns(sinks), 0, 'l');

% test growth
sol = optimizeCbModel(R3,'max','one');
fprintf('Growth rate after blocking sinks in the inward direction: %f.\n', sol.f)

% find active sink and demand reactions and show their flux
% sol.v contains the fluxes
DMs = contains(R3.rxns,'DM_'); % boolean vector 
sinks = contains(R3.rxns,'sink_');
DMs_and_sinks = DMs + sinks; % boolean vector
DMs_and_sinks_fluxes_inw = find(sol.v < 0 & DMs_and_sinks);

table(R3.rxns(DMs_and_sinks_fluxes_inw), sol.v(DMs_and_sinks_fluxes_inw), R3.rxnNames(DMs_and_sinks_fluxes_inw),'VariableNames',{'ID', 'Flux', 'Name'})

%% Show that the triglycerides sink alone solves the problem

% block inward flux for all sink reactions
R3 = changeRxnBounds(R3,'sink_tag_hs[c]', -1000, 'l');

% test growth
sol = optimizeCbModel(R3,'max','one');
fprintf('Growth rate after allowing triglycerides: %f.\n', sol.f)

% find active sink and demand reactions and show their flux
% sol.v contains the fluxes
DMs = contains(R3.rxns,'DM_'); % boolean vector 
sinks = contains(R3.rxns,'sink_');
DMs_and_sinks = DMs + sinks; % boolean vector
DMs_and_sinks_fluxes_inw = find(sol.v < 0 & DMs_and_sinks);

table(R3.rxns(DMs_and_sinks_fluxes_inw), sol.v(DMs_and_sinks_fluxes_inw), R3.rxnNames(DMs_and_sinks_fluxes_inw),'VariableNames',{'ID', 'Flux', 'Name'})

%% Block all sinks and demand reactions in all directions and add Rtotal synthesis
% block sinks and DM_s in all directions
R3 = changeRxnBounds(R3, R3.rxns(DMs), 0, 'b');
R3 = changeRxnBounds(R3, R3.rxns(sinks), 0, 'b');

% Add Rtotal synthesis
[R3, rxnIDexists] = addReaction(R3, 'Rtotal_synth', 'metaboliteList',{'stcoa[c]', 'pmtcoa[c]', 'odecoa[c]', 'lneldccoa[c]', 'Rtotalcoa[c]'},...
    'stoichCoeffList',[-1 -1 -1 -1 4], ...
    'reversible',true);

% test growth
sol = optimizeCbModel(R3,'max','one');
fprintf('Growth rate after blocking sinks and demand reactions and adding Rtotal synthesis: %f.\n', sol.f)

% find active sink and demand reactions and show their flux
% sol.v contains the fluxes
DMs = contains(R3.rxns,'DM_'); % boolean vector 
sinks = contains(R3.rxns,'sink_');
DMs_and_sinks = DMs + sinks; % boolean vector
DMs_and_sinks_fluxes_inw = find(sol.v < 0 & DMs_and_sinks);

table(R3.rxns(DMs_and_sinks_fluxes_inw), sol.v(DMs_and_sinks_fluxes_inw), R3.rxnNames(DMs_and_sinks_fluxes_inw),'VariableNames',{'ID', 'Flux', 'Name'})


%% Combine Rtotal species
% define metabolite IDs that need to change
map_to_Rtotal_c = {'Rtotal2[c]','Rtotal3[c]'};
map_to_Rtotal_e = {'Rtotal2[e]','Rtotal3[e]'};
map_to_Rtotalcoa_c = {'Rtotal2coa[c]','Rtotal3coa[c]'};
change = {'Rtotal3coa[m]','Rtotal3crn[c]','Rtotal3crn[m]'};

% rename some metabolite IDs
for met = change
    indx = findMetIDs(R3,met);
    R3.mets(indx) = strrep(met,'3','');
end

% Replace stoichiometries to new species
for met = map_to_Rtotal_c
    row = findMetIDs(R3,met);
    if row == 0
        break
    end
    newrow = findMetIDs(R3,'Rtotal[c]');
    
    active_reactions = find(R3.S(row,:));
    stoich = R3.S(row, active_reactions);
    R3.S(row, active_reactions) = 0;
    R3.S(newrow, active_reactions) = stoich;
end

for met = map_to_Rtotal_e
    row = findMetIDs(R3,met);
    newrow = findMetIDs(R3,'Rtotal[e]');
    
    active_reactions = find(R3.S(row,:));
    stoich = R3.S(row, active_reactions);
    R3.S(row, active_reactions) = 0;
    R3.S(newrow, active_reactions) = stoich;
end

for met = map_to_Rtotalcoa_c
    row = findMetIDs(R3,met);
    newrow = findMetIDs(R3,'Rtotalcoa[c]');
    
    active_reactions = find(R3.S(row,:));
    stoich = R3.S(row, active_reactions);
    R3.S(row, active_reactions) = 0;
    R3.S(newrow, active_reactions) = stoich;
end

% test growth
sol = optimizeCbModel(R3,'max','one');
fprintf('Growth rate after equation Rtotal species: %f.\n', sol.f)

% find active sink and demand reactions and show their flux
% sol.v contains the fluxes
DMs = contains(R3.rxns,'DM_'); % boolean vector 
sinks = contains(R3.rxns,'sink_');
DMs_and_sinks = DMs + sinks; % boolean vector
DMs_and_sinks_fluxes_inw = find(sol.v < 0 & DMs_and_sinks);

table(R3.rxns(DMs_and_sinks_fluxes_inw), sol.v(DMs_and_sinks_fluxes_inw), R3.rxnNames(DMs_and_sinks_fluxes_inw),'VariableNames',{'ID', 'Flux', 'Name'})
