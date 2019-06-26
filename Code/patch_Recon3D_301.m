%% Settings and setup
format compact

initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

%% Load Recon3D
Recon3DModel = readCbModel('../Models/Recon3DModel_301.mat');

% make sure the set full biomass as objective
Recon3DModel.c = zeros(size(Recon3DModel.c));
Recon3DModel.c(strcmp(Recon3DModel.rxns,'biomass_reaction')) = 1;
sol = optimizeCbModel(Recon3DModel,'max');
fprintf('Growth rate on the vanilla R3 with biomass as the objective: %f.\n',sol.f);

%% Patch Recon3D: combine Rtotal species + add synthesis step
% define metabolite IDs that need to change
map_to_Rtotal_c = {'Rtotal2[c]','Rtotal3[c]'};
map_to_Rtotal_e = {'Rtotal2[e]','Rtotal3[e]'};
map_to_Rtotalcoa_c = {'Rtotal2coa[c]','Rtotal3coa[c]'};
change = {'Rtotal3coa[m]','Rtotal3crn[c]','Rtotal3crn[m]'};

% rename some metabolite IDs
for met = change
    indx = findMetIDs(Recon3DModel,met);
    Recon3DModel.mets(indx) = strrep(met,'3','');
end

% Replace stoichiometries to new species
for met = map_to_Rtotal_c
    row = findMetIDs(Recon3DModel,met);
    if row == 0
        break
    end
    newrow = findMetIDs(Recon3DModel,'Rtotal[c]');
    
    active_reactions = find(Recon3DModel.S(row,:));
    stoich = Recon3DModel.S(row, active_reactions);
    Recon3DModel.S(row, active_reactions) = 0;
    Recon3DModel.S(newrow, active_reactions) = stoich;
end

for met = map_to_Rtotal_e
    row = findMetIDs(Recon3DModel,met);
    newrow = findMetIDs(Recon3DModel,'Rtotal[e]');
    
    active_reactions = find(Recon3DModel.S(row,:));
    stoich = Recon3DModel.S(row, active_reactions);
    Recon3DModel.S(row, active_reactions) = 0;
    Recon3DModel.S(newrow, active_reactions) = stoich;
end

for met = map_to_Rtotalcoa_c
    row = findMetIDs(Recon3DModel,met);
    newrow = findMetIDs(Recon3DModel,'Rtotalcoa[c]');
    
    active_reactions = find(Recon3DModel.S(row,:));
    stoich = Recon3DModel.S(row, active_reactions);
    Recon3DModel.S(row, active_reactions) = 0;
    Recon3DModel.S(newrow, active_reactions) = stoich;
end

% Add Rtotal synthesis
[Recon3DModel, rxnIDexists] = addReaction(Recon3DModel, 'Rtotal_synth', 'metaboliteList',{'stcoa[c]', 'pmtcoa[c]', 'odecoa[c]', 'lneldccoa[c]', 'Rtotalcoa[c]'},...
    'stoichCoeffList',[-1 -1 -1 -1 4], ...
    'reversible',true);

writeCbModel(Recon3DModel,'mat','../Models/Recon3DModel_301_patched.mat')