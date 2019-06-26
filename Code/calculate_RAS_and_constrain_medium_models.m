%% Settings and setup
format compact

initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

%% Read transcriptomics data
% genes were mapped from gene name to Entrez in Python with the myGene module
% df = tdfread('Data/microarray_data_with_entrez_genes.csv'); % Microarray
df = tdfread('../Data/RNAseq_data_with_entrez_genes.csv'); % RNA-seq

% store the cell line columns individually
PLC = df.PLC;
Huh7 = df.Huh7;

datasetnames = {'PLC','Huh7'};
datasets = {PLC, Huh7};

%% Calculate RAS scores for all reactions
R3 = readCbModel('../Models/medium_only/Recon3DModel_301_patched_M1.mat'); % load one of our patched medium models

RASmatrix = zeros(length(R3.rxns),length(datasets));
for i = 1:length(datasets)
    data = datasets{i};
    transcripts = data;
    
    % loop over all reactions to define the RAS
    for j=1:length(R3.rxns)
        if(strcmp(R3.rules(j),'')) % no genes linked to this reaction
            Coeff = NaN;
        else
            rules = R3.rules{j};
            Coeff = getScore(rules, transcripts); 
        end
        RASmatrix(j,i) = Coeff; 
    end
end

% How many reactions will get a restricted bound?
length(find(~isnan(RASmatrix(:,1))))

% export the data for analysis in Python
save('../Data/RAS_scores.mat','RASmatrix','-mat')

%% Plot maximal biomass flux vs. NNR factor on all 12 models
% which NNR factors to try?
scales = logspace(3,-2,50); % logarithmically from 1 to 0.01
final_NNR_index = 40; % was 13
final_NNR = scales(final_NNR_index);

% store results in a matrix
% columns for the NNR factor, rows for the 12 models
% values will be the maximal biomass flux
results = zeros(12,length(scales));

% loop over the two datasets
model_count = 1;
for j = 1:length(datasets)
    % get RAS scores for this cell line
    cell_line = datasetnames{j};
    RAScolumn = RASmatrix(:,j); % PLC, Huh7

    % loop over the 6 medium models
    for medium = {'1','2','3','4','5','6'}
        % load model
        boundedModel_orig = readCbModel(['../Models/medium_only/Recon3DModel_301_patched_M' char(medium) '.mat']);

        % loop over NNR factor
        NNR_count = 1;
        for NNR = scales
            boundedModel = boundedModel_orig;
            maxScore = 1000;

            % loop over all rxns and adjust bounds
            for i=1:length(boundedModel.rxns)
                rxn_val = RAScolumn(i);

                if isnan(rxn_val) % no genes linked to reaction
                    % check if bounds are 1000 (this excludes specifically set
                    % exchange reactions
                    if boundedModel.lb(i) == -1000
                        boundedModel.lb(i) = -1*maxScore;
                    end
                    if boundedModel.ub(i) == 1000
                        boundedModel.ub(i) = maxScore;
                    end
                else % nonzero value
                    if boundedModel.lb(i) == -1000
                        boundedModel.lb(i) = -1*NNR*rxn_val;
                    end
                    if boundedModel.ub(i) == 1000
                        boundedModel.ub(i) = NNR*rxn_val;
                    end
                end
            end

            % store maximal biomass flux
            sol = optimizeCbModel(boundedModel,'max');
            results(model_count, NNR_count) = sol.obj;

            % save models for the final NNR factor
            if NNR == final_NNR 
                writeCbModel(boundedModel, 'mat', ['../Models/Recon3DModel_301_patched_M' char(medium), '_', cell_line, '_NNR=', num2str(NNR),'.mat']);
            end
            
            NNR_count = NNR_count + 1;
            
        end
        
        model_count = model_count + 1;
    end
end

% plot results
% PLC
figure('pos',[10 10 1000 400],'DefaultAxesFontSize',16)
semilogx(scales,results(1:6,:),'--') % PLC
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
semilogx(scales,results(7:12,:),':') % Huh7
xlim([0.01,500])
xlabel('NNS')
ylabel('Maximal biomass flux')
legend('PLC M1','PLC M2','PLC M3','PLC M4','PLC M5','PLC M6',...
        'Huh7 M1','Huh7 M2','Huh7 M3','Huh7 M4','Huh7 M5','Huh7 M6')
line([final_NNR, final_NNR],[0,1.4],'Color','black','LineStyle','--','HandleVisibility','off')
hold off;

fig = gcf; fig.PaperPositionMode = 'auto'; % on-screen size
export_fig('../Figures/biomass_vs_NNR.png','-png','-r200','-p0.01')

% biomass bar plot
figure('pos',[10 10 1000 400],'DefaultAxesFontSize',16)
bar(results(:,final_NNR_index),'k')
xticklabels({'PLC M1','PLC M2','PLC M3','PLC M4','PLC M5','PLC M6',...
                'Huh7 M1','Huh7 M2','Huh7 M3','Huh7 M4','Huh7 M5','Huh7 M6'})
xtickangle(45)
ylabel('Maximal biomass flux')
export_fig('../Figures/biomass_at_final_NNR.png','-png','-r200','-p0.01')

%% Compare various objectives across all models
objectives = {'EX_glc_D[e]','EX_gln_L[e]','EX_o2[e]','EX_phe_L[e]','EX_lac_L[e]','EX_pyr[e]','ATPM','EX_co2[e]'};

% store results in a matrix: rows for models, columns for objectives
results_min = zeros(12, length(objectives));
results_max = zeros(12, length(objectives));
model_count = 1;
for cell_line = {'PLC','Huh7'}
    for medium = {'1','2','3','4','5','6'}
        temp_model_orig = readCbModel(['../Models/Recon3DModel_301_patched_M' char(medium), '_', char(cell_line), '_NNR=', num2str(scales(final_NNR_index)),'.mat']);                    
                                
        r_count = 1;
        for r = objectives

            if strcmp(r,'ATPM')
               [temp_model, rxnIDexists] = addReaction(temp_model_orig, 'ATPM', 'metaboliteList',{'adp[c]', 'pi[c]', 'h2o[c]','h[c]', 'atp[c]'},...
                                    'stoichCoeffList',[1 1 -1 1 -1], ...
                                    'reversible',false); 
            else
                temp_model = temp_model_orig;
            end
            
            if strcmp(r,'EX_co2[e]')
                % get maximal biomass flux
                sol = optimizeCbModel(temp_model,'max');
                
                % enfore the maximal biomass flux
                temp_model.lb(strcmp(temp_model.rxns,'biomass_reaction')) = sol.f;
                
                % minimize and maximize the sum of CO2 and bicarbonate exchange
                temp_model.c = zeros(size(temp_model.c));
                temp_model.c(strcmp(temp_model.rxns,'EX_co2[e]')) = 1;
                temp_model.c(strcmp(temp_model.rxns,'EX_hco3[e]')) = 1;
                sol = optimizeCbModel(temp_model,'min');
                minFlux = sol.f;
                sol = optimizeCbModel(temp_model,'max');
                maxFlux = sol.f;
            else
                opt = 100; % force maximal biomass flux
                [minFlux, maxFlux] = fluxVariability(temp_model, opt, 'max', r);
            end

            results_min(model_count, r_count) = minFlux;
            results_max(model_count, r_count) = maxFlux;

            r_count = r_count + 1;
        end

        model_count = model_count + 1;
    end
end

% look at results for things that are close to boundaries
results_min
results_max

% plot
set(0,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
titles = {'Glucose','Glutamine','Oxygen','Phenylalanine','Lactate','Pyruvate','ATP','CO_2 + bicarbonate'};
figure('pos',[10 10 1400 1000],'DefaultAxesFontSize',16)
plot_count = 1;
for i = 1:length(objectives) % loop over models
    subplot(4,2,plot_count)
    plot([1:12; 1:12],[(results_min(:,i))'; results_max(:,i)'], 'LineWidth',5,'Color','k')
    hold on;
    for j = 1:12 % loop over models for current objective i and if range is small show a marker
        if abs(results_min(j,i) - results_max(j,i)) < 1e-1
            plot(j, results_max(j,i),'*k')
        end
    end
    xlim([0,13]) % for extra space on the RHS
    yrange = max(results_max(:,i)) - min(results_min(:,i));
    ylim([min(0, min(results_min(:,i))) max(0, max(results_max(:,i)))])
    xticks(1:12)
    xticklabels({'PLC M1','PLC M2','PLC M3','PLC M4','PLC M5','PLC M6',...
                'Huh7 M1','Huh7 M2','Huh7 M3','Huh7 M4','Huh7 M5','Huh7 M6'})
    xtickangle(45)
    ylabel('Flux [mM/h]')
    title(titles(i))
    
    plot_count = plot_count + 1;
end

fig = gcf; fig.PaperPositionMode = 'auto'; % on-screen size
export_fig(strcat('../Figures/FVA_input_output.png'),'-png','-r200','-p0.01')

%% WITHOUT CHOLESTEROL: NO OXYGEN NEEDED?
objectives = {'EX_o2[e]'};

% store results in a matrix: rows for models, columns for objectives
results_min = zeros(12, length(objectives));
results_max = zeros(12, length(objectives));
model_count = 1;
for cell_line = {'PLC','Huh7'}
    for medium = {'1','2','3','4','5','6'}
        temp_model_orig = readCbModel(['../Models/Recon3DModel_301_patched_M' char(medium), '_', char(cell_line), '_NNR=', num2str(scales(final_NNR_index)),'.mat']);                    
                       
        % set cholesterol in biomass reaction to zero
        r_indx = find(strcmp(temp_model.rxns,'biomass_reaction'));
        met_indx = find(strcmp(temp_model.mets,'chsterol[c]'));
        temp_model_orig.S(met_indx, r_indx) = 0;
        
        r_count = 1;
        for r = objectives

            temp_model = temp_model_orig;

            opt = 100; % force maximal biomass flux
            [minFlux, maxFlux] = fluxVariability(temp_model, opt, 'max', r);

            results_min(model_count, r_count) = minFlux;
            results_max(model_count, r_count) = maxFlux;

            r_count = r_count + 1;
        end

        model_count = model_count + 1;
    end
end

% look at results for things that are close to boundaries
results_min
results_max

%% Respiration reactions at maximal growth
objectives = {'ATPS4mi','NADH2_u10mi','CYOR_u10mi','CYOOm2i','CYOOm3i','PDHm'};

% store results in a matrix: rows for models, columns for objectives
% store results in a matrix: rows for models, columns for objectives
results_min = zeros(12, length(objectives));
results_max = zeros(12, length(objectives));
model_count = 1;
for cell_line = {'PLC','Huh7'}
    for medium = {'1','2','3','4','5','6'}
        temp_model_orig = readCbModel(['../Models/Recon3DModel_301_patched_M' char(medium), '_', char(cell_line), '_NNR=', num2str(scales(final_NNR_index)),'.mat']);                                            
       
        % without forcing biomass
        r_count = 1;
        for r = objectives
            
            temp_model = temp_model_orig;
            
            opt = 100; % force maximal biomass flux
            [minFlux, maxFlux] = fluxVariability(temp_model, opt, 'max', r);

            results_min(model_count, r_count) = minFlux;
            results_max(model_count, r_count) = maxFlux;

            r_count = r_count + 1;
        end

        model_count = model_count + 1;
    end
end

% look at results for things that are close to boundaries
results_min
results_max

% plot
set(0,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
titles = {'ATP synthase (ATPS4mi)','NADH dehydrogenase (NADH2\_u10mi)','ubiquinol-cytochrome c reductase (CYOR\_u10mi)',...
    'cytochrome c oxidase (CYOOm2i)','cytochrome c oxidase (CYOOm3i)','Pyruvate dehydrogenase (PDHm)'};
figure('pos',[10 10 1400 600],'DefaultAxesFontSize',16)
plot_count = 1;
for i = 1:length(objectives) % loop over models
    subplot(3,2,plot_count)
    plot([1:12; 1:12],[(results_min(:,i))'; results_max(:,i)'], 'LineWidth',5,'Color','k')
    hold on;
    for j = 1:12 % loop over models for current objective i and if range is small show a marker
        if abs(results_min(j,i) - results_max(j,i)) < 1e-2
            plot(j, results_max(j,i),'*k')
        end
    end
    xlim([0,13]) % for extra space on the RHS
    yrange = max(results_max(:,i)) - min(results_min(:,i));
    ylim([min(0, min(results_min(:,i))) max(0, max(results_max(:,i)))])
    xticks(1:12)
    xticklabels({'PLC M1','PLC M2','PLC M3','PLC M4','PLC M5','PLC M6',...
                'Huh7 M1','Huh7 M2','Huh7 M3','Huh7 M4','Huh7 M5','Huh7 M6'})
    xtickangle(45)
    ylabel('Flux [mM/h]')
    title(titles(i))
    
    plot_count = plot_count + 1;
end

fig = gcf; fig.PaperPositionMode = 'auto'; % on-screen size
export_fig(strcat('../Figures/FVA_respiration.png'),'-png','-r200','-p0.01')

%% Write all models to SBML
for cell_line = {'PLC','Huh7'}
    for medium = {'1','2','3','4','5','6'}
        modelName = ['../Models/Recon3DModel_301_patched_M' char(medium), '_', char(cell_line), '_NNR=', num2str(scales(final_NNR_index)),'.mat'];
        m = readCbModel(modelName);
    writeCbModel(m, 'sbml', modelName);
    end
end