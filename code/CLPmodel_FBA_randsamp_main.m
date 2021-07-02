initCobraToolbox
setRavenSolver('cobra')
changeCobraSolver('gurobi','LP')

% input:
% human-GEM v1.6.0
% animal GEMs
% "spp_GEMs_for_sims.csv" containing GEMs with data
% "exchange_input_human.csv" containing open exchange reactions 
% "Simian_genera.csv" of higher primates (Simians)
% function "getCLPmodel"
% function "randomSampling_vertGEM'
% 
% output is: a collection of CLP models (animalCLPmodels)
% random sampling results
% CLP model max ATP


%% load GEMs
load('Human-GEM_1_6_0.mat')
ihuman.id = '9606';

% cd to folder containing full animal GEMs
files = dir(fullfile('./', '*.mat'));

% Load GEM list
fid = fopen('spp_GEMs_for_sims.csv');
dataFile = textscan(fid,'%s %f32 %s %s %f32 %f32 %f32 %f32','Delimiter',',','HeaderLines',1);
GEMlist_SciName = dataFile{1};
GEMlist_taxonID = dataFile{2};
GEMlist_category = dataFile{4};
GEMlist_TL = dataFile{5};
GEMlist_percent_dietprot = dataFile{6};
GEMlist_percent_dietcarb = dataFile{7};
GEMlist_percent_dietfat = dataFile{8};
fclose(fid);


animalGEMs = {};
animalGEMs = [animalGEMs,ihuman];

% load all GEMs in list & collect into structure 
for j = 2:length(GEMlist_taxonID) %1 is ihuman
    filename =  [num2str(GEMlist_taxonID(j)) '.mat'];
    load(filename);  
    animalGEMs = [animalGEMs,model];
end

animalGEMs = animalGEMs(:);
disp('animalGEMs loaded')

%% get CLPmodels

% Load constraints file (human CLP)
fid = fopen('exchange_input_human.csv');
dataFile        = textscan(fid,'%s %f32 %f32','Delimiter',',','HeaderLines',1);
exchange_mets   = dataFile{1};
exchange_lb     = dataFile{2};
exchange_ub     = dataFile{3};
fclose(fid);

% Load list of Simians
fid = fopen('Simian_genera.csv');
dataFile        = textscan(fid,'%s');
Simians         = dataFile{1};
fclose(fid);

animalCLPmodels = {};
GEMlist_IDx_converted=[];
GEMlist_IDx_failed=[];



for i = 1:length(animalGEMs)
   disp(['converting ' num2str(i) ' of ' num2str(length(animalGEMs)) ' GEMs'])
   model = animalGEMs{i}; 
   open_exchange_mets = exchange_mets;
   open_exchange_lb = exchange_lb;
   open_exchange_ub = exchange_ub;
   
   %diet components
   k = find(ismember(open_exchange_mets, 'glucose'));
   open_exchange_lb(k) = -GEMlist_percent_dietcarb(i) / 180.16 * 1000;
   open_exchange_ub(k) = -GEMlist_percent_dietcarb(i) / 180.16 * 1000;
   
   k = find(ismember(open_exchange_mets, 'fatty acid-uptake pool'));
   open_exchange_lb(k) = -GEMlist_percent_dietfat(i) / 280.63 * 1000;
   open_exchange_ub(k) = -GEMlist_percent_dietfat(i) / 280.63 * 1000;
   
   k = find(ismember(open_exchange_mets, 'protein-uptake pool'));
   open_exchange_lb(k) = -GEMlist_percent_dietprot(i) / 111.06 * 1000;
   open_exchange_ub(k) = -GEMlist_percent_dietprot(i) / 111.06 * 1000;
   
   %nitrogenous waste
   if ismember(GEMlist_category(i), 'Fishes')
       open_exchange_lb(find(ismember(open_exchange_mets, 'urea'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'urea'))) = 1000;
       open_exchange_lb(find(ismember(open_exchange_mets, 'NH3'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'NH3'))) = 1000;
       
       open_exchange_lb(find(ismember(open_exchange_mets, 'urate'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'urate'))) = 0;
       open_exchange_lb(find(ismember(open_exchange_mets, 'allantoin'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'allantoin'))) = 0;
   end
   
   if ismember(GEMlist_category(i), 'Birds') | ismember(GEMlist_category(i), 'Reptiles')
       open_exchange_lb(find(ismember(open_exchange_mets, 'urate'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'urate'))) = 1000;
       
       open_exchange_lb(find(ismember(open_exchange_mets, 'NH3'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'NH3'))) = 0;
       open_exchange_lb(find(ismember(open_exchange_mets, 'urea'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'urea'))) = 0;
       open_exchange_lb(find(ismember(open_exchange_mets, 'allantoin'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'allantoin'))) = 0;
   end
   
   if ismember(GEMlist_category(i), 'Cartilaginous fishes') | ismember(GEMlist_category(i), 'Amphibians')
       open_exchange_lb(find(ismember(open_exchange_mets, 'urea'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'urea'))) = 1000;
       
       open_exchange_lb(find(ismember(open_exchange_mets, 'NH3'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'NH3'))) = 0;
       open_exchange_lb(find(ismember(open_exchange_mets, 'urate'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'urate'))) = 0;
       open_exchange_lb(find(ismember(open_exchange_mets, 'allantoin'))) = 0;
       open_exchange_ub(find(ismember(open_exchange_mets, 'allantoin'))) = 0;
   end
   
   if ismember(GEMlist_category(i), 'Mammals')
       mamm_genus = strsplit(GEMlist_SciName{i}, ' ');
       mamm_genus = mamm_genus(1);
       if sum(ismember(Simians, mamm_genus)) == 0
           open_exchange_lb(find(ismember(open_exchange_mets, 'urea'))) = 0;
           open_exchange_ub(find(ismember(open_exchange_mets, 'urea'))) = 1000;
           open_exchange_lb(find(ismember(open_exchange_mets, 'allantoin'))) = 0;
           open_exchange_ub(find(ismember(open_exchange_mets, 'allantoin'))) = 1000;
           
           open_exchange_lb(find(ismember(open_exchange_mets, 'NH3'))) = 0;
           open_exchange_ub(find(ismember(open_exchange_mets, 'NH3'))) = 0;
           open_exchange_lb(find(ismember(open_exchange_mets, 'urate'))) = 0;
           open_exchange_ub(find(ismember(open_exchange_mets, 'urate'))) = 0;
       end
   end
   
   % fix Opossum (13616) 
   if GEMlist_taxonID(i) == 13616
       Opossum_addback = {'HMR_3881';'HMR_4216';'HMR_4184'};

        for Opossum_i = 1:length(Opossum_addback) %
            %grab metabolites in the rxn
            query_rxn = Opossum_addback{Opossum_i}; %
            ihuman_rxn_idx = find(ismember(ihuman.rxns, query_rxn));
            query_met_idx = find(ihuman.S(:,ihuman_rxn_idx));

            %grab metabolite IDs: first check that all metabolites are there, if not add it
            model_met_idx = [];
            for Opossum_j = 1:length(query_met_idx)
               query_met_j = ihuman.mets{query_met_idx(Opossum_j)}; 
               query_metName_j = ihuman.metNames{query_met_idx(Opossum_j)}; 
               model_met_j = find(ismember(model.mets,query_met_j));
               if isempty(model_met_j)
                   model.mets(end+1) = {query_met_j};
                   model.metNames(end+1) = {query_metName_j};
                   model.b(end+1) = ihuman.b(query_met_idx(Opossum_j));
                   model_met_idx = [model_met_idx; length(model.mets)];
               else
                   model_met_idx = [model_met_idx; model_met_j];
               end
            end

            %add the reaction
            model.S(model_met_idx,end+1) = ihuman.S(query_met_idx,ihuman_rxn_idx);
            model.rxns(end+1) = {query_rxn};
            model.ub(end+1) = ihuman.ub(ihuman_rxn_idx);
            model.lb(end+1) = ihuman.lb(ihuman_rxn_idx);
            model.c(end+1) = ihuman.c(ihuman_rxn_idx);

            model.rev(end+1) = ihuman.rev(ihuman_rxn_idx);
            model.rxnNames(end+1) = ihuman.rxnNames(ihuman_rxn_idx);
            model.grRules(end+1) = ihuman.grRules(ihuman_rxn_idx);

            rxn_gene_idx = find(ihuman.rxnGeneMat(ihuman_rxn_idx,:));
            rxn_gene = ihuman.genes(rxn_gene_idx);
            rxnGeneMat_idx = find(ismember(model.genes,rxn_gene));

            if isempty(rxnGeneMat_idx)
                rxnGeneMat_idx = length(model.genes)+1;
                model.genes(rxnGeneMat_idx) = ihuman.genes(rxn_gene_idx);
            end

            model.rxnGeneMat(end+1,rxnGeneMat_idx) = 1;
            model.subSystems(end+1) = ihuman.subSystems(ihuman_rxn_idx);
            model.eccodes(end+1) = ihuman.eccodes(ihuman_rxn_idx);
            model.rxnReferences(end+1) = ihuman.rxnReferences(ihuman_rxn_idx);
            model.rxnConfidenceScores(end+1) = ihuman.rxnConfidenceScores(ihuman_rxn_idx);

        end

   end
   
   % fix Tasmanian Devil (9305)
   if GEMlist_taxonID(i) == 9305
       TasDevil_addback = {'HMR_4216';'HMR_4406';'HMR_8022';'HMR_4184'};

        for TasDevil_i = 1:length(TasDevil_addback) %
            %grab metabolites in the rxn
            query_rxn = TasDevil_addback{TasDevil_i}; %
            ihuman_rxn_idx = find(ismember(ihuman.rxns, query_rxn));
            query_met_idx = find(ihuman.S(:,ihuman_rxn_idx));

            %grab metabolite IDs: first check that all metabolites are there, if not add it
            model_met_idx = [];
            for TasDevil_j = 1:length(query_met_idx)
               query_met_j = ihuman.mets{query_met_idx(TasDevil_j)}; 
               query_metName_j = ihuman.metNames{query_met_idx(TasDevil_j)}; 
               model_met_j = find(ismember(model.mets,query_met_j));
               if isempty(model_met_j)
                   model.mets(end+1) = {query_met_j};
                   model.metNames(end+1) = {query_metName_j};
                   model.b(end+1) = ihuman.b(query_met_idx(TasDevil_j));
                   model_met_idx = [model_met_idx; length(model.mets)];
               else
                   model_met_idx = [model_met_idx; model_met_j];
               end
            end

            %add the reaction
            model.S(model_met_idx,end+1) = ihuman.S(query_met_idx,ihuman_rxn_idx);
            model.rxns(end+1) = {query_rxn};
            model.ub(end+1) = ihuman.ub(ihuman_rxn_idx);
            model.lb(end+1) = ihuman.lb(ihuman_rxn_idx);
            model.c(end+1) = ihuman.c(ihuman_rxn_idx);

            model.rev(end+1) = ihuman.rev(ihuman_rxn_idx);
            model.rxnNames(end+1) = ihuman.rxnNames(ihuman_rxn_idx);
            model.grRules(end+1) = ihuman.grRules(ihuman_rxn_idx);

            rxn_gene_idx = find(ihuman.rxnGeneMat(ihuman_rxn_idx,:));
            rxn_gene = ihuman.genes(rxn_gene_idx);
            rxnGeneMat_idx_all = [];
            for TasDevil_k = 1:length(rxn_gene)
                rxnGeneMat_idx = find(ismember(model.genes,rxn_gene(TasDevil_k)));
                if isempty(rxnGeneMat_idx)
                    rxnGeneMat_idx = length(model.genes)+1;
                    model.genes(rxnGeneMat_idx) = rxn_gene(TasDevil_k);
                end
                rxnGeneMat_idx_all = [rxnGeneMat_idx_all;rxnGeneMat_idx];
            end
            model.rxnGeneMat(end+1,rxnGeneMat_idx_all) = 1;
            model.subSystems(end+1) = ihuman.subSystems(ihuman_rxn_idx);
            model.eccodes(end+1) = ihuman.eccodes(ihuman_rxn_idx);
            model.rxnReferences(end+1) = ihuman.rxnReferences(ihuman_rxn_idx);
            model.rxnConfidenceScores(end+1) = ihuman.rxnConfidenceScores(ihuman_rxn_idx);

        end
   end
   
   %make all model ID = taxonID
   model.id = GEMlist_taxonID(i);
   
   % main 
   try
      CLPmodel = getCLPmodel(model, open_exchange_mets, open_exchange_lb, open_exchange_ub);
      animalCLPmodels = [animalCLPmodels,CLPmodel]; 
      GEMlist_IDx_converted = [GEMlist_IDx_converted;i];
   catch
      disp(['GEM number ' num2str(i) ' failed to convert']) 
      GEMlist_IDx_failed = [GEMlist_IDx_failed;i];
   end
end

animalCLPmodels = animalCLPmodels(:);

disp('CLPmodels constructed')


% for the constructed CLPmodels, move on to calculate max ATP


%% max ATP

CLPmodel_max_ATP=[];

for i = 1:length(animalCLPmodels)
   model = animalCLPmodels{i};
   model.c(:) = 0;
   model.c(find(ismember(model.rxns, 'EX_atp[e]'))) = 1; 
   max_ATP_sol = optimizeCbModel(model);
   CLPmodel_max_ATP = [CLPmodel_max_ATP;max_ATP_sol.f];
end

max_ATP_results = table(GEMlist_taxonID(GEMlist_IDx_converted),CLPmodel_max_ATP);


disp('max ATP done')


%% max ATP with AMP allowance


CLPmodel_maxATP_AMP2g=[];
CLPmodel_maxATP_AMP5g=[];
CLPmodel_maxATP_AMP10g=[];
CLPmodel_maxATP_AMPunlimited=[];

for i = 1:length(animalCLPmodels)
    disp(i);
   model = animalCLPmodels{i};
   model.c(:) = 0;
   model.c(find(ismember(model.rxns, 'EX_atp[e]'))) = 1; 
   
   model.ub(find(ismember(model.rxns, 'HMR_9262'))) = 0;
   model.lb(find(ismember(model.rxns, 'HMR_9262'))) = -5.76;
   max_ATP_sol = optimizeCbModel(model);
   CLPmodel_maxATP_AMP2g = [CLPmodel_maxATP_AMP2g;max_ATP_sol.f];
   
   model.ub(find(ismember(model.rxns, 'HMR_9262'))) = 0;
   model.lb(find(ismember(model.rxns, 'HMR_9262'))) = -14.4;
   max_ATP_sol = optimizeCbModel(model);
   CLPmodel_maxATP_AMP5g = [CLPmodel_maxATP_AMP5g;max_ATP_sol.f];
   
   model.ub(find(ismember(model.rxns, 'HMR_9262'))) = 0;
   model.lb(find(ismember(model.rxns, 'HMR_9262'))) = -28.8;
   max_ATP_sol = optimizeCbModel(model);
   CLPmodel_maxATP_AMP10g = [CLPmodel_maxATP_AMP10g;max_ATP_sol.f];
   
   model.ub(find(ismember(model.rxns, 'HMR_9262'))) = 0;
   model.lb(find(ismember(model.rxns, 'HMR_9262'))) = -1000;
   max_ATP_sol = optimizeCbModel(model);
   CLPmodel_maxATP_AMPunlimited = [CLPmodel_maxATP_AMPunlimited;max_ATP_sol.f];
end

maxATP_with_AMP_allowance = table(GEMlist_taxonID(GEMlist_IDx_converted), ...
    CLPmodel_max_ATP, CLPmodel_maxATP_AMP2g, ...
    CLPmodel_maxATP_AMP5g, CLPmodel_maxATP_AMP10g, CLPmodel_maxATP_AMPunlimited);
writetable(maxATP_with_AMP_allowance, 'CLPmodel_maxATP_with_AMP_allowance');


disp('max ATP with AMP allowance done')



%% random sampling:

for i = 1:length(animalCLPmodels)
    disp(['random sampling ' num2str(i) ' of ' num2str(length(animalCLPmodels)) ' GEMs'])
    model = animalCLPmodels{i}; 
    %constrain ATP production to maxATP
    %this asks: during the production of maxATP (from 100g of food), what are
    %the average fluxes
    model.ub(find(ismember(model.rxns, 'EX_atp[e]'))) = CLPmodel_max_ATP(i);
    model.lb(find(ismember(model.rxns, 'EX_atp[e]'))) = CLPmodel_max_ATP(i);
    
    randomSample_sol = randomSampling_vertGEM(model,1000,true, true, false);
    
    randomSample_mean = mean(randomSample_sol,2); %mean of rows
    randomSample_rxnIDx = find(randomSample_mean~=0);
    randomSample_flux = randomSample_mean(randomSample_rxnIDx);
    %randomSample_subsystem = model.subSystems(randomSample_rxnIDx);
    randomSample_rxns = model.rxns(randomSample_rxnIDx);
    
    randomSample_result_table = table(randomSample_rxnIDx,randomSample_rxns,randomSample_flux);
    
    filename = ['randomSample_maxATPconstrained_' num2str(GEMlist_taxonID(GEMlist_IDx_converted(i)))];
    writetable(randomSample_result_table, filename);
end

disp('random sampling done')







