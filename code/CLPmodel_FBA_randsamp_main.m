initCobraToolbox
setRavenSolver('cobra')
changeCobraSolver('gurobi','LP')

% input requirements:
% human-GEM v1.6.0 
% a folder containing all full animal GEMs (30 plus human = 31 total)
% "spp_GEMs_for_sims_ver2.csv" containing GEMs with data
% "exchange_input_human.csv" containing open exchange reactions 
% "Simian_genera.csv" of higher primates (Simians)
% "CaseStudy_61622_monkey_Li_2006.csv", Golden Snub-Nosed Monkey (TaxonID: 61622) case study
% "CaseStudy_9627_fox_Needham_2014.csv", Red Fox (TaxonID: 9627) case study
% function "getCLPmodel"
% function "randomSampling_vertGEM'
% 
% output is: a collection of CLP models (animalCLPmodels)
% random sampling results
% CLP model max ATP


%% load GEMs

load('Human-GEM_1_6_0.mat')
ihuman.id = '9606';

% Load GEM list
fid = fopen('spp_GEMs_for_sims_ver2.csv');
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

%% get CLP models

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

writetable(max_ATP_results, 'max_ATP_results');

disp('max ATP done')




%% random sampling:

for i = 1:length(animalCLPmodels)
    disp(['random sampling ' num2str(i) ' of ' num2str(length(animalCLPmodels)) ' GEMs'])
    model = animalCLPmodels{i}; 
    %constrain ub of ATP production to maxATP, lb to 90% of maxATP
	%do not calculate mean
    model.ub(find(ismember(model.rxns, 'EX_atp[e]'))) = CLPmodel_max_ATP(i);
    model.lb(find(ismember(model.rxns, 'EX_atp[e]'))) = 0.9 * CLPmodel_max_ATP(i);
    
    randomSample_sol = randomSampling_vertGEM(model,1000,true, true, false);
    randomSample_rxns = model.rxns;
    randomSample_result_table = table(randomSample_rxns,randomSample_sol);
    
    filename = ['randomSample_maxATPconstrained_' num2str(GEMlist_taxonID(GEMlist_IDx_converted(i)))];
    writetable(randomSample_result_table, filename);
end

disp('random sampling done')





%% case study: Golden Snub-Nosed Monkey (TaxonID: 61622) 
% data from: https://onlinelibrary.wiley.com/doi/epdf/10.1002/ajp.20220

% Load seasonal diet data (Fig 2)
fid = fopen('CaseStudy_61622_monkey_Li_2006.csv');
dataFile        = textscan(fid,'%f32 %s %f32 %f32 %f32 %f32 %f32 %f32 %f32',...
    'Delimiter',',','HeaderLines',1);
monkey_month     = dataFile{2};
monkey_dietcarb     = dataFile{7};
monkey_dietfat     = dataFile{8};
monkey_dietprot     = dataFile{9};
fclose(fid);

model = animalCLPmodels{28}; 

if model.id ~= 61622
    EM='Wrong model! Get model for golden-snub-nosed monkey, TaxonID 61622';
    dispEM(EM);
end

%sanity check: re-calculate max ATP 
monkey_ATP=[];
model.c(:) = 0;
model.c(find(ismember(model.rxns, 'EX_atp[e]'))) = 1; 
max_ATP_sol = optimizeCbModel(model);
%max_ATP_sol.f;


%find all exchange rxns 
all_exchange_rxns = getExchangeRxns(model); %this is rxn IDs
all_exchange_rxns = find(contains(model.rxns, all_exchange_rxns)); %rxn indices
all_exchange_metIDx = zeros(length(all_exchange_rxns),1);
for i = 1:length(all_exchange_rxns)
	ex_met = find(model.S(:,all_exchange_rxns(i)));
	all_exchange_metIDx(i,1) = ex_met;
end
all_exchange_mets = model.metNames(all_exchange_metIDx);

%find dietcarb, fat, and prot exchange rxns
dietcarb = all_exchange_rxns(find(ismember(all_exchange_mets, 'glucose')));
dietfat = all_exchange_rxns(find(ismember(all_exchange_mets, 'fatty acid-uptake pool')));
dietprot = all_exchange_rxns(find(ismember(all_exchange_mets, 'protein-uptake pool')));


for i = 1:length(monkey_month)
    model.ub(dietcarb) = -monkey_dietcarb(i) /180.16 * 1000;
    model.lb(dietcarb) = -monkey_dietcarb(i) /180.16 * 1000;
    model.ub(dietfat) = -monkey_dietfat(i) /280.63 * 1000;
    model.lb(dietfat) = -monkey_dietfat(i) /280.63 * 1000;
    model.ub(dietprot) = -monkey_dietprot(i) /111.06 * 1000;
    model.lb(dietprot) = -monkey_dietprot(i) /111.06 * 1000;
    
    max_ATP_sol = optimizeCbModel(model);
    monkey_ATP = [monkey_ATP;max_ATP_sol.f];
    
end


monkey_results_table = table(monkey_month,monkey_ATP);
writetable(monkey_results_table, 'casestudy_monkey_results');

disp('monkey case study done')



%% case study: Red Fox (TaxonID: 9627)
% data from: https://link.springer.com/article/10.1007/s13364-014-0188-7

% Load seasonal diet data (Fig 2)
fid = fopen('CaseStudy_9627_fox_Needham_2014.csv');
dataFile        = textscan(fid,'%s %f32 %f32 %f32 %f32 %f32 %f32 %f32',...
    'Delimiter',',','HeaderLines',1);
fox_season     = dataFile{1};
fox_dietcarb     = dataFile{6};
fox_dietfat     = dataFile{7};
fox_dietprot     = dataFile{8};
fclose(fid);

model = animalCLPmodels{8}; 

if model.id ~= 9627
    EM='Wrong model! Get model for red fox, TaxonID 9627';
    dispEM(EM);
end

%sanity check: re-calculate max ATP
fox_ATP=[];
model.c(:) = 0;
model.c(find(ismember(model.rxns, 'EX_atp[e]'))) = 1; 
max_ATP_sol = optimizeCbModel(model);
%max_ATP_sol.f

%find all exchange rxns 
all_exchange_rxns = getExchangeRxns(model); %this is rxn IDs
all_exchange_rxns = find(contains(model.rxns, all_exchange_rxns)); %rxn indices
all_exchange_metIDx = zeros(length(all_exchange_rxns),1);
for i = 1:length(all_exchange_rxns)
	ex_met = find(model.S(:,all_exchange_rxns(i)));
	all_exchange_metIDx(i,1) = ex_met;
end
all_exchange_mets = model.metNames(all_exchange_metIDx);

%find dietcarb, fat, and prot exchange rxns
dietcarb = all_exchange_rxns(find(ismember(all_exchange_mets, 'glucose')));
dietfat = all_exchange_rxns(find(ismember(all_exchange_mets, 'fatty acid-uptake pool')));
dietprot = all_exchange_rxns(find(ismember(all_exchange_mets, 'protein-uptake pool')));


for i = 1:length(fox_season)
    model.ub(dietcarb) = -fox_dietcarb(i) /180.16 * 1000;
    model.lb(dietcarb) = -fox_dietcarb(i) /180.16 * 1000;
    model.ub(dietfat) = -fox_dietfat(i) /280.63 * 1000;
    model.lb(dietfat) = -fox_dietfat(i) /280.63 * 1000;
    model.ub(dietprot) = -fox_dietprot(i) /111.06 * 1000;
    model.lb(dietprot) = -fox_dietprot(i) /111.06 * 1000;
    
    max_ATP_sol = optimizeCbModel(model);
    fox_ATP = [fox_ATP;max_ATP_sol.f];
    
end


fox_results_table = table(fox_season,fox_ATP);
writetable(fox_results_table, 'casestudy_fox_results');

disp('fox case study done')












