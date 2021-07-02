function CLPmodel = getCLPmodel(model, open_exchange_mets, open_exchange_lb, open_exchange_ub)
% getCLPmodel
% CLPmodel allows exchange of diet as CLP (carb, lipid, and protein) compositions
%
% If open_exchange_mets, open_exchange_lb, and open_exchange_ub are given, output is 
% constrained by these bounds, AND reactions that carry zero flux when maximizing growth
% and when maximizing ATP production are constrained to lb==ub==0. This CLPmodel is 
% then ready for random sampling.
%
%   model               a model structure
%   open_exchange_mets  exchange reaction indices that should stay open
%   open_exchange_lb    lb of exchange reaction indices that should stay open
%   open_exchange_ub    ub of exchange reaction indices that should stay open
%
% Rosemary Yu, 2021-04-26

if nargin<2
    open_exchange_mets = [];
	open_exchange_lb = [];
	open_exchange_ub = [];
end

%Delete model.rxnFrom and model.metFrom (only in HumanGEM):
if isfield(model,'rxnFrom')
    model = rmfield(model,'rxnFrom');
end

if isfield(model,'metFrom')
    model = rmfield(model,'metFrom');
end

%Check that the model is feasible as it is given
sol=optimizeCbModel(model);
if isempty(sol.x)
    EM='The model has no feasible solution';
    dispEM(EM);
end

%step 1: add biomass_CLP reaction: first copy the existing biomass reaction
biomass = find(contains(model.rxns, 'biomass'));
biomass = biomass(end);
model.S(:,end+1) = model.S(:,biomass);
model.rxns(end+1) = {'biomass_CLP'};
model.lb(end+1) = 0;
model.ub(end+1) = 1000;
model.rev(end+1) = 0;
model.c(end+1) = 1;
model.rxnNames(end+1) = {'Generic biomass reaction for CLP model'};
model.grRules(end+1) = model.grRules(biomass);
model.rxnGeneMat(end+1,:) = model.rxnGeneMat(biomass,:);
model.subSystems(end+1) = model.subSystems(biomass);
model.eccodes(end+1) = model.eccodes(biomass);
model.rxnReferences(end+1) = model.rxnReferences(biomass);
model.rxnConfidenceScores(end+1) = model.rxnConfidenceScores(biomass);

%remove DNA, RNA, cofactor pool, and metabolite pool from biomass_CLP
biomass_CLP = find(contains(model.rxns, 'biomass_CLP'));
biomass_DNA = find(ismember(model.metNames,'DNA')); %there are two
biomass_DNA = biomass_DNA(find(ismember(model.metComps(biomass_DNA), 8))); %metComps 8 is nuclear
biomass_RNA = find(ismember(model.metNames,'RNA'));
biomass_cof = find(ismember(model.metNames,'cofactor_pool_biomass'));
biomass_met = find(ismember(model.metNames,'metabolite_pool_biomass'));

model.S(biomass_DNA, biomass_CLP) = 0;
model.S(biomass_RNA, biomass_CLP) = 0;
model.S(biomass_cof, biomass_CLP) = 0;
model.S(biomass_met, biomass_CLP) = 0;

%step 2: add protein-uptake pool exchange reaction
model.S(end+1,:) = 0;
model.S(:,end+1) = 0;
model.S(end,end) = -1;

biomass_prot = find(ismember(model.metNames,'protein_pool_biomass'));
model.rxns(end+1) = {'EX_M10013[c]_CLP'}; %M10013 is metID of protein pool
model.mets(end+1) = model.mets(biomass_prot);
model.lb(end+1) = -1000;
model.ub(end+1) = 0;
model.rev(end+1) = 0;
model.c(end+1) = 0;
model.b(end+1) = model.b(biomass_prot);
model.rxnNames(end+1) = {'protein-uptake pool exchange reaction for CLP model'};
model.grRules(end+1) = model.grRules(biomass);
model.rxnGeneMat(end+1,:) = model.rxnGeneMat(biomass,:);
model.subSystems(end+1) = model.subSystems(biomass);
model.eccodes(end+1) = model.eccodes(biomass);
model.metNames(end+1) = {'protein-uptake pool'};
model.metComps(end+1) = model.metComps(biomass_prot);
model.inchis(end+1) = model.inchis(biomass_prot);
model.metFormulas(end+1) = model.metFormulas(biomass_prot);
model.rxnReferences(end+1) = model.rxnReferences(biomass);
model.rxnConfidenceScores(end+1) = model.rxnConfidenceScores(biomass);
model.metCharges(end+1) = model.metCharges(biomass_prot);

% step 3: add protein-uptake pool breakdown reaction
protein_uptake = find(contains(model.rxns, 'EX_M10013[c]_CLP'));
model.S(:,end+1) = model.S(:,protein_uptake);

%breakdown to amino acids (stoichiometry is identical to 'Protein pool for biomass reaction')
amino_acids       = {'alanine';	... 
                    'arginine';	... 
                    'asparagine'; ... 
                    'aspartate';	... 
                    'cysteine';	... 
                    'glutamate';	... 
                    'glutamine';	... 
                    'glycine';	... 
                    'histidine';	... 
                    'isoleucine';	... 
                    'leucine';	... 
                    'lysine';	... 
                    'methionine';	... 
                    'phenylalanine';	... 
                    'proline';	... 
                    'serine';	... 
                    'threonine';	... 
                    'tryptophan';	... 
                    'tyrosine';	... 
                    'valine'};    
AA_IDx = find(ismember(model.metNames,amino_acids)); 
AA_IDx = AA_IDx(find(ismember(model.metComps(AA_IDx), 4))); %metComps 4 is cytoplasmic
AA_IDx([6 7]) = AA_IDx([7 6]); %glutamine and glutamate positions are swapped when sorted on 3-letter codes
protein_pool = find(ismember(model.rxnNames, 'Protein pool for biomass reaction'));
protein_pool_mets = find(model.S(:,protein_pool));
protein_pool_stoichiometry = model.S(protein_pool_mets,protein_pool);
protein_pool_mets = model.metNames(protein_pool_mets(find(protein_pool_stoichiometry > 0)));
protein_pool_stoichiometry = protein_pool_stoichiometry(protein_pool_stoichiometry > 0);
protein_pool_sort = table(protein_pool_mets, protein_pool_stoichiometry);
protein_pool_sort = sortrows(protein_pool_sort);
AA_stoichiometry = table2array(protein_pool_sort(2:21,2)); 
model.S(AA_IDx,end) = AA_stoichiometry;

model.rxns(end+1) = {'protein_breakdown_CLP'};
model.lb(end+1) = 0;
model.ub(end+1) = 1000;
model.rev(end+1) = 0;
model.c(end+1) = 0;
model.rxnNames(end+1) = {'protein-uptake pool breakdown reaction for CLP model'};
model.grRules(end+1) = model.grRules(protein_uptake);
model.rxnGeneMat(end+1,:) = model.rxnGeneMat(protein_uptake,:);
model.subSystems(end+1) = model.subSystems(protein_uptake);
model.eccodes(end+1) = model.eccodes(protein_uptake);
model.rxnReferences(end+1) = model.rxnReferences(protein_uptake);
model.rxnConfidenceScores(end+1) = model.rxnConfidenceScores(protein_uptake);

% step 4: if open_exchange_mets etc are given, proceed to constrain CLPmodel
% (if not given, output CLPmodel is without CLP-specific constraints)
if ~isempty(open_exchange_mets)
    %close all exchange reactions in model
	all_exchange_rxns = getExchangeRxns(model); %this is rxn IDs
	all_exchange_rxns = find(contains(model.rxns, all_exchange_rxns)); %rxn indices
	model.ub(all_exchange_rxns) = 0;
	model.lb(all_exchange_rxns) = 0;
	
	%find reaction indices of open_exchange_mets
	all_exchange_metIDx = zeros(length(all_exchange_rxns),1);
	for i = 1:length(all_exchange_rxns)
		ex_met = find(model.S(:,all_exchange_rxns(i)));
		all_exchange_metIDx(i,1) = ex_met;
	end
	all_exchange_mets = model.metNames(all_exchange_metIDx);
	
	open_exchange_IDx = zeros(length(open_exchange_mets),1);
	for j = 1:length(open_exchange_mets)
		query_met = find(ismember(all_exchange_mets, open_exchange_mets(j)));
		open_exchange_IDx(j,1) = all_exchange_rxns(query_met);
	end
	
	%constrain open_exchange_IDx by the given lb and ub
	model.ub(open_exchange_IDx) = open_exchange_ub;
	model.lb(open_exchange_IDx) = open_exchange_lb;
	
	%step 5: QC: maxmize biomass and ATP
	model.c(:) = 0;
	model.c(open_exchange_IDx(find(ismember(open_exchange_mets, 'biomass')))) = 1; 
	max_biomass_sol = optimizeCbModel(model);
	if isempty(max_biomass_sol.x)
		EM='The CLPmodel has no feasible solution';
		dispEM(EM);
	end
	
	model.c(:) = 0;
	model.c(open_exchange_IDx(find(ismember(open_exchange_mets, 'ATP')))) = 1; 
	max_ATP_sol = optimizeCbModel(model);
	if isempty(max_ATP_sol.x)
		EM='The CLPmodel has no feasible solution';
		dispEM(EM);
    end
   
    %Step 6: constrain ATP production to max
    %model.ub(open_exchange_IDx(find(ismember(open_exchange_mets, 'ATP')))) = max_ATP_sol.f;
    %model.lb(open_exchange_IDx(find(ismember(open_exchange_mets, 'ATP')))) = max_ATP_sol.f;
    
    %close reactions that cannot carry flux
    %close_zeroflux_rxns = [];
    %counter=0;
    %for k = 1:length(model.rxns)
    %    model.c(:) = 0;
    %    model.c(k) = 1; 
    %    close_rxn_sol = optimizeCbModel(model);
    %    if close_rxn_sol.f == 0
    %        model.c(k) = -1; 
    %        close_rxn_sol_rev = optimizeCbModel(model);
    %        if close_rxn_sol_rev.f == 0
    %           close_zeroflux_rxns = [close_zeroflux_rxns;k]; 
    %           counter = counter+1;
    %        end
    %    end
    %    if rem(k,10)==0
    %       disp(['finished constraining ' num2str(k) '/' num2str(length(model.rxns)) ' reactions']) 
    %    end
    %end
    %model.ub(close_zeroflux_rxns) = 0;
	%model.lb(close_zeroflux_rxns) = 0;
    
    %re-constrain open_exchange_IDx by the given lb and ub
	%model.ub(open_exchange_IDx) = open_exchange_ub;
	%model.lb(open_exchange_IDx) = open_exchange_lb;
    
	%re-set objective function to maximize ATP
	%model.c(:) = 0;
	%model.c(open_exchange_IDx(find(ismember(open_exchange_mets, 'ATP')))) = 1; 
	
end

%output CLPmodel ready for random sampling
model.description = 'CLPmodel';
CLPmodel = model;

end