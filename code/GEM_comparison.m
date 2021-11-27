%
% compare GEMs
% 

% load table for species
T = readtable('../data/GEMs/Table_S1.csv');

% load models into a structure GEMs
%clear GEMs
GEMs.id = T.taxonID;
for i=1:length(GEMs.id)
    if GEMs.id(i) == 9606
        filename = '../data/GEMs/HumanGEM_1_6.mat';
        load(filename);
        model = ihuman;
    else
        filename = ['../data/GEMs/' num2str(GEMs.id(i)) '.mat'];
        load(filename);
    end
    GEMs.model{i,1} = model;
    GEMs.rxnNum{i,1}  = length(model.rxns);
    GEMs.metNum{i,1}  = length(model.mets);
    GEMs.geneNum{i,1} = length(model.genes);
end
GEMs.category = T.Category;

% model comparison
compStruct = compareMultipleModels(GEMs.model,false,true,GEMs.category);
% based on this overview analysis, 4 figures can be made:
% 1. Structural similarity, 2. tSNE plot, 3. Subsystem coverage, and full


%% tSNE
%proj_coords = tsne(double(compStruct.reactions.matrix'),'Algorithm','exact','Distance','hamming','NumDimensions',3,'Perplexity',5);
proj_coords = tsne(double(compStruct.reactions.matrix'),'Distance','hamming','NumDimensions',3);
axis_labels = {'tSNE 1';'tSNE 2';'tSNE 3'};

% prepare groupVector
[groupNames,~,groupVector] = unique(GEMs.category);

% plot structure comparison results
figure; hold on;
color_data = groupVector;
color_palette = lines(length(groupNames));
colormap(color_palette);

scatter3(proj_coords(:,1),proj_coords(:,2),proj_coords(:,3),120,color_data,'filled');

% use text label to present GEO id, gender and month
text(proj_coords(:,1), proj_coords(:,2)+5, proj_coords(:,3)+5, T.CommonName);

view(135,25);  
grid on
xlabel(axis_labels{1}); ylabel(axis_labels{2}); zlabel(axis_labels{3});
set(gca,'FontSize',15,'LineWidth',1.25);
title('Model Structural Comparison','FontSize',16,'FontWeight','bold')

% add legend
for i = 1:length(groupNames)
    h(i) = scatter3([],[],[],120,color_palette(i,:),'filled');
end
[l, hobj, hout, mout] = legend(h,groupNames);
M = findobj(hobj,'type','patch');
set(M,'MarkerSize',11); 


%% check deviation in subsystems
subCoverage = (compStruct.subsystems.matrix - mean(compStruct.subsystems.matrix,2))./mean(compStruct.subsystems.matrix,2)*100;
subs_ind = any(abs(subCoverage) > 25,2);   % filter subsystems to plot
subs_ind(1) = 0;   
cmap = custom_cmap('redblue');
genHeatMap(subCoverage(subs_ind,:), 'colNames', T.CommonName, 'rowNames', compStruct.subsystems.ID(subs_ind),...
        'clusterDim','both','clusterDist','euclidean','colorMap', cmap,...
        'colorBounds', [-100, 100],'gridColor', 'k');
set(gca,'FontSize',15);
xtickangle(45)


%% comapre rxn content using clustergram
cg = clustergram(compStruct.structComp, 'Symmetric', false,...
    'Colormap', flipud(bone), 'RowLabels', T.CommonName,...
    'DisplayRatio', [0.1500 0.1500],...
    'ColumnLabels', T.CommonName);
plot(cg);
set(cg,'LineWidth',1);
xtickangle(45)
set(gca,'FontSize',15);
set(gca,'LineWidth',1);


%% scatter plot for model components
figure; hold on;
box on
sz = 100;
range = 1:31;
scatter(range, cell2mat(GEMs.rxnNum), sz, 'r', '^', 'filled');
scatter(range, cell2mat(GEMs.metNum), sz, 'b', '', 'filled');
scatter(range, cell2mat(GEMs.geneNum), sz, 'm', 'h', 'filled');
xtickangle(45)

set(gca,'FontSize',16,'FontName','Helvetica');
xlim([0 32]);
ylabel('# GEM components','FontSize',18,'FontName','Helvetica');
xlabel('Species','FontSize',18,'FontName','Helvetica');
set(gca,'xtick',range,'xticklabel',T.CommonName)
% prepare legend
[l, hobj, hout, mout] = legend({'Reaction';'Metabolite';'Gene'},'FontSize',16,'FontName','Helvetica','location','ne');
M = findobj(hobj,'type','patch');
set(M,'MarkerSize',16);  


