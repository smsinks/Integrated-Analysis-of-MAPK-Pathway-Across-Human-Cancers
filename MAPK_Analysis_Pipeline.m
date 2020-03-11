%% Comprehensive Molecular Characterisation of MAPK Pathways in Cancer

% With the recent conclusion of the data generation phase of The Cancer
% Genome Atlas (TCGA), there is opportunity for systematic analyses of the
% entire TCGA pan-cancer cohort, including analyses focusing on specific
% oncogenic pathways. The aim of our study was to take a comprehensively
% examine the entire PI3K/AKT/mTOR pathway and its components in the over
% 10,000 human cancers and 32 cancer types profiled by TCGA, these data
% involving multiple molecular profiling platforms, including proteomics.

clc ; clear ; close all ;

% check that the achilles and CCLE datasets have been downloaded and added 
% to the folder
if ~exist('Achilles_sample_info_19Q4.xlsx','file')
    error(['Please download and add to the folder the following files',...
        ' from https://depmap.org/portal/download/ :',...)
        ' CCLE_expression.csv, CCLE_mutations.csv',....
        ', and Achilles_gene_effect_CRISPR_19Q4.csv']) 
end

% check if the new files have been added to the matlab path
if ~exist('setup_env.m','file')
    error(['Please download the MATLAB Connectivity Map toolbox from',...
        ' https://github.com/cmap/cmapM and follow the installation',...
        'instructions to instal'])
end
    
%% Get All the Genes Involved in MAPK Signalling

mapkGenes = readtable('MAPKgenes.xlsx') ;
mapkGenes = unique(mapkGenes) ;

% ============== Get Data for My Genes from cBioPortal ==============

% Either load the cBioPortal data or get from cBioPortal If the data does
% not exist
if exist('mapkDatawithAllMutations.mat','file')
    fprintf('\n Loading cBioPortal MAPK pathway Data \n')
    % this data has matched results load mapkData.mat

    % this dataset contains all the mutations
    load mapkDatawithAllMutations.mat
    
%     error('Get the more than 40000 data from cBioPortal')
else
    fprintf(['\nGetting cBioPortal Data from http:/',...
        '/www.cbioportal.org/public-portal\n'])
    
    % set the bias of the: because mutation have more data. I will get all
    % data that has mutation data here
    mutBias = true;
    
    % get the data from the cBioPortal online repository
    [mutations, cancerStudies, cnaData, mrna, clinicalData ] = ...
        getcBioPortalDataAllStudies(unique(mapkGenes.Gene) , mutBias) ;
    
    % change the CancerStudy to Upper
    mutations.CancerStudy = upper(mutations.CancerStudy);
    cnaData.CancerStudy = upper(cnaData.CancerStudy);
    
    if mutBias ~= true
        % remove the duplicate sample Ids from the mutations data, copy
        % number data and the clinical data
        [~, theUnique] = unique(mutations.SampleIds);
        mutations = mutations(theUnique,:);
        [~, theUnique] = unique(cnaData.SampleIds);
        cnaData = cnaData(theUnique,:);
        
        % arrange the mutations data according to the cnaData
        [~, theLocs] = ismember(cnaData.SampleIds,mutations.SampleIds);
        mutations = mutations(theLocs, :);
        
        % throw in an assession
        assert(all(strcmp(mutations.SampleIds,cnaData.SampleIds)))
        
        % arrange the clinical data according to the mutations data
        theUnique = ismember(clinicalData.SampleIds,mutations.SampleIds);
        clinicalData = clinicalData(theUnique,:);
        [~, theLocs] = ismember(mutations.SampleIds,clinicalData.SampleIds);
        clinicalData = clinicalData(theLocs, :);
        
        % throw in an assession
        assert(all(strcmp(mutations.SampleIds,clinicalData.SampleIds)))
    end
    
    if mutBias ~= true
        % save a copy of the data
        save('mapkData.mat','mutations', 'cancerStudies', 'cnaData',...
            'mrna','clinicalData')
    else
        % save the mutation baised data
        save('mapkDatawithAllMutations.mat','mutations', ...
            'cancerStudies', 'cnaData','mrna','clinicalData')
    end
end

%% Clean Up the Data 

% some of the studies e.g THYROID and BREAST do not have oncotree
% annotations IDs therefore change these
changeCodes = {'LUAD','LUNG';'COADREAD','CRC';'THAP','THYROID';...
    'ESCA','EGC';'BRCA','BREAST' ; } ;
[locThem,these] = ismember(mutations.CancerStudy,changeCodes(:,2));
these(these == 0) = [] ;
mutations.CancerStudy(locThem) = changeCodes(these,1) ;
cnaData.CancerStudy(locThem) = changeCodes(these,1) ;

% add match the clinical data with mutations data arrange the clinical data
% according to the mutations data
[~, theUnique ] = unique(clinicalData.SampleIds);
clinicalData = clinicalData(theUnique,:);
[these] = ismember(clinicalData.SampleIds ,mutations.SampleIds);
clinicalData= clinicalData(these, :);
[these] = ismember(mutations.SampleIds,clinicalData.SampleIds);
mutations = mutations(these, :);

% also change the cancer study from the TMB because it contains many
% different cancer. These are in the clinical data. Therefore, I replace
% the TMB with the IDs from the the clinical data oncotree convert clinical
% cancer study to categorical
clinicalData.CancerStudy = categorical(clinicalData.CancerStudy);
mutations.CancerStudy(ismember(mutations.CancerStudy,'TMB')) = ...
    clinicalData.ONCOTREE_CODE(clinicalData.CancerStudy == 'TMB') ;

% also change these TMB studies within the clinical data
clinicalData.CancerStudy(clinicalData.CancerStudy == 'TMB') = ...
    clinicalData.ONCOTREE_CODE(clinicalData.CancerStudy == 'TMB');

% check the NaN values to empty strings
for pp = 3:width(mutations)
    mutations.(pp) = strrep(mutations.(pp),'NaN','') ;
end

% remove the proteins withh less than 10 samples from the
catsIn = categories( categorical( mutations.CancerStudy)) ;
numOfCats = countcats( categorical( mutations.CancerStudy)) ;

% now get the categoricals with less than 10 samples 
less10Samples = catsIn(numOfCats < 10) ;

checkO = clinicalData(ismember(clinicalData.CancerStudy, ...
    less10Samples) , [1:5,10]) ;
checkO = removevars(checkO, 'AGE');

% get the unique values so that It easy to manually change the oncotree IDs
[~, theUnique] = unique(checkO.CancerStudy) ;
checkStudies = checkO( theUnique, [1,3,end])  ;
checkStudies.NewOncoTree = checkStudies.ONCOTREE_CODE ;

% create an array to change the values in the check studies table of the
% NewOncoTree. change the codes of following cancer studies
getCodes = {'Bladder Cancer','BLCA'; 'Breast Cancer','BRCA' ; ...
    'Colorectal Cancer','COADREAD'; 'Cancer of Unknown Primary',...
    'Unknown'; 'Esophagogastric Cancer','ESCA' ; ...
    'Glioma','GB'; 'Head and Neck Cancer','HNSC' ; ...
    'Melanoma','SKCM';'Non-Small Cell Lung Cancer','NSCLC' ; ...
    'Renal Cell Carcinoma','CCRCC'} ;
[locThem,these] = ismember(checkStudies.CANCER_TYPE, getCodes(:,1) );
checkStudies.NewOncoTree(locThem) = getCodes(these,2) ;

% now change the codes in the mutations, clinical data and copy number data
[locThem,these] = ismember(mutations.CancerStudy, ...
    checkStudies.CancerStudy);
these(these == 0) = [] ;
mutations.CancerStudy(locThem) = checkStudies.NewOncoTree(these) ;
cnaData.CancerStudy(locThem) = checkStudies.NewOncoTree(these);
clinicalData.CancerStudy(locThem) = checkStudies.NewOncoTree(these) ;

% throw in an assession
assert(all(strcmp(mutations.SampleIds,cellstr(clinicalData.SampleIds) )))

clear catsIn numOfCats less10Samples checkO theUnique locThem these ...
    getCodes checkStudies 

%% ====== get alterations of MAPK pathway and put them in one table ======

% create a table with two variable that are required by the function
% getAlterations frequency. One column should be the pathwayName and the
% other should be Genes (which contains a list of gene that are involved in
% that MAPK pathway

mapkPathways = unique(mapkGenes.Pathway);
proteinClass = unique(mapkGenes.proteinClass);

% get all the mapk genes and pathway plus transcriptoin factors
mapkPathwayAlterations = cell( length(mapkPathways),1 );
for ii = 1:length(mapkPathways)
    mapkPathwayAlterations(ii,1) = strcat( mapkPathways(ii),' with Targets');
    mapkPathwayAlterations{ii,2} = strjoin( mapkGenes.Gene( ...
        ismember(mapkGenes.Pathway,mapkPathways(ii) ) )  );
end

% get the mapk genes without transcription factors
mapkPathwayOnly = cell( length(mapkPathways),1 );
for ii = 1:length(mapkPathways)
    mapkPathwayOnly(ii,1) = mapkPathways(ii) ;
    
    % get the mapk genes plus transcription factors
    mapkOnlyGenes = mapkGenes.Gene( ...
        ismember(mapkGenes.Pathway,mapkPathways(ii) ) ) ;
    
    % remove the transcption factors
    locTFs = ismember(mapkOnlyGenes, ...
        mapkGenes.Gene(ismember( mapkGenes.proteinClass,'Target')) ) ;
    mapkOnlyGenes(locTFs) = [];
    
    % add to the table
    mapkPathwayOnly{ii,2} = strjoin(mapkOnlyGenes);
end

% all both alteration in one cell array
mapkPathwayAlterations = [mapkPathwayOnly;mapkPathwayAlterations] ;

% also add all the mapk pathway genes without the Targers to the alteration
% and add the variable names
mapkPathwayAlterations(end+1,:) = {'All MAPK pathways', ...
    strjoin(mapkGenes.Gene(~ismember( mapkGenes.proteinClass,'Target') ))};

% also add all the mapk pathway genes to the alteration  and add the
% variable names
mapkPathwayAlterations(end+1,:) = {'All MAPK with Targets', ...
    strjoin(mapkGenes.Gene) } ;
mapkPathwayAlterations = array2table(mapkPathwayAlterations,...
    'VariableNames',{'pathwayName','Genes'} );

% get all the indivial class of proteins involved in the mapk
mapkProteinsAlterations = cell( length(proteinClass),1 );
for ii = 1:length(proteinClass)
    mapkProteinsAlterations(ii,1) = proteinClass(ii) ;
    mapkProteinsAlterations{ii,2} = strjoin( mapkGenes.Gene( ...
        ismember(mapkGenes.proteinClass,proteinClass(ii) ) ) ) ;
end

mapkProteinsAlterations = array2table(mapkProteinsAlterations,...
    'VariableNames',{'pathwayName','Genes'} );

% === get alterations of each MAPK pathways and put them in one table ===

% This is a table with columns as cancer types and row as MAPK pathways

% get the genes involved in a MAPK pathwy and return only these for the
% copy number data and mutations data
fprintf('\n Getting MAPK Pathway Alterations \n')
MAPK_Alterations_Results = find_MAPK_AlterationFreq( ...
    mapkPathwayAlterations,...
    mutations) ;

% also get percentage of alteration involved in MAPK
MAPK_Proteins_Results = find_MAPK_AlterationFreq( ...
    mapkProteinsAlterations,...
    mutations) ;

% get the mean alteration of the each MAPK protein class across various
% cancer
MAPK_Mean_Protein_Alterations = array2table( ...
    mean( MAPK_Proteins_Results{:,2:end} ) ,'VariableNames',...
    MAPK_Proteins_Results.Properties.VariableNames(2:end) ) ;

MAPK_Mean_Protein_Alterations = rows2vars(MAPK_Mean_Protein_Alterations);
MAPK_Mean_Protein_Alterations.Properties.VariableNames = ...
    {'Pathway','PercentAltered'} ;

% also get the mean alterations of the Four MAPK pathways Across cancers
MAPK_Mean_Pathway_Alterations = array2table( ...
    mean( MAPK_Alterations_Results{:,2:end} ) ,'VariableNames',...
    MAPK_Alterations_Results.Properties.VariableNames(2:end) ) ;
MAPK_Mean_Pathway_Alterations = ...
    rows2vars(MAPK_Mean_Pathway_Alterations) ;
MAPK_Mean_Pathway_Alterations.Properties.VariableNames = ...
    {'Pathway','PercentAltered'} ;

fprintf('\n Getting Frequency of Alterations \n')
% get the frequency of gene alterations of all MAPK genes so that I can
% finally draw the pathway

% ====== Plot 5: Each MAPK Pathways Alterations Across Cancers ===========
eachGeneMutations  = sum( ~ismissing( mutations{:,3:end} ) ) / ...
    height(mutations) ;
eachGeneMutations = array2table(eachGeneMutations *100,'VariableNames',...
    mutations.Properties.VariableNames(3:end) );
eachGeneMutations = rows2vars(eachGeneMutations );
eachGeneMutations.Properties.VariableNames = ...
    {'Gene','PercentAltered'} ;
eachGeneMutations = sortrows(eachGeneMutations,'PercentAltered','descend');

%  ========  Finally Get Mutations of Each Gene Across Cancers =========
% get the mutations from the table by converting them to double
AcrossCancerMutations = [ mutations(:,1) ,...
    array2table( double(~ismissing( mutations{:,3:end} )) ) ];
AcrossCancerMutations.Properties.VariableNames(2:end) = ...
    mutations.Properties.VariableNames(3:end) ;

% convert to categorical so group stats can work properly
AcrossCancerMutations.CancerStudy = categorical(...
    AcrossCancerMutations.CancerStudy);

% get the groups starts
AcrossCancerMutations = grpstats(AcrossCancerMutations,...
    'CancerStudy','sum');

% divide each row by the total number of sample in each cancer study to get
% the total percentage of mutations in for each gene
for ii = 3:width(AcrossCancerMutations)
    % conver to the percentage of tumours with mutations
    AcrossCancerMutations.(ii) = ...
        (AcrossCancerMutations.(ii)./AcrossCancerMutations.GroupCount)*100;
end

fprintf('\n Getting Overall Alterations for Each Gene \n')
% finally change the variables names by remove sum from each genes
AcrossCancerMutations.Properties.VariableNames(3:end) = ...
    extractAfter(AcrossCancerMutations.Properties.VariableNames(3:end),...
    'sum_');

% set up the alteration of each gene across across 
mapkGeneMutations = innerjoin(mapkGenes, eachGeneMutations);
[~, these] = unique( mapkGeneMutations.Gene);
mapkGeneMutations = mapkGeneMutations(these, :) ;
mapkGeneMutations.Properties.VariableNames(end) = "meanAltered" ;

% find the maximim alteartion of each gene across the various tumours
eachGeneMax = rows2vars(AcrossCancerMutations(:, [1,3:end]), ...
    'VariableNamesSource','CancerStudy' );
eachGeneMax.maxAltered = max(eachGeneMax{:,2:end} ,[],2) ;
eachGeneMax.Properties.VariableNames(1) = "Gene" ;
eachGeneMax = eachGeneMax(:, {'Gene','maxAltered'}) ;

% now join this table to the one which shows mean gene alteration in each
% gene across cancers 
mapkGeneMutations = innerjoin(mapkGeneMutations, eachGeneMax)  ;

% print out the total number of samples
TotalSample =  sum(AcrossCancerMutations.GroupCount);

clear mapkPathways ii locTFs mapkPathwayOnly ...
    mapkOnlyGenes aa locThem these changeCodes theLocs theUnique ...
    mapkProteinsAlterations mapkPathwayAlterations pp

%% Visualise the results 

% get only the alteration with the targets
mapkAlter4plot = MAPK_Mean_Pathway_Alterations( ...
    contains( MAPK_Mean_Pathway_Alterations.Pathway,'With'), :) ;
mapkAlter4plot.Pathway = regexprep( mapkAlter4plot.Pathway, ...
    'With+\w*',' Pathway') ;

% plot the frequence of MAPK pathway alterations for each path
figure()
subplot(1,2,1)
mapkAlter4plot.Pathway = categorical(mapkAlter4plot.Pathway);
barh(mapkAlter4plot.Pathway, mapkAlter4plot.PercentAltered, ...
    'DisplayName','mapkAlter4plot.PercentAltered')

% adjust the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold')
xlabel('Tumours with Alterations (%) ')
title('MAPK Pathways Mutations','FontSize',16)

% plot th extent of pathway protein alteration across the MAPK pathways 
subplot(1,2,2)
MAPK_Mean_Protein_Alterations.Pathway = ...
    categorical(MAPK_Mean_Protein_Alterations.Pathway);
barh(MAPK_Mean_Protein_Alterations.Pathway, ...
    MAPK_Mean_Protein_Alterations.PercentAltered, ...
    'DisplayName','MAPK_Mean_Protein_Alterations.PercentAltered')

% adjust the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold')
xlabel('Tumours with Alterations (%) ')
title('MAPK Mutations by Encoded Protein','FontSize',16)

% produce a clustergram of the MAPK alterations across cancers. This should
% exclude the of the table columns 
cgo1 = clustergram( MAPK_Alterations_Results{:, 6:end-2}',...
    'rowlabels', ...
    regexprep( ...
    MAPK_Alterations_Results.Properties.VariableNames( 6:end-2),...
    'With+\w*',' Pathway') ,...
    'columnlabels',cellstr(MAPK_Alterations_Results.CancerStudy),...
    'colormap', jet ,'standardize','column',...
    'ColumnPDist','euclidean','Dendrogram',6 ,'Linkage','complete');

% set the font size of the clustergram 
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 35,'FontWeight','bold')


% add the mutation frequency to the top of the clustergram
% move up the dendograms
cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
    'Color',{[1 1 0],[0.6 0.6 1]});
set(cgo1,'ColumnGroupMarker',cm)

% get the purity data and plot it on top of the clustergram
cgAxes = plot(cgo1);
% set(cgAxes, 'Clim', [-2,2])

% add a heatmap on top to show sample purity
purityBar =  MAPK_Alterations_Results ;
[~, PosTargets] = ismember(cgo1.ColumnLabels , ...
    MAPK_Alterations_Results.CancerStudy);
purityBar = purityBar(PosTargets,:) ;

% throw in an assertion 
assert(all(strcmp(cgo1.ColumnLabels',...
    cellstr(purityBar.CancerStudy))) ,'Size Mismatch') ;
axes('position',[0.2423 0.735 0.5760 0.15] );
bar(purityBar.AllMAPKWithTargets','FaceColor',[0.717 0.274 1]);
ylabel('% mutated')

% adjust the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'XTickLabel',[],'XLim',[0.5 height(purityBar)+0.5] ,...
    'FontWeight','bold')

% ***********************************************************************
% Alteration within the classical pathway are mutually exclusive to those
% in the p38 and JNK pathway 
% ************************************************************************

% % cluster the gene alteration across cancer types 
% clustergram( AcrossCancerMutations{:, 3:end}',...
%     'rowlabels',AcrossCancerMutations.Properties.VariableNames(3:end),...
%     'columnlabels',cellstr(AcrossCancerMutations.CancerStudy),...
%     'colormap', redbluecmap ,'standardize','column',...
%     'ColumnPDist','euclidean','Dendrogram',6 ,'Linkage','complete');

% ================ Table 1: Tumours Being Studies ========================

% get the cancerStudy and their Codes.
cancerClass = cancerStudies(:,[1:2,end]);

% extact the cancer codes and the cancer study names
cancerClass.cancerTypeId = categorical( upper( ...
    extractBefore(cellstr(cancerClass.cancerTypeId),'_') ) );
cancerClass.name = strtrim( extractBefore(cancerClass.name,'(' ) )  ;

% add the cancer study details to the table
cancerClass = [ cancerClass, cancerStudies(:,3)] ;

% get the cell lines from the data
cellLines = clinicalData(~cellfun(@isempty, clinicalData.SAMPLE_CLASS),:);

% ================ Plot 2: MAPK Pathways in DrawIO ======================
% in DrawIO

% Also try a network in yEd
MAPKnetwork = createNetwork2(mapkGenes.Gene);

% trim the network to return only proteins that are in the initial gene
% list
genesInList = all( [ ismember(MAPKnetwork.Protein1,mapkGenes.Gene) ,...
    ismember(MAPKnetwork.Protein2,mapkGenes.Gene) ] , 2) ;
trimmedMAPKnetwork = MAPKnetwork(genesInList,:);

% remove all transcription factors from the Proteins1 so that they do not
% enfluence the network connectivity
trimmedMAPKnetwork = trimmedMAPKnetwork( ...
    ismember(trimmedMAPKnetwork.Protein1, ...
    mapkGenes.Gene( ~ismember( mapkGenes.proteinClass,'Target'))) , :) ;

% add the protein class to the trimmed network 
[~, locProtClass] = ismember( trimmedMAPKnetwork.Protein1 , ...
    mapkGenes.Gene ) ;
trimmedMAPKnetwork.proteinClass = mapkGenes.proteinClass(locProtClass) ;
trimmedMAPKnetwork.pathway = mapkGenes.Pathway(locProtClass) ;

% create a graph to edit in cytoscape
writetable(trimmedMAPKnetwork,'yED_MAPK_Network.xlsx','Sheet',1) ;
writetable(MAPKnetwork,'yED_MAPK_NetworkComplete.xlsx','Sheet',1) ;

% It turn out that the network is too complex for the automatic yEd
% Algorythm to lay it out clearly. Therefore the manual method using drawIO
% will be employed to create the beautiful figure 1!!!!!

% save some supplementary files 
writetable( mapkGenes ,'Supplementary File 1.xlsx', ...
    'Sheet','MAPK Pathway Genes') ;
writetable( cancerClass ,'Supplementary File 1.xlsx', ...
    'Sheet','Cancer Studies') ;
writetable( mapkGeneMutations,'Supplementary File 1.xlsx',...
    'Sheet', 'Specific MAPK Gene Mutations');
writetable( AcrossCancerMutations , 'Supplementary File 1.xlsx',...
    'Sheet', 'Cancer-MAPK Gene Mutations');
writetable( MAPK_Proteins_Results , 'Supplementary File 1.xlsx',...
    'Sheet', 'Cancer-Protein Type Mutations');
writetable( mapkAlter4plot , 'Supplementary File 1.xlsx',...
    'Sheet', 'Cancer-MAPK pathways Mutations');
writetable( MAPK_Alterations_Results , 'Supplementary File 1.xlsx',...
    'Sheet', 'MAPK Alteration Across Cancers');

clear purityBar mapkAlter4plot cgo1 locProtClass

%% Correlation of Pathway Mutations to Patient Survival

% Get the correlation between each of the three mapk kinase pathway
% alteration and the treatment outcomes in term of the overall survival and
% disease free pregression of patients of each cancer with or without
% alteration in each of the four MAPK pathway

% throw in an assession
assert(all(strcmp(mutations.SampleIds,cellstr(clinicalData.SampleIds) )))

% put together the copy number data and the mutations data to find tumours
% with alteration in the MAPK pathways at both mutations and copy number
% changes
mutatedMAPKgenes = mutations ;
mutatedMAPKgenes{:,3:end} = num2cell( ...
    double( ~cellfun(@isempty,mutatedMAPKgenes{:,3:end}) )  ) ;

% now find the samples that have alterations in mapk genes get the low
% metabolism group: Get for new analysis!!!
mapkAlteredIds = double( sum(cell2mat(mutatedMAPKgenes{:,3:end}),2 ) > 0 );

% convert the logical into a cell and add it as a variable to the clinical
% data table
mapkStatus = ...
    strrep( (cellstr(num2str(mapkAlteredIds))),'0','Not mutated');
mapkStatus = strrep(mapkStatus,'1','Mutated');
try % you dont want to get an error after the variable has been added to the table
    clinicalData = addvars(clinicalData,mapkStatus,'Before',...
        'ABNORMAL_LYMPHOCYTE_PERCENT') ;
catch
end

% Get the Overall Survival Data and the Disease Free Survival Data and
% delete the missing low from the data
OsData = [clinicalData.OS_MONTHS, ...
    clinicalData.OS_STATUS ,num2cell(mapkAlteredIds)] ;
OsData(any(cellfun(@isempty,OsData),2),:) = [] ;

fprintf('\nNumber of Tumour with no MAPK alterations = %d \n',...
    sum(cell2mat(OsData(:,3)) ))
fprintf('\nNumber of Tumours with MAPK alterations = %d \n',...
    sum(~cell2mat(OsData(:,3)) ))

% Anonotate the groups to either high metabolism or low metabolism and
% perferom K-M Survival Analysis
groups = cellstr( num2str(cell2mat(OsData(:,3)) ) ) ;
groups = strrep(groups,'0','Not mutated') ;
groups = strrep(groups,'1','mutated') ;

% ============= PERFORM OVERALL SURVIVAL ANALYSIS =================
[~, ~, stats] = MatSurv( str2double(OsData(:,1)), ...
    lower(OsData(:,2)) , groups, 'NoRiskTable',true) ;

% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
    'TickDir','out')
ylevel = 0.65 ;
annotationData = str2double(OsData(:,1));
for jj = 1:numel(unique(groups))
    text(max(annotationData)/3,ylevel, ...
        sprintf('Median OS: (%s) = %g\n', ...
        stats.GroupNames{jj},stats.MedianSurvivalTime(jj)),...
        'FontWeight','bold','FontSize',10)
    ylevel = ylevel+0.04 ;
end
title('\bf Overall Survival','FontSize',16)

% ============= PERFORM DISEASE FREE SURVIVAL ANALYSIS =================
% get only the rows of the data that have disease free data
dfsClinicalData = clinicalData(~cellfun( ...
    @isempty,clinicalData.DFS_STATUS),:);
dfsmapkAlteredIds = mapkAlteredIds(~cellfun( ...
    @isempty,clinicalData.DFS_STATUS))  ;

% create for the data for disease free survival: first make the disease
% free data compatible with the matSurv function
dfsData = [dfsClinicalData.DFS_MONTHS, ...
    strrep(regexprep(dfsClinicalData.DFS_STATUS,'/(\w+)',''),...
    'Recurred', 'Relapsed'),num2cell(dfsmapkAlteredIds)] ;
dfsData(any(cellfun(@isempty,dfsData),2),:) = [] ;

fprintf('\nNumber of Tumour with no MAPK alterations = %d \n',...
    sum(cell2mat(dfsData(:,3)) ))
fprintf('\nNumber of Tumours with MAPK alterations = %d \n',...
    sum(~cell2mat(dfsData(:,3)) ))

% Anonotate the groups to either high metabolism or low metabolism and
% perferom K-M Survival Analysis
groups = cellstr( num2str(cell2mat(dfsData(:,3)) ) ) ;
groups = strrep(groups,'0','Not mutated') ;
groups = strrep(groups,'1','mutated') ;

[~, ~, statsDFS] = MatSurv( str2double(dfsData(:,1)), ...
    dfsData(:,2), groups, 'NoRiskTable',true) ;

% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
    'TickDir','out')
ylevel = 0.65 ;
annotationData = str2double(dfsData(:,1));
for jj = 1:numel(unique(groups))
    text(max(annotationData)/3,ylevel, ...
        sprintf('Median DFS: (%s) = %g\n', ...
        statsDFS.GroupNames{jj},statsDFS.MedianSurvivalTime(jj)),...
        'FontWeight','bold','FontSize',10)
    ylevel = ylevel+0.04 ;
end
title('\bf Disease Free Survival','FontSize',16)

% clear some variables
clear ylevel annotationData groups OsData dfsData metabolicStatus ...
    lowMetaLoc lowMetabolism statsDFS stats ans jj heatData ...
    participants tempMutations dfsmapkAlteredIds dfsClinicalData ...
    mapkStatus mapkAlteredIds mutatedMAPKgenes mapkAlter4plot ...
    MAPKnetwork

%% Find the Impact of Each pathway or Protein Class on Clinical Outcomes

% here we using the proteinClass and pathway from the data
outcomeVar = mapkGenes.Properties.VariableNames(3:end);
clinicalOutcomes = ["overallSurvival","diseaseFreeSurvival"];

% loop over the disease outcomes
for pp = 1:length(clinicalOutcomes)
    % loop over the the pathways and proteinsclasses
    for kk = 1:length(outcomeVar)
        fprintf('\n Running suvival analysis for %s \n', outcomeVar{kk})
        % covert the pathways to categorical for faster indexing
        mapkGenes.(outcomeVar{kk}) =categorical(mapkGenes.(outcomeVar{kk}));
        
        % get the genes altered for each pathway
        proteinORpathways = unique(mapkGenes.(outcomeVar{kk})) ;
        
        % preallocated the appendTable for the samples with mutations in
        % each of the MAPK signalling pathways
        alterationSummary = zeros(height(mutations), ...
            length(proteinORpathways));
        
        for ii = 1:length(proteinORpathways)
            % get the mutations in each pathway using the alterations of
            % the genes for each of those pathways
            curPathwayTable = mutations{ :, ismember( ...
                mutations.Properties.VariableNames, ...
                mapkGenes.Gene( ...
                mapkGenes.(outcomeVar{kk}) == proteinORpathways(ii) ) ) };
            
            % now find the samples that have alterations in mapk genes get
            % the low metabolism group: Get for new analysis!!!
            alterationSummary(:,ii) = ...
                any(~cellfun(@isempty,curPathwayTable),2);
            
        end
        
        % now names the pathways that are altered for each samples
        osMutations = cell(height(mutations), 1);
        for ii = 1:length(proteinORpathways)
            % find the locatinos of the current pathays have mutations
            locMuts = alterationSummary(:,ii) == true ;
            
            % put the pathways in the cell array
            osMutations(locMuts) = cellstr( proteinORpathways(ii) );
        end
        
        % find the missing cell array section and name them as not mutated
        osMutations(cellfun(@isempty,osMutations) ) =...
            {'not mutated'};
        
        % annotated the samples with multiple mapk pathway alterations
        osMutations( sum(alterationSummary, 2) > 1) = {'multiple'} ;
        
        % perform the overall survival for these groups of pathways get the
        % clinical data for the HM cancers
        try
            if strcmp(outcomeVar{kk},'Pathway')
                clinicalData = addvars(clinicalData,osMutations,...
                    'NewVariableNames',{'alteredMAPKpathway'},'After',...
                    'CANCER_TYPE');
                % convert to categorical
                clinicalData.alteredMAPKpathway = ...
                    categorical(clinicalData.alteredMAPKpathway) ;
            elseif strcmp(outcomeVar{kk}, 'proteinClass')
                clinicalData = addvars(clinicalData,osMutations,...
                    'NewVariableNames',{'alteredMAPKproteins'},'After',...
                    'CANCER_TYPE');
                % convert to categorical
                clinicalData.alteredMAPKproteins = ...
                    categorical(clinicalData.alteredMAPKproteins) ;
            end
        catch
            fprintf('\nalteredMAPKpathway in already in the table\n\n')
        end
        
        % Get the Overall Survival Data and the Disease Free Survival Data
        % and delete the missing low from the data
        if strcmp(outcomeVar{kk},'Pathway')
            % convert to categorical
            clinicalData.alteredMAPKpathway = ...
                categorical(clinicalData.alteredMAPKpathway) ;
            
            % Get the Overall Survival Data and the Disease Free Survival
            % Data and delete the missing low from the data
            if strcmp(clinicalOutcomes(pp),'overallSurvival')
                OsData = [clinicalData.OS_MONTHS, ...
                    clinicalData.OS_STATUS ,  ...
                    cellstr(clinicalData.alteredMAPKpathway)  ] ;
                OsData(any(cellfun(@isempty,OsData),2),:) = [] ;
            elseif strcmp(clinicalOutcomes(pp),'diseaseFreeSurvival')
                % =======PERFORM DISEASE FREE SURVIVAL ANALYSIS ==========
                % get only the rows of the data that have disease free data
                dfsClinicalData = clinicalData(~cellfun( ...
                    @isempty,clinicalData.DFS_STATUS),:);
                
                % make the outcome of only two types
                dfsClinicalData.DFS_STATUS = strrep(regexprep( ...
                    dfsClinicalData.DFS_STATUS,'/(\w+)',''),...
                    'Recurred', 'Relapsed');
                % finally get the dfs data
                OsData = [dfsClinicalData.DFS_MONTHS, ...
                    dfsClinicalData.DFS_STATUS,  ...
                    cellstr(dfsClinicalData.alteredMAPKpathway)  ] ;
                OsData(any(cellfun(@isempty,OsData),2),:) = [] ;
            end
            
            % get the number of tumours that belong to each cluster
            summary(clinicalData.alteredMAPKpathway)
            
            % specific the protein names to be used on the plots
            proteinNames = categories(clinicalData.alteredMAPKpathway);
            
        elseif strcmp(outcomeVar{kk}, 'proteinClass')
            if strcmp(clinicalOutcomes(pp),'overallSurvival')
                OsData = [clinicalData.OS_MONTHS, ...
                    clinicalData.OS_STATUS , ...
                    cellstr(clinicalData.alteredMAPKproteins)  ] ;
                OsData(any(cellfun(@isempty,OsData),2),:) = [] ;
            elseif strcmp(clinicalOutcomes(pp),'diseaseFreeSurvival')
                % ======= PERFORM DISEASE FREE SURVIVAL ANALYSIS ========
                % get only the rows of the data that have disease free data
                dfsClinicalData = clinicalData(~cellfun( ...
                    @isempty,clinicalData.DFS_STATUS),:);
                
                % make the outcome of only two types
                dfsClinicalData.DFS_STATUS = strrep(regexprep( ...
                    dfsClinicalData.DFS_STATUS,'/(\w+)',''),...
                    'Recurred', 'Relapsed');
                
                % finally get the dfs data
                OsData = [dfsClinicalData.DFS_MONTHS, ...
                    dfsClinicalData.DFS_STATUS,  ...
                    cellstr(dfsClinicalData.alteredMAPKproteins)  ] ;
                OsData(any(cellfun(@isempty,OsData),2),:) = [] ;
            end
            
            % get the number of tumours that belong to each cluster
            summary(clinicalData.alteredMAPKproteins)
            
            % specific the protein names to be used on the plots
            proteinNames = categories(clinicalData.alteredMAPKproteins);
        end
        
        % Anonotate the groups to either high metabolism or low metabolism
        % and perferom K-M Survival Analysis
        groups = OsData(:,3) ;
        if strcmp(outcomeVar{kk},'Pathway') && ...
                strcmp(clinicalOutcomes(pp),'overallSurvival')
            NoriskTableVar = true ; % make this false to pruduce the risk table
        elseif strcmp(outcomeVar{kk}, 'proteinClass')
            NoriskTableVar = true ;
        end
        
        % ============= PERFORM OVERALL SURVIVAL ANALYSIS =================
        [~, ~, stats] = MatSurv( str2double(OsData(:,1)), ...
            lower(OsData(:,2)) ,OsData(:,3), ...
            'NoRiskTable',NoriskTableVar ,'PairWiseP',true ,...
            'setLegendWithNumbers' , true) ;
        
        % change the groups names in the pairwise stats table to make them
        % easily readble
        if strcmp(clinicalOutcomes(pp),'overallSurvival')
            if strcmp(outcomeVar{kk}, 'proteinClass')
                osProteinsPairwise = struct2table( stats.ParwiseStats ) ;
                osProteinsPairwise.GroupNames = stats.ParwiseName ;
                osProteinsPairwise = ...
                    sortrows(osProteinsPairwise,'p_MC','ascend');
            elseif strcmp(outcomeVar{kk},'Pathway')
                osPathwayPairwise = struct2table( stats.ParwiseStats );
                osPathwayPairwise.GroupNames = stats.ParwiseName ;
                osPathwayPairwise = ...
                    sortrows(osPathwayPairwise,'p_MC','ascend');
            end
        elseif strcmp(clinicalOutcomes(pp),'diseaseFreeSurvival')
            if strcmp(outcomeVar{kk}, 'proteinClass')
                dfsProteinsPairwise = struct2table( stats.ParwiseStats );
                dfsProteinsPairwise.GroupNames = stats.ParwiseName ;
                dfsProteinsPairwise = ...
                    sortrows(dfsProteinsPairwise,'p_MC','ascend');
            elseif strcmp(outcomeVar{kk},'Pathway')
                dfsPathwayPairwise = struct2table( stats.ParwiseStats );
                dfsPathwayPairwise.GroupNames = stats.ParwiseName ;
                dfsPathwayPairwise = ...
                    sortrows(dfsPathwayPairwise,'p_MC','ascend');
            end
        end
        
        % annotate the plot
        set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
            'TickDir','out')
       
        % get the data to use for adding annotations to the plots          
        % add the title to the plot
        if pp == 1
            if strcmp(outcomeVar{kk},'Pathway')
                title(['\bf OS:',...
                    ' Mutations in Specific MAPK Pathways Genes'], ...
                'FontSize',15)
            elseif strcmp(outcomeVar{kk},'proteinClass')
                title(['\bf OS:',...
                    ' Types of MAPK Encoding Genes'],'FontSize',15)
            end
        elseif pp == 2        
            % add the title to the plot
            if strcmp(outcomeVar{kk},'Pathway')
                title(['\bf DFS:',...
                    ' Mutations in Specific MAPK Pathways Genes'], ...
                    'FontSize',15)
            elseif strcmp(outcomeVar{kk},'proteinClass')
                title(['\bf DFS:',...
                    ' Types of MAPK Encoding Genes'],...
                    'FontSize',15)
            end         
        end
        
        % create pie chart showing the alteration by get the labels to
        % use in the pie chart
        if pp == 1 && strcmp(outcomeVar{kk},'Pathway') 
            pieLabels = clinicalData(:, {'alteredMAPKpathway','AGE'} );
            pieLabels.AGE = str2double(pieLabels.AGE) ;
            pieLabels = grpstats(pieLabels,'alteredMAPKpathway') ;
            
            % create the string variable for the pie chart labels
            finalLabels = strcat(cellstr(pieLabels.alteredMAPKpathway),...
                strcat( strcat(':XX',num2str(pieLabels.GroupCount) ) ,...
                'XXTumours') ) ;
            
            finalLabels = strrep( strrep(finalLabels,' ', '') ,'XX',' ') ;
            
            % specificy the parts of the pie chart to explode
            explode = pieLabels.GroupCount / sum(pieLabels.GroupCount) ...
                < 0.05 ;
            
            % finally produce the pie chart
            figure()
            pie(pieLabels.GroupCount ,explode, finalLabels)
            
        elseif pp == 1 && strcmp(outcomeVar{kk},'proteinClass')
            pieLabels = clinicalData(:, {'alteredMAPKproteins','AGE'} );
            pieLabels.AGE = str2double(pieLabels.AGE) ;
            pieLabels = grpstats(pieLabels,'alteredMAPKproteins') ;
            
            % create the string variable for the pie chart labels
            finalLabels = strcat(cellstr(pieLabels.alteredMAPKproteins),...
                strcat( strcat(':XX',num2str(pieLabels.GroupCount) ) ,...
                'XXTumours') ) ;
            
            finalLabels = strrep( strrep(finalLabels,' ', '') ,'XX',' ') ;
            
            % specificy the parts of the pie chart to explode
            explode = pieLabels.GroupCount / sum(pieLabels.GroupCount) ...
                < 0.05 ;
            
            % finally produce the pie chart
            figure()
            pie(pieLabels.GroupCount ,explode, finalLabels)
        end      
    end
end

% save some results to excel 
writetable( osPathwayPairwise,'Supplementary File 2.xlsx', ...
    'Sheet','Pathways - OS Pairwise Comp')
writetable( osProteinsPairwise,'Supplementary File 2.xlsx', ...
    'Sheet','Proteins - OS Pairwise Comp')
writetable( dfsPathwayPairwise,'Supplementary File 2.xlsx', ...
    'Sheet','Pathways - DFS Pairwise Comp')
writetable( dfsProteinsPairwise,'Supplementary File 2.xlsx', ...
    'Sheet','Proteins - DFS Pairwise Comp')

clear locMuts jj ii ylevel pathwayNames groups OsData osMAPKmutations ...
    proteinNames stats osProteinMutations alterationSummary ...
    proteinClasses NoriskTableVar yMove FontS curPathwayTable kk ...
    genesInList mapkAlteredIds osMutations dfsClinicalData pp ...
    clinicalOutcomes outcomeVar explode finalLabels pieLabels

%% save the clinical data to used for plotting the data in tableau

fprintf('\n Saving the clinical Data to Excel \n')
% writetable(clinicalData, 'mapkClinicalData.xlsx')

%% Use Achilles Data To Assess Upstream proteins have more Impact on Cell

% ========== viablity after crisper gene pertubations ==============
% load the achilles data for crisper pertubation and the sample info
fprintf('\n Loading CRISPR data and setting up analysis \n')
crispr = readtable('Achilles_gene_effect_CRISPR_19Q4.csv');
sampleInfo = readtable('Achilles_sample_info_19Q4.csv');

% change the cell line ID to those of the common names from the DepMap Ids
crispr.Properties.VariableNames(1) = sampleInfo.Properties.VariableNames(1);
crispr = innerjoin( sampleInfo(:,1:2), crispr) ;

% clean up the variable names and table
crispr(:,{'DepMap_ID'}) = [] ;
[uniqueNames, uniqueLocations] = ...
    unique( extractBefore(crispr.Properties.VariableNames,'_') ) ;
crispr = crispr(:, uniqueLocations);
crispr.Properties.VariableNames = uniqueNames ;
crispr = movevars(crispr,{'stripped'}, 'Before',1);
crispr.Properties.VariableNames(1) = "cell_line" ;
for ii = 2:width(crispr)
    if iscell(crispr.(ii))
        crispr.(ii) = str2double(crispr.(ii)) ;
    end
end

% remove non cancer cell lines from the crispr data 
fprintf('\n The number of all cell lines = %d\n', height(crispr) ) 
crispr(ismember( crispr.cell_line, ...
    sampleInfo.stripped_cell_line_name( ismember( ...
    sampleInfo.disease ,{'Fibroblast','Engineered'}) ) ) , : )  = [];
fprintf('\n The number of all CANCER cell lines = %d\n', height(crispr) ) 

% Dependent Cell Line: A cell line is considered dependent if it has a
% probability of dependency greater than 0.5.

% Dependency Score: Outcome from DEMETER2 or CERES. A lower score means
% that a gene is more likely to be dependent in a given cell line. A score
% of 0 is equivalent to a gene that is not essential whereas a score of -1
% corresponds to the median of all common essential genes.

% Identify pan-dependent genes as those for whom 90 of cell lines rank the
% gene above a given dependency cutoff. The cutoff is determined from the
% central minimum in a histogram of gene ranks in their 90th percentile
% least dependent line

% save a copy of the crispr data to be use later
crisprAll = crispr ;

% ********************************************************************
% MAPK genes more essential than other genes? ARE the oncogenesmore
% essential than tumour suppressor genes
% ********************************************************************

% check on average if the cancer cell lines are more depend on mapk
% signalling compared to other pathways

nonMapkCrispr = crispr(:, ...
    ~ismember(crispr.Properties.VariableNames ,mapkGenes.Gene)) ;
crispr = [crispr(:,1) crispr(:, ...
    ismember(crispr.Properties.VariableNames ,mapkGenes.Gene)) ] ;

% get the list of oncogene and tumour suppressor genes
% The source databases include:
% The COSMIC consensus cancer gene database
% Uniport Knowledgebase
% TSGene database
% ONCOGene database
% load the files
cosmic = readtable('cancer_gene_censensus.csv');
cosmicOGs = cosmic.GeneSymbol(contains(cosmic.RoleInCancer,'oncogene'));
cosmicTSG = cosmic.GeneSymbol(contains(cosmic.RoleInCancer,'TSG'));
uniprot_OGs = readtable('uniprot_protooncogene.xlsx');
uniprot_TSGs = readtable('uniprot_TSGs.xlsx');
TSGenes = readtable('Human_TSGs.txt');
ONCOgenes = readtable('ONCOgenes_human.txt');

% we put together the oncogenes and TSGs
oncogenes = unique( [cosmicOGs; ...
    extractBefore(uniprot_OGs.EntryName,'_') ; ONCOgenes.OncogeneName ] );
TSGs = unique( [cosmicTSG ; extractBefore(uniprot_TSGs.EntryName,'_'); ...
    TSGenes.Var2(2:end) ] );

fprintf('\n Determining the dependence scores for genes \n')

%% get the mean of the other groups and the mapk pathway for both oncogenes
% and tumour suppressor genes and put them into a new table
meanDependencies = addvars( crispr(:,1), ...
    nanmean( crispr {:, ismember( ...
    crispr.Properties.VariableNames,oncogenes) } ,2) , ...
    nanmean( crispr {:,ismember( ...
    crispr.Properties.VariableNames,TSGs) } ,2) ,...
    nanmean( nonMapkCrispr {:, ismember( ...
    nonMapkCrispr.Properties.VariableNames, oncogenes) } ,2) , ...
    nanmean( nonMapkCrispr{:, ismember( ...
    nonMapkCrispr.Properties.VariableNames,TSGs) } ,2) ,...
    'NewVariableNames',{'mapkOG','mapkTSG','otherOG','otherTSG'}) ;
% 
% % old codes below 
% anovaDepedencies = [meanDependencies.mapkOG ; meanDependencies.mapkTSG; ...
%     meanDependencies.otherOG ; meanDependencies.otherTSG ] ;
% anovaGroups = repelem({'mapkOG';'mapkTSG';'otherOG';'otherTSG'}, ...
%     height(meanDependencies) ,1) ;

% ******************  TRY SOMETHING ELSE *******************
% ######################################################################

mapkOG = crispr {:, ismember( ...
    crispr.Properties.VariableNames,oncogenes) } ;
mapkTSG = crispr {:,ismember( ...
    crispr.Properties.VariableNames,TSGs) } ;
otherOG = nonMapkCrispr {:, ismember( ...
    nonMapkCrispr.Properties.VariableNames, oncogenes) } ;
otherTSG = nonMapkCrispr{:, ismember( ...
    nonMapkCrispr.Properties.VariableNames,TSGs) } ;

% new code 
% make the dependence scores into one array for the ttest 
anovaDepedencies = [mapkOG(:) ; mapkTSG(:) ; otherOG(:) ;otherTSG(:) ] ;
anovaGroups = [ repelem({'mapkOG'}, length( mapkOG(:) ),1)  ; ...
    repelem({'mapkTSG'}, length( mapkTSG(:) ) ,1) ; ...
    repelem({'otherOG'},length( otherOG(:) ) ,1) ; ...
    repelem({'otherTSG'}, length( otherTSG(:) ) ,1) ] ;

% ######################################################################

% perform the anova test
[p,tbl,fStats] = anova1(anovaDepedencies , anovaGroups,'off') ;

% perform a multiple comparision test
multiResults = multcompare(fStats , 'Display','off','CType','bonferroni');

% save the anova results
mapkVsOtherDependAnovaResults.pValue = p;
mapkVsOtherDependAnovaResults.tbl = tbl;
mapkVsOtherDependAnovaResults.fStats = fStats;
mapkVsOtherDependAnovaResults.MultiCompare = multiResults;

% set the color to the box plots
rng(6);
color = rand(4,3) ;

% plot the data
figure()
boxplot( filloutliers(anovaDepedencies ,'linear'), anovaGroups ,...
    'Color', flipud(color) ,'Symbol','r+') ;

% set some figure properties and add title ot the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold','XTick',1:4, 'XTickLabel',...
    {'MAPK OGs','MAPK TSGs','Other OGs','Other TSGs'},'YLim',[-0.9, 1.1])
ylabel('Dependancy Score')
title('Gene Dependance','FontSize',16)

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)

% set the color of the box plots
h4 = findobj(gca,'Tag','Box') ;
for jj=1:length(h4)
    patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
        color(jj,:),'FaceAlpha', .8 ,'LineStyle','-');
end

% get the groups from the plot first after getting the
% signifcant results only
multiResults = multiResults(multiResults(:,end) < 0.05,: ) ;
groups =  multiResults(:,1:2) ;

% plot the significance value to the plot
sigstar({groups(1,:),groups(2,:),groups(3,:),groups(4,:), ...
    groups(5,:)} , multiResults(:,end))  ;

% *********************************************************************
% What cancer type has a highest MAPK dependences score
% *********************************************************************

sampleInfo.Properties.VariableNames{2} = 'cell_line';
meanDependencies = innerjoin(...
    sampleInfo(:,{'disease','cell_line'}), meanDependencies);
meanDependencies.disease = categorical(meanDependencies.disease);
meanDependencies = removevars( meanDependencies, {'cell_line'}) ;

% get the group means
cancerDependecies = grpstats(meanDependencies, 'disease', 'mean') ;

clear meanCrisprOther uniqueNames uniqueLocations ...
    meanMapk pValue tstats ii jj h4 otherCrispr cosmic cosmicOGs ...
    cosmicTSG uniprot_OGs uniprot_TSGs TSGenes ONCOgenes groups ...
    multiResults colors p tbl fStats anovaDepedencies ...
    meanDependencies anovaGroups color

%% Get the Dependency Between Each Type

% *********************************************************************
% What cancer cell lines have the highest MAPK OG and TSG score
% *********************************************************************

% add the cancer types and clean up the variables just for oncogenes
mostAffectedCCLs = innerjoin(sampleInfo(:,{'disease','cell_line'}), ...
    crisprAll(:, [true , ...
    ismember(crisprAll.Properties.VariableNames(2:end) , ...
    oncogenes) ] ) );
mostAffectedCCLs.disease = categorical(mostAffectedCCLs.disease);

% get the cancer types
cancerTypes = categories(mostAffectedCCLs.disease);

% specific the pathways
pathways = categories(mapkGenes.Pathway);

% the first loo uses all each MAPK pathway genes vs Other MAPK pathway
% genes for each cancer types Whereas the second loop uses all other
% oncogenes vs the MAPK genes to show that the MAPK are more essential that
% other genes combined
comparisonTypes = {'mapkVsmapk','mapkVsOthers'} ;
for ll = 1:length(comparisonTypes)
    
    % loop over each mapk pathway
    for kk = 1:length(pathways)
        
        % prelloacte the ttest table and avarege dependence score for all
        % the mapk alteration in different cancers
        mapkTtestDependenceScore = unique( mostAffectedCCLs(:,1) ) ;
        
        % then over cancer cell lines of each cancer type
        for ii = 1:length(cancerTypes)
            
            % find the data for the current study and other cancers study
            mapkStudyScores = mostAffectedCCLs( ...
                mostAffectedCCLs.disease == cancerTypes{ii} , :) ;
            
            % get the dependence scores for MAPK genes and the other genes
            mapkScores = nanmedian( mapkStudyScores{ :, ...
                ismember(...
                mapkStudyScores.Properties.VariableNames ,...
                mapkGenes.Gene(mapkGenes.Pathway == pathways(kk) ) ) },2);
            
            %  *********** TRY SOMETHING ELSE *****************
            mapkScores = mapkStudyScores{ :, ...
                ismember(...
                mapkStudyScores.Properties.VariableNames ,...
                mapkGenes.Gene(mapkGenes.Pathway == pathways(kk) ) ) } ;
            
            mapkScores = filloutliers( mapkScores(:) ,'linear') ;
            
            % get the other scores to use for the comparison in the test
            if strcmp(comparisonTypes(ll), 'mapkVsmapk')
                allOtherScores = nanmedian(mostAffectedCCLs{:, ...
                    ismember(mostAffectedCCLs.Properties.VariableNames ,...
                    mapkGenes.Gene(mapkGenes.Pathway == pathways(kk) ))},2);
                
                %  *********** TRY SOMETHING ELSE *****************
                allOtherScores = mostAffectedCCLs{:, ...
                    ismember(mostAffectedCCLs.Properties.VariableNames ,...
                    mapkGenes.Gene(mapkGenes.Pathway == pathways(kk) ))} ;
                                         
                allOtherScores =  filloutliers( allOtherScores(:) ,...
                    'linear' );
                
            elseif strcmp(comparisonTypes(ll), 'mapkVsOthers')
                % versus all other oncogenes within the cells
                allOtherScores = nanmedian( crisprAll{:, ...
                    [false ,false, ~ismember(...
                    crisprAll.Properties.VariableNames(3:end),...
                    mapkGenes.Gene) & ismember( ...
                    crisprAll.Properties.VariableNames(3:end), ...
                    oncogenes)]} , 2);
                          
                %  *********** TRY SOMETHING ELSE *****************
                allOtherScores = crisprAll{:, ...
                    [false ,false, ~ismember(...
                    crisprAll.Properties.VariableNames(3:end),...
                    mapkGenes.Gene) & ismember( ...
                    crisprAll.Properties.VariableNames(3:end), ...
                    oncogenes)]} ;
                
                 allOtherScores = filloutliers( allOtherScores(:) ,...
                     'linear') ;
            end
            
            % perform the ttest for the comparison for the current crisper
            % data
            [~,curPvalue,CI,stats] = ttest2( allOtherScores ,mapkScores,...
                'Vartype','unequal') ;
            
            % create a table of the drug response to the cell lines to the
            % anticancer drugs
            meanOtherScores = nanmean(allOtherScores ,'All') ;
            meanMapkScores = nanmean(mapkScores ,'All') ;
            
            % add the test results to the append
            mapkTtestDependenceScore(ii,2:7) = ...
                num2cell([meanMapkScores , meanOtherScores...
                CI(1), CI(2), stats.tstat, curPvalue]);
            
            % add the names of the values to the table
            if ii == 1
                % make the table variable different depending on the types
                % of analyis either MAPK vs MAPK or MAPK vs Others
                if strcmp(comparisonTypes(ll), 'mapkVsmapk')
                    mapkTtestDependenceScore.Properties.VariableNames(2:end) = ...
                        {'meanCancerMAPK','meanAllCancersMAPKGenes', ...
                        'lowerBound','upperBound','tValue','pValue'} ;
                else
                    mapkTtestDependenceScore.Properties.VariableNames(2:end) = ...
                        {'meanCancerMAPK','meanAllOtherGenes',...
                        'lowerBound','upperBound','tValue','pValue'} ;              
                end
                    
                % get the groups to plot
                figureGroups = [ ...
                    repelem( {'All'},length(allOtherScores), 1) ; ...
                    repelem(cancerTypes(ii),length(mapkScores), 1 ) ] ;
                figureGroups = categorical(figureGroups);
                curPlotData = [allOtherScores ; mapkScores ] ;
            else
                % get the groups again and the data to plot
                figureGroups  = [ figureGroups ; ...
                    repelem(cancerTypes(ii),length(mapkScores), 1 ) ];
                curPlotData  = [curPlotData ; mapkScores] ;
            end
            
        end
        
        % sort the table according to the level is sgnificance  and the
        % cancer cell lines with nan pvalues because they did not have
        % enough data
        mapkTtestDependenceScore  = sortrows(mapkTtestDependenceScore, ...
            'tValue','ascend') ;
%         nanCellLines = mapkTtestDependenceScore.disease( ...
%             isnan(mapkTtestDependenceScore.pValue) ) ;
%         mapkTtestDependenceScore( ...
%             isnan(mapkTtestDependenceScore.pValue),:) = [] ;
        
        % add the current pathway to the table 
        mapkTtestDependenceScore = addvars( ...
            mapkTtestDependenceScore, repmat(pathways(kk) , ...
            height(mapkTtestDependenceScore),1 ), 'After','disease',...
            'NewVariableNames','Pathway')  ;
        
        % divide the eachSample with the actual number of samples so that I
        % could properly the number of samples per cancer cell line
        % add the number of samples to the figures groups
        eachSamples = [categories(figureGroups), ...
            string(countcats(figureGroups) )] ;    
        
        % also use the crispr data to count the groups
        crisprEachSample = [ ...
            categories( categorical( mostAffectedCCLs.disease) ) ,...
            string(countcats(categorical(mostAffectedCCLs.disease))) ] ;
        
        assert( all( strcmp( eachSamples(2:end, 1) , ...
            crisprEachSample(:,1))) , 'Allignment Mismatch') 
        
        % add the groups to th each sample groups 
        crisprEachSample = [ [ "All", ...
             string(  sum( str2double(  crisprEachSample(:,2) ))) ] ; ...
              crisprEachSample ] ;
         eachSamples = crisprEachSample ;
         
        % also clean up the data so that the row with nan values are not
        % present
%         curPlotData(ismember(figureGroups,nanCellLines ) ) = [] ;
%         figureGroups(ismember(figureGroups,nanCellLines )) = [];
%         figureGroups = removecats(figureGroups, cellstr(nanCellLines));
        
        % add the number of samples to the row labels
        [~, locSamples] = ismember(cellstr( figureGroups),  ...
            eachSamples(:,1))  ;
        figureGroupsNew = strcat(cellstr( figureGroups), ...
            strcat( strcat( 'XX', eachSamples(locSamples,2) ) ,'YY') ) ;
         figureGroupsNew = categorical( ...
            strrep( strrep( figureGroupsNew, 'XX','(') ,'YY',')') ) ;
        
        % arrange the group according to the level of significance  
        % add to growing table
        if kk == 1
            cancerTypeMAPKdependence = mapkTtestDependenceScore ;
        else
            cancerTypeMAPKdependence = [ cancerTypeMAPKdependence; ...
                mapkTtestDependenceScore ];
        end
        
        % plot the data to figure
        if kk == 1
            figure();
            clf
            set(gcf,'position',[30, 120, 1200, 600]);
            hold on
        end
        % plot them in a loop
        subplot(1,4, kk)
        boxplot( filloutliers(curPlotData,'linear') ,figureGroupsNew, ...
            'Orientation','Horizontal','OutlierSize', 2);
        
        % set some figure properties and add title ot the figure
        set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
            'FontWeight','bold')
        title( string(pathways(kk)),'FontSize',11 ,'FontWeight','bold')
        
        % remove the ytick label for kk > 1
        ax = gca;
        if kk > 1
            % remove the y axis
            % ax.YAxis.Visible = 'off';
            set(gca,'YTickLabel',[])
        end
        
        % set the line width of the box plots
        set(findobj(gca,'type','line'),'linew',1)
                
        % specific the color of the plots using the p values of the test
        plotColors = [0.1 0.9 0.1;  0.1940 0.3840 0.9560 ; ...
            0.9900 0.2250 0.0980; 0.7 0.7 0.7] ;
        
        % score add color to the box plots get the unique groups
        uniqueGroups = flipud( categories(figureGroups) );
        h4 = findobj(gca,'Tag','Box') ;
        for jj = 1:length(h4)
            
            % get the pvalue  and t values
            curPvalue = mapkTtestDependenceScore.pValue( ....
                mapkTtestDependenceScore.disease == uniqueGroups(jj) ) ;
            curTvalue = mapkTtestDependenceScore.tValue( ....
                mapkTtestDependenceScore.disease == uniqueGroups(jj) ) ;
            
            % seletect the colors
            if strcmp( uniqueGroups(jj) , 'All')
                color = plotColors(1,:)  ;
            elseif curTvalue > 0 && curPvalue < 0.05
                % for cancer cell that dependent on MAPK signalling
                color = plotColors(2,:) ;
            elseif curTvalue < 0 && curPvalue < 0.05
                % for cancer cell that dependent are not MAPK signalling         
                color = plotColors(3,:) ;
            else
                % for non significant comparisions
                color = plotColors(4,:) ;
            end
            patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
                color,'FaceAlpha', 0.8) % 'LineWidth',1);
        end
        
    end
    
    % add a super title
    if strcmp(comparisonTypes(ll), 'mapkVsmapk')
        % sgtitle('\bf MAPK Genes vs all MAPK Genes','FontSize',16)
        annotation('textbox',[0.40 0.93 0.57 0.08],...
        'String', 'Between Cancer Type MAPK pathways Depedency',...
        'FitBoxToText','on','FontWeight','bold','FontSize',16,...
        'EdgeColor','none')
    elseif strcmp(comparisonTypes(ll), 'mapkVsOthers')
        % sgtitle('\bf MAPK Genes vs All Other Genes','FontSize',16)
        annotation('textbox',[0.40 0.93 0.57 0.08],...
        'String','MAPK Genes vs All Other Genes',...
        'FitBoxToText','on','FontWeight','bold','FontSize',16,...
        'EdgeColor','none')
    end
    
    annotation('textbox',[0.46 0.017 0.57 0.06],...
        'String','Dependancy Score',...
        'FitBoxToText','on','FontWeight','bold','FontSize',12,...
        'EdgeColor','none')
    
    % add a legend to figure  using the create legend function
    % First specify the legend labels and yPoint and xStart
    legendLabels = {'All Cell Lines';'Increased Dependency';...
        'Reduced Dependency'; 'Neutral'} ;
    yPoint = 0.85 ; xStart = 0.89 ; myLgdTitle = 'Cell Fitness';
    createLegendInternal(yPoint, xStart, legendLabels , plotColors ,myLgdTitle)

    hold off % of the figure
    
    % save the data to supplementary files
    if strcmp(comparisonTypes(ll), 'mapkVsmapk')
%         cancerTypeMAPKdependence = sortrows(cancerTypeMAPKdependence ,...
%             'tValue','ascend') ;
        writetable(cancerTypeMAPKdependence,'Supplementary File 3.xlsx',...
            'Sheet','MAPK vs MAPK Genes Dependence')
    elseif strcmp(comparisonTypes(ll), 'mapkVsOthers')   
%         cancerTypeMAPKdependence = sortrows(cancerTypeMAPKdependence ,...
%             'tValue','ascend') ;
        writetable(cancerTypeMAPKdependence,'Supplementary File 3.xlsx',...
            'Sheet','MAPK vs Other Genes Dependence') ;     
    end
    
end

clear eachSamples  locSamples figureGroupsNew legendLabels  yPoint ...
    xStart legendLabels plotColors myLgdTitle crisprEachSample

%% ====== do this for the protein classess of the MAPK

% add the cancer types and clean up the variables just for oncogenes
mostAffectedCCLs = innerjoin(...
    sampleInfo(:,{'disease','cell_line'}),crisprAll );
mostAffectedCCLs.disease = categorical(mostAffectedCCLs.disease);

% get the cancer types
cancerTypes = categories(mostAffectedCCLs.disease);

% the first loo uses all each MAPK pathway genes vs Other MAPK pathway
% genes for each cancer types Whereas the second loop uses all other
% oncogenes vs the MAPK genes to show that the MAPK are more essential that
% other genes combined
comparisonTypes = {'mapkVsmapk','mapkVsOthers'} ;
for ll = 1:length(comparisonTypes)
    
    % loop over each mapk pathway
    for kk = 1:length(proteinClass)
        
        % prelloacte the ttest table and avarege dependence score for all
        % the mapk alteration in different cancers
        mapkTtestDependenceScore = unique( mostAffectedCCLs(:,1) ) ;
        
        % then over cancer cell lines of each cancer type
        for ii = 1:length(cancerTypes)
            
            % find the data for the current study and other cancers study
            mapkStudyScores = mostAffectedCCLs( ...
                mostAffectedCCLs.disease == cancerTypes{ii} , :) ;
            
            % get the dependence scores for MAPK genes and the other genes
            mapkScores =  mapkStudyScores{ :, ...
                ismember(...
                mapkStudyScores.Properties.VariableNames ,...
                mapkGenes.Gene( ...
                mapkGenes.proteinClass == proteinClass(kk))) } ;
             mapkScores = fillmissing( ...
                 filloutliers( mapkScores(:) ,'linear') ,'linear') ;
            
            % get the other scores to use for the comparison in the test
            if strcmp(comparisonTypes(ll), 'mapkVsmapk')
                allOtherScores = mostAffectedCCLs{:, ...
                    ismember(mapkStudyScores.Properties.VariableNames ,...
                    mapkGenes.Gene( ...
                    mapkGenes.proteinClass == proteinClass(kk) ))} ;
                allOtherScores = fillmissing( ...
                    filloutliers(allOtherScores(:), 'linear') ,'linear');
                
            elseif strcmp(comparisonTypes(ll), 'mapkVsOthers')
                % versus all other oncogenes within the cells
                allOtherScores =  crisprAll{:, ...
                    [false ,false, ~ismember(...
                    crisprAll.Properties.VariableNames(3:end),...
                    mapkGenes.Gene) & ismember( ...
                    crisprAll.Properties.VariableNames(3:end),...
                    oncogenes)]};
                allOtherScores = fillmissing( ...
                    filloutliers(allOtherScores(:), 'linear') ,'linear');
                
%                     allOtherScores = nanmedian( crisprAll{:, ...
%                     [false ,false, ~ismember(...
%                     crisprAll.Properties.VariableNames(3:end),...
%                     mapkGenes.Gene)]} , 2);
            end
            
            % perform the ttest for the comparison for the current crisper
            % data
            [~,curPvalue,CI,stats] = ttest2( allOtherScores ,mapkScores,...
                'Vartype','unequal') ;
            
            % create a table of the drug response to the cell lines to the
            % anticancer drugs
            meanOtherScores = nanmean(allOtherScores ,'All') ;
            meanMapkScores = nanmean(mapkScores ,'All') ;
            
            % add the test results to the append
            mapkTtestDependenceScore(ii,2:7) = ...
                num2cell([meanMapkScores , meanOtherScores...
                CI(1), CI(2), stats.tstat, curPvalue]);
                   
            % add the names of the values to the table
            if ii == 1
                % set different variable names
                if strcmp( comparisonTypes{ll}, 'mapkVsmapk' )
                    mapkTtestDependenceScore.Properties.VariableNames(2:end) = ...
                        {'meanCancerMAPK','meanAllCancersMAPKGenes', ...
                        'lowerBound','upperBound','tValue','pValue'} ;
                elseif strcmp( comparisonTypes{ll}, 'mapkVsOthers' )
                    mapkTtestDependenceScore.Properties.VariableNames(2:end) = ...
                        {'meanCancerMAPK','meanAllOtherGenes',...
                        'lowerBound','upperBound','tValue','pValue'} ;
                end
                                    
                % get the groups to plot
                figureGroups = [ ...
                    repelem( {'All'},length(allOtherScores), 1) ; ...
                    repelem(cancerTypes(ii),length(mapkScores), 1 ) ] ;
                figureGroups = categorical(figureGroups);
                curPlotData = [allOtherScores ; mapkScores ] ;
            else
                % get the groups again and the data to plot
                figureGroups  = [ figureGroups ; ...
                    repelem(cancerTypes(ii),length(mapkScores), 1 ) ];
                curPlotData  = [curPlotData ; mapkScores] ;
            end
            
        end
        
        % sort the table according to the level is sgnificance  and the
        % cancer cell lines with nan pvalues because they did not have
        % enough data
        mapkTtestDependenceScore  = sortrows(mapkTtestDependenceScore, ...
            'tValue','ascend') ;
        
        % add the current pathway to the table
        mapkTtestDependenceScore = addvars( ...
            mapkTtestDependenceScore, repmat(proteinClass(kk) , ...
            height(mapkTtestDependenceScore),1 ), 'After','disease',...
            'NewVariableNames','Pathway')  ;
        
%         nanCellLines = mapkTtestDependenceScore.disease( ...
%             isnan(mapkTtestDependenceScore.pValue) ) ;
        mapkTtestDependenceScore( ...
            isnan(mapkTtestDependenceScore.pValue),:) = [] ;
        
        % also clean up the data so that the row with nan values are not
        % present
%         curPlotData(ismember(figureGroups,nanCellLines ) ) = [] ;
%         figureGroups(ismember(figureGroups,nanCellLines )) = [];
%         figureGroups = removecats(figureGroups, cellstr(nanCellLines));

        % divide the eachSample with the actual number of samples so that I
        % could properly the number of samples per cancer cell line
        % add the number of samples to the figures groups
        eachSamples = [categories(figureGroups), ...
            string(countcats(figureGroups) )] ;    
        
        % also use the crispr data to count the groups
        crisprEachSample = [ ...
            categories( categorical( mostAffectedCCLs.disease) ) ,...
            string(countcats(categorical(mostAffectedCCLs.disease))) ] ;
        
        assert( all( strcmp( eachSamples(2:end, 1) , ...
            crisprEachSample(:,1))) , 'Allignment Mismatch')
        
        % add the groups to th each sample groups
        crisprEachSample = [ [ "All", ...
            string(  sum( str2double(  crisprEachSample(:,2) ))) ] ; ...
            crisprEachSample ] ;
        eachSamples = crisprEachSample ;
        
        % add the number of samples to the row labels
        [~, locSamples] = ismember(cellstr( figureGroups),  ...
            eachSamples(:,1))  ;
        figureGroupsNew = strcat(cellstr( figureGroups), ...
            strcat( strcat( 'XX', eachSamples(locSamples,2) ) ,'YY') ) ;
        figureGroupsNew = categorical( ...
            strrep( strrep( figureGroupsNew, 'XX','(') ,'YY',')') ) ;
        
        % add to growing table
        if kk == 1
            cancerTypeProteinClassdependence = mapkTtestDependenceScore ;
        else
            cancerTypeProteinClassdependence = ....
                [ cancerTypeProteinClassdependence; ...
                mapkTtestDependenceScore ];
        end
        
        % plot the data to figure
        if kk == 1
            figure();
            clf
            set(gcf,'position',[30, 120, 1200, 600]);
            hold on
        end
        % plot them in a loop
        subplot(2,5,kk)
        boxplot( filloutliers(curPlotData,'linear'), figureGroupsNew ,...
            'Orientation','Horizontal','OutlierSize', 2) ;
        
        % set some figure properties and add title ot the figure
        set(gca,'LineWidth',1,'FontSize',8 ,'Box','off',...
            'FontWeight','bold')
        title( string(proteinClass(kk)),'FontSize',9)
        
        % remove the ytick label for kk > 1
        ax = gca;
        if kk ~= 1
            if kk ~= 6
                % remove the y axis
                % ax.YAxis.Visible = 'off';
                set(gca,'YTickLabel',[])
            end
        end
        
        % set the line width of the box plots
        set(findobj(gca,'type','line'),'linew',1)
        
        % add a reference line to the plot refEnd = ylim; line([0,0],[0,
        % refEnd(2)],'LineStyle','--','Color','k','LineWidth',1)
        
        % specific the color of the plots using the p values of the test
        plotColors = [0.1 0.9 0.1;  0.1940 0.3840 0.9560 ; ...
            0.9900 0.2250 0.0980; 0.7 0.7 0.7] ;
        
        % specific the color of the plots using the p values of the test
        % score add color to the box plots get the unique groups
        uniqueGroups = flipud( categories(figureGroups ) );
        h4 = findobj(gca,'Tag','Box') ;
        for jj = 1:length(h4)
            
            % get the pvalue  and t values
            curPvalue = mapkTtestDependenceScore.pValue( ....
                mapkTtestDependenceScore.disease == uniqueGroups(jj) ) ;
            curTvalue = mapkTtestDependenceScore.tValue( ....
                mapkTtestDependenceScore.disease == uniqueGroups(jj) ) ;
            
            % seletect the colors
            if strcmp( uniqueGroups(jj) , 'All')
                color = plotColors(1,:)  ;
            elseif curTvalue > 0 && curPvalue < 0.05
                % for cancer cell that dependent on MAPK signalling
                color = plotColors(2,:) ;
            elseif curTvalue < 0 && curPvalue < 0.05
                % for cancer cell that dependent are not MAPK signalling         
                color = plotColors(3,:) ;
            else
                % for non significant comparisions
                color = plotColors(4,:) ;
            end
            patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
                color,'FaceAlpha', 0.8 ) % ,'LineWidth',1;
        end
        
    end
    
    % add a super title
    if strcmp(comparisonTypes(ll), 'mapkVsmapk')
        % sgtitle('\bf MAPK Genes vs all MAPK Genes','FontSize',16)
        annotation('textbox',[0.40 0.93 0.57 0.08],...
        'String', 'Between Cancer Type MAPK Pathway Gene Class Depedency',...
        'FitBoxToText','on','FontWeight','bold','FontSize',16,...
        'EdgeColor','none')
    elseif strcmp(comparisonTypes(ll), 'mapkVsOthers')
        % sgtitle('\bf MAPK Genes vs All Other Genes','FontSize',16)
        annotation('textbox',[0.40 0.93 0.57 0.08],...
        'String','MAPK Pathway Genes Classess vs All Other Genes',...
        'FitBoxToText','on','FontWeight','bold','FontSize',16,...
        'EdgeColor','none')
    end
    
    annotation('textbox',[0.46 0.017 0.57 0.06],...
        'String','Dependancy Score',...
        'FitBoxToText','on','FontWeight','bold','FontSize',12,...
        'EdgeColor','none')
    
    % add a legend to figure  using the create legend function
    % First specify the legend labels and yPoint and xStart
    legendLabels = {'All Cell Lines';'Dependent';...
        'non-Dependent'; 'Neutral'} ;
    yPoint = 0.85 ; xStart = 0.89 ; myLgdTitle = 'Cell Fitness';
    createLegendInternal(yPoint, xStart, legendLabels , plotColors ,myLgdTitle)

    hold off % of the figure
    
    % save the data to supplementary files
    if strcmp(comparisonTypes(ll), 'mapkVsmapk')
        writetable(cancerTypeProteinClassdependence, ...
            'Supplementary File 3.xlsx',...
            'Sheet','Gene Class- MAPK vs MAPK Dep')
    elseif strcmp(comparisonTypes(ll), 'mapkVsOthers')   
        writetable(cancerTypeProteinClassdependence, ...
            'Supplementary File 3.xlsx',...
            'Sheet','Gene Class- MAPK vs Other Genes') ;     
    end
    
end

clear h4 color plotData figureGroups ii meanOtherScores CI stats ...
    curPvalue allOtherScores mapkStudyScores mapkScores cancerTypes ...
    mostAffectedCCLs nanCellLines meanOtherScores  meanMapkScores ...
    uniqueGroups mapkTtestDependenceScore ii jj kk curTvalue ...
    curPlotData comparisonTypes nonMapkCrispr

%% Correlation Between Cell Line Dependency Scores

% *********************************************************************
% Do cell lines of the same cancer type have similar dependency scores
% *********************************************************************

% get the clustering data and the missing data
clusterDataDep = innerjoin(sampleInfo(:,{'disease','cell_line'}),crispr );
clusterDataDep = fillmissing(clusterDataDep ,'linear','DataVariables',...
    clusterDataDep.Properties.VariableNames(3:end)) ;

% plot the clutergram of dependency scores for tumour types
% corData = corr(clusterDataDep{:,3:end});
cgo = clustergram( clusterDataDep{:,3:end},...
    'columnlabels',clusterDataDep.Properties.VariableNames(3:end) ,...
    'rowlabels', clusterDataDep.cell_line,...
    'colormap',redbluecmap ,'standardize','row',...
    'ColumnPDist','cosine','Dendrogram',6 ,'Linkage','average');

% % get a copy of the clustergram object and plot it for a sigle clustergram
% cgoCopy = cgo ;
% % get the colour bar on the plot
% cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
%     'Color',{[1 1 0],[0.6 0.6 1]});
% set(cgoCopy,'ColumnGroupMarker',cm ,'RowLabels',clusterDataDep.cell_line);
%
% % make the figure very big to include the full row labels
% cgAxes = plot(cgoCopy);
% set(cgAxes) ;

% get the mapk pathway to which the genes belong and add that to op of the
% clustergram

% ===== get the location of the data to produce the multiple plots ====
% there is a small bug where the row label seem to change by themselves
[~, locCols] = ismember(cgo.ColumnLabels, ...
    clusterDataDep.Properties.VariableNames(3:end) ) ;
[~, locRows] = ismember(cgo.RowLabels, clusterDataDep.cell_line) ;
heatData = clusterDataDep{:,3:end} ;
heatData =  heatData(locRows,locCols);

%% ==================== Produce Multiple Plots ====================
figure();
% clf
% set(gcf,'position',[100,50,1000,600]);
% the first number is how far the figure will be from the x-axis and the
% seceond number is now far the figure will be from the y-axis. The third
% number is was far the figure will run across the figure bar and the last
% number is far it will displaced allow the y-axis

axes('position',[0.20, 0.15, 0.74, 0.56]);
heatmap(cgo.ColumnLabels, flip(cgo.RowLabels), heatData,...
    'Colormap',redbluecmap,'ColorbarVisible','off',...
    'GridVisible','off' ,'ColorLimits',[-1 0.8] ,'FontColor','k' ,...
    'FontSize',8 , ...
    'XLabel',sprintf('%d MAPK Pathway Genes', width(clusterDataDep) ) ,...
    'YLabel', sprintf('%d Cell Lines', height(clusterDataDep) ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET THE COLORS %%%%%%%%%%%%%%%%%%%%%%%%%%
%  the JCO palette
classesColors.proteinClass = [
    0                   0.450980392156863   0.760784313725490
    0.937254901960784   0.752941176470588                   0
    0.525490196078431   0.525490196078431   0.525490196078431
    0.803921568627451   0.325490196078431   0.298039215686275
    0.478431372549020   0.650980392156863   0.862745098039216
    0                   0.235294117647059   0.403921568627451
    0.560784313725490   0.466666666666667                   0
    0.231372549019608   0.231372549019608   0.231372549019608
    0.654901960784314   0.188235294117647   0.188235294117647
    0.290196078431373   0.411764705882353   0.564705882352941];
classesColors.proteinClass = classesColors.proteinClass(...
    1:numel(unique(mapkGenes.proteinClass)) , :) ;

% the nature palette
classesColors.pathway = [
    0.901960784313726   0.294117647058824   0.207843137254902
    0.301960784313725   0.733333333333333   0.835294117647059
    0                   0.627450980392157   0.529411764705882
    0.235294117647059   0.329411764705882   0.533333333333333
    0.952941176470588   0.607843137254902   0.498039215686275
    0.517647058823529   0.568627450980392   0.705882352941177
    0.568627450980392   0.819607843137255   0.760784313725490
    0.862745098039216                   0                   0
    0.494117647058824   0.380392156862745   0.282352941176471
    0.690196078431373   0.611764705882353   0.521568627450980
    ];

classesColors.pathway = classesColors.pathway( ...
    1:numel(unique(mapkGenes.Pathway)) ,:) ;

% the science color palette
classColors.SciencePalette = [
    0.231372549019608   0.286274509803922   0.572549019607843
    0.933333333333333                   0                   0
    0                   0.545098039215686   0.270588235294118
    0.388235294117647   0.094117647058824   0.474509803921569
    0.372549019607843   0.333333333333333   0.607843137254902
    0.733333333333333                   0   0.129411764705882
    0                   0.509803921568627   0.501960784313725
    0.635294117647059                   0   0.337254901960784
    0.501960784313725   0.505882352941176   0.501960784313725
    0.105882352941176   0.098039215686275   0.098039215686275
    ];

classesColors.proteinClass = classColors.SciencePalette ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the alterations rate of mutations and copy number changees in each
% tumour type
% change the arrangement in the alterationRate
% first get the MAPK genes that are present in the DeepMap database
[~, locGenes] = ismember( cgo.ColumnLabels ,mapkGenes.Gene);
mapkGenes = mapkGenes(locGenes, :) ;

% throw in an assertion
assert(all( strcmp(cgo.ColumnLabels' ,mapkGenes.Gene)))

% add a heatmap on top to show the supervised classifier output
axes('position',[0.20, 0.715, 0.74,0.025]);
ultraBars( double(mapkGenes.Pathway)', classesColors.pathway ) ;
% add the name of the genes to the left of heatmap
dim = [0.1, 0.715,0.1,0.025];
annotation('textbox',dim,'String','Pathway',...
    'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
    'HorizontalAlignment','right','FontWeight','bold',...
    'VerticalAlignment','middle');
hold off

% add the class the protein to the chart
axes('position',[0.20, 0.742, 0.74,0.025]);
ultraBars(double(mapkGenes.proteinClass)', classesColors.proteinClass );
% add the name of the bar to the left
dim = [0.1, 0.742, 0.1,0.025];
annotation('textbox',dim,'String','MAPK Gene Class',...
    'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
    'HorizontalAlignment','right','FontWeight','bold',...
    'VerticalAlignment','middle');
hold off

% % now get the data of the cancer study and change the arrangement in the
% cancer studies
barData = clusterDataDep(:,3:end) ;
[~, locCols] = ismember(cgo.ColumnLabels, barData.Properties.VariableNames );
barData = barData(:,locCols) ;

% throw in an assertion
assert(all(strcmp(cgo.ColumnLabels',barData.Properties.VariableNames')))

% get the values above a give threshold of the essentialies scores
barData = table2array(barData) <  -0.50;
barData = sum(barData, 1);

axes('position',[0.2, 0.77, 0.74, 0.025]);
heatmap(barData,'Colormap', parula ,'GridVisible','off' ,...
    'FontColor','none','ColorbarVisible','off','ColorLimits',[-1 75]);
% add the name of the bar to the left
dim = [0.1, 0.77, 0.1, 0.025];
annotation('textbox',dim,'String','Dependecy Score',...
    'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
    'HorizontalAlignment','right','FontWeight','bold',...
    'VerticalAlignment','middle');

% add the legend to the figure
% add a legend to figure  using the create legend function
% First specify the legend labels and yPoint and xStart
legendLabels = categories( mapkGenes.Pathway) ;

% specify where to plot the legend and the font size and the width of the
% rectangle and the text box
yPoint = 0.70 ; xStart = 0.05 ; myLgdTitle = 'Pathway'; fontSizes = [8,10];
rectAndTextBox = [0.01 ,0.07] ;

% create the legend
createLegendInternal(yPoint, xStart, legendLabels , classesColors.pathway , ...
    myLgdTitle , fontSizes ,rectAndTextBox)

% add the legend to the figure
% add a legend to figure  using the create legend function
% First specify the legend labels and yPoint and xStart
legendLabels = categories( mapkGenes.proteinClass) ;
% replace Guanosine Exchange Factor with GEF and Dual Specificy Phosphatase
% withh DUSP
legendLabels = strrep(legendLabels, 'DualSpecificityPhosphatase','DUSP');
legendLabels = strrep(legendLabels, 'GuanosineExchangeFactor','GEF');

% specify where to plot the legend and the font size and the width of the
% rectangle and the text box
yPoint = 0.55 ; xStart = 0.05 ; myLgdTitle = 'MAPK Gene Class';
fontSizes = [8,10]; rectAndTextBox = [0.01 ,0.07] ;

% create the legend
createLegendInternal(yPoint, xStart, legendLabels , classesColors.proteinClass,...
    myLgdTitle ,fontSizes ,rectAndTextBox)

% replot the heatmap because there seems to be a problem
axes('position',[0.20, 0.15, 0.74, 0.56]);
heatmap(cgo.ColumnLabels, flip(cgo.RowLabels), heatData,...
    'Colormap',redbluecmap,'ColorbarVisible','off',...
    'GridVisible','off' ,'ColorLimits',[-1 0.8] ,'FontColor','k' ,...
    'FontSize',8 , ...
    'XLabel',sprintf('%d MAPK Pathway Genes', width(clusterDataDep) ) ,...
    'YLabel', sprintf('%d Cell Lines', height(clusterDataDep) ) );

clear locX locY colorM cgAxes cgoCopy clusterGeneDep nonMapkCrisp ...
    fontSizes rectAndTextBox myLgdTitle xStart yPoint legendLabels ...
    ii ll locCols locGenes LocThem headData dim X these ans

%% Correlation betweeen Essentiality and Other data

fprintf('\n Loading methylation data for the cell lines \n')
% load the methylation data and clean it up
CCLEmethyl = readtable('CClE_methylation_processed.csv');

% change the variable names and remove the methyl data with TT because
% these are duplicates
CCLEmethyl(:, contains(CCLEmethyl.Properties.VariableNames,'TT_')) = [] ;
CCLEmethyl.Properties.VariableNames = ...
    regexprep(CCLEmethyl.Properties.VariableNames, '\_+\w*','') ;

for ii = 2:width(CCLEmethyl)
    if iscell( CCLEmethyl.(ii))
        CCLEmethyl.(ii) = str2double(CCLEmethyl.(ii)) ;
    end
end

% get the avarage of methylaiton and transpose the methylation table
% CCLEmethylAverage = CCLEmethyl(:, 1:2);
CCLEmethyl = removevars(CCLEmethyl,'avg');
CCLEmethyl.cluster = categorical(CCLEmethyl.cluster);
CCLEmethyl = groupsummary(CCLEmethyl ,'cluster','mean') ;
CCLEmethyl = removevars(CCLEmethyl, 'GroupCount');
CCLEmethyl = rows2vars( CCLEmethyl ,'VariableNamesSource','cluster');
CCLEmethyl.Properties.VariableNames(1) = "cell_line" ;
CCLEmethyl.cell_line = extractAfter(CCLEmethyl.cell_line,'_');

% =========== load and process the CCLE mutation data ===============
fprintf('\n Processing CCLE mutations data \n')
CCLEmutations = readtable('CCLE_mutations.csv');

% return only the MAPK genes so that processMAF files runs faster
CCLEmutations = CCLEmutations( ismember( ...
    CCLEmutations.Hugo_Symbol , mapkGenes.Gene) , :) ;

% the data turn out to be a mutation annotation file that needs to be
% processed into a table 
CCLEmutations = processMAF(CCLEmutations) ;

% change the DepMap Cell Ids to the cell line IDs
CCLEmutations.Properties.VariableNames(1) = "DepMap_ID" ;
CCLEmutations = innerjoin(sampleInfo(:,[1,2]),CCLEmutations );

% get only the cell line that are present in the CRISPR data 
CCLEmutations = CCLEmutations(...
    ismember(CCLEmutations.cell_line, crispr.cell_line) , :) ;

% make the methylation data the sample length as the crispr data
% get the missing genes and add them to the methylation data using NaN
these = ismember(crispr.Properties.VariableNames , ...
    CCLEmutations.Properties.VariableNames) ;
missingGenes = crispr.Properties.VariableNames(~these);

% add the NaN
theMissing = array2table( ...
    NaN(height(CCLEmutations), length(missingGenes) ) ,...
    'VariableNames', missingGenes) ;
CCLEmutations = [ CCLEmutations , theMissing ] ;

% arrange the mutations data in the same way order as the samples cluster
[~, LocThem ] = ismember(cgo.ColumnLabels , ...
    CCLEmutations.Properties.VariableNames) ;
LocThem(LocThem == 0) = [] ;
CCLEmutations = [CCLEmutations(:,[1,2]), CCLEmutations(:,LocThem)] ;

% throw in an assertion
assert( all( strcmp( CCLEmutations.Properties.VariableNames(3:end)',...
    cgo.ColumnLabels') ) )

% ============== Load and process the copy number data ==============
fprintf('\n Processing CCLE copy number data \n')

CCLEcna = readtable('CCLE_gene_cn.csv') ;

% change the variable names in the table 
varNames = matlab.lang.makeUniqueStrings( [{'DepMap_ID'}, ...
    extractBefore(CCLEcna.Properties.VariableNames(2:end), '_') ] ); 
CCLEcna.Properties.VariableNames = varNames  ;
 
% add the covertional sample names to the table and return only the cell
% line that have crispr data
CCLEcna = innerjoin(sampleInfo(:,[1,2]),CCLEcna );
CCLEcna = CCLEcna( ...
    ismember(CCLEcna.cell_line, crispr.cell_line), :)  ;

% arrange the mutations data in the same way order as the samples cluster
[~, LocThem ] = ismember(cgo.ColumnLabels , ...
    CCLEcna.Properties.VariableNames) ;
LocThem(LocThem == 0) = [] ;
CCLEcna = [CCLEcna(:,[1,2]), CCLEcna(:,LocThem)] ;

% throw in an assertion
assert( all( strcmp( CCLEcna.Properties.VariableNames(3:end)',...
    cgo.ColumnLabels') ) )

% ================= let now also load the mrna data ==================
fprintf('\n Loading mrna trascription data for the cell lines \n')
if ~exist('CCLE_mrna_processed.csv','file')
    CCLEmrna = readtable('CCLE_expression.csv');
    
    % clean up the data
    [~, these] =  unique(regexprep(CCLEmrna.Properties.VariableNames,...
        '\_+\w*','') ,'stable');
    CCLEmrna = CCLEmrna(:, these) ;
    CCLEmrna.Properties.VariableNames = ...
        regexprep(CCLEmrna.Properties.VariableNames, '\_+\w*','') ;
    
    % save a copy of the other mrnas
    CCLEothermrna = CCLEmrna(:, ...
        ~ismember(CCLEmrna.Properties.VariableNames, mapkGenes.Gene ))  ;
    
    % return only the MAPK genes
    CCLEmrna = [ CCLEmrna(:,1), CCLEmrna(:, ...
        ismember(CCLEmrna.Properties.VariableNames,mapkGenes.Gene )) ] ;
    
    % change the cell line ID to those of the common names from the DepMap
    % Ids
    CCLEmrna.Properties.VariableNames(1) = ...
        sampleInfo.Properties.VariableNames(1);
    CCLEmrna = innerjoin( sampleInfo(:,1:2), CCLEmrna) ;
    CCLEmrna = CCLEmrna(:,2:end) ;
    CCLEmrna.Properties.VariableNames(1) = "cell_line" ;
    
    % do the same for all the mrna data
    CCLEothermrna.Properties.VariableNames(1) = ...
        sampleInfo.Properties.VariableNames(1);
    CCLEothermrna = innerjoin( sampleInfo(:,1:2), CCLEothermrna) ;
    CCLEothermrna = CCLEothermrna(:,2:end) ;
    CCLEothermrna.Properties.VariableNames(1) = "cell_line" ;
    
    % save to excel
    writetable(CCLEmrna,'CCLE_mrna_processed.csv');
    writetable(CCLEothermrna,'CCLE_mrna_nonMAPK.csv')
else
    % load the mrna data that was processed
    CCLEmrna = readtable('CCLE_mrna_processed.csv') ;
end

% make the mrna and methylation data the same length and order as the
% crispr data
[~, LocThem ] = ismember(cgo.ColumnLabels  , ...
    CCLEmrna.Properties.VariableNames) ;
LocThem(LocThem == 0) = [] ;
CCLEmrna = [CCLEmrna(:,1), CCLEmrna(:, LocThem) ] ;

% also get only the cell line that are present in the crispr
[~, LocThem ] = ismember(flip(cgo.RowLabels) , CCLEmrna.cell_line) ;
LocThem(LocThem == 0) = [] ;
CCLEmrna = CCLEmrna(LocThem,:) ;

% throw in an assertion
assert( all( strcmp( CCLEmrna.Properties.VariableNames(2:end)',...
    cgo.ColumnLabels') ) )

% make the methylation data the sample length as the crispr data
% get the missing genes and add them to the methylation data using NaN
these = ismember(crispr.Properties.VariableNames , ...
    CCLEmethyl.Properties.VariableNames) ;
missingGenes = crispr.Properties.VariableNames(~these);

% add the NaN
theMissing = array2table( ...
    NaN(height(CCLEmethyl), length(missingGenes) ) ,...
    'VariableNames', missingGenes) ;
CCLEmethyl = [ CCLEmethyl , theMissing ] ;

% make them the same length as the other genes
[~, LocThem ] = ismember(cgo.ColumnLabels , ...
    CCLEmethyl.Properties.VariableNames) ;
LocThem(LocThem == 0) = [] ;
CCLEmethyl = [CCLEmethyl(:,1), CCLEmethyl(:, LocThem) ];

% also get only the cell line that are present in the crispr
[~, LocThem ] = ismember(flip(cgo.RowLabels) , CCLEmethyl.cell_line) ;
LocThem(LocThem == 0) = [] ;
CCLEmethyl = CCLEmethyl(LocThem,:) ;

% throw in an assertion
assert( all( strcmp( CCLEmethyl.Properties.VariableNames(2:end)',...
    cgo.ColumnLabels') ) )

% ==================== Produce Multiple Plots ====================
figure();
clf
set(gcf,'position',[100,50,1000,600]);
% the first number is how far the figure will be from the x-axis and the
% seceond number is now far the figure will be from the y-axis. The third
% number is was far the figure will run across the figure bar and the last
% number is far it will displaced allow the y-axis

axes('position',[0.15, 0.10, 0.74, 0.20]);
heatmap(cgo.ColumnLabels, flip(cgo.RowLabels), heatData,...
    'Colormap',redbluecmap,'ColorbarVisible','off',...
    'GridVisible','off' ,'ColorLimits',[-1 0.8] ,'FontSize',8 );

mrnaHeatData = filloutliers(CCLEmrna{:,2:end} ,'linear');
axes('position',[0.15, 0.31, 0.74, 0.20]);
heatmap(mrnaHeatData ,'Colormap',redgreencmap,'GridVisible','off' ,...
    'FontColor','none','ColorbarVisible','off','ColorLimits',[-10 20]);

axes('position',[0.15, 0.52, 0.74, 0.20]);
heatmap(CCLEmethyl{:,2:end},'Colormap',redbluecmap,'GridVisible','off',...
    'FontColor','none','ColorbarVisible','off' ,...
    'MissingDataColor',[0.8 0.8 0.8]);

axes('position',[0.15, 0.73, 0.74, 0.20]);
colors = redblue(20) ;
heatmap( zscore(CCLEcna{:,3:end}) ,'Colormap', colors ,...
    'GridVisible','off','FontColor','none','ColorbarVisible','off' ,...
    'MissingDataColor',[0.8 0.8 0.8] ,'ColorLimits',[-1 1]);

% axes('position',[0.15, 0.65, 0.74, 0.13]);
% CCLEmutationNum = double( ~ismissing( CCLEmutations(:, 3:end) ) ) ;
% heatmap( CCLEmutationNum ,'Colormap',redgreencmap,'GridVisible','off',...
%     'FontColor','none','ColorbarVisible','off' ,...
%     'MissingDataColor',[0.8 0.8 0.8] ,'ColorLimits',[-1 1]);

% create a textbox
% Create textbox
textToPut = {sprintf('CRISPR Fitness: %d Cell lines', height(crispr) ), ...
    sprintf('mRNA: %d Cell Lines', height(CCLEmrna) ), ...
    sprintf('DNA methylation: %d Cell Lines', height(CCLEmethyl) ), ...
    sprintf('CNV: %d Cell Lines', height(CCLEcna) ) } ;
textPos = [0.20, 0.40, 0.60, 0.80] ;
for ii = 1:length(textToPut)
    annotation('textbox',[0.052 textPos(ii) 0.1 0.0355],...
        'VerticalAlignment','middle','String', textToPut{ii} , ...
        'FontWeight','bold','FontSize',12,'FontName','Helvetica Neue',...
        'EdgeColor','none','HorizontalAlignment','right');
end

% Create textbox
annotation('textbox',...
    [0.285678733 0.94205734665 0.62879638009 0.04564907097],...
    'String',{['Relation Between Gene Essentialiy and Other',...
    ' Molecular Profiling Dataset' ]},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'EdgeColor',[1 1 1]);

% *********************************************************************
% Is there a correlation between the gene expression signature and the
% gene dependence scores
% *********************************************************************
% get the row-wise correlation of between dependence scores and the
% expression of the genes

% make the data the same size
[~, locRows] = ismember(CCLEmrna.cell_line,clusterDataDep.cell_line)  ;
[~, locColumns] = ismember( ...
    CCLEmrna.Properties.VariableNames , ...
    clusterDataDep.Properties.VariableNames );

% remove the 0s for the indexing arrays and the get the data
locRows(locRows == 0 ) = [];
locColumns(locColumns == 0) = [] ;
corrDataDep = clusterDataDep(locRows,:) ;
corrDataDep  = [corrDataDep(:,1) , corrDataDep(:, locColumns) ]  ;

% throw in an assertion
assert( all (strcmp(corrDataDep.cell_line , CCLEmrna.cell_line)))
assert(all (strcmp ( corrDataDep.Properties.VariableNames(2:end)' ,...
    CCLEmrna.Properties.VariableNames') ) )

% now get the data plot the correlation
matchedDepData = corrDataDep(:, 2:end);
corrDataDep = corrDataDep(:, 3:end) ;
geneEssentalityValue = sum( corrDataDep{:,:} < -0.5  );
[~, locMax] = maxk( geneEssentalityValue , 1) ;
topEssentialMAPKgene = corrDataDep.Properties.VariableNames{locMax} ;
% [~, locMin] = mink( geneEssentalityValue , 1) ;
% leastEssentialMAPKgene = corrDataDep.Properties.VariableNames{locMin} ;

% plot the data for MYC genes
xVar = corrDataDep.(topEssentialMAPKgene)  ;
yVar = CCLEmrna.(topEssentialMAPKgene) ;

figure()
scatter(xVar, yVar,40,'filled','MarkerFaceColor',[0.3010 0.7450 0.9330] )
hold on
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold')

% add a reference line to the plot
X = [ones(length(xVar),1) xVar] ;
b1 = X\yVar    ;     yCalc1 = X*b1 ;
plot(xVar,yCalc1,'LineWidth',2.5)
% calculate the pearson's linear correation
[r2 , pR ] = corr(xVar,yVar, 'Type','Pearson');
% legend({'points','Linear'},'Location','Best')
text(0.4, 0.9, strcat( sprintf("R = %0.2f, P ", r2),...
    convertPValue2SuperScript(pR)) ,'FontSize',14,'FontWeight','bold',...
    'Units','normalized')

% add figure labels
xlabel('Dependency Score')
ylabel('mRNA Expression')
title('MYC: mRNA-Depedency Score Correlation','FontSize',14)

%% do the same for all genes and dependece scores

% return only the oncogenes for the achilles data and the mrna expression
% data and measure the correlation 
corrDataDepOGs = corrDataDep(:, ismember( ...
    corrDataDep.Properties.VariableNames, oncogenes) );
CCLEmrnaOGs = CCLEmrna(:, ismember( ...
    CCLEmrna.Properties.VariableNames , oncogenes) ) ;
    
% throw in an assertion 
assert( all( strcmp( CCLEmrnaOGs.Properties.VariableNames' ,...
   corrDataDepOGs.Properties.VariableNames') ) )

% get the data correlation
corrDataDepOGs = corrDataDepOGs{:,:} ;
corrMrna = CCLEmrnaOGs{:,:} ;
rowMeanDep = corrDataDepOGs (:)  ;
rowMeanMrna = corrMrna(:) ;

% set the color

figure()
scatter(rowMeanDep, rowMeanMrna,25, 'filled' ,...
    'MarkerFaceColor',[0.4660 0.6740 0.1880] )
hold on
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold','YLim',[-0.5 14])

% add a reference line to the plot
X = [ones(length(rowMeanDep),1) rowMeanDep] ;
b1 = X\rowMeanMrna    ;     yCalc1 = X*b1 ;
plot(rowMeanDep,yCalc1,'LineWidth',2)
% legend({'points','Linear'},'Location','Best')
% calculate the pearson's linear correation
[r2 , pR ] = corr(rowMeanDep,rowMeanMrna, 'Type','Pearson', ...
    'Rows','complete') ;
if pR == 0
    text(0.4, 0.9, strcat( sprintf("R = %0.2f, P ", r2),...
        convertPValue2SuperScript(pR)) ,'FontSize',14,'FontWeight','bold',...
        'Units','normalized')
else
    text(0.4, 0.9, strcat( sprintf("R = %0.2f, P = ", r2),...
        convertPValue2SuperScript(pR)) ,'FontSize',14,'FontWeight','bold',...
        'Units','normalized')
end

% add figure labels
xlabel('Dependency Score')
ylabel('mRNA Expression (TPM)')
title('MAPK Oncogenes','FontSize',14)

clear xVar yVar x b1 corrDataDepReg rowMeanDep rowMeanMrna yCalc1 ...
    geneEssentalityValue locMax locMin locColumns locRows r2 pR ...
    corrMrna CCLEmethyl CCLEmutationNum CCLEcna corrDataDepOGs ...
    CCLEmrnaOGs

%% Plot a Histogram of Dependence Verse Expression 

% ***********************************************************************
% To further explore this class, we directly correlated gene expression and
% gene dependence.
% ***********************************************************************

% get a matched dataset for the crispr data and the mrna data and arrange
% the cell line in the same order

% load the expresssion data 
fprintf('\n Loading the mRNA expression data \n')
cclemrna1 = readtable('CCLE_mrna_processed.csv');
cclemrna2 = readtable('CCLE_mrna_nonMAPK.csv');
CCLEmrnaAll = innerjoin( cclemrna1, cclemrna2) ;

fprintf('\nFinding correlation between gene fitness and expression \n')
matchCrisprData = crisprAll( ... 
    ismember(crisprAll.cell_line, CCLEmrnaAll.cell_line) , : )  ;
[~, locThem] = ismember(matchCrisprData.cell_line, CCLEmrnaAll.cell_line) ;
CCLEmrnaAll = CCLEmrnaAll(locThem, :) ;

% also do the same for the genes 
[~, themA, themB ] = intersect( ...
    matchCrisprData.Properties.VariableNames ,...
    CCLEmrnaAll.Properties.VariableNames ,'stable') ;
CCLEmrnaAll = CCLEmrnaAll(:, themB) ;
matchCrisprData = matchCrisprData(:, themA);

% throw in an assertion 
assert(all(strcmp(matchCrisprData.cell_line, CCLEmrnaAll.cell_line)))
assert(all(strcmp(  matchCrisprData.Properties.VariableNames',...
    CCLEmrnaAll.Properties.VariableNames')) ) 

% preallocate
crisprmrnaCorrAll = matchCrisprData.Properties.VariableNames(2:end)' ;

% get the correlation
for ii = 1:length(crisprmrnaCorrAll)
    % get the current gene cripsr data and mrna expression data
    curGene = crisprmrnaCorrAll{ii} ;
    
    % check if th gene is present in the expression dataset if not then
    % continue
    if ~any( ismember(CCLEmrnaAll.Properties.VariableNames,curGene) )
        crisprmrnaCorrAll(ii,[2,3]) = num2cell([NaN,NaN]) ;
        continue
    else     
        % get the gene data 
        curDepGene = matchCrisprData.(curGene);
        curmRNA = CCLEmrnaAll.(curGene) ;
        
        % find the Pearson's correlation scores
        [r2 , pR ] = corr(curDepGene,curmRNA, 'Type','Pearson',...
            'Rows','complete') ;
        
        % put them in the matrix 
        crisprmrnaCorrAll(ii,[2,3]) = num2cell( [r2 , pR ]);
        
    end
end

% convert to a table and clean up the variables without p values and also
% add the essentiality labels
crisprmrnaCorrAll = array2table(crisprmrnaCorrAll , 'VariableNames' ,...
    {'HugoSymbol','Correlation','pValue'} ) ;

% load the fitness data from achilles and clean up the table 
achilleEssentialGenes = readtable('Achilles_common_essentials.csv');
achilleEssentialGenes = [achilleEssentialGenes;
    achilleEssentialGenes.Properties.VariableNames ] ;
achilleEssentialGenes.Properties.VariableNames = {'HugoSymbol','Entrez'} ;

% check that the genes in crisprmrnaCorrAll are arrange in the sample order
% as those in the matchCrisprData
assert(all(strcmp( matchCrisprData.Properties.VariableNames(2:end)' ,...
    crisprmrnaCorrAll.HugoSymbol) ) ) 

% do the same for the mrna data 
[~ ,locThem] = ismember(crisprmrnaCorrAll.HugoSymbol, ...
    CCLEmrnaAll.Properties.VariableNames) ;
locThem = locThem(locThem ~= 0);
CCLEmrnaAllMatched = CCLEmrnaAll(:, locThem) ;

% do the same for the crispr data 
[~ ,locThem] = ismember(CCLEmrnaAllMatched.Properties.VariableNames, ...
    crisprmrnaCorrAll.HugoSymbol) ;
locThem = locThem(locThem ~= 0);
crisprmrnaCorrAll = crisprmrnaCorrAll(locThem,:) ;

% throw in an assertion
assert(all(strcmp( CCLEmrnaAllMatched.Properties.VariableNames' ,...
    crisprmrnaCorrAll.HugoSymbol) ) ) 

% add the gene fitness rank of the genes let me also add the crispr gene
% score to the table
locThem = ismember(crisprmrnaCorrAll.HugoSymbol, ...
    achilleEssentialGenes.HugoSymbol) ;
crisprmrnaCorrAll = addvars(crisprmrnaCorrAll ,locThem , ...
    sum( matchCrisprData{:,2:end} < -0.5 )' ,...
    mean(CCLEmrnaAllMatched{:, :})' ,...
    'NewVariableNames',{'EssentialityClass','DependentCCL','meanExpr'})  ; 
crisprmrnaCorrAll.Correlation = cell2mat(crisprmrnaCorrAll.Correlation);
crisprmrnaCorrAll.pValue = cell2mat(crisprmrnaCorrAll.pValue);
crisprmrnaCorrAll( isnan( crisprmrnaCorrAll.pValue ), :) = [] ;
crisprmrnaCorrAll = sortrows(crisprmrnaCorrAll, 'pValue','ascend');

% =========================== plot the histogram ========================

% % the first one for the essential genes
% figure()
% histogram(crisprmrnaCorrAll.Correlation( ...
%     crisprmrnaCorrAll.EssentialityClass == true) )
% hold on 
% histogram(crisprmrnaCorrAll.Correlation( ...
%     crisprmrnaCorrAll.EssentialityClass == false) )
% ylabel('Count')
% xlabel('Self mRNA vs Dependency Correlation')
% set(gca,'LineWidth',1.5,'Box','off','FontWeight','bold')
% legend({'Essential','Others'},'Location','Best') 
% hold off
% 
% % the second one for the high ranking genes genes
% figure()
% histogram(crisprmrnaCorrAll.Correlation( ...
%     crisprmrnaCorrAll.DependentCCL >= 50 ) )
% hold on 
% histogram(crisprmrnaCorrAll.Correlation( ...
%     crisprmrnaCorrAll.DependentCCL < 50) )
% ylabel('Count')
% xlabel('Self mRNA vs Dependency Correlation')
% set(gca,'LineWidth',1.5,'Box','off','FontWeight','bold')
% legend({'High Ranking','Others'},'Location','Best') 
% hold off

%% =====  the first one for the essential genes ==========
figure()
subplot(1,4,1)
% the essential genes
histogram(crisprmrnaCorrAll.Correlation( ...
    crisprmrnaCorrAll.Correlation <= -0.3), 'FaceColor',[0.1 0.1 1])
hold on 
% the inbetween genes
histogram(crisprmrnaCorrAll.Correlation( ...
    crisprmrnaCorrAll.Correlation > -0.3 &  ...
    crisprmrnaCorrAll.Correlation < 0.3) ,300, ...
    'FaceColor',[0.7 0.7 0.7] ,'EdgeColor',[0.7 0.7 0.7])
% the tumour suppressor genes
histogram(crisprmrnaCorrAll.Correlation( ...
    crisprmrnaCorrAll.Correlation >= 0.3), 'FaceColor',[1 0.1 0.1])

ylabel('Count')
xlabel('Self mRNA vs Dependency Correlation')
set(gca,'LineWidth',1.5,'Box','off','FontSize',12,'FontWeight','bold')

% add a legend to the plot
legendLabels = {'Over Expressed Oncogenes','Others','Essential/TSGs'} ;
plotColors = [0.1 0.1 1; 0.7 0.7 0.7 ; 1 0.1 0.1 ] ;
createLegendInternal(0.85, 0.15, legendLabels , plotColors,...
    'Groups' , [12 10] ,[0.018 ,0.12])
hold off

% == plot the scatter plot for the most negative and positive correlation =  

% do the same for the second gene
myGenes = {'IRF4','SOX10','LTO1'}; % 'LTO1' 'CDKN1A'

for ii = 1:length(myGenes)
    xVar = matchCrisprData.(myGenes{ii})  ;
    yVar = CCLEmrnaAll.(myGenes{ii})   ;
    
    % change the facecolor 
    if ii ~= 3
        myColor = [ 0.1 0.1 1 ];
    else
        myColor = [1 0.1 0.1 ] ;
    end
    
    subplot(1,length(myGenes)+1,ii+1)
    scatter(xVar, yVar,30,'filled','MarkerFaceColor', myColor )
    hold on
    set(gca,'FontSize',12,'LineWidth',1.5,'FontSize',12 ,'Box','off',...
        'FontWeight','bold','YLim', [ -1 10])
    
    % add a reference line and the correlation coefficent
    addReferenceLineToPlot(xVar, yVar)
    
    % add figure labels
    xlabel('Dependency Score')
    ylabel('mRNA Expression')
    title(sprintf('%s',myGenes{ii}) ,'FontSize',16)
end


% save to excel 
writetable( crisprmrnaCorrAll, 'Supplementary File 3.xlsx', ...
    'Sheet','Fitness-Expression Correlation')

clear r2 pR curDepGene curmRNA curGene matchCrisprData CCLEmrnaAll ...
    cclemrna1 cclemrna2 myColor xVar yVar matchCrisprData ...
    plotColors crisprmrnaCorrAll myGenes legendLabels plotColors ...
    locThem achilleEssentialGenes r2 pR ii curGene curmRNA ...
    cclemrna1 cclemrna2 themA themB

% What thresholds should I use to decide if a gene is really having a
% significant effect on a cell line? Although it depends on the risk of
% false positives you’re willing to tolerate, for most applications a
% cutoff of 0.5 in gene dependency probability or greater makes sense. For
% gene effect, a score less than -0.5 represents noteworthy depletion in
% most cell lines.

%% Plot the Gene Ontology Biological Processess

% I used Enricher to query the GO Biological Processess that are
% overrepresented using two list of genes: 1: those which show a positive
% correlation coeffienct that is greater than 0.3 between the CRISPR
% fitness and the self-mRNA expression and those which have a correlation
% coefficent less than -0.3 i.e., negative

% load the enrichr results'
negCrisprGenes = readtable(...
            'GO_Biological_Process_2018_CRISPR_negative_correlation.txt') ;
negCrisprGenes = sortrows(negCrisprGenes,'AdjustedP_value','ascend') ;
posCrisprGenes = readtable(...
            'GO_Biological_Process_2018_CRISPR_positive Correlation') ;
posCrisprGenes = sortrows(posCrisprGenes,'AdjustedP_value','ascend') ;

% specificy the number of rows to plot
numOfrows = 15 ;
% plot the bar graphs
figure()
% plot for the tumour suppressor genes
barh( -log(negCrisprGenes.P_value(1:numOfrows)) ,'r',...
    'FaceAlpha',0.5 , 'EdgeAlpha',0.5 ,'EdgeColor','r');
set(gca,'YDir','reverse','Box','off','YTick',[],'LineWidth',1 ,...
    'YLim',[0.5 numOfrows+0.5] ,'XLim',[14 28],'FontWeight', 'bold', ...
    'XColor',[1 1 1] ,'XTick',[])
ylabel('Ontology Terms')
title('Essential /TSG Enrichment','FontSize',14 ,'FontWeight', 'bold')
for ii = 1:numOfrows
   text(14.1,ii, negCrisprGenes.Term(ii),'FontSize' ,12,'FontWeight','bold')
end

% plot for the oncogenes and overexpressed transcription factors
figure()
barh( -log(posCrisprGenes.P_value(1:numOfrows)) ,'b',...
    'FaceAlpha',0.3 ,'EdgeAlpha',0.3, 'EdgeColor','b');
set(gca,'YDir','reverse','Box','off','YTick',[],'LineWidth',1 ,...
    'YLim',[0.5 numOfrows+0.5] ,...
    'XLim', [0 -log(posCrisprGenes.P_value(1))+3 ], ...
    'FontWeight', 'bold', 'XColor',[1 1 1] ,'XTick',[])
ylabel('Ontology Terms')
title('Over Expressed Oncogenes / Transcription Factors', ...
    'FontSize',14 ,'FontWeight', 'bold')
for ii = 1:numOfrows
   text(0.1,ii, posCrisprGenes.Term(ii),'FontSize' ,12,'FontWeight','bold')
end

% save the result tot he suppllementary file 
writetable(negCrisprGenes, 'Supplementary File 3.xlsx','Sheet',...
    'GO Bio Process CRISPR NEG corr') ;
writetable(posCrisprGenes,'Supplementary File 3.xlsx','Sheet', ...
    'GO Bio Process CRISPR POS corr') ;

%% Are gene essentialities associated with and mutations 

% Show that pancreatic cancer cell line are highly dependent on KRAS which
% is mutated in most tumours
krasCrispr = clusterDataDep(:,{'disease','cell_line','KRAS'} )  ;
krasExpr = innerjoin( sampleInfo(:, {'disease','cell_line'}) , ...
    CCLEmrna( :,{'cell_line','KRAS'} ) ) ;

% also add the mutation data to the expression data 
krasMut = CCLEmutations(:,{'cell_line','KRAS'} )  ;
krasMut.Properties.VariableNames(2) = "KRASmut" ;

% add that to the kras mRNA expression data 
krasExpr = innerjoin(krasExpr, krasMut) ;

% make sure that the two dataset have the same cell line and arrange them
% in the same order between the mRNA expression data and the CRISPR data
[~, iA, iB] = intersect(krasCrispr.cell_line, krasExpr.cell_line) ;
krasCrispr = krasCrispr(iA, :) ;
krasExpr = krasExpr(iB, :) ;

% throw in an assertion 
assert( all( strcmp( krasCrispr.cell_line, krasExpr.cell_line) ) )

% set the plot colors between pancreatic cancer and other cancers and also
% between those with mutations in the KRAS genes and those with no
% mutations in the kras genes

% set the colors 
colors1 = [0.392 0.831 0.074 ; 0.1 0 1] ;

% ================= produce the scatter plots ==================
figure() 
clf ;
axes('position',[0.15, 0.33, 0.74, 0.60]);
gscatter( krasCrispr.KRAS , krasExpr.KRAS, ...
    ismember(krasCrispr.disease,'Pancreatic Cancer') , colors1, '..',30)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('KRAS Dependency Score') ; 
ylabel('KRAS mRNA Expression') ; 
title('Cell Fitness: Pancreatic vs Other CLL','FontSize',14',...
    'FontWeight','bold')
hold on 
% add a reference line and the correlation coefficent 
addReferenceLineToPlot(krasCrispr.KRAS, krasExpr.KRAS)
legend({'Others CCL','Pancreatic CLL'} , 'Location','NorthEast')

% ============ add box plot at the bottom of the chart =============

% get the limit of the scatter plot to use those as the limit of the bar
% graph 
theXlim = get(gca,'XLim') ;

% plot the data
axes('position',[0.15, 0.1, 0.74, 0.15]);
boxplot( ...
    krasCrispr.KRAS , ismember(krasCrispr.disease,'Pancreatic Cancer'),...
    'Orientation','horizontal','Width',0.7 ,'Color', colors1 ,...
    'Symbol','r+' ) ;

% set some figure properties and add title ot the figure
set(gca,'Box','off','YColor',[1 1 1], 'XColor',[1 1 1] ,'XLim', theXlim)

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)

% add color to the box plots
colors2 = flipud(colors1) ;
h4 = findobj(gca,'Tag','Box') ;
for jj=1:length(h4)
    patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
        colors2(jj,:),'FaceAlpha', .8 ,'LineWidth',1);
end

% perform a ttest and add the p value to the plot 
[~, pValue ] =  ttest2(krasCrispr.KRAS , ...
    ismember(krasCrispr.disease,'Pancreatic Cancer'),...
    'Vartype','unequal') ;
text( 0.4, 0.1 , convertPValue2SuperScript(pValue) ,...
    'Units','normalized', 'FontWeight','bold','FontSize',12)
hold off

% ======================== subplot 2 =============================
figure()
clf ;
axes('position',[0.15, 0.33, 0.74, 0.60]);
colors1 = [1 0.43 0.1649 ; 0.0745 0.623 1] ;
gscatter( krasCrispr.KRAS , krasExpr.KRAS, ...
    ~cellfun(@isempty ,krasExpr.KRASmut) , colors1, '..',30) ;
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('KRAS Dependency Score') ; 
ylabel('KRAS mRNA Expression') ; 
title('Cell Fitness: KRAS Mutation Status','FontSize',14',...
    'FontWeight','bold')
hold on 
% add a reference line and the correlation coefficent 
addReferenceLineToPlot(krasCrispr.KRAS, krasExpr.KRAS)
legend({'not mutated','KRAS mutated'} , 'Location','NorthEast')
hold off

% ============ add box plot at the bottom of the chart =============

% get the limit of the scatter plot to use those as the limit of the bar
% graph 
theXlim = get(gca,'XLim') ;

% plot the data
axes('position',[0.15, 0.1, 0.74, 0.15]);
boxplot( ...
    krasCrispr.KRAS , ~cellfun(@isempty ,krasExpr.KRASmut) ,...
    'Orientation','horizontal','Width',0.7 ,'Color', colors1 ,...
    'Symbol','r+') ;

% set some figure properties and add title ot the figure
set(gca,'Box','off','YColor',[1 1 1], 'XColor',[1 1 1] ,'XLim',theXlim)

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)

% add color to the box plots
colors2 = flipud(colors1) ;
h4 = findobj(gca,'Tag','Box') ;
for jj=1:length(h4)
    patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
        colors2(jj,:),'FaceAlpha', .8 ,'LineWidth',1);
end

% perform a ttest and add the p value to the plot 
[~, pValue ] =  ttest2(krasCrispr.KRAS , ...
    ~cellfun(@isempty ,krasExpr.KRASmut) ,...
    'Vartype','unequal') ;
text( 0.4, 0.1 ,sprintf('p = %0.3e',pValue) ,'Units','normalized',...
    'FontWeight','bold','FontSize',12)
hold off

% plot some pathways to show the connection between KRAS, MYC, BCAR1 - the
% three genes to which pancreatic cancer cell lines are most dependent on
% from the MAPK kinase pathway

clear iA iB krasExpr krasCrispr colors1 krasMut ovCancer colors2 ...
    pValue theXlim CCLEmutationMun

%% Check Gene Essentiality for Other genes and cancer

% show that those cancer cell line with 100% p53 mutations increase their
% fitness once p53 is silenced

% define the input parameters for SKIN cancer BRAF mutations
checkGene = 'BRAF' ;
cancerCheck = 'Skin Cancer' ;

% the pvalues for the left plot then the t stats for the left plot and the
% t stats for the right plot
[pValueBRAF, tStatsBRAF, stats2BRAF ]= getMutationExprCrisprAssociation(...
    checkGene, clusterDataDep, ...
    sampleInfo , CCLEmutations , CCLEmrna , cancerCheck) ;

% define the input parameters for SKIN cancer BRAF mutations
checkGene = 'MAPK1' ;
cancerCheck = 'Skin Cancer' ;
getMutationExprCrisprAssociation( checkGene, clusterDataDep, ...
    sampleInfo , CCLEmutations , CCLEmrna , cancerCheck)

% define the input parameters for pancreatic cancer KRAS mutations
checkGene = 'KRAS' ;
cancerCheck = 'Pancreatic Cancer' ;
getMutationExprCrisprAssociation( checkGene, clusterDataDep, ...
    sampleInfo , CCLEmutations , CCLEmrna , cancerCheck)

% lung cancer
checkGene = 'NRAS' ;
cancerCheck = 'Skin Cancer' ;
getMutationExprCrisprAssociation( checkGene, clusterDataDep, ...
    sampleInfo , CCLEmutations , CCLEmrna , cancerCheck)

% Mutated oncogenes (e.g., NRAS, BRAF, and KRAS) are among the most robust
% dependencies observed in DRIVE (Figure 2A). As expected, BRAF mutation
% and dependence are mainly observed in colon, thyroid, and melanoma lines,
% while NRAS mutation and dependence is most prevalent
% in melanoma. KRAS mutation and dependence occurs in colon,
% pancreatic, and lung lineages. Consistent with TCGA (Cancer Genome
% Atlas Network, 2012) and emerging clinical data (Mayer et al., 2017),

%%
% checkGene = 'BRAF' ;
% cancerCheck = 'Breast Cancer' ;
% getMutationExprCrisprAssociation(checkGene, clusterDataDep, ...
%     sampleInfo , CCLEmutations , CCLEmrna , cancerCheck)

%% Relationship Between Gene Essentility and Copy Number Change

% let get the MAPK pathway genes most affect by CNV 
% also get percentage of alteration involved in MAPK
cnaAlterations = cnaData ;
cnaAlterations( cellfun(@isempty, cnaAlterations.CancerStudy) , :)  = [] ;
cnaAlterations.CancerStudy = categorical(cnaAlterations.CancerStudy) ;
cnaAlterations.SampleIds = [] ;

% get all values that are not equal to 0
cnaAlterations{:,2:end} = cnaAlterations{:,2:end} ~= 0 ;

% calculate the frequency of each gene CNA 
cnaAlterations = grpstats(cnaAlterations, 'CancerStudy', ...
    'sum', 'DataVars', cnaAlterations.Properties.VariableNames(2:end) ) ;
cnaAlterations.Properties.VariableNames(3:end) = ...
    extractAfter(cnaAlterations.Properties.VariableNames(3:end), ...
    '_') ;

% now get the percetange
cnaAlterations{:,3:end} = ( cnaAlterations{:,3:end}./ ...
    cnaAlterations.GroupCount  ) * 100;

% get the genes most affected by copy number changes 
cnaMostAlteredGenes = array2table( max(cnaAlterations{:,3:end},[],1), ...
    'VariableNames',cnaAlterations.Properties.VariableNames(3:end) )  ;
cnaMostAlteredGenes = rows2vars(cnaMostAlteredGenes) ;

% add variable names and sort the rows 
cnaMostAlteredGenes.Properties.VariableNames = ...
    {'HugoSymbol','MaxAlteration'};

% also add th percentange alteration to the table 
cnaMostAlteredGenes = addvars( cnaMostAlteredGenes , ...
    mean(cnaAlterations{:,3:end})' , 'After','HugoSymbol' ,...
    'NewVariableNames','Percentage') ;

cnaMostAlteredGenes = sortrows(cnaMostAlteredGenes,...
    'Percentage','descend');

writetable(cnaMostAlteredGenes ,'Supplementary File 3.xlsx', ...
    'Sheet','CNA Alterations');

%% ======== get the ccle copy number data and clean up the data =======
fprintf('\n load CCLE copy number data \n') 

CCLEcna = readtable('CCLE_cnvOldData.txt');
CCLEcna(:, {'Var2','Var3'} ) = [] ;
CCLEcna.Properties.VariableNames = matlab.lang.makeUniqueStrings( ...
    ['HugoSymbol', CCLEcna{1, 2:end} ] ) ;
CCLEcna( 1:3 , :) = [] ;

% return only the MAPK genes and transpose the table
CCLEcna = CCLEcna( ismember( CCLEcna.HugoSymbol , mapkGenes.Gene), :) ;
CCLEcna = rows2vars(CCLEcna ,'VariableNamesSource','HugoSymbol') ;
CCLEcna.Properties.VariableNames(1) = "cell_line" ;

% convert the data to double from cell 
for ii = 2:width(CCLEcna)
    CCLEcna.(ii) = str2double(CCLEcna.(ii)) ;  
end

%% check the correlation and dependency 

checkGene = 'FLNA' ;
cancerCheck = 'Ovarian Cancer' ;
getMutationExprCrisprAssociation(checkGene, clusterDataDep, ...
    sampleInfo , CCLEcna , CCLEmrna , cancerCheck) ;

checkGene = 'DUSP9' ;
cancerCheck = 'Breast Cancer' ;
getMutationExprCrisprAssociation(checkGene, clusterDataDep, ...
    sampleInfo , CCLEcna , CCLEmrna , cancerCheck) ;

%% Can we predict the genes essentiality based on the pan cancer expression

% % tranpose the expression and dependence score dataset
% machineLearningDep = rows2vars(matchedDepData,...
%     'VariableNamesSource','cell_line') ;
% machineLearningMrna = rows2vars(CCLEmrna, ...
%     'VariableNamesSource','cell_line') ;
%
% % change the first variable names to genes
% machineLearningDep.Properties.VariableNames(1) = "HugoSymbol" ;
% machineLearningMrna.Properties.VariableNames(1) = "HugoSymbol" ;
%
% % replace the missing values
% machineLearningDep = fillmissing( machineLearningDep ,'nearest', ...
%     'DataVariables', machineLearningDep.Properties.VariableNames(2:end) );
% machineLearningMrna= fillmissing( machineLearningMrna ,'nearest', ...
%     'DataVariables', machineLearningMrna.Properties.VariableNames(2:end) );
%
% % throw in an assertion
% assert(all( strcmp( machineLearningDep.HugoSymbol , ...
%     machineLearningMrna.HugoSymbol)))
%
% % get machine learning data for one cell lines
% machinePredData = addvars( machineLearningMrna, machineLearningDep.(6) ,...
%     'After','HugoSymbol','NewVariableNames',...
%     strcat('pred_',machineLearningDep.Properties.VariableNames(6)  ) );

clear machineLearningDep machineLearningMrna

%% *********************************************************************
% Are the most essential highly expressed than those genes that non
% essential across the cancer cell lines
% *********************************************************************


% % now get the data of the cancer study and change the arrangement in the
% cancer studies
barData = clusterDataDep(:,3:end) ;
[~, locCols] = ismember(cgo.ColumnLabels, barData.Properties.VariableNames );
barData = barData(:,locCols) ;

% throw in an assertion
assert(all(strcmp(cgo.ColumnLabels',barData.Properties.VariableNames')))

% get the values above a give threshold of the essentialies scores
barData = table2array(barData) <  -0.50;
barData = sum(barData, 1);

% get data of the same length and only for the genes on which the cell
% lines are highly dependence on
cutOffDep = 100 ;

% compare the expression of highly essential genes to those that are not
% essential
mrnaHeatData = CCLEmrna{:,2:end} ;
essentialExpr = nanmean( mrnaHeatData(:, barData > cutOffDep) ,2) ;
othersExpr = nanmean( mrnaHeatData(:, barData <= cutOffDep) ,2) ;

% create any array for the groups and compare the mean mrna expression
% between the essential genes and the non essential genes
[~,pValue,~,~] = ttest2(essentialExpr,othersExpr,'Vartype','unequal') ;

% create an array of groups
essentialGroups = ...
    [repmat({'Core Essential Genes'},length(essentialExpr),1) ;...
    repmat({'Other Genes'}, length(othersExpr),1) ] ;

% plot the data
figure()
boxplot([essentialExpr;othersExpr] , categorical(essentialGroups) ) ;
hold on

% set some figure properties and add title ot the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold')
ylabel('mRNA Transcript Levels (log TPM)')
% title('Gene Transcription','FontSize',16)

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)

% add color to the box plots
color = [0.9900 0.2250 0.0980;0.1940 0.3840 0.9560 ] ;
h4 = findobj(gca,'Tag','Box') ;
for jj=1:length(h4)
    patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
        color(jj,:),'FaceAlpha', .8 ,'LineWidth',1);
end
text(0.45, 0.9,sprintf('p = %0.1f', pValue) , ...
    'FontWeight','bold','FontSize',14 ,'Unit','normalize')

hold off

clear jj color essentialExpr pValue CI tStats cutOffDep ...
    essentialGroups othersExpr essentialExpr mrnaHeatData LocThem ...
    theMissing these missingGenes locRows locColumns corrDataDep ...
    curGene ii leastEssentialMAPKgene geneEssentalityValue

%% Consistentely Expressed Highly Expressed Gene

% *********************************************************************
% Are other genes that are consistently up regulated across cancer cell
% line more likely to be essential than those not up regulated
% *********************************************************************

% load all the mrna expression data and change the names of the cell lines
% to match those in the crispr data
fprintf('\n Loading all the mrna expression data \n')
CCLEothermrna = readtable('CCLE_mrna_nonMAPK.csv') ;

% get the intersect between the two datasets
[~, locA, locB ] = intersect(...
    CCLEothermrna.Properties.VariableNames,...
    crisprAll.Properties.VariableNames ,'stable') ;

% get the clustering data
CCLEothermrnaClust = CCLEothermrna(:, locA) ;
crisprAllClust = crisprAll(:, locB) ;

% throw in an assertion
assert(all(strcmp(CCLEothermrnaClust.Properties.VariableNames', ...
    crisprAllClust.Properties.VariableNames') )) ;

% filling the missing values
fprintf('\n Filling the missing values \n')
crisprAllClust = fillmissing(crisprAllClust,'linear','DataVariables',...
    crisprAllClust.Properties.VariableNames(2:end)) ;

% cluster the data
if ~exist('cgoAllAchilles.mat','file')
    fprintf('\n Clustering the large dataset \n')
    cgoAll = clustergram( crisprAllClust{:,2:end}, ...
        'RowLabels',crisprAllClust.cell_line, ...
        'ColumnLabels', crisprAllClust.Properties.VariableNames(2:end),...
        'colormap',redbluecmap ,...
        'standardize','row','ColumnPDist','cosine','Linkage','average');
    
    % save the clustergram object
    fprintf('\n Saving the clustergram object \n')
    save('cgoAllAchilles.mat','cgoAll') ;
    fprintf('\n Clustergram has been saved \n')
else % just the save clustergram object
    fprintf('loading the clustergram data \n')
    load 'cgoAllAchilles.mat'
end

%% Show Relationship Between Fitness and Expression for all Tissues

% there is a small bug where the row label seem to change by themselves
fprintf('\n Producing the heatmap - will take a while \n')
[~, locCols] = ismember(cgoAll.ColumnLabels, ...
    crisprAllClust.Properties.VariableNames(2:end) ) ;
[~, locRows] = ismember(cgoAll.RowLabels, crisprAllClust.cell_line) ;
heatDataAll = crisprAllClust{:,2:end} ;
heatDataAll =  heatDataAll(locRows,locCols);

% make the mrna and methylation data the same length and order as the
% crispr data
[~, LocThem ] = ismember(cgoAll.ColumnLabels  , ...
    CCLEothermrna.Properties.VariableNames) ;
LocThem(LocThem == 0) = [] ;
CCLEothermrnaclust = [CCLEothermrna(:,1), CCLEothermrna(:,LocThem)];

% also get only the cell line that are present in the crispr
[~, LocThem ] = ismember(flip(cgoAll.RowLabels), ...
    CCLEothermrnaclust.cell_line);
LocThem(LocThem == 0) = [] ;
CCLEothermrnaclust = CCLEothermrnaclust(LocThem,:) ;

% throw in an assertion
assert( all( strcmp( ...
    CCLEothermrnaclust.Properties.VariableNames(2:end)',...
    cgoAll.ColumnLabels') ) )

% also get the evenly process tcga tumour data and the gtex data
% this data contains both normal tissue, cancer and adjucent tissues
Gtex = readtable('Normal_GTeX.txt','HeaderLines',4);
Gtex.GeneID = [] ;

% get only the genes present in other tissues
Gtex = Gtex(ismember(Gtex.GeneName, cgoAll.ColumnLabels),  :) ;

% transpose the table and add the missing genes to the table
Gtex = rows2vars(Gtex,'VariableNamesSource','GeneName') ;
Gtex.Properties.VariableNames(1) = "Tissue" ;

% make the Gtex data the sample length as the clustergram all data
% get the missing genes and add them to the Gtex data using NaN
these = ismember(cgoAll.ColumnLabels,Gtex.Properties.VariableNames(2:end));
missingGenes = cgoAll.ColumnLabels(~these);

% add the NaN
theMissing = array2table( NaN(height(Gtex), length(missingGenes) ) ,...
    'VariableNames', missingGenes) ;
Gtex = [ Gtex , theMissing ] ;

% now get the genes that are present in the mrna data and gtex data
Gtex = [Gtex(:,1),  Gtex(:, ismember(...
    Gtex.Properties.VariableNames , cgoAll.ColumnLabels) ) ] ;

% throw in an assertion
assert( width(Gtex(:,2:end)) == length(cgoAll.ColumnLabels) )

% arrange the genes in the same order and the column labels
[~, locThem ] = ismember(...
    cgoAll.ColumnLabels ,Gtex.Properties.VariableNames) ;
locThem(locThem == 0)  = [] ;
Gtex =[Gtex(:,1) , Gtex(:, locThem) ];

% throw in an assertion
assert(all( strcmp( Gtex.Properties.VariableNames(2:end)' ,...
    cgoAll.ColumnLabels') ))

% load the TCGA mrna data for tumours using the TCGA ids that are present
% in the cancer studies and set the reduce size to false because I want to
% return all the genes all the cancer studies
if ~exist('mrnaTCGAHeatmapData.mat','file')
    % get all the mrna expression data for the TCGA tumours
    reduceSize = false ;
    tcgaCancerTypes = extractBefore( cellstr(cancerStudies.cancerTypeId(...
        contains(cancerStudies.name,'TCGA') ) ) ,'_')  ;
    mrnaTCGA = getMedianTranscritomicsData(tcgaCancerTypes, reduceSize) ;
    
    % arrange the TCGA mrna expression data in the same order as that of
    % the essentiallity scores from the achilles project arrange the genes
    % in the same order and the column labels
    [~, locThem ] = ismember(...
        cgoAll.ColumnLabels , mrnaTCGA.Properties.VariableNames) ;
    locThem(locThem == 0)  = [] ;
    mrnaHeatData =[mrnaTCGA(:,1) , mrnaTCGA(:, locThem) ];
    
    %     % get the median expression for each TCGA tumour
    %     % ***I can produce this plot without filtering on powerful computer***
    %     mrnaHeatData.CancerStudy = categorical(mrnaHeatData.CancerStudy) ;
    %     fprintf('\n Finding the mean expression from the tall array \n')
    %     fprintf('\n This will take a while if its the first run \n')
    %     mrnaHeatData = grpstats(mrnaHeatData,'CancerStudy','mean') ;
    %     fprintf('\n Gathering tall array data \n')
    %     mrnaHeatData = gather(mrnaHeatData) ;
    
    % save the data so that next time I dont have to wait for a long time
    % for this analysis to run
    save('mrnaTCGAHeatmapData.mat','mrnaHeatData')
    
    % and remove the mrnaTCGA data to free up some space for processing the
    % dataset
    clear mrnaTCGA
else % if the data has already been processed
    fprintf('\n Loading the TCGA mrna expression data \n')
    load('mrnaTCGAHeatmapData.mat')
end

% make the TCGA data the sample length as the clustergram all data
% get the missing genes and add them to the TCGA data using NaN
these = ismember(cgoAll.ColumnLabels, ...
    mrnaHeatData.Properties.VariableNames(2:end));
missingGenes = cgoAll.ColumnLabels(~these);

% add the NaN
theMissing = array2table( NaN(height(mrnaHeatData), length(missingGenes) ) ,...
    'VariableNames', missingGenes) ;
mrnaHeatData = [ mrnaHeatData , theMissing ] ;

% and arrange the genes in the same order as the clutergram
[~, these] = ismember(cgoAll.ColumnLabels, ...
    mrnaHeatData.Properties.VariableNames);
mrnaHeatData = [ mrnaHeatData(:,1) , mrnaHeatData(:, these) ] ;

% throw in an assertion
assert(all( strcmp( mrnaHeatData.Properties.VariableNames(2:end)' ,...
    cgoAll.ColumnLabels') ))

% ==================== Produce Multiple Plots ====================
figure();
% clf ; set(gcf,'position',[200,100,1200,600]);
% the first number is how far the figure will be from the x-axis and the
% seceond number is now far the figure will be from the y-axis. The third
% number is was far the figure will run across the figure bar and the last
% number is far it will displaced allow the y-axis

axes('position',[0.15, 0.10, 0.74, 0.18]);
heatmap(heatDataAll,'Colormap',redbluecmap,'ColorbarVisible','off',...
    'FontColor','none','GridVisible','off' ,'ColorLimits',[-1 0.8]);

mrnaCCLEHeatDataAll = filloutliers(CCLEothermrnaclust{:,2:end} ,'linear');
axes('position',[0.15, 0.30, 0.74, 0.18]);
heatmap(mrnaCCLEHeatDataAll ,'Colormap',redgreencmap,'GridVisible','off',...
    'FontColor','none','ColorbarVisible','off','ColorLimits',[-10 18]);

axes('position',[0.15, 0.50, 0.74, 0.18]);
GtexHeatData = filloutliers(Gtex{:,2:end} ,'linear');
heatmap(GtexHeatData,'Colormap',redbluecmap,'GridVisible','off',...
    'FontColor','none','ColorbarVisible','off',...
    'MissingDataColor',[0.8 0.8 0.8] ,'ColorLimits',[-10 30] , ...
    'Title',...
    'Relation Between Transcript Expression and Gene Essentiality');

axes('position',[0.15, 0.70, 0.74, 0.18]);
heatmap(mrnaHeatData{:,2:end},'Colormap',parula,'GridVisible','off',...
    'FontColor','none','ColorbarVisible','off',...
    'MissingDataColor',[0.8 0.8 0.8] ,'ColorLimits',[-200 800] , ...
    'Title',...
    'Relationship Between mRNA Transcript Expression and Gene Essentiality');

% annoate the plot
% Create doublearrow
annotation('doublearrow',[0.152941176470588 0.889592760180996],...
    [0.072 0.072],'LineWidth',3);

% Create textbox
annotation('textbox',...
    [0.416384615384616 0.062 0.200809954751131 0.0299572039942939],...
    'String','          18023 Genes',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[1 1 1]);

% create a textbox
% Create textbox
textToPut = { ...
    sprintf('Achilles CRISPR Fitness: %d Cell lines',size(heatDataAll,1)), ...
    sprintf('mRNA: %d CCLE Cell Lines',size(mrnaCCLEHeatDataAll,1)),...
    sprintf('mRNA: %d GTEx Normal Tissues', height(Gtex) ) , ...
    sprintf('mRNA: %d TCGA Cancer Studies', size(mrnaHeatData,1)) } ;

textPos = [0.17, 0.37, 0.57, 0.77] ;
for ii = 1:length(textToPut)
    annotation('textbox',[0.0500 textPos(ii) 0.1000 0.0395],...
        'VerticalAlignment','middle','String', textToPut{ii} , ...
        'FontWeight','bold','FontSize',12,'FontName','Helvetica Neue',...
        'EdgeColor','none','VerticalAlignment','cap',...
        'HorizontalAlignment','right');
end

% Create textbox
annotation('textbox',...
    [0.285678733 0.89205734665 0.50879638009 0.04564907097],...
    'String', ...
    {'Relationship Between mRNA Transcript Expression and Gene Essentiality'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'EdgeColor',[1 1 1]);

fprintf('\n Rendering the heatmap \n')

% remove some variables
clear locX locY locThem CCLEothermrnaclust mrnaCCLEHeatDataAll ...
    GtexClust GtexHeatData mrnaColor textToPut textPos

%% Find the relationship between protein expression and gene essentiality

% *********************************************************************
% Are the essential genes more high expresssed at protein level
% *********************************************************************
% HPAproteins = readtable('normal_tissue_protein.csv') ;

% not enough proteomics data

%% Which MAPK pathways is more essential across cancer types

% get the essentiality for each pathway
pathways = categories(mapkGenes.Pathway);
pathwayCrispr = crispr(:,1) ;
for ii = 1:length(pathways)
    % get the crispr values for each pathway genes
    curPathwayCrispr = crispr{:, ismember( ...
        crispr.Properties.VariableNames ,...
        mapkGenes.Gene(mapkGenes.Pathway == pathways(ii) ) ) & ...
        ismember(crispr.Properties.VariableNames , oncogenes) } ;
    
    % add those to the table
    pathwayCrispr = addvars(pathwayCrispr, nanmean(curPathwayCrispr,2) ,...
        'NewVariableNames', pathways(ii) );
end

% get the essentiality for each mapk protein type
proteinClass = categories(mapkGenes.proteinClass);
proteinClassCrispr = crispr(:,1) ;
for ii = 1:length(proteinClass)
    % get the crispr values for each pathway genes
    curPathwayCrispr = crispr{:, ismember( ...
        crispr.Properties.VariableNames ,...
        mapkGenes.Gene(mapkGenes.proteinClass == proteinClass(ii) ) )...
        & ismember(crispr.Properties.VariableNames , oncogenes) } ;
    
    % add those to the table
    proteinClassCrispr = addvars(proteinClassCrispr, ...
        nanmean(curPathwayCrispr,2),'NewVariableNames', proteinClass(ii));
end

% ***********************************************************************
% which among the mapk gene is most essential to cancer cell lines
% ***********************************************************************

mostEssential = crispr{:,2:end};
for ii = 1:size(mostEssential,2)
    % get the values that are two standard deviations away
    % here I need to get 61 strongly selective (essential) MAPK instances
    mostEssential(:,ii) = mostEssential(:,ii) <  -0.50;
end
% save a copy of most essential for the next stage
mostEssentialOG = mostEssential;

% use this to find the most essential genes
mostEssential = [crispr(:,1) , array2table( mostEssential, ...
    'VariableNames',crispr.Properties.VariableNames(2:end) ) ] ;

% get the most essential gene
mostEssentialGene = table(crispr.Properties.VariableNames(2:end)', ...
    sum( mostEssential{:,2:end})','VariableNames',{'HugoSymbol','Score'});
mostEssentialGene = sortrows(mostEssentialGene,'Score','descend');
mostEssentialGene.HugoSymbol = categorical(mostEssentialGene.HugoSymbol);

% save the file to excel 
writetable(mostEssentialGene,'Supplementary File 3.xlsx','Sheet',...
    'Dependency of MAPK Genes')

% MAPK = Strongly Selective
% MYC = Common Essential
% data = readtable('Achilles_common_essentials.csv');
% head(data)
% data = data(ismember(data.AAAS ,mapkGenes.Gene) ,:)

% remove the gene with 0 scores
mostEssentialGene(mostEssentialGene.Score < 50,:) = [] ;
mostEssentialGene.HugoSymbol = categorical( ...
    cellstr( mostEssentialGene.HugoSymbol) , ...
    flipud( cellstr(mostEssentialGene.HugoSymbol)) );

% plot the data
figure()
barh(mostEssentialGene.HugoSymbol, mostEssentialGene.Score )
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold')
xlabel('Number of Dependent Cell Lines')
ylabel('MAPK Pathways Gene')
title('Most Essential MAPK Pathways Genes','FontSize',16)

% ***********************************************************************
% Which is the most cell cancer type that is dependant on mapk signalling
% ***********************************************************************

% get the dependency score for each cancer type
mostEssential = mostEssentialOG ;
mostEssential4Cancer = addvars( crispr(:,1), sum(mostEssential,2), ...
    'NewVariableNames','CellLineScore');

% add the cancer types and clean up the variables
mostEssential4Cancer = innerjoin(...
    sampleInfo(:,{'disease','cell_line'}),mostEssential4Cancer);

% plot a figure the essentiallity in different cancer types
figure()
boxplot(mostEssential4Cancer.CellLineScore,mostEssential4Cancer.disease)
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold' ,'XTickLabelRotation',45)

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)

% add the labels
ylabel('Number of Genes')
xlabel('Cancer Types')
title('Between Cancer Type Dependence Score','FontSize',16)

% =========================== Response to drugs =====================

% ************************************************************************
% Do cell line with the higher dependence score response the most to drug
% which target the MAPK signalling pathway
% ************************************************************************

% load in the drug response data of the cell lines
fprintf('\n Reading GDSC Screened Compounds Data\n ')
gdsc1Drugs = readtable('GDSC1 drugs.csv');
gdsc2Drugs = readtable('GDSC2 drugs.csv');
gdscDrugs = vertcat(gdsc1Drugs, gdsc2Drugs);

% return only the MAPK targeting drugs
gdscDrugs.pathway_name = categorical(gdscDrugs.pathway_name);
gdscDrugs = gdscDrugs(gdscDrugs.pathway_name == 'ERK MAPK signaling' |...
   gdscDrugs.pathway_name == 'JNK and p38 signaling' ,: ) ;

% what drugs do we have
gdscDrugs.targets = categorical(gdscDrugs.targets);

figure()
histogram(gdscDrugs.targets ,'FaceColor',[0.4660 0.6740 0.1880] )
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold' ,'XTickLabelRotation',45)
ylabel('Number of Drugs')
xlabel('Drug Targets')
title('Anti ERK MAPK Signalling Drugs','FontSize',14)

fprintf('\n Reading GDSC Fitted Dose Response Curve Data\n ')
gdsc1 = readtable('GDSC1_fitted_dose_response_15Oct19.xlsx');
gdsc2 = readtable('GDSC2_fitted_dose_response_15Oct19.xlsx');
gdscDoseResponse = vertcat(gdsc1, gdsc2);

% remove some variables
clear gdsc1Drugs gdsc2Drugs gdsc1 gdsc2

% return only drugs that targets the MAPK signalling pathway
gdscDoseResponse.PATHWAY_NAME = categorical(gdscDoseResponse.PATHWAY_NAME);
gdscDoseResponse = gdscDoseResponse( ...
    gdscDoseResponse.PATHWAY_NAME == 'ERK MAPK signaling' | ...
    gdscDoseResponse.PATHWAY_NAME == 'JNK and p38 signaling' , :) ;

% also return only gdsc drug data for cell line that also have drug
% essentiallity data within the project Achilles
% unlike the GDSC cell line those in the Achilles project  do not have
% dashes in the names: do away with those
gdscDoseResponse.CELL_LINE_NAME  = ...
    strrep(gdscDoseResponse.CELL_LINE_NAME ,'-','') ;
gdscDoseResponse.Properties.VariableNames(5) = "cell_line";
gdscDoseResponse = innerjoin( mostEssential4Cancer, gdscDoseResponse) ;

% complete calculating the dependence ranks for cancer types and add those
% to the gdsc dose response table
% now remove the cell lines for the data
mostEssential4Cancer = removevars(mostEssential4Cancer, 'cell_line');
mostEssential4Cancer.disease = categorical(mostEssential4Cancer.disease);

% get the summary statistics
mostEssential4Cancer = groupsummary( ...
    mostEssential4Cancer,'disease','median');
mostEssential4Cancer.Properties.VariableNames(3) = "medianCancerTypeScore" ;
mostEssential4Cancer = sortrows(mostEssential4Cancer,...
    'medianCancerTypeScore','descend');

% add to the gdscDoseResponse data
gdscDoseResponse.disease = categorical(gdscDoseResponse.disease);
gdscDoseResponse = innerjoin(mostEssential4Cancer(:,[1,end]), ...
    gdscDoseResponse) ;

% save the results to the supplementary file 
writetable( gdscDoseResponse, 'Supplementary File 4.xlsx', 'Sheet',...
    'MAPK Drug GDSC Dose Response')

%% Find Differenes in the Dose Response of Cell lines

% ************************************************************************
% are these dependence scores related to response of these cancer cell
% lines to various map kinase inhibitors
% ************************************************************************

% ========= between cell line responses to MAPK inhibitors ========
% first classify the cell lines into two groups: those with higher MAPK
% dependence scores and those with lower MAPK dependence scores then
% compared the dose response profiles of these cancer cell line to the
% various anticancer drugs

% add the cell line rank variable to the table
gdscDoseResponse = addvars( gdscDoseResponse, ...
    gdscDoseResponse.CellLineScore > ...
    median(gdscDoseResponse.CellLineScore), 'After', 'cell_line', ...
    'NewVariableNames','CellLineBaseRank');

% ************************************************************************
% are cell line with higher depedence score significantly more responpsive
% to MAPK pathway targeting drugs
% ************************************************************************

% compare the dose response of the cell line with high and lower ranks for
% all the MAPK drugs
% perform a ttest with unequal variable assumed
[~,p,ci,stats] = ttest2( ...
    gdscDoseResponse.LN_IC50(gdscDoseResponse.CellLineBaseRank == 1),...
    gdscDoseResponse.LN_IC50(gdscDoseResponse.CellLineBaseRank == 0),...
    'Vartype','unequal') ;

% get the mean values for the cell lines with metabolic alterations
meanHighRank = mean(...
    gdscDoseResponse.LN_IC50(gdscDoseResponse.CellLineBaseRank == 1) ) ;

meanLowRank = mean(...
    gdscDoseResponse.LN_IC50(gdscDoseResponse.CellLineBaseRank == 0) );

% produce the boxplot for the comparision
figure()
boxplot(gdscDoseResponse.LN_IC50,gdscDoseResponse.CellLineBaseRank  ) ;
hold on

% set some figure properties and add title ot the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold','XTick',1:2, 'XTickLabel', {'Lower','Higher'})
ylabel('Log IC_5_0')
xlabel('Cell Line MAPK Dependence Groups')
% title('Dose Response Variations','FontSize',16)

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)

% add color to the box plots
color = [ 0.1940 0.3840 0.9560 ;0.9900 0.2250 0.0980 ] ;
h4 = findobj(gca,'Tag','Box') ;
for jj=1:length(h4)
    patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
        color(jj,:),'FaceAlpha', .8 ,'LineWidth',1);
end
text(0.5, 0.95, strcat("p = ", convertPValue2SuperScript(p) ), ...
    'FontWeight','bold','FontSize',14 ,'Units','normalize' , ...
    'HorizontalAlignment','Center')

hold off

% #######################################################################
% Those with lower MAPK dependencies are significantly less responsive to
% MAPK targeting drugs.
% #######################################################################

% add to the table
resposeResults = array2table([meanHighRank ,meanLowRank, ...
    stats.tstat, ci(1), ci(2) , p] ,'VariableNames',...
    {'meanHighRank','meanLowRank','tValue','upperBound',...
    'lowerBound','pValue'}) ;

% create a struct for these results that will later be saved to a
% supplementary file
doseResponseResults.cellLineAllMAPK = resposeResults ;


% save the results to the supplementary file 
writetable( resposeResults, 'Supplementary File 4.xlsx', 'Sheet',...
    'Response High vs Low Dep Score')

%%
% ************************************************************************
% are cell line with higher depedence score significantly more responpsive
% to MAPK pathway targeting drugs
% ************************************************************************

% convert the drug names to categorical values
gdscDoseResponse.DRUG_NAME  = categorical(gdscDoseResponse.DRUG_NAME);

% perform t-test to compare the drug of responose of cells line that have
% higher MAPK dependence scores vs those with lower MAPK dependence scores
appendTable = unique(gdscDoseResponse(:,{'DRUG_NAME'})) ;

for ii = 1:height(appendTable)
    % get the name of the drug
    curDrug = appendTable.DRUG_NAME(ii) ;
    
    % perform a ttest with unequal variable assumed
    [~,p,ci,stats] = ttest2( ...
        filloutliers( gdscDoseResponse.Z_SCORE( ...
        gdscDoseResponse.CellLineBaseRank== true & ...
        gdscDoseResponse.DRUG_NAME == curDrug), 'linear'), ...
        filloutliers( gdscDoseResponse.Z_SCORE( ...
        gdscDoseResponse.CellLineBaseRank  == false & ...
        gdscDoseResponse.DRUG_NAME == curDrug),'linear') ,...
        'Vartype','unequal') ;
    
    % get the mean values for the cell lines with metabolic alterations
    meanHighRank = mean(...
        gdscDoseResponse.Z_SCORE(...
        gdscDoseResponse.CellLineBaseRank == true & ...
        gdscDoseResponse.DRUG_NAME == curDrug)) ;
    
    meanLowerRank = mean(...
        gdscDoseResponse.Z_SCORE(...
        gdscDoseResponse.CellLineBaseRank == false & ...
        gdscDoseResponse.DRUG_NAME == curDrug)) ;
    
    % add to the table
    appendTable(ii,2:7) = num2cell([meanHighRank,meanLowerRank, ...
        stats.tstat,ci(1), ci(2) , p]);
end

% add the variable names to the table
appendTable.Properties.VariableNames(2:7) = ...
    {'meanHigherRank','meanLowerRank','tValue','lowerBound',...
    'upperBound','pValue'} ;
resposeResults = addvars( appendTable, ...
    mafdr(appendTable.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;

resposeResults = sortrows(resposeResults,'pValue','ascend');

% create a struct for these results that will later be saved to a
% supplementary file
doseResponseResults.cellLineEachMAPK = resposeResults ;

% save the results to the supplementary file 
writetable( resposeResults, 'Supplementary File 4.xlsx', 'Sheet',...
    'Between Cell Line Response')

% #######################################################################
% It turns out that indeed those cell lines with lower MAPK ranks tend to
% repond poorly to various MAPK inhibitors compared to those with higher
% MAPK dependence scores
% #######################################################################

% plot a figure for this comparison
% get the for most signifcation pathways between the high metabolism and

% add the metabolic status of the cell lines to be taken from the Drug
% respose data and join the two tables

% ====== Produce the BAR graph that Darren suggests for Drug Action ======
forDrugsPlot = gdscDoseResponse(:,...
    {'CellLineBaseRank','DRUG_NAME','Z_SCORE'}) ;
forDrugsPlot.DRUG_NAME = categorical(forDrugsPlot.DRUG_NAME) ;

% get the median values of for each drug
drugBars = grpstats(forDrugsPlot,{'CellLineBaseRank','DRUG_NAME'},...
    {'mean','sem'} ,'DataVars','Z_SCORE') ;

% arrange the drugs according to the level of significance
try 
    drugBars = innerjoin(drugBars, resposeResults(:, [1,end-1] ) );
catch
    drugBars = innerjoin(drugBars, resposeResults(:, [2,end-1] ) ) ;
end
drugBars = sortrows(drugBars, 'pValue','ascend');

% get the data that contains nonmean
meanVars = drugBars{ : , contains(drugBars.Properties.VariableNames,...
    'mean')} ;
semVars = drugBars{ : , contains(drugBars.Properties.VariableNames,...
    'sem')} ;

% reshape the variable so that they could be used to plot the grouped bar
% graphs
% reshape the variable so that they could be used to plot the grouped bar
% graphs
meanVars = reshape(meanVars, 2, length(meanVars)/2 );
semVars = reshape(semVars, 2, length(semVars)/2);

crisprColor = [0.9900 0.2250 0.0980;0.1940 0.3840 0.9560 ] ;
mostSigPath = cellstr(resposeResults.DRUG_NAME);

figure();
% clf ;set(gcf,'position',[500,200,900,500]);
hold on
hb = bar(1:length(mostSigPath),meanVars','BarWidth',1 ) ;

% change the color of the second bars graphs
% add color to the box plots
hb(1).FaceColor = crisprColor(1,:);
hb(2).FaceColor = crisprColor(2,:);

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for kk = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(kk).XData+hb(kk).XOffset;
    er = errorbar(xData, meanVars(kk,:), semVars(kk,:),'k.');
    er.Color = crisprColor(kk,:) ;
end
ylabel('IC_5_0  Z-score')
title('Between Cell Line Dose-Response')
% add the tick and ticklabel
set(gca,'XTick',[1:length(mostSigPath)] , ...
    'XTickLabel', regexprep( mostSigPath ,' (+\w*+\)','') ,...
    'XTickLabelRotation',90,'LineWidth',1,'FontSize',12, ...
    'FontWeight','bold')

% add the stars for pvalues
for kk = 1:numel(mostSigPath)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    % get the location of the first stars
    xVar = kk ;
    
    % make the start of the stars different for the variable
    if max(meanVars(:,kk)) > 1
        yVar = max(meanVars(:,kk)) + max(semVars(:,kk)) + 0.01 ;
    else
        yVar = max(meanVars(:,kk)) + max(semVars(:,kk)) +0.01 ;
    end
    
    % specify the distance between the star
    % shift = 0.10 ;
    shift = max(meanVars,[],'all')/6 ;
    
    % change the shift if the values are both negative
    if yVar < 0
        yVar = 0.02 ;
    end
    
    % get the p values and specificy the number of star to plot at that
    % points on the chart
    pvalue = resposeResults.pValue(kk) ;
    
    if pvalue < 0.001
        % plot three stars
        lineData = [xVar yVar; xVar yVar + shift; xVar yVar + 2*shift] ;
    elseif pvalue < 0.01
        % plot two stars
        lineData = [xVar, yVar; xVar yVar + shift] ;
    elseif pvalue < 0.05
        lineData = [xVar, yVar];
        % plot two stars
    else
        % plot nothing
        continue;
    end
    
    % finally produce the line plot
    line(lineData(:,1),lineData(:,2),'LineStyle','none',...
        'Marker','pentagram','MarkerSize',10,'Color','k',...
        'MarkerFaceColor','k')
end

fontSizes = [11,10];
rectAndTextBox = [0.01 ,0.12] ;
createLegendInternal(0.85, 0.75, {'Lower','Higher'} , crisprColor,...
    'MAPK Dependence', fontSizes ,rectAndTextBox )

hold off

clear hb lineData xVar yVar shift kk er xData mostSigPath semVars ...
    meanVars drugBars forDrugsPlot appendTable meanLowerRank ...
    meanHighRank p ci stats ii lgd ll jj locA locB locGenes LocThem ...
    missingGenes meanLowRank pValue these X curPathwayCrispr ...
    color corData barData annotationData rectAndTextBox fontSizes

%% compare the dose response of the cell line with high and lower ranks for

% ************************************************************************
% Is there in difference in the drug action of MAPK targeting drug between
% tumours of that have difference dependence scores
% ************************************************************************

% add a new variable to the compare between cancer with high dependence
% score and lower dependence scores
try
    gdscDoseResponse = addvars(gdscDoseResponse , ...
        gdscDoseResponse.medianCancerTypeScore > median( ...
        gdscDoseResponse.medianCancerTypeScore), 'Before','cell_line',...
        'NewVariableNames','medianCancerTypeRank') ;
catch
    fprintf('\n medianCancerTypeRank has already been added \n')
end

% compare the dose response of the cell line with high and lower ranks for
% all the MAPK drugs
% perform a ttest with unequal variable assumed
[~,p,ci,stats] = ttest2( ...
    gdscDoseResponse.Z_SCORE(gdscDoseResponse.medianCancerTypeRank == 1),...
    gdscDoseResponse.Z_SCORE(gdscDoseResponse.medianCancerTypeRank == 0),...
    'Vartype','unequal') ;

% get the mean values for the cell lines with metabolic alterations
meanHighRank = mean(...
    gdscDoseResponse.Z_SCORE(gdscDoseResponse.medianCancerTypeRank == 1) );

meanLowRank = mean(...
    gdscDoseResponse.Z_SCORE(gdscDoseResponse.medianCancerTypeRank == 0) );

% produce the boxplot for the comparision
figure()
boxplot(gdscDoseResponse.Z_SCORE,gdscDoseResponse.medianCancerTypeRank  ) ;
hold on

% set some figure properties and add title ot the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold','XTick',1:2, 'XTickLabel', ...
    {'Lower Dependence','Higher Dependence'})
ylabel('IC_5_0 Z-score')
xlabel('Cell Line Groups')
title('Between Cancer Dose Response Variations','FontSize',14)

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)

% add color to the box plots
color = [ 0.1940 0.3840 0.9560 ;0.9900 0.2250 0.0980 ] ;
h4 = findobj(gca,'Tag','Box') ;
for jj=1:length(h4)
    patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
        color(jj,:),'FaceAlpha', .8 ,'LineWidth',1);
end
text(0.5, 0.95, strcat("p = ", convertPValue2SuperScript(p) ), ...
    'FontWeight','bold','FontSize',14 ,'Units','normalize' , ...
    'HorizontalAlignment','Center')

hold off

% #######################################################################
% Those with lower MAPK dependencies are significantly less responsive to
% MAPK targeting drugs.
% #######################################################################

% add to the table
resposeResults = array2table([meanHighRank ,meanLowRank, ...
    stats.tstat, ci(1), ci(2) , p] ,'VariableNames',...
    {'meanHighRank','meanLowRank','tValue','upperBound',...
    'lowerBound','pValue'}) ;

% create a struct for these results that will later be saved to a
% supplementary file
doseResponseResults.CancerAllMAPK = resposeResults ;

% save the results to the supplementary file 
writetable( resposeResults, 'Supplementary File 4.xlsx', 'Sheet',...
    'Between Cancer Drug Response')

%% between cancer dose response differences for each MAPK drugs

% ************************************************************************
% are cell line with higher depedence score significantly more responpsive
% to MAPK pathway targeting drugs
% ************************************************************************

% perform t-test to compare the drug of responose of cells line that have
% higher MAPK dependence scores vs those with lower MAPK dependence scores
appendTable = unique(gdscDoseResponse(:,[14,15])) ;

for ii = 1:height(appendTable)
    % get the name of the drug
    curDrug = appendTable.DRUG_NAME(ii) ;
    
    % perform a ttest with unequal variable assumed
    [~,p,ci,stats] = ttest2( ...
        filloutliers( gdscDoseResponse.Z_SCORE( ...
        gdscDoseResponse.medianCancerTypeRank == true & ...
        gdscDoseResponse.DRUG_NAME == curDrug), 'linear'), ...
        filloutliers( gdscDoseResponse.Z_SCORE( ...
        gdscDoseResponse.medianCancerTypeRank  == false & ...
        gdscDoseResponse.DRUG_NAME == curDrug),'linear') ,...
        'Vartype','unequal') ;
    
    % get the mean values for the cell lines with metabolic alterations
    meanHighRank = mean(...
        gdscDoseResponse.Z_SCORE(...
        gdscDoseResponse.medianCancerTypeRank == true & ...
        gdscDoseResponse.DRUG_NAME == curDrug)) ;
    
    meanLowerRank = mean(...
        gdscDoseResponse.Z_SCORE(...
        gdscDoseResponse.medianCancerTypeRank == false & ...
        gdscDoseResponse.DRUG_NAME == curDrug)) ;
    
    % add to the table
    appendTable(ii,3:8) = num2cell([meanHighRank,meanLowerRank, ...
        stats.tstat,ci(1), ci(2) , p]);
end

% add the variable names to the table
appendTable.Properties.VariableNames(3:8) = ...
    {'meanHigherRank','meanLowerRank','tValue','lowerBound',...
    'upperBound','pValue'} ;
resposeResults = addvars( appendTable, ...
    mafdr(appendTable.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;

resposeResults = sortrows(resposeResults,'pValue','ascend');

% create a struct for these results that will later be saved to a
% supplementary file
doseResponseResults.CancerEachMAPK = resposeResults ;

% save the results to the supplementary file 
writetable( resposeResults, 'Supplementary File 4.xlsx', 'Sheet',...
    'Each Cancer Drug Response')

% #######################################################################
% It turns out that indeed those cell lines with lower MAPK ranks tend to
% repond poorly to various MAPK inhibitors compared to those with higher
% MAPK dependence scores
% #######################################################################

% plot a figure for this comparison
% get the for most signifcation pathways between the high metabolism and

% add the metabolic status of the cell lines to be taken from the Drug
% respose data and join the two tables

% ====== Produce the BAR graph that Darren suggests for Drug Action ======
forDrugsPlot = gdscDoseResponse(:,...
    {'medianCancerTypeRank','DRUG_NAME','Z_SCORE'}) ;
forDrugsPlot.DRUG_NAME = categorical(forDrugsPlot.DRUG_NAME) ;

% get the median values of for each drug
drugBars = grpstats(forDrugsPlot,{'medianCancerTypeRank','DRUG_NAME'},...
    {'mean','sem'} ,'DataVars','Z_SCORE') ;

% arrange the drugs according to the level of significance
drugBars = innerjoin(drugBars, resposeResults(:, [1,end-1] ) );
drugBars = sortrows(drugBars, 'pValue','ascend');

% get the data that contains nonmean
meanVars = drugBars{ : , contains(drugBars.Properties.VariableNames,...
    'mean')} ;
semVars = drugBars{ : , contains(drugBars.Properties.VariableNames,...
    'sem')} ;

% reshape the variable so that they could be used to plot the grouped bar
% graphs
meanVars = reshape(meanVars, 2, length(meanVars)/2 );
semVars = reshape(semVars, 2, length(semVars)/2);

crisprColor = [0.9900 0.2250 0.0980;0.1940 0.3840 0.9560 ] ;
mostSigPath = cellstr(resposeResults.DRUG_NAME);

figure();
% clf ;set(gcf,'position',[500,200,900,500]);
hold on
hb = bar(1:length(mostSigPath),meanVars','BarWidth',1) ;

% change the color of the second bars graphs
% add color to the box plots
hb(1).FaceColor = crisprColor(1,:);
hb(2).FaceColor = crisprColor(2,:);

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for kk = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(kk).XData+hb(kk).XOffset;
    er = errorbar(xData, meanVars(kk,:), semVars(kk,:),'k.');
    er.Color = crisprColor(kk,:) ;
end
ylabel('IC_5_0  Z-score')
title('Between Cancer Dose-Response')
% add the tick and ticklabel
set(gca,'XTick', [1:length(mostSigPath)] , ...
    'XTickLabel', regexprep( mostSigPath ,' (+\w*+\)','') ,...
    'XTickLabelRotation',90,'LineWidth',1,'FontSize',12, ...
    'FontWeight','bold')

% add the stars for pvalues
for kk = 1:numel(mostSigPath)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    % get the location of the first stars
    xVar = kk ;
    
    % make the start of the stars different for the variable
    if max(meanVars(:,kk)) > 1
        yVar = max(meanVars(:,kk)) + max(semVars(:,kk)) + 0.01 ;
    else
        yVar = max(meanVars(:,kk)) + max(semVars(:,kk)) +0.01 ;
    end
    
    % specify the distance between the star
    % shift = 0.10 ;
    shift = max(meanVars,[],'all')/6 ;
    
    % change the shift if the values are both negative
    if yVar < 0
        yVar = 0.02 ;
    end
    
    % get the p values and specificy the number of star to plot at that
    % points on the chart
    pvalue = resposeResults.pValue(kk) ;
    
    if pvalue < 0.001
        % plot three stars
        lineData = [xVar yVar; xVar yVar + shift; xVar yVar + 2*shift] ;
    elseif pvalue < 0.01
        % plot two stars
        lineData = [xVar, yVar; xVar yVar + shift] ;
    elseif pvalue < 0.05
        lineData = [xVar, yVar];
        % plot two stars
    else
        % plot nothing
        continue;
    end
    
    % finally produce the line plot
    line(lineData(:,1),lineData(:,2),'LineStyle','none',...
        'Marker','pentagram','MarkerSize',10,'Color','k',...
        'MarkerFaceColor','k')
end

% add a legend to the figure
fontSizes = [11,10];
rectAndTextBox = [0.01 ,0.12] ;
createLegendInternal(0.85, 0.75, {'Lower','Higher'} , crisprColor,...
    'MAPK Dependence', fontSizes ,rectAndTextBox )

hold off

clear hb lineData xVar yVar shift kk er xData mostSigPath semVars ...
    meanVars drugBars forDrugsPlot appendTable meanLowerRank ...
    meanHighRank p ci stats ii lgd ll jj locA locB locGenes LocThem ...
    missingGenes meanLowRank pValue these X curPathwayCrispr ...
    color corData barData annotationData heatData heatDataAll pValue

% HAVE A REPLOT OF THE MAPK PATHWAY WITH EACH GENE COLORS BY THE NUMBER OF
% CANCER CELL LINES WHICH ARE DEPENDANT ON IT. THE ROW BY THE DEPENDANCE
% SCORE OF EACH TYPE OF PROTEIN AND THE COLUMNS BY THE DEPENDANCE OF THE
% PATHWAY
%% GDSC Response of cell line that are: 

% 1) Have a mutated oncogenes verses those with no mutations 
% 2) Highly dependent on an oncogene vs those that are not
fprintf('\n Checking the Dose response between CRISPR genes \n')

% perform t-test to compare the drug of responose of cells line that have
% higher dependency on one genes vs those with lower MAPK dependence scores
gdscDrugs = unique(gdscDoseResponse(:,[13,14])) ;
appendTable = array2table( crispr.Properties.VariableNames(2:end)' , ...
    'VariableName',{'HugoSymbol'} ) ;
gdscCrisprGenesDiff = [] ;

for ii = 1:height(gdscDrugs)
    % get the name of the drug
    curDrug = gdscDrugs.DRUG_NAME(ii) ;
    
    fprintf('\n Running loop for drug %s number %d of %d \n', ...
       string(curDrug) , ii, height(gdscDrugs) ) 
    
    % make a copy for the append table during this loop iteration
    appendTableCopy = appendTable ;
    
    % loop across all the cancer cell lines
    for jj = 2:width(crispr)
        % get the current cell lines that have a high dependency on the the
        depCellLines = crispr.cell_line( crispr.(jj) < -0.5 ) ;
        
        % check that there atleast more than 5 cell lines are present in
        % each group
        if length(depCellLines) < 5 || ...
                width(crispr) - length(depCellLines) < 5
            continue
        end
        
        % perform a ttest with unequal variable assumed
        [~,p,ci,stats] = ttest2( ...
            filloutliers( gdscDoseResponse.Z_SCORE( ...
            ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug), 'linear'), ...
            filloutliers( gdscDoseResponse.Z_SCORE( ...
            ~ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug),'linear') ,...
            'Vartype','unequal') ;
        
        % get the mean values for the cell lines with metabolic alterations
        meanHighRank = mean(...
            gdscDoseResponse.Z_SCORE(...
            ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug)) ;
        
        meanLowerRank = mean(...
            gdscDoseResponse.Z_SCORE(...
            ~ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug)) ;
        
        % make a copy of the append table and add values to it
        appendTableCopy(jj-1,2:8) = [ cellstr(curDrug) , ...
            num2cell([meanHighRank,meanLowerRank,...
            stats.tstat,ci(1), ci(2) , p]) ];
        
    end
    
    % delete the empty rows of the append table copy no durg
    appendTableCopy( cellfun(@isempty, appendTableCopy{:,2}) , :)  = [] ;
    
    % add to the final table
    gdscCrisprGenesDiff = [gdscCrisprGenesDiff; appendTableCopy] ;
end

% add the variable names and sort the rows
gdscCrisprGenesDiff.Properties.VariableNames(2:8) = ...
    {'Drug','meanDoseRespHigherDep','meanDoseRespLowDep','tStat', ...
    'lowerBound','upperBound','pValue'} ;
gdscCrisprGenesDiff = addvars( gdscCrisprGenesDiff, ...
    mafdr( gdscCrisprGenesDiff.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;
gdscCrisprGenesDiff = movevars(gdscCrisprGenesDiff ,'Drug','Before', ...
    'HugoSymbol') ;
gdscCrisprGenesDiff(isnan(gdscCrisprGenesDiff.pValue), :) = [] ;
gdscCrisprGenesDiff = unique(gdscCrisprGenesDiff) ;
gdscCrisprGenesDiff = sortrows(gdscCrisprGenesDiff,'pValue','ascend');
gdscCrisprGenesDiff.HugoSymbol = ...
    categorical(gdscCrisprGenesDiff.HugoSymbol);
gdscCrisprGenesDiff.Drug = categorical(gdscCrisprGenesDiff.Drug);

% save the results to the supplementary file 
writetable( gdscCrisprGenesDiff, 'Supplementary File 4.xlsx', 'Sheet',...
    'CRISPR Gene Drug Response')

%% Find the difference in the dose Response for Mutants

fprintf('\n Checking the Dose response between mutants and non-mutants \n')
% perform t-test to compare the drug of responose of cells line that have
% mutation vs those with no mutations
gdscDrugs = unique(gdscDoseResponse(:,[13,14])) ;
appendTable = array2table( CCLEmutations.Properties.VariableNames(3:end)' , ...
    'VariableName',{'HugoSymbol'} ) ;
gdscMutGenesDiff = [] ;

for ii = 1:height(gdscDrugs)
    % get the name of the drug
    curDrug = gdscDrugs.DRUG_NAME(ii) ;
    
    % make a copy for the append table during this loop iteration
    appendTableCopy = appendTable ;
    
    % loop across all the cancer cell lines
    for jj = 2:width(crispr)
        % get the current cell lines that have a high dependency on the the
        
        % it turn out that a few genes have NaN mutation so I take care of
        % this 
        if isnumeric(CCLEmutations.(jj))
            continue
        else
            depCellLines = CCLEmutations.cell_line( ...
                ~cellfun(@isempty, CCLEmutations.(jj) )  ) ;
        end
        
        % check that there atleast more than 5 cell lines are present in
        % each group
        if length(depCellLines) < 5 || ...
                width(CCLEmutations) - length(depCellLines) < 6
            continue
        end
        
        % perform a ttest with unequal variable assumed
        [~,p,ci,stats] = ttest2( ...
            filloutliers( gdscDoseResponse.Z_SCORE( ...
            ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug), 'linear'), ...
            filloutliers( gdscDoseResponse.Z_SCORE( ...
            ~ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug),'linear') ,...
            'Vartype','unequal') ;
        
        % get the mean values for the cell lines with metabolic alterations
        meanHighRank = mean(...
            gdscDoseResponse.Z_SCORE(...
            ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug)) ;
        
        meanLowerRank = mean(...
            gdscDoseResponse.Z_SCORE(...
            ~ismember(gdscDoseResponse.cell_line,depCellLines) & ...
            gdscDoseResponse.DRUG_NAME == curDrug)) ;
        
        % make a copy of the append table and add values to it
        appendTableCopy(jj-1,2:8) = [ cellstr(curDrug) , ...
            num2cell([meanHighRank,meanLowerRank,...
            stats.tstat,ci(1), ci(2) , p]) ];
        
    end
    
    % delete the empty rows of the append table copy no durg
    appendTableCopy( cellfun(@isempty, appendTableCopy{:,2}) , :)  = [] ;
    
    % add to the final table
    gdscMutGenesDiff = [gdscMutGenesDiff; appendTableCopy] ;
end

% add the variable names and sort the rows
gdscMutGenesDiff.Properties.VariableNames(2:8) = ...
    {'Drug','meanDoseRespHigherDep','meanDoseRespLowDep','tStat', ...
    'lowerBound','upperBound','pValue'} ;
gdscMutGenesDiff = addvars( gdscMutGenesDiff, ...
    mafdr( gdscMutGenesDiff.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;
gdscMutGenesDiff = movevars(gdscMutGenesDiff ,'Drug','Before', ...
    'HugoSymbol') ;
gdscMutGenesDiff(isnan(gdscMutGenesDiff.pValue), :) = [] ;
gdscMutGenesDiff = unique(gdscMutGenesDiff) ;
gdscMutGenesDiff = sortrows(gdscMutGenesDiff,'pValue','ascend');
gdscMutGenesDiff.HugoSymbol = categorical(gdscMutGenesDiff.HugoSymbol);
gdscMutGenesDiff.Drug = categorical(gdscMutGenesDiff.Drug);

% save the results to the supplementary file 
writetable( gdscMutGenesDiff, 'Supplementary File 4.xlsx', 'Sheet',...
    'Mutations Gene Drug Response')

%% ============== plot the figure for these comparisions =================

% ==================== THE INTERGRATED PLOT ==================

% specify the data that will be used to make the plot
currentToPlot = categorical(["CRISPR","Mutations"]) ;

% make a copy of the data that will be used to produe the plot
for pp = 1:length(currentToPlot)
    % get the required dataset
    if currentToPlot(pp) == "CRISPR"
        gdscCrisprGenesPlot = gdscCrisprGenesDiff( ...
            gdscCrisprGenesDiff.FDR < 0.05 ,{'Drug','HugoSymbol',...
            'tStat','pValue','FDR'} )  ;
    elseif currentToPlot(pp) == "Mutations"
        % make a copy of the data that will be used to produe the plot
        gdscCrisprGenesPlot = gdscMutGenesDiff( ...
            gdscCrisprGenesDiff.FDR < 0.05 ,{'Drug','HugoSymbol',...
            'tStat','pValue','FDR'} )  ;
    end
    
    % get only the significant results based on the FDR and also convert
    % the FDR value to either 2 for lower response or 3 for higher response
    % and 1 for not signicant This will the data compatible with the
    % oncoPrint function
    gdscCrisprGenesPlot.Value =  ( gdscCrisprGenesPlot.tStat > 0) + 2 ;
    
    % prellocate some variables
    plotTable = ones( length(unique( gdscCrisprGenesPlot.Drug )) , ...
        length( unique(gdscCrisprGenesPlot.HugoSymbol) ) ) ;
    uniqueDrugs = unique( gdscCrisprGenesPlot.Drug ) ;
    uniqueGenes = unique( gdscCrisprGenesPlot.HugoSymbol ) ;
    
    for ii = 1:length(uniqueDrugs)
        % get the data for one drug
        curData = gdscCrisprGenesPlot( ...
            ismember(gdscCrisprGenesPlot.Drug, uniqueDrugs(ii) ) , ...
            {'Drug','HugoSymbol','Value'} ) ;
        
        % find the position of the genes in the plot table
        [~,locGenes] = ismember(curData.HugoSymbol,uniqueGenes ) ;
        
        % add to the growing table
        plotTable(ii,locGenes) = curData.Value ;
    end
    
    % convert to a tables
    plotTable = addvars( array2table( plotTable , ...
        'VariableNames' , cellstr(uniqueGenes) ) , uniqueDrugs, ...
        'Before', 1 ,'NewVariableNames', {'Drugs'} ) ;
    
    % find the dependence that are associated with only resistance or
    % sensitive to MAPK inhibitors
    sensitiveOnly = plotTable( :, [true, all( plotTable{:, 2:end} ~= 3)]);
    sensitiveOnly = addvars( array2table( ...
        sensitiveOnly.Properties.VariableNames(2:end)','VariableNames', ...
        "HugoSymbol" ) , sum(sensitiveOnly{:, 2:end} > 1)', ...
        sum(sensitiveOnly{:, 2:end} > 1)'/height(sensitiveOnly)*100 , ...
        'NewVariableNames', {'numberOfDrugs','PercentageOfDrugs'} ) ;
    sensitiveOnly = sortrows(sensitiveOnly,'PercentageOfDrugs','descend')
    
    resistantOnly = plotTable( :, [true, all( plotTable{:, 2:end} ~= 2)]);
    resistantOnly = addvars( array2table( ...
        resistantOnly.Properties.VariableNames(2:end)','VariableNames', ...
        "HugoSymbol" ) , sum(resistantOnly{:, 2:end} > 1)', ...
        sum(resistantOnly{:, 2:end} > 1)'/height(resistantOnly)*100 , ...
        'NewVariableNames', {'numberOfDrugs','PercentageOfDrugs'} ) ;
    resistantOnly = sortrows(resistantOnly,'PercentageOfDrugs','descend')
    
    % specify the colors for the plot responseColors = [0.9 0.9 0.9 ;
    % 0.1000 0.2000 0.6000; 0.944 0.337 0.1112];
    responseColors = [0.9 0.9 0.9 ; 0.1000 0.6000 0.2000; 0.944 0.337 0.1112];
    
    % plot the data on an oncoprint
    if currentToPlot(pp) == "CRISPR"
        global plotTime
        plotTime = 2;
        % PLOT THE COMPLEX HEAT MAP
        axisLast = oncoPrintInternal( plotTable{:, 2:end}, responseColors,...
            cellstr(plotTable.Drugs) , true) ;
    elseif currentToPlot(pp) == "Mutations"
        global plotTime
        plotTime = 3;
        % PLOT THE COMPLEX HEAT MAP
        axisLast = oncoPrintInternal( plotTable{:, 2:end}, responseColors,...
            cellstr(plotTable.Drugs) , true) ;
    end
    
    % create an array that will be used to color the text in the plots
    colorTextArray = rows2vars(plotTable,'VariableNamesSource','Drugs');
    colorTextArray.Properties.VariableNames(1) = "HugoSymbol" ;
    colorTextArray.color(all(colorTextArray{:, 2:end} ~= 3 ,2)) = 1 ;
    colorTextArray.color(all(colorTextArray{:, 2:end} ~=2 ,2)) = -1 ;
    colorTextArray.color( all(colorTextArray{:, 2:end} == 2 ,2) & ...
        all(colorTextArray{:, 2:end} == 3 ,2)) = 0 ;
    colorTextArray = colorTextArray(:,{'HugoSymbol','color'}) ;
    responseColorsText = [0.1000 0.6000 0.2000; 0.1 0.1 0.1; ...
        0.944 0.337 0.11] ;
    % thow in an assertion
    assert( all( strcmp( colorTextArray.HugoSymbol,cellstr(uniqueGenes))))
    
    % add the gene names to the top of the plot demacate the x values into
    % two equal part based on teh number genes
    if currentToPlot(pp) == "CRISPR"
        increaseby = 0.0202 + 0.0202/10;
        xFirst = 0.011 ;
    elseif currentToPlot(pp) == "Mutations"
        increaseby = 0.0102 + 0.0104/10;
        xFirst = 0.0055 ;
    end
    
    % add the genes and colors
    for jj = 1:length(uniqueGenes)
        
        % specificy the gene label colors
        if colorTextArray.color(jj) == 1
            curColor = responseColorsText(1,:) ;
        elseif colorTextArray.color(jj) == 0
            curColor = responseColorsText(2,:) ;
        elseif colorTextArray.color(jj) == -1
            curColor = responseColorsText(3,:) ;
        end
        % add the gene names to the top of the plot
        text( xFirst, 6, string(uniqueGenes(jj)),...
            'FontSize', 10, 'HorizontalAlignment', 'right', ...
            'FontWeight', 'bold','VerticalAlignment', 'middle' ,...
            'Rotation',90 ,'Units','Normalize','Color' ,curColor );
        xFirst = xFirst + increaseby ;
    end
    
    % add the gene class to the top of the plot
    
    % get the MAPK genes and arrange them in the same order as the plot
    [~, locThem] = ismember(uniqueGenes ,mapkGenes.Gene ) ;
    topMapkGenes = mapkGenes(locThem, :) ;
    
    % throw in an assertion
    assert( all( strcmp( topMapkGenes.Gene , cellstr(uniqueGenes) ) ) )
    
    % add the number of genes that are dependent of these cell lines and
    % also the dependence score of mean
    
    % =============== prellocate some variable =========== 
    % the dimesions boxplot
    yInitial = 0.8 ;
    dims = [0.03, yInitial, 0.1, 0.025 ] ;
    % the axis position and change y position and names of the bars
    axPos = axisLast; axPos(2) = yInitial ;
    if currentToPlot(pp) == "CRISPR"
        barNames = {'Cascade','Gene Class','Dependence'} ;
    elseif currentToPlot(pp) == "Mutations"
        barNames = {'Cascade','Gene Class','Mutations'} ;
    end
    increaseBy = 0.03;
    
    % add a heatmap on top to show the supervised classifier output
    for ii = 1:length(barNames)
        axes('position', axPos);
        switch ii
            case 1
                ultraBars( double(topMapkGenes.Pathway)', ...
                    classesColors.pathway ) ;
            case 2
                ultraBars(double(topMapkGenes.proteinClass)', ...
                    classesColors.proteinClass );
            case 3
                
                % % now get the data of the cancer study and change the
                % arrangement in the cancer studies
                barData = clusterDataDep(:,3:end) ;
                [~, locCols] = ismember(topMapkGenes.Gene, ...
                    barData.Properties.VariableNames );
                barData = barData(:,locCols) ;
                
                % throw in an assertion
                assert(all(strcmp(topMapkGenes.Gene, ...
                    barData.Properties.VariableNames')))
                
                % get the values above a give threshold of the essentialies
                % scores
                barData = table2array(barData) <  -0.50;
                barData = sum(barData, 1);
                heatmap(barData,'Colormap', parula ,'GridVisible','off',...
                    'FontColor','none','ColorbarVisible','off', ...
                    'ColorLimits',[-1 75]);
                
        end
        
        % add the name of the genes to the left of heatmap
        annotation('textbox', dims ,'String', barNames{ii} ,...
            'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
            'HorizontalAlignment','right','FontWeight','bold',...
            'VerticalAlignment','middle');
        
        % increase the y axis values
        axPos(2) = axPos(2) + increaseBy ;
        dims(2) = dims(2) + increaseBy ;
    end
    % save the data to excel
    if currentToPlot(pp) == "CRISPR"
        writetable(sensitiveOnly, 'Supplementary File 4.xlsx', 'Sheet',...
            'CRISPR Drug Sensitive Genes')
        writetable(resistantOnly, 'Supplementary File 4.xlsx', 'Sheet',...
            'CRISPR Drug Resistance Genes')
        writetable(gdscCrisprGenesDiff,'Supplementary File 4.xlsx', ...
            'Sheet','CRISPR Drug Sensitive Assoc')
    elseif currentToPlot(pp) == "Mutations"
        writetable(sensitiveOnly, 'Supplementary File 4.xlsx', 'Sheet',...
            'Mutation Drug Sensitive Genes')
        writetable(resistantOnly, 'Supplementary File 4.xlsx', 'Sheet',...
            'Mutation Drug Resistance Genes')
        writetable(gdscMutGenesDiff, 'Supplementary File 4.xlsx', ...
            'Sheet','Mutation Drug Sensitive Assoc')
        
    end
end

% ===================== add a legend to the figure =================

% specify where to plot the legend and the font size and the width of the
% rectangle and the text box
figure()
yPoint = 0.90 ; xStart = 0.5 ; myLgdTitle = 'Pathway'; fontSizes = [8,10];
rectAndTextBox = [0.01 ,0.10] ;
legendLabels = categories( mapkGenes.Pathway);
% create the legend
createLegendInternal(yPoint, xStart, legendLabels , ...
    classesColors.pathway , myLgdTitle , fontSizes ,rectAndTextBox)

% add the legend to the figure
% add a legend to figure  using the create legend function
% First specify the legend labels and yPoint and xStart
legendLabels = categories( mapkGenes.proteinClass) ;
% replace Guanosine Exchange Factor with GEF and Dual Specificy Phosphatase
% withh DUSP
legendLabels = strrep(legendLabels, 'DualSpecificityPhosphatase','DUSP');
legendLabels = strrep(legendLabels, 'GuanosineExchangeFactor','GEF');

% specify where to plot the legend and the font size and the width of the
% rectangle and the text box
yPoint = 0.75 ; myLgdTitle = 'MAPK Pathway Gene Class';

% create the legend
createLegendInternal(yPoint, xStart, legendLabels , ...
    classesColors.proteinClass,myLgdTitle ,fontSizes ,rectAndTextBox)

% do the legend of the headmap colors
yPoint = 0.45 ; 
responseColorsText = [0.1000 0.6000 0.2000; 0.8 0.8 0.8; ...
    0.944 0.337 0.11] ;

createLegendInternal(yPoint, xStart, ...
    {'More Sensitive','No Difference','More Resistant'} , ...
    responseColorsText,'Drug Response', fontSizes  )

clear hb lineData xVar yVar shift kk er xData mostSigPath semVars ...
    meanVars drugBars forDrugsPlot appendTable meanLowerRank ...
    meanHighRank p ci stats ii lgd ll jj locA locB locGenes LocThem ...
    missingGenes meanLowRank pValue these X curPathwayCrispr ...
    color corData barData annotationData heatData heatDataAll pValue ...
    appendTableCopy locGenes uniqueDrugs uniqueGenes curData plotTime ...
    topMapkGenes axPos dims increaseBy  barData locCols topMapkGenes ...
    yInitial barNames locThem xFirst curColor colorTextArray ...
    responseColorsText responseColors currentToPlot pp

%% IF inhibition does not to kill the cell, is this deal to divergence of signals

% through the signaling pathway... what if we inhibit two proteins does
% this really kill the cell? LINCS has combination treatment -- we use
% those

% load the proteomics data from the links project
% GCT and GCTx files can be read in the same way.
% We'll use the same two files throughout this tutorial.
gct_file_locDrugs = fullfile(['/Users/sinkala/Documents/MATLAB/MAPK ',...
    'Pan-cancer Analysis/LINCS proteomics'], 'LINCS_P100_all_data.gct');


% LINCS proteomics/GSE101406_Broad_LINCS_P100_Level3_QCNORM_n1684x96.gctx'
% contains data of the proteomics response of 6 cancer cell lines to drug
% pertubation - I am still not sure what the columns stand for
p100DrugsData = cmapm.Pipeline.parse_gctx(gct_file_locDrugs);

% GCT data representation
% GCT and GCTx files are both represented in memory as structures.
disp(class(p100DrugsData));

% =================== Layout of the GCT structure ======================
% The GCT structure comprises of the following fields:
%
% * mat : Numeric data matrix [RxC]
% * rid : Cell array of row ids
% * rhd : Cell array of row annotation fieldnames
% * rdict : Dictionary of row annotation fieldnames to indices
% * rdesc : Cell array of row annotations
% * cid : Cell array of column ids
% * chd : Cell array of column annotation fieldnames
% * cdict : Dictionary of column annotation fieldnames to indices
% * cdesc : Cell array of column annotations

% display the p100 drug response of the cell lines
disp(p100DrugsData);

% get the data into three table:
% 1 for the drugs, 2 for the drug description and 3, for protein
% description
p100DrugCorrData = addvars( array2table( p100DrugsData.mat ,...
    'VariableNames',p100DrugsData.cid ) , p100DrugsData.rid, ...
    'Before', 1,'NewVariableNames','drugName') ;

p100RowDesc = array2table( p100DrugsData.rdesc ,'VariableNames',...
    p100DrugsData.rhd) ;
p100ColDesc = array2table( p100DrugsData.cdesc ,'VariableNames',...
    p100DrugsData.chd) ;

% return only only the cell line that are present within the GDSC database
% for the and find out if there are difference in the response of cell line
% that respond better to MAPK drug and those that respond better
theMAPKdrugs = ismember( p100RowDesc.moa , ...
    {'p38 MAPK inhibitor','JNK inhibitor','RAF inhibitor',...
    'MEK inhibitor' } ) ;
p100RowDesc = p100RowDesc(theMAPKdrugs, :);
p100DrugCorrData = p100DrugCorrData(theMAPKdrugs, :) ;

% replace the drugs with the correct names
p100DrugCorrData.Properties.VariableNames = strrep( ...
    p100DrugCorrData.Properties.VariableNames,'PD0325901','PD-0325901');
p100DrugCorrData.Properties.VariableNames = strrep( ...
    p100DrugCorrData.Properties.VariableNames,'SP600125',...
    'Pyrazolanthrone');

% do the same for the columns
theMAPKdrugs = ismember( p100ColDesc.moa , ...
    {'p38 MAPK inhibitor','JNK inhibitor','RAF inhibitor',...
    'MEK inhibitor' } ) ;
p100ColDesc = p100ColDesc(theMAPKdrugs,:) ;
p100DrugCorrData = [ p100DrugCorrData(:,1),  ...
    p100DrugCorrData(:, ...
    contains(p100DrugCorrData.Properties.VariableNames, ...
    p100ColDesc.pert_iname ,'IgnoreCase',true) )  ];

% throw in an assertions
assert( all( size(p100DrugCorrData) == ...
    [height(p100RowDesc), height(p100ColDesc)+1 ]) ,'Size mismatch')

% since p100ColDesc has repeating values, these can not be used to find
% theh actual location of the samples in a clustergram. Therefore, I insert
% a new column with just number from 1 : height of p100ColDesc
p100ColDesc.NameCellID = strrep( strcat( ...
    strcat(p100ColDesc.cell_id, ':') ,...
    p100ColDesc.name)  ,'-','');

% make sure that the rows are arranged in the same order as the columns in
p100DrugCorrData.Properties.VariableNames(1) = "pert_id" ;
p100DrugCorrData = innerjoin( p100RowDesc(:,...
    {'pert_id','pert_iname'}), p100DrugCorrData) ;
p100DrugCorrData = removevars(p100DrugCorrData,"pert_id");

% also replace the dashed in the variable names
p100DrugCorrData.Properties.VariableNames(2:end) = ...
    strrep(p100DrugCorrData.Properties.VariableNames(2:end), '-','') ;

% create a correlation matrix
% correl = corrcoef(p100DrugCorrData{:,2:end});
p100heatData = p100DrugCorrData{:,2:end} ;

% provide an option for redblue coloring
classesColors.redblue20 = redbluecmap ; % parula ; %redblue(200); %

fprintf('\n clustering the p100 LINCS data \n')
corrCGOp100 = clustergram(p100heatData,...
    'Colormap',classesColors.redblue20, ...
    'RowLabels', p100DrugCorrData.pert_iname , ...
    'ColumnLabels',p100DrugCorrData.Properties.VariableNames(:,2:end),...
    'ColumnPDist','euclidean',...
    'Linkage','complete','Dendrogram',8.6);  % correlation

% anotate the clustergram
addTitle(corrCGOp100,'Phophoprotein Drug Response')
% addYLabel(corrCGOp100,'Drugs')

% move the dendograms up to allow space to put the groups
cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
    'Color',{[1 1 0],[0.6 0.6 1]});
set(corrCGOp100,'ColumnGroupMarker',cm ,'ColumnLabelsRotate',90)
% set(corrCGOp100,'RowGroupMarker',cm )

% arrange the columns description in the same order as the sample clustered
[~,locCols] = ismember(corrCGOp100.ColumnLabels',p100ColDesc.NameCellID);
p100ColDesc = p100ColDesc(locCols,:);
assert( all(strcmp(corrCGOp100.ColumnLabels', p100ColDesc.NameCellID)) )

% do the same with the rows labels
[~,locRows] = ismember(corrCGOp100.RowLabels,p100RowDesc.pert_iname);
p100RowDesc = p100RowDesc(locRows,:);
assert( all(strcmp(corrCGOp100.RowLabels, p100RowDesc.pert_iname)) )

% convert to categorical so they could be used for the ploting using
% ultraBars
p100ColDesc.cell_id = categorical(p100ColDesc.cell_id);
p100ColDesc.moa = categorical(p100ColDesc.moa);
p100RowDesc.moa = categorical(p100RowDesc.moa);

% define the colors
% create colors for al the groups that are to be used in the top bar charts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classesColors.moa = [0.49,0.18,0.56; 01.00,0.41,0.16;...
    0.4947 0.9691 0.3367; 0.3484 0.7023 0.8022] ; %
classesColors.cell_id = rand(numel(unique(p100ColDesc.cell_id)),3) ;
% rng(1)
% classesColors.moa = rand(numel(unique(p100ColDesc.moa)),3) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First cluster the data to get the groups that are clustered that will be
% in the heatmap
plot(corrCGOp100) ;
hold on
axes('position',[0.2423 0.730 0.5760 0.0250] );
ultraBars( double(p100ColDesc.moa)', classesColors.moa, ...
    'MOA')

% % create a dendogram for the heatmap that i will produce
% axes('position',[0.196,0.80,0.608,0.15]);
% treeData = pdist(ProtExprKmeans','Correlation');
% tree = linkage(treeData,'complete');
% leafOrder = optimalleaforder(tree,pdist(ProtExprKmeans'));
% H = dendrogram(tree, 0,'Reorder',leafOrder,'ColorThreshold','default');
% set(H,'LineWidth',1.5)
% set(gca,'YTickLabel',[],'Visible','off','XDir','reverse');

% define the y point to start the plot
yPoint = 0.123809523809524 ;
yPointText = 0.115 + 0.0857/2 ;

% add the groups to the right of the clustergram
for ii = 1:height(p100RowDesc)
    % find the mode of action for the current drug and a bar to the right
    % of the with the corresponding colors
    switch p100RowDesc.moa(ii)
        case 'JNK inhibitor'
            curColor = classesColors.moa(1,:) ;
        case 'MEK inhibitor'
            curColor = classesColors.moa(2,:) ;
        case 'RAF inhibitor'
            curColor = classesColors.moa(3,:) ;
        case 'p38 MAPK inhibitor'
            curColor = classesColors.moa(4,:) ;
    end
    
    % add the legend text
    annotation('textbox',[0.83 yPointText 0.12 0.0230],...
        'String',p100RowDesc.pert_iname{ii},'FontSize',6,...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]);
    
    % add the legend color
    annotation('rectangle',[0.821 yPoint 0.010 0.0857],...
        'EdgeColor',curColor, 'FaceColor',curColor);
    
    % move the y point down
    yPoint = yPoint + 0.0857 ;
    yPointText = yPointText + 0.0857 ;
end

% add a legend to the figure
% Create textbox
annotation('textbox',[0.1023 0.863 0.16 0.0250],...
    'String','Drug Mode of Action','FontSize',8,...
    'FontName','Helvetica Neue','FitBoxToText','off',...
    'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]);

% specificy the y values starts and mode of actions for the drugs
yPoint = 0.830 ;
moa = categories(p100ColDesc.moa);

for ii = 1:length(classesColors.moa)
    % add the legend color
    annotation('rectangle',[0.1023 yPoint 0.018 0.023],...
        'EdgeColor',classesColors.moa(ii,:), ...
        'FaceColor',classesColors.moa(ii,:));
    
    % add the legend text
    annotation('textbox',[0.118 yPoint 0.12 0.0230],...
        'String',moa{ii},'FontSize',6,...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]);
    
    % move the y point down
    yPoint = yPoint - 0.03 ;
end

hold off

clear theMAPKdrugs p100DrugsData gct_file_locDrugs H leafOrder ...
    tree yInitial increaseby jj boxPos boxNumber boxColors dim ...
    barData2 xEndPos ySize toAddFirst barNames locCols locRows correl ...
    yPoint ii moa yPointText

% save the p100 data to excel 
writetable( p100DrugCorrData, 'Supplementary File 5.xlsx','Sheet', ...
    'p100 Drug Response')

%% Using only Cell lines with Data in the GDSC and Achilles Project

% now return only cell line that are present within the GDSC and Achilles
% project
theMAPKdrugs = ismember( p100ColDesc.cell_id , ...
    gdscDoseResponse.cell_line) | ismember( p100ColDesc.cell_id , ...
    crispr.cell_line) ;
p100ColDesc = p100ColDesc(theMAPKdrugs,:) ;

% return the cell line present in the column description for the
% p100DrugCorrData also
theMAPKdrugs = ismember(...
    p100DrugCorrData.Properties.VariableNames,...
    p100ColDesc.NameCellID );
p100DrugCorrData = [p100DrugCorrData(:,1),p100DrugCorrData(:,theMAPKdrugs)];

% throw in an assertions
assert( all( size(p100DrugCorrData) == ...
    [height(p100RowDesc), height(p100ColDesc)+1 ]) ,'Size mismatch')
assert(all(ismember(p100DrugCorrData.Properties.VariableNames(2:end)', ...
    p100ColDesc.NameCellID ) ) )

% get the clustering data
p100heatData = p100DrugCorrData{:,2:end} ;

fprintf('\n clustering the p100 LINCS data \n')
corrCGOp100 = clustergram(p100heatData,...
    'Colormap',classesColors.redblue20, ...
    'RowLabels', p100DrugCorrData.pert_iname , ...
    'ColumnLabels',p100DrugCorrData.Properties.VariableNames(:,2:end),...
    'ColumnPDist','seuclidean',...
    'Linkage','complete','Dendrogram',8.6);  % correlation

% anotate the clustergram
addTitle(corrCGOp100,'Phophoprotein Drug Response')
% addYLabel(corrCGOp100,'Drugs')

% move the dendograms up to allow space to put the groups
cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
    'Color',{[1 1 0],[0.6 0.6 1]});
set(corrCGOp100,'ColumnGroupMarker',cm ,'ColumnLabelsRotate',90)
% set(corrCGOp100,'RowGroupMarker',cm )

% arrange the columns description in the same order as the sample clustered
[~,locCols] = ismember(corrCGOp100.ColumnLabels',p100ColDesc.NameCellID);
p100ColDesc = p100ColDesc(locCols,:);
assert( all(strcmp(corrCGOp100.ColumnLabels', p100ColDesc.NameCellID)) )

% do the same with the rows labels
[~,locRows] = ismember(corrCGOp100.RowLabels,p100RowDesc.pert_iname);
p100RowDesc = p100RowDesc(locRows,:);
assert( all(strcmp(corrCGOp100.RowLabels, p100RowDesc.pert_iname)) )

% convert to categorical so they could be used for the ploting using
% ultraBars
p100ColDesc.cell_id = categorical(p100ColDesc.cell_id);
p100ColDesc.moa = categorical(p100ColDesc.moa);
p100RowDesc.moa = categorical(p100RowDesc.moa);

% define the colors
% create colors for al the groups that are to be used in the top bar charts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classesColors.moa = [0.49,0.18,0.56; 01.00,0.41,0.16;...
    0.4947 0.9691 0.3367; 0.3484 0.7023 0.8022] ; %
classesColors.cell_id = rand(numel(unique(p100ColDesc.cell_id)),3) ;
% rng(1)
% classesColors.moa = rand(numel(unique(p100ColDesc.moa)),3) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First cluster the data to get the groups that are clustered that will be
% in the heatmap
plot(corrCGOp100) ;
hold on
axes('position',[0.2423 0.730 0.5760 0.0250] );
ultraBars( double(p100ColDesc.moa)', classesColors.moa, ...
    'MOA')

% % create a dendogram for the heatmap that i will produce
% axes('position',[0.196,0.80,0.608,0.15]);
% treeData = pdist(ProtExprKmeans','Correlation');
% tree = linkage(treeData,'complete');
% leafOrder = optimalleaforder(tree,pdist(ProtExprKmeans'));
% H = dendrogram(tree, 0,'Reorder',leafOrder,'ColorThreshold','default');
% set(H,'LineWidth',1.5)
% set(gca,'YTickLabel',[],'Visible','off','XDir','reverse');

% define the y point to start the plot
yPoint = 0.123809523809524 ;
yIncrease = (0.725 - yPoint)/length(p100RowDesc.pert_iname);
yPointText = 0.115 + 0.0857/2 ;

% add the groups to the right of the clustergram
for ii = 1:height(p100RowDesc)
    % find the mode of action for the current drug and a bar to the right
    % of the with the corresponding colors
    switch p100RowDesc.moa(ii)
        case 'JNK inhibitor'
            curColor = classesColors.moa(1,:) ;
        case 'MEK inhibitor'
            curColor = classesColors.moa(2,:) ;
        case 'RAF inhibitor'
            curColor = classesColors.moa(3,:) ;
        case 'p38 MAPK inhibitor'
            curColor = classesColors.moa(4,:) ;
    end
    
    % add the legend text
    annotation('textbox',[0.83 yPointText 0.12 0.0230],...
        'String',p100RowDesc.pert_iname{ii},'FontSize',10,...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]);
    
    % add the legend color
    annotation('rectangle',[0.821 yPoint 0.010 yIncrease],...
        'EdgeColor',curColor, 'FaceColor',curColor);
    
    % move the y point down
    yPoint = yPoint + yIncrease ;
    yPointText = yPointText + yIncrease ;
end

% add a legend to the figure
% Create textbox
annotation('textbox',[0.1023 0.863 0.16 0.0250],...
    'String','Drug Mode of Action','FontSize',8,...
    'FontName','Helvetica Neue','FitBoxToText','off',...
    'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]);

% specificy the y values starts and mode of actions for the drugs
yPoint = 0.830 ;
moa = categories(p100ColDesc.moa);

for ii = 1:length(classesColors.moa)
    % add the legend color
    annotation('rectangle',[0.1023 yPoint 0.018 0.023],...
        'EdgeColor',classesColors.moa(ii,:), ...
        'FaceColor',classesColors.moa(ii,:));
    
    % add the legend text
    annotation('textbox',[0.118 yPoint 0.12 0.0230],...
        'String',moa{ii},'FontSize',6,...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]);
    
    % move the y point down
    yPoint = yPoint - 0.03 ;
end

hold off

clear theMAPKdrugs p100DrugsData gct_file_locDrugs H leafOrder ...
    tree yInitial increaseby jj boxPos boxNumber boxColors dim ...
    barData2 xEndPos ySize toAddFirst barNames locCols locRows correl ...
    yPoint ii moa yPointText varNames yIncrease

%% ======== INTEGRATIVE ANALYSIS OF ALL THE PERTUBATION DATASETS =======

% #######################################################################
% IT TURNS OUT THAT THE CELL LINE HAVE A HIGHLY CORRELATED RESPONSE TO MEK
% INHIBITORS AND VARIED RESPONSES TO JAK INHIBITORS -- can we use p100 data
% , L1000 data, GDSC dose response and Achilles gene dependencies to
% understand why this is the case why this is the case.

% ALSO THE RESPONSE OF THE CELL LINE TO MEK AND RAF INHIBITORS IS NOT
% HIGHLY CORREALTED IN CERTAIN INSTANCE -- can we use the above datasets to
% illuminate why this is the case?
% #######################################################################

% It turns out the GDSC and/or CCLE do not have all the required drugs
% therefore, I can not answer this question concretely. It should be
% interesting to find out why this is the case in features

%% Correlate the response of these cell line to p100 responses

% load the new lincs drug data
% load the proteomics data from the links project
% GCT and GCTx files can be read in the same way.
% We'll use the same two files throughout this tutorial.
gct_file_locDrugs = fullfile(['/Users/sinkala/Documents/MATLAB/MAPK ',...
    'Pan-cancer Analysis/LINCS proteomics'], 'LINCS_P100_all_data.gct');

% LINCS proteomics/GSE101406_Broad_LINCS_P100_Level3_QCNORM_n1684x96.gctx'
% contains data of the proteomics response of 6 cancer cell lines to drug
% pertubation - I am still not sure what the columns stand for
p100DrugsData = cmapm.Pipeline.parse_gctx(gct_file_locDrugs);

% get the data into three table:
% 1 for the drugs, 2 for the drug description and 3, for protein
% description
p100ColDesc2 = array2table( p100DrugsData.cdesc ,'VariableNames',...
    p100DrugsData.chd) ;

% do the same for the columns
theMAPKdrugs = ismember( p100ColDesc2.moa , ...
    {'p38 MAPK inhibitor','JAK inhibitor','RAF inhibitor',...
    'MEK inhibitor' } ) ;
p100ColDesc2 = p100ColDesc2(theMAPKdrugs,:) ;

% remove the astrocytes from the data 
p100ColDesc2 = p100ColDesc2( ...
    ~ismember( p100ColDesc2.cell_id , {'Astrocyte','NPC'}) , :) ;

% return only the dose resposes of these p100 cancer cell line and the p100
% drugs from the GDSC database
lincsGdscResponse =  gdscDoseResponse( ...
    ismember(gdscDoseResponse.DRUG_NAME, p100ColDesc2.name) & ...
    ismember(gdscDoseResponse.cell_line, p100ColDesc2.cell_id), :) ;

% get the response of the cell line from the CCLE and clean up the data
lincsCcleResponse = readtable('CCLE Drug Sensitivity Data.xlsx');
lincsCcleResponse.Properties.VariableNames = ...
    lincsCcleResponse{2,:} ;
lincsCcleResponse(1:2,:) = [] ;
lincsCcleResponse.Properties.VariableNames(1:2) = ...
    {'CancerType','cell_line'} ;

% return only the dose resposes of these p100 cancer cell line and the p100
% drugs from the CCLE database
lincsCcleResponse = lincsCcleResponse( ...
    ismember(lincsCcleResponse.Compound, p100ColDesc.name) & ...
    ismember( strrep( strrep(lincsCcleResponse.cell_line,'-',''), ...
      ' ',''), p100ColDesc.cell_id), :);

% change the variable names
lincsCcleResponse.Properties.VariableNames{5} = 'Doses_uM';
lincsCcleResponse.Properties.VariableNames{6} = 'ActivityDataMedian';
lincsCcleResponse.Properties.VariableNames{7} = 'ActivitySD';
lincsCcleResponse.Properties.VariableNames{11} = 'IC50';

% plot the dose response curves for the two drugs
figure()
for ii = 1:height(lincsCcleResponse)
    curDose = str2double( split( lincsCcleResponse.Doses_uM(ii) ,',') );
    curResponse = str2double( split( ...
        lincsCcleResponse.ActivityDataMedian(ii),',' ) );
    ActivitySD = str2double( split( lincsCcleResponse.ActivitySD(ii),','));
    
    % use the dose response function to produce the plo
    hold on
    plot(log(curDose),curResponse,'LineWidth',2 ,...
        'MarkerSize',40,'Marker','.') ;
    
%     % add the error bars to the plot
%     er = errorbar(log(curDose), curResponse, log(ActivitySD) );
end
% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',14,'Box','off','FontWeight','bold')
xlabel('Log Dose (uM)') ;
ylabel('Activity Response')
lgdd = strcat( lincsCcleResponse.cell_line , strcat( ...
    strcat(": IC50 = ",lincsCcleResponse.IC50) ,'(µM)') ) ;
legend(lgdd','Location','best') ;
title('\bf Response to PD-0325901','FontSize',16)
hold off

% plot the response of MCF7 and A549 cell line to Selumetinib from the GDSC
% database (either a dose response curve or just the avarage response

% get only the two cell line data. Because Selumetinib is repeated. Get
% only one example of it
lincsGdscCommon = unique( lincsGdscResponse( ...
    ismember( lincsGdscResponse.cell_line ,...
    unique(p100ColDesc2.cell_id) ), :)) ;
lincsGdscCommon = lincsGdscCommon(2:end, :)  ;
lincsGdscCommon.cell_line = categorical(lincsGdscCommon.cell_line);

% return only the unique rows 
[~, locUnique ] = unique(lincsGdscCommon.cell_line) ;
lincsGdscCommon = lincsGdscCommon(locUnique, :) ;

% produce a bar graph
figure();
bar(lincsGdscCommon.cell_line, lincsGdscCommon.LN_IC50 ,...
    'BarWidth',0.7,'FaceColor',[0.5 0.7 0.2]) ;
ylabel('Log IC_5_0 ')
title('Response to Selumetinib','FontSize',15)
% add the tick and ticklabel
set(gca,'LineWidth',1,'FontWeight','bold','Box','off','YLim',[-2 3.5])

% *** the RAW data does not have drug IDS or dose response profiles ****
% load the GDSC dose response values and return only the two cell line
% data = readtable('GDSC1_public_raw_data_15Oct19.csv');
% data = data(ismember(data.CELL_LINE_NAME,lincsGdscTwo.cell_line), :) ;
%
% % also only return NA and my drug concentrations
% data2 = data(ismember(data.DRUG_ID,{'Selumetinib','PD-0325901'} ), :);

% ************************************************************************
% Here just as there is a high correlation between the p100 response of
% MCF7 to Selumetinib and PD-0325901 the MCF7 cell line are resistant to
% both of these MEK1 inhibitors. Also the A549 cell are hypersensitive to
% both Selumetinib and PD-0325901. Need to check the response p100 response
% of the MCF7 and A549 cell line
% ************************************************************************

% get the p100 response of these cell line to to PD-0325901 from the level
% 3 data lincs data
% load the proteomics data from the links project
gct_file_location = fullfile(['/Users/sinkala/Documents/MATLAB/MAPK ',...
    'Pan-cancer Analysis/LINCS proteomics'],  ...
    'GSE101406_Broad_LINCS_P100_Level3_QCNORM_n1684x96.gctx');

% LINCS proteomics/GSE101406_Broad_LINCS_P100_Level3_QCNORM_n1684x96.gctx'
% contains data of the proteomics response of 6 cancer cell lines to drug
% pertubation - I am still not sure what the columns stand for
p100Level3All = cmapm.Pipeline.parse_gctx(gct_file_location);

% GCT data representation
% GCT and GCTx files are both represented in memory as structures.
disp(class(p100Level3All));
disp(p100Level3All);

% get the data into three table:
% 1 for the drugs, 2 for the drug description and 3, for protein
% description
p100Level3 = addvars( array2table( p100Level3All.mat , 'VariableNames',...
    p100Level3All.cid ) , p100Level3All.rid, 'Before', 1, ...
    'NewVariableNames','proteins') ;

% load the pertabugen info from the lincs project and return only the two
% cell line that I want to check for and the two drugs plus the untreated
% DMSO to be used for baseline
p100L3Drugs = readtable('GSE101406_Broad_LINCS_P100_inst_info.txt');
p100L3Drugs.pert_iname = categorical(p100L3Drugs.pert_iname);
p100L3Drugs = p100L3Drugs( ...
    ismember(p100L3Drugs.cell_id, lincsGdscCommon.cell_line) & ...
    ismember(p100L3Drugs.pert_iname, ...
    [ unique(p100ColDesc.pert_iname) ;{'DMSO'} ]   ) , :) ;

% now only get the columns of the p100Level3 data that has the sample names
% contained in the p100L3Drugs
p100Level3 = [p100Level3(:,1) , p100Level3(:, ismember(...
    p100Level3.Properties.VariableNames, p100L3Drugs.inst_id  ) ) ] ;

% fill the missing variable in the dataset
p100Level3 = fillmissing(p100Level3, 'linear','DataVariables', ...
    p100Level3.Properties.VariableNames(2:end) ) ;

% now plot the correlation between the response of A549 cell lines to
% selumetinib and PD-0325901, and also the response of MCF7 cell line to
% the same two drugs
intergrateCellLines = unique(p100ColDesc2.cell_id) ;
figure()
for ii = 1:length(intergrateCellLines)-1
    % get the sample Ids for which the cell line where treated with a
    % particular drugs
    selumetinibIDs = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'selumetinib' ) ) ;
    PD03ids = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'PD-0325901' )) ;
    
    % skip the loop if the dataset is empty
    if isempty( PD03ids) || isempty(selumetinibIDs)
        continue
    end
    
    % now get the p100 response data of the cellline to these drugs
    selResponse = p100Level3{ :, ismember(...
        p100Level3.Properties.VariableNames , selumetinibIDs) } ;
    PD03Response = p100Level3{:, ismember(...
        p100Level3.Properties.VariableNames , PD03ids) } ;
    
    % convert the response into an arrange
    selResponse = double( selResponse(:) ) ;
    PD03Response = double( PD03Response(:) ) ;
    
    % produce the plots 
    subplot(1, length(intergrateCellLines)-2 , ii)
    scatter(selResponse, PD03Response,40,'filled',...
        'MarkerFaceColor',[0.3010 0.7450 0.9330] )
    hold on
    set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
        'FontWeight','bold')
    
    % add a reference line to the plot
    X = [ones(length(selResponse),1) selResponse] ;
    b1 = X\PD03Response    ;     yCalc1 = X*b1 ;
    plot(selResponse,yCalc1,'LineWidth',1.5)
    % calculate the pearson's linear correation
    [r2 , pR ] = corr(selResponse,PD03Response, 'Type','Pearson');    
    if pR == 0
        text(0.5, 0.95, strcat( sprintf("R = %0.2f, P ", r2), ...
            convertPValue2SuperScript(pR) ) ,'FontWeight','bold', ...
            'FontSize',12 ,'Units','normalize', ...
            'HorizontalAlignment','Center')
    else
        text(0.5, 0.95, strcat( sprintf("R = %0.2f, P = ", r2), ...
            convertPValue2SuperScript(pR) ) ,'FontWeight','bold', ...
            'FontSize',12 ,'Units','normalize', ...
            'HorizontalAlignment','Center')
    end
   
    % add figure labels
    xlabel('Response to selumetinib')
    ylabel('Response to PD-0325901')
    title(sprintf('p100 Response of %s Cell Line', ...
        intergrateCellLines{ii}),'FontSize',14)   
end
hold off

% now plot the response to two cell line to the drugs using a head map
% these should exclude the treatement by DMSO
treatmentIDs = p100L3Drugs.inst_id( ...
       ismember(p100L3Drugs.pert_iname , ...
       [ unique(p100ColDesc.pert_iname) ;{'DMSO'} ] )) ;
txP100L3 = [p100Level3(:,1) , p100Level3(:, ismember(...
        p100Level3.Properties.VariableNames , treatmentIDs) ) ] ;   
    
% check thhe protein names so they are easy to visualise 
txP100L3.proteins = extractBefore( extractAfter(txP100L3.proteins ,'_'),...
    '[') ;
txP100L3.proteins = strcat( ...
    strcat(extractBefore (txP100L3.proteins, '_') ,'-'), ...
    extractBetween( txP100L3.proteins ,'_','_') ) ;

% remove the duplicates
txP100L3.proteins = strrep( ...
    matlab.lang.makeUniqueStrings(txP100L3.proteins) ,'_','-')  ;

% i need to change the names in the heatmap therefore I can those to drug
% type and cell used 
p100L3Drugs.CellLineDrug =  strrep( matlab.lang.makeUniqueStrings ( ...
    strcat(p100L3Drugs.cell_id ,cellstr(p100L3Drugs.pert_iname) ) ) , ...
    '_','-');

% get the cell line that have names and use those to change the heatmap
% naming
p100L3DrugsHeat = p100L3Drugs(p100L3Drugs.pert_iname ~= 'DMSO', :) ;
[~, locIds ] = ismember(p100L3DrugsHeat.inst_id , ...
    txP100L3.Properties.VariableNames ) ;
locIds(locIds == 0) = [] ;
txP100L3 = [ txP100L3(:,1)  , txP100L3(:, locIds) ] ;

% throw in an assertion 
assert( all( strcmp( p100L3DrugsHeat.inst_id , ...
    txP100L3.Properties.VariableNames(2:end)' ) ) ) 

figure()
hAll = heatmap(txP100L3.proteins , p100L3DrugsHeat.CellLineDrug, ...
    zscore( txP100L3{:,2:end}'),'Colormap',parula,'ColorbarVisible','on',...
    'GridVisible','on');   
title(hAll, 'p100 Response of the Cell Line') ;

clustergram( txP100L3{:,2:end}' ,'RowLabels', ...
    p100L3DrugsHeat.CellLineDrug , 'ColumnLabels', txP100L3.proteins,...
    'ColumnPDist','seuclidean','ColorMap', redgreencmap , ...
    'Linkage','complete','Standardize','column',...
    'Dendrogram', 15.6) ;

clustergram( txP100L3{:,2:end} ,'RowLabels', ...
   txP100L3.proteins , 'ColumnLabels', p100L3DrugsHeat.CellLineDrug  ,...
    'ColumnPDist','seuclidean','ColorMap', redgreencmap , ...
    'Linkage','complete','Standardize','row',...
    'Dendrogram', 15.6) ;
    
% return only the MAPK proteins and produce another heatmap 
txP100L3MAPKonly = txP100L3( contains( txP100L3.proteins  ,...
    mapkGenes.protein ,'IgnoreCase',true) , :) ;

figure()
hAll = heatmap(txP100L3MAPKonly.proteins , p100L3DrugsHeat.CellLineDrug,...
    txP100L3MAPKonly{:,2:end}','Colormap',parula, ...
    'ColorbarVisible','on','GridVisible','on');   
title(hAll, 'p100 MAPK Response of the Cell Line') ;

%%
% ************************************************************************
% SUPRISINGLY THE RESPONSE OF THESE CELL TWO CELL AT P100 LEVELS SEEM TO
% BE VERY SIMILAR UUMMMMM
% ************************************************************************

% let me check the rate of the change for protein when the cell line are
% treated with DMSO or with drugs 

% specifiy the colors for the bards 
barColors = [0.5 0.9 0.2; 0.9 0.5 0.2 ; 0.2 0.5 0.9; 0.4 0.5 0.2 ;...
    0.2 0.5 0.9; 0.4 0.5 0.2] ;
p100Level3.proteins = categorical(p100Level3.proteins);

% preallacate the mean changes values
meanChangeP100 = [];
varNames = cell(1,1);

intergrateCellLines = unique(p100ColDesc2.cell_id) ;
figure()
for ii = 1:2 % length(intergrateCellLines)
    
    % get the sample Ids for which the cell line where treated with a
    % particular drugs
    selumetinibIDs = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'selumetinib' ) ) ;
    PD03ids = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'PD-0325901' )) ;
    
    % get the response of the cell line to DMSO
    DMSOids = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'DMSO' ) ) ;
    
    % now get the p100 response data of the cellline to these drugs
    selResponse = p100Level3{ :, ismember(...
        p100Level3.Properties.VariableNames , selumetinibIDs) } ;
    PD03Response = p100Level3{:, ismember(...
        p100Level3.Properties.VariableNames , PD03ids) } ;
    DMSOResponse = p100Level3{:, ismember(...
        p100Level3.Properties.VariableNames ,  DMSOids) } ;
    
    % save the mean response of the cell lines to the two drugs 
    selResponseMean = mean(selResponse,2) - mean(DMSOResponse,2) ;
    PD03ResponseMean = mean(PD03Response,2) - mean(DMSOResponse,2) ;
    
    % specific the point of the plot
    if ii == 1
        jj = ii+1 ;
        
        % save the mean values into a table
        varNames(1) = strcat(intergrateCellLines(ii),'_selumetinib') ;
        varNames(2) = strcat(intergrateCellLines(ii),'_PD-0325901') ;
        meanChangeP100 = array2table( ...
            [selResponseMean, PD03ResponseMean],'VariableNames',varNames);  
    else
        jj = ii+2 ;
        
        % save the mean values into a table
        varNames(1) = strcat(intergrateCellLines(ii),'_selumetinib') ;
        varNames(2) = strcat(intergrateCellLines(ii),'_PD-0325901') ;
        meanChangeP100 = addvars(meanChangeP100, ...
            selResponseMean, PD03ResponseMean, ...
            'NewVariableNames',varNames);
    end
    
        % EXIT THE LOOP JUST IN CASE 
    % #############################################################
    if ii > 4 
        continue;
    end
    % #############################################################
    
    % plot a bar graph for selumetinib
    if ii == 1
        subplot(4,1,ii)
         bar(p100Level3.proteins, selResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(ii, :) )
    
%     Change the color for a particular bar by setting the FaceColor property to 'flat'. Then change the corresponding row in the CData matrix to the new RGB triplet. For example, change the color of the second bar.
%     
%     b = bar(1:10,'FaceColor','flat');
%     b.CData(2,:) = [0 0.8 0.8];

    else
        subplot(4,1,jj-1)
         bar(p100Level3.proteins, selResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(jj, :) )
    end 
   
    set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
        'FontWeight','bold','XTickLabel',[])
    ylabel('Mean Change')
    title(sprintf('%s Reponse to %s', ...
        intergrateCellLines{ii}) ,'FontSize',13)
    
    % define the groups bar names 
    protein4Plot = categorical( ...
        matlab.lang.makeUniqueStrings (extractBefore ( ...
        extractAfter ( extractBefore(cellstr(p100Level3.proteins ) , ...
        '[') ,'_') , '_')  ) ) ;
    
    % plot a bar graph for selumetinib
    subplot(4,1,jj)
    
    %  THIS KEEP FAILING SO I PUT A CONDITION WILL EDIT LATER
    
    % **********************************************************
    try 
    bar(protein4Plot, PD03ResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(ii+1, :) )
    catch
    end
    % **********************************************************
    
    
    ylabel('Mean Change') 
    
    % put xaxis labels and and Tick labels on the last plots
    if jj == 4
        xlabel('MAPK Proteins')
        set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
            'FontWeight','bold')
    else
        set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
            'FontWeight','bold','XTickLabel',[])
    end
    title( sprintf('%s Reponse to ',...
        intergrateCellLines{ii}) ,'FontSize',13)
    
    % also perform a ttest to compare the rate of change before treatment
    % and after treatment 

end

% add the proteins to the table 
meanChangeP100 = addvars( meanChangeP100, p100Level3.proteins, ...
    'Before',1 ,'NewVariableNames','proteins') ;

%% Do the plot for all the drugs

% preallacate the mean for all the drugs present in the p100 changes values
meanChangeP100All = [];
intergrateCellLines = unique(p100ColDesc2.cell_id) ;
intergratedDrugs = categories(p100L3Drugs.pert_iname) ;

% remove DMSO from the drugs and A357 from the response data because it was
% only treated with on drug
intergratedDrugs(ismember( intergratedDrugs, ...
    {'DMSO','CC-401','losmapimod','pyrazolanthrone','vemurafenib'} ...
    )) = [] ;

for ii = 1:length(intergrateCellLines)
    for jj = 1:length(intergratedDrugs )
        
        % get the sample Ids for which the cell line where treated with a
        % particular drugs
        selumetinibIDs = p100L3Drugs.inst_id( ...
            ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
            ismember(p100L3Drugs.pert_iname , intergratedDrugs{jj} ) ) ;
        
        % get the response of the cell line to DMSO
        DMSOids = p100L3Drugs.inst_id( ...
            ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
            ismember(p100L3Drugs.pert_iname , 'DMSO' ) ) ;
        
        % now get the p100 response data of the cellline to these drugs
        selResponse = p100Level3{ :, ismember(...
            p100Level3.Properties.VariableNames , selumetinibIDs) } ;
        DMSOResponse = p100Level3{:, ismember(...
            p100Level3.Properties.VariableNames ,  DMSOids) } ;
        
        % save the mean response of the cell lines to the two drugs
        selResponseMean = nanmean(selResponse,2) - nanmean(DMSOResponse,2);
            
        % specific the point of the plot
        if ii == 1
            % save the mean values into a table
            varNames = sprintf('%s_%s', ...
                intergrateCellLines{ii}, intergratedDrugs{jj} ) ;
            meanChangeP100All = array2table( ...
                selResponseMean,'VariableNames',cellstr(varNames) ) ;
        else
            % save the mean values into a table
            varNames = sprintf('%s_%s', ...
                intergrateCellLines{ii}, intergratedDrugs{jj} ) ;
            meanChangeP100All = addvars(meanChangeP100All, ...
                selResponseMean, ...
                'NewVariableNames', cellstr(varNames) );
        end
        
    end
end

% delete the columns with all NaNs that do not have the drugs
meanChangeP100All( :, all( isnan( meanChangeP100All{:,:} ) , 1) ) = [] ;

% add the proteins to the table 
meanChangeP100All = addvars( meanChangeP100All, p100Level3.proteins, ...
    'Before',1 ,'NewVariableNames','proteins') ;

% there seems to be a bug in my code therefore the combine the all mean
% change when the mean change 
meanChangeP100 = innerjoin( meanChangeP100, meanChangeP100All) ;

% check thhe protein names so they are easy to visualise 
meanChangeP100.proteins = cellstr(meanChangeP100.proteins) ;
meanChangeP100.proteins = extractBefore( extractAfter( ...
    meanChangeP100.proteins ,'_'), '[') ;
meanChangeP100.proteins = strcat( ...
    strcat(extractBefore (meanChangeP100.proteins, '_') ,'-'), ...
    extractBetween( meanChangeP100.proteins ,'_','_') ) ;

% remove the duplicates
meanChangeP100.proteins = strrep( ...
    matlab.lang.makeUniqueStrings(meanChangeP100.proteins) ,'_','-')  ;

%% Produce the heatmap of the mean change values
figure()
hAll = heatmap( meanChangeP100.proteins , ...
    strrep( meanChangeP100.Properties.VariableNames(2:end),'_','-'), ...
    meanChangeP100{:,2:end}', 'Colormap',jet, 'ColorbarVisible','on', ...
    'GridVisible','on');   
title(hAll, 'p100 MAPK Response of the Cell Line') ;

% ==================== plot only the MAPK genes now ================= 
% get the MAPK only genes
treatmentIDs = p100L3Drugs.inst_id( ...
       ismember(p100L3Drugs.pert_iname , ...
       {'selumetinib','PD-0325901','DMSO'}));
mapkP100L3 = [p100Level3(:,1) , p100Level3(:, ismember(...
        p100Level3.Properties.VariableNames , treatmentIDs) ) ] ;   
    
% convert to categorical so that I can use extract functions
mapkP100L3.proteins = cellstr( mapkP100L3.proteins) ;
    
% check thhe protein names so they are easy to visualise 
mapkP100L3.proteins = extractBefore( extractAfter( ...
    mapkP100L3.proteins ,'_'), '[') ;
mapkP100L3.proteins = strcat( ...
    strcat(extractBefore (mapkP100L3.proteins, '_') ,'-'), ...
    extractBetween( mapkP100L3.proteins ,'_','_') ) ;

% remove the duplicates
mapkP100L3.proteins = strrep( ...
    matlab.lang.makeUniqueStrings(mapkP100L3.proteins) ,'_','-')  ;

% i need to change the names in the heatmap therefore I can those to drug
% type and cell used 
p100L3Drugs.CellLineDrug =  strrep( matlab.lang.makeUniqueStrings ( ...
    strcat(cellstr(p100L3Drugs.cell_id), ...
    cellstr(p100L3Drugs.pert_iname) ) ) , ...
    '_','-');

diffP100L3MAPKonly = mapkP100L3( contains( cellstr(mapkP100L3.proteins),...
    mapkGenes.protein ,'IgnoreCase',true) , :) ;
diffP100L3MAPKonly.proteins = categorical(...
    cellstr(diffP100L3MAPKonly.proteins) );

intergrateCellLines = {'MCF7','A549'} ;
figure()
for ii = length(intergrateCellLines)
    % get the sample Ids for which the cell line where treated with a
    % particular drugs
    selumetinibIDs = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'selumetinib' ) ) ;
    PD03ids = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'PD-0325901' )) ;
    
    % get the response of the cell line to DMSO
    DMSOids = p100L3Drugs.inst_id( ...
        ismember(p100L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(p100L3Drugs.pert_iname , 'DMSO' ) ) ;
    
    % now get the p100 response data of the cellline to these drugs
    selResponse = diffP100L3MAPKonly{ :, ismember(...
        diffP100L3MAPKonly.Properties.VariableNames , selumetinibIDs) } ;
    PD03Response = diffP100L3MAPKonly{:, ismember(...
        diffP100L3MAPKonly.Properties.VariableNames , PD03ids) } ;
    DMSOResponse = diffP100L3MAPKonly{:, ismember(...
        diffP100L3MAPKonly.Properties.VariableNames ,  DMSOids) } ;
    
    % save the mean response of the cell lines to the two drugs 
    selResponseMean = mean(selResponse,2) - mean(DMSOResponse,2) ;
    PD03ResponseMean = mean(PD03Response,2) - mean(DMSOResponse,2) ;
    
    % specific the point of the plot
    if ii == 1
        jj = ii+1 ;
    else
        jj = ii+2 ;
    end 
    
    % plot a bar graph for selumetinib
    if ii == 1
        subplot(4,1,ii)
         bar(diffP100L3MAPKonly.proteins, selResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(ii, :) )
    else
        subplot(4,1,jj-1)
         bar(diffP100L3MAPKonly.proteins, selResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(jj, :) )
    end 
   
    set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
        'FontWeight','bold','XTickLabel',[])
    ylabel('Mean Change')
    title(sprintf('%s Reponse to selumetinib', ...
        intergrateCellLines{ii}) ,'FontSize',16)
    
    % plot a bar graph for selumetinib
    subplot(4,1,jj)
    bar(diffP100L3MAPKonly.proteins, PD03ResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(ii+1, :) )
    ylabel('Mean Change') 
    
    % put xaxis labels and and Tick labels on the last plots
    if jj == 4
        xlabel('MAPK Proteins')
        set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
            'FontWeight','bold')
    else
        set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
            'FontWeight','bold','XTickLabel',[])
    end
    title(sprintf('%s Reponse to selumetinib', ...
        intergrateCellLines{ii}) ,'FontSize',16)
    
    % also perform a ttest to compare the rate of change before treatment
    % and after treatment 

end

% ***********************************************************************
% What is the function of these proteins that tend to be different between
% the proteins that I have intendified
% ***********************************************************************

% let have a pca plot showing where the cells start off and then converge
% after treatment to a common theme. I think that this is what is
% happenning here since there is a very strong correlation between the
% responses of cell to drug pertubation 

% filter out the proteins that do not vary between the treatment conditions
% obtain expression measurements
ProtExpression = p100Level3{:,2:end};
genesClustering = p100Level3{:,1}; % obtain the genes

% remove nan values
nanIndices = any(isnan(ProtExpression),2);
ProtExpression(nanIndices,:) = [];
genesClustering(nanIndices) = [];
numel(genesClustering)

% Gene profiling experiments typically include genes that exhibit little
% variation in their profile and are generally not of interest. These genes
% are commonly removed from the data.

% Mask = genevarfilter(Data) calculates the variance for each gene
% expression profile in Data and returns Mask, which identifies the gene
% expression profiles with a variance less than the 10th percentile. Mask
% is a logical vector with one element for each row in Data. The elements
% of Mask corresponding to rows with a variance greater than the threshold
% have a value of 1, and those with a variance less than the threshold are
% 0.

for times = 1:4 % the number of filter times
    mask = genevarfilter(ProtExpression);
    
    ProtExpression = ProtExpression(mask,:);
    genesClustering = genesClustering(mask);
    numel(genesClustering)
    
    % filter out genes below a certain fold change threshold
    [~,ProtExpression,genesClustering] = ...
        genelowvalfilter(ProtExpression,genesClustering,'absval',log2(2));
    numel(genesClustering)
    
    % filter genes below a certain percentile: VERY POWERFUL discriminant
    [~,ProtExpression,genesClustering] = ...
        geneentropyfilter(ProtExpression,genesClustering,'prctile',20);
    numel(genesClustering)
end

% get the groups of thhe dataset from the drug informations by first
% arranging the two tables in the same order
[~, locIds ] = ismember( p100Level3.Properties.VariableNames(2:end), ...
    p100L3Drugs.inst_id );
p100L3Drugs.inst_id = p100L3Drugs.inst_id(locIds) ;

% throw in an assertion 
assert( all( strcmp( p100L3Drugs.inst_id , ...
    p100Level3.Properties.VariableNames(2:end)' ) ) ) 

% specify the groups the from the line and there treatment 
p100L3Drugs.groups = categorical( ...
    regexprep(p100L3Drugs.CellLineDrug,'-+\w*','') );

% specificy the colours
rng(18)
p100L3Drugs.pert_iname = categorical(cellstr(p100L3Drugs.pert_iname)) ;
colors1 = rand(length(categories(p100L3Drugs.pert_iname) ), 4) ;
rng(2)
p100L3Drugs.cell_id = categorical(p100L3Drugs.cell_id);
colors2 = rand(length((p100L3Drugs.cell_id) ), 4) ;

rng default % for reproducibility
yValues = tsne(ProtExpression','Algorithm','exact', ...
    'Distance','seuclidean')  ;

figure()
subplot(1,2,1)
gscatter(yValues(:,1),yValues(:,2), p100L3Drugs.pert_iname,colors1,'..',40)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('tSNE-1') ; 
ylabel('tSNE-2') ; 
title('p100 Response to Treatment','FontSize',14','FontWeight','bold')

subplot(1,2,2)
gscatter(yValues(:,1),yValues(:,2), p100L3Drugs.cell_id,colors2,'..',40)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('tSNE-1') ; 
ylabel('tSNE-2') ; 
title('p100 Response to Treatment','FontSize',14','FontWeight','bold')

% add the legend and change the marker size 
% [~, objh] = legend(lgdGroups,'location','Best',...
%     'Fontsize', 12);
% %// set font size as desired
% objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
% set(objhl, 'Markersize', 40); %// set marker size as desired

clear p100DrugsData gct_file_locDrugs theMAPKdrugs ActivitySD ...
    lgdd lincsGdscResponse X r2 pR selResponse ...
    PD03Response PD03ids selumetinibIDs intergrateDrugs ...
    intergrateCellLines PD03ResponseMean selResponseMean ...
    treatmentIDs txP100L3 p100L3DrugsHeat hAll txP100L3MAPKonly ...
    DMSOids DMSOResponse jj mapkP100L3  selResponseMean ...
    PD03ResponseMean selResponse DMSOids PD03ids diffP100L3MAPKonly ...
    treatmentIDs p100L3Drugs intergrateCellLines varNames ...
    meanChangeP100 yCalc1 xVar ii b1 ax lgdGroups ProtExpression ...
    objh yValues mask nanIndices samples times locIds colors1 colors2 ...
    geneConvert protein4Plot locUnique

%% LINCS L1000 analysis

% load the L1000 data 
% load the proteomics data from the links project
% GCT and GCTx files can be read in the same way.
% We'll use the same two files throughout this tutorial.
gct_file_L100loc = fullfile(['/Users/sinkala/Documents/MATLAB/MAPK ',...
    'Pan-cancer Analysis/LINCS proteomics'],...
    'GSE101406_Broad_LINCS_L1000_Level4_ZSPCINF_mlr12k_n1667x12328.gctx');

% LINCS proteomics/GSE101406_Broad_LINCS_P100_Level3_QCNORM_n1684x96.gctx'
% contains data of the proteomics response of 6 cancer cell lines to drug
% pertubation - I am still not sure what the columns stand for
L1000data = cmapm.Pipeline.parse_gctx(gct_file_L100loc);
disp(L1000data)

% ************************************************************************
% Here just as there is a high correlation between the p100 response of
% MCF7 to Selumetinib and PD-0325901 the MCF7 cell line are resistant to
% both of these MEK1 inhibitors. Also the A549 cell are hypersensitive to
% both Selumetinib and PD-0325901. Need to check the response p100 response
% of the MCF7 and A549 cell line
% ************************************************************************

% get the p100 response of these cell line to to PD-0325901 from the level
% 3 data lincs data
% load the proteomics data from the links project
gct_file_L100Level3 = fullfile(['/Users/sinkala/Documents/MATLAB/MAPK ',...
    'Pan-cancer Analysis/LINCS proteomics'],  ...
    'GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx');

% contains data of the proteomics response of 6 cancer cell lines to drug
% pertubation - I am still not sure what the columns stand for
l1000Level3All = cmapm.Pipeline.parse_gctx(gct_file_L100Level3);

% GCT data representation
% GCT and GCTx files are both represented in memory as structures.
disp(class(l1000Level3All));
disp(l1000Level3All);

% get the data into three table:
% 1 for the drugs, 2 for the drug description and 3, for protein
% description
l1000Level3 = addvars( array2table( l1000Level3All.mat , 'VariableNames',...
    l1000Level3All.cid ) , l1000Level3All.rid, 'Before', 1, ...
    'NewVariableNames','pr_gene_id') ;

% convert the entrez symbols to hugo genes symbol 
% get the gene IDs from ensembl using perl API or just upload the
% downloaded data and return only the protein coding genes
% system('perl webExample.pl getMyIds.xml') % perl apif0-
geneConvert = readtable('GSE101406_Broad_LINCS_L1000_gene_info.txt');
geneConvert.Properties.VariableNames(2) = "HugoSymbol" ;
l1000Level3.pr_gene_id = str2double(l1000Level3.pr_gene_id) ;

% join the two table and remove the theh entrez gene symbol and those rows
% with missing genes names 
l1000Level3 = innerjoin( geneConvert(:,[1,2]) ,l1000Level3) ;
l1000Level3 = removevars(l1000Level3,'pr_gene_id');
l1000Level3(cellfun(@isempty,l1000Level3.HugoSymbol), :) = [] ;

% load the pertabugen info from the lincs project and return only the two
% cell line that I want to check for and the two drugs plus the untreated
% DMSO to be used for baseline
l1000L3Drugs = readtable('GSE101406_Broad_LINCS_L1000_inst_info.txt');
l1000L3Drugs.pert_iname = categorical(l1000L3Drugs.pert_iname);
l1000L3Drugs = l1000L3Drugs( ...
    ismember(l1000L3Drugs.cell_id, lincsGdscCommon.cell_line) & ...
    ismember(l1000L3Drugs.pert_iname, ...
    {'selumetinib','PD-0325901','DMSO'}) , :) ;

% now only get the columns of the p100Level3 data that has the sample names
% contained in the p100L3Drugs
l1000Level3 = [l1000Level3(:,1) , l1000Level3(:, ismember(...
    l1000Level3.Properties.VariableNames, l1000L3Drugs.inst_id  ) ) ] ;

% fill the missing variable in the dataset
l1000Level3 = fillmissing(l1000Level3, 'linear','DataVariables', ...
    l1000Level3.Properties.VariableNames(2:end) ) ;

% now plot the correlation between the response of A549 cell lines to
% selumetinib and PD-0325901, and also the response of MCF7 cell line to
% the same two drugs
intergrateCellLines = {'MCF7','A549'} ;
figure()
for ii = 1:length(intergrateCellLines)
    % get the sample Ids for which the cell line where treated with a
    % particular drugs
    selumetinibIDs = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'selumetinib' ) ) ;
    PD03ids = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'PD-0325901' )) ;
    
    % now get the p100 response data of the cellline to these drugs
    selResponse = l1000Level3{ :, ismember(...
        l1000Level3.Properties.VariableNames , selumetinibIDs) } ;
    PD03Response = l1000Level3{:, ismember(...
        l1000Level3.Properties.VariableNames , PD03ids) } ;
    
    % convert the response into an arrange
    selResponse = double( selResponse(:) ) ;
    PD03Response = double( PD03Response(:) ) ;
    
    % produce the plots 
    subplot(1,2,ii)
    scatter(selResponse, PD03Response,40,'filled',...
        'MarkerFaceColor',[0.3010 0.7450 0.9330] )
    hold on
    set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
        'FontWeight','bold')
    
    % add a reference line to the plot
    X = [ones(length(selResponse),1) selResponse] ;
    b1 = X\PD03Response    ;     yCalc1 = X*b1 ;
    plot(selResponse,yCalc1,'LineWidth',1.5)
    % calculate the pearson's linear correation
    [r2 , pR ] = corr(selResponse,PD03Response, 'Type','Pearson');
    text(0.5, 0.95, strcat( sprintf("R = %0.2f, P = ", r2),...
        convertPValue2SuperScript(pR) ) ,'FontSize',12, ...
        'FontWeight','bold','Units','Normalize', ...
        'HorizontalAlignment','Center')
    
    % add figure labels
    xlabel('Response to selumetinib')
    ylabel('Response to PD-0325901')
    title(sprintf('L1000 Response of %s Cell Line', ...
        intergrateCellLines{ii}),'FontSize',14)   
end
hold off

% now plot the response to two cell line to the drugs using a head map
% these should exclude the treatement by DMSO
treatmentIDs = l1000L3Drugs.inst_id( ...
       ismember(l1000L3Drugs.pert_iname , {'selumetinib','PD-0325901'}));
txL1000L3 = [l1000Level3(:,1) , l1000Level3(:, ismember(...
        l1000Level3.Properties.VariableNames , treatmentIDs) ) ] ;   

% i need to change the names in the heatmap therefore I can those to drug
% type and cell used 
l1000L3Drugs.CellLineDrug =  strrep( matlab.lang.makeUniqueStrings ( ...
    strcat(l1000L3Drugs.cell_id ,cellstr(l1000L3Drugs.pert_iname) ) ) , ...
    '_','-');

% get the cell line that have names and use those to change the heatmap
% naming 
p100L3DrugsHeat = l1000L3Drugs(l1000L3Drugs.pert_iname ~= 'DMSO', :) ;
[~, locIds ] = ismember( txL1000L3.Properties.VariableNames(2:end) ,...
    p100L3DrugsHeat.inst_id ) ;
p100L3DrugsHeat.inst_id = p100L3DrugsHeat.inst_id(locIds) ;

% throw in an assertion 
assert( all( strcmp( p100L3DrugsHeat.inst_id , ...
    txL1000L3.Properties.VariableNames(2:end)' ) ) ) 

% return only the MAPK genes response after treatment
txL1000L3 = txL1000L3( ismember( txL1000L3.HugoSymbol ,...
    mapkGenes.protein), :) ;

figure()
hAll = heatmap(txL1000L3.HugoSymbol , p100L3DrugsHeat.CellLineDrug, ...
    txL1000L3{:,2:end}','Colormap',parula,'ColorbarVisible','on',...
    'GridVisible','on');   
title(hAll, 'L1000 Response of the Cell Line') ;

clustergram(txL1000L3{:,2:end} ,'RowLabels', ...
    txL1000L3.HugoSymbol , 'ColumnLabels', p100L3DrugsHeat.CellLineDrug,...
    'ColumnPDist','seuclidean','ColorMap', redgreencmap , ...
    'Linkage','complete','Standardize','row',...
    'Dendrogram', 15.6) ;
    
%% ==== return only only the essential genes for both cell line =====

% get the Achilles results of the two cell line and clean them up 
crisprTwo = crispr(ismember(crispr.cell_line, lincsGdscCommon.cell_line), :);
crisprTwo = rows2vars(crisprTwo,'VariableNamesSource','cell_line') ;
crisprTwo.Properties.VariableNames(1) = "HugoSymbol" ;

% return only the essential genes those with scores of less than -0.50 in
% each cell line 
crisprTwo = crisprTwo(crisprTwo.MCF7 < -0.5 | crisprTwo.A549 < -0.5 , :)  ;

txL1000L3MAPKonly = txL1000L3( contains( txL1000L3.HugoSymbol  ,...
   crisprTwo.HugoSymbol,'IgnoreCase',true) , :) ;

figure()
hAll = heatmap(txL1000L3MAPKonly.HugoSymbol  , ...
    p100L3DrugsHeat.CellLineDrug,txL1000L3MAPKonly{:,2:end}', ...
    'Colormap',parula,'ColorbarVisible','on','GridVisible','on');   
title(hAll, 'L1000 Essential MAPK Gene Response') ;

% ************************************************************************
% SUPRISINGLY THE RESPONSE OF THESE CELL TWO CELL AT P100 LEVELS SEEM TO
% BE VERY SIMILAR UUMMMMM
% ************************************************************************

% let me check the rate of the change for protein when the cell line are
% treated with DMSO or with drugs 

% specifiy the colors for the bards 
barColors = [0.5 0.9 0.2; 0.9 0.5 0.2 ; 0.2 0.5 0.9; 0.4 0.5 0.2 ; ...
   0.5 0.3 0.7; 0.3 0.8 0.6] ;
l1000Level3.HugoSymbol = categorical(l1000Level3.HugoSymbol);

% get only the mapk genes again 
txL1000L3 = l1000Level3(ismember( l1000Level3.HugoSymbol, ...
    mapkGenes.protein), :)  ;

% preallacate the mean changes values
meanChangeL1000 = [];
varNames = cell(1,1) ;

% intergrateCellLines = {'MCF7','A549','A375'} ;
intergrateCellLines = {'MCF7','A549'} ;
% figure()
for ii = 1:length(intergrateCellLines)
    % get the sample Ids for which the cell line where treated with a
    % particular drugs
    selumetinibIDs = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'selumetinib' ) ) ;
    PD03ids = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'PD-0325901' )) ;
    
    % get the response of the cell line to DMSO
    DMSOids = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'DMSO' ) ) ;
    
    % now get the p100 response data of the cellline to these drugs
    selResponse = txL1000L3{ :, ismember(...
       txL1000L3.Properties.VariableNames , selumetinibIDs) } ;
    PD03Response = txL1000L3{:, ismember(...
        txL1000L3.Properties.VariableNames , PD03ids) } ;
    DMSOResponse = txL1000L3{:, ismember(...
        txL1000L3.Properties.VariableNames ,  DMSOids) } ;
    
    % save the mean response of the cell lines to the two drugs 
    selResponseMean = mean(selResponse,2) - mean(DMSOResponse,2) ;
    PD03ResponseMean = mean(PD03Response,2) - mean(DMSOResponse,2) ;
    
    % specific the point of the plot
    if ii == 1
%         jj = ii+1 ;
        
        % save the mean values into a table
        varNames(1) = strcat(intergrateCellLines(ii),'_selumetinib') ;
        varNames(2) = strcat(intergrateCellLines(ii),'_PD-0325901') ;
        meanChangeL1000 = array2table( ...
            [selResponseMean, PD03ResponseMean],'VariableNames',varNames);  
    else
%         jj = ii+2 ;
        
        % save the mean values into a table
        varNames(1) = strcat(intergrateCellLines(ii),'_selumetinib') ;
        varNames(2) = strcat(intergrateCellLines(ii),'_PD-0325901') ;
        meanChangeL1000 = addvars(meanChangeL1000, ...
            selResponseMean, PD03ResponseMean, ...
            'NewVariableNames',varNames);
    end
    
%     % EXIT THE LOOP JUST IN CASE 
%     % #############################################################
%     if ii > 4 || jj > 4
%         continue;
%     end
%     % #############################################################
%     
%     % convert to categorical allow for plotting 
%     txL1000L3.HugoSymbol = categorical(txL1000L3.HugoSymbol)  ;
%     
%     % plot a bar graph for selumetinib
%     if ii == 1
%         subplot(4,1,ii)
%          bar(txL1000L3.HugoSymbol, selResponseMean ,...
%         'BarWidth',0.7,'FaceColor', barColors(ii, :) )
%     else
%         subplot(4,1,jj-1)
%          bar(txL1000L3.HugoSymbol , selResponseMean ,...
%         'BarWidth',0.7,'FaceColor', barColors(jj, :) )
%     end 
%    
%     set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
%         'FontWeight','bold','XTickLabel',[])
%     ylabel('Mean Change')
%     title(sprintf('%s Reponse to selumetinib', ...
%         intergrateCellLines{ii}) ,'FontSize',12)
%     
%     % plot a bar graph for selumetinib
%     subplot(4,1,jj)
%     bar(txL1000L3.HugoSymbol , PD03ResponseMean ,...
%         'BarWidth',0.7,'FaceColor', barColors(ii+1, :) )
%     ylabel('Mean Change') 
%     
%     % put xaxis labels and and Tick labels on the last plots
%     if jj == 4
%         xlabel('MAPK Proteins')
%         set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
%             'FontWeight','bold')
%     else
%         set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
%             'FontWeight','bold','XTickLabel',[])
%     end
%     title(sprintf('%s Reponse to selumetinib', ...
%         intergrateCellLines{ii}) ,'FontSize',12)
%     
%     % also perform a ttest to compare the rate of change before treatment
%     % and after treatment 

end

% add the gene names to the data
meanChangeL1000  = addvars(  meanChangeL1000 , ...
    txL1000L3.HugoSymbol ,'Before',1 , 'NewVariableNames','HugoSymbol');
 
% ************************************************************************
% get the reponse of all the cell lines 
% preallacate the mean for all the drugs present in the p100 changes values
meanChangeP100All = [];
% intergrateCellLines = {'MCF7','A549','A375'} ;
intergrateCellLines = {'MCF7','A549'} ;
intergratedDrugs = categories(l1000L3Drugs.pert_iname) ;

% remove DMSO from the drugs and A357 from the response data because it was
% only treated with on drug
intergratedDrugs(~ismember( intergratedDrugs, ...
    {'selumetinib','PD-0325901' } )) = [] ;
  
for ii = 1:length(intergrateCellLines)
    for jj = 1:length(intergratedDrugs )
        
        % get the sample Ids for which the cell line where treated with a
        % particular drugs
        selumetinibIDs = l1000L3Drugs.inst_id( ...
            ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
            ismember(l1000L3Drugs.pert_iname , intergratedDrugs{jj} ) ) ;
        
        % get the response of the cell line to DMSO
        DMSOids = l1000L3Drugs.inst_id( ...
            ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
            ismember(l1000L3Drugs.pert_iname , 'DMSO' ) ) ;
        
        % now get the L1000 response data of the cellline to these drugs
        selResponse = txL1000L3{ :, ismember(...
            txL1000L3.Properties.VariableNames , selumetinibIDs) } ;
        DMSOResponse = txL1000L3{:, ismember(...
            txL1000L3.Properties.VariableNames ,  DMSOids) } ;
        
        % save the mean response of the cell lines to the two drugs
        selResponseMean = nanmean(selResponse,2) - nanmean(DMSOResponse,2);
            
        % specific the point of the plot
        if ii == 1
            % save the mean values into a table
            varNames = sprintf('%s_%s', ...
                intergrateCellLines{ii}, intergratedDrugs{jj} ) ;
            meanChangeL1000All = array2table( ...
                selResponseMean,'VariableNames',cellstr(varNames) ) ;
        else
            % save the mean values into a table
            varNames = sprintf('%s_%s', ...
                intergrateCellLines{ii}, intergratedDrugs{jj} ) ;
            meanChangeL1000All = addvars(meanChangeL1000All, ...
                selResponseMean, ...
                'NewVariableNames', cellstr(varNames) );
        end
        
    end
end

% delete the columns with all NaNs that do not have the drugs
meanChangeL1000All( :, all( isnan( meanChangeL1000All{:,:} ) , 1) ) = [] ;

% add the proteins to the table 
meanChangeL1000All= addvars( meanChangeL1000All, txL1000L3.HugoSymbol, ...
    'Before',1 , 'NewVariableNames','HugoSymbol');

% there seems to be a bug in my code therefore the combine the all mean
% change when the mean change 
meanChangeL1000 = innerjoin( meanChangeL1000, meanChangeL1000All) ;

% remove the duplicates
meanChangeL1000.HugoSymbol = strrep( ...
    matlab.lang.makeUniqueStrings( ...
    cellstr( meanChangeL1000.HugoSymbol) ),'_','-')  ;

% ************************************************************************

% Produce the heatmap of the mean change values
figure()
hAll = heatmap(meanChangeL1000.HugoSymbol  , ...
    strrep( meanChangeL1000.Properties.VariableNames(2:end), ... 
    '_','-') , zscore( meanChangeL1000{:,2:end} )', ...
    'Colormap',jet, 'ColorbarVisible','on','GridVisible','on');   
title(hAll, 'L1000 MAPK Response of the Cell Line') ;

% ==================== plot only the MAPK genes now ================= 
% get the MAPK only genes
treatmentIDs = l1000L3Drugs.inst_id( ...
       ismember(l1000L3Drugs.pert_iname , ...
       {'selumetinib','PD-0325901','DMSO'}));
mapkL1000L3 = [txL1000L3(:,1) , txL1000L3(:, ismember(...
        txL1000L3.Properties.VariableNames , treatmentIDs) ) ] ;   
    
% convert to categorical so that I can use extract functions
mapkL1000L3.HugoSymbol = cellstr( mapkL1000L3.HugoSymbol ) ;

% i need to change the names in the heatmap therefore I can those to drug
% type and cell used 
l1000L3Drugs.CellLineDrug =  strrep( matlab.lang.makeUniqueStrings ( ...
    strcat(l1000L3Drugs.cell_id ,cellstr(l1000L3Drugs.pert_iname) ) ) , ...
    '_','-');

diffL1000L3MAPKonly = mapkL1000L3( contains(...
    cellstr(mapkL1000L3.HugoSymbol),mapkGenes.protein ,'IgnoreCase',true),:);
diffL1000L3MAPKonly.proteins = categorical(...
    cellstr(diffL1000L3MAPKonly.HugoSymbol) );

% intergrateCellLines = {'MCF7','A549'} ;
figure()
for ii = 1:length(intergrateCellLines)
    % get the sample Ids for which the cell line where treated with a
    % particular drugs
    selumetinibIDs = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'selumetinib' ) ) ;
    PD03ids = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'PD-0325901' )) ;
    
    % get the response of the cell line to DMSO
    DMSOids = l1000L3Drugs.inst_id( ...
        ismember(l1000L3Drugs.cell_id, intergrateCellLines(ii) ) & ...
        ismember(l1000L3Drugs.pert_iname , 'DMSO' ) ) ;
    
    % now get the p100 response data of the cellline to these drugs
    selResponse = diffL1000L3MAPKonly{ :, ismember(...
        diffL1000L3MAPKonly.Properties.VariableNames , selumetinibIDs) } ;
    PD03Response = diffL1000L3MAPKonly{:, ismember(...
        diffL1000L3MAPKonly.Properties.VariableNames , PD03ids) } ;
    DMSOResponse = diffL1000L3MAPKonly{:, ismember(...
        diffL1000L3MAPKonly.Properties.VariableNames ,  DMSOids) } ;
    
    % save the mean response of the cell lines to the two drugs 
    selResponseMean = mean(selResponse,2) - mean(DMSOResponse,2) ;
    PD03ResponseMean = mean(PD03Response,2) - mean(DMSOResponse,2) ;
    
    % convert to categorical to allow for plotting
    diffL1000L3MAPKonly.HugoSymbol = categorical( ...
        diffL1000L3MAPKonly.HugoSymbol) ;
    
    % specific the point of the plot
    if ii == 1
        jj = ii+1 ;
    else
        jj = ii+2 ;
    end 
    
    % skip the loop when the value is greater than 5
%     if jj > 5 
%         continue
%     end
    
    % plot a bar graph for selumetinib
    if ii == 1
        subplot(6,1,ii)
         bar(diffL1000L3MAPKonly.HugoSymbol , selResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(ii, :) )
    else
        subplot(4,1,jj-1)
         bar(diffL1000L3MAPKonly.HugoSymbol , selResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(jj, :) )
    end 
   
    set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
        'FontWeight','bold','XTickLabel',[])
    ylabel('Mean Change')
    title(sprintf('%s Reponse to selumetinib', ...
        intergrateCellLines{ii}) ,'FontSize',16)
    
    % plot a bar graph for selumetinib
    subplot(6,1,jj)
    bar(diffL1000L3MAPKonly.HugoSymbol , PD03ResponseMean ,...
        'BarWidth',0.7,'FaceColor', barColors(ii+1, :) )
    ylabel('Mean Change') 
    
    % put xaxis labels and and Tick labels on the last plots
    if jj == 4
        xlabel('MAPK Proteins')
        set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
            'FontWeight','bold')
    else
        set(gca,'FontSize',12,'LineWidth',1,'FontSize',10 ,'Box','off',...
            'FontWeight','bold','XTickLabel',[])
    end
    title(sprintf('%s Reponse to selumetinib', ...
        intergrateCellLines{ii}) ,'FontSize',16)
    
    % also perform a ttest to compare the rate of change before treatment
    % and after treatment 

end

% ***********************************************************************
% What is the function of these proteins that tend to be different between
% the proteins that I have intendified
% ***********************************************************************

% let have a pca plot showing where the cells start off and then converge
% after treatment to a common theme. I think that this is what is
% happenning here since there is a very strong correlation between the
% responses of cell to drug pertubation 

% filter out the proteins that do not vary between the treatment conditions
% obtain expression measurements
ProtExpression = l1000Level3{:,2:end};
genesClustering = l1000Level3{:,1}; % obtain the genes

% remove nan values
nanIndices = any(isnan(ProtExpression),2);
ProtExpression(nanIndices,:) = [];
genesClustering(nanIndices) = [];
numel(genesClustering)

% Gene profiling experiments typically include genes that exhibit little
% variation in their profile and are generally not of interest. These genes
% are commonly removed from the data.

% Mask = genevarfilter(Data) calculates the variance for each gene
% expression profile in Data and returns Mask, which identifies the gene
% expression profiles with a variance less than the 10th percentile. Mask
% is a logical vector with one element for each row in Data. The elements
% of Mask corresponding to rows with a variance greater than the threshold
% have a value of 1, and those with a variance less than the threshold are
% 0.

for times = 1:10 % the number of filter times
    mask = genevarfilter(ProtExpression);
    
    ProtExpression = ProtExpression(mask,:);
    genesClustering = genesClustering(mask);
    numel(genesClustering)
    
    % filter out genes below a certain fold change threshold
    [~,ProtExpression,genesClustering] = ...
        genelowvalfilter(ProtExpression,genesClustering,'absval',log2(3));
    numel(genesClustering)
    
    % filter genes below a certain percentile: VERY POWERFUL discriminant
    [~,ProtExpression,genesClustering] = ...
        geneentropyfilter(ProtExpression,genesClustering,'prctile',20);
    numel(genesClustering)
end

% get the groups of thhe dataset from the drug informations by first
% arranging the two tables in the same order
[~, locIds ] = ismember( l1000Level3.Properties.VariableNames(2:end), ...
    l1000L3Drugs.inst_id );
l1000L3Drugs.inst_id = l1000L3Drugs.inst_id(locIds) ;

% throw in an assertion 
assert( all( strcmp( l1000L3Drugs.inst_id , ...
    l1000Level3.Properties.VariableNames(2:end)' ) ) ) 

% specify the groups the from the line and there treatment 
l1000L3Drugs.groups = categorical( ...
    regexprep(l1000L3Drugs.CellLineDrug,'-+\w*','') );

% specificy the colours
rng(18)
l1000L3Drugs.pert_iname = categorical(cellstr(l1000L3Drugs.pert_iname)) ;
colors1 = rand(length(categories(l1000L3Drugs.pert_iname) ), 3) ;
rng(5)
l1000L3Drugs.cell_id = categorical(l1000L3Drugs.cell_id);
colors2 = rand(length((l1000L3Drugs.cell_id) ), 3) ;

rng default % for reproducibility
yValues = tsne(ProtExpression','Algorithm','exact', ...
    'Distance','correlation')  ;

figure()
subplot(1,2,1)
gscatter(yValues(:,1),yValues(:,2), l1000L3Drugs.pert_iname,colors1,'..',40)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('tSNE-1') ; 
ylabel('tSNE-2') ; 
title('L1000 Response to Treatment','FontSize',14','FontWeight','bold')

subplot(1,2,2)
gscatter(yValues(:,1),yValues(:,2), l1000L3Drugs.cell_id,colors2,'..',40)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('tSNE-1') ; 
ylabel('tSNE-2') ; 
title('p100 Response to Treatment','FontSize',14','FontWeight','bold')

clear p100DrugsData gct_file_locDrugs theMAPKdrugs ActivitySD ...
    lgdd lincsGdscResponse X r2 pR selResponse ...
    PD03Response PD03ids selumetinibIDs intergrateDrugs ...
    intergrateCellLines PD03ResponseMean selResponseMean ...
    treatmentIDs txP100L3 p100L3DrugsHeat hAll txP100L3MAPKonly ...
    DMSOids DMSOResponse jj mapkP100L3  selResponseMean ...
    PD03ResponseMean selResponse DMSOids PD03ids diffP100L3MAPKonly ...
    treatmentIDs p100L3Drugs intergrateCellLines varNames ...
    meanChangeP100 yCalc1 xVar ii b1 ax lgdGroups ProtExpression ...
    objh yValues mask nanIndices samples times locIds colors1 colors2 ...
    geneConvert curDose curColor plotColors PosTargets ...
    pvalue theMissing txL1000L3 txL1000L3MAPKonly ...
    gct_file_location gct_file_L100loc gct_file_L1000Level3  ...
    curResponse colors checkGene cancerCheck barColors ans 


%% Check the dependence of KRAS mutants 

% From Project Drive: In lung cancer, we observed a substantial number of
% KRAS mutant cell lines in which KRAS is dispensable for growth. Instead,
% these KRAS mutant lung lines are susceptible to NFE2L2 (NRF2) and SMARCA2
% (BRM) knock- down correlated with loss-of-function mutations in KEAP1
% and/or low SMARCA4 (BRG1), respectively (Figure 2F). This has important
% implications for treating KRAS mutant cancers with MAPK pathway
% inhibitors, as co-occurring mutation is predicted to lead to de novo
% resistance.

% ************************************************************************
% Does this mean that when the KRAS mutant cell line are treatment with
% MAPK inhibitors, the levels of  NFE2L2 (NRF2) and SMARCA2 (BRM) go down? 
% ************************************************************************

%% MAPK Pathway alteration across tumours of different population

% use ICGC data to show that alteration within the MAPK are the same or may
% differ across tumour types based on the patients demographics

% for liver cancer

% for breast cancer


%% What drugs can we use to treat cancers with Lower MAPK dependence scores?


%% Internal Functions

function addReferenceLineToPlot(xVariable, yVariable)
    % This function add a reference line to the plot and the R square value

    % preform a regression calculation and add it to the plot
    X = [ones(length(xVariable),1) xVariable] ;
    b1 = X\yVariable    ;     yCalc1 = X*b1 ;
    plot(xVariable,yCalc1,'LineWidth',1.5,'Color','k')
    
    % calculate the pearson's linear correation  and add it to the plot
    [r2 , pR ] = corr(xVariable,yVariable, 'Type','Pearson',...
        'Rows','complete');
    if pR < 0.0001
        text(0.4, 0.9 , strcat( sprintf( ...
            "R = %0.2f, P = ", r2), convertPValue2SuperScript(pR)  ), ...
            'Units','normalized','FontSize',11,'FontWeight','bold')
    else
        text(0.4, 0.9 ,sprintf('R = %0.2f, P = %0.4f', r2, pR), ...
            'Units','normalized','FontSize',11,'FontWeight','bold')
    end
end


% *********************** end of function ************************
% ****************************************************************



function [pValue , stats, stats2] = getMutationExprCrisprAssociation( ...
    checkGene, crisprData, ...
    sampleInfoData , mutationData , mrnaData , checkCancer)

    % Show that pancreatic cancer cell line are highly dependent on KRAS
    % which is mutated in most tumours
    geneCrispr = crisprData(:,{'disease','cell_line',checkGene} )  ;
    geneExpr = innerjoin( sampleInfoData(:, {'disease','cell_line'}) , ...
        mrnaData( :,{'cell_line',checkGene} ) ) ;

    % also add the mutation data to the expression data
    geneMut = mutationData(:,{'cell_line',checkGene} )  ;
    geneMut.Properties.VariableNames{2} = sprintf('%smut',checkGene) ;

    % add that to the kras mRNA expression data
    geneExpr = innerjoin(geneExpr, geneMut) ;

    % make sure that the two dataset have the same cell line and arrange
    % them in the same order between the mRNA expression data and the
    % CRISPR data
    [~, iA, iB] = intersect(geneCrispr.cell_line, geneExpr.cell_line) ;
    geneCrispr = geneCrispr(iA, :) ;
    geneExpr = geneExpr(iB, :) ;

    % throw in an assertion
    assert( all( strcmp( geneCrispr.cell_line, geneExpr.cell_line) ) )

    % set the plot colors between pancreatic cancer and other cancers and
    % also between those with mutations in the KRAS genes and those with no
    % mutations in the kras genes

    % set the colors
    colors1 = [0.392 0.831 0.074 ; 0.1 0 1] ;

    % ================= produce the scatter plots ==================
    figure()
    clf ;
    axes('position',[0.15, 0.33, 0.30, 0.60]);
    gscatter( geneCrispr.(checkGene) , geneExpr.(checkGene), ...
        ismember(geneCrispr.disease, checkCancer) , colors1, '..',20)
    set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
    xlabel( sprintf('%s Dependence', checkGene) ) ;
    ylabel( sprintf('%s mRNA Expression',checkGene) ) ;
    title(sprintf('%s : %s vs Other CLL' ,checkGene, checkCancer),...
        'FontSize',14','FontWeight','bold')
    hold on
    % add a reference line and the correlation coefficent
    addReferenceLineToPlot(geneCrispr.(checkGene), geneExpr.(checkGene) )
    legend({'Others CCL',sprintf('%s CLL',checkCancer)} ,...
        'Location','Best')

    % ============ add box plot at the bottom of the chart =============

    % get the limit of the scatter plot to use those as the limit of the
    % bar graph
    theXlim = get(gca,'XLim') ;

    % plot the data
    axes('position',[0.15, 0.1, 0.30, 0.15]);
    boxplot( ...
        geneCrispr.(checkGene) , ...
        ismember(geneCrispr.disease, checkCancer),...
        'Orientation','horizontal','Width',0.7 ,'Color', colors1 ,...
        'Symbol','r+' ) ;

    % set some figure properties and add title ot the figure
    set(gca,'Box','off','YColor',[1 1 1], 'XColor',[1 1 1] ,'XLim',theXlim)

    % set the line width of the box plots
    set(findobj(gca,'type','line'),'linew',2)

    % add color to the box plots
    colors2 = flipud(colors1) ;
    h4 = findobj(gca,'Tag','Box') ;
    for jj=1:length(h4)
        patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
            colors2(jj,:),'FaceAlpha', .8 ,'LineWidth',1);
    end

    % perform a ttest and add the p value to the plot
    [~, pValue,~, stats] =  ttest2( geneCrispr.(checkGene) , ...
        ismember(geneCrispr.disease, checkCancer),...
        'Vartype','unequal') ;
    
    % put the p value on the plot depending on its value
    if pValue < 0.0001
    text( 0.35, 0.1 , convertPValue2SuperScript(pValue),...
        'Units','normalized','FontWeight','bold','FontSize',12)
    else
       text( 0.35, 0.1 ,sprintf('p = %0.4f',pValue),...
           'Units','normalized','FontWeight','bold','FontSize',12) 
    end
    hold off

    % ******************************************************************
    % ======================== subplot 2 =============================
    axes('position',[0.50, 0.33, 0.30, 0.60]);
    
    % specificy the groups for the mutations data and the copy number
    % change
    if iscell( geneExpr.(sprintf('%smut',checkGene)) ) 
        alterationsGroups = ...
            ~cellfun( @isempty ,geneExpr.(sprintf('%smut',checkGene)) ) ;
        
        weHave = 'mutations' ;
        
        % specifiy the colors for mutations
        colors1 = [1 0.43 0.1649 ; 0.0745 0.623 1] ;
         
    elseif isnumeric( geneExpr.(sprintf('%smut',checkGene)) ) 
        alterationsGroups = geneExpr.(sprintf('%smut',checkGene))  ;
        
        weHave = 'CNA' ;
        
        % specify the color for CNA
         colors1 = [1 0.43 0.1649 ; 0.0745 0.623 1 ; 0.49,0.18,0.56] ;
    else
        error('Check that the mutations or CNA data correctly formatted')
    end
    
    % now plot the data 
    gscatter( geneCrispr.(checkGene) , geneExpr.(checkGene) , ...
         alterationsGroups  ,colors1, '..',20) ;
     
    set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
    xlabel( sprintf('%s Dependence', checkGene) ) ;
    ylabel( sprintf('%s mRNA Expression',checkGene) ) ;
    
    % have different titles for mutations for copy number data
    if strcmp( weHave , 'mutations')
        title(sprintf('Cell Fitness: %s Mutation Status',checkGene),...
            'FontSize',14','FontWeight','bold')
        hold on
        % add a reference line and the correlation coefficent
        addReferenceLineToPlot(geneCrispr.(checkGene),...
            geneExpr.(checkGene) )
        legend({'not mutated',sprintf('%s mutated',checkGene )} , ...
            'Location','NorthEast')
        hold off
    elseif strcmp( weHave , 'CNA')
        title(sprintf('Cell Fitness: %s CNV Status',checkGene),...
            'FontSize',14','FontWeight','bold')
        hold on
        % add a reference line and the correlation coefficent
        addReferenceLineToPlot(geneCrispr.(checkGene),...
            geneExpr.(checkGene) )
        legend({ sprintf('%s Deleted',checkGene ) ,'not altered', ...
             sprintf('%s Amplified',checkGene ) } , ...
            'Location','NorthEast')
        hold off
    end


    % ============ add box plot at the bottom of the chart =============

    % get the limit of the scatter plot to use those as the limit of the bar
    % graph
    theXlim = get(gca,'XLim') ;

    % plot the data
    axes('position',[0.50, 0.1, 0.30, 0.15]);
    boxplot( ...
        geneCrispr.(checkGene) , alterationsGroups  ,...
        'Orientation','horizontal','Width',0.7 ,'Color', colors1 ,...
        'Symbol','r+') ;

    % set some figure properties and add title ot the figure
    set(gca,'Box','off','YColor',[1 1 1], 'XColor',[1 1 1] ,'XLim',theXlim)

    % set the line width of the box plots
    set(findobj(gca,'type','line'),'linew',2)

    % add color to the box plots
    colors2 = flipud(colors1) ;
    h4 = findobj(gca,'Tag','Box') ;
    for jj=1:length(h4)
        patch(get(h4(jj),'XData'),get(h4(jj),'YData'),...
            colors2(jj,:),'FaceAlpha', .8 ,'LineWidth',1);
    end

    % perform a ttest and add the p value to the plot
    % have different titles for mutations for copy number data
    if strcmp( weHave , 'mutations')
        [~, pValue,~, stats2 ] =  ttest2(geneCrispr.(checkGene) , ...
            alterationsGroups  ,'Vartype','unequal') ;
    elseif strcmp( weHave , 'CNA')
         pValue =  anova1(geneCrispr.(checkGene), alterationsGroups,'off');
    end
     
    
    % put the p values
    if pValue < 0.0001
        text( 0.35, 0.1 , convertPValue2SuperScript(pValue),...
            'Units','normalized','FontWeight','bold','FontSize',12)
        hold off
    else
        text( 0.35, 0.1 ,sprintf('p = %0.4f',pValue) ,...
            'Units','normalized','FontWeight','bold','FontSize',12)
        hold off   
    end

end


% *********************** end of function ************************
% ****************************************************************

% ====================== another function =================



function createLegendInternal(yPoint, xStart, legendLabels , plotColors,...
    myLgdTitle , fontSizes ,rectAndTextBox)

    % specificy the y values starts and mode of actions for the drugs
    % yPoint = 0.830 ; xStart = 0.1023 ;
    xStartText = xStart + 0.01 ;
    yPointTitle = yPoint + 0.03 ;

    % specify the font size to be used in the plot
    if ~exist('fontSizes','var')
        fontSizes = [10, 12] ;
    end

    % specifiy the rectangle and text length
    if ~exist('rectAndTextBox','var')
        rectAndTextBox = [0.018 ,0.12] ;
    end

    % check for errors
    if ~isnumeric(yPoint) || ~isnumeric(xStart)
        error('Both yPoint and xStarts should be numeric values')
    elseif yPoint > 1 || xStart > 1
        error('Both yPoint and xStarts should be less than 1')
    elseif ~isnumeric(plotColors)
        error('plot Color should be numeric')
    end

    if size(plotColors,1) ~= length(legendLabels)
        error('There should be a color for each legend names')
    end

    if iscategorical( legendLabels)
        legendLabels = categories(legendLabels);
    end

    for ii = 1:length(legendLabels)
        % add the legend color
        annotation('rectangle',[xStart yPoint rectAndTextBox(1) 0.023],...
            'EdgeColor', plotColors(ii,:), ...
            'FaceColor', plotColors(ii,:));

        % add the legend text
        annotation('textbox',[xStartText yPoint rectAndTextBox(2) 0.0230],...
            'String',legendLabels{ii},'FontSize',fontSizes(1),...
            'FontName','Helvetica Neue','FitBoxToText','off',...
            'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
            'VerticalAlignment','middle','FontWeight','normal')

        % move the y point down
        yPoint = yPoint - 0.03 ;
    end

    % add the title
    annotation('textbox',[xStart yPointTitle rectAndTextBox(2) 0.0230],...
        'String', myLgdTitle,'FontSize',fontSizes(2),...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
        'VerticalAlignment','middle','FontWeight','bold',...
        'HorizontalAlignment','left');

end


% *********************** end of function ************************
% ****************************************************************

%  ================== another function  ====================


function c = redblue(m)
%   REDBLUE Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG,
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b];

end


% *********************** end of function ************************
% ****************************************************************

%  ================== another function  ====================


function varargout=sigstar(groups,stats,nosort)
    % sigstar - Add significance stars to bar charts, boxplots, line charts, etc,
    %
    % H = sigstar(groups,stats,nsort)
    %
    % Purpose
    % Add stars and lines highlighting significant differences between pairs of groups. 
    % The user specifies the groups and associated p-values. The function handles much of 
    % the placement and drawing of the highlighting. Stars are drawn according to:
    %   * represents p<=0.05
    %  ** represents p<=1E-2
    % *** represents p<=1E-3
    %
    %
    % Inputs
    % groups - a cell array defining the pairs of groups to compare. Groups defined 
    %          either as pairs of scalars indicating locations along the X axis or as 
    %          strings corresponding to X-tick labels. Groups can be a mixture of both 
    %          definition types.
    % stats -  a vector of p-values the same length as groups. If empty or missing it's 
    %          assumed to be a vector of 0.05s the same length as groups. Nans are treated
    %          as indicating non-significance.
    % nsort -  optional, 0 by default. If 1, then significance markers are plotted in 
    %          the order found in groups. If 0, then they're sorted by the length of the 
    %          bar.
    %
    % Outputs
    % H - optionally return handles for significance highlights. Each row is a different
    %     highlight bar. The first column is the line. The second column is the text (stars).
    %     
    %
    % Examples
    % 1. 
    % bar([5,2,1.5])
    % sigstar({[1,2], [1,3]})
    %
    % 2. 
    % bar([5,2,1.5])
    % sigstar({[2,3],[1,2], [1,3]},[nan,0.05,0.05])
    %
    % 3.  **DOESN'T WORK IN 2014b**
    % R=randn(30,2);
    % R(:,1)=R(:,1)+3;
    % boxplot(R)
    % set(gca,'XTick',1:2,'XTickLabel',{'A','B'})
    % H=sigstar({{'A','B'}},0.01);
    % ylim([-3,6.5])
    % set(H,'color','r')
    %
    % 4. Note the difference in the order with which we define the groups in the 
    %    following two cases. 
    % x=[1,2,3,2,1];
    % subplot(1,2,1)
    % bar(x)
    % sigstar({[1,2], [2,3], [4,5]})
    % subplot(1,2,2)
    % bar(x)
    % sigstar({[2,3],[1,2], [4,5]})
    %
    % ALSO SEE: demo_sigstar
    %
    % KNOWN ISSUES:
    % 1. Algorithm for identifying whether significance bar will overlap with 
    %    existing plot elements may not work in some cases (see line 277)
    % 2. Bars may not look good on exported graphics with small page sizes.
    %    Simply increasing the width and height of the graph with the 
    %    PaperPosition property of the current figure should fix things.
    %
    % Rob Campbell - CSHL 2013

    %Input argument error checking

    %If the user entered just one group pair and forgot to wrap it in a cell array 
    %then we'll go easy on them and wrap it here rather then generate an error
    if ~iscell(groups) & length(groups)==2
        groups={groups};
    end

    if nargin<2 
        stats=repmat(0.05,1,length(groups));
    end
    if isempty(stats)
        stats=repmat(0.05,1,length(groups));
    end
    if nargin<3
        nosort=0;
    end


    %Check the inputs are of the right sort
    if ~iscell(groups)
        error('groups must be a cell array')
    end

    if ~isvector(stats)
        error('stats must be a vector')
    end

    if length(stats)~=length(groups)
        error('groups and stats must be the same length')
    end

    %Each member of the cell array groups may be one of three things:
    %1. A pair of indices.
    %2. A pair of strings (in cell array) referring to X-Tick labels
    %3. A cell array containing one index and one string
    %
    % For our function to run, we will need to convert all of these into pairs of
    % indices. Here we loop through groups and do this. 

    xlocs=nan(length(groups),2); %matrix that will store the indices 
    xtl=get(gca,'XTickLabel');  

    for ii=1:length(groups)
        grp=groups{ii};

        if isnumeric(grp)
            xlocs(ii,:)=grp; %Just store the indices if they're the right format already

        elseif iscell(grp) %Handle string pairs or string/index pairs

            if isstr(grp{1})
                a=strmatch(grp{1},xtl);
            elseif isnumeric(grp{1})
                a=grp{1};
            end
            if isstr(grp{2})
                b=strmatch(grp{2},xtl);
            elseif isnumeric(grp{2})
                b=grp{2};
            end

            xlocs(ii,:)=[a,b];
        end

        %Ensure that the first column is always smaller number than the second
        xlocs(ii,:)=sort(xlocs(ii,:));

    end

    %If there are any NaNs we have messed up. 
    if any(isnan(xlocs(:)))
        error('Some groups were not found')
    end

    %Optionally sort sig bars from shortest to longest so we plot the shorter ones first
    %in the loop below. Usually this will result in the neatest plot. If we waned to 
    %optimise the order the sig bars are plotted to produce the neatest plot, then this 
    %is where we'd do it. Not really worth the effort, though, as few plots are complicated
    %enough to need this and the user can define the order very easily at the command line. 
    if ~nosort
        [~,ind]=sort(xlocs(:,2)-xlocs(:,1),'ascend');
        xlocs=xlocs(ind,:);groups=groups(ind);
        stats=stats(ind);
    end



    %-----------------------------------------------------
    %Add the sig bar lines and asterisks 
    holdstate=ishold;
    hold on

    H=ones(length(groups),2); %The handles will be stored here

    y=ylim;
    yd=myRange(y)*0.05; %separate sig bars vertically by 5% 

    for ii=1:length(groups)
        thisY=findMinY(xlocs(ii,:))+yd;
        H(ii,:)=makeSignificanceBar(xlocs(ii,:),thisY,stats(ii));
    end
    %-----------------------------------------------------




    %Now we can add the little downward ticks on the ends of each line. We are
    %being extra cautious and leaving this it to the end just in case the y limits
    %of the graph have changed as we add the highlights. The ticks are set as a
    %proportion of the y axis range and we want them all to be the same the same
    %for all bars.
    yd=myRange(ylim)*0.01; %Ticks are 1% of the y axis range
    for ii=1:length(groups)
        y=get(H(ii,1),'YData');
        y(1)=y(1)-yd;
        y(4)=y(4)-yd;   
        set(H(ii,1),'YData',y)
    end




    %Be neat and return hold state to whatever it was before we started
    if ~holdstate
        hold off
    elseif holdstate
        hold on
    end


    %Optionally return the handles to the plotted significance bars (first column of H)
    %and asterisks (second column of H).
    if nargout>0
        varargout{1}=H;
    end


end %close sigstar



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Internal functions

function H=makeSignificanceBar(x,y,p)
    %makeSignificanceBar produces the bar and defines how many asterisks we get for a 
    %given p-value


    if p<=1E-3
        stars='***'; 
    elseif p<=1E-2
        stars='**';
    elseif p<=0.05
        stars='*';
    elseif isnan(p)
        stars='n.s.';
    else
        stars='';
    end
            
    x=repmat(x,2,1);
    y=repmat(y,4,1);

    H(1)=plot(x(:),y,'-k','LineWidth',1.0,'Tag','sigstar_bar');

    %Increase offset between line and text if we will print "n.s."
    %instead of a star. 
    if ~isnan(p)
        offset=0.005;
    else
        offset=0.02;
    end

    starY=mean(y)+myRange(ylim)*offset;
    H(2)=text(mean(x(:)),starY,stars,...
        'HorizontalAlignment','Center',...
        'BackGroundColor','none',...
        'Tag','sigstar_stars','FontSize',20);

    Y=ylim;
    if Y(2)<starY
        ylim([Y(1),starY+myRange(Y)*0.05])
    end


end %close makeSignificanceBar



function Y=findMinY(x)
    % The significance bar needs to be plotted a reasonable distance above all the data points
    % found over a particular range of X values. So we need to find these data and calculat the 
    % the minimum y value needed to clear all the plotted data present over this given range of 
    % x values. 
    %
    % This version of the function is a fix from Evan Remington
    oldXLim = get(gca,'XLim');
    oldYLim = get(gca,'YLim');

    axis(gca,'tight')
    set(gca,'xlim',x) %Matlab automatically re-tightens y-axis

    yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
    Y = max(yLim);

    axis(gca,'normal')
    set(gca,'XLim',oldXLim,'YLim',oldYLim)

end %close findMinY


function rng=myRange(x)
    %replacement for stats toolbox range function
    rng = max(x) - min(x);
end %close myRange


% *********************** end of function ************************
% ****************************************************************

% ======================= another function =========================
function pSuperScript = convertPValue2SuperScript(p)
    % converts the number to scientific superscript for printing on a
    % figure
    pS = num2str(p) ;

    % get the first number
    firstNumbers = extractBefore(pS,'e') ;

    % check if there is a decimal place. then only get the first 4 numbers
    if contains( firstNumbers  ,'.')
        firstNumbers = firstNumbers(1:4) ;
    end
    
    % get the correctly formated p value
    pSuperScript = sprintf('%s x 10^{%d}', firstNumbers, ...
        str2double(extractAfter(pS, 'e') )) ;
    
    % if the p value is large
    if p > 0.0001
       pSuperScript = sprintf('%0.4f', p) ;
    elseif p == 0
         pSuperScript = sprintf('< 1 x 10^{%d}', -300) ;
    end

end

% *********************** end of function ************************
% ****************************************************************

% ======================= another function =========================


% This function processes TCGA data for pancancer studies into single files
% for each dataset for all studies


function [mutations, cancerStudies, cnaData, mrna, clinicalData ] = ...
    getcBioPortalDataAllStudies(myGenes,mutBias)

% the api seem to be non functional at the moment. So I use a text filed
% that i download from cBioportal to run my analysis

try
    % Get Data from cBioPortal
    % Set web API URL (excluding 'webservice.do', trailing slash optional)
    cgdsURL = 'http://www.cbioportal.org/';
    % Get list of available cancer types
    cancerStudies = getcancerstudies(cgdsURL);
    cancerStudies = struct2table(cancerStudies);
catch
    cancerStudies = readtable('cancerStudies.txt');
    cancerStudies.Properties.VariableNames(1) = {'cancerTypeId'};
end

% Get the list cancer studies codes & return only Cancer studies that have
% both mutations and copy number alterations data
toKeepStudies = false(height(cancerStudies),1) ;
for kk = 1:height(cancerStudies)
    fprintf('\n Checking Genetic Profiles for %s \n',...
    cancerStudies.name{kk} )
     
    % check if the study is a TCGA provisional then do away with it
    if contains(cancerStudies.name{kk},'Provisional', ...
            'IgnoreCase',true) 
        fprintf('\n This is a Provisional Study - Excluded \n')
        fprintf('\n ======================================================= \n\n')
        continue
    end
    
    % remove Non pan cancer TCGA studies because they form duplicate
    % studies in the data
    % check if the study is a TCGA provisional then do away with it
    if contains(cancerStudies.name{kk},'TCGA', 'IgnoreCase',true)
        if ~contains(cancerStudies.name{kk},'PanCancer Atlas', ...
                'IgnoreCase',true)
            fprintf('\n This is a Non Pan-Cancer Study - Excluded \n')
            fprintf(['\n =======================================',...
            '================ \n\n'])
            continue
        end
    end
    
    % remove CCLE studies because they are cell line data this removes
    if contains(cancerStudies.name{kk},...
            {'Cell Line','NCI-60','xenograft','Pediatric','Histio',...
            'MIXED','Summit'},'IgnoreCase',true )
        fprintf('\n These are CCLE or NCI data - Excluded \n')
        fprintf(['\n =======================================',...
            '================ \n\n'])
        continue
    end

%     % remove CCLE studies because they are cell line data
%     if contains(cancerStudies.name{kk},...
%             {'Cell Line','NCI-60','xenograft','Pediatric','Histio',...
%              },'IgnoreCase',true )
%         fprintf('\n These are CCLE or NCI data - Excluded \n')
%         fprintf(['\n =======================================',...
%             '================ \n\n'])
%         continue
%     end
%     
    % get the data from cBioPortal 
    geneticProfiles = getgeneticprofiles(cgdsURL, ...
        cancerStudies.cancerTypeId{kk}); 
    
    % use the mutation biase to get all the studies with mutation
    if mutBias == true
        % find out if the data contains mutations
        mutsInData = sum( contains(geneticProfiles.geneticProfileName, ...
            {'Mutations'},'IgnoreCase',true ) ) ;
        mutsComp = 1 ;
    else
        % find out if the data contains mutations
        mutsInData = sum( contains(geneticProfiles.geneticProfileName, ...
            {'Mutations','copy-number'},'IgnoreCase',true ) ) ;
        mutsComp = 2;
    end
    
    % now get the data
    if mutsInData >= mutsComp % if the data has copy number data and mutations
        toKeepStudies(kk) = true;
        fprintf('\n PASSED!! The study has mutations and copy number data \n')
    else
        fprintf('\n The study has NO mutations and copy number data \n')
    end
    fprintf('\n ======================================================= \n\n')
end
cancerStudies = cancerStudies( toKeepStudies ,:) ;


% ============ Process the Mutations Data into a Single Table ===========


mutations = [] ; % this also helps to clear the present Mutations data
cancerMuts = zeros(height(cancerStudies),2); % to allocate the cancer mutations
fprintf('\n\n')
for ii = 1:height(cancerStudies)
    % read the cancer study from the text files one after the other. The
    % cancer Ids are found in the cBioPortal table "cancerStudies"
    % get the case list for current cancer study and then get the copy
    % number case List from the struct
    caseLists = getcaselists(cgdsURL, cancerStudies.cancerTypeId{ii} );
   
    % specify a different criterial for the TCGA studies
    if contains(cancerStudies.cancerTypeId{ii},'tcga','IgnoreCase',true)
        caseListId = caseLists.caseListId( ...
            contains(caseLists.caseListId,'complete','IgnoreCase',true) ) ;
    else
        caseListId = caseLists.caseListId( ...
            contains(caseLists.caseListId,'sequenced','IgnoreCase',true) ) ;
    end
    
    % if the study has no mutation data continue  
    if isempty(caseListId)
        continue
        % make sure that the caseListId is not empty
%         assert(~isempty(caseListId))
    end
    % make sure that the caseListId is not empty
    assert(~isempty(caseListId))
    
    % get the available genetic profile for the data and then get the copy
    % number data from the struct
    geneticProfiles = getgeneticprofiles(cgdsURL, ...
        cancerStudies.cancerTypeId{ii}) ;
    
    % specity different criteria for TCGA and other studies
    if contains(cancerStudies.cancerTypeId{ii},'tcga','IgnoreCase',true)
        geneticProfileId = geneticProfiles.geneticProfileId( ...
            contains(geneticProfiles.geneticProfileId ,'mutations') ) ;
    else
        geneticProfileId = geneticProfiles.geneticProfileId( ...
            contains(geneticProfiles.geneticProfileId ,'mutations',...
            'IgnoreCase',true)) ;
        % sometime the cna data in the other studies is also given as
        % gistic
        if isempty(geneticProfileId)
            geneticProfileId = geneticProfiles.geneticProfileId( ...
                contains(geneticProfiles.geneticProfileId ,'mutations') ) ;
        end
    end
    
    % check that there is no mutation data and skip the study
    if isempty(geneticProfileId)
        continue
    end
    % make sure that the genetic profile is not empty
    assert(~isempty(geneticProfileId))
   
    % now get the copy number data from cBioPortal
    tempMut = getprofiledata(cgdsURL, caseListId{1}, ...
        geneticProfileId{1}, myGenes, false) ;
    
    % convert the copy number data to a table
    tempMUTtable = array2table(tempMut.data' ,'VariableNames',...
        matlab.lang.makeValidName(tempMut.common') ) ;
    tempMut = [tempMut.caseId, tempMUTtable] ;

    fprintf('\n\n Now Processing Mutation data for: %s Study # %d\n',...
        cancerStudies.cancerTypeId{ii},ii)
    
    % get the genes which are found in metabolic pathways. Also remove the
    % Entrez Gene Ids and then transponse the table so that the genes are
    % not on top
    tempMut.Properties.VariableNames(1) = {'SampleIds'};
    tempMut = addvars(tempMut,...
        upper(repmat(extractBefore( cancerStudies.cancerTypeId(ii),'_'),...
        height(tempMut),1)) ,'Before','SampleIds','NewVariableNames',...
        'CancerStudy') ;
    
    % check the NaN values to empty strings
    for pp = 2:width(tempMut)
        tempMut.(pp) = strrep(tempMut.(pp),'NaN','') ;   
    end
    
    % remove the duplicate rows of samples Ids
    [~, nonDuplicates] = unique(tempMut.SampleIds);
    tempMut = tempMut(nonDuplicates,:);
    
    % add the total number of samples to cancer mutatins
    cancerMuts(ii,1) = height(tempMut); 
    
    % this is different from the copy number data as not all the metabolic
    % pathway genes are mutated. Therefore I need to add the genes that are
    % not mutated to the table create a table of height = that of tempMut
    % and length = to that of the unmutated genes
    missingGenes = setdiff(myGenes,tempMut.Properties.VariableNames(3:end));
    dummyMut = array2table( cell(height(tempMut),length(missingGenes)) );
    dummyMut.Properties.VariableNames = strrep(missingGenes,'-','_') ;
    tempMut = [tempMut,dummyMut] ;
    
    % add to the growing table
    switch ii
        case ii == 1
            mutations = vertcat(mutations,tempMut) ;
        otherwise
            %     [C,ia,ib] = intersect(___) also returns index vectors ia
            %     and ib using any of the previous syntaxes. Generally, C =
            %     A(ia) and C = B(ib)
            [~,ia,ib] = intersect(mutations.Properties.VariableNames, ...
                tempMut.Properties.VariableNames,'stable') ;
            mutations = mutations(:,ia);
            tempMut = tempMut(:,ib);
            mutations = vertcat(mutations, tempMut) ;
    end
end

% add the total number of samples to the clinical data
cancerStudies = addvars(cancerStudies, cancerMuts(:,1),...
    'NewVariableNames',{'Samples'}) ;

% remove the duplicate rows of samples Ids
[~, nonDuplicates] = unique(mutations.SampleIds);
mutations = mutations(nonDuplicates,:);

% clear some variables
clear tempMut ia ib ii cancerTypeId dummyMut missingGenes



% ====== Get the Gene Expression Profile of the Metabolic Genes ======



mrna = [] ; % this also helps to clear the present Mutations data
% once the loop goes past this point then toAdd will be true
toAddFirst = true;

for ii = 1:height(cancerStudies)
    % read the cancer study from the text files one after the other. The
    % cancer Ids are found in the cBioPortal table "cancerStudies"
    cancerTypeId = string ( extractBefore( ...
        cancerStudies.cancerTypeId(ii), '_tcga')) ;
    try
    tempRNA = readtable( strcat(cancerTypeId,...
        '_tcga_data_RNA_Seq_v2_expression_median') );
    catch % not all studies are TCGA
        continue
    end
    fprintf('Now Processing mRNA Data for: %s Study # %d\n',cancerTypeId,ii)
    fprintf('  The total number of genes is %d \n',...
        length( unique(tempRNA.Hugo_Symbol) ) )
    
    % get the genes which are found in metabolic pathways.
    tempRNA = tempRNA( ismember(tempRNA.Hugo_Symbol, myGenes) ,: ) ;
    fprintf('  Remaining genes is %d\n\n',...
        length( unique(tempRNA.Hugo_Symbol) ) )

    try
        tempRNA{:,3:end} = strrep(tempRNA{:,2:end},'NA','NaN');
    catch
    end
    % change the cell containts to double
    for jj = 2:width(tempRNA)
        if iscell(tempRNA.(jj))
            tempRNA.(jj) = str2double(tempRNA.(jj)) ;
        end
    end
    % transpose the table
    tempRNA = rows2vars(tempRNA,'VariableNamesSource','Hugo_Symbol') ;
    tempRNA.Properties.VariableNames(1) = {'SampleIds'};
    
    % add cancer study
    tempRNA = addvars(tempRNA,repmat(cancerTypeId,height(tempRNA),1) ,...
        'Before','SampleIds','NewVariableNames','CancerStudy') ;
    
    % this is different from the copy number data as not all the metabolic
    % pathway genes are mutated. Therefore I need to add the genes that are
    % not mutated to the table create a table of height = that of tempMut
    % and length = to that of the unmutated genes
    missingGenes =setdiff(myGenes,tempRNA.Properties.VariableNames(3:end));
    dummyRNA = array2table( cell(height(tempRNA),length(missingGenes)) );
    dummyRNA.Properties.VariableNames = strrep(missingGenes,'-','_') ;
    tempRNA = [tempRNA,dummyRNA] ;
    
    % add to the growing table
    switch toAddFirst
        case true
            mrna = vertcat(mrna,tempRNA) ;
            
            % change the value to toAdFirts
            toAddFirst = false;
        otherwise
            %     [C,ia,ib] = intersect(___) also returns index vectors ia
            %     and ib using any of the previous syntaxes. Generally, C =
            %     A(ia) and C = B(ib)
            [~,ia,ib] = intersect(mrna.Properties.VariableNames, ...
                tempRNA.Properties.VariableNames,'stable') ;
            mrna = mrna(:,ia);
            tempRNA = tempRNA(:,ib);
            mrna = vertcat(mrna, tempRNA) ;
    end
end


% =========  Process the Copy Number Data into One Table ============


% This fetches mutation, copy number alterations and expression data for
% all the genes (allGenes) for each patient ID

cnaData = [] ; % this also helps to clear the present CNAdata
% once the loop goes past this point then toAdd will be true
toAddFirst = true;

for ii = 1:height(cancerStudies)
    % read the cancer study from the text files one after the other. The
    % cancer Ids are found in the cBioPortal table "cancerStudies"
    
    % get the case list for current cancer study and then get the copy
    % number case List from the struct
    caseLists = getcaselists(cgdsURL, cancerStudies.cancerTypeId{ii} );
   
    % specify a different criterial for the TCGA studies
    if contains(cancerStudies.cancerTypeId{ii},'tcga','IgnoreCase',true)
        caseListId = caseLists.caseListId( ...
            contains(caseLists.caseListId,'complete','IgnoreCase',true) ) ;
    else
        caseListId = caseLists.caseListId( ...
            contains(caseLists.caseListId,'cna','IgnoreCase',true) ) ;
    end
    
    if isempty(caseListId)
        continue
    end
    
    % get the available genetic profile for the data and then get the copy
    % number data from the struct
    geneticProfiles = getgeneticprofiles(cgdsURL, ...
        cancerStudies.cancerTypeId{ii}) ;
    
    % specity different criteria for TCGA and other studies
    if contains(cancerStudies.cancerTypeId{ii},'tcga','IgnoreCase',true)
        geneticProfileId = geneticProfiles.geneticProfileId( ...
            contains(geneticProfiles.geneticProfileId ,'gistic') ) ;
    else
        geneticProfileId = geneticProfiles.geneticProfileId( ...
            contains(geneticProfiles.geneticProfileId ,'cna',...
            'IgnoreCase',true)) ;
        % sometime the cna data in the other studies is also given as
        % gistic
        if isempty(geneticProfileId)
            geneticProfileId = geneticProfiles.geneticProfileId( ...
                contains(geneticProfiles.geneticProfileId ,'gistic') ) ;
        end
    end
    
    if isempty(geneticProfileId)
        continue
        % make sure that the genetic profile is not empty
%         assert(~isempty(geneticProfileId))
    end
   
   
    % now get the copy number data from cBioPortal
    tempCNA = getprofiledata(cgdsURL, caseListId{1}, ...
        geneticProfileId{1}, myGenes, false) ;
    
    % convert the copy number data to a table
    tempCNAtable = array2table(tempCNA.data' ,'VariableNames',...
        matlab.lang.makeValidName(tempCNA.common') ) ;
    tempCNA = [tempCNA.caseId, tempCNAtable] ;

    fprintf('\n\n Now Processing Copy NumberData for: %s Study # %d\n',...
        cancerStudies.cancerTypeId{ii},ii)
    
    % get the genes which are found in metabolic pathways. Also remove the
    % Entrez Gene Ids and then transponse the table so that the genes are
    % not on top
    tempCNA.Properties.VariableNames(1) = {'SampleIds'};
    tempCNA = addvars(tempCNA,...
        upper(repmat(extractBefore( cancerStudies.cancerTypeId(ii),'_'),...
        height(tempCNA),1)) ,'Before','SampleIds','NewVariableNames',...
        'CancerStudy') ;
    
    % check that the copy number data contains more than 80% of the MAPK
    % genes that are being obtained
    % get the total number of missing genes and remove the sample with more
    % then 20% of the genes missing
    theMissing = sum( isnan(str2double(tempCNA{:,2:end}) ), 2)...
        /length(myGenes);
    tempCNA(theMissing > 0.2, :) = [];
    
    % do away with the cancer study if the number of sample after removing
    % the sample with missing data is less than 20 
    if height(tempCNA) < 40 
        continue 
    end
    
    % check that the data are not decimal points when I have to do away
    % with in my table 
    if any ( any( rem( str2double(tempCNA{:,3:end}) ,1) ) )
        continue 
    end
    
    % find the common genes between the two dataset and then combine the
    % table.
    switch toAddFirst
        case true
            % if the data is a cell array convert to double
            if iscell(tempCNA.(3) )
                for jj = 3:width(tempCNA)
                    tempCNA.(jj) = str2double(tempCNA.(jj)) ;
                end
            end
            cnaData = vertcat(cnaData,tempCNA) ;
            
            % change toAdd to false
            toAddFirst = false ;
        otherwise
            %     [C,ia,ib] = intersect(___) also returns index vectors ia
            %     and ib using any of the previous syntaxes. Generally, C =
            %     A(ia) and C = B(ib)
            [~,ia,ib] = intersect(cnaData.Properties.VariableNames, ...
                tempCNA.Properties.VariableNames,'stable') ;
            cnaData = cnaData(:,ia);
            tempCNA = tempCNA(:,ib);
            
            % sometimes the table variables are cell array and not double
            % therefore I have to take that into account
            try
                cnaData = vertcat(cnaData, tempCNA) ;
            catch
                fprintf('The Data %d is of Cell Type\n',ii)
                for jj = 3:width(tempCNA)
                    tempCNA.(jj) = str2double(tempCNA.(jj)) ;
                end
                cnaData = vertcat(cnaData, tempCNA) ;
            end
    end
end

% remove the 1 and - 1 or - 0
for ii = 3:width(cnaData)
    cnaData.(ii)(cnaData.(ii) == -1|cnaData.(ii) == 1|cnaData.(ii) == -0) = 0  ;
end

clear tempCNA jj cancerTypeId ia ib



% ====  Process The Clinical and Sample Data into a Single Table ========



clinicalData = [] ;
fprintf('\n\n')
for ii = 1:height(cancerStudies)
    
    % get the current cancer study
    cancerTypeId = extractBefore(cancerStudies.cancerTypeId{ii},'_') ;
    
    fprintf('\nNow Processing Clinical Data for: %s Study # %d\n',...
        cancerTypeId,ii)
    
    % get the case list for current cancer study and then get the copy
    % number case List from the struct
    caseLists = getcaselists(cgdsURL, cancerStudies.cancerTypeId{ii} );
   
    % specify a different criterial for the TCGA studies
    if contains(cancerStudies.cancerTypeId{ii},'tcga','IgnoreCase',true)
        caseListId = caseLists.caseListId( ...
            contains(caseLists.caseListId,'complete','IgnoreCase',true) ) ;
    else
        caseListId = caseLists.caseListId( ...
            contains(caseLists.caseListId,'cna','IgnoreCase',true) ) ;
        if isempty(caseListId)
            caseListId = caseLists.caseListId( ...
            contains(caseLists.caseListId,'sequenced','IgnoreCase',true) ) ;
        end
    end
    % make sure that the caseListId is not empty
    assert(~isempty(caseListId))
    
    tempTable = getclinicaldata(cgdsURL, caseListId{1});
    appendTable = array2table([tempTable.caseId, tempTable.data]) ;
    appendTable.Properties.VariableNames = ['SampleIds';...
        matlab.lang.makeValidName(tempTable.clinVariable)] ;
    
    % add cancer study
    appendTable = addvars(appendTable,...
        upper( repmat( cellstr(cancerTypeId), height(appendTable),1) ),...
        'Before','SampleIds','NewVariableNames','CancerStudy') ;
    
    % this is different from the copy number data as not all the
    % clinical data . Therefore I need to add the genes that are
    % have the same names the table create a table of
    if ii ~=1
        missingClinicals = setdiff(...
            clinicalData.Properties.VariableNames, ...
            appendTable.Properties.VariableNames) ;
        dummyClins1 = array2table( cell(height(appendTable),...
            length(missingClinicals)) );
        dummyClins1.Properties.VariableNames = missingClinicals ;
        appendTable = [appendTable,dummyClins1] ;
        
        % also add the dummy clinical data to the growing table
        missingClinicals2 = setdiff(appendTable.Properties.VariableNames,...
            clinicalData.Properties.VariableNames) ;
        dummyClins2 = array2table( cell(height(clinicalData),...
            length(missingClinicals2)) );
        dummyClins2.Properties.VariableNames = missingClinicals2 ;
        clinicalData = [clinicalData,dummyClins2] ;
    end
    
    % add the current clinical data to the overall clinical data
    switch ii
        case ii == 1
            clinicalData = vertcat(clinicalData,appendTable) ;
        otherwise           
            % sometimes the table variables are cell array and not double
            % therefore I have to take that into account
            clinicalData = vertcat(clinicalData, appendTable) ;
    end
end


% ======== Get the common samples for each all genetic profiles =========


% extract the sample Ids before the -01. This is what is making the
% difference. First replace all end of the name with _01
% mrna.SampleIds = strcat( extractBefore(mrna.SampleIds,13), '_01');
% % cnaData.SampleIds = strcat( extractBefore(cnaData.SampleIds,13), '_01');
% % mutations.SampleIds = strcat( extractBefore(mutations.SampleIds,13),'_01');
% clinicalData.SampleIds = strrep( strcat( extractBefore(...
%     clinicalData.SampleIds,13),'_01') , '-', '_');

if mutBias ~= true
    commonIds = intersect(cnaData.SampleIds,mutations.SampleIds);
    fprintf('\nThe number of Pancancer Patients is %d \n',length(commonIds) )
    
    % retain only the row with matching IDs mrna = mrna(...
    % ismember(mrna.SampleIds, commonIds), :) ;
    cnaData = cnaData( ismember(cnaData.SampleIds, commonIds), :) ;
    mutations = mutations( ismember(mutations.SampleIds, commonIds), :);
    % clinicalData = clinicalData(ismember(clinicalData.SampleIds,...
    %     commonIds),:);
end

% convert the cancerType to categories
cancerStudies.cancerTypeId = categorical(cancerStudies.cancerTypeId);


end % end of function


% *********************** end of function ************************
% ****************************************************************

% ======================= another function =========================

% This function returns the alteration frequency for mutations and copy
% number data

function pathwayAlterations = find_MAPK_AlterationFreq(...
    metabolicPathways,mutations,cnaData)
% This is a table with columns as cancer types and row as metabolic
% pathways

% get the genes involved in a metabolic pathwy and return only these for
% the copy number data and mutations data

for ii = 1:height(metabolicPathways)
    
    % get the genes involved
    pathwayGenes = split(metabolicPathways.Genes(ii));
    pathwayMuts = mutations(:,[true, false, ...
        ismember(mutations.Properties.VariableNames(3:end), pathwayGenes)]);
    
    % process the copy number add
    if nargin == 3
        pathwayCopyNumber = cnaData(:,[true, false, ...
            ismember(cnaData.Properties.VariableNames(3:end),...
            pathwayGenes)] );
    end
    
    % get the mutatations in each samples
    pathwayMuts = addvars( pathwayMuts(:,1), ...
        double(any(~cellfun(@isempty,pathwayMuts{:,2:end}),2) ) , ...
        'NewVariableNames','Overall') ;
    pathwayMuts.CancerStudy = categorical(pathwayMuts.CancerStudy);
    
    % combine the two tables: ONLY IF there are also copy number
    % alterations
    if nargin == 3
        if any(ismember(cnaData.Properties.VariableNames(3:end),...
                pathwayGenes))
            
            % also get alterations for the copy number data. Sometimes the
            % data is cell if iam deleteion with GDSC data
            try % for TCGA and CCLE data
                pathwayCopyNumber2 = addvars( pathwayCopyNumber(:,1), ...
                    double( any(pathwayCopyNumber{:,2:end}, 2) ) , ...
                    'NewVariableNames','Overall') ;
            catch % for GDSC data
                pathwayCopyNumber2 = addvars( pathwayCopyNumber(:,1), ...
                    double(any(...
                    ~cellfun(@isempty,pathwayCopyNumber{:,2:end}),2) ),...
                    'NewVariableNames','Overall') ;
                
            end
            pathwayCopyNumber2.CancerStudy = ...
                categorical(pathwayCopyNumber2.CancerStudy);
            
            % add to the total mutations table
            pathwayMuts.Overall =  double( ...
                any([pathwayMuts.Overall,pathwayCopyNumber2.Overall] ,2) ) ;
        end
    end
    % convert the zeroes to NaN to make group stats easiler to do
    pathwayMuts.Overall(pathwayMuts.Overall == 0) = NaN ;
    
    % now get the group stats for the two tables: first get the stats for
    % the table without zeros and then with zeroes
    pathwayMuts = grpstats(pathwayMuts,'CancerStudy','numel') ;
%     pathwayCopyNumber =
%     grpstats(pathwayCopyNumber,'CancerStudy','numel');
    
    % create a table that has the over alterations percentage for both
    % mutations and copy number data
    if ii == 1
        pathwayAlterations = addvars(pathwayMuts(:,1), ...
            round( ...
            pathwayMuts.numel_Overall./pathwayMuts.GroupCount,3)*100,...
            'NewVariableNames', ...
            matlab.lang.makeValidName(metabolicPathways.pathwayName(ii)) );
    else % % join the two tables
        tempAlterations = addvars(pathwayMuts(:,1), ...
            round( ...
            pathwayMuts.numel_Overall./pathwayMuts.GroupCount,3)*100,...
            'NewVariableNames', ...
            matlab.lang.makeValidName(metabolicPathways.pathwayName(ii)) );
        
        pathwayAlterations = innerjoin(pathwayAlterations, tempAlterations);
    end
end

end


% *********************** end of function ************************
% ****************************************************************

% ====================== another function =========================

% get the all genes in are in the 3 pathways from the ucsc pathway

function mTOR = createNetwork2(mTOR)

    fprintf('\n Reading Super Pathway Data \n')
    ucscPathway = readtable('My_MAPK_Super_Pathway.xlsx');

    ogmTOR = mTOR;
    mTOR = ucscPathway( contains(ucscPathway.Protein1,mTOR) & ...
        contains(ucscPathway.Protein2, mTOR) , :) ;

    % ogCellcycle = cellCycle ; cellCycle = ucscPathway(
    % contains(ucscPathway.Protein1,cellCycle) & ...
    %     contains(ucscPathway.Protein2,cellCycle) , :);
    %
    % ogPathwaysInCancer = pathwaysInCancer; pathwaysInCancer =
    % ucscPathway( contains(ucscPathway.Protein1,...
    %     pathwaysInCancer) &
    %     contains(ucscPathway.Protein2,pathwaysInCancer),:);

    % ========== These results look very good so far!!!!!!! ============
    % clear locK2 locK pathways GOKegg1 GOKegg2 GoKeggC1Terms GoKeggC2Terms
    % ...
    %     GoBioC1Terms GoBioC2Terms GObio1 GObio2 to_go combinedScores ...
    %     goTerm termEnd ucscPathway

    % ================== Now create the mTOR pathaway ==============

    % I will have to start with a simpler (smaller) graph create a graph
    % from edge: first delete the selfloop nodes
    selfLoop = strcmpi(mTOR.Protein1 , mTOR.Protein2) ;
    mTOR(selfLoop,:) = [] ;

    % create table for yED that also contains the origanal proteins and
    % also add the differential gene expression log p-value to be used to
    % colour the notes
    ogProteins = contains(mTOR.Protein1 ,ogmTOR) ;
    mTOR.NodeInGO = double(ogProteins) ;

    % remove that bad arrow from the data
    bad = contains(mTOR.Interaction,'-t>');
    mTOR.Interaction(bad,1) = {'->t'};

    % now create a graph
    mtorGraph = digraph(mTOR.Protein1 , mTOR.Protein2);
    mtorGraph.Edges.Interaction = mTOR.Interaction ;

    % plot the graph
    figure()
    hMTOR = plot(mtorGraph,'layout','force','usegravity',false,...
        'MarkerSize',10,'ArrowSize', 6,'EdgeAlpha',0.80 ,...
        'LineWidth', 0.5000);
    set(gca,'FontSize',12,'FontWeight','bold','visible', 'off')
    title('mTOR-Network','FontSize',18)
    hold on
    % get the nodes that have edge for interactions from biogrid
    allInters =  unique(mtorGraph.Edges.Interaction) ;
    for ii = 1:length(allInters)
        cur_inter = allInters(ii,1) ;
        locsG = contains(mtorGraph.Edges.Interaction,cur_inter);
        [sOut,tOut] = findedge(mtorGraph);
        allEdges = [sOut,tOut];
        % check = mtorGraph.Edges(locsG,:) ;
        subGraph = allEdges(locsG,:) ;
        subGraph = reshape( subGraph',1,[]) ;
        % if the interaction is just protein-protein
        if strcmp(cur_inter,'->i')
            highlight(hMTOR,subGraph,'EdgeColor',[0.73 0.49 0.43], ...
                'LineWidth',1.5 ,'LineStyle','--') % ,'ArrowPosition',1)
        elseif strcmp(cur_inter,'->p')
            highlight(hMTOR,subGraph,'EdgeColor','b','LineWidth',2)
        elseif strcmp(cur_inter,'-a>')
            highlight(hMTOR,subGraph,'EdgeColor',[0.32 0.79 0.43],...
                'LineWidth',2)
        elseif strcmp(cur_inter,'-a|')
            highlight(hMTOR,subGraph,'EdgeColor','r','LineWidth',2)
        else
            highlight(hMTOR,subGraph,'EdgeColor',[0.5 0.5 0.5],...
                'LineWidth',1.5)
        end

    end
    hold off
end

% *********************** end of function ************************
% ****************************************************************

% ======================= another function =========================

% plot the mutation profile of the samples on a grid
function exitPos = oncoPrintInternal(inData, colors, rowNames, ...
    horizonBars, nextAxis)
% Input:
% inData: a matrix and vector of plotting data
% colors: colors for each unique value of inData
% rowNames: a cell array for name for each row
% horizonBar: specifies whether to add a horizontal bar to each plot or not

% get the number of unique numbers and remove the 0 which means not plot
uniqueVars = unique(inData) ;

% create the legend variable if the inData is a row vector
if size(inData,1) == 1
    lgdVar = split( num2str(uniqueVars) )';
end

% get the number of colours to plot
if ~exist('colors','var') || isempty(colors)
    colors = rand(length(uniqueVars), 3);
end

% check the are is a row name for each row in inData
if size(inData,1) > 1 % only valid for matrix data
    if size(inData,1) ~= size(rowNames,1)
        error('row names should contain a name for each row in plotting data')
    end
end

% check if the orizontal bars are required on the plot
if ~exist('horizonBars','var')
    horizonBars = false;
end

% check for color errors
if length(uniqueVars) > size(colors,1)
    error('A color must be specified for each value in Data')
end

% make a plot of multiple bar chart of top of each other
% add the subtype plots to the heatmap using a loop
% initialise some variables
global plotTime
figure(plotTime*100)

yInitial = 0.15; yPosSaved = yInitial;
ySize = 0.025; increaseby = 0.0225; % 0.44 0.4
xEndPos = 0.7 ; % 0.7750
% position of the percentage and bargraph
xEndPos1 = xEndPos+0.10;
xEndPos2 = xEndPos1+0.05;
% loop over the rows and ascend by column
for jj = 1:size(inData,1) % begin plotting from the bottem
    % define the next axis for the follow up plots
    if ~exist('nextAxis','var') || isempty(nextAxis)
        axes('position',[0.1300,yInitial,xEndPos,ySize]);
    elseif exist('nextAxis','var') && jj > 1
        axes('position',[0.1300,yInitial,xEndPos,ySize]);
    else
        axes('position',nextAxis);
        yInitial = nextAxis(2) ;
        ySize = nextAxis(4) ; xEndPos = nextAxis(3) ;
        yPosSaved = yInitial ;
    end
    for ii = 1:numel(uniqueVars)
        plotData = double(ismember(inData(jj,:),uniqueVars(ii) )) ;
        bar(plotData,'FaceColor',colors(ii,:),'EdgeColor',[1 1 1] ,...
            'BarWidth',0.9) ;
        hold on
        
        % add the name of the genes to the left of heatmap
        if exist('rowNames','var')
            dim = [0.02 yInitial 0.11 increaseby];
            annotation('textbox',dim,'String',rowNames{jj},...
                'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
                'HorizontalAlignment','right','FontWeight','bold',...
                'VerticalAlignment','middle');
        end
        % add the percentage of mutated genes to the right of heatmap
        dim = [xEndPos1 yInitial 0.05 increaseby];
        percMut = sum(inData(jj,:) > 1)/length(inData(jj,:))*100;
        annotation('textbox',dim,'String',...
            strcat(num2str(round(percMut)),'%'),...
            'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
            'HorizontalAlignment','right','FontWeight','bold',...
            'VerticalAlignment','middle');
    end
    % change the plot properties
    set(gca,'GridColor',[1 1 1], 'XLim', [0.5  size(inData,2)+0.5],...
        'XColor',[1 1 1] ,'YColor',[1 1 1],'YTickLabel',[],...
        'XTickLabel',[],'FontWeight','bold','YTick',[],'XTick',[])
    % increase the value to change the colors and plot positions
    yInitial = yInitial + increaseby;
end
% add a grey box to the plot usign the annotation function
% dim = [0.1300, yPosSaved, xEndPos, increaseby*size(inData,1)];
% annotation('rectangle',dim ,'Color',[0.5, 0.5, 0.5])
% hold off

% plot the horizontal bars if they are required
% prellocate the bar data size
barhData = zeros(size(inData,1),numel(uniqueVars)) ;
if horizonBars == true
    axes('position',[xEndPos2, yPosSaved, 0.10, ...
        increaseby*size(inData,1) ]);
    for kk = 2:numel(uniqueVars)
        barhData(:,kk) = sum(inData == uniqueVars(kk),2) ;
    end
    bar1 = barh(barhData,'stacked','BarWidth',0.85) ;
    % make sure there are no colors and spaces between the axis and the
    % first and last bar
    set(gca,'GridColor',[1 1 1], 'YLim', [0.5 size(inData,1)+0.52 ], ...
        'XColor',[0.3 0.3 0.3] ,'YColor',[1 1 1],'FontSize',10,...
        'YTick',[],'FontWeight','bold','Box','off','TickDir', 'out',...
        'LineWidth',1,'XAxisLocation','origin')
    
    % annoate the bar graph
    for ii = 1:size(barhData,2)
        set(bar1(ii),'FaceColor',colors(ii,:))
    end   
end

% add a vertical bar graph to the plot
axes1 = axes('position',[0.1300,yPosSaved-0.137,xEndPos,0.13]);
hold(axes1,'on');
barhData = zeros(numel(uniqueVars),size(inData,2)) ;
for kk = 2:numel(uniqueVars)
    barhData(kk,:) = sum(inData == uniqueVars(kk),1) ;
end
bar1 = bar(barhData','stacked','BarWidth',0.85) ;
% make sure there are no colors and spaces between the axis and the
% first and last bar
set(gca,'GridColor',[1 1 1], 'XLim', [0.5 size(inData,2)+0.52 ],...
    'XColor',[1 1 1] ,'YColor',[0.3 0.3 0.3],'FontSize',10,...
    'XTick',[],'FontWeight','bold','Box','off','TickDir', 'out',...
    'LineWidth',1,'XAxisLocation','origin','YDir','reverse',...
    'YAxisLocation','right')

% annoate the bar graph
for ii = 1:size(barhData,1)
    set(bar1(ii),'FaceColor',colors(ii,:))
end

% save the exit positions
exitPos = [0.1300,yInitial,xEndPos, ySize] ;

end

% *********************** end of function ************************

% ****************************************************************

% ==================== Another Internal Function ======================

% This function return mrna expressiond data from the cBioportal given a
% list of cancers

function mrnaData = getMedianTranscritomicsData(cancerTypes, reduceSize)

% redcue the size of the genes to only a few if reduce size is true
if reduceSize == true
    % get only the landmark genes from the tall array
    landMarkGenes = readtable('landmarkGenes.txt');
end
    
for ii = 1:length(cancerTypes)
    % read the cancer study from the text files one after the other. The
    % cancer Ids are found in the cBioPortal table "cancerStudies"
    cancerTypeId = lower(cancerTypes{ii});
    fprintf('Now processing mRNA data for: %s Study # %d\n',cancerTypeId,ii)
    tempRNA = readtable( strcat(cancerTypeId,...
        '_tcga_data_RNA_Seq_v2_expression_median') ) ;
    fprintf('  The total number of genes is %d \n',...
        length( unique(tempRNA.Hugo_Symbol) ) )
    
    % remove the entrenz gene symbol and convert the cell array data in the
    % table to double and also remove the genes with Hugo gene symbols from
    % the data 
    tempRNA.Entrez_Gene_Id = [] ;
    tempRNA( ismissing(tempRNA(:,1)) , :) = [] ;
    
    % redcue the size of the genes to only a few if reduce size is true
    if reduceSize == true
        % get only the landmark genes
        tempRNA = tempRNA(ismember(tempRNA.Hugo_Symbol,...
            landMarkGenes.Symbol ) , :) ;
    end
    
    try
        tempRNA{:,3:end} = strrep(tempRNA{:,2:end},'NA','NaN');
    catch
    end
    % change the cell containts to double
    for jj = 2:width(tempRNA)
        if iscell(tempRNA.(jj))
            tempRNA.(jj) = str2double(tempRNA.(jj)) ;
        end
    end
    % transpose the table
    tempRNA = rows2vars(tempRNA,'VariableNamesSource','Hugo_Symbol') ;
    tempRNA.Properties.VariableNames(1) = {'SampleIds'};
    
    % add cancer study
%     tempRNA =  addvars(tempRNA,repmat(cellstr( ...
%         cancerTypeId),height(tempRNA),1) ,...
%         'Before','SampleIds','NewVariableNames','CancerStudy') ;
    
    % #############################################################
    % this is only valid for the median expression levels 
    tempRNA.SampleIds = [] ;
    geneNames = tempRNA.Properties.VariableNames(2:end) ;
    tempRNA = array2table( mean(tempRNA{:,2:end} ,1)  , ...
        'VariableNames',geneNames) ;
    % add the cancer study
    tempRNA =  addvars(tempRNA,repmat(cellstr( ...
        cancerTypeId),height(tempRNA),1) ,...
        'Before',1,'NewVariableNames','CancerStudy') ;
    
    % ###############################################################
    
    % add to the growing table
    switch ii
        case ii == 1
            mrnaData = tempRNA ;
        otherwise
            %     [C,ia,ib] = intersect(___) also returns index vectors ia
            %     and ib using any of the previous syntaxes. Generally, C =
            %     A(ia) and C = B(ib)
            [~,ia,ib] = intersect(mrnaData.Properties.VariableNames, ...
                tempRNA.Properties.VariableNames,'stable') ;
            mrnaData = mrnaData(:,ia);
            tempRNA = tempRNA(:,ib);
            mrnaData = vertcat(mrnaData, tempRNA) ;
    end
end

end



