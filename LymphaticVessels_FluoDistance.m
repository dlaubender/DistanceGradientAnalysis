function LymphaticVessels_FluoDistance

% This function will calculate a euclidean distance mask and extract
% average fluorescence values as a function of distance. It will let you
% select the images from the second channel (relative to the bloodvessel
% channel) to apply the distance mask to.
%
% Requirements:
% - Drawn and saved Lymphatic ROIs (see function DrawSave_LymphaticROIs.mat)
%
% Parameters to set:
% - Pixel2Microns:      Pixel to micron conversion factor
% - ExpLabelLength:     Set this depending on image labeling between
%                       channels (length from end till common label)
% 
% Outputs:
% - Figure:             Graph plotting normalized distance dependent 
%                       fluorecence (saved as .fig, .svg and .pdf)
% - Data Matrices:      Distance dependent fluorescent values are saved as 
%                       .mat and .xlsx files
%
%
% David Laubender 2020


%% parameter switchboard

% Set Pixel to micron resolution
Pixel2Microns = 1/1.5;

% Define number of characters used for differentiating different imaging 
% channels (including ".tif" --> e.g. here: ...C0.tif and ...C1.tif 
% --> last 5 characters)
ExpLabelLength = 5;

% name experimental conditions
Name_StainingOI = inputdlg('Define name for staining of interest (no spaces)','Name experimental conditions');
Name_IsoControl = inputdlg('Define name for control (no spaces)','Name experimental conditions');

% define number of staining protocols: staining of interest and negative control
Nr_StainingProtocols = 2;

% preallocate
DistanceStatistics = [];
for m = 1:Nr_StainingProtocols            
    DistanceStatistics(m).RawDistancePixelFs = [];
    DistanceStatistics(m).NormDistancePixelFs = [];
    DistanceStatistics(m).AllNormDistancePixelFs = [];
    DistanceStatistics(m).AllRawDistancePixelFs = [];
    DistanceStatistics(m).DistanceMeans = [];
    DistanceStatistics(m).DistanceSTDs = [];
    DistanceStatistics(m).DistanceSEMs = [];
    DistanceStatistics(m).GlobalDistanceMeans = [];
    DistanceStatistics(m).GlobalDistanceSTDs = [];
    DistanceStatistics(m).GlobalDistanceSEMs = [];
    DistanceStatistics(m).normGlobalDistanceMeans = [];
    DistanceStatistics(m).normGlobalDistanceSTDs = [];
    DistanceStatistics(m).normGlobalDistanceSEMs = [];
end

tic

for k = 1:Nr_StainingProtocols         
    %% Load data by selection

    % select all Interstitium images
    if k == 1
        [FileNames,FilePath] = uigetfile('.tif','Select staining of interest Interstitium images','MultiSelect','on');
    else
        [FileNames,FilePath] = uigetfile('.tif','Select control Interstitium images','MultiSelect','on');                
    end
    Nr_files = length(FileNames);

    % change directory to filepath
    cd(FilePath)

    % preallocate
    AB_Interstitium = cell(Nr_files,1);
    LymphMasks = cell(Nr_files,1);
%     AB_Lymph = cell(Nr_files,1);

    % load Interstitium images and Lymph ROIs with Pixel Conversion Factors
    for i = 1:Nr_files
        FileName_temp = FileNames{i};
        AB_Interstitium{i} = double(imread(FileName_temp));
        LymphMasks{i} = load(['ROI_Mask_' FileName_temp(1:end-ExpLabelLength) '0.mat']);

        clear FileName_temp
        clear ImInfo_temp
    end


    %% calculate and plot Distance Mask from Lymphatic ROI
    
    % determine longest distance from Lymph vessel from all pictures for preallocation
    for h = 1:Nr_files
        DistanceLength = length(unique(round(bwdist(LymphMasks{h,1}.ROI_Mask,'euclidean')*Pixel2Microns)));
        if h == 1
            MaxDistance = DistanceLength;
        elseif DistanceLength > MaxDistance
            MaxDistance = DistanceLength;
        end
    end
    
    % determine longest distance from Lymph vessel from all pictures for preallocation
    PixelFrequency = zeros(MaxDistance,1);
    for h = 1:Nr_files
        [~,~,PixelFrequency_ic] = unique(round(bwdist(LymphMasks{h,1}.ROI_Mask,'euclidean')*Pixel2Microns));
        PixelFrequency_temp = accumarray(PixelFrequency_ic,1);
        PixelFrequency(1:length(PixelFrequency_temp)) = PixelFrequency(1:length(PixelFrequency_temp)) + PixelFrequency_temp;
    end
    MaxPixelFrequency = max(max(PixelFrequency));
    
    % preallocate
    DistanceStatistics(k).AllNormDistancePixelFs = NaN(MaxPixelFrequency,MaxDistance);
    DistanceStatistics(k).AllRawDistancePixelFs = NaN(MaxPixelFrequency,MaxDistance);
    if k == 1
        Excel_StainingOI = NaN(MaxDistance,Nr_files);
    elseif k == 2
        Excel_IsoControl = NaN(MaxDistance,Nr_files);
    end
    
    for j = 1:Nr_files

        % create distance mask around Lymph Vessel mask (round to nearest decimal)
        DistanceMask = round(bwdist(LymphMasks{j,1}.ROI_Mask,'euclidean')*Pixel2Microns);


        %% extract distance dependent fluorescence values (normalized to F(max)) of Interstitium (relative from Lymph Vessel)

        % preallocate
        Fluorescence_Dist = cell(1,length(unique(DistanceMask))-1);
        MEAN_FDist = NaN(1,length(unique(DistanceMask))-1);
        STD_FDist = NaN(1,length(unique(DistanceMask))-1);
        SEM_FDist = NaN(1,length(unique(DistanceMask))-1);

        % find F(max) for normalization
        Fmax = max(max(AB_Interstitium{j}(~LymphMasks{j,1}.ROI_Mask)));
       
        % extract distance dependent fluorescence values & stats
        Mask_Distances = unique(DistanceMask);
        for i = 1:length(Mask_Distances)
            ID = Mask_Distances(i);
            if ID > 0
                DistanceIDs = (DistanceMask == ID);
                Fluorescence_Dist{ID} = AB_Interstitium{j}(DistanceIDs);
                MEAN_FDist(ID) = nanmean(Fluorescence_Dist{ID});
                STD_FDist(ID) = nanstd(Fluorescence_Dist{ID});
                SEM_FDist(ID) = nanstd(Fluorescence_Dist{ID})/sqrt(length(Fluorescence_Dist{ID}));
                clear DistanceIDs

                % save variables
                DistanceStatistics(k).RawDistancePixelFs{j,ID} = Fluorescence_Dist{ID};
                DistanceStatistics(k).NormDistancePixelFs{j,ID} = Fluorescence_Dist{ID}/Fmax;
                if j == 1
                    DistanceStatistics(k).AllNormDistancePixelFs(1:length(Fluorescence_Dist{ID}),ID) = Fluorescence_Dist{ID}/Fmax;
                    DistanceStatistics(k).AllRawDistancePixelFs(1:length(Fluorescence_Dist{ID}),ID) = Fluorescence_Dist{ID};
                else
                    First_NaN_position = find(isnan(DistanceStatistics(k).AllNormDistancePixelFs(:,ID)) == 1, 1,'first');
                    DistanceStatistics(k).AllNormDistancePixelFs(First_NaN_position:First_NaN_position+length(Fluorescence_Dist{ID})-1,ID) = Fluorescence_Dist{ID}/Fmax;
                    DistanceStatistics(k).AllRawDistancePixelFs(First_NaN_position:First_NaN_position+length(Fluorescence_Dist{ID})-1,ID) = Fluorescence_Dist{ID};
                    clear First_NaN_position
                end
                DistanceStatistics(k).DistanceMeans(j,ID) = MEAN_FDist(ID);
                DistanceStatistics(k).DistanceSTDs(j,ID) = STD_FDist(ID);
                DistanceStatistics(k).DistanceSEMs(j,ID) = SEM_FDist(ID);
                if k == 1
                    Excel_StainingOI(ID,j) = MEAN_FDist(ID);
                elseif k == 2
                    Excel_IsoControl(ID,j) = MEAN_FDist(ID);
                end
            end
            clear ID
        end   
    end

    % calculate Means, STDs and SEMs for each distance over all pictures (repetitions)
    for i = 1:length(DistanceStatistics(k).AllNormDistancePixelFs(1,:))
        DistanceStatistics(k).GlobalDistanceMeans(i) = nanmean(DistanceStatistics(k).AllRawDistancePixelFs(:,i));
        DistanceStatistics(k).GlobalDistanceSTDs(i) = nanstd(DistanceStatistics(k).AllRawDistancePixelFs(:,i));
        DistanceStatistics(k).GlobalDistanceSEMs(i) = nanstd(DistanceStatistics(k).AllRawDistancePixelFs(:,i))/sqrt(sum(~isnan(DistanceStatistics(k).AllRawDistancePixelFs(:,i))));
    end
end


%% Normalize Means, STDs and SEMs to Means (non-control) for plotting

for m = 1:Nr_StainingProtocols
    if m == 1
        DistanceStatistics(m).normGlobalDistanceMeans = DistanceStatistics(m).GlobalDistanceMeans/max(DistanceStatistics(m).GlobalDistanceMeans);
        DistanceStatistics(m).normGlobalDistanceSTDs = DistanceStatistics(m).GlobalDistanceSTDs/max(DistanceStatistics(m).GlobalDistanceMeans);
        DistanceStatistics(m).normGlobalDistanceSEMs = DistanceStatistics(m).GlobalDistanceSEMs/max(DistanceStatistics(m).GlobalDistanceMeans);
    else
        DistanceStatistics(m).normGlobalDistanceMeans = DistanceStatistics(m).GlobalDistanceMeans/max(DistanceStatistics(m-1).GlobalDistanceMeans);
        DistanceStatistics(m).normGlobalDistanceSTDs = DistanceStatistics(m).GlobalDistanceSTDs/max(DistanceStatistics(m-1).GlobalDistanceMeans);
        DistanceStatistics(m).normGlobalDistanceSEMs = DistanceStatistics(m).GlobalDistanceSEMs/max(DistanceStatistics(m-1).GlobalDistanceMeans);        
    end
end

    
%% plot distance dependent fluorescence values (normalized to F(max)) of Interstitium (relative from Lymph Vessel)

% plot mean F with SEM 
f2 = figure('Name','Distance Dependent F','color','w');
hold on
errorbar((1:length(DistanceStatistics(1).normGlobalDistanceMeans)),DistanceStatistics(1).normGlobalDistanceMeans,DistanceStatistics(1).normGlobalDistanceSEMs,'.','MarkerSize',15)
errorbar((1:length(DistanceStatistics(2).normGlobalDistanceMeans)),DistanceStatistics(2).normGlobalDistanceMeans,DistanceStatistics(2).normGlobalDistanceSEMs,'.','MarkerSize',15)
title('Distance Dependent Fluorescence');
xlabel('Distance from Lymphvessel [um]');
ylabel('Normalized Fluorescence Intensity');
legend('Staining of interest','iso control')
xlim([0 100]);
ylim_curr = get(gca,'ylim'); 
ylim_curr(2) = 1;
set(gca,'ylim',ylim_curr)
hold off

toc


%% save plot and data matrix

% save plot
saveas(f2,['DistanceDependentF_' char(Name_StainingOI) '.fig']);
saveas(f2,['DistanceDependentF_' char(Name_StainingOI) '.svg']);
saveas(f2,['DistanceDependentF_' char(Name_StainingOI) '.pdf']);

% save data matrix
save(['DistanceDependentF_ExcelMeanMatrices_' char(Name_StainingOI) '.mat'],'Excel_StainingOI','Excel_IsoControl')
save(['DistanceDependentF_DataMatrix_' char(Name_StainingOI) '.mat'],'DistanceStatistics','-v7.3')
xlswrite(['DistanceDependentF_' char(Name_StainingOI) '.xlsx'],Excel_StainingOI);
xlswrite(['DistanceDependentF_' char(Name_IsoControl) '.xlsx'],Excel_IsoControl);



