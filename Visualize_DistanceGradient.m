function Visualize_DistanceGradient

% This function allows to select and load in ROIed Lymphatic AB images and  
% overlays a color-coded distance gradient from the ROI borders.
%
% Requirements:
% - Drawn and saved Lymphatic ROIs (see function DrawSave_LymphaticROIs.mat)
%
% Parameters to set:
% - Pixel2Microns:      Pixel to micron conversion factor
% - ExpLabelLength:     Set this depending on image labeling between
%                       channels (length from end till common label)
%
% Output:
% - Distance gradient image
%
% David Laubender 2020


%% parameter switchboard

% Set Pixel to micron resolution
Pixel2Microns = 1/1.5;

% Define number of characters used for differentiating different imaging 
% channels (including ".tif" --> e.g. here: ...C0.tif and ...C1.tif 
% --> last 5 characters)
ExpLabelLength = 5;

        
%% select and load all relevant data

% select all Cell/Interstitium images
[FileNames,FilePath] = uigetfile('.tif','Select Cell images','MultiSelect','on');
Nr_files = length(FileNames);

% change directory to filepath
cd(FilePath)

% preallocate
LymphMasks = cell(Nr_files,1);
Cell_image = cell(Nr_files,1);
Lymph_image = cell(Nr_files,1);

% load Cell/Interstitium images and Lymph ROIs with Pixel Conversion Factors
for i = 1:Nr_files
    FileName_temp = FileNames{i};
%         ImInfo_temp = imfinfo(FileName_temp);
%         Pixel2Microns(i) = 1/ImInfo_temp.XResolution;
    LymphMasks{i} = load(['ROI_Mask_' FileName_temp(1:end-ExpLabelLength) '0.mat']);
    Cell_image{i} = imread([FileName_temp(1:end-ExpLabelLength) '0.tif']);
    clear FileName_temp
    clear ImInfo_temp
end


%% calculate and plot distance gradient from Lymphatic ROI for each selected image

for j = 1:Nr_files

    % create distance mask around Lymph Vessel mask (round to nearest decimal)
    DistanceMask = round(bwdist(LymphMasks{j,1}.ROI_Mask,'euclidean')*Pixel2Microns);
    
    % find pixels of Lymph Vessel Mask and Perimeter
    [ROI_IDs_Y,ROI_IDs_X] = find(LymphMasks{j,1}.ROI_Mask == 1);
    Lymph_mask_perim = bwperim(LymphMasks{j,1}.ROI_Mask);
    [ROIperim_IDs_Y,ROIperim_IDs_X] = find(Lymph_mask_perim);

    % plot distance mask
    f1 = figure;
    imagesc(DistanceMask)
    hold on
    colormap(gca,'cool');
    c = colorbar;
    c.Label.String = 'Distance to nearest LymphVessel [um]';
    plot(ROI_IDs_X,ROI_IDs_Y,'.k')
    plot(ROIperim_IDs_X,ROIperim_IDs_Y,'.w','MarkerSize',1.5)
    title('Distance Mask')
    axis square
    axis off
    hold off

    % save figure
    saveas(f1,[char(FileNames(j)) '_ROI_DistanceMask.fig']);
    saveas(f1,[char(FileNames(j)) '_ROI_DistanceMask.svg']);
    saveas(f1,[char(FileNames(j)) '_ROI_DistanceMask.pdf']);
end





