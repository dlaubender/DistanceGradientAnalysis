function Overlay_20MicronContourLines

% This function allows to select and load in ROIed Lymphatic AB images and  
% overlays color-coded 20,40,60,80 and 100um contour lines from the ROI 
% borders.
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
% - Contour overlay images (Lymphatic and Cellular stain)
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
    Lymph_image{i} = imread([FileName_temp(1:end-ExpLabelLength) '0.tif']);
    Cell_image{i} = imread([FileName_temp(1:end-ExpLabelLength) '1.tif']);
    clear FileName_temp
    clear ImInfo_temp
end


%% calculate and plot 20 microns contour lines from Lymphatic ROI for each selected image

for j = 1:Nr_files

    % create distance mask around Lymph Vessel mask (round to nearest decimal)
    DistanceMask = round(bwdist(LymphMasks{j,1}.ROI_Mask,'euclidean')*Pixel2Microns);

    % extract contour lines in 20 micron steps up to 100 microns
    [CL20_Y,CL20_X] = find(DistanceMask == 20);
    [CL40_Y,CL40_X] = find(DistanceMask == 40);
    [CL60_Y,CL60_X] = find(DistanceMask == 60);
    [CL80_Y,CL80_X] = find(DistanceMask == 80);
    [CL100_Y,CL100_X] = find(DistanceMask == 100);

    % find pixels of ROI and Perimeter
    Lymph_mask_perim = bwperim(LymphMasks{j,1}.ROI_Mask);
    [ROIperim_IDs_Y,ROIperim_IDs_X] = find(Lymph_mask_perim);
    
    % define colormap for colormap indexing
    Colormap = cool(5);

    % plot cell picture and overlay contour lines
    f1 = figure;
    imshow(Cell_image{j});
    hold on
    plot(ROIperim_IDs_X,ROIperim_IDs_Y,'.w','MarkerSize',1.5)
    plot(CL20_X,CL20_Y,'.','Color',[Colormap(1,:)],'MarkerSize',1.5)
    plot(CL40_X,CL40_Y,'.','Color',[Colormap(2,:)],'MarkerSize',1.5)
    plot(CL60_X,CL60_Y,'.','Color',[Colormap(3,:)],'MarkerSize',1.5)
    plot(CL80_X,CL80_Y,'.','Color',[Colormap(4,:)],'MarkerSize',1.5)
    plot(CL100_X,CL100_Y,'.','Color',[Colormap(5,:)],'MarkerSize',1.5)
%     title('20 micron contour lines')
    axis off
    hold off

    % save figure
    saveas(f1,[char(FileNames(j)) '_Cell_ROI_20umContour.fig']);
    export_fig(f1,[char(FileNames(j)) '_Cell_ROI_20umContour.png']);
%     saveas(f1,[char(FileNames(j)) '_Cell_ROI_20umContour.pdf']);
    
    
    % plot lymph picture and overlay contour lines
    f2 = figure;
    imshow(Lymph_image{j});
    hold on
    plot(ROIperim_IDs_X,ROIperim_IDs_Y,'.w','MarkerSize',1.5)
    plot(CL20_X,CL20_Y,'.','Color',[Colormap(1,:)],'MarkerSize',1.5)
    plot(CL40_X,CL40_Y,'.','Color',[Colormap(2,:)],'MarkerSize',1.5)
    plot(CL60_X,CL60_Y,'.','Color',[Colormap(3,:)],'MarkerSize',1.5)
    plot(CL80_X,CL80_Y,'.','Color',[Colormap(4,:)],'MarkerSize',1.5)
    plot(CL100_X,CL100_Y,'.','Color',[Colormap(5,:)],'MarkerSize',1.5)
%     title('20 micron contour lines')
    axis off
    hold off

    % save figure
    saveas(f2,[char(FileNames(j)) '_Lymph_ROI_20umContour.fig']);
    export_fig(f2,[char(FileNames(j)) '_Lymph_ROI_20umContour.png']);
%     export_fig(f2,[char(FileNames(j)) '_Lymph_ROI_20umContour.pdf']);
end





