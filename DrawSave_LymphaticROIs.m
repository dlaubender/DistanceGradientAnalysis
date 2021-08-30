function DrawSave_LymphaticROIs

% This function allows for freehand drawing of ROIs on top of the Lymphatic
% AB image and saves the results as a binary mask.
%
% David Laubender 2020



close all
clear all

% select all tifs to be ROI-ed
[FileNames,~] = uigetfile('.tif','Select all images to be ROI-ed','MultiSelect','on');    

% enter number of picture to be ROI-ed
Nr_Reps = length(FileNames);
clear FileNames


% draw and extract ROIs for each picture
for i = 1:Nr_Reps

    % select & load Lymphatic AB tiff
    [FileName,FilePath] = uigetfile({'.tif'},'Select an image');
    AB_Lymph = imread(strcat(FilePath,FileName));

    % define number of necessary ROIs to be drawn based on picture
    imshow(AB_Lymph)
    answer = inputdlg('Define number of necessary ROIs','Lymphatic Sample');
    Nr_ROIs = str2num(answer{1});
    clear answer
    close gcf

    % preallocate ROI structure
    ROIs = cell(1,Nr_ROIs);

    % draw ROIs for binary mask on lymphatic antibody picture
    for i = 1:Nr_ROIs
        imshow(AB_Lymph)
        ROI = drawfreehand();
        ROIs{i} = ROI.createMask();
    end

    % combine all ROIs
    for j = 1:Nr_ROIs
        if j == 1
            ROI_Mask = ROIs{j};
        else
            ROI_Mask = ROI_Mask + ROIs{j};
        end
    end

    % binarize in case of overlap
    ROI_Mask(ROI_Mask > 1) = 1;
    ROI_Mask = logical(ROI_Mask);

    % Save ROI Mask for later processing
    cd(FilePath)
    save(['ROI_Mask_' FileName(1:end-4) '.mat'],'ROI_Mask');
    close gcf
end