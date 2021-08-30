# DistanceGradientAnalysis
### This repository includes the essential code for ROI-dependent distance fluorescence gradient analysis for the publication Holtkamp et al. (2021)

These scripts are not written as a general toolbox, but are tailored to the publication's requirements and raw-data labeling logic.

/////////////////////////////////////////////////////

David Laubender 2020

written in Matlab R2018b

Associated publication: Holtkamp et al. (2021)

////////////////////////////////////////////////////

# Analysis logic:

First, name image files per condition distincly with the channel-identifier in the end. In this case, two spectrally distinct antibodies were used. 
The first (channel 0 --> "C0") stained against lymphatic vessels and the second (channel 1 --> "C1") against cells of interest. 
Lymphatic vessels were used to draw ROI masks, which were then applied to the second channel respectively to analyze the location of the cells of interest 
as a function of distance.

#### 1) DrawSave_LymphaticROIs.m
#### 2) LymphaticVessels_FluoDistance.m


///////////////////////////////////////////////////

# Additional scripts:

These two scripts allow to visualize the ROI-dependent distance gradients either smoothly or overlayed on the raw images as 20um contour distance lines.

#### 3) Visualize_DistanceGradient.m
#### 4) Overlay_20MicronContourLines.m


