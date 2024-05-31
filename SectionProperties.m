%Alexander Walsh 5/30/2024
clear
format compact

%This script is used to find the section properties of an arbitrary shape
%(moments of inertia, area, centroid). To use it:
%1. Create a blank white png image RGB(255,255,255) using a tool such as
%MS paint or any image editing program
%2. Draw your section in black RGB(0,0,0)
%3. Add a horizontal scale line anywhere on the image in red RGB(255,0,0)
%Any other colors are ignored, only pure white, black, and red will work.

%Define your inputs as show in the following example:

%{
    %Define the image file
    image = imread("ExampleSectionImage01.png");
    %Define the length of the red horizontal scale line
    scaleLength = 3;
    %If you want to find the moments around some point other than the
    geometric centroid of the shape, set forceCentroid as true and input
    the coordinates you want to use in inches (otherwise set it to false).
    forceCentroid = false;
    xCentroid = 0;
    yCentroid = 0;
%}

%============================== USER INPUT ================================

image = imread("ExampleSectionImage01.png");
scaleLength = 3;

forceCentroid = false;
xCentroid = 0;
yCentroid = 0;

%==========================================================================
dimensions = size(image);
width = dimensions(2); height = dimensions(1);
numPx = width*height;

%iterate through each pixel of the image
xObjData = zeros(1,numPx);
yObjData = zeros(1,numPx);
xScaleData = zeros(1,numPx);
yScaleData = zeros(1,numPx);
bii = 1;
rii = 1;
firstLoop=true;
h = waitbar(0,'...');
for xx = 1:width
    for yy = 1:height
        %check the color of the current pixel
        localPixel = image(yy,xx,:);
        localPixel = reshape(localPixel,1,3);
        if localPixel==[255,255,255]
            %White, do nothing
        elseif localPixel==[0,0,0]
            %Black, add to object data array
            xObjData(bii) = xx;
            yObjData(bii) = yy;
            bii = bii+1;
        elseif localPixel==[255,0,0]
            %Red, add to scale data array
            xScaleData(rii) = xx;
            yScaleData(rii) = yy;
            rii = rii+1;
        end
    end
    %fprintf("%.1f %%\n",xx/width*100)
    if firstLoop
        firstLoop=false;
    end
    waitbar(xx/width,h,"Calculating")
end
delete(h)
%get rid of unused indices in the matrices
xObjData(bii:end)=[];
yObjData(bii:end)=[];
xScaleData(rii:end)=[];
yScaleData(rii:end)=[];

%find the image scale per Px and area of each pixel
scaleLengthPx = max(xScaleData)-min(xScaleData);
imageScale = scaleLength/scaleLengthPx;%in/px
pixelArea = imageScale.^2;

%Calculate the dimensions of the object
objectArea = length(xObjData)*pixelArea;
objectHeight = imageScale*(max(yObjData)-min(yObjData));
objectWidth = imageScale*(max(xObjData)-min(xObjData));
%Calculate the centroid of the object
xWeighted = (xObjData*imageScale)*pixelArea;
yWeighted = (yObjData*imageScale)*pixelArea;
if forceCentroid == false
    xCentroid = sum(xWeighted)/objectArea;%relative to image border
    xCenterPixel = floor((xCentroid/imageScale));
    yCentroid = sum(yWeighted)/objectArea;%relative to image border
    yCenterPixel = floor((yCentroid/imageScale));
    xCentroid = xCentroid-(min(xObjData)*imageScale);%relative to bottom of obj
    yCentroid = yCentroid-(min(yObjData)*imageScale);%relative to left of obj
end

xPositions = (xObjData-min(xObjData))*imageScale;
yPositions = (yObjData-min(yObjData))*imageScale;

Ixx = sum(pixelArea*((yPositions-yCentroid).^2));
Iyy = sum(pixelArea*((xPositions-xCentroid).^2));

%output data

disp("  Section Property Results:")
fprintf("Area: %.3f sq.in\nIxx: %.3f qu.in\nIyy: %.3f qu.in\n",objectArea,Ixx,Iyy)
fprintf("Width: %.3f in\nHeight: %.3f in\n",objectWidth,objectHeight)
fprintf("Centroid X: %.3f in\nCentroid Y: %.3f in\n\n",xCentroid,yCentroid)
