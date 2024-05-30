%Alexander Walsh 5/30/2024
clear
format compact

%This is a program to determine the stability allowables of a curved plate
%made of composite material loaded under uniaxial compression and/or shear.
%Define your laminate and loading below as shown in the following example:

%{
    %define the layup schedule (degrees) from top to bottom
    layup = [0,45,90,-45,0];
    %define layer thicknesses (in)
    thickness = [0.1,0.1,0.1,0.1,0.1];
    %define the material for each layer (index in the materialProperties.xlsx file)
    material = [1,1,1,1,1];
    %panel dimensions a = unloaded side, b = loaded side arc length (both in inches, always do it this way even if a<b).
    a = 8;
    b = 4;
    %loading Nc and Ns are the compressive and shear running loads, lb/in
    Nc = 1;
    Ns = 0;
    %define the edge constraints "C" = clamped, "SS" = simply supported
    edges = "SS"
    %define the panel radius (inches)
    r = 10;
    
%}


%============================== USER INPUT ================================

layup = [45,-45,0,90,0,90,0,-45,45];
thickness = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01];
material = [1,1,1,1,1,1,1,1,1];
a = 36;
b = 14;
Nc = 60;
Ns = 10;
edges = "SS";
r = 5;

%==========================================================================

if (length(thickness)~=length(material))||(length(layup)~=length(material))||(length(layup)~=length(thickness))
    error("Check your inputs. The layup, thickness, material, and load matrices must be the same length");
end

if (edges~="SS")&&(edges~="C")
    error("Invalid edge constraints. Edges must either be clamped 'C' or simply supported 'S'");
end

materialData = readtable("materialProperties.xlsx");
nLayers = length(layup);

R = [1,0,0;0,1,0;0,0,2];
%generate T, T inverse, Q, S matrices for each layer
T = zeros(3,3,nLayers);
Tinv = zeros(3,3,nLayers);
Q = zeros(3,3,nLayers);
S = zeros(3,3,nLayers);
Qbar = zeros(3,3,nLayers);
layup = deg2rad(layup);
%calculate the qbar matrix (see CSH section 2-14)
for ii = 1:nLayers
    theta = layup(ii);
    E1 = materialData.E1(material(ii));
    E2 = materialData.E2(material(ii));
    G12 = materialData.G12(material(ii));
    v12 = materialData.v12(material(ii));
    v21 = (E2/E1)*v12;
    T(:,:,ii) = [cos(theta)^2,sin(theta)^2,2*sin(theta)*cos(theta);...
                sin(theta)^2,cos(theta)^2,-2*sin(theta)*cos(theta);...
                -1*sin(theta)*cos(theta),sin(theta)*cos(theta),(cos(theta)^2)-(sin(theta)^2)];
    Tinv(:,:,ii) = inv(T(:,:,ii));
    Q(:,:,ii) = [(E1/(1-v12*v21)),(v12*E2)/(1-v12*v21),0;...
                (v12*E2)/(1-v12*v21),(E2/(1-v12*v21)),0;...
                0,0,G12];
    Qbar(:,:,ii) = (T(:,:,ii)^-1)*Q(:,:,ii)*R*T(:,:,ii)*(R^-1);
    S(:,:,ii) = [1/E1,-v12/E1,0;-v12/E1,1/E2,0;0,0,1/G12];
end

%get the midpoint of each layer
tlam = sum(thickness);
zmid = tlam/2;
zbottom = cumsum(thickness);
ztop = zbottom-thickness;
zbark = mean([zbottom;ztop]);%relative to upper surface
zbark = zbark-zmid;%relative to the midplane of the laminate
%create the A, B, and D matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
for ii = 1:nLayers
    tk = thickness(ii);
    A = A+Qbar(:,:,ii)*tk;
    B = B+Qbar(:,:,ii)*tk*zbark(ii);
    D = D+Qbar(:,:,ii)*(tk*(zbark(ii)^2)+((tk^3)/(12)));
end
ABBD = [A,B;B,D];
D = mean([D(1,1),D(1,2),D(2,2),D(3,3)]);%stiffness parameter

%find the effective poisson's ratio for the whole panel
veff = 0;
for ii = 1:nLayers
    veff = veff+materialData.v12(material(ii));
end
veff = veff/nLayers;

%find geometric parameter Z
Z = ((b^2)/(r*tlam))*sqrt(1-(veff^2));

if(Z>1000)
    error("Geometric parameter Z is too large. This could be caused by: Panel dimension b too large, radius too small, thickness too small, or poisson's ratio too large.")
    %If this error occurs, don't just delete this line and try again
    %anyway. The plots in section 8-9,10,11 do not have Kc or Ks values for
    %Z above 1000 and something is wrong with your geometry.
end
%Get the Ks and Kc values for various edge length ratio cases and
%constraints
%NOTE: The original jpg images must be used. If they're changed then this
%section of the code will not be calibrated and the selected K values will
%be completely wrong.
%Kc
im1 = imread("CSH fig 8-3-1.jpg");
im2 = imread("CSH fig 8-3-2.jpg");
im3 = imread("CSH fig 8-3-3.jpg");
im4 = imread("CSH fig 8-3-4.jpg");
im5 = imread("CSH fig 8-3-5.jpg");

xCoef = log10(Z)/4;
topPx = xCoef*(2.687e3-269)+269;
bottomPx = xCoef*(2.708e3-260)+260;
topPy = xCoef*(107-86)+86;
bottomPy = xCoef*(1.886e3-1.889e3)+1.889e3;

fig = figure;
imshow(im1);
titleString = "Select the y-value along the red line associated with your r/t and ("+string(edges)+") edge constraints (Your r/t ="+string(r/tlam)+")";
title(titleString);
hold on
plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
[~,y] = ginput(1);
yCoef1 = (y-1.889e3)/(86-1.889e3);
yCoef2 = (y-1.886e3)/(107-1.886e3);
Kc = 1000^(mean([yCoef1,yCoef2]));
close(fig);

Ncr_c = (Kc*pi*pi*D)/(b^2);

%Ks
xCoef = log10(Z)/3;
if(a>=b)&&(edges=="SS")
    topLeft = [475,80];
    topRight = [2.2335e3,86];
    bottomLeft = [427.5,1.883e3];
    bottomRight = [2.2065e3,1.913e3];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im2);
    titleString = "Select the y-value along the red line associated with your a/b. (Your a/b ="+string(a/b)+")";
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Ks = 1000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_s = (Ks*pi*pi*D)/(b^2);
end

%a and b switched here
if(a<b)&&(edges=="SS")
    topLeft = [467,43.5];
    topRight = [2159,130.5];
    bottomLeft = [419,1.864e3];
    bottomRight = [2150,1.8735e3];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im3);
    titleString = "Select the y-value along the red line associated with your a/b. (Your a/b ="+string(b/a)+")";
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Ks = 1000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_s = (Ks*pi*pi*D)/(a^2);
end

if(a>=b)&&(edges=="C")
    topLeft = [242,103.5];
    topRight = [2114,127.5];
    bottomLeft = [233,1.9215e3];
    bottomRight = [2096,1935.5];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im4);
    titleString = "Select the y-value along the red line associated with your a/b. (Your a/b ="+string(a/b)+")";
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Ks = 1000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_s = (Ks*pi*pi*D)/(b^2);
end

%a and b switched here
if(a<b)&&(edges=="C")
    topLeft = [277.758,94.8824];
    topRight = [2178.1,97.8797];
    bottomLeft = [280.7556,1938.3];
    bottomRight = [2172.1,1917.3];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im5);
    titleString = "Select the y-value along the red line associated with your a/b. (Your a/b ="+string(b/a)+")";
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Ks = 1000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_s = (Ks*pi*pi*D)/(a^2);
end

Rc = Nc/Ncr_c;
Rs = Ns/Ncr_s;

MS = (2/(Rc+sqrt((Rc^2)+4*(Rs^2))))-1;

disp("  Curved Panel Results:")
fprintf("Compression Allowable: %.3f\nShear Allowable: %.3f\nMS: %.3f\n",Ncr_c,Ncr_s,MS)