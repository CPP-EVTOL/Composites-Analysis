%Alexander Walsh 5/30/2024
clear
format compact

%This is a program to determine the stability allowables of a cylinder
%loaded under uniaxial compression, torsion, and bending.
%Define your laminate and loading below as shown in the following example:

%{
    %define the layup schedule (degrees) from top to bottom
    layup = [0,45,90,-45,0];
    %define layer thicknesses (in)
    thickness = [0.1,0.1,0.1,0.1,0.1];
    %define the material for each layer (index in the materialProperties.xlsx file)
    material = [1,1,1,1,1];
    %cylinder dimensions L = length, r = laminate midplane radius (both in inches)
    L = 36
    r = 5;
    %loading P = compression load (lb), M = moment (in*lb), Nxy
    = torsion running load (lb/in);
    P = 1;
    M = 10;
    Nxy = 1;
%}


%============================== USER INPUT ================================

layup = [0,90];
thickness = [0.02,0.02];
material = [1,1];
L = 36;
r = 8;
P = 220;
M = 0;
Nxy = 0;

%==========================================================================

if (length(thickness)~=length(material))||(length(layup)~=length(material))||(length(layup)~=length(thickness))
    error("Check your inputs. The layup, thickness, material, and load matrices must be the same length");
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
Z = ((L^2)/(r*tlam))*sqrt(1-(veff^2));

if(Z>(10^5))
    error("Geometric parameter Z is too large. This could be caused by: Cylinder length too large, radius too small, thickness too small, or poisson's ratio too large.")
    %If this error occurs, don't just delete this line and try again
    %anyway. The plots in section 8-18,19,21 do not have Kc or Kt values for
    %Z above 10^5 and something is wrong with your geometry.
end

if(Z<1)
    error("Geometric parameter Z is too small. This could be caused by: Cylinder length too small, radius too large, thickness too large, or poisson's ratio too small.")
    %If this error occurs, don't just delete this line and try again
    %anyway. The plots in section 8-18,19,21 do not have Kc or Kt values for
    %Z below 1 and something is wrong with your geometry.
end

if((r/tlam)<100)
    disp("CAUTION: r/t is too small for accurate compression results to be calculated.")
    Kc = 1;
end

%Get the Kc and Kt values for various r/t values and edge constraints
%NOTE: The original jpg images must be used. If they're changed then this
%section of the code will not be calibrated and the selected K values will
%be completely wrong.
im1 = imread("CSH fig 8-4-1.jpg");
im2 = imread("CSH fig 8-4-2.jpg");
im3 = imread("CSH fig 8-4-3.jpg");
im4 = imread("CSH fig 8-4-4.jpg");
im5 = imread("CSH fig 8-4-6.jpg");

%Kc
xCoef = log10(Z)/5;
titleString = "Select the y-value along the red line associated with your edge constraint (r/t = "+string(r/tlam)+")";
if((r/tlam)>=100)&&((r/tlam)<=500)
    topLeft = [248,274];
    topRight = [2099,298];
    bottomLeft = [212,1672];
    bottomRight = [2090,1690];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im1);
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Kc = 10000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_c = (Kc*pi*pi*D)/(L^2);
end

if((r/tlam)>500)&&((r/tlam)<=1000)
    topLeft = [271.45,286.85];
    topRight = [2193.5,337.83];
    bottomLeft = [211.48,1813.1];
    bottomRight = [2184.5,1849.1];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im2);
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Kc = 10000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_c = (Kc*pi*pi*D)/(L^2);
end

if((r/tlam)>1000)&&((r/tlam)<=2000)
    topLeft = [265,302];
    topRight = [2062,314];
    bottomLeft = [259,1781];
    bottomRight = [2086,1769];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im3);
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Kc = 10000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_c = (Kc*pi*pi*D)/(L^2);
end

if((r/tlam)>2000)
    topLeft = [301.87,316.79];
    topRight = [2344,331.78];
    bottomLeft = [283.88,1951.1];
    bottomRight = [2341,1918];
    topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
    bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
    topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
    bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

    fig = figure;
    imshow(im4);
    title(titleString);
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
    yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
    Kc = 10000^(mean([yCoef1,yCoef2]));
    close(fig)

    Ncr_c = (Kc*pi*pi*D)/(L^2);
end

%kt (torsion)
topLeft = [191.5,65];
topRight = [2.6255e3,89];
bottomLeft = [171.5,1585];
bottomRight = [2619.5,1579];
topPx = xCoef*(topRight(1)-topLeft(1))+topLeft(1);
bottomPx = xCoef*(bottomRight(1)-bottomLeft(1))+bottomLeft(1);
topPy = xCoef*(topRight(2)-topLeft(2))+topLeft(2);
bottomPy = xCoef*(bottomRight(2)-bottomLeft(2))+bottomLeft(2);

fig = figure;
imshow(im5);
titleString = "Select the y-value along the red line associated with your edge constraint";
title(titleString);
hold on
plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
[~,y] = ginput(1);
yCoef1 = (y-bottomLeft(2))/(topLeft(2)-bottomLeft(2));
yCoef2 = (y-bottomRight(2))/(topRight(2)-bottomRight(2));
Kt = 1000^(mean([yCoef1,yCoef2]));
close(fig)

Ncr_t = (Kt*pi*pi*D)/(L^2);

%get maximum running load from compression and moment
Nmax = (P/(2*pi*r))+(M/(pi*r*r));

Rc = Nmax/Ncr_c;
Rt = Nxy/Ncr_t;

MS = (2/(Rc+sqrt((Rc^2)+4*(Rt^2))))-1;

disp("  Cylinder Results:")
fprintf("Max compression runnung load: %.3f lb/in\nCompression Allowable: %.3f lb/in\nTorsion Allowable: %.3f lb/in\nMS: %.3f\n\n",Nmax,Ncr_c,Ncr_t,MS)