%Alexander Walsh 5/29/2024
clear
format compact

%This is a program to determine the stability allowables of a flat plate
%made of composite material loaded under uniaxial compression and/or shear.
%If your laminate is specially orthotropic, use the equations in the CSH
%section 8-6 instead for a more accurate result.
%Define your laminate and loading below as shown in the following example: 

%{
    %define the layup schedule (degrees) from top to bottom
    layup = [0,45,90,-45,0];
    %define layer thicknesses (in)
    thickness = [0.1,0.1,0.1,0.1,0.1];
    %define the material for each layer (index in the materialProperties.xlsx file)
    material = [1,1,1,1,1]; (Uni Glass/Epoxy=1, Uni Boron/Epoxy=2, etc.)
    %panel dimensions a = unloaded side, b = loaded side (both in inches)
    a = 8;
    b = 4;
    %loading Nc and Ns are the compressive and shear running loads, lb/in
    Nc = 1;
    Ns = 0;
%}


%============================== USER INPUT ================================

layup = [0,0,0,0,0];
thickness = [0.02,0.02,0.02,0.02,0.02];
material = [3,3,3,3,3];
a = 28;
b = 6;
Nc = 10;
Ns = 200;

%==========================================================================

if (length(thickness)~=length(material))||(length(layup)~=length(material))||(length(layup)~=length(thickness))
    error("Check your inputs. The layup, thickness, material, and load matrices must be the same length");
end

if (a/b<1)
    error("a/b is not within a valid range (it's less than 0.5)")
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

%Get the Ks and Kc values for various edge length ratio cases
if ((a/b)>=5)&&((a/b)<25)
    edges = input("Are the unloaded edges clamped or simply supported? Respond with either C or SS\n",'s');
    if edges == "SS"
        Ks = 5.35+4*(b/a)^2;
        Kc = 4;
    else
        Ks = 9.5;
        Kc = 6.98;
    end
end

if ((a/b)>=25)
    edges = input("Are the unloaded edges clamped or simply supported? Respond with either C or SS\n",'s');
    if edges == "SS"
        Ks = 5.35;
        Kc = 4;
    else
        Ks = 692.331;
        Kc = 0;
    end
end

if ((a/b)<5)
    %NOTE: The original jpg images must be used. If they're changed then this
    %section of the code will not be calibrated and the selected K values will
    %be completely wrong.
    im1 = imread("CSH fig 8-2-6.jpg");
    im2 = imread("CSH fig 8-2-7.jpg");
    
    xCoef = (a/b)/5;
    topPx = xCoef*(2.1095e3-118.4465)+118.4465;
    bottomPx = xCoef*(2.1295e3-94.4574)+94.4574;
    topPy = xCoef*(78.2241-86.2205)+86.2205;
    bottomPy = xCoef*(3.2328e3-3.2608e3)+3.2608e3;
    %get Kc from CSH fig. 8.2.6
    fig = figure;
    imshow(im1);
    title("Select the y-value along the red line associated with your edge constraints (case A or C)")
    xlabel("CSH figure 8.2.6")
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-3.2608e3)/(86.2205-3.2608e3);
    yCoef2 = (y-3.2328e3)/(78.2241-3.2328e3);
    Kc = mean([yCoef1,yCoef2])*16;
    close(fig)
    %get Ks from fig. 8.2.7
    topPx = xCoef*(2.0835e3-112.5)+112.5;
    bottomPx = xCoef*(2.1135e3-148.5)+148.5;
    topPy = xCoef*(80-80)+80;
    bottomPy = xCoef*(2009-2.045e3)+2.045e3;
    fig = figure;
    imshow(im2);
    title("Select the y-value alond the red line associated with your edge constraints")
    xlabel("CSH figure 8.2.7")
    hold on
    plot([bottomPx,topPx],[bottomPy,topPy],'r--','LineWidth',2)
    [~,y] = ginput(1);
    yCoef1 = (y-2.045e3)/(80-2.045e3);
    yCoef2 = (y-2009)/(80-2009);
    Ks = mean([yCoef1,yCoef2])*10+5;
    close(fig)
end

Ncr_c = (Kc*pi*pi*D)/(b^2);%compression allowable running load
Ncr_s = (Ks*pi*pi*D)/(b^2);%shear allowable running load
Rc = Nc/Ncr_c;
Rs = Ns/Ncr_s;
MS = (2/(Rc+sqrt((Rc^2)+4*(Rs^2))))-1;

disp("  Flat Panel results:")
fprintf("Compression allowable: %.3f lb/in\nShear allowable: %.3flb/in\nMS: %.3f\n\n",Ncr_c,Ncr_s,MS)