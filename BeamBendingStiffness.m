%Alexander Walsh 5/30/2024
clear
format compact

%This is a script for getting the equivalent "EI" for a composite laminate.
%The "EI" parameter is used in the beam bending canned formulas at the back
%of the structures books (APPX.D) to find deflection/slope of beams with  
%various loadings and supports. If you are looking for the stresses and
%strains associated with a bending moment, use CLTbasic.m along with the  
%moment calculated in Appendix. D instead. This script also calculates the
%moment and deflection for the most common beam loading case (canteliever
%beam with a downward point load at one end) using case 1a.
%NOTE: It is only accurate for symmetrical laminates.

%Define your laminate and loading below as shown in the following example: 

%{
    %define the layup schedule (degrees) from top to bottom
    layup = [0,45,90,-45,0];
    %define layer thicknesses (in)
    thickness = [0.1,0.1,0.1,0.1,0.1];
    %define the material for each layer (index in the materialProperties.xlsx file)
    material = [1,1,1,1,1];
    %define the width of the beam (inches)
    b = 3;
    %define point load (lb)
    P = 20;
    %define beam length (inches)
    L = 10;
%}

%============================== USER INPUT ================================
%NOTE: Layers are numbered from the top down

layup = [45,-45,0,90,0,90,0,-45,45];
thickness = [0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02];
material = [1,1,1,1,1,1,1,1,1];
b = 3;
P = 20;
L = 10;

%==========================================================================

%Yes, the code here is mostly identical to CLTbasic.m and the other CLT
%based programs. Finding the ABBD matrix is the first step in any laminate
%calculation.

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

EIeq = b*D(1,1);

Mx = P*L;
ymax = -(P*(L^3))/(3*EIeq);
thetamax = -(P*(L^2))/(2*EIeq);

disp("  Beam Stiffness Results:");
fprintf("Equivalent EI: %.3f lb*sq.in\nMx: %.3f in*lb\nMax deflection: %.3f in\nMax angle: %.3f deg\n\n",EIeq,Mx,ymax,rad2deg(thetamax));