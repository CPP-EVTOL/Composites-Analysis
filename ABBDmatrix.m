%Alexander Walsh 5/29/2024
clear
format compact

%This is a classical lamination thory program that determines the margins
%of safety for each layer of a laminate. It does NOT consider stability,
%only stress. Crippling and buckling must be analyzed separately.
%Define your laminate and loading below as shown in the following example: 

%{
    %define the layup schedule (degrees) from top to bottom
    layup = [0,45,90,-45,0];
    %define layer thicknesses (in)
    thickness = [0.1,0.1,0.1,0.1,0.1
    %define the material for each layer (index in the materialProperties.xlsx file)
    material = [1,1,1,1,1];
    %define the loading for the laminate [Nx,Ny,Nxy,Mx,My,Mxy] where runing
    loads are in lb/in and running moments are in lb. See CSH fig. 4.4.2
    load = [0,0,0,0,0,0];
%}

%============================== USER INPUT ================================
%NOTE: Layers are numbered from the top down

layup = [45,-45,0,-45,45];
thickness = [0.1,0.1,0.1,0.1,0.1];
material = [3,3,3,3,3];
load = [100000,0,0,0,0,0];

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

%display which couplings are present (arbitrary values, just to get an idea)
c1 = mean(abs([A(1,3),A(2,3)]));
c2 = mean(abs([B(1,1),B(1,2),B(2,2)]));
c3 = mean(abs([B(1,3),B(2,3),B(3,3)]));
c4 = mean(abs([D(1,3),D(2,3)]));
cmax = max([c1,c2,c3,c4]);
c1 = c1/cmax; c2 = c2/cmax; c3 = c3/cmax; c4 = c4/cmax;
disp("  Coupling present in the laminate:")
fprintf("Shear/Axial:     "); printSlash(c1);
fprintf("Bending/Axial:   "); printSlash(c2);
fprintf("Axial/Torsion:   "); printSlash(c3);
fprintf("Torsion/Bending: "); printSlash(c4);

%calculate strains at the midplane
midplaneStrains = (ABBD^-1)*(load');
%calculate the strains and stresses for each layer
layerStrainsBody = zeros(3,nLayers);
layerStressesBody = zeros(3,nLayers);
layerStressesPrincipal = zeros(3,nLayers);
layerStrainsPrincipal = zeros(3,nLayers);
for ii = 1:nLayers
    layerStrainsBody(:,ii) = midplaneStrains(1:3)+zbark(ii)*midplaneStrains(4:6);%strains of each layer in the body axis (x and y directions)
    layerStressesBody(:,ii) = Qbar(:,:,ii)*layerStrainsBody(:,ii);%stresses in each layer in the body axis (x and y directions)
    layerStressesPrincipal(:,ii) = T(:,:,ii)*layerStressesBody(:,ii);%stresses in each layer in the principal axes (1 and 2 direction)
    layerStrainsPrincipal(:,ii) = S(:,:,ii)*layerStressesPrincipal(:,ii);%strains in each layer in the principal axes (1 and 2 direction)
end

%get MS's for each layer using the Hoffman Failure Criterion
MS = zeros(1,nLayers);
for ii = 1:nLayers
    %applied load
    f1 = layerStressesPrincipal(1,ii);
    f2 = layerStressesPrincipal(2,ii);
    f12 = layerStressesPrincipal(3,ii);
    %allowable (with adjusted sign, see CSH 2-20)
    F1t = materialData.F1t(material(ii))*sign(f1);
    F1c = materialData.F1c(material(ii))*sign(f1);
    F2t = materialData.F2t(material(ii))*sign(f2);
    F2c = materialData.F2c(material(ii))*sign(f2);
    F12 = materialData.F12(material(ii));
    %calculate MS
    term1 = -((f1^2)/(F1c*F1t));
    term2 = (f1*f2)/(F1c*F1t);
    term3 = -((f2^2)/(F2c*F2t));
    term4 = f1*(F1c+F1t)/(F1c*F1t);
    term5 = f2*(F2c+F2t)/(F2c*F2t);
    term6 = (f12^2)/(F12^2);
    MS(ii) = (1/(term1+term2+term3+term4+term5+term6))-1;
end

%print margins of safety
disp("  MS for each layer:")
for ii = 1:nLayers
    fprintf("Layer %.0f: %.3f\n",ii,MS(ii))
end

%accepts a number 0-1, prints 0 to 10 slashes followed by spaces inside []
function printSlash(input)
    n = round(10*input);
    fprintf("[");
    for ii = 1:n
        fprintf("/");
    end
    for ii = n:10
        fprintf(" ");
    end
    fprintf("]\n");
end