%Alexander Walsh 6/3/2024
clear
format compact

%This script is for preliminary sizing for a laminate
%under either uniaxial load, biaxial loading of the same sign, or biaxial
%loading of opposite signs (shear). The percentages of 0, 90, and 45 degree
%plies are calculated using John Smith's 10% rule iteratively in reverse.
%See CSH section 5-7 for more information on the 10% rule.

%{
                                HOW TO USE:

1. Input the loading on your panel as a running load in the "loads"
variable. It should be a 1x2 matrix with [Nx, Ny]. If you have a shear
load, it must be converted to a biaxial load of opposite signs first.

2. Input the thickness of each lamina (layer) in the "thickness" variable
in inches.

3. Input the material as an index of the materialProperties.xlsx
spreadsheet into the "material" variable.


%}

%============================== USER INPUT ================================

loads = [-23,20];
thickness = 0.02;
material = 1;

%==========================================================================

if((loads(2)==0)||(loads(1)==0))
    loadCase = 1;
elseif(sign(loads(1))~=sign(loads(2)))
    loadCase = 3;
elseif(sign(loads(1))==sign(loads(2)))
    loadCase = 2;
else
    error("Invalid load case.")
end

materialData = readtable("materialProperties.xlsx");


switch loadCase
    case 1
        %
    case 2
        %
    case 3
        %
    otherwise
        %
end