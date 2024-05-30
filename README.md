What is this?

This is a collection of scripts and resources for analyzing composite structures. The equations and methods used are from Dr. Coburn's Composite Strength Handbook 2nd edition (abbreviated CSH throughout).
The purpose of these scripts are to allow people who haven't taken the composites class to do some basic analysis of parts with known layups and loadings. I've summarized each script or file below:

materialProperties.xlsx

This contains all of the material data for different composite materials. It is referenced by each MATLAB script as well. A description of the different data points and their associated units is located in materialPropertiesInfo.txt. When testing new materials to get allowables, you don't need every one of these points. The important ones are: E1, E2, G12, F1t, F1c, F2t, F2c, and v12.

CLTbasic.m

This is a script that uses "Classical Lamination Theory" to find the strains and stresses of each lamina in the laminate. It also outputs margins of safety for each lamina as well as a non-numerical visualization of the different couplings that are present in the laminate. Various combinations of lamina orientations result in coupling between axial loads, torsion, bending and shear. Here are some examples:
[0,90,0] does not contain coupling.
[0,45,0] contains shear/axial coupling.
[0,0,90] contains bending/axial coupling.
[-45,45,90] contains bending/axial, axial/torsion, and some torsion/bending coupling.
As an input, this script takes the layup schedule (orientation of each lamina), the thickness of each lamina, the material of each lamina, and the loading on the laminate.

FlatPanelStability.m

This is a script that is used to find the buckling allowable of a flat panel. It works for panels in uniaxial compression and/or shear and with the unloaded edges either both clamped or both simply supported. Depending on the panel dimensions given, you may have to manually select a Kc and Ks value from figures 8.2.6 and 8.2.7. When prompted, click the y value on the plot corresponding to your edge constraint (see the edge cases on the figure, find the line associated with either case A or C, and click where it intersects with the red dotted line). You also may be prompted to identify the edge cases by typing "SS" or "C" for simply supported and clamped respectively. Note: The reason that edge cases B, D, and E are not supported in this script is because there is no line given for determining Ks. You can calculate these cases by hand for uniaxial compression alone relatively easily by following 8-3 thru 8-5 and using the D (stiffness parameter) value from this script.
