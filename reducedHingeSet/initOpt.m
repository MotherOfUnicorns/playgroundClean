function opt=initOpt(varargin)
%initOpt() initialize the simulations with default options. To change 
%specific options use initOpt('option1',value1,'option2',value2,...).
%
%The most important options are:
%1) 'inputType': 'preDefined', 'individual' or 'user'.
%2) 'template': any string indicating the specific template used to generate 
%   the architected material. 
%   -> When 'preDefined' is used, possible values are '#1','#2',...'#28' as 
%   defined in Fig. 3 and the supplementary information Fig. S6. Other 
%   values relate directly to the figures in the paper: Fig1a, Fig1b, 
%   Fig1c, Fig2f, Fig2g, Fig4a, Fig6#a, Fig6#b, Fig6#c, Fig6#d.
%   -> When 'individual' is used, possible values are the names of the
%   platonic solids, archimedean solids and some prisms, including: 
%   'tetrahedron', 'cube', 'octahedron', truncated tetrahedron',
%   'cuboctahedron', 'truncated cube', 'truncated octahedron',
%   'rhombicuboctahedron', 'truncated cuboctahedron', 'triangular prism',
%   'hexagonal prism', 'octagonal prism' and 'decagonal prism'.
%   -> When 'user' is used, a string needs to be provided of a script file 
%   (without the .m extension). An example is provided as 'userInputFile',
%   which also provides more information on defining your own user file.
%3) 'plot': 'info', 'result' or 'modes'. 'info' should be used when setting 
%   up a new template, while 'result' will run the mode analysis and 
%   determine the deformation modes. 'modes' is a summarized version of
%   'result'.
%
%See the intiOpt() file for a few other (less important) options. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GEOMETRY OPTIONS DEFAULT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input type that defines from which database to choose from
%Options: 'preDefined', 'inidividual' or 'user'
opt.inputType='preDefined';
%Name of the template
opt.template='#18';
%Turn periodic boundary conditions on or off
%Options: 'on' or 'off'
opt.periodic='on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT AND MOVIE OPTIONS DEFAULT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type of results to display         
%Options: 'info' or 'result'
opt.plot='result';
%Number of frames to show the modes the modes
%Integer value > 0
opt.frames=30;
%Scales the amplitude of the modes
%Value between 0 and 1: 0 no deformation, 1 some hinge angles will be zero
opt.scale=0.75;
%Save all results in a gif file
%Options: 'on' or 'off'
opt.saveMovie='off';
%When using periodic boundary condtions, number of unit cells to show
%Integer value >0
opt.plotPer=2;
%Plot options default view
opt.AZ=-54;
opt.EL=14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extrusion length
%Value >0
opt.Lextrude=1;
%Random offset for hinge stiffness
%Any value, but typically equal to zero or 1e-7
opt.perturbStiff=1e-7; %small value to seperate modes with same frequency.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UPDATE OPTIONS BASED ON INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:2:nargin
    opt.(varargin{i}) = varargin{i+1};
end
