topdir = '~/projects/single_cell_replication_timing/';

addpath([topdir '/processing'])
addpath([topdir '/analysis'])
addpath([topdir '/visualization'])

% Load sample list
samples = {'GM12878', 'GM12891', 'GM12892', 'H1', 'H7', 'H9', 'HCT116', 'RKO', 'MCF7', ...
    'GM18507', 'COLO829'};

cell_line_names = samples;
cell_line_names{7} = 'HCT-116';
cell_line_names{9} = 'MCF-7';
cell_line_names{14} = 'COLO-829';

% Set global graphic parameters
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigureColor', 'white') % figure
set(groot, 'defaultAxesUnits', 'centimeters', 'defaultAxesFontSizeMode', 'manual', ...
    'defaultAxesFontName', 'Arial', ...
    'defaultAxesFontSize', 6.75, ...
    'defaultAxesLabelFontSizeMultiplier', 1.037, 'defaultAxesTitleFontSizeMultiplier', 1.037, ...
    'defaultAxesBox', 'on', 'defaultAxesTickLength', [0 0], 'defaultAxesNextPlot', 'add', ...
    'defaultAxesToolbarVisible', 'off') % axes
set(groot, 'defaultLegendFontSizeMode', 'manual', 'defaultLegendFontSize', 7, ...
    'defaultLegendBox', 'on') % legend

% Color scheme
g1_dark = '#969696';
g1_light = '#DCDCDC';
s_dark = '#00441B';
s_light = {'#58BC91', '#4EB3D3', '#FA9FB5'};
