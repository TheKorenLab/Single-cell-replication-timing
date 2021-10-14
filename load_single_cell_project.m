topdir = '~/projects/single_cell_replication_timing/';

addpath([topdir '/processing'])
addpath([topdir '/analysis'])
addpath([topdir '/visualization'])

% Load sample list
samples = {'GM12878', 'GM12891', 'GM12892', 'H1', 'H7', 'H9', 'HCT116', 'RKO', 'MCF7', 'GM18507'};

% Set global graphics parameters
set(groot, 'defaultFigureUnits', 'inches', 'defaultFigureColor', 'white') % figure
set(groot, 'defaultAxesFontName', 'Arial', 'defaultAxesFontSize', 9, ...
    'defaultAxesLabelFontSizeMultiplier', 1.1, ...
    'defaultAxesBox', 'on', 'defaultAxesTickLength', [0 0], 'defaultAxesNextPlot', 'add', ...
    'defaultAxesToolbarVisible', 'off') % axes
set(groot, 'defaultLegendFontSize', 11, 'defaultLegendFontSizeMode', 'manual') % legend

% Color scheme
g1_dark = '#969696';
g1_light = '#DCDCDC';
s_dark = '#00441B';
s_light = {'#58BC91', '#4EB3D3', '#FA9FB5'};
