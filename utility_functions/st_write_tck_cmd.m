function cmd = st_write_tck_cmd(fact_input,tck_out,cfg)
% fact_input: an image file (nifti, mif etc) with 3XN volumes. Each volume
% holds a single peak orientation to be used for tractography. See MRTrix
% 3.0 documentation for more details.
% tck_out: Output file name.
% 
% Optional arguments are given as cfg fields:
% cfg.fg: a fiber group, or a fiber group tck file
% cfg.plane: 0 sagittal; 1: coronal; 2: axial.
% cfg.mode: 1,2,3
% etc.
% See "mrview --help" for all options
% Example:
%     cfg = [];
%     cfg.maxlength = '5000';
%     cfg.step = '0.2';
%     cfg.angle = '60';
%     cfg.seed_grid_per_voxel = [seed_img ' 5']
%     cfg.select = '10000';
%     cfg.mask = mask_img;

%% Create the command
cmd = ['tckgen -algorithm FACT '];

fNames = fields(cfg);


if isfield(cfg,'include')
    for rI = 1:length(cfg.include)
        if any(strcmp(fNames,'include'))
                cmd = [cmd, '-include ', cfg.include{rI}];
        end
    end
    fNames(strcmp(fNames,'include')) = [];   
end

for fI = 1:length(fNames)
    cmd = [cmd, ' -', fNames{fI} ' ' cfg.(fNames{fI})];
end

% Only grab the image at the end of the process
% if ~isfield(cfg.tractography,'tsf_load') % For per-streamline color, don't capture the image until the user changes the colormap manually
cmd = [cmd ' ' fact_input ' ' tck_out ' -force'];