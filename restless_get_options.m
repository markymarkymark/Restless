% -------------------------------------------------------------
function options = restless_get_options(optfile)

options = [];
if (nargin < 1), optfile = ''; end
    
% --- Set hardcoded default option values ---
if (ispc())
    defopts.scratchpath = [getenv('HOMEDRIVE') getenv('HOMEPATH') '\Documents\DATA\restless\scratch'];
    defopts.RTtoppath = 'Z:\RTexport_current'; 
else
    defopts.scratchpath = [getenv('HOME') '/data/restless/scratch'];
    if (isdir('/mnt/rtexport/RTexport_current')), defopts.RTtoppath = '/mnt/rtexport/RTexport_current';
    else                                          defopts.RTtoppath = '/mnt/rtexport'; end
end

defopts.wildcard            = '*.dcm';
defopts.empty_scratch       = 1;           % delete contents of scratch folder on start (need this for next file check to work!)
defopts.filesort            = '';          % sort list of dicoms (options are 'date', 'name', or '')
defopts.sleeptime           = 0.100;       % how long to sleep between checks for new Dicoms
defopts.maxreps             = 1000;        % program will stop after this many volumes are read
defopts.make_niftis         = 1;           % must be set for moco to work
defopts.do_moco             = 1;
defopts.prime_moco          = 1;           % prime moco with pre-acquired fmri "single-rep"
defopts.moco_primequality   = 0.1;         % if priming moco, use better quality
defopts.moco_quality        = 0.9;         % if not able to prime, then need to run w/o lengthy priming step
defopts.moco_thresh         = 0.25;        % visual level for too much motion
defopts.moco_only           = 1;           % only monitor motion, don't bother making resliced images
defopts.moco_rel            = 1;           % compute relative motion params
defopts.plot_rms            = 1;           % plot RMS metric of motion
defopts.plot_moco           = 1;
defopts.plot_ymax           = 0.4;         % initial motion plot ranges (rescales as needed)
defopts.plot_xmax           = 100;         
defopts.liveprint           = 1;
defopts.profile             = 0;

% --- restore options from file (if provided) ---
if (~isempty(optfile))
    if (exist(optfile,'file') ~= 0)
        s       = load(optfile);
        options = s.options;
    else
        fprintf(2,'ERROR: options file %s does not exist\n',optfile);
        return
    end
    
    % --- Fill in any missing options with defaults ---
    options.scratchpath         = checkstruct(options,'scratchpath',        defopts.scratchpath);
    options.RTtoppath           = checkstruct(options,'RTtoppath',          defopts.RTtoppath);
    options.wildcard            = checkstruct(options,'wildcard',           defopts.wildcard);
    options.empty_scratch       = checkstruct(options,'empty_scratch',      defopts.empty_scratch);
    options.filesort            = checkstruct(options,'filesort',           defopts.filesort);
    options.sleeptime           = checkstruct(options,'sleeptime',          defopts.sleeptime);
    options.maxreps             = checkstruct(options,'maxreps',            defopts.maxreps);
    options.make_niftis         = checkstruct(options,'make_niftis',        defopts.make_niftis);
    options.do_moco             = checkstruct(options,'do_moco',            defopts.do_moco);
    options.prime_moco          = checkstruct(options,'prime_moco',         defopts.prime_moco);
    options.moco_primequality   = checkstruct(options,'moco_primequality',  defopts.moco_primequality);
    options.moco_quality        = checkstruct(options,'moco_quality',       defopts.moco_quality);
    options.moco_thresh         = checkstruct(options,'moco_thresh',        defopts.moco_thresh);
    options.moco_only           = checkstruct(options,'moco_only',          defopts.moco_only);
    options.moco_rel            = checkstruct(options,'moco_rel',           defopts.moco_rel);
    options.plot_rms            = checkstruct(options,'plot_rms',           defopts.plot_rms);
    options.plot_moco           = checkstruct(options,'plot_moco',          defopts.plot_moco);
    options.plot_ymax           = checkstruct(options,'plot_ymax',          defopts.plot_ymax);
    options.plot_xmax           = checkstruct(options,'plot_xmax',          defopts.plot_xmax);
    options.liveprint           = checkstruct(options,'liveprint',          defopts.liveprint);
    options.profile             = checkstruct(options,'profile',            defopts.profile);
   
% --- no options file give, use straight defaults ---
else
    options = defopts;
end
return