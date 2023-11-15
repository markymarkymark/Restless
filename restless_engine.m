function status = restless_engine(mode, options)
% --------------------------------------------------------------------
% Display motion in realtime from Siemens RT Export folder
%
% MElliott 3/2015
% --------------------------------------------------------------------

% --- static global variable ---
global restless_plotzoom

% --- static local vars ---
persistent TR timeout_time count RTpath figh lastfile nwait vol0 tstart last_mat center
persistent mocopars rmsmotion plotpars tnow nwaits filenames
persistent xrange yrange

status    = 0;                 % error
stopcolor = [0.8 0.8 0.8];
runcolor  = [1 0.6 0.6];

switch (mode)
    % --------------------------------------------
    % --- MODE = 0 Initialize (called once)
    % --------------------------------------------
    case 0
        
        % --- check paths ---
        if (~isdir(options.RTtoppath)), fprintf(2,'ERROR: RT Export folder %s does not exist!\n',options.RTtoppath); return; end
        if (~isdir(options.scratchpath)), fprintf(2,'ERROR: Scratch folder %s does not exist!\n',options.scratchpath); return; end
        fprintf(1,'\n\n---------------------------------------------------------------\n');
        fprintf('Finding most recent RT Export subfolder (may take a few moments)...\n');
        subjfolder = get_next_file(options.RTtoppath,'*','','date');
        RTpath = [options.RTtoppath filesep() subjfolder];
        fprintf(1,'Using RT export folder: %s\n',RTpath);
        
        % --- empty scratch folder ---
        if (options.empty_scratch)
            recycle('off');
            delete([options.scratchpath filesep() '*.nii']);
            delete([options.scratchpath filesep() '*mcpar.txt']);
        end
        
        % --- arrays to track timings, etc, ... ---
        mocopars  = zeros(options.maxreps,6);
        rmsmotion = zeros(options.maxreps);
        plotpars  = zeros(options.maxreps,4);
        tnow      = zeros(options.maxreps,1);
        nwaits    = zeros(options.maxreps,1);
        filenames = cell(options.maxreps,1);
        last_mat  = eye(4); % identity matrix
        
        % --- set up for moco ---
        if (options.do_moco)
            vol0 = []; % indicate that moco priming must still be done
        end
        
        % --- set up plot window for displaying moco results ---
        if (options.do_moco && options.plot_moco)
            set(0,'units','pixels');
            screen = get(0,'screensize');
            %mp = get(0, 'MonitorPositions');    % gets sizes of all displays
            %ndisplays = size(mp,1);
            %screen = mp(ndisplays,:);           % get size of last display
            width = round(0.75*screen(3));
            height = round(0.50*screen(4));
            x0 = screen(3)/2 - width/2;         % centered left-right
            y0 = screen(4) - height - 50;       % at top with room for border
            figh=figure(1);
            clf(figh);
            set(figh, 'menubar', 'none');
            callstr = 'set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))) ; uiresume ' ;
            set(figh,'position',[x0 y0 width height],'name','RestLess Plot','keypressfcn',callstr,'windowstyle','normal','numbertitle','off','userdata',-1);
            set(figh,'color',stopcolor);
            
            hold off
            plot(1:options.maxreps,ones(options.maxreps,1)*options.moco_thresh,'--m')
            restless_plotzoom = 1;
            xrange = [1 options.plot_xmax];
            yrange = [0 options.plot_ymax];
            axis([xrange yrange/restless_plotzoom]);
            title('Real Time Motion Params');
            xlabel('measurement #')
            if (options.moco_rel), ylab = 'relative ';
            else                   ylab = ''; end
            if (options.plot_rms), ylabel([ylab 'RMS disp (mm)']);
            else                   ylabel([ylab '|translations| (mm)']); end
            drawnow
        end
        
        % --- Get started, return file with latest timestamp ---
        [filename,stat] = get_next_file(RTpath,options.wildcard,'',options.filesort);
        if (~stat), return; end
        
        % --- use latest Dicom to prime moco system ---
        if (options.do_moco && options.prime_moco)
            errmsg  = '';
            if (isempty(filename))
                errmsg = 'No reference dicom file found.';
            else
                infile    = [RTpath filesep() filename];
                [~,froot] = fileparts(infile);
                outfile   = [options.scratchpath filesep() froot '.nii'];
                hdr       = ME_spm8_dicom_headers4(infile);
                if (isempty(hdr))
                    errmsg = 'Could not read dicom file.';
                elseif (isempty(strfind(hdr{1}.ImageType,'MOSAIC')))
                    errmsg = 'Dicom is not from a mosaiced EPI scan.';
                else
                    tmp = ME_spm8_dicom_convert4(hdr,outfile);
                    if (isempty(tmp))
                        errmsg = 'Error converting dicom to nifti.';
                    else
                        fprintf(1,'Priming moco system with %s (quality = %g).\n',filename,options.moco_primequality);
                        vol0   = spm_vol(outfile);
                        nifti_coreg(vol0,[],'rigidbody_fast','',options.moco_primequality);
                        center = nifti_get_center(vol0); % need for rmsmotion calc
                    end
                end
            end
            if (~isempty(errmsg))
                vol0 = []; 
                fprintf(2,'\nWARNING: %s Could not prime moco system in advance of run.\n',errmsg);
            end
        end
        
        % --- Initialize variables for repeated calls with MODE = 1 ---
        lastfile = filename;
        count   = 1;
        nwait   = 0;
        timeout_time = 300;                     % wait time for first file to arrive (sec)
        fprintf(1,'Ready and waiting for scan to begin...\n');
        if (options.do_moco && options.plot_moco), set(figh,'color',runcolor); end
        drawnow
        status = 1; % success
        
        % -----------------------------------------------------------------------
        % --- MODE = 1 repeatedly called to check for new dicoms and update moco
        % -----------------------------------------------------------------------
    case 1
        
        % --- Quit if max # of files have been found ---
        if (count > options.maxreps)
            fprintf(1,'\n\nMaxreps = %1d reached. Done.\n',options.maxreps);
            status = 2;  % this means successful stop
            return;
        end
        
        % --- Set how long we'll wait for a new file to appear ---
        if (count == 2), timeout_time = 3*TR; end
        maxwait = fix(timeout_time/options.sleeptime);
        
        % --- See if a new file is ready ---
        [filename,stat] = get_next_file(RTpath,options.wildcard,lastfile,options.filesort);
        if (~stat), return; end
        
        % --- No new file yet, return and wait for next time to check ---
        if (isempty(filename))
            nwait = nwait + 1;
            if (nwait >= maxwait)
                fprintf(2,'ERROR: no new file found before %g second timeout occured (%1d waits).\n',timeout_time, nwait);
            else
                status = 1;
            end
            return
        end
        
        % --- Found the next file ---
        if (count == 1), tstart = tic(); end
        
        % --- Record stats ---
        tnow(count)      = toc(tstart);
        nwaits(count)    = nwait;
        filenames{count} = filename;
        if (options.liveprint)
            fprintf(1,'%3d) %s %10.3f (%5d) %s',count,filenames{count},tnow(count),nwaits(count));
        end
        
        % --- Read Dicom ---
        infile    = [RTpath filesep() filename];
        [~,froot] = fileparts(infile);
        outfile   = [options.scratchpath filesep() froot '.nii'];
        hdr       = ME_spm8_dicom_headers4(infile);
        if (isempty(hdr)),                                fprintf(2,'ERROR: Could not read new dicom file.\n'); return; end
        if (isempty(strfind(hdr{1}.ImageType,'MOSAIC'))), fprintf(2,'ERROR: Dicom is not from a mosaiced EPI scan.\n'); return; end
        if (count == 1), TR = hdr{1}.RepetitionTime/1000; end
        
        % --- Convert Dicom to NIFTI ---
        if (options.make_niftis)
            tmp = ME_spm8_dicom_convert4(hdr,outfile);
            if (isempty(tmp)),                            fprintf(2,'ERROR: Error converting dicom to nifti.\n'); return; end
            vol = spm_vol(outfile);
        end
        
        % --- Moco ---
        if (options.do_moco)
            % --- Moco not primed, so prime and that's all ---
            if (isempty(vol0))
                vol0 = vol;
                nifti_coreg(vol0,[],'rigidbody_fast','',options.moco_quality);
                center = nifti_get_center(vol0); % need for rmsmotion calc
                
                % --- compute motion of new volume ---
            else
                if (options.moco_only)                                              % compute motion mat only, no reslice
                    routfile = '';
                else
                    routfile = [options.scratchpath filesep() 'r_' froot '.nii'];   % compute motion AND reslice
                end
                mat = nifti_coreg(vol0,vol,'rigidbody_fast',routfile);              % moco calc
                
                % --- Record motion estimates ---
                rmsmotion(count)  = rmsdiff(last_mat,mat,center);   % Jenkinson RMS motion metric
                pars              = spm_imatrix(mat/last_mat);      % convert to translations and rotations
                mocopars(count,:) = -pars(1:6);                     % take "-" to match SPM plots convention
                
                % --- Write moco params to scratch folder as well ---
                if (~options.moco_only)
                    routfile = [options.scratchpath filesep() 'r_' froot '_mcpar.txt'];
                    fp=fopen(routfile,'w'); fprintf(fp,'%g %g %g %g %g %g \n',-pars(1:6)); fclose(fp);
                end
                
                % --- record motion as relative to previous volume
                if (options.moco_rel), last_mat = mat; end
            end
            if (options.liveprint)
                fprintf(' %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',mocopars(count,1:3),mocopars(count,4:6)*180/pi);
            end
        else
            if (options.liveprint)
                fprintf(1,'\n'); % finish  the print statement above
            end
        end
        
        % --- Plot moco ---
        if (options.do_moco && options.plot_moco && count > 1)
            
            % --- Get latest moco plot data ---
            plotpars(count,1:3) = abs(mocopars(count,1:3));
            plotpars(count,4)   = rmsmotion(count);
            % dmax = max(plotpars(count,:));
            % if (dmax > options.moco_thresh), beep; end
           
            % --- plot it ---
            figure(figh);
            hold on
            if (count > xrange(2)), xrange(2) = 2 * xrange(2); end      % expand x-range as needed
            new_range = [xrange yrange/restless_plotzoom];              % possibly new plot range
            old_range = axis();                                         % current plot window range
            if (new_range(2) ~= old_range(2) || new_range(4) ~= old_range(4))
                axis(new_range);
            end
            if (options.plot_rms)
                plot(count-1:count,plotpars(count-1:count,4),'b','LineWidth',2)
            else
                plot(count-1:count,plotpars(count-1:count,1),'b','LineWidth',2)
                plot(count-1:count,plotpars(count-1:count,2),'g','LineWidth',2)
                plot(count-1:count,plotpars(count-1:count,3),'r','LineWidth',2)
            end
            hold off
            drawnow
        end
        
        % --- Reset for next file to be found ---
        nwait = 0;
        count = count + 1;
        lastfile = filename;
        
        status = 1; % this means success, return and come back for more
        
    % --------------------------------------------
    % --- MODE = 2 End (called once)
    % --------------------------------------------
    case 2
        if (options.do_moco && options.plot_moco), set(figh,'color',stopcolor); end
        
        if (~options.liveprint)
            for i=1:count-1
                fprintf(1,'%3d) %s %10.3f (%5d) %s',i,filenames{i},tnow(i),nwaits(i));
                if (options.do_moco)
                    fprintf(' %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',mocopars(i,1:3),mocopars(i,4:6)*180/pi);
                else
                    fprintf(1,'\n'); % finish  the print statement above
                end
            end
        end
        
        status = 1; % success
end

return



