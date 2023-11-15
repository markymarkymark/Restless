function [filename,stat] = get_next_file(folder,wildcard,lastfile,sortoption)
% Return next file in a folder
%
% inputs:
%   lastfile - if '' (or []), then return last file in filelist
%              else return next file from filelist after "lastfile"
%
%   sortoption = 'date' - sort filelist by modification date
%   sortoption = 'name' - sort filelist alphabetically by filename
%                          default is to take filelist as dir() returns it
%
% outputs:
%   filename - next file, '' if none.
%   stat     - 1, success
%              0, error

% --- use much faster Unix system() method if unix ---
if (isunix() && 0), [filename,stat] = get_next_file_fast(folder,wildcard,lastfile,sortoption); return; end

filename = '';
stat     = 0;  % error

% --- get directory listing ---
sstr   = [folder filesep() wildcard];
files  = dir(sstr);
nfiles = numel(files);
%if (nfiles == 0), fprintf(2,'ERROR: no files matched %s\n',sstr); return; end
if (nfiles == 0), stat = 1; return; end % no new file, return empty string

if (nfiles > 1)
    switch sortoption
        % --- sort found files by date, oldest to newest ---
        case 'date'
            dates    = extractfield(files,'datenum');
            [~,ix] = sort(dates);
            files    = files(ix);
            
        % --- sort found files by name, A to Z ---
        case 'name'
            names    = extractfield(files,'name');
            [~,ix] = sort(names);
            files    = files(ix);
            
        % --- no sort ---
        otherwise
    end
end

% --- find first file in list that comes after lastfile ---
if (~isempty(lastfile))
    i = nfiles;
    while(i > 0 && ~isequal(lastfile,files(i).name)), i = i -1; end                                     % fastest to count up from the end of the list
    if     (i == 0), fprintf(2,'ERROR: lastfile = "%s" not found in "%s\n"',lastfile,sstr); return      % no match for lastfile!
    elseif (i == nfiles);                                                                               % lastfile is last in filelist, return empty string
    else   filename = files(i+1).name;                                                                  % return next file after lastfile
    end
    
% --- No lastfile provided, return last file from list ---
else
    filename = files(nfiles).name;
end

stat = 1;
end

