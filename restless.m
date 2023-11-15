% to compile
%
%  mcc -m restless -a spm_dicom_dict.mat 

function varargout = restless(varargin)
% RESTLESS M-file for restless.fig
%      RESTLESS, by itself, creates a new RESTLESS or raises the existing
%      singleton*.
%
%      H = RESTLESS returns the handle to a new RESTLESS or the handle to
%      the existing singleton*.
%
%      RESTLESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESTLESS.M with the given input arguments.
%
%      RESTLESS('Property','Value',...) creates a new RESTLESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before restless_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to restless_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help restless

% Last Modified by GUIDE v2.5 07-Dec-2017 15:17:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @restless_OpeningFcn, ...
                   'gui_OutputFcn',  @restless_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before restless is made visible.
function restless_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to restless (see VARARGIN)

% Choose default command line output for restless
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes restless wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- initialize GUI objects ---
init_gui(handles);
set(hObject,'name','RestLess'); % get cute
movegui(hObject,'southeast');
set(hObject,'CloseRequestFcn',@restless_closegui); % handle user close of GUI
return


% --- Outputs from this function are returned to the command line.
function varargout = restless_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in CHECK_domoco.
function CHECK_domoco_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_domoco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_domoco


% --- Executes on button press in CHECK_mocoonly.
function CHECK_mocoonly_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_mocoonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_mocoonly


% --- Executes on button press in CHECK_mocorel.
function CHECK_mocorel_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_mocorel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_mocorel


% --- Executes on button press in CHECK_plotmoco.
function CHECK_plotmoco_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_plotmoco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_plotmoco


% --- Executes on button press in CHECK_primemoco.
function CHECK_primemoco_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_primemoco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_primemoco



function EDIT_primequal_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_primequal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_primequal as text
%        str2double(get(hObject,'String')) returns contents of EDIT_primequal as a double

options = read_options(handles);
update_gui(handles,options);
return

% --- Executes during object creation, after setting all properties.
function EDIT_primequal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_primequal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDIT_mocothresh_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_mocothresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_mocothresh as text
%        str2double(get(hObject,'String')) returns contents of EDIT_mocothresh as a double

options = read_options(handles);
update_gui(handles,options);

% --- Executes during object creation, after setting all properties.
function EDIT_mocothresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_mocothresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDIT_maxreps_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_maxreps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_maxreps as text
%        str2double(get(hObject,'String')) returns contents of EDIT_maxreps as a double

options = read_options(handles);
update_gui(handles,options);

% --- Executes during object creation, after setting all properties.
function EDIT_maxreps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_maxreps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDIT_ymax_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_ymax as text
%        str2double(get(hObject,'String')) returns contents of EDIT_ymax as a double

options = read_options(handles);
update_gui(handles,options);

% --- Executes during object creation, after setting all properties.
function EDIT_ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDIT_xmax_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_xmax as text
%        str2double(get(hObject,'String')) returns contents of EDIT_xmax as a double

options = read_options(handles);
update_gui(handles,options);

% --- Executes during object creation, after setting all properties.
function EDIT_xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BUTTON_choose1.
function BUTTON_choose1_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_choose1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options = read_options(handles);
path = uigetdir(options.RTtoppath,'Select the Real Time Dicom Export folder');
if (~ischar(path)), return; end
options.RTtoppath = path;
update_gui(handles,options);
return

% --- Executes on button press in BUTTON_choose2.
function BUTTON_choose2_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_choose2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options = read_options(handles);
path = uigetdir(options.scratchpath,'Select the scratch folder');
if (~ischar(path)), return; end
options.scratchpath = path;
update_gui(handles,options);
return

% --- Executes on button press in CHECK_emptyscratch.
function CHECK_emptyscratch_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_emptyscratch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_emptyscratch


% --- Executes on button press in CHECK_makenifti.
function CHECK_makenifti_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_makenifti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_makenifti


% --- Executes on button press in CHECK_liveprint.
function CHECK_liveprint_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_liveprint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_liveprint



function EDIT_sleeptime_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_sleeptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_sleeptime as text
%        str2double(get(hObject,'String')) returns contents of EDIT_sleeptime as a double

options = read_options(handles);
update_gui(handles,options);

% --- Executes during object creation, after setting all properties.
function EDIT_sleeptime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_sleeptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDIT_wildcard_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_wildcard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_wildcard as text
%        str2double(get(hObject,'String')) returns contents of EDIT_wildcard as a double


% --- Executes during object creation, after setting all properties.
function EDIT_wildcard_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_wildcard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CHECK_plotrms.
function CHECK_plotrms_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_plotrms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_plotrms



function EDIT_mocoqual_Callback(hObject, eventdata, handles)
% hObject    handle to EDIT_mocoqual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDIT_mocoqual as text
%        str2double(get(hObject,'String')) returns contents of EDIT_mocoqual as a double

options = read_options(handles);
update_gui(handles,options);
return

% --- Executes during object creation, after setting all properties.
function EDIT_mocoqual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDIT_mocoqual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RADIO_sortname.
function RADIO_sortname_Callback(hObject, eventdata, handles)
% hObject    handle to RADIO_sortname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RADIO_sortname


% --- Executes on button press in RADIO_sortdate.
function RADIO_sortdate_Callback(hObject, eventdata, handles)
% hObject    handle to RADIO_sortdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RADIO_sortdate


% --- Executes on button press in RADIO_sortnone.
function RADIO_sortnone_Callback(hObject, eventdata, handles)
% hObject    handle to RADIO_sortnone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RADIO_sortnone


% --- Executes on button press in BUTTON_go.
function BUTTON_go_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options = read_options(handles);
str     = get(hObject,'String');

if (isequal(str,'Start'))
    
    % --- create a timer object that will be the background process ---
    t = timer();
    t.Period = options.sleeptime;
%    t.ExecutionMode = 'fixedRate';
    t.ExecutionMode = 'fixedSpacing';
    t.StartFcn = {@restless_startfunc, options, handles}; 
    t.TimerFcn = {@restless_timerfunc, options}; 
    t.StopFcn  = {@restless_stopfunc,  options, handles}; 
    handles.timer = t;
    guidata(hObject, handles);
      
    start(handles.timer);

elseif (isequal(str,'Stop'))
    
    stop(handles.timer);
end
    
return

% --- Executes on button press in BUTTON_savedef.
function BUTTON_savedef_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_savedef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isequal(questdlg('Are you sure you want to change the default options?','Save Defaults','yes','no','no'),'no')), return; end
file    = get_defaults_file();
options = read_options(handles);
save(file,'options');
return


% --- Executes on button press in BUTTON_save.
function BUTTON_save_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = fileparts(get_defaults_file());
[file,path] = uiputfile([path filesep() '*.mat'],'Save current options');
if (~ischar(file)), return; end
options = read_options(handles);
save([path filesep() file],'options');
return

% --- Executes on button press in BUTTON_load.
function BUTTON_load_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = fileparts(get_defaults_file());
[file,path] = uigetfile([path filesep() '*.mat'],'Choose a saved options file');
if (~ischar(file)), return; end
options = restless_get_options([path filesep() file]);
if (~isempty(options)), update_gui(handles,options); end
return

% --- Executes on button press in BUTTON_zoomin.
function BUTTON_zoomin_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global restless_plotzoom

restless_plotzoom = restless_plotzoom * 2;

return

% --- Executes on button press in BUTTON_zoomout.
function BUTTON_zoomout_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global restless_plotzoom

restless_plotzoom = restless_plotzoom / 2;

return
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% User-written routines that interect with the GUI
% (above here is autogenerated code from the GUIDE app)
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


% -----------------------------------------------------------------------
function init_gui(handles)

global restless_plotzoom

options = restless_get_options(get_defaults_file()); % get options from defaults file
if (isempty(options)) % get hard coded default options
    options = restless_get_options(); 
end

update_gui(handles,options);

restless_plotzoom = 1;
return

% -----------------------------------------------------------------------
function update_gui(handles,options)

set(handles.TEXT_RTtoppath,'string',options.RTtoppath);
set(handles.TEXT_scratchpath,'string',options.scratchpath);

set(handles.CHECK_emptyscratch,          'Value',options.empty_scratch);
set(handles.CHECK_makenifti,          'Value',options.make_niftis);
set(handles.CHECK_liveprint,          'Value',options.liveprint);
set(handles.CHECK_domoco,          'Value',options.do_moco);
set(handles.CHECK_mocoonly,          'Value',options.moco_only);
set(handles.CHECK_mocorel,          'Value',options.moco_rel);
set(handles.CHECK_plotmoco,          'Value',options.plot_moco);
set(handles.CHECK_plotrms,          'Value',options.plot_rms);
set(handles.CHECK_primemoco,          'Value',options.prime_moco);

set(handles.EDIT_wildcard,'string',options.wildcard);
set(handles.EDIT_sleeptime,'string',sprintf('%6.3f',options.sleeptime));
set(handles.EDIT_mocothresh,'string',sprintf('%6.3f',options.moco_thresh));
set(handles.EDIT_maxreps,'string',sprintf('%1d',options.maxreps));
set(handles.EDIT_ymax,'string',sprintf('%6.3f',options.plot_ymax));
set(handles.EDIT_xmax,'string',sprintf('%1d',options.plot_xmax));
set(handles.EDIT_mocoqual,'string',sprintf('%6.3f',options.moco_quality));
set(handles.EDIT_primequal,'string',sprintf('%6.3f',options.moco_primequality));

switch (options.filesort)
    case 'date'
        set(handles.RADIO_sortdate, 'Value',1);
    case 'name'
        set(handles.RADIO_sortname, 'Value',1);
    otherwise
        set(handles.RADIO_sortnone, 'Value',1);
end
return

% -----------------------------------------------------------------------
function options = read_options(handles)

options.RTtoppath = get(handles.TEXT_RTtoppath,'string');
options.scratchpath = get(handles.TEXT_scratchpath,'string');

options.empty_scratch = get(handles.CHECK_emptyscratch,          'Value');
options.make_niftis = get(handles.CHECK_makenifti,          'Value');
options.liveprint = get(handles.CHECK_liveprint,          'Value');
options.do_moco = get(handles.CHECK_domoco,          'Value');
options.moco_only = get(handles.CHECK_mocoonly,          'Value');
options.moco_rel = get(handles.CHECK_mocorel,          'Value');
options.plot_moco = get(handles.CHECK_plotmoco,          'Value');
options.plot_rms = get(handles.CHECK_plotrms,          'Value');
options.prime_moco = get(handles.CHECK_primemoco,          'Value');

options.wildcard = get(handles.EDIT_wildcard,'string');
options.sleeptime = str2double(get(handles.EDIT_sleeptime,'string'));
options.moco_thresh = str2double(get(handles.EDIT_mocothresh,'string'));
options.maxreps = round(str2double(get(handles.EDIT_maxreps,'string')));
options.plot_ymax = str2double(get(handles.EDIT_ymax,'string'));
options.plot_xmax = round(str2double(get(handles.EDIT_xmax,'string')));
options.moco_quality = str2double(get(handles.EDIT_mocoqual,'string'));
options.moco_primequality = str2double(get(handles.EDIT_primequal,'string'));

if (get(handles.RADIO_sortdate,'Value'))
    options.filesort = 'date';
elseif (get(handles.RADIO_sortname,'Value'))
    options.filesort = 'name';
else  
    options.filesort = '';
end

return

% -----------------------------------------------------------------------
function file = get_defaults_file()

persistent path

if (isempty(path))
    if (ispc())
        home = [getenv('HOMEDRIVE') getenv('HOMEPATH') '\Documents'];
    else
        home = getenv('HOME');
    end
    if (isdir([home filesep() 'data' filesep() 'restless']))
        path = [home filesep() 'data' filesep() 'restless'];
    elseif (isdir([home filesep() 'data' filesep() 'RTfmri']))
        path = [home filesep() 'data' filesep() 'RTfmri'];
    elseif (isdir(home))
        path = home;
    else
        path = pwd();
    end
end

%path = fileparts(mfilename('fullpath'));
%path = fileparts(which('restless.m'));
file = [path filesep() 'restless_defaults.mat'];
return


% -----------------------------------------------------------------------
function restless_closegui(src,callbackdata)

if (ishghandle(1)), close(1); end;  % close the plot window
closereq();                         % run normal GUI close function

return


% -----------------------------------------------------------------------
function restless_startfunc(tobj, event, options, handles)

%disp('Starting')
set(handles.BUTTON_go, 'BackgroundColor',[0.8 0 0]);
set(handles.BUTTON_go, 'String','Stop');
drawnow
status = restless_engine(0, options);
if (status ~= 1), stop(tobj); end
return

% -----------------------------------------------------------------------
function restless_timerfunc(tobj, event, options)

%disp('running')
status = restless_engine(1, options);
if (status ~= 1), stop(tobj); end
return

% -----------------------------------------------------------------------
function restless_stopfunc(tobj, event, options, handles)

%disp('Stoping and deleting timer object')
delete(tobj);
restless_engine(2, options);
set(handles.BUTTON_go, 'BackgroundColor',[0 0.8 0]);
set(handles.BUTTON_go, 'String','Start');
drawnow
return
