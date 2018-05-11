function varargout = PMA_analyzer_v3(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PMA_analyzer_v3_OpeningFcn, ...
                   'gui_OutputFcn',  @PMA_analyzer_v3_OutputFcn, ...
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
function PMA_analyzer_v3_OpeningFcn(hObject, eventdata, handles, varargin)


handles.output = hObject;

guidata(hObject, handles);
initialize(handles,1);
addlistener(handles.slider_windowPMA,'Value','PostSet',@listener_slider_windowPMA_Callback);
addlistener(handles.frame_index,'Value','PostSet',@listener_frame_index_Callback);

function varargout = PMA_analyzer_v3_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%  PMA WINDOW  %%%%%%%%%%%%%%%%%%%%%%%%%%
function windowPMA_CreateFcn(hObject, eventdata, handles)
function windowPMA_ButtonDownFcn(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%  BACKGROUND  %%%%%%%%%%%%%%%%%%%%%%%%%%
function background_WindowButtonMotionFcn(hObject, eventdata, handles)

%JWK left click, right click, scroll
function background_WindowButtonDownFcn(hObject, eventdata, handles)
global pma;

cursor_PMA = get(handles.windowPMA, 'CurrentPoint');
cursor_x = fix(cursor_PMA(1,1)); cursor_y =fix(cursor_PMA(1,2));

if get(handles.frame_index,'value') == 0 || cursor_x <= 3 || ...
        cursor_x >= pma.width-3 || cursor_y <= 3 || cursor_y >= pma.height-3
    return
end

mbutton = get(hObject, 'SelectionType');

disp(mbutton)

if strcmpi(mbutton, 'normal')
    index = istherespot(cursor_x,cursor_y,handles);
   
    if index > 0 %jwk mouse is over a peak get selected peak
        set_tmpspot(index);
        make_tmpspotframe()
        draw_frame()
        draw_trace()
    else         %jwk mouse is not over a peak make current peak
        set_newtmpspot(cursor_x,cursor_y)
        make_tmpspotframe()
        draw_frame()
        draw_trace()
    end
end

% --- Executes on scroll wheel click while the figure is in focus.
function background_WindowScrollWheelFcn(hObject, eventdata, handles)
global pma;
if eventdata.VerticalScrollCount < 0
        move_slider(-1)
elseif eventdata.VerticalScrollCount > 0
        move_slider(1)
end

%JWK_FIX 'a' docking 's' hemifusion 'd'fusion 'z' delete
function background_WindowKeyPressFcn(hObject, eventdata, handles)
global g_tmpspot;
global g_spots
global MAX_NUM_SPOTS;

spot_exist = 1; 
if g_tmpspot.index < 1 || g_tmpspot.index > MAX_NUM_SPOTS
   spot_exist = 0; 
end
%JWK_INCOMPLETE
display(eventdata.Key)
switch eventdata.Key
    case 'a'
        if spot_exist ==  1
        else
            if g_tmpspot.x(1) == 0
                return
            else
            disp('ADD')
            add_tmpspots(0)
            draw('normal')
            end
        end
    case 'q'
        if spot_exist ==  1
        else
            if g_tmpspot.x(1) == 0
                return
            else
                disp('ADD')
                add_tmpspots(1)
                draw('normal')
            end
        end
    case 's'
        if spot_exist ==  0
        else
           set_fusiontime(g_tmpspot.index,1)
           draw('normal')
        end
    case 'd'
        if spot_exist ==  0
        else
           set_fusiontime(g_tmpspot.index,2)
           draw('normal')
        end
    case 'z'
        if spot_exist == 1
            disp('DEL')
            display('this is current index')
            del_tmpspots()
        else
        end 
                  
end
% jwk_later cannot overlap
if g_tmpspot.x > 0
    if strcmp(eventdata.Key, 'leftarrow')
        move_tmpspot('l')
    elseif strcmp(eventdata.Key, 'rightarrow')
        move_tmpspot('r')
    elseif strcmp(eventdata.Key, 'uparrow')
        move_tmpspot('u')
    elseif strcmp(eventdata.Key, 'downarrow')
        move_tmpspot('d')
    end
    
end
function background_KeyPressFcn(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE  %%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize(handles,is_initial)
global ghandles;    ghandles = handles;
initialize_etc();
if is_initial ~= 1
    initialize_tmpspot(); initialize_spots(); 
    initialize_figure(handles); initialize_coeff(handles);
end

function initialize_figure(handles)
% clear figures
cla(handles.windowPMA);
cla(handles.windowTRACE)
cla(handles.windowHIST);
% set figure axis
set(handles.windowPMA,  'xTick', [], 'yTick', []);
colormap(handles.windowPMA,gray(128))

% initialize edit string

function initialize_tmpspot()
global g_tmpspot;
global pma;
global MAX_NUM_SPOTS

g_tmpspot   = struct('index',-1,'x',zeros(1,2),'y',zeros(1,2),...
            'donor',zeros(1,pma.len),'acceptor',zeros(1,pma.len),'etime',ones(3,1)*pma.len);
        
g_tmpspot.index     = -1;
g_tmpspot.x         = zeros(1,2);
g_tmpspot.y         = zeros(1,2);

g_tmpspot.donor     = zeros(1,pma.len);
g_tmpspot.acceptor  = zeros(1,pma.len);
% jwk_change
% make_tmpspotframe(); %this clears tmpspotframe
% draw_frame()
function initialize_spots()
global g_spots
global pma
global MAX_NUM_SPOTS

g_spots     = struct('total',0,'x',zeros(MAX_NUM_SPOTS,2),'y',zeros(MAX_NUM_SPOTS,2),...
            'donor',zeros(MAX_NUM_SPOTS,pma.len),'acceptor',zeros(MAX_NUM_SPOTS,pma.len),'etime',ones(MAX_NUM_SPOTS,3)*pma.len);
        
g_spots.total     = 0;
g_spots.x         = zeros(MAX_NUM_SPOTS,2);
g_spots.y         = zeros(MAX_NUM_SPOTS,2);

g_spots.donor     = zeros(MAX_NUM_SPOTS,pma.len);
g_spots.acceptor  = zeros(MAX_NUM_SPOTS,pma.len);
% jwk_change
% make_maskframe()
% draw_frame()
function initialize_etc()
global MAX_NUM_SPOTS
global NUM_BIN
global g_maskframe
global g_tmpspotframe
global IS_MEAN
global dir_script

dir_script = mfilename('fullpath'); dir_script = dir_script(1:end-numel(mfilename));
g_maskframe         = zeros(512,512);
g_tmpspotframe      = zeros(512,512);
% MAX_NUM_SPOTS       = 10;
MAX_NUM_SPOTS       = 5000;

IS_MEAN         = 0;
NUM_BIN         = 5; %actually number of frames in bin

function initialize_coeff(handles)
global FRAME_MAX

%value
default_num_bin     = 50;
default_FRAME_MAX    = 100;

%set
set(handles.edit_NUMBIN, 'value', default_num_bin);
set(handles.edit_NUMBIN, 'string',default_num_bin);
set(handles.edit_FRAMEMAX, 'value', default_FRAME_MAX);
set(handles.edit_FRAMEMAX, 'String', default_FRAME_MAX);

%get
FRAME_MAX = str2num(get(handles.edit_FRAMEMAX, 'string'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUTTON CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function button_pma()
function btn_loadPMA_Callback(hObject, eventdata, handles)
% clear global
[filename, pathname] = uigetfile( ...
{  '*.pma','pma-files (*.pma)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file', ...
   'MultiSelect', 'on');
if filename == 0 
    return
end
cd(pathname)
load_pmafile(pathname,filename(1:end-4),handles);
initialize(handles,0);
draw('initial');
function btn_loadPMA_ButtonDownFcn(hObject, eventdata, handles)

function btn_loadSPOTS_Callback(hObject, eventdata, handles)
global pma
global g_spots
name_gspot = [pma.dir,'new_',pma.name,'_gspots.mat']
A = load(name_gspot);
g_spots = A.g_spots;
g_spots.etime
draw('normal')

function edit_CONTRAST_Callback(hObject, eventdata, handles)
str   = get(handles.edit_CONTRAST, 'string');
contrast   = str2num(str);
set(handles.windowPMA,  'xTick', [], 'yTick', []);
colormap(handles.windowPMA,gray(contrast))

function edit_NUMBIN_Callback(hObject, eventdata, handles)
global NUM_BIN
str   = get(handles.edit_NUMBIN, 'string');
NUM_BIN   = str2num(str);
draw('normal')
function edit_FRAMEMAX_Callback(hObject, eventdata, handles)
global FRAME_MAX
str   = get(handles.edit_FRAMEMAX, 'string');
FRAME_MAX   = str2num(str);
draw('normal')
function btn_save_Callback(hObject, eventdata, handles)
disp('SAVE PKS, TRACE')
write_variables()
% write_pks()
% write_trace()
write_hist()
% write_config()
function frame_index_Callback(hObject, eventdata, handles)
frame_index = str2double(get(handles.frame_index, 'String'));
set(handles.slider_windowPMA, 'Value', frame_index);
draw_frame();

%added by bhawk 20180508
function listener_frame_index_Callback(hObject, eventdata)
global ghandles
frame_index = str2double(get(ghandles.frame_index, 'String'));
if frame_index > 2500
    disp('frame_index listener triggered')
    set(ghandles.background, 'Color','red')
else
    set(ghandles.background, 'Color','white')
end

function listener_slider_windowPMA_Callback(hObject, eventdata)
global ghandles
frame_index = round(get(ghandles.slider_windowPMA,'value'));
set(ghandles.frame_index, 'value', frame_index);
set(ghandles.frame_index, 'string', frame_index);
draw('slider');
function edit_CONTRAST_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_FRAMEMAX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_NUMBIN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_endFRAME_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_windowPMA_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function frame_index_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function static_pmalen_CreateFcn(hObject, eventdata, handles)
function clear_global()
clear pma
clear g_spot
clear MAX_NUM_SPOTS
clear ghandles
clear g_tmpspot
clear MAX_NUM_SPOTS
clear NUM_BIN
clear g_maskframe
clear g_tmpspotframe
clear IS_MEAN
clear dir_script
function slider_windowPMA_Callback(hObject, eventdata, handles)

%%% PRIVATE %%%
function load_pmafile(pathname,filename,handles)
global pma
downstairs = 1;
if downstairs == 1
    fullpath = [pathname,filename,'.pma'];
    field_len     = 'len'; field_width   = 'width';  field_height  = 'height';
    field_area    = 'area';field_frames  = 'frames'; field_name    = 'name';
    field_dir     = 'dir'; field_savename= 'savename';
    
    fileinfo    = dir(fullpath);    filesize    = fileinfo.bytes;
    fid         = fopen(fullpath);  
    pma_width   = fread(fid, 1, 'int16')
    pma_height  = fread(fid, 1, 'int16') 
    pma_len     = floor((filesize-4)/(pma_width*pma_height))
    pma_frames  = zeros(pma_height,pma_width,pma_len);
else
    fullpath = [pathname,filename,'.pma'];
    field_len     = 'len'; field_width   = 'width';  field_height  = 'height';
    field_area    = 'area';field_frames  = 'frames'; field_name    = 'name';
    field_dir     = 'dir'; field_savename= 'savename';
    
    fileinfo    = dir(fullpath);    filesize    = fileinfo.bytes;
    fid         = fopen(fullpath);  pma_width   = fread(fid, 1, 'int16');
    pma_height  = fread(fid, 1, 'int16'); pma_len     = (filesize-4)/(pma_width*pma_height);
    pma_frames  = zeros(pma_height,pma_width,pma_len);
    
end
pma = struct(field_len,pma_len,field_width,pma_width,field_height,pma_height,field_area,pma_width*pma_height,...
    field_frames,pma_frames,field_name,filename,field_dir,pathname,field_savename,'');

set(handles.static_filename,'string',pma.name);

h = waitbar(0, 'Loading frames from pma file...');
for i=1:pma.len
    pma.frames(:,:,i) = fread(fid, [512,512], 'uint8');
    pma.frames(:,:,i) = rot90(pma.frames(:,:,i));
    waitbar(i/pma_len);
end
close(h)
fclose(fid);

set(handles.static_pmalen, 'value', pma.len);
set(handles.static_pmalen, 'string', pma.len);
set(handles.frame_index, 'value', 1);
set(handles.frame_index, 'string', 1);

% %%%%%%%%%%%%%%% slider setup %%%%%%%%%%%%%%%%%%%%%%%%
set(handles.slider_windowPMA, 'Min', 1);
set(handles.slider_windowPMA, 'Max', pma.len);
set(handles.slider_windowPMA, 'SliderStep', [1,20]/(pma.len-1));
set(handles.slider_windowPMA, 'value', 1);

function write_config()
global ghandles
global pma

start_frame     = str2num(get(ghandles.edit_startFRAME, 'string'));
end_frame       = str2num(get(ghandles.edit_endFRAME, 'string'));

donor_level     = str2num(get(ghandles.edit_DonorLevel, 'string'));
acceptor_level  = str2num(get(ghandles.edit_AcceptorLevel, 'string'));
leakage         = str2num(get(ghandles.edit_leakage, 'string'));
correction_factor   = str2num(get(ghandles.edit_Correction, 'string'));
num_bin             = str2num(get(ghandles.edit_NUMBIN, 'string'));
fid = fopen([pma.dir,'new_',pma.name,'_config.txt'],'wt');
fprintf(fid,'start frame   :: %f \n',start_frame);
fprintf(fid,'end frame     :: %f \n',end_frame);
fprintf(fid,'donor    level:: %f \n',donor_level);
fprintf(fid,'acceptor level:: %f \n',acceptor_level);
fprintf(fid,'leakage       :: %f \n',leakage);
fprintf(fid,'correction    :: %f \n',correction_factor);
fprintf(fid,'num bin       :: %f \n',num_bin);

fclose(fid);

function write_trace()
%jwk trace write
global g_spots;
global pma;
global MAX_NUM_SPOTS

%jwk get information
fid = fopen([pma.dir,'new_',pma.name,'.traces'], 'w');
fwrite(fid,pma.len,'int32');
fwrite(fid,g_spots.total*2,'int16');
disp('WRITE TRACE ::')
disp('pma.len g_spots.total')
disp([pma.len g_spots.total])

tmp = zeros(g_spots.total, pma.len);
spot_index = find(g_spots.x(:,1) > 0);
trace_index  = 0;
if ~isempty(spot_index)
    for i = 1: g_spots.total
        trace_index = trace_index + 1;
        tmp(trace_index,:) = g_spots.donor(i,:);
        trace_index = trace_index + 1;
        tmp(trace_index,:) = g_spots.acceptor(i,:);
    end
end

index       = (1:g_spots.total*pma.len*2);
fdata       = zeros(g_spots.total*pma.len*2,1,'int16');
fdata(index)= tmp(index);
fwrite(fid,fdata,'int16');
fclose(fid);  
function write_variables()
global pma;
global g_spots
name_gspot = [pma.dir,'new_',pma.name,'_gspots.mat']
save(name_gspot,'g_spots');

function write_hist()
global g_spots;
global pma;
global NUM_BIN;


% time of fused vesicles "docked hemifused fused"


fuse_index = find(g_spots.etime(:,3) < pma.len);
dbg_str = sprintf('SAVING number spot_index:%d',numel(fuse_index)); disp(dbg_str);
time_index  = 0;
if ~isempty(fuse_index)
save_name = [pma.dir,'new_',pma.name,'_fused_vesicles'];
save_name = sprintf('%s_%d.time',save_name,numel(fuse_index));
fid = fopen(save_name,'wt');
    for i = fuse_index'
        time_index = time_index + 1;
        g_spots.etime(i,4) = g_spots.etime(i,3)-g_spots.etime(i,1);
        if g_spots.etime(i,4) == 0
            g_spots.etime(i,5) = 1;
        end
        fprintf(fid,'%d,\t   %d,\t   %d,\t   %d\n',time_index,g_spots.etime(i,1),g_spots.etime(i,2),g_spots.etime(i,3));
    end
    fclose(fid);
end

%added 20180511
fusion_instant = sum(g_spots.etime(:,5));

% Diff btn "hemi-dock  fused-hemi  fused-docked"
fid = fopen([pma.dir,'new_',pma.name,'_fused_diff.time'],'wt');

spot_index = find(g_spots.etime(:,3) < pma.len);
dbg_str = sprintf('SAVING number spot_index:%d',numel(spot_index)); disp(dbg_str);
time_index  = 0;
if ~isempty(spot_index)
    for i = spot_index'
        time_index = time_index + 1;
        fprintf(fid,'%d,\t   %d,\t   %d,\t   %d\n',time_index,g_spots.etime(i,2)-g_spots.etime(i,1),g_spots.etime(i,3)-g_spots.etime(i,2), g_spots.etime(i,3)-g_spots.etime(i,1));
    end
end
fclose(fid);

% time of docked vesicles "dock  fused-hemi  fused-docked"


docked_index = find(g_spots.etime(:,1) < pma.len);
dbg_str = sprintf('SAVING number spot_index:%d',numel(docked_index)); disp(dbg_str);
time_index  = 0;
if ~isempty(docked_index)
save_name = [pma.dir,'new_',pma.name,'_docked_vesicles'];
save_name = sprintf('%s_%d.time',save_name,numel(docked_index));
fid = fopen(save_name,'wt');
    for i = docked_index'
        time_index = time_index + 1;
        fprintf(fid,'%d,\t   %d,\t   %d,\t   %d\n',time_index,g_spots.etime(i,1),g_spots.etime(i,2),g_spots.etime(i,3));
    end
end
fclose(fid);

%added 20180511 write summary file for events: dock fusion and instant
%fusion

save_name = [pma.dir,'new_',pma.name,'_summary.txt'];
fid = fopen(save_name,'wt');
fprintf(fid,'Docked:\t\t\t%d\n',length(docked_index));
fprintf(fid,'Fused:\t\t\t%d\n',length(fuse_index));
fprintf(fid,'Instant fusion:\t\t%d\n',fusion_instant);
fclose(fid);

% docked_index   = find(g_spots.etime(:,1) ~= pma.len);
save_name = [pma.dir,'new_',pma.name,'_hist_fusion'];
fused_index    = find(g_spots.etime(:,3) ~= pma.len);
if ~isempty(fused_index)
    fuse_time  = g_spots.etime(fused_index,3)-g_spots.etime(fused_index,1);
    max_fuse_time = ceil(max(fuse_time)/NUM_BIN)*NUM_BIN;
    xbins = 0:NUM_BIN:max_fuse_time;
    xbins = xbins + NUM_BIN/2;
    [hist_data b] = hist(fuse_time,xbins);
    xlswrite(save_name,hist_data'/numel(fused_index));
end

save_name = [pma.dir,'new_',pma.name,'_hist_docking'];
dock_time      = g_spots.etime(g_spots.etime(:,1) ~= pma.len,1);
max_dock_time = ceil(max(dock_time)/NUM_BIN)*NUM_BIN;
xbins = 0:NUM_BIN:max_dock_time;
[hist_data b] = hist(dock_time,xbins);
cum_hist_data = cumsum(hist_data);
xlswrite(save_name,cum_hist_data');


function move_slider(direction)
global ghandles
global pma
frame_index = str2double(get(ghandles.frame_index, 'String'));
if direction < 0
    frame_index = frame_index -1;
    if frame_index > 0
        set(ghandles.frame_index, 'value', frame_index);
        set(ghandles.frame_index, 'string', frame_index);
        draw('slider');
    else
    end
else
    frame_index = frame_index + 1;
    if frame_index <= pma.len
        set(ghandles.frame_index, 'value', frame_index);
        set(ghandles.frame_index, 'string', frame_index);
        draw('slider');
    else
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     SPOTS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% PUBLIC %%%%%%%
function spots_pma()
% function add_spots_pks(x,y,donor,acceptor)
function add_spots_pks(x,y)
global g_spots
global MAX_NUM_SPOTS
global pma
x_d = x(1:2:end); x_a = x(2:2:end);
y_d = y(1:2:end); y_a = y(2:2:end);

num_spots = numel(x_d);

dbg_str = sprintf('size x_d: %d ',num_spots);
disp(dbg_str);

if num_spots > MAX_NUM_SPOTS
    disp('PKS is has more spots than MAX_NUM_SPOTS')
    num_spots = MAX_NUM_SPOTS;
end

donor       = zeros(num_spots,pma.len);
acceptor    = zeros(num_spots,pma.len);

dbg_str = sprintf('size acceptor:'); disp(dbg_str);
disp(size(acceptor))
dbg_str = sprintf('size donor:'); disp(dbg_str);
disp(size(donor))

circle = [ ...       
        0.0037    0.0743    0.2019    0.0743    0.0037    ; ...
        0.0123    0.2466    0.6703    0.2466    0.0123    ; ...
        0.0183    0.3679    1.0000    0.3679    0.0183    ; ...
        0.0123    0.2466    0.6703    0.2466    0.0123    ; ...
        0.0037    0.0743    0.2019    0.0743    0.0037    ; ...     
];

radius = 2;
x_d = round(x_d);     x_a = round(x_a);
y_a = 512*ones(num_spots,1)-round(y_a);
y_d = 512*ones(num_spots,1)-round(y_d);

dbg_str = sprintf('size x_d: %d ',numel(x_d)); disp(dbg_str);
dbg_str = sprintf('size y_d: %d ',numel(y_d)); disp(dbg_str);
dbg_str = sprintf('size x_a: %d ',numel(x_a)); disp(dbg_str);
dbg_str = sprintf('size y_a: %d ',numel(y_a)); disp(dbg_str);


%get donor accpetor trace from PMA
for j = 1: num_spots
    for i = 1:pma.len
        donor(j,i) = ...
            sum(sum(pma.frames(y_d(j)-radius:y_d(j)+radius,x_d(j)-radius:x_d(j)+radius,i).*circle));
    end
    
    %JWK_INCOMPLETE
    for i = 1:pma.len
        acceptor(j,i) = ...
            sum(sum(pma.frames(y_a(j)-radius:y_a(j)+radius,x_a(j)-radius:x_a(j)+radius,i).*circle));
    end
end
g_spots.x(1:num_spots,1)        = x_d;
g_spots.y(1:num_spots,1)        = y_d;
g_spots.x(1:num_spots,2)        = x_a;
g_spots.y(1:num_spots,2)        = y_a;
g_spots.total                   = num_spots;
g_spots.donor(1:num_spots,:)     = donor;
g_spots.acceptor(1:num_spots,:)  = acceptor;
dbg_str = sprintf('size x_d: %d ',numel(g_spots.x(1:num_spots,1))); disp(dbg_str);
dbg_str = sprintf('size y_d: %d ',numel(g_spots.y(1:num_spots,1))); disp(dbg_str);
dbg_str = sprintf('size x_a: %d ',numel(g_spots.x(1:num_spots,2))); disp(dbg_str);
dbg_str = sprintf('size y_a: %d ',numel(g_spots.y(1:num_spots,2))); disp(dbg_str);

[tmp index ] = sort(g_spots.y(1:num_spots,1));

g_spots.x(1:num_spots,1)             =g_spots.x(index,1);
g_spots.x(1:num_spots,2)             =g_spots.x(index,2);
g_spots.y(1:num_spots,1)             =g_spots.y(index,1);
g_spots.y(1:num_spots,2)             =g_spots.y(index,2);
g_spots.donor(1:num_spots,:)       =g_spots.donor(index,:);
g_spots.acceptor(1:num_spots,:)    =g_spots.acceptor(index,:);

initialize_tmpspot()

%type 0 is add current frame index type 1 is add closest spike previous to
%current frame index
function add_tmpspots(type)
global g_spots;
global g_tmpspot;
global MAX_NUM_SPOTS;
global pma
global ghandles
current_index = get_currentframeindex()
current_frame = pma.frames(:,:,current_index);

donor_frame = current_frame(g_tmpspot.x(1)-1:g_tmpspot.x(1)+1,g_tmpspot.y(1)-1:g_tmpspot.y(1)+1);
accep_frame = current_frame(g_tmpspot.x(2)-1:g_tmpspot.x(2)+1,g_tmpspot.y(2)-1:g_tmpspot.y(2)+1);

if type == 1 
    search_frame = smooth(g_tmpspot.donor(1:current_index));
    A = search_frame(1:end-3);
    B = search_frame(4:end);
    C = B-A;
    [max_diff max_diff_index]   = max(C);
    current_index = max_diff_index + 2;
end
%JWK_INCOMPLETE calibrate
% [maxD,ind] = max(donor_frame(:));
% [Dx,Dy] = ind2sub(size(donor_frame),ind)
% 
% [maxA,ind] = max(accep_frame(:));
% [Ax,Ay] = ind2sub(size(accep_frame),ind)
% 
% Dx = Dx -2; Dy = Dy -2;
% Ax = Ax -2; Ay = Ay -2;

Dx = 0; Dy = 0; Ax = 0; Ay = 0;

%JWK_incomplete
%check if any available spots
index = get_insertspotindex(g_tmpspot.y(1));
display(index)
if (index > 0) && index <= MAX_NUM_SPOTS &&g_spots.total < MAX_NUM_SPOTS
   g_spots.x(index,1)           = round(g_tmpspot.x(1))-Dx;
   g_spots.y(index,1)           = round(g_tmpspot.y(1))-Dy;
   g_spots.x(index,2)           = round(g_tmpspot.x(2))-Ax;
   g_spots.y(index,2)           = round(g_tmpspot.y(2))-Ay;
   g_spots.donor(index,:)       = g_tmpspot.donor;
   g_spots.acceptor(index,:)    = g_tmpspot.acceptor;
   g_spots.total                = g_spots.total+1;
   g_spots.etime(index,:)       = g_tmpspot.etime;
   g_spots.etime(index,1)       = current_index;
   
    % JWK_INCOMPLETE initialize_tmpspot();
    set_tmpspot(index)
%     make_maskframe();
%     make_tmpspotframe();
%     draw_frame();
else

end
function del_tmpspots()
global g_tmpspot;
global MAX_NUM_SPOTS;
display('index')
g_tmpspot.index
%check if selected
if (g_tmpspot.index > 0)||(g_tmpspot.index <= MAX_NUM_SPOTS) 
    initialize_spots_index(g_tmpspot.index);
%     initialize_tmpspot();
    make_maskframe();
    make_tmpspotframe();
    draw_frame();
else
    disp('NOT EXISTING SPOT')
end
function index = get_insertspotindex(y)
global g_spots
global MAX_NUM_SPOTS
global pma
index = -1;
if g_spots.total < MAX_NUM_SPOTS
   tmp = find(g_spots.etime(:,1)== pma.len);
   index = tmp(1);
end

%JWK_INCOMPLETE
function move_tmpspot(arrow)
global g_tmpspot
global g_spots
global pma
radius = 3;
curr_x = g_tmpspot.x(1); curr_y = g_tmpspot.y(1);
switch arrow
    case 'l'
        if (curr_x - 1) < radius
            return
        end
        set_newtmpspot(curr_x -1,curr_y)
    case 'r'
        if (curr_x + 1) > (pma.width - radius)
            return
        end
        set_newtmpspot(curr_x+1,curr_y)
    case 'd'
        if (curr_y + 1) > (pma.height - radius)
            return
        end
        set_newtmpspot(curr_x,curr_y+1)
    case 'u'
        if (curr_y - 1) < radius
            return
        end
        set_newtmpspot(curr_x,curr_y-1)       
end
%JWK_INCOMPLETE

function set_newtmpspot(x,y)
global g_tmpspot;
global pma;
%jwk make gaussian
%JWK_INCOMPLETE


circle = [ ...       
        0.0037    0.0743    0.2019    0.0743    0.0037    ; ...
        0.0123    0.2466    0.6703    0.2466    0.0123    ; ...
        0.0183    0.3679    1.0000    0.3679    0.0183    ; ...
        0.0123    0.2466    0.6703    0.2466    0.0123    ; ...
        0.0037    0.0743    0.2019    0.0743    0.0037    ; ...     
];

radius = 2;

disp('set new tmpspot x y')
disp([x y])
mapped = make_mapped_coordinate(x,y);    


if x > 256
    g_tmpspot.index     = -1;
    disp('ACCEPTOR')
    g_tmpspot.x(1)      = floor(mapped(1));
    g_tmpspot.y(1)      = floor(mapped(2));
    
    g_tmpspot.x(2)      = x;
    g_tmpspot.y(2)      = y;
    
    g_tmpspot.etime(1)  = get_currentframeindex();
    g_tmpspot.etime(2)  = pma.len;
    g_tmpspot.etime(3)  = pma.len;
    disp('x_1 x_2 y_1 y_2')
    disp([g_tmpspot.x(1) g_tmpspot.x(2) g_tmpspot.y(1) g_tmpspot.y(2)])
else
    disp('DONOR')
    g_tmpspot.index     = -1;
    g_tmpspot.x(1)      = x;
    g_tmpspot.y(1)      = y;
    
    g_tmpspot.x(2)      = floor(mapped(1));
    g_tmpspot.y(2)      = floor(mapped(2));
    g_tmpspot.etime(1)  = get_currentframeindex();
    g_tmpspot.etime(2)  = pma.len;
    g_tmpspot.etime(3)  = pma.len;
    disp('x_1 x_2 y_1 y_2')
    disp([g_tmpspot.x(1) g_tmpspot.x(2) g_tmpspot.y(1) g_tmpspot.y(2)])
end
%JWK_INCOMPLETE
for i = 1:pma.len
    g_tmpspot.donor(i) = ...
        sum(sum(pma.frames(g_tmpspot.y(1)-radius:g_tmpspot.y(1)+radius,g_tmpspot.x(1)-radius:g_tmpspot.x(1)+radius,i).*circle));
end

%JWK_INCOMPLETE
for i = 1:pma.len
    g_tmpspot.acceptor(i) = ...
        sum(sum(pma.frames(g_tmpspot.y(2)-radius:g_tmpspot.y(2)+radius,g_tmpspot.x(2)-radius:g_tmpspot.x(2)+radius,i).*circle));
end

    

make_tmpspotframe()
draw_frame()
draw_trace()
function set_tmpspot(i)
global g_tmpspot;
global g_spots;
global ghandles
global pma

g_tmpspot.index        = i;
g_tmpspot.x(1)         = g_spots.x(i,1);
g_tmpspot.x(2)         = g_spots.x(i,2);
g_tmpspot.y(1)         = g_spots.y(i,1);
g_tmpspot.y(2)         = g_spots.y(i,2);
g_tmpspot.donor        = g_spots.donor(i,:);
g_tmpspot.acceptor     = g_spots.acceptor(i,:);
g_tmpspot.etime        = g_spots.etime(i,:);

make_tmpspotframe()
draw_frame()
draw_trace()

%%%%%%% PRIVATE %%%%%%%
function mapped = make_mapped_coordinate(x,y)

P = zeros(4,4);
Q = zeros(4,4);

if x > 256
   %selected acceptor 
    x = x - 256-1;
    y = y-1;

    P(:) = [ ...
        -266871000000 ...
        -4023970000000 ...
        19141500000 ...
        -23513000 ...
        94802400000000 ...
        129600000000 ...
        -568250000 ...
        686717 ...
        57233800000 ...
        -1188760000 ...
        5179620 ...
        -6333.17...
        -163853000 ...
        3212910 ...
        -14029.4 ...
        17.3716 ...
        ]
    Q(:) = [ ...
        -11732800000 ...
        98463700000000 ...
        8387000000 ...
        -11872800 ...
        -74637800000 ...
        43426000000 ...
        -237806000 ...
        329130 ...
        9311360000 ...
        -354871000 ...
        1897220 ...
        -2613.1 ...
        -27862000 ...
        891458 ...
        -4743.65 ...
        6.52014000000000...
        ]

    
    
    input = [ 1 x x^2 x^3; y x*y y*x^2 y*x^3; y^2 y^2*x y^2*x^2 y^2*x^3;y^3 y^3*x y^3*x^2 y^3*x^3 ];
    mapped_x = sum(sum(input .*P / 100000000000000))+1;
    mapped_y = sum(sum(input .*Q / 100000000000000))+1;
    
else
   %selected donor
    x = x-1;
    y = y-1;
    
    P(:) = [ ...
        300526000000    4004180000000   -19059000000    23503900 ...
        105141000000000 -129242000000   567943000       -689949 ...
        -56638500000     1188130000      -5191790        6382.570000 ...
        162327000       -3215840        14081           -17.52190000]
    
    Q(:) = [ ...
        15600800000 ...
        101479000000000 ...
        -8083900000 ...
        11484600 ...
        760930000000 ...
        -42481900000 ...
        233099000 ...
        -323524 ...
        -9410120000 ...
        348777000 ...
        -1868090 ...
        2578.84 ...
        28065300 ...
        -879026 ...
        4686.480000 ...
        -6.453270000 ...
        ]

    input = [ 1 x x^2 x^3; y x*y y*x^2 y*x^3; y^2 y^2*x y^2*x^2 y^2*x^3;y^3 y^3*x y^3*x^2 y^3*x^3 ];
    mapped_x = sum(sum(input .*P / 100000000000000))+256+1;
    mapped_y = sum(sum(input .*Q / 100000000000000))+1;

end

mapped   = [mapped_x mapped_y];
function index = istherespot(x,y,handles)
global MAX_NUM_SPOTS
global g_spots
index = 0;
current_frame = get_currentframeindex();
for i=1:MAX_NUM_SPOTS     
    if (abs(g_spots.y(i,1)-y)<2 && abs(g_spots.x(i,1)-x)<2) || (abs(g_spots.y(i,2)-y)<2 && abs(g_spots.x(i,2)-x)<2)
        if current_frame > g_spots.etime(i,1) && current_frame < g_spots.etime(i,3) 
            index = i;
            display('this is occupied')
            return
        end
    end
end
index = 0;

function initialize_spots_index(index)
global g_spots
global pma


g_spots.total                           = g_spots.total-1;
cursor_x = g_spots.x(index,1); cursor_y = g_spots.y(index,1);

g_spots.x(index,:)        = [-1 -1];
g_spots.y(index,:)        = [-1 -1];
g_spots.donor(index,:)    = zeros(1,pma.len);
g_spots.acceptor(index,:) = zeros(1,pma.len);
g_spots.etime(index,:)    = pma.len*ones(1,3);

% set_tmpspot(index)
set_newtmpspot(cursor_x,cursor_y)
make_maskframe()
draw_frame()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      DRAW      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% PUBLIC %%%%%%%
function draw_pma()
function draw(type)
if strcmp(type,'slider')
    draw_frame(); make_maskframe(); make_tmpspotframe();
elseif strcmp(type, 'initial')
    draw_frame();
else
    make_maskframe(); make_tmpspotframe();
    draw_frame(); draw_trace(); draw_hist();
end
function draw_hist()
draw_hist_fusion()
draw_hist_docking()

function draw_hist_fusion()
global ghandles
global g_spots
global pma
global NUM_BIN
axes(ghandles.windowHIST);
% docked_index   = find(g_spots.etime(:,1) ~= pma.len);
fused_index    = find(g_spots.etime(:,3) ~= pma.len);

if ~isempty(fused_index)
    
    fuse_time  = g_spots.etime(fused_index,3)-g_spots.etime(fused_index,1);
    
    max_fuse_time = ceil(max(fuse_time)/NUM_BIN)*NUM_BIN;
    if max_fuse_time == 0
        max_fuse_time = 5;
    end
    xbins = 0:NUM_BIN:max_fuse_time;
    xbins = xbins + NUM_BIN/2;
    [hist_data b] = hist(fuse_time,xbins);
    bar(b,hist_data);
    axis([0 max_fuse_time 0 max(hist_data)*1.2])
    hold on;
    time_const = median(fuse_time);
    line([time_const time_const], [0 max(hist_data)*1.2],'LineWidth', 1, 'Color', 'r');

    hold off;
    grid on;
end

function draw_hist_docking()
global ghandles
global g_spots
global pma
global NUM_BIN
axes(ghandles.windowHIST2);
dock_time      = g_spots.etime(g_spots.etime(:,1) ~= pma.len,1);

max_dock_time = ceil(max(dock_time)/NUM_BIN)*NUM_BIN;
xbins = 0:NUM_BIN:max_dock_time;


[hist_data b] = hist(dock_time,xbins);

cum_hist_data = cumsum(hist_data);


plot(b,cum_hist_data);
axis([0 max_dock_time 0 max(cum_hist_data)*1.2])
% bar(b,cum_hist_data);
grid on;


function draw_trace()
global ghandles
global g_tmpspot
global g_spots
global FRET
global FRET_MIN
global FRAME_MAX
cla(ghandles.windowTRACE);
axes(ghandles.windowTRACE);

plot_acceptor   = g_tmpspot.acceptor;
plot_donor      = g_tmpspot.donor;
hold on 
% axis([0 pma.len 0 max(g_tmpspot.donor)*1.2])

current_frame = get_currentframeindex();
max_trace = max([max(plot_acceptor) max(plot_donor)]);
% if g_tmpspot.etime(1)> 0
line([g_tmpspot.etime(1) g_tmpspot.etime(1)], [0 max_trace],'LineWidth', 4, 'Color', 'k'); hold on;
line([g_tmpspot.etime(2) g_tmpspot.etime(2)], [0 max_trace],'LineWidth', 4, 'Color', 'b'); hold on;
line([g_tmpspot.etime(3) g_tmpspot.etime(3)], [0 max_trace],'LineWidth', 4, 'Color', 'r'); hold on;
plot(plot_acceptor, 'Parent', ghandles.windowTRACE,'LineWidth', 1, 'Color', 'r');
plot(plot_donor, 'Parent', ghandles.windowTRACE,'LineWidth', 1, 'Color', 'g'); hold off;
% end

%JWK NOW
xlim(ghandles.windowTRACE,[1 numel(g_tmpspot.acceptor)]);
ylim(ghandles.windowTRACE,[0  max(max(g_tmpspot.donor,g_tmpspot.acceptor))*1.2]);
set(ghandles.windowTRACE, 'yTick', 0:25:ceil(max(max(g_tmpspot.donor,g_tmpspot.acceptor))*1.2) );
hold off

function draw_frame()
global ghandles
global g_maskframe
global g_tmpspotframe
       
frame_index = str2double(get(ghandles.frame_index, 'String'));
%frame_index = get(ghandles.frame_index,'Value');
frame = get_frame(frame_index);
frame(g_maskframe == 1) =  50;
frame(g_maskframe == 5) = 150;

frame(g_tmpspotframe == 1) = 150;
frame(g_tmpspotframe == 5) = 200;

image(frame, 'Parent', ghandles.windowPMA);
set(ghandles.windowPMA, 'xTick', 0, 'yTick', 0);

%%%%%%% PRIVATE %%%%%%%
function current_index = get_currentframeindex()
global ghandles
current_index = str2double(get(ghandles.frame_index, 'String'));
function set_fusiontime(index,fusion_type)            
global g_spots
global g_tmpspot
global ghandles
current_index = str2double(get(ghandles.frame_index, 'String'));
if fusion_type == 1
    g_spots.etime(index,2) = current_index;
    g_tmpspot.etime(2) = current_index;
elseif fusion_type == 2
    g_spots.etime(index,3) = current_index;
    g_tmpspot.etime(3) = current_index;
elseif fusion_type == 3
    g_spots.etime(index,2) = current_index;
    g_spots.etime(index,3) = current_index;
    g_tmpspot.etime(2) = current_index;
    g_tmpspot.etime(3) = current_index;
end

function frame = get_frame(index)
%%%%%%%%%%%%%%% loading frames from pma file %%%%%%%%
global pma;
frame = pma.frames(:,:,index);
function cal_fret()
global g_spots
global FRET
global pma
global MAX_NUM_SPOTS


plot_donor    = g_spots.donor(1:g_spots.total,:)-g_spots.donor_level;
plot_acceptor = g_spots.acceptor(1:g_spots.total,:)-g_spots.acceptor_level;

FRET.value    = zeros(MAX_NUM_SPOTS,pma.len);
total         = plot_donor + plot_acceptor;
FRET.value(1:g_spots.total,:)  = plot_acceptor ./ total;
g_spots.fret = FRET.value;

function make_maskframe()
global g_spots;
global g_maskframe;
global pma;
radius = 4;

circle_D = ... 
  [ 1 1 1 0 0 0 1 1 1;...
    1 1 0 0 0 0 0 1 1;...
    1 0 0 0 0 0 0 0 1;...
    0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0;...
    1 0 0 0 0 0 0 0 1;...
    1 1 0 0 0 0 0 1 1;...
    1 1 1 0 0 0 1 1 1];

current_frame = get_currentframeindex();
g_maskframe = zeros(pma.height,pma.width);
assignin('base','g_maskframe_previous',g_maskframe)
%jwk_fix
% spot_index = find(g_spots.x(:,1) > 0);

docked_index   = find(g_spots.etime(:,1) <= current_frame);
fused_index = find(g_spots.etime(:,3) < current_frame);
[x y] = size(g_spots.etime);
A = zeros(x,1);
A(docked_index) = 1;
A(fused_index) = 0;
spot_index = find(A>0);

%jwk 20170614
if ~isempty(spot_index)
    for i = spot_index'
        if ((round(g_spots.y(i,2))-radius) * (round(g_spots.y(i,2))+radius) * (round(g_spots.x(i,2))-radius) * (round(g_spots.x(i,2))+radius)) < 0
            display('asdfasdfsadfsdafsdafasdfsdafsdafsdafdfsaasdf')
            spot_index(i)
        else
        g_maskframe(round(g_spots.y(i,2))-radius : round(g_spots.y(i,2))+radius, round(g_spots.x(i,2))-radius : round(g_spots.x(i,2))+radius) ...
            = g_maskframe(round(g_spots.y(i,2))-radius : round(g_spots.y(i,2))+radius, round(g_spots.x(i,2))-radius : round(g_spots.x(i,2))+radius) |circle_D;
        end
    end
    g_maskframe = g_maskframe * 5;
    for i = spot_index'
        if ((round(g_spots.y(i,1))-radius) * (round(g_spots.y(i,1))+radius) * (round(g_spots.x(i,1))-radius) * (round(g_spots.x(i,1))+radius)) < 0
            display('LEFTQ!!!!!!')
            i
            (round(g_spots.y(i,1))-radius)  
            (round(g_spots.y(i,1))+radius)  
            (round(g_spots.x(i,1))-radius)  
            (round(g_spots.x(i,1))+radius)
        else
            g_maskframe(round(g_spots.y(i,1))-radius : round(g_spots.y(i,1))+radius, round(g_spots.x(i,1))-radius : round(g_spots.x(i,1))+radius) ...
                = g_maskframe(round(g_spots.y(i,1))-radius : round(g_spots.y(i,1))+radius, round(g_spots.x(i,1))-radius : round(g_spots.x(i,1))+radius) |circle_D;
        end
    end
end
function make_tmpspotframe()
global g_tmpspotframe;
global g_tmpspot
global pma;


x = g_tmpspot.x(1);
y = g_tmpspot.y(1);


%JWK_INCOMPLETE must check if mapped x y are in boundary
if x < 1
    g_tmpspotframe = zeros(pma.height,pma.width);
    return
end

radius = 3;
circle = ... 
  [ 0 0 1 1 1 0 0;...
    0 1 0 0 0 1 0;...
    1 0 0 0 0 0 1;...
    1 0 0 0 0 0 1;...
    1 0 0 0 0 0 1;...
    0 1 0 0 0 1 0
    0 0 1 1 1 0 0];

g_tmpspotframe = zeros(pma.height,pma.width);
g_tmpspotframe(round(y)-radius : round(y)+radius, round(x)-radius : round(x)+radius) = circle;
x = g_tmpspot.x(2);
y = g_tmpspot.y(2);
g_tmpspotframe(round(y)-radius : round(y)+radius, round(x)-radius : round(x)+radius) = 5 * circle;

function slider_windowCONTRAST_Callback(hObject, eventdata, handles)
function slider_windowCONTRAST_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in chkGrid.
function chkGrid_Callback(hObject, eventdata, handles)
% hObject    handle to chkGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Grid = [handles.gridline_horizontal_bottom, handles.gridline_horizontal_top, handles.gridline_vertical];
status = get(handles.chkGrid, 'Value');
if status == 1
    set(Grid, 'Visible', 'On')
end
if status == 0
    set(Grid, 'Visible', 'Off')
end    
        
% --- Executes on button press in RB_gridline_vertical.
function RB_gridline_vertical_Callback(hObject, eventdata, handles)
% hObject    handle to RB_gridline_vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

status = get(handles.RB_gridline_vertical, 'Value');
if status == 1
    set(handles.gridline_vertical, 'Visible', 'On')
end
if status == 0
    set(handles.gridline_vertical, 'Visible', 'Off')
end    


% --- Executes on button press in RB_gridline_horizontal_bottom.
function RB_gridline_horizontal_bottom_Callback(hObject, eventdata, handles)
% hObject    handle to RB_gridline_horizontal_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

status = get(handles.RB_gridline_horizontal_bottom, 'Value');
if status == 1
    set(handles.gridline_horizontal_bottom, 'Visible', 'On')
end
if status == 0
    set(handles.gridline_horizontal_bottom, 'Visible', 'Off')
end 


% --- Executes on button press in RB_gridline_horizontal_top.
function RB_gridline_horizontal_top_Callback(hObject, eventdata, handles)
% hObject    handle to RB_gridline_horizontal_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

status = get(handles.RB_gridline_horizontal_top, 'Value');
if status == 1
    set(handles.gridline_horizontal_top, 'Visible', 'On')
end
if status == 0
    set(handles.gridline_horizontal_top, 'Visible', 'Off')
end 
