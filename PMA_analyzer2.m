% PMA_analyzer2
% Joseph Zambreno, Iowa State University
% 2/10/18

function varargout = PMA_analyzer2(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PMA_analyzer2_OpeningFcn, ...
                   'gui_OutputFcn',  @PMA_analyzer2_OutputFcn, ...
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


% --- Executes just before PMA_analyzer2 is made visible.
function PMA_analyzer2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PMA_analyzer2 (see VARARGIN)

% Choose default command line output for PMA_analyzer2
handles.output = hObject;

handles.timer = timer(...
    'ExecutionMode', 'fixedRate', ...       % Run timer repeatedly.
    'Period', 0.03, ...                        % Attempt a real-time display.
    'TimerFcn', {@timer_callback}); % Specify callback function.

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PMA_analyzer2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global pmaData;
global fileData;
global windowData;
global appState;
global ghandles;
ghandles = handles;

appState.pma_loaded = 0;
appState.events_analyzed = 0;
appState.display_contrast = 127;
appState.display_threshold = 12;
appState.black_and_white = 0;
appState.display_frame = 1;
appState.current_point = [512, 512];
appState.saved_zoom = zeros(2,2);
appState.prospective_zoom = zeros(1,2);
appState.prospective_clickcount = zeros(1,2);
appState.filter_min = 1000;
appState.filter_max = 2500;

appState.eventtypes = {'unknown', 'docking', 'instant fusion', 'delayed fusion'}; 
set(ghandles.prospective1_eventtype_pulldown, 'String', appState.eventtypes);
set(ghandles.saved1_eventtype_pulldown, 'String', appState.eventtypes);
set(ghandles.saved2_eventtype_pulldown, 'String', appState.eventtypes);


% Create a callback function to track the pma_axes current position
set(ghandles.figure1, 'WindowButtonMotionFcn', @mouseMove);

% Create a callback function for click events
set(ghandles.figure1, 'WindowButtonDownFcn', @mouseClick);


% Create a close callback
set(ghandles.figure1, 'CloseRequestFcn',@close_request);

% Configure play/pause button
appState.play_pause = 1; % Start off paused, I think
set(ghandles.play_pause,'Enable','off');
set(ghandles.play_pause,'String',char(hex2dec('25BA')));
set(ghandles.play_pause,'Enable','on');

% Configure current point and filename text
set(ghandles.current_point,'String','');
set(ghandles.filename_text,'String','');

% --- Outputs from this function are returned to the command line.
function varargout = PMA_analyzer2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
return;


function mouseClick(obj,evt)
    global ghandles;
    global appState;
    global pmaData;
    global eventData;
    
    if (appState.pma_loaded == 0)
        return;
    end;

    % Differentiate between single and double-click events
    persistent chk
    persistent clickPoint
    if isempty(chk)
      chk = 1;
      clickPoint{1} = get(ghandles.pma_axes,'CurrentPoint');
      clickPoint{2} = get(ghandles.saved1_axes1,'CurrentPoint');
      clickPoint{3} = get(ghandles.saved1_axes2,'CurrentPoint');
      clickPoint{4} = get(ghandles.saved2_axes1,'CurrentPoint');
      clickPoint{5} = get(ghandles.saved2_axes2,'CurrentPoint');
      clickPoint{6} = get(ghandles.prospective1_axes1,'CurrentPoint');
      clickPoint{7} = get(ghandles.prospective1_axes2,'CurrentPoint');   
      pause(0.25); %Add a delay to distinguish single click from a double click
      if chk == 1
        clickType = 'Single';          
        chk = [];
      end
    else
      chk = [];
      clickType = 'Double';
    end
    
    % This ignores the second click of a double-click
    if ~exist('clickType','var')
        return;
    end;
    
    % Check if we've clicked inside the pma_axes
    currentPoint = clickPoint{1};
    if ((currentPoint(1,1) >= 0) & (currentPoint(1,1) < pmaData.width) & ...
             (currentPoint(1,2) >= 0) & (currentPoint(1,2) < pmaData.height))

         % Add a prospective event at that location      
        candidate_event.start_frame = appState.display_frame;
        candidate_event.end_frame = appState.display_frame+20;
        if (candidate_event.end_frame > pmaData.numframes)
            candidate_event.end_frame = pmaData.numframes;
        end;
        j = round(currentPoint(1,1))+1;
        i = round(currentPoint(1,2))+1;
        candidate_event.center = [i, j];
        neighborhood_radius = 4;
        neighborhood_radius = saturate_radius(neighborhood_radius, j, i, pmaData.width, pmaData.height);
        candidate_event.radii = neighborhood_radius*ones(candidate_event.end_frame-candidate_event.start_frame+1,1);
        candidate_event.neighborhood_intensity = 0;
        candidate_event.source = 'User';
        candidate_event.eventtype = 1;

        for k3=1:pmaData.numframes
          candidate_event.neighborhood_intensity(k3) = mean(mean(pmaData.frames(i-neighborhood_radius:i+neighborhood_radius,j-neighborhood_radius:j+neighborhood_radius,k3)));
        end;
        
        if (eventData.num_prospective_events == 0)
            eventData.prospective_events = candidate_event;
        else
            % Just insert at the end for now
            eventData.prospective_events = [eventData.prospective_events candidate_event];
        end;
        eventData.num_prospective_events = eventData.num_prospective_events + 1;

        % We want to always see the most recent event anyway, right?
        appState.current_prospective_event = eventData.num_prospective_events;

        display_prospective_frame();
        display_pma_frame();
        return;
    end;

    
    
    
    if (appState.events_analyzed == 0)
        return;
    end;
    
    % Only zoom on double-click events
    if strcmp(clickType,'Double')
        currentPoint = get(ghandles.saved1_axes1,'CurrentPoint');
        XLim = get(ghandles.saved1_axes1,'XLim');
        YLim = get(ghandles.saved1_axes1,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.saved_zoom(1,1) = 1 - appState.saved_zoom(1,1);
            display_saved_frame();
            return;
        end;

        currentPoint = get(ghandles.saved1_axes2,'CurrentPoint');
        XLim = get(ghandles.saved1_axes2,'XLim');
        YLim = get(ghandles.saved1_axes2,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
             (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.saved_zoom(1,2) = 1 - appState.saved_zoom(1,2);
            display_saved_frame();
            return;
        end;

        currentPoint = get(ghandles.saved2_axes1,'CurrentPoint');
        XLim = get(ghandles.saved2_axes1,'XLim');
        YLim = get(ghandles.saved2_axes1,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.saved_zoom(2,1) = 1 - appState.saved_zoom(2,1);
            display_saved_frame();
            return;
        end;

        currentPoint = get(ghandles.saved2_axes2,'CurrentPoint');
        XLim = get(ghandles.saved2_axes2,'XLim');
        YLim = get(ghandles.saved2_axes2,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
             (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.saved_zoom(2,2) = 1 - appState.saved_zoom(2,2);
            display_saved_frame();
            return;
        end;
    else
        currentPoint = get(ghandles.saved1_axes1,'CurrentPoint');
        current_start_frame = eventData.candidate_events(appState.current_saved_event).start_frame;       
        XLim = get(ghandles.saved1_axes1,'XLim');
        YLim = get(ghandles.saved1_axes1,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.display_frame = current_start_frame;                     

            set(ghandles.current_frame_slider,'Value',appState.display_frame);
            set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
            display_pma_frame();
            return;
        end;

        currentPoint = get(ghandles.saved1_axes2,'CurrentPoint');
        current_start_frame = eventData.candidate_events(appState.current_saved_event).start_frame;       
        XLim = get(ghandles.saved1_axes2,'XLim');
        YLim = get(ghandles.saved1_axes2,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
             (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.display_frame = current_start_frame;           
            set(ghandles.current_frame_slider,'Value',appState.display_frame);
            set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
            display_pma_frame();
            return;
            
        end;

        currentPoint = get(ghandles.saved2_axes1,'CurrentPoint');
        current_start_frame = eventData.candidate_events(appState.current_saved_event+1).start_frame;       
        XLim = get(ghandles.saved2_axes1,'XLim');
        YLim = get(ghandles.saved2_axes1,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.display_frame = current_start_frame;           
            set(ghandles.current_frame_slider,'Value',appState.display_frame);
            set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
            display_pma_frame();
            return;
            
        end;

        currentPoint = get(ghandles.saved2_axes2,'CurrentPoint');
        current_start_frame = eventData.candidate_events(appState.current_saved_event+1).start_frame;
        XLim = get(ghandles.saved2_axes2,'XLim');
        YLim = get(ghandles.saved2_axes2,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
             (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.display_frame = current_start_frame;
            set(ghandles.current_frame_slider,'Value',appState.display_frame);
            set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
            display_pma_frame();
            return;
            
        end;        
        
        
        
    end;
    
    if (eventData.num_prospective_events == 0)
        return;
    end;

    % Only zoom the prospective events on a double-click. 
    % Otherwise, adjust the prospective_event startframe and endframe
    if strcmp(clickType,'Double')
        currentPoint = get(ghandles.prospective1_axes1,'CurrentPoint');
        XLim = get(ghandles.prospective1_axes1,'XLim');
        YLim = get(ghandles.prospective1_axes1,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.prospective_zoom(1) = 1 - appState.prospective_zoom(1);
            display_prospective_frame();
            return;
            
        end;

        currentPoint = get(ghandles.prospective1_axes2,'CurrentPoint');
        XLim = get(ghandles.prospective1_axes2,'XLim');
        YLim = get(ghandles.prospective1_axes2,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            appState.prospective_zoom(2) = 1 - appState.prospective_zoom(2);
            display_prospective_frame();
            return;
            
        end;
    
    % Single click - move start / end points
    else
        currentPoint = clickPoint{6};
        XLim = get(ghandles.prospective1_axes1,'XLim');
        YLim = get(ghandles.prospective1_axes1,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            current_start_frame = eventData.prospective_events(appState.current_prospective_event).start_frame;
            current_end_frame = eventData.prospective_events(appState.current_prospective_event).end_frame;

            % Move start point
            if (appState.prospective_clickcount(1) == 0)
                if (round(currentPoint(1,1)) < current_end_frame)
                    eventData.prospective_events(appState.current_prospective_event).start_frame = round(currentPoint(1,1));
                    appState.prospective_clickcount(1) = 1 - appState.prospective_clickcount(1);
                    display_prospective_frame();
                end;
            % Move end point
            else
                if (round(currentPoint(1,1)) > current_start_frame)
                    eventData.prospective_events(appState.current_prospective_event).end_frame = round(currentPoint(1,1));
                    appState.prospective_clickcount(1) = 1 - appState.prospective_clickcount(1);
                    display_prospective_frame();
                end;         
            end;
            % Have to resize the radii array
            eventData.prospective_events(appState.current_prospective_event).radii = eventData.prospective_events(appState.current_prospective_event).radii(1)*ones(eventData.prospective_events(appState.current_prospective_event).end_frame - eventData.prospective_events(appState.current_prospective_event).start_frame+1,1);
            
            % No matter what, move the pma_display to the current selected
            % event
            current_start_frame = eventData.prospective_events(appState.current_prospective_event).start_frame;
            appState.display_frame = current_start_frame;
            set(ghandles.current_frame_slider,'Value',appState.display_frame);
            set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
            display_pma_frame();
            return;        
            
        end;

        currentPoint = clickPoint{7};
        XLim = get(ghandles.prospective1_axes2,'XLim');
        YLim = get(ghandles.prospective1_axes2,'YLim');
        if ((currentPoint(1,1) >= XLim(1)) & (currentPoint(1,1) < XLim(2)) & ...
                 (currentPoint(1,2) >= YLim(1)) & (currentPoint(1,2) < YLim(2)))
            current_start_frame = eventData.prospective_events(appState.current_prospective_event).start_frame;
            current_end_frame = eventData.prospective_events(appState.current_prospective_event).end_frame;

            % Move start point
            if (appState.prospective_clickcount(1) == 0)
                if (round(currentPoint(1,1)) < current_end_frame)
                    eventData.prospective_events(appState.current_prospective_event).start_frame = round(currentPoint(1,1));
                    appState.prospective_clickcount(1) = 1 - appState.prospective_clickcount(1);
                    display_prospective_frame();
                end;
            % Move end point
            else
                if (round(currentPoint(1,1)) > current_start_frame)
                    eventData.prospective_events(appState.current_prospective_event).end_frame = round(currentPoint(1,1));
                    appState.prospective_clickcount(1) = 1 - appState.prospective_clickcount(1);
                    display_prospective_frame();
                end;
            
            end;
            % Have to resize the radii array
            eventData.prospective_events(appState.current_prospective_event).radii = eventData.prospective_events(appState.current_prospective_event).radii(1)*ones(eventData.prospective_events(appState.current_prospective_event).end_frame - eventData.prospective_events(appState.current_prospective_event).start_frame+1,1);
            
            % No matter what, move the pma_display to the current selected
            % event
            current_start_frame = eventData.prospective_events(appState.current_prospective_event).start_frame;
            appState.display_frame = current_start_frame;
            set(ghandles.current_frame_slider,'Value',appState.display_frame);
            set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
            display_pma_frame();
            return;         
            
         end;
        
    end;
        
    
    
    
function mouseMove(varargin) 
    global ghandles;
    global appState;
    global pmaData;
    
    if (appState.pma_loaded == 0)
        return;
    end;
        
     currentPoint = get(ghandles.pma_axes,'CurrentPoint');

     % Only update the text if we are inside the axes we care about
     if ((currentPoint(1,1) >= 0) & (currentPoint(1,1) < pmaData.width) & ...
             (currentPoint(1,2) >= 0) & (currentPoint(1,2) < pmaData.height))
         set(ghandles.current_point,'String',sprintf('Current Point: [%d, %d]',...
             round(currentPoint(1,2))+1, round(currentPoint(1,1))+1));
     else
         set(ghandles.current_point,'String', '');
     end;

     
    if (appState.events_analyzed == 0)
        return;
    end;

    
return;


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
    global fileData;
    global pmaData;
    global ghandles;
    global appState;
    global eventData;

    if (appState.events_analyzed == 0)
        return;
    end;

    fileData.CSVfileName = [fileData.shortfileName '.csv'];
    filename = fullfile(fileData.pathName, fileData.CSVfileName);

    h3 = waitbar(0,'Saving event 1','Name','Saving analysis data...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
        'WindowStyle', 'modal');
    setappdata(h3,'canceling',0)      

    fid = fopen(filename, 'w');
    
    fprintf(fid, '%s\n', datestr(datetime('now')));
    M1 = {'EventID', 'CenterX', 'CenterY', 'FrameStart', 'FrameEnd', 'minRadius', 'maxRadius', 'Source', 'Type'}; 
    for i=1:length(M1)-1
        fprintf(fid, '%s, ', M1{i});
    end;
    fprintf(fid, '%s\n', M1{length(M1)});
    
    for i=1:eventData.num_candidate_events
        if getappdata(h3,'canceling')            
            break
        end      
               
        if (mod(i,10) == 0)
            waitbar(i/eventData.num_candidate_events,h3,sprintf('Saving event %d of %d', i, eventData.num_candidate_events));
        end;

        candidate_event = eventData.candidate_events(i);

        fprintf(fid, '%d, %d, %d, %d, %d, %d, %d, %s, %s\n', i, ...
                                    candidate_event.center(1), ...
                                    candidate_event.center(2), ...
                                    candidate_event.start_frame, ...
                                    candidate_event.end_frame, ...
                                    min(candidate_event.radii), ...
                                    max(candidate_event.radii), ...
                                    candidate_event.source, ...
                                    appState.eventtypes{candidate_event.eventtype});
    
    end;
    fclose(fid);        

    
    fileData.MATfileName = [fileData.shortfileName '.mat'];
    filename = fullfile(fileData.pathName, fileData.MATfileName);
    waitbar(1,h3,sprintf('Saving workspace variables...'));
    warning('off', 'MATLAB:Figure:FigureSavedToMATFile');
    save(filename, '-v7.3', 'fileData', 'ghandles', 'appState', 'eventData');
    
    delete(h3);    


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fileData;
global pmaData;
global ghandles;
global appState;
global eventData;

[fileData.fileName,fileData.pathName] = uigetfile('*.pma','Select the PMA data file');
fileData.shortfileName = fileData.fileName(1:length(fileData.fileName)-4);
filename = fullfile(fileData.pathName, fileData.fileName);
if (~isequal(fileData.fileName,0) & ~isempty(dir(filename)))

    appState.filter_min = str2num(get(ghandles.filter_edit1,'String'));
    appState.filter_max = str2num(get(ghandles.filter_edit2,'String'));
    
    load_pmafile_with_progress_bar(filename);
    
    
    % now that we know the number of frames, clean up the filter range
    if (appState.filter_max > pmaData.numframes)
        appState.filter_max = pmaData.numframes;
    end;
    if (appState.filter_min >= appState.filter_max)
        appState.filter_min = appState.filter_max - 1;
    end;
    set(ghandles.filter_edit1,'String',sprintf('%d', appState.filter_min));
    set(ghandles.filter_edit2,'String',sprintf('%d', appState.filter_max));

    set(ghandles.filter_edit1,'Enable','off');
    set(ghandles.filter_edit2,'Enable','off');
    
    appState.display_frame = appState.filter_min;    
    appState.frame_range = appState.filter_max-appState.filter_min+1;
    
    set(ghandles.total_frames,'String',sprintf('/ %d', pmaData.numframes));
    set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
    set(ghandles.contrast_threshold_edit,'String',sprintf('%d', appState.display_contrast));
    
    set(ghandles.current_frame_slider,'Value',appState.display_frame);
    set(ghandles.current_frame_slider,'Min', 1);
    set(ghandles.current_frame_slider,'Max', pmaData.numframes);
    set(ghandles.current_frame_slider,'SliderStep',[1/pmaData.numframes 10/pmaData.numframes]);

    set(ghandles.filename_text,'String',fileData.shortfileName);
    
    display_pma_frame();
    
    analyze_events_with_progress_bar();    
    appState.current_saved_event = 1;
    appState.current_prospective_event = 0;
    display_saved_frame();
    
end

%---------------------------------------------------------------------
function timer_callback(~,~)
    global appState;
    global pmaData;
    global ghandles;

    
    % If we're running, make sure to update
    if (appState.play_pause == 0)
        display_pma_frame();

        % Stop running if we've reached the end
        appState.display_frame = appState.display_frame + 1;
        if (appState.display_frame > appState.filter_max)
            appState.display_frame = appState.filter_max;            
            set(ghandles.play_pause,'Enable','off');
            set(ghandles.play_pause,'String',char(hex2dec('25C0'))); % Rewind button
            set(ghandles.play_pause,'Enable','on');
            appState.play_pause = 1;
            stop(ghandles.timer);
        end;

        % Test to see how updating the frame count + scroll looks
        set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
        set(ghandles.current_frame_slider,'Value',appState.display_frame);
        
    end;
        
%---------------------------------------------------------------------    
function display_prospective_frame()    
    global ghandles;
    global pmaData;
    global appState;
    global eventData;

    if (appState.events_analyzed == 0)
        return;
    end;

    % Display the prospective event
    visibility = 'off';
    if (eventData.num_prospective_events >= 1)
        visibility = 'on';
    end;

    set(ghandles.prospective1_axes1,'visible',visibility);
    set(get(ghandles.prospective1_axes1,'children'),'visible',visibility)   
    set(ghandles.prospective1_axes2,'visible',visibility);
    set(get(ghandles.prospective1_axes2,'children'),'visible',visibility)   
    set(ghandles.prospective1_eventnum_text,'visible',visibility);
    set(ghandles.prospective1_location_text,'visible',visibility);
    set(ghandles.prospective1_frame_text,'visible',visibility);
    set(ghandles.prospective1_radius_text,'visible',visibility);
    set(ghandles.prospective1_source_text,'visible',visibility);
    set(ghandles.prospective1_eventtype_pulldown,'visible',visibility);
    set(ghandles.prospective1_addevent_pushbutton,'visible',visibility);
    set(ghandles.prospective1_addevent_text,'visible',visibility);        
    set(ghandles.prospective1_ignoreevent_pushbutton,'visible',visibility);
    set(ghandles.prospective1_ignoreevent_text,'visible',visibility);        

    if (eventData.num_prospective_events == 1)
        visibility = 'off';
    end;
    set(ghandles.prospective_slider,'visible',visibility);        

    
    if (eventData.num_prospective_events >= 1)
        prospective1_event = appState.current_prospective_event;
        candidate_event = eventData.prospective_events(prospective1_event);
        
        set(ghandles.prospective1_eventnum_text,'String',sprintf('Event: %d of %d', prospective1_event, eventData.num_prospective_events));
        set(ghandles.prospective1_location_text,'String',sprintf('Location: [%d, %d]', candidate_event.center(1), candidate_event.center(2)));
        set(ghandles.prospective1_frame_text,'String',sprintf('Frames: [%d, %d]', candidate_event.start_frame, candidate_event.end_frame));
        set(ghandles.prospective1_radius_text,'String',sprintf('Radius: [%d - %d]', min(candidate_event.radii), max(candidate_event.radii)));
        set(ghandles.prospective1_source_text,'String',['Source: ', candidate_event.source]);
        set(ghandles.prospective1_eventtype_pulldown,'Value', candidate_event.eventtype);      

        set(ghandles.prospective_slider,'Max', eventData.num_prospective_events);
        set(ghandles.prospective_slider,'Value',eventData.num_prospective_events-prospective1_event+1);
        set(ghandles.prospective_slider,'Min', 1);
        set(ghandles.prospective_slider,'SliderStep',[1/eventData.num_prospective_events 2/eventData.num_prospective_events]);     
        
        axes(ghandles.prospective1_axes1);
        cla;
        hold on;
        set(ghandles.prospective1_axes1,'fontsize',6);
        title('Center Intesity');
        xdata = 1:pmaData.numframes;
        plot_data = permute(pmaData.frames(candidate_event.center(1), candidate_event.center(2),1:pmaData.numframes),[3 2 1]);
        if (appState.prospective_zoom(1) == 1)
            title('Center Intensity (zoomed)');
            xdata = candidate_event.start_frame:candidate_event.end_frame;
            plot_data2 = plot_data(candidate_event.start_frame:candidate_event.end_frame);
            plot_data = plot_data2;
        end;
        plot(xdata, plot_data);
        plot([candidate_event.start_frame candidate_event.start_frame], [min(plot_data) max(plot_data)], 'g'); 
        plot([candidate_event.end_frame candidate_event.end_frame], [min(plot_data) max(plot_data)], 'r'); 

        
        axes(ghandles.prospective1_axes2);
        cla;
        hold on;
        set(ghandles.prospective1_axes2,'fontsize',6);
        title('Neighborhood Intensity');
        xdata = 1:pmaData.numframes;
        plot_data = candidate_event.neighborhood_intensity;
        if (appState.prospective_zoom(2) == 1)
            title('Neighborhood Intensity (zoomed)');
            xdata = candidate_event.start_frame:candidate_event.end_frame;
            plot_data2 = plot_data(candidate_event.start_frame:candidate_event.end_frame);
            plot_data = plot_data2;
        end;
        plot(xdata, plot_data);
        plot([candidate_event.start_frame candidate_event.start_frame], [min(plot_data) max(plot_data)], 'g'); 
        plot([candidate_event.end_frame candidate_event.end_frame], [min(plot_data) max(plot_data)], 'r'); 

    drawnow;        
    end;    
    
    
%---------------------------------------------------------------------    
function display_saved_frame()    
    global ghandles;
    global pmaData;
    global appState;
    global eventData;

    if (appState.events_analyzed == 0)
        return;
    end;
   
    % Display the first saved event
    visibility = 'off';
    if (eventData.num_candidate_events >= 1)
        visibility = 'on';
    end;
    set(ghandles.saved1_axes1,'visible',visibility);
    set(get(ghandles.saved1_axes1,'children'),'visible',visibility)   
    set(ghandles.saved1_axes2,'visible',visibility);
    set(get(ghandles.saved1_axes2,'children'),'visible',visibility)   
    set(ghandles.saved1_eventnum_text,'visible',visibility);
    set(ghandles.saved1_location_text,'visible',visibility);
    set(ghandles.saved1_frame_text,'visible',visibility);
    set(ghandles.saved1_radius_text,'visible',visibility);
    set(ghandles.saved1_source_text,'visible',visibility);
    set(ghandles.saved1_eventtype_pulldown,'visible',visibility);
    set(ghandles.saved1_removeevent_pushbutton,'visible',visibility);
    set(ghandles.saved1_removeevent_text,'visible',visibility);        
    
    % Display the second saved event and scroll bar
    visibility = 'off';
    if (eventData.num_candidate_events >= 2)
        visibility = 'on';
    end;

    set(ghandles.saved2_axes1,'visible',visibility);
    set(get(ghandles.saved2_axes1,'children'),'visible',visibility)   
    set(ghandles.saved2_axes2,'visible',visibility);
    set(get(ghandles.saved2_axes2,'children'),'visible',visibility)   
    set(ghandles.saved2_eventnum_text,'visible',visibility);
    set(ghandles.saved2_location_text,'visible',visibility);
    set(ghandles.saved2_frame_text,'visible',visibility);
    set(ghandles.saved2_radius_text,'visible',visibility);
    set(ghandles.saved2_source_text,'visible',visibility);
    set(ghandles.saved2_eventtype_pulldown,'visible',visibility);
    set(ghandles.saved2_removeevent_pushbutton,'visible',visibility);
    set(ghandles.saved2_removeevent_text,'visible',visibility);        
    set(ghandles.saved_slider,'visible',visibility);        
    
       
    if (eventData.num_candidate_events >= 1)
        saved1_event = appState.current_saved_event;
        candidate_event = eventData.candidate_events(saved1_event);
        
        set(ghandles.saved1_eventnum_text,'String',sprintf('Event: %d of %d', saved1_event, eventData.num_candidate_events));
        set(ghandles.saved1_location_text,'String',sprintf('Location: [%d, %d]', candidate_event.center(1), candidate_event.center(2)));
        set(ghandles.saved1_frame_text,'String',sprintf('Frames: [%d, %d]', candidate_event.start_frame, candidate_event.end_frame));
        set(ghandles.saved1_radius_text,'String',sprintf('Radius: [%d - %d]', min(candidate_event.radii), max(candidate_event.radii)));
        set(ghandles.saved1_source_text,'String',['Source: ', candidate_event.source]);
        set(ghandles.saved1_eventtype_pulldown,'Value', candidate_event.eventtype);      

        axes(ghandles.saved1_axes1);
        cla;
        hold on;
        set(ghandles.saved1_axes1,'fontsize',6);
        title('Center Intesity');
        xdata = 1:pmaData.numframes;
        plot_data = permute(pmaData.frames(candidate_event.center(1), candidate_event.center(2),1:pmaData.numframes),[3 2 1]);
        if (appState.saved_zoom(1,1) == 1)
            title('Center Intensity (zoomed)');
            xdata = candidate_event.start_frame:candidate_event.end_frame;
            plot_data2 = plot_data(candidate_event.start_frame:candidate_event.end_frame);
            plot_data = plot_data2;
        end;
        plot(xdata, plot_data);
        plot([candidate_event.start_frame candidate_event.start_frame], [min(plot_data) max(plot_data)], 'k'); 
        plot([candidate_event.end_frame candidate_event.end_frame], [min(plot_data) max(plot_data)], 'k'); 

        
        axes(ghandles.saved1_axes2);
        cla;
        hold on;
        set(ghandles.saved1_axes2,'fontsize',6);
        title('Neighborhood Intensity');
        xdata = 1:pmaData.numframes;
        plot_data = candidate_event.neighborhood_intensity;
        if (appState.saved_zoom(1,2) == 1)
            title('Neighborhood Intensity (zoomed)');
            xdata = candidate_event.start_frame:candidate_event.end_frame;
            plot_data2 = plot_data(candidate_event.start_frame:candidate_event.end_frame);
            plot_data = plot_data2;
        end;
        plot(xdata, plot_data);
        plot([candidate_event.start_frame candidate_event.start_frame], [min(plot_data) max(plot_data)], 'k'); 
        plot([candidate_event.end_frame candidate_event.end_frame], [min(plot_data) max(plot_data)], 'k'); 

        
    end;

    if (eventData.num_candidate_events >= 2)
        saved2_event = saved1_event + 1;
        candidate_event = eventData.candidate_events(saved2_event);
        
        set(ghandles.saved2_eventnum_text,'String',sprintf('Event: %d of %d', saved2_event, eventData.num_candidate_events));
        set(ghandles.saved2_location_text,'String',sprintf('Location: [%d, %d]', candidate_event.center(1), candidate_event.center(2)));
        set(ghandles.saved2_frame_text,'String',sprintf('Frames: [%d, %d]', candidate_event.start_frame, candidate_event.end_frame));
        set(ghandles.saved2_radius_text,'String',sprintf('Radius: [%d - %d]', min(candidate_event.radii), max(candidate_event.radii)));
        set(ghandles.saved2_source_text,'String',['Source: ', candidate_event.source]);
        set(ghandles.saved2_eventtype_pulldown,'Value', candidate_event.eventtype);      

        set(ghandles.saved_slider,'Value',eventData.num_candidate_events-saved1_event);
        set(ghandles.saved_slider,'Min', 1);
        set(ghandles.saved_slider,'Max', eventData.num_candidate_events-1);
        set(ghandles.saved_slider,'SliderStep',[1/eventData.num_candidate_events 2/eventData.num_candidate_events]);

        axes(ghandles.saved2_axes1);
        cla;
        hold on;
        set(ghandles.saved2_axes1,'fontsize',6);
        title('Center Intensity');       
        xdata = 1:pmaData.numframes;
        plot_data = permute(pmaData.frames(candidate_event.center(1), candidate_event.center(2),1:pmaData.numframes),[3 2 1]);
        if (appState.saved_zoom(2,1) == 1)
            title('Center Intensity (zoomed)');
            xdata = candidate_event.start_frame:candidate_event.end_frame;
            plot_data2 = plot_data(candidate_event.start_frame:candidate_event.end_frame);
            plot_data = plot_data2;
        end;
        plot(xdata, plot_data);
        plot([candidate_event.start_frame candidate_event.start_frame], [min(plot_data) max(plot_data)], 'k'); 
        plot([candidate_event.end_frame candidate_event.end_frame], [min(plot_data) max(plot_data)], 'k'); 

        
        axes(ghandles.saved2_axes2);
        cla;
        hold on;
        set(ghandles.saved2_axes2,'fontsize',6);
        title('Neighborhood Intensity');
        xdata = 1:pmaData.numframes;
        plot_data = candidate_event.neighborhood_intensity;
        if (appState.saved_zoom(2,2) == 1)
            title('Neighborhood Intensity (zoomed)');
            xdata = candidate_event.start_frame:candidate_event.end_frame;
            plot_data2 = plot_data(candidate_event.start_frame:candidate_event.end_frame);
            plot_data = plot_data2;
        end;
        plot(xdata, plot_data);
        plot([candidate_event.start_frame candidate_event.start_frame], [min(plot_data) max(plot_data)], 'k'); 
        plot([candidate_event.end_frame candidate_event.end_frame], [min(plot_data) max(plot_data)], 'k'); 
    
        
        
    end;
    
    
    drawnow;
    return;


%---------------------------------------------------------------------
function display_pma_frame()
    global ghandles;
    global pmaData;
    global appState;
    global eventData;
    axes(ghandles.pma_axes) 

    if (appState.pma_loaded == 0)
        return;
    end;

    if (appState.black_and_white == 0)    
        imshow(pmaData.frames(:,:,appState.display_frame), [0 appState.display_contrast]);
    else
        bwframe = pmaData.frames(:,:,appState.display_frame) > appState.display_threshold;;        
        imshow(bwframe)
    end;

    % Annotate with candidate events
    if (appState.events_analyzed)
        for j=1:eventData.num_candidate_events
            if ((eventData.candidate_events(j).start_frame <= appState.display_frame) & (eventData.candidate_events(j).end_frame) >= appState.display_frame)
               viscircles(flip(eventData.candidate_events(j).center), eventData.candidate_events(j).radii(appState.display_frame-eventData.candidate_events(j).start_frame+1),'LineWidth',1,'LineStyle','--');
            end;
        end;

        for j=1:eventData.num_prospective_events
            if ((eventData.prospective_events(j).start_frame <= appState.display_frame) & (eventData.prospective_events(j).end_frame) >= appState.display_frame)
                viscircles(flip(eventData.prospective_events(j).center), eventData.prospective_events(j).radii(appState.display_frame-eventData.prospective_events(j).start_frame+1),'LineWidth',1,'LineStyle','--','EdgeColor','B');
            end;
        end;
        
        
    end;
    
    drawnow;
    
    
    return;

    

%---------------------------------------------------------------------
function analyze_events_with_progress_bar()

    global appState;
    global pmaData;
    global eventData;
    
    appState.events_analyzed = 0;

    h2 = waitbar(0,'Analyzing frame 1','Name','Analyzing PMA data...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
        'WindowStyle', 'modal');
    setappdata(h2,'canceling',0)      
    
    width = pmaData.width;
    height = pmaData.height;
%    numframes = pmaData.numframes;
    numframes = appState.frame_range;

    num_candidate_events = 0;
    frame_counter = 1;
    for i=appState.filter_min:appState.filter_max
        if getappdata(h2,'canceling')            
            appState.events_analyzed = 0;
            break
        end        
        if (mod(frame_counter,10) == 0)
            waitbar(frame_counter/numframes,h2,sprintf('Analyzing frame %4d of %4d (%d candidate centers)', frame_counter, numframes, num_candidate_events));
        end;
        
        frame_counter = frame_counter + 1;
        
        centers = zeros(width, height);
        radii = zeros(width, height);
    
        frame = pmaData.frames(:,:,i);
        bwframe = frame > 12; % TO-DO: Configurable analysis
        bwframe3 = bwmorph(bwframe,'close');    
        stats = regionprops(bwframe3,'Centroid','MajorAxisLength','MinorAxisLength');

        candidate_centers = round(cat(1, stats.Centroid));
        candidate_MajorAxisLengths = cat(1, stats.MajorAxisLength);
        candidate_MinorAxisLengths = cat(1, stats.MinorAxisLength);
        candidate_radii = mean([candidate_MajorAxisLengths candidate_MinorAxisLengths],2)/2; 

        for j=1:length(candidate_centers)
        
            if (candidate_radii(j) > 3) % TO-DO: Configurable analysis
                centers(candidate_centers(j,2), candidate_centers(j,1)) = 1;
                radii(candidate_centers(j,2), candidate_centers(j,1)) = candidate_radii(j);
            end;
    
        end;
    
        sparse_centers{i} = sparse(centers);
        sparse_radii{i} = sparse(radii);

        num_candidate_events = num_candidate_events + sum(sum(centers));
    end;
    
    % Analysis needs to look to the future even if filter_max < numframes
    if (appState.filter_max < pmaData.numframes)
        sparse_matrix(pmaData.width,pmaData.height) = 1;
        for i=appState.filter_max+1:pmaData.numframes;
            sparse_centers{i} = sparse(sparse_matrix);
            sparse_radii{i} = sparse(sparse_matrix);
        end;
    end;
    
    num_candidate_events = 0;

    frame_counter = 1;
    for i2=appState.filter_min:appState.filter_max

        if getappdata(h2,'canceling')            
            appState.events_analyzed = 0;
            break
        end        
        if (mod(frame_counter,10) == 0)
            waitbar(frame_counter/numframes,h2,sprintf('Second pass: frame %4d of %4d (%d candidate events)', frame_counter, numframes, num_candidate_events));
        end;
        
        frame_counter = frame_counter + 1;
        
        [rows,cols,radii] = find(sparse_radii{i2});
    
        for iterator=1:length(radii) 
        
            i = rows(iterator);
            j = cols(iterator);
            radius = radii(iterator);

            tradius = floor(1.5*radius); %TO-DO: Configurable analysis
            tradius = saturate_radius(tradius, j, i, width, height);

            hit_count = 0;
            miss_count = 0;
            clear candidate_radii;
            % Current approach: for each center/radius, see how many
            % consecutive overlapping (in X/Y) centers you can find. 
            k = 1;
            while (miss_count < 3) % TO-DO: Configurable analysis
                candidate_radii(k) = tradius;
                
                if (i2+k>pmaData.numframes)
                    break;
                end;

                % Did we find another center within the tradius? Update the
                % count (and the tradius. Currently we take the most recent
                % radius but could argue for finding the max). 
                if (sum(sum(full(sparse_centers{i2+k}(i-tradius:i+tradius,j-tradius:j+tradius)))) ~= 0)
                    hit_count = hit_count + 1;
                    miss_count = 0;
                    tradius2 = max(max(full(sparse_radii{i2+k}(i-tradius:i+tradius,j-tradius:j+tradius))));
                    tradius2 = floor(1.5*tradius2); %TO-DO: Configurable analysis
                    tradius = saturate_radius(tradius2, j, i, width, height);                   
                else
                    miss_count = miss_count + 1;
                end;
                
                k = k + 1;
            end;
            
            k = k - 1;
            
            % If our hit_count is "better" than our miss_count (+ other
            % conditions), consider it a candidate
            if (hit_count >= 8) % TO-DO: Configurable analysis

                candidate_event.start_frame = i2;
                candidate_event.end_frame = i2+k-1;
                candidate_event.center = [i, j];
                candidate_event.radii = candidate_radii;
                candidate_event.neighborhood_intensity = 0;
                candidate_event.source = 'Auto';
                candidate_event.eventtype = 1;
                
                

                % Get rid of overlapping candidates              
                valid_candidate = 1;

                for k2=1:num_candidate_events

                    % This is a bit conservative. If the radii ever overlap
                    % we consider it an overlapping event.
                    existing_max_radius = max(candidate_events(k2).radii);
                    candidate_max_radius = max(candidate_event.radii);
                    tradius = max(existing_max_radius, candidate_max_radius);
                                      
                    if ((candidate_events(k2).start_frame <= candidate_event.start_frame) & (candidate_events(k2).end_frame >= candidate_event.start_frame) ...
                            & (abs(candidate_events(k2).center(1)-candidate_event.center(1)) < tradius) & (abs(candidate_events(k2).center(2)-candidate_event.center(2)) < tradius))
                        valid_candidate = 0;
                        break;
                    end;
                end;

                if (valid_candidate == 1)
                    %disp(['possible canddiate at [', num2str(i), ', ', num2str(j), '], frames: [', num2str(candidate_event.start_frame), ', ', num2str(candidate_event.end_frame), '], radii: ', num2str(candidate_event.radii), ']']);               

                    % Calculate neighborhood intensity
                    for k3=1:pmaData.numframes
                        neighborhood_radius = 4;
                        neighborhood_radius = saturate_radius(neighborhood_radius, j, i, width, height);  
                        candidate_event.neighborhood_intensity(k3) = mean(mean(pmaData.frames(i-neighborhood_radius:i+neighborhood_radius,j-neighborhood_radius:j+neighborhood_radius,k3)));
                    end;
                        
                    %for k3=1:length(candidate_event.radii)
                    %    candidate_event.neighborhood_intensity(k3) = mean(mean(pmaData.frames(i-candidate_event.radii(k3):i+candidate_event.radii(k3),j-candidate_event.radii(k3):j+candidate_event.radii(k3),candidate_event.start_frame+k3-1)));
                    %end;
                    
                    num_candidate_events = num_candidate_events + 1;
                    candidate_events(num_candidate_events) = candidate_event;
                end;
                
            end;    
        end;
    end;    
    
    
    delete(h2);    
    
    if (num_candidate_events > 0)
        appState.events_analyzed = 1;
        eventData.candidate_events = candidate_events;
        eventData.num_candidate_events = num_candidate_events;
    end;

    eventData.prospective_events = 0;
    eventData.num_prospective_events = 0;
    
%---------------------------------------------------------------------
% Helper function to make sure we never go out of bound when looking
% at pixel neighborhoods. This is actually a little tricky
function tradius2 = saturate_radius(tradius, center_i, center_j, width, height)

    % If the center_i <= the radius, we have to shrink the radius
    if (center_i <= tradius)
       tradius = center_i-1;
    end;
    
    % If center_j is still <= the radius, shrink it. We are guaranteed to
    % already by safe on the i dimension.
    if (center_j <= tradius)
       tradius = center_j-1;
    end;
        
    % Same approach for overlapping the end of a row/column
    if (center_i+tradius >= width)
       tradius = width-center_i;
    end;
    if (center_j+tradius >= height)
       tradius = height-center_j;
    end;

    tradius2 = tradius;

    
%---------------------------------------------------------------------
function load_pmafile_with_progress_bar(filename)

    global appState;
    global fileData;
    global pmaData;
    appState.pma_loaded = 0;

    waitbar_text = strrep(fileData.fileName,'_','\_');
    
    h = waitbar(0,['Reading ', waitbar_text],'Name','Loading PMA frames...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
        'WindowStyle', 'modal');
    setappdata(h,'canceling',0)      

    
    
    fileinfo    = dir(filename);
    filesize    = fileinfo.bytes;
    fid         = fopen(filename);  
    pmaData.width   = fread(fid, 1, 'int16');
    pmaData.height  = fread(fid, 1, 'int16'); 
    pmaData.numframes     = floor((filesize-4)/(pmaData.width*pmaData.height));
    pmaData.frames  = zeros(pmaData.height,pmaData.width,pmaData.numframes);
    
    appState.pma_loaded = 1;

    pma_len = pmaData.numframes;
    for i=1:pma_len
        if getappdata(h,'canceling')            
            appState.pma_loaded = 0;
            break
        end
        
        if (mod(i,10) == 0)
            waitbar(i/pma_len,h,sprintf('Reading frame %4d of %4d',i, pma_len));
        end;
            
        pmaData.frames(:,:,i) = fread(fid, [pmaData.width, pmaData.height], 'uint8');
        pmaData.frames(:,:,i) = rot90(pmaData.frames(:,:,i));
    end
    delete(h);
    fclose(fid);

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(hObject,'Position'); %Pos(3 and 4) are width and height




% --- Executes on slider movement.
function current_frame_slider_Callback(hObject, eventdata, handles)
% hObject    handle to current_frame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    global appState;
    global ghandles;

    appState.display_frame = round(get(hObject,'Value'));
    set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
    display_pma_frame();

    
    
% --- Executes during object creation, after setting all properties.
function current_frame_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_frame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in black_and_white.
function black_and_white_Callback(hObject, eventdata, handles)
% hObject    handle to black_and_white (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of black_and_white
    global appState;
    global ghandles;
    appState.black_and_white = get(hObject,'Value');
    if (appState.black_and_white == 0)
        set(ghandles.contrast_threshold_edit,'String',appState.display_contrast);
        set(ghandles.contrast_threshold_text,'String','Contrast [1-255]: ');
    else
        set(ghandles.contrast_threshold_edit,'String',appState.display_threshold);
        set(ghandles.contrast_threshold_text,'String','Threshold [1-255]: ');
    
    end;
    display_pma_frame();
    

function contrast_threshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contrast_threshold_edit as text
%        str2double(get(hObject,'String')) returns contents of contrast_threshold_edit as a double
    global appState;
    global ghandles;

    % Check to see if in range
    value = str2num(get(hObject,'String'));
    if (value < 1)
        value = 1;
    end;
    if (value > 255)
        value = 255;
    end;    
    set(hObject,'String',num2str(value));
    
    if (appState.black_and_white == 0)        
        appState.display_contrast = value;
    else
        appState.display_threshold = value;
    end;
    
    display_pma_frame();



% --- Executes during object creation, after setting all properties.
function contrast_threshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function current_frame_edit_Callback(hObject, eventdata, handles)
% hObject    handle to current_frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global appState;
    global ghandles;
    global pmaData;

    display_frame = str2num(get(hObject,'String'));
    if ((display_frame < 1) | (display_frame > pmaData.numframes) | (mod(display_frame,1) ~= 0))
        set(hObject,'String', num2str(appState.display_frame));
    else
        appState.display_frame = display_frame;
        set(ghandles.current_frame_slider,'Value',appState.display_frame);
        display_pma_frame();
    end;
        



% --- Executes during object creation, after setting all properties.
function current_frame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in play_pause.
function play_pause_Callback(hObject, eventdata, handles)
% hObject    handle to play_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global appState;
global pmaData;
global ghandles;

if (appState.pma_loaded == 0)
    return;
end;

% Pause button has been hit
if (appState.play_pause == 0)
    set(hObject,'String', char(hex2dec('25BA'))); % Play button
    set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
    set(ghandles.current_frame_slider,'Value',appState.display_frame);    
    appState.play_pause = 1;
    stop(ghandles.timer);

% Play (or rewind) button has been hit
elseif (appState.play_pause == 1)

    if (appState.display_frame >= pmaData.numframes) % Likely rewind condition
        appState.display_frame = 1;
        set(hObject,'String', char(hex2dec('25BA'))); % Play button
        set(ghandles.current_frame_edit,'String',sprintf('%d', appState.display_frame));
        set(ghandles.current_frame_slider,'Value',appState.display_frame);            
    else
        set(hObject,'String', [char(10074), char(10074)]); % Pause button
        appState.play_pause = 0;
        start(ghandles.timer);
    end;
end;


function close_request(varargin)
    global ghandles;

try
	listOfTimers = timerfindall; % List all timers, just for info.
	% Get handle to the one timer that we should have.
	if isempty(listOfTimers)
		% Exit if there is no timer to turn off.
		return;
	end
	handleToTimer = ghandles.timer;    
	% Stop that timer.
	stop(handleToTimer);
	% Delete all timers from memory.
	listOfTimers = timerfindall;
	if ~isempty(listOfTimers)
		delete(listOfTimers(:));
	end
catch ME
	errorMessage = sprintf('Error in StopTimer().\nThe error reported by MATLAB is:\n\n%s', ME.message);
	fprintf('%s\n', errorMessage);
	WarnUser(errorMessage);
end

    delete(ghandles.figure1);

return;


% --- Executes on selection change in prospective1_eventtype_pulldown.
function prospective1_eventtype_pulldown_Callback(hObject, eventdata, handles)
% hObject    handle to prospective1_eventtype_pulldown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global appState;
    global eventData;
    eventData.prospective_events(appState.current_prospective_event).eventtype = get(hObject,'Value');
    display_prospective_frame();



% --- Executes during object creation, after setting all properties.
function prospective1_eventtype_pulldown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prospective1_eventtype_pulldown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in saved1_eventtype_pulldown.
function saved1_eventtype_pulldown_Callback(hObject, eventdata, handles)
% hObject    handle to saved1_eventtype_pulldown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global appState;
    global eventData;
    eventData.candidate_events(appState.current_saved_event).eventtype = get(hObject,'Value');
    display_saved_frame();

% --- Executes during object creation, after setting all properties.
function saved1_eventtype_pulldown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saved1_eventtype_pulldown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saved1_removeevent_pushbutton.
function saved1_removeevent_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saved1_removeevent_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global appState;
    global eventData;

    saved1_event = appState.current_saved_event;
    eventData.candidate_events(saved1_event) = [];
    eventData.num_candidate_events = eventData.num_candidate_events - 1;
    if (saved1_event >= eventData.num_candidate_events)
        appState.current_saved_event = eventData.num_candidate_events - 1;
    end;

    display_saved_frame();

% --- Executes on button press in saved2_removeevent_pushbutton.
function saved2_removeevent_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saved1_removeevent_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global appState;
    global eventData;

    saved2_event = appState.current_saved_event+1;
    eventData.candidate_events(saved2_event) = [];
    eventData.num_candidate_events = eventData.num_candidate_events - 1;
    if (saved2_event >= eventData.num_candidate_events)
        appState.current_saved_event = eventData.num_candidate_events - 1;
    end;

    display_saved_frame();

    
    
% --- Executes on selection change in saved2_eventtype_pulldown.
function saved2_eventtype_pulldown_Callback(hObject, eventdata, handles)
% hObject    handle to saved2_eventtype_pulldown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global appState;
    global eventData;
    eventData.candidate_events(appState.current_saved_event+1).eventtype = get(hObject,'Value');
    display_saved_frame();


% --- Executes during object creation, after setting all properties.
function saved2_eventtype_pulldown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saved2_eventtype_pulldown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function saved_slider_Callback(hObject, eventdata, handles)
% hObject    handle to saved_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global appState;
    global ghandles;
    global eventData;

    appState.current_saved_event = eventData.num_candidate_events - round(get(hObject,'Value'));
    if (appState.current_saved_event == eventData.num_candidate_events)
        appState.current_saved_event = appState.current_saved_event - 1;
    end;
    display_saved_frame();
    

% --- Executes during object creation, after setting all properties.
function saved_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saved_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in prospective1_addevent_pushbutton.
function prospective1_addevent_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to prospective1_addevent_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global appState;
    global eventData;

    prospective1_event = appState.current_prospective_event;
    current_event = eventData.prospective_events(prospective1_event);
    eventData.prospective_events(prospective1_event) = [];
    eventData.num_prospective_events = eventData.num_prospective_events - 1;
    if (prospective1_event > eventData.num_prospective_events)
        appState.current_prospective_event = eventData.num_prospective_events; 
    end;

    % Add it to the end. I think this is ok.
    eventData.candidate_events = [eventData.candidate_events current_event];
    eventData.num_candidate_events = eventData.num_candidate_events + 1;

    display_pma_frame();
    display_prospective_frame();
    display_saved_frame();


% --- Executes on button press in prospective1_ignoreevent_pushbutton.
function prospective1_ignoreevent_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to prospective1_ignoreevent_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global appState;
    global eventData;

    prospective1_event = appState.current_prospective_event;
    eventData.prospective_events(prospective1_event) = [];
    eventData.num_prospective_events = eventData.num_prospective_events - 1;
    if (prospective1_event > eventData.num_prospective_events)
        appState.current_prospective_event = eventData.num_prospective_events ;
    end;

    display_pma_frame();
    display_prospective_frame();



% --- Executes on slider movement.
function prospective_slider_Callback(hObject, eventdata, handles)
% hObject    handle to prospective_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global appState;
    global ghandles;
    global eventData;

    round(get(hObject,'Value'));
    appState.current_prospective_event = eventData.num_prospective_events - round(get(hObject,'Value'))+1;
%    if (appState.current_prospective_event == eventData.num_prospective_events)
%        appState.current_prospective_event = appState.current_prospective_event - 1;
%    end;
    display_prospective_frame();


% --- Executes during object creation, after setting all properties.
function prospective_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prospective_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function filter_edit1_Callback(hObject, eventdata, handles)
% hObject    handle to filter_edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global appState;
    global ghandles;

    filter_min = str2num(get(hObject,'String'));
    if ((filter_min >= appState.filter_max) | (filter_min < 1) | (mod(filter_min,1) ~= 0))
        set(hObject,'String',appState.filter_min);
        return;
    else
        appState.filter_min = filter_min;
        display_pma_frame();
        display_prospective_frame();
        display_saved_frame();
    end;



% --- Executes during object creation, after setting all properties.
function filter_edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_edit2_Callback(hObject, eventdata, handles)
% hObject    handle to filter_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global appState;
    global ghandles;
    
    filter_max = str2num(get(hObject,'String'));
    if ((filter_max <= appState.filter_min) | (filter_max > 9999) | (mod(filter_max,1) ~= 0))
        set(hObject,'String',appState.filter_max);
        return;
    else
        appState.filter_max = filter_max;
        display_pma_frame();
        display_prospective_frame();
        display_saved_frame();
    end;


% --- Executes during object creation, after setting all properties.
function filter_edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
