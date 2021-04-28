function varargout = p4(varargin) %this is the setup code for the figure and initializes the callbacks

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @p4_OpeningFcn, ...
                   'gui_OutputFcn',  @p4_OutputFcn, ...
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


% Static values
global Fmin % Minimun frequency
Fmin = 50;
global Fmax % Maximum frequency
Fmax = 600;
global threshold_amp  % Minimum amplitude for a peak in order for autocorrelation to be called - the higher number the less peaks it will accept
threshold_amp = 0.2;
global harmonic_deviation % Max deviation from fundemental frequency in order to be harmonic (in %)
harmonic_deviation = 0.2;
global threshold_ratio % Max ratio in heights of autocorr allowed between fundamental and harmonics
threshold_ratio = 0.4;
global f0_target_sample_time % sample time from the previous f0 to be used to decide the next f0
f0_target_sample_time = 0.08;
global target_tol  % Max change from previous f0 to target f0 (in%)
target_tol = 0.3;
global initFs % Sampling frequency (Hz)
initFs = 44000;
global nbits % Bits per sample (can look up audio bit depth for details)(this can be reduced to 8 to save computational power
nbits = 16;
global precision % Save as single precision vector
precision = 'single';
global nchans  % Single audio channel (mono)
nchans = 1;
global initscale_factor %factor of how much to scale when compressing the signal (1 keeps it - 0 makes it completely pitch accurate) 
initscale_factor = 0.15;

% --- Executes just before visuals - setup all of the handles for the
% software and the figure
function p4_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to project (see VARARGIN)

% Choose default command line output for project
handles.output = hObject;
set(handles.axes1,'XTickLabel',{});
set(handles.axes1,'YTickLabel',{});

% Axes2 is just for mouse click input
set(handles.axes2,'XTick',[]);
set(handles.axes2,'YTick',[]);
set(handles.axes2,'XTickLabel',{});
set(handles.axes2,'YTickLabel',{});
set(handles.axes2,'ButtonDownFcn','p4(''mouseclick'',gcbo,[],guidata(gcbo))');


% Set defaults for the program - retrived from set_default function
handles = set_defaults(handles);

% Update handles structure
guidata(hObject, handles);

%setup of the default handles for the software
function handles = set_defaults(handles)

% Scale setup 
handles.scale_options.indices = [];    
handles.scale_options.freqs = [];
handles.scale_options.notes = {};
handles.scale_options.fund_index = [];

% recording object and playback object
handles.record_obj = [];
handles.play_obj = [];

% Initialize structure for saving/manipulating sound data
handles.sound.A = [];
handles.sound.A_corrected = [];
handles.sound.Fs = [];
handles.sound.t = [];
handles.sound.f0 = [];
handles.sound.f0_save = [];
handles.sound.f0_corrected = [];
handles.sound.f0_corrected_save = [];
handles.sound.tcalc = [];
handles.sound.selected_points = logical([]);
handles.sound.selected_points_save = logical([]);

% Initialize status reporting
handles.status.isrecording = false;
handles.status.isplaying = false;
handles.status.iscompress = false;
handles.status.issnap = false;
handles.status.isview = false;
handles.status.ismousezoom = false;
handles.status.Xlim = [];
handles.status.Ylim = [];


[handles.scale_options.indices,...
    handles.scale_options.freqs,...
    handles.scale_options.notes] = ...
    get_scale('Cmajor');

handles = update_GUI(handles);
set(handles.original_button,'value',1);
set(handles.plot_selection,'value',2);

%setup in order to update the plot - did so it's always on pitch to make
%code shorter and we don't need waveform
function handles = update_plot(handles)

val = get(handles.plot_selection,'value');
axes(handles.axes1);

if val == 2
    if ~isempty(handles.sound.tcalc) && ~isempty(handles.sound.f0)
        hplot = [];
        f0 = handles.sound.f0;
        f0(f0 < 1) = NaN;
        if all(isnan(f0))
            cla
            hwarn = warndlg('No pitch detected with current settings!','WARNING');
            waitfor(hwarn);
        else
            hplot(1) = semilogy(handles.axes1,handles.sound.tcalc,f0,'.-','color',[0.7 0.7 0.7]);
            if ~isempty(handles.sound.f0_corrected)
                hold on
                f02 = handles.sound.f0_corrected;
                f02(f02 < 1) = NaN;
                hplot(2) = semilogy(handles.axes1,handles.sound.tcalc,f02,'.-','color','b');
                
                if any(handles.sound.selected_points)
                    f03 = f02;
                    f03(~handles.sound.selected_points) = NaN;
                    semilogy(handles.axes1,...
                        handles.sound.tcalc,f03,'.-','color','r');
                end
                hold off
            end
            grid on
            Ylim = [min(f0)/1.2 max(f0)*1.2];
            set(handles.axes1,'Ylim',Ylim);
            set(handles.axes1,'Ytick',handles.scale_options.freqs);
            set(handles.axes1,'YtickLabel',handles.scale_options.notes);
            if isempty(handles.status.Xlim)
                set(handles.axes1,'Xlim',[min(handles.sound.t) max(handles.sound.t)]);
            else
                set(handles.axes1,'Xlim',handles.status.Xlim);
            end
        end
    else
        cla
    end
end

% Invisible axes for getting mouse information
set(handles.axes2,'Xlim',get(handles.axes1,'Xlim'));
set(handles.axes2,'Ylim',get(handles.axes1,'Ylim'));
set(handles.axes2,'Yscale',get(handles.axes1,'Yscale'));
set(handles.axes2,'ButtonDownFcn','p4(''mouseclick'',gcbo,[],guidata(gcbo))');
set(handles.axes2,'Color','none');
axes(handles.axes2);

clear f0 f02 f03

%function for the mode the program is in and how the buttons should be
%available or not
function handles = update_GUI(handles)

if handles.status.isrecording
    set(handles.play_button,'enable','off');
    set(handles.record_button,'enable','off');
    set(handles.stop_button,'enable','on');
elseif handles.status.isplaying
    set(handles.play_button,'enable','off');
    set(handles.record_button,'enable','off');
    set(handles.stop_button,'enable','off');
else
    if ~isempty(handles.sound.A)
        set(handles.play_button,'enable','on');
    else
        set(handles.play_button,'enable','off');
    end
    set(handles.record_button,'enable','on');
    set(handles.stop_button,'enable','off');
end
    
% --- Executes on button press in record_button.
function record_button_Callback(hObject, ~, handles)
% hObject    handle to record_button (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
global initFs
global nbits
global nchans
try
    handles.record_obj = audiorecorder(initFs,...
        nbits,nchans);
    record(handles.record_obj);

    handles.status.isrecording = true;

    handles = update_GUI(handles);
    guidata(hObject, handles);

end

% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, ~, handles)
% hObject    handle to stop_button (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
global initFs
global precision

if handles.status.isrecording
    stop(handles.record_obj);
    handles.status.isrecording = false;
    
    % Get the waveform
    handles.sound.A = getaudiodata(handles.record_obj,...
        precision);
    
    % Save sampling frequency
    handles.sound.Fs = initFs;
    
    
    % Calculate time vector
    handles.sound.t = (0:length(handles.sound.A)-1)./initFs;
    
    
    
    handles = run_pitch_detection(handles);
    handles.status.isview = true;
    handles = update_plot(handles);
    handles = update_GUI(handles);
    guidata(hObject, handles);
    
end

% --- Executes on button press in play_button.
function play_button_Callback(hObject, ~, handles)
% hObject    handle to play_button (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Decide which to play
val = get(handles.original_button,'val');
if val == 1
    if isempty(handles.sound.A)
        return
    end
    A = handles.sound.A;
else
    if isempty(handles.sound.A_corrected)
        handles = run_pitch_correction(handles);
    end
    A = handles.sound.A_corrected;
end



handles.status.isplaying = true;
handles = update_GUI(handles);
handles = update_plot(handles);

% Get the current plot limits
Ylim = get(handles.axes1,'Ylim');
if isempty(handles.status.Xlim)
    tmin = min(handles.sound.t);
    tmax = max(handles.sound.t);
else
    [what,Imin] = min(abs(handles.status.Xlim(1)-handles.sound.t));
    [what,Imax] = min(abs(handles.status.Xlim(2)-handles.sound.t));
    tmin = handles.sound.t(Imin);
    tmax = handles.sound.t(Imax);
    A = A(Imin:Imax);
end

% Draw the vertical line to scan across
axes(handles.axes1);
hline = line(tmin.*[1 1],[Ylim],'color','r');

pause(0.1);

% Start playing the wave file
p = audioplayer(A, handles.sound.Fs); 
play(p);
tic;
while 1
    tnow = toc + tmin;
    if tnow > tmax
        break
    end
    % Scan the line
    set(hline,'XData',[tnow tnow]);
    pause(0.02);
end
set(hline,'XData',[tmin tmin]);

% Done
handles.status.isplaying = false;
handles = update_GUI(handles);

guidata(hObject, handles);

clear A;
axes(handles.axes2);


% --- Executes on selection change in plot_selection. //waveform plot
function plot_selection_Callback(hObject, ~, handles)
% hObject    handle to plot_selection (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A_corrected)
    handles = run_pitch_correction(handles);
end
handles = update_plot(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
% sets the colors of the background
function plot_selection_CreateFcn(hObject, ~, ~)
% hObject    handle to plot_selection (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in original_button.
function original_button_Callback(hObject, ~, handles)
% handles    structure with handles and user data (see GUIDATA)

set(handles.original_button,'val',1);
set(handles.modified_button,'val',0);


% --- Executes on button press in modified_button.
function modified_button_Callback(hObject, eventdata, handles)
% hObject    handle to modified_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of modified_button

set(handles.original_button,'val',0);
set(handles.modified_button,'val',1);




% --- Executes on button press in nearest pitch button.
function nearest_pitch_button_Callback(hObject, ~, handles)

if isempty(handles.sound.f0_corrected) || ~any(handles.sound.selected_points) ||...
        isempty(handles.scale_options.freqs);
    return
end

pitch_direction = get(gcbo,'String');

% Pull out the frequency and time of the selected points
I = find(handles.sound.selected_points);
sp = handles.sound.selected_points(I(1):I(end));
f0 = handles.sound.f0_corrected(I(1):I(end));
t = handles.sound.tcalc(I(1):I(end));
t = t-t(1);

% Come up with a envelope for compressing that goes from 1 down to the
% desired compression factor over the desired time interval
if length(t) < 2
    where = 1;
elseif 0 > t(2)
    [what,where] = min(abs(t - 0));
else
    where = 1;
end

tmp = f0(where:end);
mean_f0 = 10^mean(log10(tmp(sp(where:end))));
if strcmp(pitch_direction,'Down')
    Iu = find(handles.scale_options.freqs < 0.99*mean_f0);
    mean_new = max(handles.scale_options.freqs(Iu));
elseif strcmp(pitch_direction,'Up')
    Iu = find(handles.scale_options.freqs > 1.01*mean_f0);
    mean_new = min(handles.scale_options.freqs(Iu));
end
factor = mean_new/mean_f0.*ones(size(t));
if where > 1
    factor(1:where) = linspace(1,factor(end),where);
end

% Scale the deviations around the mean in a log sense
f0 = f0.*factor;
handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
handles.sound.f0_save(end+1,:) = handles.sound.f0;
handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;
handles.sound.f0_corrected(handles.sound.selected_points) = f0(sp);

handles.sound.A_corrected = [];

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

clear t f0 what where factor


% --- Executes on button press in apply_compress_button.
function cleaner_button_Callback(hObject, ~, handles)
global initscale_factor
                  
if isempty(handles.sound.f0_corrected) || ~any(handles.sound.A)
    return
end

% get frequency and correlating time of the points selected
f0 = handles.sound.f0_corrected(handles.sound.selected_points);
I = find(handles.sound.selected_points);
t = handles.sound.tcalc(I(1):I(end));
t = t-t(1);

% compress the signal down to the desired factor which has been selected in
% initscale_factor
factor = ones(size(t)).*initscale_factor;
if 0 > t(2)
    [what,where] = min(abs(t - 0));
    factor(1:where) = linspace(1,initscale_factor,where);
end
factor = factor(handles.sound.selected_points(I(1):I(end)));

% find the deviations and put them in logarithmic form and take the mean of
% this in order to a
f0 = 10.^(mean(log10(f0)) + factor.*(log10(f0)-mean(log10(f0))));
handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
handles.sound.f0_save(end+1,:) = handles.sound.f0;
handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;
handles.sound.f0_corrected(handles.sound.selected_points) = f0;

handles.sound.A_corrected = [];

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

clear t f0 what where factor



function handles = run_pitch_detection(handles)
global Fmin 
global Fmax 
global threshold_amp  
global harmonic_deviation 
global threshold_ratio 
global target_tol  

if isempty(handles.sound.A) %ensures if there's no signal it doesn't run the pitch detection but stops.
    return
end

plot_on = false;

block = 2*round(0.04*handles.sound.Fs/2);  % takes a block of the signal to find the pitch (this can be changed according to what makes sense)
step = block/2; % how many steps in each block, this is how many values from peak to peak

N = floor((length(handles.sound.A)-block)/step); % Number of frequency computations to be done according to the size of signal and block/step amount

% Initialize variables for storing results
f0 = zeros(1,N);                        % making a 1D vector which is the size of the maximum amount of results, fills it up with 0 for now
tcalc = zeros(1,N);                     % same as before, making a 1D vector that stores the time value of calculated results size of max computations
f0_target = [];                         % to know what block of signal will be the next one computed

f0_samples = floor...                   % figures out how many f0 samples there will be, rounded down towards nearest whole number.
    (0.08/(step/handles.sound.Fs));

% makes waitbar for visual purpouses 
hwait = waitbar(0,'Determining pitch...');
set(hwait,'Name','Please Wait');
waitbar_count = 0;
waitbar_update = round(N/10);
pause(0.001);


I = 1;  %initialises I which is used to know how many blocks have been through, to calculate the following block etc. 
for n = 1:N         %starts a for loop which runs from 0 and until it has gone through every computation as declared with N before
    
    % Update waitbar
    if waitbar_count > waitbar_update
        waitbar(n/N,hwait);
        waitbar_count = 0;
    end
    waitbar_count = waitbar_count + 1;
    
    % extracts the block of signal that will be processed in this
    % itteration
    Atemp = handles.sound.A(I:I+block-1); %makes matrix with the signal values that goes from I -> the end of this block being processed -1
    ttemp = handles.sound.t(I:I+block-1); %same here just with time instead of the signal
    I = I+step;
    
    % call the find_f0_timedomain function to find the frequency -
    % important that n isn't confused with N.
    [f0(n),acorr,Nshifts,Tshifts] = ...
        find_f0_timedomain(Atemp,handles.sound.Fs,Fmin,Fmax,...
        threshold_amp,threshold_ratio,harmonic_deviation,...
        f0_target,target_tol);
    tcalc(n) = median(ttemp);
    
    if plot_on
        figure(2)
        subplot(2,1,1),
        plot(ttemp,Atemp)
        title(sprintf('Frequency: %0.1f Hz',f0(n)))
        subplot(2,1,2),
        plot(Tshifts,acorr)
        pause
    end
    
    % If the previous and current are invalid, then make the last one invalid
    % as well
    if n > 2
        if f0(n-2) < 1 && f0(n) < 1
            f0(n-1) = 0;
        end
    end
            
    % check previous samples to decide target frequency for next
    % computation in order for the signal to be continous.
    if n >= f0_samples
        f0_chunk = f0(n-f0_samples+1:n);
        f0_target = f0_chunk(f0_chunk>0);
        if ~isempty(f0_target)
            f0_target = mean(f0_target);
        end
    end
    
    
end

if ishandle(hwait)
    close(hwait);
end

% Save calculation
handles.sound.A_corrected = [];
handles.sound.f0 = f0;
handles.sound.f0_save = [];
handles.sound.f0_corrected = f0;
handles.sound.f0_corrected_save= [];
handles.sound.tcalc = tcalc;
handles.sound.selected_points = false(size(f0));
handles.sound.selected_points_save= false(size(f0));
handles.sound.block = block;                  
handles.sound.step = step;
    
function handles = run_pitch_correction(handles)
global Fmin 
global initscale_factor

if isempty(handles.sound.A)
    return
end

plot_on = false;

A = handles.sound.A;
t = handles.sound.t;
f0 = handles.sound.f0;
f02 = handles.sound.f0_corrected;
Fs = handles.sound.Fs;
scale_factor = zeros(size(f0));
A_corrected = zeros(size(A));
block = handles.sound.block;      %block is 0,04            
step = handles.sound.step;
N = length(f0);

max_acorr_shift = 0;
max_acorr_amp = 0;  

% Waitbar for visual purpouses
hwait = waitbar(0,'Correcting pitch...');
set(hwait,'Name','Please Wait');
waitbar_count = 0;
waitbar_update = round(N/10);
pause(0.001);

I = 1;
for n = 1:N
    
    % Update waitbar
    if waitbar_count > waitbar_update
        waitbar(n/N,hwait);
        waitbar_count = 0;
    end
    waitbar_count = waitbar_count + 1;

    % Pull out a chunk of the signal
    Atemp = A(I:I+block-1); %same as in pitch detection
    ttemp = t(I:I+block-1);
       
    if f0(n) > 0 && ~isnan(f0(n)) %if statement to ensure the frequency isn't empty
        
        % find the scale_factor for the given itteration of n.
        scale_factor(n) = f02(n)/f0(n);
    
        % does the interpolation to get new datapoints of the signal.
        tp = mean(ttemp) + (ttemp-mean(ttemp)).*scale_factor(n); % subtracts the mean of ttemp and adds the scale_factor whereafter it adds the mean of ttemp to get new data points of scaled values
        Ainterp = interp1(ttemp,Atemp ,tp)'; %does liniar interpolation (interp1(sample point,corresponding value, query point)
        Ivalid = find(~isnan(Ainterp)); %finds all nonvalid values in Ainterp
        Ainterp(isnan(Ainterp)) = 0; %puts all nonvalid values to 0
        
        Nperiod = ceil(1/(f02(n)/Fs)); %makes the N for this period, divides the corrected frequency with the sampling rate and then divides with 1 and rounds it to nearest positive number.
    else
        % if there is no change in shift run this
        scale_factor(n) = 1;
        Ainterp = Atemp;  
        Ivalid = find(~isnan(Ainterp));
        Nperiod = ceil(1/(Fmin/Fs));
    end
    
    if n == 1       %first itteration of the for loop assign Ainterp the A_corrected values in the form of a matrix indexed from I->I+block-1
        A_corrected(I:I+block-1) = Ainterp;
        
    else
        
        % Gets a period from the new waveform of the interpolated values
        Achunk = Ainterp(Ivalid(1):Ivalid(1)+Nperiod-1); %makes a matrix which goes from first invalid->invalid+Nperiod-1
        factor = sum(abs(Achunk))*2; %set's all values to abselute values and takes the sum *2
        
        if all(scale_factor(n-1:n) == 1) %if all numbers are ==1 then skip the next part.

        else

            %setup variables for autocorrelation
            max_acorr_amp = 0; 
            max_acorr_shift = 0;

            for Nshift = -round(Nperiod/2)+ (1:Nperiod) %make Nshift an matrix based on the Nperiod values. 

                % Calculate makeshift autocorrelation
                acorr = 1-sum(abs(Achunk - A_corrected((I:I+Nperiod-1) + Nshift + Ivalid(1))))./factor;

                if acorr > max_acorr_amp %if the autocorrelation is over 0
                    max_acorr_amp = acorr; %set the max ammplitude to correlated value
                    max_acorr_shift = Nshift;
                end
            end
        end
        
        [what,where] = min(abs(Achunk - ...             %finds the minimum values in correlated values
            A_corrected((I:I+Nperiod-1) + max_acorr_shift + Ivalid(1))));
        
        if plot_on %plot stuff
            Asave = A_corrected((I:I+Nperiod-1) + max_acorr_shift + Ivalid(1));
        end
        %takes the correlated values from where->length of the array and
        %replace the empty values in Ainterp with the corrected values from
        %where->end to make the corrected signal complete. 
        A_corrected((I+where-1:I+length(Ivalid)-1) + max_acorr_shift + Ivalid(1)) = ... 
            Ainterp(Ivalid(where:end));
        
        if plot_on %plot stuff
            tt = 1:length(Achunk);
            figure(1);
            plot(tt,Achunk,tt,Asave,tt,A_corrected((I:I+Nperiod-1) + max_acorr_shift + Ivalid(1)))
            title(sprintf('Acorr = %0.3f (Nshift = %i)',max_acorr_amp,max_acorr_shift))
            hold on
            plot(tt(where),Achunk(where),'*')
            hold off
            pause
        end

        

    end
    
    I = I+step; %counts I with the step for the next itteration
    
end

%reset values 
handles.sound.A_corrected = A_corrected;
initscale_factor = scale_factor;

clear f0 f02 scale factor A A_corrected

if ishandle(hwait)
    close(hwait);
end



% Function for making the mouse work on the plot
function mouseclick(hObject, ~, handles)


if isempty(handles.sound.f0)
    return
end

if handles.status.ismousezoom
    
  
    
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                  % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    xminmax = [p1(1) p1(1)+offset(1)];
    yminmax = [p1(2) p1(2)+offset(2)];
    
    if diff(xminmax) < diff(handles.sound.t(1:2))
        dt = diff(handles.sound.t(1:2));
        xminmax = mean(xminmax) + [-dt/2 dt/2];
    end

    handles.status.Xlim = xminmax;
    handles.status.ismousezoom = false;
    handles.status.isview = true;
    
    handles = update_plot(handles);
    handles = update_GUI(handles);
    
    guidata(hObject, handles);
    set(gcf,'Pointer','arrow');
    
    return;

elseif get(handles.plot_selection,'value') ~= 1
    clicktype = get(gcf,'SelectionType');
    if strcmp(clicktype,'normal') || strcmp(clicktype,'alt')

        point1 = get(gca,'CurrentPoint');    % button down detected
        finalRect = rbbox;                  % return figure units
        point2 = get(gca,'CurrentPoint');    % button up detected
        point1 = point1(1,1:2);              % extract x and y
        point2 = point2(1,1:2);
        p1 = min(point1,point2);             % calculate locations
        offset = abs(point1-point2);         % and dimensions
        xminmax = [p1(1) p1(1)+offset(1)];
        yminmax = [p1(2) p1(2)+offset(2)];

        if strcmp(clicktype,'normal')
            handles.sound.selected_points(:) = false;
        end

        handles.sound.selected_points(...
            handles.sound.f0_corrected > yminmax(1) & ...
            handles.sound.f0_corrected < yminmax(2) & ...
            handles.sound.tcalc > xminmax(1) & ...
            handles.sound.tcalc < xminmax(2) & ...
            handles.sound.f0_corrected > 0) = true;

    elseif strcmp(clicktype,'extend')

        % Save for undoing
        handles.sound.f0_save(end+1,:) = handles.sound.f0;
        handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
        handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;

        % Get initial click point
        point1 = get(gca,'CurrentPoint');    % button down detected
        point1 = point1(1,1:2);
        handles.point1 = point1;
        handles.single_point = false;

        % If no points have previously been selected, select the closest point
        if all(~handles.sound.selected_points)
            a = axis;
            tnorm = (handles.sound.tcalc - a(1))/(a(2)-a(1));
            fnorm = (handles.sound.f0_corrected - a(3))/(a(4)-a(3));
            xnorm = (point1(1) - a(1))/(a(2)-a(1));
            ynorm = (point1(2) - a(3))/(a(4)-a(3));
            [what,where] = min((tnorm-xnorm).^2+(fnorm-ynorm).^2);
            handles.sound.selected_points(where) = true;
            handles.single_point = true;
        end

        guidata(hObject, handles);
        pause(0.05);


        % Set function to track mouse position
        handles.tstart = now;
        guidata(hObject, handles);
        set(gcf,'WindowButtonMotionFcn','p4(''wbmcb'',gcbo,[],guidata(gcbo))');
        set(gcf,'WindowButtonUpFcn','p4(''wbucb'',gcbo,[],guidata(gcbo))');

        handles = update_GUI(handles);

    end
end

guidata(hObject, handles);
handles = update_plot(handles);

% ensures no keyboard input
function nokeyresponse(hObject, ~, ~)

function [f0,acorr,Nshifts,Tshifts] = find_f0_timedomain(A,Fs,Fmin,Fmax,...
    threshold_amp,threshold_ratio,harmonic_deviation,target_f0,target_tol)

% This function calculates the fundamental frequency of a signal in the
% time domain

% OUTPUTS:
%
%   f0:         The interpolated fundamental freqency (Hz)
%   
%   acorr:      The vector of "autocorrelation" values
%
%   Nshifts:    The vector of shift indices associated with the values in
%               acorr
%
%   Tshifts:    The vector of time shifts associated with the values in
%               acorr
%
% INPUTS:
%   
%   A:          The input signal A(t), sampled at an interval Ts
%
%   Fs:         The sample frequency (Hz)
%
%   Fmin:       (Optional) The minimum expected frequency (Hz)
%
%   Fmax:       (Optional) The maximum expected freqency (Hz)
% 
%   threshold_amp: (Optional) The min autocorr amplitude considered peak
%   
%   threshold_ratio: (Optional) The min ratio between max autocorr peak and
%                       a given peak to be considered as a possible peak
%
%   harmonic deviation: (Optional) Tolerance for looking for another peak
%                       at half the frequency of the highest
%                       autocorrelation peak
%
%   target_f0:  (Optional) The probable target frequency
%
%   target_tol: (Optional) The max error between target_f0 and current freq
%                       (0.1 = 10%)
%

if nargin < 9 || isempty(target_tol)
    target_tol = 0.1;
end
if nargin < 8 || isempty(target_f0)
    target_f0 = [];
end
if nargin < 7 || isempty(harmonic_deviation)
    harmonic_deviation = 0.03;
end
if nargin < 6 || isempty(threshold_ratio)
    threshold_ratio = 0.7;
end
if nargin < 5 || isempty(threshold_amp)
    threshold_amp = 0.3;
end
if nargin < 4 || isempty(Fmax)
    Fmax = 800;
end
if nargin < 3 || isempty(Fmin)
    Fmin = 40;
end

f0 = 0;
acorr = [];
Nshifts = [];
Tshifts = [];

% Calculate index of min and max shift
Nshift_min = floor(1/(Fmax/Fs));
Nshift_max = ceil(1/(Fmin/Fs));
if Nshift_max+Nshift_min > length(A) || Nshift_min >= Nshift_max
    disp('Error in find_f0_timedomain.m: Chunk size must be larger!')
    return
end

% Pull out the chunk of the signal and calculate a scale factor
% The scale factor is the probable maximum value
Achunk = A(1:Nshift_max);
scale_factor = sum(abs(Achunk))*2;

% Calculate a vector of all shifts
Nshifts = Nshift_min:Nshift_max;
Tshifts = Nshifts / Fs;


% Shift and calculate a makeshift autocorrelation
acorr = zeros(1,length(Nshifts));
index = 1;
for Nshift = Nshifts
    
    % Break out if we are shifting the chunk beyond the end of the vector A
    if Nshift_max + Nshift > length(A)
        Nshifts = Nshifts(1:index-1);
        Tshifts = Tshifts(1:index-1);
        acorr = acorr(1:index-1);
        break
    end
    
    % Calculate makeshift autocorrelation
    acorr(index) = 1-sum(abs(Achunk - A([1:Nshift_max] + Nshift)))./scale_factor; 
    index = index + 1;
    
end

% Try to make the autocorrelation "level" and with the minimum at zero
acorr = acorr - polyval(polyfit(Nshifts,acorr,1),Nshifts);
acorr = acorr - min(acorr);

% Find a list of all of the peaks above the threshold
[maxX,maxY,Imax,maxX_fit,peakX,peakY] = peakfind(Tshifts,acorr,5,1/Fmax,[],'',false);
% Keep only those harmonics that are within a certain height of the largest
% peak and whose height is above the threshold amplitude

% Make a counter to keep track of which peaks remain after each criteria is
% applied below...
Ipeak = 1:length(maxX_fit);

% Keep those whose heights are above the threshold
Ipeak = Ipeak(maxY > threshold_ratio*max(maxY) & maxY > threshold_amp);
maxX_fit = maxX_fit(Ipeak);

% If only 1 or zero peaks remain, return
N = length(Ipeak);

if N < 2
    return
end
    
% Keep only those harmonics whose frequency is less than a certain
% deviation around a multiple of the fundamental
[maxX_fit,Ix] = sort(maxX_fit);
Ipeak = Ipeak(Ix);
possible_f0 = zeros(1,N);
for n = 1:N-1
    if any(abs(1 - maxX_fit(n+1:N)./(2*maxX_fit(n))) < harmonic_deviation);
        if isempty(target_f0)
            f0 = 1./maxX_fit(n);
            Ipeak = Ipeak(n);
            break
        else
            possible_f0(n) = 1./maxX_fit(n);
        end
    end
end

% Now keep only the harmonic that is closest to the desired harmonic
if isempty(target_f0)

else
    f0_error = abs(possible_f0./target_f0 - 1);
    [what,where] = min(f0_error);
    if what < target_tol
        f0 = possible_f0(where);
        Ipeak = Ipeak(where);
    end
end

% Now "fine tune" the frequency around the peak using a parabolic fit
if f0 > 0
    p = polyfit(peakX(Ipeak,:),peakY(Ipeak,:),2);
    Xfit = -p(2)/2/p(1);
    f0 = 1/Xfit;
end



function [maxX,maxY,Imax,maxX_fit,peakX,peakY] = peakfind(X,Y,Nmax,Xsep,Ymin,sortmethod,peakfit)


% This function finds the maxima of the given function, and
% returns the first Nvalues, sorted in decending order
%
% INPUTS:
%
% X = a vector of X-values for the signal
%
% Y = a fector of Y-values for the signal
%
% Nmax = maximum number of values to return
%
% Xsep = the minimum acceptable separation between "peaks"
% 
% Ymin = the minimum Y value to return
%
% sortmethod = how to sort the peaks...
%               'maxY' = Sort in decending order starting with maximum Y-valued peak
%               'minX' = Sort in ascending order starting with minimum X-valued peak
%
% peakfit = set to 'true' to do a parabolic fit and refine the maxX values
%
% OUTPUTS:
%
% maxX = the X-coordinates of all of the N peaks
%
% maxY = the Y-heights of all of the N peaks
%
% Imax = indices of all of the N peaks
%
% maxX_fit = a more refined list of the X-coordinates based on a parabolic
%           fit to each peak
%
% peakX = a Nx5 matrix with each row containing the 5 X-values around
%           each peak
%
% peakY = a Nx5 matrix with each row containing the 5 Y-values around each
%           peak

if nargin < 7 || isempty(peakfit)
    peakfit = true;
end
if nargin < 6 || isempty(sortmethod)
    sortmethod = 'maxY';
end
if nargin < 5
    Ymin = [];
end
if nargin < 4 || isempty(Xsep)
    Xsep = 0;
end
if nargin < 3 || isempty(Nmax)
    Nmax = 10;
end


% Differentiate, and find change in sign
Ydiff = diff(Y);
Imax1 = find(Ydiff(1:end-1) > 0 & Ydiff(2:end) < 0);
Imax2 = find(Ydiff(1:end-2) > 0 & Ydiff(2:end-1) == 0 & Ydiff(3:end) < 0);
Imax3 = find(Ydiff(1:end-3) > 0 & Ydiff(2:end-2) == 0 & Ydiff(3:end-1) == 0 & Ydiff(4:end) < 0);

% Concatenate all of the possible peaks
maxX = [X(Imax1+1) X(Imax2+1) X(Imax3+2)];
maxY = [Y(Imax1+1) Y(Imax2+1) Y(Imax3+2)];
Imax = [(Imax1+1) (Imax2+1) (Imax3+2)];

% Re-sort if desired
if strcmp(sortmethod,'minX')
    [maxX,IX] = sort(maxX,'ascend');
    maxY = maxY(IX);
    Imax = Imax(IX);
elseif strcmp(sortmethod,'maxY')
    [maxY,IY] = sort(maxY,'descend');
    maxX = maxX(IY);
    Imax = Imax(IY);
end

% Remove any peaks that are below the min peak value allowed
if ~isempty(Ymin)
    Irem = find(maxY < Ymin);
    maxY(Irem) = [];
    maxX(Irem) = [];
    Imax(Irem) = [];
end

% Remove any peaks that are closer than the minimum allowed peak seperation
% in X
if Xsep > 0
    Iremove = zeros(1,length(maxX));
    I = 0;
    for n = 2:length(maxX)
        if abs(maxX(n) - maxX(n-1)) < Xsep
            maxX(n) = maxX(n-1);
            maxY(n) = maxY(n-1);
            I = I+1;
            Iremove(I) = n;
        end
    end
    Iremove = Iremove(1:I);
    maxX(Iremove) = [];
    maxY(Iremove) = [];
    Imax(Iremove) = [];
end

% Remove any peaks beyond the max number to return
if length(maxX) > Nmax
    maxX = maxX(1:Nmax);
    maxY = maxY(1:Nmax);
    Imax = Imax(1:Nmax);
end

maxX_fit = maxX;
peakX = zeros(length(Imax),5);
peakY = zeros(length(Imax),5);
% Do a 2nd-order poly fit around the peak to pinpoint the peak location
for n = 1:length(Imax);
    if Imax(n) > 3 && Imax(n) < length(X)-2
        I = Imax(n)+[-2:2];
        if peakfit
            p = polyfit(X(I),Y(I),2);
            maxX_fit(n) = -p(2)/2/p(1);
        end
        peakX(n,1:5) = X(I);
        peakY(n,1:5) = Y(I);
    end
end

function [indices,freqs,notes] = get_scale(~)

indices = [];

indices = [indices-12 indices indices+12 indices+24 indices+36];
indices = indices(indices >=0 & indices <=48);

[indices,freqs,notes] = get_note_matrix(indices);

function [indices,freqs,notes] = get_note_matrix(note)

if nargin < 1
    note = '';
end

indices = (0:48)';
freqs = 55.*2.^(indices./12);
notes = {
    'A1'
    'A#/Bb1'
    'B1'
    'C1'
    'C#/Db1'
    'D1'
    'D#/Eb1'
    'E1'
    'F1'
    'F#/Gb1'
    'G1'
    'G#/Ab1'
    'A2'
    'A#/Bb2'
    'B2'
    'C2'
    'C#/Db2'
    'D2'
    'D#/Eb2'
    'E2'
    'F2'
    'F#/Gb2'
    'G2'
    'G#/Ab2'
    'A3'
    'A#/Bb3'
    'B3'
    'C3'
    'C#/Db3'
    'D3'
    'D#/Eb3'
    'E3'
    'F3'
    'F#/Gb3'
    'G3'
    'G#/Ab3'
    'A4'
    'A#/Bb4'
    'B4'
    'C4'
    'C#/Db4'
    'D4'
    'D#/Eb4'
    'E4'
    'F4'
    'F#/Gb4'
    'G4'
    'G#/Ab4'
    'A5'
    };

if ~isempty(note)
    if ischar(note)
        I = find(strcmp(notes,note));
        if ~isempty(I)
            indices = indices(I);
            freqs = freqs(I);
            notes = notes(I);
        else
            indices = [];
            freqs = [];
            notes = '';
        end
    else
        note = note(note>=0 & note<=48);
        indices = note;
        freqs = freqs(note+1);
        notes = notes(note+1);
    end
end
