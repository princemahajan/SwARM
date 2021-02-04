function varargout = SwARM(varargin)
% SWARM MATLAB code for SwARM.fig
%      SWARM, by itself, creates a new SWARM or raises the existing
%      singleton*.
%
%      H = SWARM returns the handle to a new SWARM or the handle to
%      the existing singleton*.
%
%      SWARM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SWARM.M with the given input arguments.
%
%      SWARM('Property','Value',...) creates a new SWARM or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SwARM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SwARM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SwARM

% Last Modified by GUIDE v2.5 01-Apr-2018 17:46:06

% Begin initialization code - DO NOT EDIT

% Author: Bharat Mahajan (https://github.com/princemahajan)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SwARM_OpeningFcn, ...
                   'gui_OutputFcn',  @SwARM_OutputFcn, ...
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

% --- Executes just before SwARM is made visible.
function SwARM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SwARM (see VARARGIN)

% Choose default command line output for SwARM
handles.output = hObject;

% Display Logo

% Initialize with the saved settings

% SwARM Configuration file. It MUST EXIST WITH PROPER FORMAT
handles.SwarmConfFile = 'SwARMConf.mat';
handles.SwarmConfStruct = 'SwARMConf';

% Load SwARM Settings Data Structure
load(handles.SwarmConfFile,handles.SwarmConfStruct);

% save configuration in GUIDE handle
handles.SwarmConf = SwARMConf;

% Init Swarm Run structure
%handles.SwarmRun.MJD_EPOCH = 2400000.5;
handles.SwarmRun.GMAT_MJD_EPOCH = 2430000.0; % See GMAT documentation
handles.SwarmRun.ResultsAvailable = false; % no results at this time

% Initially turn off formation ICs update status
handles.SwarmRun.DepModified = false;

% GMAT
handles.SwarmRun.GMATChStatesFile = strcat(SwARMConf.GMAT.ExePath, '\\output\\ChiefStates.txt'); 
handles.SwarmRun.GMATDepStatesFile = strcat(SwARMConf.GMAT.ExePath, '\\output\\DeputyStates.txt'); 
handles.SwarmRun.GMATChElemFile = strcat(SwARMConf.GMAT.ExePath, '\\output\\ChiefOE.txt');
handles.SwarmRun.GMATDepElemFile = strcat(SwARMConf.GMAT.ExePath, '\\output\\DepOE.txt');

handles.SwarmRun.GMATFuncFile = strcat(SwARMConf.GMAT.ExePath, '\\userfunctions\\gmat\\GetCDStates.gmf');

% Initualize the GUI with saved settings
handles = LoadGUI(handles);

% all buttons and plots disable except Init
set(handles.button_Prop,'Enable','off');
set(handles.button_Save,'Enable','off');
set(handles.button_Results,'Enable','off');
set(handles.button_UpdateFF,'Enable','off');
set(handles.AbsMenu,'Enable','off');
set(handles.RelMenu,'Enable','off');

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes SwARM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SwARM_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function text_coefffile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_coefffile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_coefffile_Callback(hObject, eventdata, handles)
% hObject    handle to text_coefffile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_coefffile as text
%        str2double(get(hObject,'String')) returns contents of text_coefffile as a double



% % --- Executes during object creation, after setting all properties.
% function volume_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to volume (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



% function volume_Callback(hObject, eventdata, handles)
% % hObject    handle to volume (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of volume as text
% %        str2double(get(hObject,'String')) returns contents of volume as a double
% volume = str2double(get(hObject, 'String'));
% if isnan(volume)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end

% % Save the new volume value
% handles.metricdata.volume = volume;
% guidata(hObject,handles)

% --- Executes on button press in button_Init.
function button_Init_Callback(hObject, eventdata, handles)
% hObject    handle to button_Init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read all the settings and generate initial conditions and other stuff for
% numerical as well as analytical propagation

% Disable all the buttons until initialization is done
set(handles.button_Init,'Enable','off');
set(handles.button_Prop,'Enable','off');
set(handles.button_Save,'Enable','off');
set(handles.button_Results,'Enable','off');
set(handles.button_UpdateFF,'Enable','off');


% Disable Results section
set(handles.AbsMenu,'Enable','off');
set(handles.RelMenu,'Enable','off');

% No results available until Generate Results button is used
handles.SwarmRun.ResultsAvailable = false;

% set progress bar 0 %
set(handles.progressbar,'Value',0);
drawnow;

% Read the GUI
SwarmConf = ReadGUI(handles);

% Use the new Swarm configuration from now on
handles.SwarmConf = SwarmConf;

% turn off this flag
handles.SwarmRun.DepModified = false;

% Swarm run data structure
lunit = handles.SwarmConf.GravModel.Re ;
tunit = sqrt(handles.SwarmConf.GravModel.Re^3/handles.SwarmConf.GravModel.mu);
vunit = lunit/tunit;
mu = 1;
Re = 1;

% Greenwich Mean Sidereal Time and Earth rotational speed at given UTC JD
JDUTCt0 = handles.SwarmConf.Time.MJDUTCt0 + handles.SwarmRun.GMAT_MJD_EPOCH;
[GMST0, we_d] = JD2GMST( JDUTCt0 );
we = we_d/1*tunit;

% Unnormalize C and S coefficients
[Clm,Slm] = DenormCS(handles.SwarmConf.GravCoeffFN, max([2, handles.SwarmConf.GravModel.IC.n, ...
                                                         handles.SwarmConf.GravModel.AST.n,...
                                                         handles.SwarmConf.GravModel.NP.n]));
C20 = Clm(3,1);

% set progress bar 25 %
set(handles.progressbar,'Value',0.25);
drawnow;

% Compute Chief initial osculating states
sma = SwarmConf.ChiefIC.SMA/lunit;
ecc = SwarmConf.ChiefIC.ECC;
inc = deg2rad(SwarmConf.ChiefIC.INC);
raan = deg2rad(SwarmConf.ChiefIC.RAAN);
aop = deg2rad(SwarmConf.ChiefIC.AOP);
ma = deg2rad(SwarmConf.ChiefIC.M);

ChiefEq0 = [sma; ma+raan+aop; tan(inc/2)*cos(raan); tan(inc/2)*sin(raan);ecc*cos(aop+raan); ecc*sin(aop+raan)];

JacobianOn = 0;
if strcmpi(handles.SwarmConf.STMType, 'None') ~= 1
    JacobianOn = 1;
end

% Compute Absolute Initial Mean States for Chief
[EqCMX0,~, iD0,~] = EqMean2Osc2(ChiefEq0, GMST0, mu, Re, we, handles.SwarmConf.GravModel.IC.n, handles.SwarmConf.GravModel.IC.m, ...
                            Clm, Slm, handles.SwarmConf.tol, handles.SwarmConf.quadtol, ...
                            true, handles.SwarmConf.J2MEXON, JacobianOn, false, handles.SwarmConf.GravModel.AST.QuadTesseralsOn );

% Chief Initial States (osculating)
ChX0 = Eqn2RV(ChiefEq0, mu, handles.SwarmConf.tol, false);

% set progress bar 50%
set(handles.progressbar,'Value',0.5);
drawnow;

% Deputy Initial Conditions
DepX0 = 0;
if JacobianOn == 1
    rhox = SwarmConf.RM.RHOX/lunit;
    rhoy = SwarmConf.RM.RHOY/lunit;
    rhoz = SwarmConf.RM.RHOZ/lunit;
    alphax = deg2rad(SwarmConf.RM.ALPHAX);
    alphaz = deg2rad(SwarmConf.RM.ALPHAZ);
    
    dMEqn0 = FormationdEqn(rhox, rhoy, rhoz, alphax, alphaz, EqCMX0,...
                handles.SwarmConf.RM.DriftCond, mu, Re, -Clm(2:end,1), handles.SwarmConf.tol);    

    [EqDOX0,~, ~,~] = EqMean2Osc2(EqCMX0 + dMEqn0, GMST0, mu, Re, we, handles.SwarmConf.GravModel.IC.n, handles.SwarmConf.GravModel.IC.m, ...
                            Clm, Slm, handles.SwarmConf.tol, handles.SwarmConf.quadtol, ...
                            false, handles.SwarmConf.J2MEXON, false, false, handles.SwarmConf.GravModel.AST.QuadTesseralsOn );
    DepX0 = Eqn2RV(EqDOX0, mu, handles.SwarmConf.tol, false);
end

% Build time vector
tp = 2*pi*sqrt(sma^3/mu);
t_data = handles.SwarmConf.Time.NumOrbits*tp;
numpoints = handles.SwarmConf.Time.PointsPerOrbit*handles.SwarmConf.Time.NumOrbits;
StartTime = handles.SwarmConf.Time.DataStartTime/tunit;

if StartTime > 0
    tspan = [0,linspace(StartTime, StartTime + t_data, numpoints)];
else
    tspan = linspace(0,t_data, numpoints);
end

% Numerical Simulation
if strcmpi(handles.SwarmConf.NumPropType, 'None') ~= 1

    % Set up GMAT
    if strcmpi(handles.SwarmConf.NumPropType, 'GMAT') == 1    
    
        % GMAT Command
        handles.SwarmRun.GMATcmd = sprintf('%s --run --exit "%s"',...
                                    strcat(handles.SwarmConf.GMAT.ExePath,'\\bin\\GMAT'), ...
                                    handles.SwarmConf.GMAT.ScriptFN);

        % Save parameters for GMAT simulation
        GPropTime = tspan(end)*tunit;
        GChiefOE = ChX0*lunit;
        GChiefOE(4:6) = GChiefOE(4:6)/tunit;
        if JacobianOn ~= 1
            GDepOE = GChiefOE;
        else
            GDepOE = DepX0*lunit;
            GDepOE(4:6) = GDepOE(4:6)/tunit;
        end
        
        % Create GMAT function and save it
        GMATFuncFID = fopen(handles.SwarmRun.GMATFuncFile,'w');
        if GMATFuncFID == -1
            disp(strcat('SwARM: check GMAT function folder path ', handles.SwarmRun.GMATFuncFile));
        end
        
        fprintf(GMATFuncFID,strcat(' function [ChX, DepX, TOF] = GetCDStates \n\n Create Array ChX[6,1] DepX[6,1]; \n\n',...
            ' Create Variable TOF;\n\n BeginMissionSequence \n\n',...
            ' ChX(1,1) = %3.12f; \n ChX(2,1) = %3.12f; \n ChX(3,1) = %3.12f;',...
            '\n ChX(4,1) = %3.12f; \n ChX(5,1) = %3.12f; \n ChX(6,1) = %3.12f; \n\n',...
            ' DepX(1,1) = %3.12f; \n DepX(2,1) = %3.12f; \n DepX(3,1) = %3.12f; \n',...
            ' DepX(4,1) = %3.12f; \n DepX(5,1) = %3.12f; \n DepX(6,1) = %3.12f; \n\n',...
            ' TOF = %3.12f; \n'),GChiefOE(1),GChiefOE(2),GChiefOE(3),GChiefOE(4),GChiefOE(5),GChiefOE(6),...
            GDepOE(1),GDepOE(2),GDepOE(3),GDepOE(4),GDepOE(5),GDepOE(6),GPropTime);
        
        % close GMAT function file
        fclose(GMATFuncFID);
        
    end
end

% set progress bar 75%
set(handles.progressbar,'Value',0.75);
drawnow;

% Compute GMST for future computations
GMSTarr = zeros(length(tspan),1);
for ctr = 1:length(tspan)
    [GMST, ~] = JD2GMST( JDUTCt0 + (tspan(ctr)-tspan(1))*tunit/86400);
    GMSTarr(ctr) = GMST;
end

% Compute Initial geometric transformation for STM
iSig0 = 0;
if JacobianOn == 1

    % Geometric Transformation for propagation
    Sig0 = DiffEq2Curv(ChiefEq0, GMST0, handles.SwarmConf.GravModel.IC.n, handles.SwarmConf.GravModel.IC.m, ...
                        true, mu, Re, Clm, Slm, handles.SwarmConf.tol);

    iSig0 = Sig0^-1;
    
end

% Update Swarm Run
handles.SwarmRun.lunit = lunit;
handles.SwarmRun.tunit = tunit;
handles.SwarmRun.vunit = vunit;
handles.SwarmRun.mu = mu;
handles.SwarmRun.Re = Re;
handles.SwarmRun.we = we;
handles.SwarmRun.GMST0 = GMST0;
handles.SwarmRun.tspan = tspan;
handles.SwarmRun.Clm = Clm;
handles.SwarmRun.Slm = Slm;
handles.SwarmRun.C20 = C20;
handles.SwarmRun.JDUTCt0 = JDUTCt0;
handles.SwarmRun.GMSTarr = GMSTarr;

handles.SwarmRun.EqCOX0 = ChiefEq0;
handles.SwarmRun.EqCMX0 = EqCMX0;
handles.SwarmRun.ChX0 = ChX0;
handles.SwarmRun.DepX0 = DepX0;
handles.SwarmRun.iD0 = iD0;
handles.SwarmRun.iSig0 = iSig0;


disp('SwARM: Init Completed...');

% set progress bar 100 %
set(handles.progressbar,'Value',1);

% Reenable all the buttons except results and update button
set(handles.button_Init,'Enable','on');
set(handles.button_Prop,'Enable','on');
set(handles.button_Save,'Enable','on');
set(handles.button_Results,'Enable','off');
set(handles.button_UpdateFF,'Enable','off');

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_Prop.
function button_Prop_Callback(hObject, eventdata, handles)
% hObject    handle to button_Prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable all the buttons until propagation is done
set(handles.button_Init,'Enable','off');
set(handles.button_Prop,'Enable','off');
set(handles.button_Save,'Enable','off');
set(handles.button_Results,'Enable','off');
set(handles.button_UpdateFF,'Enable','off');

% Disable Results section
set(handles.AbsMenu,'Enable','off');
set(handles.RelMenu,'Enable','off');

% set progressbar to 0%
set(handles.progressbar,'Value',0);
drawnow;

tspan = handles.SwarmRun.tspan;
nPoints = length(tspan);

if strcmpi(handles.SwarmConf.STMType,'None') ~= 1
    STM_on = 1;
else
    STM_on = 0;
end

tic;

% Numerical Propagation - GMAT or MATLAB if needed
if strcmpi(handles.SwarmConf.NumPropType,'GMAT') == 1
    
    % Run GMAT
    status = system(handles.SwarmRun.GMATcmd);
    disp(['SwARM: GMAT call status = ' num2str(status)]);
    
    % Load Chief and Deputy States and nondimensionalize
    ChStates = dlmread(handles.SwarmRun.GMATChStatesFile);
    ChElems = dlmread(handles.SwarmRun.GMATChElemFile);

    GT = ChStates(:,1)/handles.SwarmRun.tunit;
    GORV = [ChStates(:,2:4)/handles.SwarmRun.lunit, ChStates(:,5:7)/handles.SwarmRun.vunit]; 
    GChElem = [ChElems(:,2)/handles.SwarmRun.lunit,ChElems(:,3),deg2rad(ChElems(:,4:7))];
    
    % Interpolate states at specified times using cubic splines
    ORV = zeros(nPoints,6);
    NumCOE = ORV;
    for ctr = 1:nPoints
        ORV(ctr,:) = spline(GT',GORV',handles.SwarmRun.tspan(ctr))';
        NumCOE(ctr,:) = spline(GT',GChElem',handles.SwarmRun.tspan(ctr))';
    end

    % Deputy Numerical states
    ODRV = 0;
    NumDOE = 0;
    if STM_on == 1
        DepStates = dlmread(handles.SwarmRun.GMATDepStatesFile);
        DepElems = dlmread(handles.SwarmRun.GMATDepElemFile);
        
        GODRV = [DepStates(:,2:4)/handles.SwarmRun.lunit, DepStates(:,5:7)/handles.SwarmRun.vunit];
        GDepElem = [DepElems(:,2)/handles.SwarmRun.lunit,DepElems(:,3),deg2rad(DepElems(:,4:7))];
        
        ODRV = zeros(nPoints,6);
        NumDOE = ODRV;
        for ctr = 1:nPoints
            ODRV(ctr,:) = spline(GT',GODRV',handles.SwarmRun.tspan(ctr))';
            NumDOE(ctr,:) = spline(GT',GDepElem',handles.SwarmRun.tspan(ctr))';
        end
    end
   
    handles.SwarmRun.NumCRV = ORV;
    handles.SwarmRun.NumCOE = NumCOE;
    handles.SwarmRun.NumDRV = ODRV;
    handles.SwarmRun.NumDOE = NumDOE;

elseif strcmpi(handles.SwarmConf.NumPropType,'MATLAB') == 1
    
    % Run Matlab Numerical Propagator
    optn = odeset('RelTol',handles.SwarmConf.tol,'AbsTol',handles.SwarmConf.tol*1e-3);
    tspan = handles.SwarmRun.tspan;
    n = handles.SwarmConf.GravModel.NP.n;
    m = handles.SwarmConf.GravModel.NP.m;
    
    % Chief
    [~,ORV] = ode113(@(T, ORV) GravAcc(T, ORV, tspan(1), handles.SwarmRun.GMST0,handles.SwarmRun.we, ...
                                n, m, handles.SwarmRun.mu,handles.SwarmRun.Re,handles.SwarmRun.Clm,handles.SwarmRun.Slm),...
                                tspan,handles.SwarmRun.ChX0,optn);
                            
    % compute Chief classical elements
    nPoints = length(handles.SwarmRun.tspan);
    NumCOE = zeros(nPoints,6);
    for ctr = 1:nPoints
        oe = OrbitElem(ORV(ctr,1:3)',ORV(ctr,4:6)',handles.SwarmRun.mu);
        NumCOE(ctr,:) = [oe(2:6),oe(end)];
    end
    
    ODRV = 0;
    NumDOE = 0;
    if STM_on == 1
        % Deputy
        [~,ODRV] = ode113(@(T, ODRV) GravAcc(T, ODRV, tspan(1), handles.SwarmRun.GMST0,handles.SwarmRun.we, ...
                                n, m, handles.SwarmRun.mu,handles.SwarmRun.Re,handles.SwarmRun.Clm,handles.SwarmRun.Slm),...
                                tspan,handles.SwarmRun.DepX0,optn);
        NumDOE = zeros(nPoints,6);
        for ctr = 1:nPoints
            oe = OrbitElem(ODRV(ctr,1:3)',ODRV(ctr,4:6)',handles.SwarmRun.mu);
            NumDOE(ctr,:) = [oe(2:6),oe(end)];
        end
    end
    handles.SwarmRun.NumCRV = ORV;
    handles.SwarmRun.NumCOE = NumCOE;
    handles.SwarmRun.NumDRV = ODRV;
    handles.SwarmRun.NumDOE = NumDOE;
end

np_time = toc;
set(handles.text_nptime,'String',np_time);

disp('SwARM: Numerical Propagation done...');

% set progress bar 25 %
set(handles.progressbar,'Value',0.25);
drawnow;

tic;

% Propagate Mean Elements
EqCMX = zeros(nPoints, 6);
mPhi = zeros(6,6,nPoints);

Jcoeff = -handles.SwarmRun.Clm(2:end,1);
for ctr = 1:nPoints
    [OEm, phi] = EqZonalMeanProp(tspan(ctr)-tspan(1), handles.SwarmRun.EqCMX0, ...
                    handles.SwarmConf.GravModel.AST.n, handles.SwarmRun.mu, handles.SwarmRun.Re, Jcoeff);
    EqCMX(ctr,:) = OEm';
    mPhi(:,:,ctr) = phi;
end

mean_time = toc;
set(handles.text_meantime,'String',mean_time);

disp('SwARM: Chief Mean Propagation Done...');

% set progress bar 35 %
set(handles.progressbar,'Value',0.35);
drawnow;

tic;

% propagate chief osculating elements
EqCOX = zeros(nPoints, 6);
Darr = zeros(6,6,nPoints);

USE_MEX = handles.SwarmConf.J2MEXON;

% update progress 4 times
pbupdate = fix(nPoints/4);
for ctr = 1:nPoints
    [EqCOX(ctr,:),~,D,~] = EqMean2Osc2(EqCMX(ctr,:)', handles.SwarmRun.GMSTarr(ctr), handles.SwarmRun.mu, handles.SwarmRun.Re, ...
                            handles.SwarmRun.we, handles.SwarmConf.GravModel.AST.n,handles.SwarmConf.GravModel.AST.m,...
                            handles.SwarmRun.Clm, handles.SwarmRun.Slm, handles.SwarmConf.tol, handles.SwarmConf.quadtol,...
                            false, USE_MEX, STM_on, false , handles.SwarmConf.GravModel.AST.QuadTesseralsOn);
    Darr(:,:,ctr) = D;
    
    if mod(ctr,pbupdate) == 0
        disp(['Chief M2O: ' num2str(ctr)]);
        % set progress bar
        set(handles.progressbar,'Value',ctr/nPoints*50/100 + 0.35);
        drawnow;
    end
end

m2o_time = toc;
set(handles.text_m2otime,'String',m2o_time);

disp('SwARM: Chief Absolute Motion Propagation Done...' );

% Generate Relative Motion STM for Deputy
STM = 0;
if STM_on == true

    STM = zeros(6,6,nPoints);

    for ctr = 1:nPoints
        % Differential Element STM
        STM(:,:,ctr) = Darr(:,:,ctr)*mPhi(:,:,ctr)*handles.SwarmRun.iD0;
        
        if strcmpi(handles.SwarmConf.STMType,'Relative States') == 1
            % Relative states STM
            Sig = DiffEq2Curv(EqCOX(ctr,:)', handles.SwarmRun.GMSTarr(ctr), ...
                handles.SwarmConf.GravModel.AST.n, handles.SwarmConf.GravModel.AST.m, true,...
                handles.SwarmRun.mu, handles.SwarmRun.Re, handles.SwarmRun.Clm, handles.SwarmRun.Slm, handles.SwarmConf.tol);
            STM(:,:,ctr) = Sig*STM(:,:,ctr)*handles.SwarmRun.iSig0;
        end
    end
    
    disp('SwARM: Relative Motion STM Generated...' );
    
end

% Save Data
handles.SwarmRun.EqCMX = EqCMX;
handles.SwarmRun.EqCOX = EqCOX;
handles.SwarmRun.STM = STM;

% set progress bar 100 %
set(handles.progressbar,'Value',1);
drawnow;


% Reenable all the buttons now
set(handles.button_Init,'Enable','on');
set(handles.button_Prop,'Enable','on');
set(handles.button_Save,'Enable','on');
set(handles.button_Results,'Enable','on');
set(handles.button_UpdateFF,'Enable','on');

% Update handles structure
guidata(hObject, handles);


% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if (hObject == handles.english)
%     set(handles.text4, 'String', 'lb/cu.in');
%     set(handles.text5, 'String', 'cu.in');
%     set(handles.text6, 'String', 'lb');
% else
%     set(handles.text4, 'String', 'kg/cu.m');
%     set(handles.text5, 'String', 'cu.m');
%     set(handles.text6, 'String', 'kg');
% end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the button_Prop flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to button_Prop the data.

disp('**********************************');
disp('************* SwARM **************');
disp('**********************************');



function text_degreeIC_Callback(hObject, eventdata, handles)
% hObject    handle to text_degreeIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_degreeIC as text
%        str2double(get(hObject,'String')) returns contents of text_degreeIC as a double


% --- Executes during object creation, after setting all properties.
function text_degreeIC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_degreeIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_degree_Callback(hObject, eventdata, handles)
% hObject    handle to text_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_degree as text
%        str2double(get(hObject,'String')) returns contents of text_degree as a double


% --- Executes during object creation, after setting all properties.
function text_degree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_np_gmat.
function radio_np_gmat_Callback(hObject, eventdata, handles)
% hObject    handle to radio_np_gmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_np_gmat


% --- Executes on button press in radio_np_matlab.
function radio_np_matlab_Callback(hObject, eventdata, handles)
% hObject    handle to radio_np_matlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_np_matlab



function text_scriptfile_Callback(hObject, eventdata, handles)
% hObject    handle to text_scriptfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_scriptfile as text
%        str2double(get(hObject,'String')) returns contents of text_scriptfile as a double


% --- Executes during object creation, after setting all properties.
function text_scriptfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_scriptfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_gmatpath_Callback(hObject, eventdata, handles)
% hObject    handle to text_gmatpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_gmatpath as text
%        str2double(get(hObject,'String')) returns contents of text_gmatpath as a double


% --- Executes during object creation, after setting all properties.
function text_gmatpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_gmatpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_stm_de.
function radio_stm_de_Callback(hObject, eventdata, handles)
% hObject    handle to radio_stm_de (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_stm_de


% --- Executes on button press in radio_stm_relstate.
function radio_stm_relstate_Callback(hObject, eventdata, handles)
% hObject    handle to radio_stm_relstate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_stm_relstate


% --- Executes on button press in english.
function english_Callback(hObject, eventdata, handles)
% hObject    handle to english (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of english


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over english.
function english_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to english (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on english and none of its controls.
function english_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to english (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_rhox_Callback(hObject, eventdata, handles)
% hObject    handle to text_rhox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_rhox as text
%        str2double(get(hObject,'String')) returns contents of text_rhox as a double


% --- Executes during object creation, after setting all properties.
function text_rhox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_rhox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_rhoy_Callback(hObject, eventdata, handles)
% hObject    handle to text_rhoy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_rhoy as text
%        str2double(get(hObject,'String')) returns contents of text_rhoy as a double


% --- Executes during object creation, after setting all properties.
function text_rhoy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_rhoy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function text_rhoz_Callback(hObject, eventdata, handles)
% hObject    handle to text_rhoz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_rhoz as text
%        str2double(get(hObject,'String')) returns contents of text_rhoz as a double


% --- Executes during object creation, after setting all properties.
function text_rhoz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_rhoz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_alphax_Callback(hObject, eventdata, handles)
% hObject    handle to text_alphax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_alphax as text
%        str2double(get(hObject,'String')) returns contents of text_alphax as a double


% --- Executes during object creation, after setting all properties.
function text_alphax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_alphax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_alphaz_Callback(hObject, eventdata, handles)
% hObject    handle to text_alphaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_alphaz as text
%        str2double(get(hObject,'String')) returns contents of text_alphaz as a double


% --- Executes during object creation, after setting all properties.
function text_alphaz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_alphaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_sma_Callback(hObject, eventdata, handles)
% hObject    handle to text_sma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_sma as text
%        str2double(get(hObject,'String')) returns contents of text_sma as a double


% --- Executes during object creation, after setting all properties.
function text_sma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_sma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_ecc_Callback(hObject, eventdata, handles)
% hObject    handle to text_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_ecc as text
%        str2double(get(hObject,'String')) returns contents of text_ecc as a double


% --- Executes during object creation, after setting all properties.
function text_ecc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_inc_Callback(hObject, eventdata, handles)
% hObject    handle to text_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_inc as text
%        str2double(get(hObject,'String')) returns contents of text_inc as a double


% --- Executes during object creation, after setting all properties.
function text_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_raan_Callback(hObject, eventdata, handles)
% hObject    handle to text_raan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_raan as text
%        str2double(get(hObject,'String')) returns contents of text_raan as a double


% --- Executes during object creation, after setting all properties.
function text_raan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_raan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_aop_Callback(hObject, eventdata, handles)
% hObject    handle to text_aop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_aop as text
%        str2double(get(hObject,'String')) returns contents of text_aop as a double


% --- Executes during object creation, after setting all properties.
function text_aop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_aop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_M_Callback(hObject, eventdata, handles)
% hObject    handle to text_M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_M as text
%        str2double(get(hObject,'String')) returns contents of text_M as a double


% --- Executes during object creation, after setting all properties.
function text_M_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_MJDUTCt0_Callback(hObject, eventdata, handles)
% hObject    handle to text_MJDUTCt0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_MJDUTCt0 as text
%        str2double(get(hObject,'String')) returns contents of text_MJDUTCt0 as a double


% --- Executes during object creation, after setting all properties.
function text_MJDUTCt0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_MJDUTCt0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_datat0_Callback(hObject, eventdata, handles)
% hObject    handle to text_datat0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_datat0 as text
%        str2double(get(hObject,'String')) returns contents of text_datat0 as a double


% --- Executes during object creation, after setting all properties.
function text_datat0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_datat0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_numorbits_Callback(hObject, eventdata, handles)
% hObject    handle to text_numorbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_numorbits as text
%        str2double(get(hObject,'String')) returns contents of text_numorbits as a double


% --- Executes during object creation, after setting all properties.
function text_numorbits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_numorbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_pointsperorbit_Callback(hObject, eventdata, handles)
% hObject    handle to text_pointsperorbit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_pointsperorbit as text
%        str2double(get(hObject,'String')) returns contents of text_pointsperorbit as a double


% --- Executes during object creation, after setting all properties.
function text_pointsperorbit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_pointsperorbit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_Save.
function button_Save_Callback(hObject, eventdata, handles)
% hObject    handle to button_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.SwarmRun.ResultsAvailable == true
    
    % save run results
    SwARMRun = handles.SwarmRun;
    
%     lunit = handles.SwarmRun.lunit;
%     tunit = handles.SwarmRun.tunit;
%     tspan = handles.SwarmRun.tspan;
%     
%     iSig0 = handles.SwarmRun.iSig0;
%     MeanOE = handles.SwarmRun.AnalMCOE;
%     OscOE = handles.SwarmRun.AnalCOE;
%     STM = handles.SwarmRun.STM;
%     ChiefRV = handles.SwarmRun.AnalCRV;
%     RelRV = handles.SwarmRun.AnalCurvRV;
%     
%     save(handles.SwarmConf.ResultsFN,'lunit','tunit','tspan','iSig0','MeanOE','OscOE','STM','ChiefRV','RelRV');
    save(handles.SwarmConf.ResultsFN,'SwARMRun');

end

% Save configurations on disk
SwARMConf = handles.SwarmConf;
save(handles.SwarmConfFile,handles.SwarmConfStruct);

disp('SwARM results saved...');

% --- Executes on selection change in AbsMenu.
function AbsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to AbsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AbsMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AbsMenu
absplot = get(handles.AbsMenu,'Value');
handles = AbsPlotData(absplot, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function AbsMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AbsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RelMenu.
function RelMenu_Callback(hObject, eventdata, handles)
% hObject    handle to RelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RelMenu

% Display Relative Results plots based on the selection in pop menu
relplot = get(handles.RelMenu,'Value');
handles = RelPlotData(relplot,handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_Results.
function button_Results_Callback(hObject, eventdata, handles)
% hObject    handle to button_Results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable all the buttons until results are generated
set(handles.button_Init,'Enable','off');
set(handles.button_Prop,'Enable','off');
set(handles.button_Save,'Enable','off');
set(handles.button_Results,'Enable','off');
set(handles.button_UpdateFF,'Enable','off');

% set progress bar 0 %
set(handles.progressbar,'Value',0);
drawnow;

tspan = handles.SwarmRun.tspan;
nPoints = length(tspan);

% Generate Chief States from Theory
AnalCRV = zeros(nPoints,6);
AnalCOE = zeros(nPoints,6);
AnalMCOE = zeros(nPoints,6);
for ctr = 1:nPoints
    AnalCRV(ctr,:) = Eqn2RV(handles.SwarmRun.EqCOX(ctr,:)' ,handles.SwarmRun.mu, handles.SwarmConf.tol, false)';
    oe = OrbitElem(AnalCRV(ctr,1:3)', AnalCRV(ctr,4:6)',handles.SwarmRun.mu);
    AnalCOE(ctr,:) = [oe(2:6),oe(end)];
    % generate Chief mean orbital elements
    AnalMCRV = Eqn2RV(handles.SwarmRun.EqCMX(ctr,:)' ,handles.SwarmRun.mu, handles.SwarmConf.tol, false)';
    oe = OrbitElem(AnalMCRV(1:3)', AnalMCRV(4:6)',handles.SwarmRun.mu);
    AnalMCOE(ctr,:) = [oe(2:6),oe(end)];
end

% set progress bar 25 %
set(handles.progressbar,'Value',0.25);
drawnow;

% Generate Numerically propagated Relative Motion
NumCurvRV = zeros(nPoints,6);
if strcmpi(handles.SwarmConf.NumPropType, 'None') ~= 1 && strcmpi(handles.SwarmConf.STMType, 'None') ~= 1
   
    % Compute Hill's frame Relative states
    for ctr = 1:nPoints
        NumCurvRV(ctr,1:6) = ECI2Curv( handles.SwarmRun.NumCRV(ctr,:)', ...
            handles.SwarmRun.NumDRV(ctr,:)', false )';
    end
    
end

% set progress bar 50 %
set(handles.progressbar,'Value',0.5);
drawnow;

% Generate Deputy States from Theory
AnalCurvRV = 0;
DiffCOE = 0;
if strcmpi(handles.SwarmConf.STMType, 'None') ~= 1
    
    rhox = handles.SwarmConf.RM.RHOX/handles.SwarmRun.lunit;
    rhoy = handles.SwarmConf.RM.RHOY/handles.SwarmRun.lunit;
    rhoz = handles.SwarmConf.RM.RHOZ/handles.SwarmRun.lunit;
    alphax = deg2rad(handles.SwarmConf.RM.ALPHAX);
    alphaz = deg2rad(handles.SwarmConf.RM.ALPHAZ);
    
    % differential equinoctial elements
    dMEqn0 = FormationdEqn(rhox, rhoy, rhoz, alphax, alphaz,handles.SwarmRun.EqCMX0,...
                handles.SwarmConf.RM.DriftCond, handles.SwarmRun.mu, handles.SwarmRun.Re, ...
                -handles.SwarmRun.Clm(2:end,1), handles.SwarmConf.tol);    

    [EqDOX0,~, ~,~] = EqMean2Osc2(handles.SwarmRun.EqCMX0 + dMEqn0, handles.SwarmRun.GMST0,...
                            handles.SwarmRun.mu, handles.SwarmRun.Re, handles.SwarmRun.we, ...
                            handles.SwarmConf.GravModel.IC.n, handles.SwarmConf.GravModel.IC.m, ...
                            handles.SwarmRun.Clm, handles.SwarmRun.Slm, handles.SwarmConf.tol, handles.SwarmConf.quadtol, ...
                            false, handles.SwarmConf.J2MEXON, false, false, handles.SwarmConf.GravModel.AST.QuadTesseralsOn );
            
    % Compute differential true longitude
    EqCOX0 = handles.SwarmRun.EqCOX0;
    psic0 = Mean2TrueLong(EqCOX0(2), EqCOX0(5), EqCOX0(6), handles.SwarmConf.tol);
    psid0 = Mean2TrueLong(EqDOX0(2), EqDOX0(5), EqDOX0(6), handles.SwarmConf.tol);
    dPsi0 = atan((tan(psid0)-tan(psic0))/(1 + tan(psid0)*tan(psic0)));
    
    % compute Deputy ICs
    dEqn0 = EqDOX0 - EqCOX0; 
    if strcmpi(handles.SwarmConf.STMType, 'Relative States') == 1
        DiffX0 = handles.SwarmRun.iSig0^-1*[dEqn0(1); dEqn0(2); dEqn0(3:6,1)];
    else
        DiffX0 = dEqn0;
    end
    
    % set progress bar 75 %
    set(handles.progressbar,'Value',0.75);
    drawnow;
    
    % Propagate Deputy States
    AnalCurvRV = zeros(nPoints,6);
    DiffCOE = AnalCurvRV;
    for ctr = 1:nPoints
        % propagate relative motion
        DiffX = handles.SwarmRun.STM(:,:,ctr)*DiffX0;
        
        % if differential elements STM then compute curvilienar states by
        % nonlinear transformation
        if strcmpi(handles.SwarmConf.STMType,'Differential Elements') == 1
            EqDOX = handles.SwarmRun.EqCOX(ctr,:) + DiffX';
            AnalDRV = Eqn2RV(EqDOX ,handles.SwarmRun.mu, handles.SwarmConf.tol, false)';
            oe = OrbitElem(AnalDRV(1:3), AnalDRV(4:6),handles.SwarmRun.mu);
            AnalDOE = [oe(2:6),oe(end)];
            DiffCOE(ctr,:) = AnalDOE - AnalCOE(ctr,:);
            AnalCurvRV(ctr,:) = ECI2Curv( AnalCRV(ctr,:)', AnalDRV', false );
        else
            AnalDRV = ECI2Curv( AnalCRV(ctr,:)', DiffX, true );
            oe = OrbitElem(AnalDRV(1:3), AnalDRV(4:6),handles.SwarmRun.mu);
            AnalDOE = [oe(2:6),oe(end)];
            DiffCOE(ctr,:) = AnalDOE - AnalCOE(ctr,:);
            AnalCurvRV(ctr,:) = DiffX';
        end
    end
end

% store results for store on disk
handles.SwarmRun.AnalCRV = AnalCRV;
handles.SwarmRun.AnalMCOE = AnalMCOE;
handles.SwarmRun.AnalCOE = AnalCOE;
handles.SwarmRun.NumCurvRV = NumCurvRV;
handles.SwarmRun.AnalCurvRV = AnalCurvRV;
handles.SwarmRun.AnalDiffCOE = DiffCOE;

% Display Absolute Results plots based on the selection in pop menu
set(handles.AbsMenu,'Enable','on');
absplot = get(handles.AbsMenu,'Value');
handles = AbsPlotData(absplot,handles);

if strcmpi(handles.SwarmConf.STMType, 'None') ~= 1
    % Display Relative Results plots based on the selection in pop menu
    set(handles.RelMenu,'Enable','on');
    relplot = get(handles.RelMenu,'Value');
    handles = RelPlotData(relplot,handles);
end

% disp('SwARM: Results generated...');
% Results available to save in mat file now
handles.SwarmRun.ResultsAvailable = true;

% set progress bar 100 %
set(handles.progressbar,'Value',1);
drawnow;

% Renable all the buttons now
set(handles.button_Init,'Enable','on');
set(handles.button_Prop,'Enable','on');
set(handles.button_Save,'Enable','on');
set(handles.button_Results,'Enable','on');
set(handles.button_UpdateFF,'Enable','on');

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton6.
function pushbutton6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function text_driftcond_Callback(hObject, eventdata, handles)
% hObject    handle to text_driftcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_driftcond as text
%        str2double(get(hObject,'String')) returns contents of text_driftcond as a double


% --- Executes during object creation, after setting all properties.
function text_driftcond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_driftcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_tol_Callback(hObject, eventdata, handles)
% hObject    handle to text_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_tol as text
%        str2double(get(hObject,'String')) returns contents of text_tol as a double


% --- Executes during object creation, after setting all properties.
function text_tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_quadtol_Callback(hObject, eventdata, handles)
% hObject    handle to text_quadtol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_quadtol as text
%        str2double(get(hObject,'String')) returns contents of text_quadtol as a double


% --- Executes during object creation, after setting all properties.
function text_quadtol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_quadtol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg1 = '\fontsize{11} \color[rgb]{.45,.31,.16} \bf SwARM: Software for Absolute and Relative Motion (Satellite Theories)';
msg2 = {'\fontsize{10} \color[rgb]{.45,.31,.16} \it Author: Bharat Mahajan', 'Contact: prince@tamu.edu, bharat.mahajan@hotmail.com',...
        'Copyright (C) 2018 by Bharat Mahajan'};
msg3 = {strcat('\fontsize{10}     SwARM is a second-order analytic propagator for absolute and relative motion of satellites in the',...
       ' presence of gravitational harmonic perturbations, with the oblateness coefficient considered as the first-order and the rest',...
       ' of the harmonics as the second-order perturbations. The oblateness'' secular effects are added up to O(3) and periodic effects',...
       ' up to O(2). For the rest of the zonal harmonics, secular and short-periodic effects are included up to O(2) and long-periodic ',...
       ' up to O(1). For tesseral harmonics, only the short-periodic effects have been included up to O(2). The theory is completely',...
       ' closed-form in the eccentricity with the tesseral effects computed using numerical quadrature. There is an option to use eccentricity',...
       ' expansions for computing the tesseral effects using a series with terms included up to O(7) of the eccentricity.'),...
       strcat('     A formation can be simulated by specifying the formation parameters that determines the mean elements of a Deputy satellite.',...
       ' The relative orbit is propagated using a STM. In case of the Relative State STM, the relative states in the curvilinear frame',...
       ' are generated using a first-order transformation between the mean and true longitudes of Chief. For a given Chief orbit, the STM',...
       ' needs to be generated only once and different formations can be propagated using STM without propagating Chief''s orbit again.',...
       ' The ability to compute numerically propagated results using GMAT or MATLAB is provided only to compare the accuracy of the analytic',...
       ' theory. All the efforts were made to minimize the time it took to implement this proof-of-concept code, so its execution speed is not',...
       ' great. Despite that keep in mind, this is a fully analytic theory, so you can use any time step.'),...
       ' ', strcat('This code is part of the supplementary material to the publication as given below.',...
       ' It is a Jagerware. That is, it is free software to use, modify and/or redistribute as far as the activities are non-commercial. ',...
       ' In any case, you should cite the original author and the publication (as given below). If you find it useful, then you owe the ',...
       ' original author a Jagerbomb!'),' ',strcat('Mahajan, B., "Absolute and Relative Motion Satellite Theories for Zonal and ',...
       ' Tesseral Gravitational Harmonics," PhD Dissertation, Texas A&M University, 2018.')};
	
d = dialog('Name','SwARM Information','Units','normalized','Position',[0.25 0.25 0.55 0.61]);
annotation(d,'textbox','Units','normalized','Position',[0 .05 1 .9],'FitBoxToText','off',...
                    'HorizontalAlignment','center','Interpreter','tex','EdgeColor','none',...
                    'String',msg1);

annotation(d,'textbox','Units','normalized','Position',[0 0 1 .9],'FitBoxToText','on',...
                    'HorizontalAlignment','center','Interpreter','tex','EdgeColor','none',...
                    'String',msg2);                
annotation(d,'textbox','Units','normalized','Position',[0 -0.15 1 .9],'FitBoxToText','off',...
                    'HorizontalAlignment','left','Interpreter','tex','EdgeColor','none',...
                    'String',msg3);                                
uicontrol(d, 'Style', 'pushbutton', 'String', 'CLOSE','Units','normalized','Position',[0.45 0 .1 .05],...
                        'Callback', 'close');    
                    
                    
                    
% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function text_mu_Callback(hObject, eventdata, handles)
% hObject    handle to text_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_mu as text
%        str2double(get(hObject,'String')) returns contents of text_mu as a double


% --- Executes during object creation, after setting all properties.
function text_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_Re_Callback(hObject, eventdata, handles)
% hObject    handle to text_Re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_Re as text
%        str2double(get(hObject,'String')) returns contents of text_Re as a double


% --- Executes during object creation, after setting all properties.
function text_Re_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_Re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_orderIC_Callback(hObject, eventdata, handles)
% hObject    handle to text_orderIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_orderIC as text
%        str2double(get(hObject,'String')) returns contents of text_orderIC as a double


% --- Executes during object creation, after setting all properties.
function text_orderIC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_orderIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to text_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_degree as text
%        str2double(get(hObject,'String')) returns contents of text_degree as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_order_Callback(hObject, eventdata, handles)
% hObject    handle to text_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_order as text
%        str2double(get(hObject,'String')) returns contents of text_order as a double


% --- Executes during object creation, after setting all properties.
function text_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_degreeNum_Callback(hObject, eventdata, handles)
% hObject    handle to text_degreeNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_degreeNum as text
%        str2double(get(hObject,'String')) returns contents of text_degreeNum as a double


% --- Executes during object creation, after setting all properties.
function text_degreeNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_degreeNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_orderNum_Callback(hObject, eventdata, handles)
% hObject    handle to text_orderNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_orderNum as text
%        str2double(get(hObject,'String')) returns contents of text_orderNum as a double


% --- Executes during object creation, after setting all properties.
function text_orderNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_orderNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function progressbar_Callback(hObject, eventdata, handles)
% hObject    handle to progressbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function progressbar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to progressbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over progressbar.
function progressbar_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to progressbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in check_mex.
function check_mex_Callback(hObject, eventdata, handles)
% hObject    handle to check_mex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_mex



function text_resultfile_Callback(hObject, eventdata, handles)
% hObject    handle to text_resultfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_resultfile as text
%        str2double(get(hObject,'String')) returns contents of text_resultfile as a double


% --- Executes during object creation, after setting all properties.
function text_resultfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_resultfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SwARM Internal Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SwarmConf = ReadGUI(handles)

% Chief Initial Conditions
SwarmConf.ChiefIC.SMA = str2double(get(handles.text_sma,'String'));
SwarmConf.ChiefIC.ECC = str2double(get(handles.text_ecc,'String'));
SwarmConf.ChiefIC.INC = str2double(get(handles.text_inc,'String'));
SwarmConf.ChiefIC.RAAN = str2double(get(handles.text_raan,'String'));
SwarmConf.ChiefIC.AOP = str2double(get(handles.text_aop,'String'));
SwarmConf.ChiefIC.M = str2double(get(handles.text_M,'String'));

% Deputy Initial Conditions
SwarmConf.RM.RHOX = str2double(get(handles.text_rhox,'String'));
SwarmConf.RM.RHOY = str2double(get(handles.text_rhoy,'String'));
SwarmConf.RM.RHOZ = str2double(get(handles.text_rhoz,'String'));
SwarmConf.RM.ALPHAX = str2double(get(handles.text_alphax,'String'));
SwarmConf.RM.ALPHAZ = str2double(get(handles.text_alphaz,'String'));
SwarmConf.RM.DriftCond = str2double(get(handles.text_driftcond,'String'));

% Other Configuration
SwarmConf.GravCoeffFN = get(handles.text_coefffile,'String');
SwarmConf.J2MEXON = get(handles.check_mex,'Value');

SwarmConf.GravModel.IC.n = str2double(get(handles.text_degreeIC,'String'));
SwarmConf.GravModel.IC.m = str2double(get(handles.text_orderIC,'String'));
SwarmConf.GravModel.AST.n = str2double(get(handles.text_degree,'String'));
SwarmConf.GravModel.AST.m = str2double(get(handles.text_order,'String'));
SwarmConf.GravModel.NP.n = str2double(get(handles.text_degreeNum,'String'));
SwarmConf.GravModel.NP.m = str2double(get(handles.text_orderNum,'String'));

% make sure drift condition degree is not greater than theory IC degree
if SwarmConf.RM.DriftCond > SwarmConf.GravModel.IC.n
    SwarmConf.RM.DriftCond = SwarmConf.GravModel.IC.n;
    set(handles.text_driftcond,'String',num2str(SwarmConf.RM.DriftCond));
end

TessSel = get(handles.pop_theory_QuadTesserals,'Value');
if TessSel == 1
    SwarmConf.GravModel.AST.QuadTesseralsOn = false;
else
    SwarmConf.GravModel.AST.QuadTesseralsOn = true;
end

% choose propagator
prop_gmat = get(handles.radio_np_gmat,'Value');
prop_matlab = get(handles.radio_np_matlab,'Value');
if prop_gmat == true
    SwarmConf.NumPropType = 'GMAT';
elseif prop_matlab == true
    SwarmConf.NumPropType = 'MATLAB';
else
    SwarmConf.NumPropType = 'None';
end

SwarmConf.ResultsFN = get(handles.text_resultfile,'String');
SwarmConf.GMAT.ScriptFN = get(handles.text_scriptfile,'String');
SwarmConf.GMAT.ExePath = get(handles.text_gmatpath,'String');

% choose STM type
stmtypede = get(handles.radio_stm_de,'Value');
stmtyperel = get(handles.radio_stm_relstate,'Value');
if stmtypede == true
    SwarmConf.STMType = 'Differential Elements';
elseif stmtyperel == true
    SwarmConf.STMType = 'Relative States';
else
    SwarmConf.STMType = 'None';
end

% Time and Tolerance
SwarmConf.tol = str2double(get(handles.text_tol,'String'));
SwarmConf.quadtol = str2double(get(handles.text_quadtol,'String'));
SwarmConf.GravModel.mu = str2double(get(handles.text_mu,'String'));
SwarmConf.GravModel.Re = str2double(get(handles.text_Re,'String'));
SwarmConf.Time.MJDUTCt0 = str2double(get(handles.text_MJDUTCt0,'String'));
SwarmConf.Time.DataStartTime = str2double(get(handles.text_datat0,'String'));
SwarmConf.Time.NumOrbits = str2double(get(handles.text_numorbits,'String'));
SwarmConf.Time.PointsPerOrbit = str2double(get(handles.text_pointsperorbit,'String'));


function handles = LoadGUI(handles)

% internal Swarm run data
SwarmConf = handles.SwarmConf;

% Chief Initial Conditions
set(handles.text_sma,'String',num2str(SwarmConf.ChiefIC.SMA, 16));
set(handles.text_ecc,'String',num2str(SwarmConf.ChiefIC.ECC,16));
set(handles.text_inc,'String',num2str(SwarmConf.ChiefIC.INC,16));
set(handles.text_raan,'String',num2str(SwarmConf.ChiefIC.RAAN,16));
set(handles.text_aop,'String',num2str(SwarmConf.ChiefIC.AOP,16));
set(handles.text_M,'String',num2str(SwarmConf.ChiefIC.M,16));

% Deputy Initial Conditions
set(handles.text_rhox,'String',num2str(SwarmConf.RM.RHOX,16));
set(handles.text_rhoy,'String',num2str(SwarmConf.RM.RHOY,16));
set(handles.text_rhoz,'String',num2str(SwarmConf.RM.RHOZ,16));
set(handles.text_alphax,'String',num2str(SwarmConf.RM.ALPHAX,16));
set(handles.text_alphaz,'String',num2str(SwarmConf.RM.ALPHAZ,16));
set(handles.text_driftcond,'String',num2str(SwarmConf.RM.DriftCond,8));

% Other Configuration
set(handles.text_coefffile,'String',SwarmConf.GravCoeffFN);
set(handles.check_mex,'Value',SwarmConf.J2MEXON);

% Gravity Model
set(handles.text_degreeIC,'String',SwarmConf.GravModel.IC.n);
set(handles.text_orderIC,'String',SwarmConf.GravModel.IC.m);
set(handles.text_degree,'String',SwarmConf.GravModel.AST.n);
set(handles.text_order,'String',SwarmConf.GravModel.AST.m);
set(handles.text_degreeNum,'String',SwarmConf.GravModel.NP.n);
set(handles.text_orderNum,'String',SwarmConf.GravModel.NP.m);

if SwarmConf.GravModel.AST.QuadTesseralsOn == true
    set(handles.pop_theory_QuadTesserals,'Value',2);
else
    set(handles.pop_theory_QuadTesserals,'Value',1);
end

% set propagator
if strcmpi(SwarmConf.NumPropType, 'GMAT') == 1
    set(handles.radio_np_gmat,'Value',1);
elseif strcmpi(SwarmConf.NumPropType, 'MATLAB') == 1
    set(handles.radio_np_matlab,'Value',1);
else
    set(handles.radio_np_none,'Value',1);
end

set(handles.text_resultfile,'String',SwarmConf.ResultsFN);
set(handles.text_scriptfile,'String',SwarmConf.GMAT.ScriptFN);
set(handles.text_gmatpath,'String',SwarmConf.GMAT.ExePath);

% set STM type
if strcmpi(SwarmConf.STMType, 'Differential Elements') == 1
    set(handles.radio_stm_de,'Value',1);
elseif strcmpi(SwarmConf.STMType, 'Relative States') == 1
    set(handles.radio_stm_relstate,'Value',1);
else
    set(handles.radio_stm_none,'Value',1);
end

% set Time and Tolerance
set(handles.text_tol,'String',num2str(SwarmConf.tol,16));
set(handles.text_quadtol,'String',num2str(SwarmConf.quadtol,16));
set(handles.text_mu,'String',num2str(SwarmConf.GravModel.mu,16));
set(handles.text_Re,'String',num2str(SwarmConf.GravModel.Re,16));
set(handles.text_MJDUTCt0,'String',num2str(SwarmConf.Time.MJDUTCt0,16));
set(handles.text_datat0,'String',num2str(SwarmConf.Time.DataStartTime,16));
set(handles.text_numorbits,'String',num2str(SwarmConf.Time.NumOrbits,16));
set(handles.text_pointsperorbit,'String',num2str(SwarmConf.Time.PointsPerOrbit,16));



function handles = AbsPlotData(absplot, handles)

% if Data start time is not at t0, then skip first data point
if handles.SwarmConf.Time.DataStartTime ~= 0
    dIn = 2;
    T = handles.SwarmRun.tspan(2:end)*handles.SwarmRun.tunit/86400;
else
    dIn = 1;
    T = handles.SwarmRun.tspan*handles.SwarmRun.tunit/86400;
end

AnalCOE = handles.SwarmRun.AnalCOE;
AnalCOE(:,1) = AnalCOE(:,1)*handles.SwarmRun.lunit;
AnalCOE(:,3:6) = rad2deg(AnalCOE(:,3:6));

AnalMCOE = handles.SwarmRun.AnalMCOE;
AnalMCOE(:,1) = AnalMCOE(:,1)*handles.SwarmRun.lunit;
AnalMCOE(:,3:6) = rad2deg(AnalMCOE(:,3:6));

AnalCRV = handles.SwarmRun.AnalCRV;
AnalCRV(:,1:3) = AnalCRV(:,1:3)*handles.SwarmRun.lunit;
AnalCRV(:,4:6) = AnalCRV(:,4:6)*handles.SwarmRun.vunit;

% if numerical propagation is enabled
if strcmpi(handles.SwarmConf.NumPropType,'None') ~= 1

    NumCOE = handles.SwarmRun.NumCOE;
    NumCOE(:,1) = NumCOE(:,1)*handles.SwarmRun.lunit;
    NumCOE(:,3:6) = rad2deg(NumCOE(:,3:6));
    
    NumCRV = handles.SwarmRun.NumCRV;
    NumCRV(:,1:3) = NumCRV(:,1:3)*handles.SwarmRun.lunit;
    NumCRV(:,4:6) = NumCRV(:,4:6)*handles.SwarmRun.vunit;
    
end

% default plot aspect ratio
handles.AbsFig.DataAspectRatioMode = 'auto';

switch(absplot)


    case {1,2,3,4,5,6}
        % Chief Classical elements
        clabels = {'a','e','i','\Omega','\omega','MA'};
        data = AnalCOE(dIn:end,absplot);
        data2 = AnalMCOE(dIn:end,absplot);
        plot(handles.AbsFig, T,data,T,data2,'.','LineWidth',2);
        
        if strcmpi(handles.SwarmConf.NumPropType,'None') ~= 1
            hold(handles.AbsFig, 'on');
            data = NumCOE(dIn:end,absplot);
            plot(handles.AbsFig, T,data,'--','LineWidth',2);
            legend(handles.AbsFig,'Theory','Mean','Num');
            hold(handles.AbsFig, 'off');
        end
        xlabel(handles.AbsFig,'Time [days]');
        ylabel(handles.AbsFig, clabels{absplot});
        
    case 7
        data = AnalCRV(dIn:end,1:3);
        plot3(handles.AbsFig, data(:,1),data(:,2),data(:,3),'LineWidth',2);
        if strcmpi(handles.SwarmConf.NumPropType,'None') ~= 1
            hold(handles.AbsFig, 'on');
            data = NumCRV(dIn:end,1:3);
            plot3(handles.AbsFig, data(:,1),data(:,2),data(:,3),'--','LineWidth',2);
            legend(handles.AbsFig,'Theory','Num');
            hold(handles.AbsFig, 'off');
        end
        xlabel(handles.AbsFig,'x');
        ylabel(handles.AbsFig,'y');
        zlabel(handles.AbsFig,'z');

        % axis equal mode
        handles.AbsFig.DataAspectRatioMode = 'manual';
        handles.AbsFig.DataAspectRatio = [1 1 1];
        handles.AbsFig.PlotBoxAspectRatio = [3 4 4];
        
    case 8
        if strcmpi(handles.SwarmConf.NumPropType,'None') ~= 1
            
            % Chief RSS Position Error
            data = NumCRV(dIn:end,1:3)-AnalCRV(dIn:end,1:3);
            data = sqrt(sum(data.^2,2));
            plot(handles.AbsFig,T,data,'LineWidth',2);
            xlabel(handles.AbsFig,'Time [days]');
            ylabel(handles.AbsFig, 'RSS Pos Err');
        else
            plot(handles.AbsFig,0,0);
            ylabel(handles.AbsFig, 'No Numerical Prop. Data');
        end
        
    case 9
        if strcmpi(handles.SwarmConf.NumPropType,'None') ~= 1
            data = NumCRV(dIn:end,4:6)-AnalCRV(dIn:end,4:6);
            data = sqrt(sum(data.^2,2));
            plot(handles.AbsFig,T,data,'LineWidth',2);
            xlabel(handles.AbsFig,'Time [days]');
            ylabel(handles.AbsFig, 'RSS Vel Err');
        else
            plot(handles.AbsFig,0,0);
            ylabel(handles.AbsFig, 'No Numerical Prop. Data');
        end
end

% Common plot settings
grid(handles.AbsFig,'on');


        
function handles = RelPlotData(relplot, handles)

% if Data start time is not at t0, then skip first data point
if handles.SwarmConf.Time.DataStartTime ~= 0
    dIn = 2;
    T = handles.SwarmRun.tspan(2:end)*handles.SwarmRun.tunit/86400;
else
    dIn = 1;
    T = handles.SwarmRun.tspan*handles.SwarmRun.tunit/86400;
end

AnalDiffCOE = handles.SwarmRun.AnalDiffCOE;
AnalDiffCOE(:,1) = AnalDiffCOE(:,1)*handles.SwarmRun.lunit;
AnalDiffCOE(:,3:6) = rad2deg(AnalDiffCOE(:,3:6));

AnalCurvRV = handles.SwarmRun.AnalCurvRV;
AnalCurvRV(:,1:3) = AnalCurvRV(:,1:3)*handles.SwarmRun.lunit;
AnalCurvRV(:,4:6) = AnalCurvRV(:,4:6)*handles.SwarmRun.vunit;

PlotErrs = true;
if strcmpi(handles.SwarmConf.NumPropType,'None') == 1 || handles.SwarmRun.DepModified == true
    % No num. data or its not valid
    PlotErrs = false;
end
    
% if numerical propagation is enabled
if PlotErrs == true

    NumCOE = handles.SwarmRun.NumCOE;
    NumCOE(:,1) = NumCOE(:,1)*handles.SwarmRun.lunit;
    NumCOE(:,3:6) = rad2deg(NumCOE(:,3:6));

    NumDOE = handles.SwarmRun.NumDOE;
    NumDOE(:,1) = NumDOE(:,1)*handles.SwarmRun.lunit;
    NumDOE(:,3:6) = rad2deg(NumDOE(:,3:6));

    NumDiffCOE = NumDOE - NumCOE;    
    
    NumCurvRV = handles.SwarmRun.NumCurvRV;
    NumCurvRV(:,1:3) = NumCurvRV(:,1:3)*handles.SwarmRun.lunit;
    NumCurvRV(:,4:6) = NumCurvRV(:,4:6)*handles.SwarmRun.vunit;

end

% default plot aspect ratio
handles.RelFig.DataAspectRatioMode = 'auto';

switch(relplot)
    
    case 1

        data = AnalCurvRV(dIn:end,1:3);
        plot3(handles.RelFig,data(:,1),data(:,2),data(:,3),'LineWidth',2);
        text(handles.RelFig,0,0,0,'Chief');
        if PlotErrs == true
            hold(handles.RelFig, 'on');
            data = NumCurvRV(dIn:end,1:3);
            plot3(handles.RelFig,data(:,1),data(:,2),data(:,3),'--','LineWidth',2);
            legend(handles.RelFig,'STM','Num.');
            hold(handles.RelFig, 'off');
        end
        xlabel(handles.RelFig,'radial');
        ylabel(handles.RelFig,'along-track');
        zlabel(handles.RelFig,'cross-track');
        % axis equal model
        handles.RelFig.DataAspectRatioMode = 'manual';
        handles.RelFig.DataAspectRatio = [1 1 1];
        handles.RelFig.PlotBoxAspectRatio = [3 4 4];

        
    case {2,3,4,5,6,7}
        clabels = {'da','de','di','d\Omega','d\omega','dMA'};
        data = AnalDiffCOE(dIn:end,relplot-1);
        plot(handles.RelFig,T,data,'LineWidth',2);
        
        if PlotErrs == true
            hold(handles.RelFig, 'on');
            data = NumDiffCOE(dIn:end,relplot-1);
            plot(handles.RelFig,T,data,'--','LineWidth',2);
            legend(handles.RelFig,'STM','Num.');
            hold(handles.RelFig, 'off');
        end
        xlabel(handles.RelFig,'Time [days]');
        ylabel(handles.RelFig, clabels{relplot-1});
                
    case 8
            data = AnalCurvRV(dIn:end,1);
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            if PlotErrs == true
                hold(handles.RelFig, 'on');
                data = NumCurvRV(dIn:end,1);
                plot(handles.RelFig,T, data(:,1),'--','LineWidth',2);
                legend(handles.RelFig,'STM','Num.');
                hold(handles.RelFig, 'off');
            end
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'radial');
        
    case 9
            data = AnalCurvRV(dIn:end,2);
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            if PlotErrs == true
                hold(handles.RelFig, 'on');
                data = NumCurvRV(dIn:end,2);
                plot(handles.RelFig,T, data(:,1),'--','LineWidth',2);
                legend(handles.RelFig,'STM','Num.');
                hold(handles.RelFig, 'off');
            end
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'along-track');
        
    case 10
            data = AnalCurvRV(dIn:end,3);
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            if PlotErrs == true
                hold(handles.RelFig, 'on');
                data = NumCurvRV(dIn:end,3);
                plot(handles.RelFig,T, data(:,1),'--','LineWidth',2);
                legend(handles.RelFig,'STM','Num.');
                hold(handles.RelFig, 'off');
            end            
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'cross-track');

    case 11
        if PlotErrs == true
            data = NumCurvRV(dIn:end,1) - AnalCurvRV(dIn:end,1);
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'radial error');
        else
            plot(handles.RelFig,0,0);
            ylabel(handles.RelFig,'No Valid Numerical Prop. Data');
        end
        
    case 12
        if PlotErrs == true
            data = NumCurvRV(dIn:end,2) - AnalCurvRV(dIn:end,2);
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'along-track error');
        else
            plot(handles.RelFig,0,0);
            ylabel(handles.RelFig,'No Valid Numerical Prop. Data');
        end
        
    case 13
        if PlotErrs == true
            data = NumCurvRV(dIn:end,3) - AnalCurvRV(dIn:end,3);
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'cross-track error');
        else
            plot(handles.RelFig,0,0);
            ylabel(handles.RelFig,'No Valid Numerical Prop. Data');
        end
        
    case 14
        if PlotErrs == true
            data = NumCurvRV(dIn:end,1:3) - AnalCurvRV(dIn:end,1:3);
            data = sqrt(sum(data.^2,2));
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'rel pos error');
        else
            plot(handles.RelFig,0,0);
            ylabel(handles.RelFig,'No Valid Numerical Prop. Data');
        end
        
    case 15
        if PlotErrs == true
            data = NumCurvRV(dIn:end,4:6) - AnalCurvRV(dIn:end,4:6);
            data = sqrt(sum(data.^2,2));
            plot(handles.RelFig,T, data(:,1),'LineWidth',2);
            xlabel(handles.RelFig,'Time [days]');
            ylabel(handles.RelFig, 'rel vel error');
        else
            plot(handles.RelFig,0,0);
            ylabel(handles.RelFig,'No Valid Numerical Prop. Data');
        end
end

% Common plot settings
grid(handles.RelFig,'on');



% Computy differential equinoctial elements for Deputy
function dEqn = FormationdEqn(rhox, rhoy, rhoz, alphax, alphaz, ChiefEq, DriftCond, mu, Re, Jcoeff, tol)

% differential Equinoctial elements for GCO
dEqn = GCODiffEquinoctial(rhox,alphax,rhoy,rhoz,alphaz,ChiefEq);

% mean drift mitigation for differential SMA for osculating elements
if DriftCond >= 2
    [ daJ2, daJn ] = ATSecDriftCond( ChiefEq, dEqn(3),dEqn(4),dEqn(5),dEqn(6),mu,Re,...
                    DriftCond, Jcoeff, tol,  true, true );
    dEqn(1) = dEqn(1) + daJ2 + daJn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SwARM Internal Functions End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in button_UpdateFF.
function button_UpdateFF_Callback(hObject, eventdata, handles)
% hObject    handle to button_UpdateFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Indicate that Deputy or formation ICs have been modified.
handles.SwarmRun.DepModified = true;

% read the current formation design parameters from GUI
handles.SwarmConf.RM.RHOX = str2double(get(handles.text_rhox,'String'));
handles.SwarmConf.RM.RHOY = str2double(get(handles.text_rhoy,'String'));
handles.SwarmConf.RM.RHOZ = str2double(get(handles.text_rhoz,'String'));
handles.SwarmConf.RM.ALPHAX = str2double(get(handles.text_alphax,'String'));
handles.SwarmConf.RM.ALPHAZ = str2double(get(handles.text_alphaz,'String'));
handles.SwarmConf.RM.DriftCond = str2double(get(handles.text_driftcond,'String'));

% make sure drift condition degree is not greater than theory IC degree
if handles.SwarmConf.RM.DriftCond > handles.SwarmConf.GravModel.IC.n
    handles.SwarmConf.RM.DriftCond = handles.SwarmConf.GravModel.IC.n;
    set(handles.text_driftcond,'String',num2str(handles.SwarmConf.RM.DriftCond));
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% create two copies of the current two figure axis
f1 = figure;
copyobj(handles.AbsFig, f1);

f2 = figure;
copyobj(handles.RelFig, f2);


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_theory_QuadTesserals.
function pop_theory_QuadTesserals_Callback(hObject, eventdata, handles)
% hObject    handle to pop_theory_QuadTesserals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_theory_QuadTesserals contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_theory_QuadTesserals


% --- Executes during object creation, after setting all properties.
function pop_theory_QuadTesserals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_theory_QuadTesserals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipushtool4_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f1 = figure;
copyobj(handles.AbsFig, f1);
ax = findobj(f1, 'type','axes');
ax.OuterPosition = [0 0 1 1];
legend;

if strcmpi(handles.SwarmConf.STMType, 'None') ~= 1

f2 = figure;
copyobj(handles.RelFig, f2);
ax = findobj(f2, 'type','axes');
ax.OuterPosition = [0 0 1 1];
legend;

end


% --- Executes on button press in push_logo.
function push_logo_Callback(hObject, eventdata, handles)
% hObject    handle to push_logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
