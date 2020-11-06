function varargout = Practice_4_Diaz_Velazquez_Bustos(varargin)
% PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS MATLAB code for Practice_4_Diaz_Velazquez_Bustos.fig
%      PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS, by itself, creates a new PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS or raises the existing
%      singleton*.
%
%      H = PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS returns the handle to a new PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS or the handle to
%      the existing singleton*.
%
%      PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS.M with the given input arguments.
%
%      PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS('Property','Value',...) creates a new PRACTICE_4_DIAZ_VELAZQUEZ_BUSTOS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Practice_4_Diaz_Velazquez_Bustos_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Practice_4_Diaz_Velazquez_Bustos_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Practice_4_Diaz_Velazquez_Bustos

% Last Modified by GUIDE v2.5 04-Nov-2018 21:44:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Practice_4_Diaz_Velazquez_Bustos_OpeningFcn, ...
                   'gui_OutputFcn',  @Practice_4_Diaz_Velazquez_Bustos_OutputFcn, ...
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


% --- Executes just before Practice_4_Diaz_Velazquez_Bustos is made visible.
function Practice_4_Diaz_Velazquez_Bustos_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Practice_4_Diaz_Velazquez_Bustos (see VARARGIN)

% Choose default command line output for Practice_4_Diaz_Velazquez_Bustos
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Practice_4_Diaz_Velazquez_Bustos wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Practice_4_Diaz_Velazquez_Bustos_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
W = str2double(get(handles.edit1,'String'));
S = str2double(get(handles.edit2,'String'));
CD0 = str2double(get(handles.edit3,'String'));
p = str2double(get(handles.edit4,'String'));
k = str2double(get(handles.edit5,'String'));
T = str2double(get(handles.edit10,'String'));
g=9.81;


vr=sqrt((2*W)/(p*S))*(k/CD0)^(1/4);
Em=(1/(2*sqrt(CD0*k)));
z= (T*Em)/W;
%%%% For MSTR%%%%
uMSTR=1;
%vMSTR=uMSTR*vr;
wMSTR=(g/vr)*sqrt(2*(z-1));
RMSTR=(vr^2*uMSTR^2)/(g*sqrt(2*(z-1)));
NMSTR=sqrt(2*z-1);
clMSTR=sqrt(CD0/k)*sqrt(2*z-1);
EMSTR=Em*(sqrt(2*z-1)/z);



%%% For STT %%%
uSTT=1/sqrt(z);
wSTT=(g/vr)*(sqrt((z^2 - 1)/z));
RSTT=(vr^2/g)*(1/sqrt(z^2 - 1));
NSTT=(1/z)*sqrt(2*z^2 - 1);
clSTT=sqrt(CD0/k)*(sqrt(2*z^2 - 1));
ESTT=Em*sqrt(2*z^2 - 1)/z^2;

%%% For Nmax %%%
uNM=sqrt(z);
wNM=(g/vr)*sqrt((z^2 - 1)/z);
RNM=(vr^2/g)*(z/sqrt(z^2 - 1));
NNM=z;
clNM=sqrt(CD0/k);
ENM=Em;
data=[uMSTR uSTT uNM; wMSTR wSTT wNM; RMSTR RSTT RNM; NMSTR NSTT NNM; clMSTR clSTT clNM; EMSTR ESTT ENM];
set(handles.uitable1,'data',data);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
W = str2double(get(handles.edit1,'String'));
S = str2double(get(handles.edit2,'String'));
CD0 = str2double(get(handles.edit3,'String'));
p = str2double(get(handles.edit4,'String'));
k = str2double(get(handles.edit5,'String'));
power = str2double(get(handles.edit8,'String'));
E = str2double(get(handles.edit9,'String'));
g=9.81;

vr=sqrt((2*W)/(p*S))*(k/CD0)^(1/4);
Em=(1/(2*sqrt(CD0*k)));
P=(E*power*1000*Em)/(vr*W);

%%%For STT%%%
uSTT=2/(3*P);
NSTT=sqrt(2*P*uSTT-uSTT^4);
clSTT=sqrt(CD0/k)*(NSTT/(uSTT)^2);
RSTT=(vr^2*uSTT^2)/(g*sqrt(NSTT^2-1));
wSTT=(g*sqrt(NSTT^2-1))/(vr*uSTT);
ESTT=clSTT/(CD0+k*(clSTT)^2);

%%For NMAX%%%
uNM=(P/2)^(1/3);
NNM=sqrt(2*P*uNM-uNM^4);
clNM=sqrt(CD0/k)*(NNM/(uNM)^2);
RNM=(vr^2*uNM^2)/(g*sqrt(NNM^2-1));
wNM=(g*sqrt(NNM^2-1))/(vr*uNM);
ENM=clNM/(CD0+k*(clNM)^2);

%%%For MSTR%%&
rts = roots([1 0 0 P -1]);
    urts = rts(imag(rts)==0);
    uMSTR=urts(urts>=0);
NMSTR=sqrt(2*P*uMSTR-uMSTR^4);
clMSTR=sqrt(CD0/k)*(NMSTR/(uMSTR)^2);
RMSTR=(vr^2*uMSTR^2)/(g*sqrt(NMSTR^2-1));
wMSTR=(g*sqrt(NMSTR^2-1))/(vr*uMSTR);
EMSTR=clMSTR/(CD0+k*(clMSTR)^2);
data=[uMSTR uSTT uNM; wMSTR wSTT wNM; RMSTR RSTT RNM; NMSTR NSTT NNM; clMSTR clSTT clNM; EMSTR ESTT ENM];
set(handles.uitable1,'data',data);




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, ~, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, ~)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, ~, handles)
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


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(~, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit1,'String','');
set(handles.edit2,'String','');
set(handles.edit3,'String','');
set(handles.edit4,'String','');
set(handles.edit5,'String','');
set(handles.edit6,'String','');
set(handles.edit7,'String','');
set(handles.edit8,'String','');
set(handles.edit9,'String','');
set(handles.edit10,'String','');




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


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
W = str2double(get(handles.edit1,'String'));
S = str2double(get(handles.edit2,'String'));
CD0 = str2double(get(handles.edit3,'String'));
p = str2double(get(handles.edit4,'String'));
k = str2double(get(handles.edit5,'String'));
clmax= str2double(get(handles.edit6,'String'));
Nlim= str2double(get(handles.edit7,'String'));
power = str2double(get(handles.edit8,'String'));
E = str2double(get(handles.edit9,'String'));
g=9.81;

vr=sqrt((2*W)/(p*S))*(k/CD0)^(1/4);
Em=(1/(2*sqrt(CD0*k)));
P=(E*power*1000*Em)/(vr*W);

%%%%RESTRICTIONS%%%


%%%For STT%%%
uSTT=2/(3*P);
NSTT=sqrt(2*P*uSTT-uSTT^4);

%%For NMAX%%%
uNM=(P/2)^(1/3);
NNM=sqrt(2*P*uNM-uNM^4);

%%%For MSTR%%&
rts = roots([1 0 0 P -1]);
    urts = rts(imag(rts)==0);
    uMSTR=urts(urts>=0);
NMSTR=sqrt(2*P*uMSTR-uMSTR^4);

if NSTT> Nlim
    NSTT=Nlim;
end

if NMSTR> Nlim
    NMSTR=Nlim;
end

if NNM> Nlim
    NNM=Nlim;
end

%%%cl%%% 
clSTT=sqrt(CD0/k)*(NSTT/(uSTT)^2);
clNM=sqrt(CD0/k)*(NNM/(uNM)^2);
clMSTR=sqrt(CD0/k)*(NMSTR/(uMSTR)^2);


     %%SST%%
RSTT=(vr^2*uSTT^2)/(g*sqrt(NSTT^2-1));
wSTT=(g*sqrt(NSTT^2-1))/(vr*uSTT);
ESTT=clSTT/(CD0+k*(clSTT)^2);


  %%nm%%
RNM=(vr^2*uNM^2)/(g*sqrt(NNM^2-1));
wNM=(g*sqrt(NNM^2-1))/(vr*uNM);
ENM=clNM/(CD0+k*(clNM)^2);


%%mstr%%
RMSTR=(vr^2*uMSTR^2)/(g*sqrt(NMSTR^2-1));
wMSTR=(g*sqrt(NMSTR^2-1))/(vr*uMSTR);
EMSTR=clMSTR/(CD0+k*(clMSTR)^2);

%%%Restriction for STT%%
if clSTT>clmax
    uSTT=((2*P)/(((k*clmax^2)/CD0)+1))^(1/3);
    NSTT=sqrt(2*P*uSTT-uSTT^4);
    if NSTT> Nlim
    NSTT=Nlim;
    end
    clSTT=clmax;
     RSTT=(vr^2*uSTT^2)/(g*sqrt(NSTT^2-1));
wSTT=(g*sqrt(NSTT^2-1))/(vr*uSTT);
ESTT=clSTT/(CD0+k*(clSTT)^2);
end

if clMSTR>clmax
    uMSTR=((2*P)/(((k*clmax^2)/CD0)+1))^(1/3);
    NMSTR=sqrt(2*P*uMSTR-uMSTR^4);
    if NMSTR> Nlim
    NMSTR=Nlim;
    end
    clMSTR=clmax;
     RMSTR=(vr^2*uMSTR^2)/(g*sqrt(NMSTR^2-1));
wMSTR=(g*sqrt(NMSTR^2-1))/(vr*uMSTR);
EMSTR=clMSTR/(CD0+k*(clMSTR)^2);
end

if clNM>clmax
    uNM=((2*P)/(((k*clmax^2)/CD0)+1))^(1/3);
    NNM=sqrt(2*P*uNM-uNM^4);
    if NNM> Nlim
    NNM=Nlim;
    end
    clNM=clmax;
     RNM=(vr^2*uNM^2)/(g*sqrt(NNM^2-1));
wNM=(g*sqrt(NNM^2-1))/(vr*uNM);
ENM=clNM/(CD0+k*(clNM)^2);
end

data=[uMSTR uSTT uNM; wMSTR wSTT wNM; RMSTR RSTT RNM; NMSTR NSTT NNM; clMSTR clSTT clNM; EMSTR ESTT ENM];
set(handles.uitable1,'data',data);


% --- Executes during object creation, after setting all properties.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
W = str2double(get(handles.edit1,'String'));
S = str2double(get(handles.edit2,'String'));
CD0 = str2double(get(handles.edit3,'String'));
p = str2double(get(handles.edit4,'String'));
k = str2double(get(handles.edit5,'String'));
clmax= str2double(get(handles.edit6,'String'));
Nlim= str2double(get(handles.edit7,'String'));
T = str2double(get(handles.edit10,'String'));
g=9.81;


vr=sqrt((2*W)/(p*S))*(k/CD0)^(1/4);
Em=(1/(2*sqrt(CD0*k)));
z= (T*Em)/W;

%%% NLIM %%%%%
NMSTR=sqrt(2*z-1);
NSTT=(1/z)*sqrt(2*z^2 - 1);
NNM=z;

if NSTT> Nlim
    NSTT=Nlim;
end

if NMSTR> Nlim
    NMSTR=Nlim;
end

if NNM> Nlim
    NNM=Nlim;
end


%%CLmax%%%

clMSTR=sqrt(CD0/k)*sqrt(2*z-1);
clSTT=sqrt(CD0/k)*(sqrt(2*z^2 - 1));
clNM=sqrt(CD0/k);

if clSTT> clmax
    clSTT=clmax;
end

if clMSTR> clmax
    clMSTR=clmax;
end

if clNM> clmax
    clNM=clmax;
end
%%%% For MSTR%%%%
uMSTR=1;
%vMSTR=uMSTR*vr;
wMSTR=(g/vr)*sqrt(2*(z-1));
RMSTR=(vr^2*uMSTR^2)/(g*sqrt(2*(z-1)));

EMSTR=Em*(sqrt(2*z-1)/z);



%%% For STT %%%
uSTT=1/sqrt(z);
wSTT=(g/vr)*(sqrt((z^2 - 1)/z));
RSTT=(vr^2/g)*(1/sqrt(z^2 - 1));


ESTT=Em*sqrt(2*z^2 - 1)/z^2;

%%% For Nmax %%%
uNM=sqrt(z);
wNM=(g/vr)*sqrt((z^2 - 1)/z);
RNM=(vr^2/g)*(z/sqrt(z^2 - 1));


ENM=Em;




data=[uMSTR uSTT uNM; wMSTR wSTT wNM; RMSTR RSTT RNM; NMSTR NSTT NNM; clMSTR clSTT clNM; EMSTR ESTT ENM];
set(handles.uitable1,'data',data);



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
W = str2double(get(handles.edit1,'String'));
S = str2double(get(handles.edit2,'String'));
CD0 = str2double(get(handles.edit3,'String'));
p = str2double(get(handles.edit4,'String'));
k = str2double(get(handles.edit5,'String'));
clmax= str2double(get(handles.edit6,'String'));
Nlim= str2double(get(handles.edit7,'String'));
T = str2double(get(handles.edit10,'String'));
g=9.81;


vr=sqrt((2*W)/(p*S))*(k/CD0)^(1/4);
Em=(1/(2*sqrt(CD0*k)));
z= (T*Em)/W;

%%% NLIM %%%%%
NMSTR=sqrt(2*z-1);
NSTT=(1/z)*sqrt(2*z^2 - 1);
NNM=z;

if NSTT> Nlim
    NSTT=Nlim;
end

if NMSTR> Nlim
    NMSTR=Nlim;
end

if NNM> Nlim
    NNM=Nlim;
end


%%CLmax%%%

clMSTR=sqrt(CD0/k)*sqrt(2*z-1);
clSTT=sqrt(CD0/k)*(sqrt(2*z^2 - 1));
clNM=sqrt(CD0/k);

if clSTT> clmax
    clSTT=clmax;
end

if clMSTR> clmax
    clMSTR=clmax;
end

if clNM> clmax
    clNM=clmax;
end
%%%% For MSTR%%%%
uMSTR=1;
%vMSTR=uMSTR*vr;
wMSTR=(g/vr)*sqrt(2*(z-1));
RMSTR=(vr^2*uMSTR^2)/(g*sqrt(2*(z-1)));

EMSTR=Em*(sqrt(2*z-1)/z);



%%% For STT %%%
uSTT=1/sqrt(z);
wSTT=(g/vr)*(sqrt((z^2 - 1)/z));
RSTT=(vr^2/g)*(1/sqrt(z^2 - 1));


ESTT=Em*sqrt(2*z^2 - 1)/z^2;

%%% For Nmax %%%
uNM=sqrt(z);
wNM=(g/vr)*sqrt((z^2 - 1)/z);
RNM=(vr^2/g)*(z/sqrt(z^2 - 1));


ENM=Em;




data=[uMSTR uSTT uNM; wMSTR wSTT wNM; RMSTR RSTT RNM; NMSTR NSTT NNM; clMSTR clSTT clNM; EMSTR ESTT ENM];
set(handles.uitable1,'data',data);
