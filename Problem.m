% Copyright 2011 2012 M.Sc. Porfirio Armando Rodríguez

% This file is part of LPApp
% LPApp is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% LPApp is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function varargout = Problem(varargin)
% PROBLEM M-file for problem.fig
%      PROBLEM, by itself, creates a new PROBLEM or raises the existing
%      singleton*.
%
%      H = PROBLEM returns the handle to a new PROBLEM or the handle to
%      the existing singleton*.
%
%      PROBLEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROBLEM.M with the given input arguments.
%
%      PROBLEM('Property','Value',...) creates a new PROBLEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Problem_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Problem_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Problem_OpeningFcn, ...
                   'gui_OutputFcn',  @Problem_OutputFcn, ...
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


% --- Executes just before problem is made visible.
function Problem_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to problem (see VARARGIN)

% Choose default command line output for problem
handles.output = hObject;

movegui(hObject,'center');
% se busca el nombre de la aplicación en la línea de llamada
LPApphandle = find(strcmp(varargin, 'LPApp'));
if ~isempty(LPApphandle)
   handles.LPApphandle = varargin{LPApphandle+1};
else
    return;
end
% Update handles structure
guidata(hObject, handles);
% se obtienen las dimensiones del problema
dimp = size(handles.LPApphandle.gui_Matrix_problem);
dim = [dimp(1)-1, dimp(2)-1];
% se rotulan las filas y columnas de las tablas
colName = cell(dim(2)+1, 1);
for i = 1:dim(2)
    colName(i) = cellstr(strcat('X',num2str(i)));
end
colName(dim(2)+1) = cellstr('Yi0');
set(handles.table_problem, 'columnname', colName);

rowName=cell(1,dim(1)+1);
for i = 1:dim(1)
    rowName(1,i) = cellstr(strcat('f',num2str(i)));
end
rowName(dim(1)+1) = cellstr('Z');
set(handles.table_problem, 'rowname', rowName);
% se establece el formato de las columnas
colFormat=cell(1,dim(2)+1);
for i = 1:(dim(2)+1)
    colFormat(1,i) = cellstr('rat');
end
set(handles.table_problem, 'columnformat', colFormat);
% se configuran las columnas como editables
colEdit = ones(1,dim(2)+1);
set(handles.table_problem, 'columneditable', (colEdit == 1));
% se despliega la tabla para el ingreso de datos
All_display = handles.LPApphandle.gui_Matrix_problem;
set(handles.table_problem, 'data', All_display);

% UIWAIT makes problem wait for user response (see UIRESUME)
% uiwait(handles.Problem);


% --- Outputs from this function are returned to the command line.
function varargout = Problem_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem Order_initial;
% Matrix_problem    guarda la especificación del problema
% Orden_initial     registra el orden de las columnas según la matriz
% identidad

% recupera la especificación dada por el usuario
All_display = get(handles.table_problem, 'data');
Matrix_problem = All_display;
% se verifica que la especificación este en forma estándar y canónica
if (~isstandardform())
    return;
end
handles.LPApphandle.Order_initial = Order_initial; % se comparte el orden inicial
handles.LPApphandle.setProblem(Matrix_problem, handles.LPApphandle); % se establece el problema actual
delete(handles.Problem);

% ---Funcion utilitaria: Verifica si se ha ingresado la formulación del
% problema de manera correcta
function iscorrect = isstandardform()
global Matrix_problem Order_initial;

dim = size(Matrix_problem);
Order_initial = zeros(1, dim(2));
% se construye la matriz identidad
identidad = eye(dim(1)-1);
% se verifica que la matriz identidad este inmersa en la especifición del
% problema
for i =1:(dim(1)-1)
    iscorrect = 0;
    for j = 1:(dim(2)-1)            
        if identidad(:, i) == Matrix_problem(1:dim(1)-1, j)
            Order_initial(j) = i;
            iscorrect = 1;
            break;
        end           
    end
    if ~iscorrect
        errordlg('La matriz no está en forma canónica.','Problema no aplicable','modal');
        return;
    end        
end
Order_initial(Order_initial == 0) = dim(1):dim(2); % se etiquetan las demás columnas
% se verifica que la solución inicial sea factible
X0 = Matrix_problem(1:(dim(1)-1), end);
if X0 < 0
    iscorrect = 0;
    errordlg('La solución inicial no es factible.','Problema no aplicable','modal');
    return;
end
% se verifica que la solución no sea degenerada
if (X0 ~= 0)~= 1 
    iscorrect = 0;
    errordlg('La solución inicial es degenerada.','Problema no aplicable','modal');
    return;
end
iscorrect = 1;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.Problem);
