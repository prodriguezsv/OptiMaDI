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

if strcmp(get(handles.LPApphandle.Simplex_transportation, 'Checked'), 'off')
    setheading(handles, char('f', 'X', 'Yi0', 'Z'));
else
    setheading(handles, char('Origen', 'Destino', 'Oferta', 'Demanda'));
end
% UIWAIT makes problem wait for user response (see UIRESUME)
% uiwait(handles.Problem);

function setheading(handles, headers)
global Dimension;

% se obtienen las dimensiones del problema
Dimension = size(handles.LPApphandle.gui_Matrix_problem);

% se rotulan las filas y columnas de las tablas
colName = cell(Dimension(2), 1);
for i = 1:Dimension(2)-1
    colName(i) = cellstr(strcat(headers(2,:),num2str(i)));
end
colName(Dimension(2)) = cellstr(headers(3,:));
set(handles.table_problem, 'columnname', colName);

rowName=cell(1,Dimension(1));
for i = 1:Dimension(1)-1
    rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
end
rowName(Dimension(1)) = cellstr(headers(4,:));
set(handles.table_problem, 'rowname', rowName);
% se establece el formato de las columnas
colFormat=cell(1,Dimension(2));
for i = 1:(Dimension(2))
    colFormat(1,i) = cellstr('rat');
end
set(handles.table_problem, 'columnformat', colFormat);
% se configuran las columnas como editables
colEdit = ones(1,Dimension(2));
set(handles.table_problem, 'columneditable', (colEdit == 1));
% se despliega la tabla para el ingreso de datos
All_display = handles.LPApphandle.gui_Matrix_problem;
set(handles.table_problem, 'data', All_display);

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
% se verifica que la especificación este en forma correcta para el método
% correspondiente
if ~iscorrectform(handles)
    return;
end
   
handles.LPApphandle.Order_initial = Order_initial; % se comparte el orden inicial
handles.LPApphandle.setProblem(Matrix_problem, handles.LPApphandle); % se establece el problema actual
delete(handles.Problem);

% ---Funcion utilitaria: Verifica si se ha ingresado la formulación del
% problema de manera correcta
function iscorrect = iscorrectform(handles)
global Matrix_problem Order_initial Dimension;

Dimension = size(Matrix_problem);
Order_initial = zeros(1, Dimension(2));
if strcmp(get(handles.LPApphandle.Simplex_transportation, 'Checked'), 'off')
    if ~iscanonicform()
        iscorrect = 0;
        return;
    end
else
    oferta = Matrix_problem(1:(Dimension(1)-1), end);
    demanda = Matrix_problem(end, 1:(Dimension(2)-1));
    if sum(oferta)~= sum(demanda)
        errordlg('El problema no está balanceado.','Problema no balanceado','modal');
        iscorrect = 0;
        return;
    end
end
if ~isfactiblesolution(handles)
    iscorrect = 0;
    return;
end
if strcmp(get(handles.LPApphandle.Simplex_dual, 'Checked'), 'off')
    if isdegeneratedsolution(handles)
        errordlg('La solución inicial es degenerada.','Posibilidad de ciclo','modal');
    end
end
iscorrect = 1;

% -------
function response = isfactiblesolution(handles)
global Matrix_problem Dimension Order_initial;

% se verifica que la solución inicial sea factible
if strcmp(get(handles.LPApphandle.Simplex, 'Checked'), 'on')    
    X0 = Matrix_problem(1:(Dimension(1)-1), end);
    if any(X0 < 0)
        errordlg('La solución inicial no es primal factible.','Método Simplex Primal no aplicable','modal');
        response = 0;
        return;
    end
elseif strcmp(get(handles.LPApphandle.Simplex_dual, 'Checked'), 'on')
    % se convierte la última en términos de Rj =cj - zj
    [Y, I] = sort(Order_initial); %#ok<ASGLU>
    Tableau_sorted = Matrix_problem(:, I);
    c = Tableau_sorted(end,:);
    z = Tableau_sorted(end, 1:(Dimension(1)-1))*Tableau_sorted(1:(Dimension(1)-1),:);

    Rj = c - z;    
    if any(Rj < 0)
        response = 0;
        errordlg('La solución inicial no es dual factible.','Método Simplex Dual no aplicable','modal');
        return;
    end
elseif strcmp(get(handles.LPApphandle.Simplex_transportation, 'Checked'), 'on')
    oferta = Matrix_problem(1:(Dimension(1)-1), end);
    demanda = Matrix_problem(end, 1:(Dimension(2)-1));
    if any(oferta < 0) || any(demanda < 0)
        response = 0;
        errordlg('La solución inicial no será factible.','Método Simplex no aplicable','modal');
        return;
    end
end
response = 1;

% ---------

function response = isdegeneratedsolution(handles)
global Matrix_problem Dimension;

% se verifica que la solución no sea degenerada
if strcmp(get(handles.LPApphandle.Simplex_transportation, 'Checked'), 'off')
    X0 = Matrix_problem(1:(Dimension(1)-1), end);
    if any(X0 == 0)
        response = 1;    
        return;
    end
else
    oferta = Matrix_problem(1:(Dimension(1)-1), end);
    demanda = Matrix_problem(end, 1:(Dimension(2)-1));
    if any(oferta == 0) || any(demanda == 0)
        response = 1;
        return;
    end
end
response = 0;
    
    
% ---------
function response = iscanonicform()
global Matrix_problem Order_initial Dimension;

% se construye la matriz identidad
identidad = eye(Dimension(1)-1);
% se verifica que la matriz identidad este inmersa en la especifición del
% problema
for i =1:(Dimension(1)-1)
    response = 0;
    for j = 1:(Dimension(2)-1)            
        if identidad(:, i) == Matrix_problem(1:Dimension(1)-1, j)
            Order_initial(j) = i;
            response = 1;
            break;
        end           
    end
    if ~response
        errordlg('La matriz no está en forma canónica.','Método no aplicable','modal');
        return;
    end        
end
Order_initial(Order_initial == 0) = Dimension(1):Dimension(2); % se etiquetan las demás columnas
response = 1;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.Problem);
