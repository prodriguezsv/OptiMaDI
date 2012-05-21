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

function varargout = LPApp(varargin)
% LPAPP M-file for LPApp.fig
%      LPAPP, by itself, creates a new LPAPP or raises the existing
%      singleton*.
%
%      H = LPAPP returns the handle to a new LPAPP or the handle to
%      the existing singleton*.
%
%      LPAPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LPAPP.M with the given input arguments.
%
%      LPAPP('Property','Value',...) creates a new LPAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LPAppv2_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LPApp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LPApp_OpeningFcn, ...
                   'gui_OutputFcn',  @LPApp_OutputFcn, ...
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

% --- Executes just before LPApp is made visible.
function LPApp_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LPAppv1 (see VARARGIN)

% Choose default command line output for LPApp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui('center');
% UIWAIT makes LPApp wait for user response (see UIRESUME)
% uiwait(handles.LPApp);

% --- Outputs from this function are returned to the command line.
function varargout = LPApp_OutputFcn(hObject, eventdata, handles) %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Función utilitaria: configura los distintos ambientes de ejecución
function set_environment(environmentname, handles)
% environmentname especifica el tipo de ambiente por configurar
% handles         estructura con manejadores y datos de usuario

if strcmp(environmentname, 'next')
    set(handles.slider_increment, 'value', 0);
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['on ';'on ';'on ';'on ';'on ';'off'];
elseif strcmp(environmentname, 'sol_multiples')
    val_switches = ['on ';'on ';'on ';'off';'on ';'on '];
    set(handles.Saveproblem, 'enable', 'on');
elseif strcmp(environmentname, 'end')
    set(handles.slider_increment, 'value', 0);
    val_switches = ['off';'off';'on ';'off';'on '; 'off'];
    set(handles.Saveproblem, 'enable', 'on');
end

gui_environment(val_switches, handles);

% --- Función utilitaria: configura el ambiente en cada paso del Método
% Simplex
function gui_environment(val_switches, handles)
% val_switches    valores para la propiedad de cada control por configurar
% handles         estructura con manejadores y datos de usuario
set(handles.popupmenu_selectvar, 'Enable', val_switches(1,:));
set(handles.slider_increment, 'Enable', val_switches(2,:));
set(handles.pushbutton_start, 'Enable', val_switches(3,:));
set(handles.pushbutton_next, 'Enable', val_switches(4,:));
set(handles.pushbutton_modify, 'Enable', val_switches(5,:));
set(handles.pushbutton_multiplesolution, 'Enable', val_switches(6,:));


% --- Executes during object creation, after setting all properties.
function table_simplexdisplay_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to table_simplexdisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in pushbutton_start.
function pushbutton_start_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem Order_initial;
% Matrix_problem    especificación del problema
% Order_initial     orden inicial de vectores columna de la tabla Simplex, el orden de la
% matrix identidad primero y luego los demás vectores

% Configura la tabla inicial del Método Simplex
handles.Order_initial = Order_initial;
setProblem(Matrix_problem, handles)


% --- Función utilitaria: obtiene las variables no básicas candidatas a
% ser seleccionadas para convertirse en básicas según el criterio de
% mejoramiento de la solución (coeficientes de costo relativo negativos)
function calc_variables(handles)
% handles       estructura con manejadores y datos de usuario
global Dimension Tableau Order_current;
% Dimension     dimensiones de la especificación del problema en términos
% del número de ecuaciones y variables
% Tableau       Tabla del Simplex actual
% Order_current orden actual de vectores columna de la tabla Simplex

if strcmp(get(handles.Simplex, 'Checked'), 'on')
    Rj = Tableau(end,1:(Dimension(2)-1)); % Recupera los coeficientes de costo reducido (Rj) de la tabla Simplex
    [Y, I] = sort(Order_current); %#ok<ASGLU> % obtiene los indices de variables ordenados de menor a mayor
    IRjneg = I(Order_current(Rj < 0));  % obtiene los indices de variables con Rj negativo
    IVar = IRjneg;
elseif strcmp(get(handles.Simplex_dual, 'Checked'), 'on')
    X0 = Tableau(1:(Dimension(1)-1), end); % Recupera los valores de las variables básicas de la tabla Simplex
    [Y, I] = sort(Order_current); %#ok<ASGLU> % obtiene los indices de variables ordenados de menor a mayor
    IX0neg = I(X0 < 0);  % obtiene los indices de variables con Rj negativo
    IVar =  IX0neg;
end

% se construye el arreglo de cadenas de las variables
dim2 = size(IVar);
Var = cell(dim2(2), 1);
for i = 1:dim2(2)
    Var(i) = cellstr(strcat('X',num2str(IVar(i))));
end

if ~isempty(Var) % en el caso que haya Rj negativos
    if strcmp(get(handles.Simplex, 'Checked'), 'on')
        [Y, p] = min(Rj(IRjneg)); %#ok<ASGLU> % obtiene el índice de variable del Rj más negativo
    else strcmp(get(handles.Simplex_dual, 'Checked'), 'on')
        [Y, p] = min(IX0neg); %#ok<ASGLU> % obtiene el índice de variable del Rj más negativo
    end
    set(handles.popupmenu_selectvar, 'string', char(Var)); % despliega las variables no básicas en la interfaz
    set_environment('next', handles); % ajusta el ambiente
    set(handles.popupmenu_selectvar, 'value', p); % se selecciona, por omisión, la variable más negativa
    calc_ratios(handles); % se calculan las razones para aplicar el criterio de factibilidad
elseif strcmp(get(handles.Simplex, 'Checked'), 'on') % verificar si hay Rj iguales a cero
    j = 1;
    IRjeqaux = I(Rj(I(Dimension(1):(Dimension(2)-1)))==0);
    dim3 = size(IRjeqaux);
    IRjeq = zeros(1, dim3(2));
    for i = Dimension(1):(Dimension(2)-1)
        if Rj(I(i)) == 0
            IRjeq(j) = I(i);
            j = j+1;
        end
    end

    dim2 = size(IRjeq);
    Var = cell(dim2(2), 1);
    for i = 1:dim2(2)
        Var(i) = cellstr(strcat('X',num2str(IRjeq(i))));
    end
    if ~isempty(Var) % en el caso que haya Rj = 0        
        set(handles.popupmenu_selectvar, 'string', char(Var));
        set_environment('sol_multiples', handles);
        set(handles.popupmenu_selectvar, 'value', 1);
        calc_ratios(handles);
    else % no hay soluciones multiples, no hay más opciones
        set(handles.popupmenu_selectvar, 'string', char('',''));
        set_environment('end', handles);
    end
end


% --- Función utilitaria: calcula las razones del criterio de factibilidad
function calc_ratios(handles)
% handles   estructura con manejadores y datos de usuario
global Tableau Dimension Order_current;
% All_display condensa los datos que se desplegarán en la interfaz

var = get(handles.popupmenu_selectvar, 'string'); % se recupera la variable no básica seleccionada
J = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2)); % se recupera el índice de variable seleccionada    
if strcmp(get(handles.Simplex, 'Checked'), 'on')
    % Se calculan las razones para la variable no básica seleccionada
    Y0 = Tableau(1:(Dimension(1)-1), end); % se recupera el vector de términos coeficientes
    Yj = Tableau(1:(Dimension(1)-1), J); % se recupera el vector columna asociada a la variable no básica
    ratios = Y0./Yj;
    All_display = zeros(Dimension(1), Dimension(2)+1);
    % se actualiza la tabla de la interfaz
    All_display(:, 1:Dimension(2)) = Tableau;
    All_display(1:(Dimension(1)-1), end) = ratios;
elseif strcmp(get(handles.Simplex_dual, 'Checked'), 'on')
    % Se calculan las razones para la variable básica seleccionada
    Fi = Tableau(Order_current(J), 1:(Dimension(2)-1)); % se recupera el vector fila asociada a la variable básica
    Rj = Tableau(end, 1:(Dimension(2)-1)); % se recupera el vector de coeficientes reducidos
    ratios = -Rj./Fi;
    All_display = zeros(Dimension(1)+1, Dimension(2));
    % se actualiza la tabla de la interfaz
    All_display(1:Dimension(1), :) = Tableau;
    All_display(end, 1:(Dimension(2)-1)) = ratios;
end
set(handles.table_simplexdisplay, 'data', All_display);  

% --- Executes on button press in pushbutton_next.
function pushbutton_next_Callback(hObject, eventdata, handles) %#ok<INUSL>
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
global Dimension;

% se recuperan las razones para aplicar el criterio de factibilidad
All_display = get(handles.table_simplexdisplay, 'data');
if strcmp(get(handles.Simplex, 'Checked'), 'on')
    ratios = All_display(1:(Dimension(1)-1),end);
    % se ejecuta el método simplex
    simplex_primal(ratios, handles);
elseif strcmp(get(handles.Simplex_dual, 'Checked'), 'on')   
    ratios = All_display(end, 1:(Dimension(2)-1));
    simplex_dual(ratios, handles);
end

calc_variables(handles);
    
% --- Se ejecuta el método simplex primal
function simplex_primal(ratios, handles)
global Tableau;

% recupera el índice de la variable no básica seleccionada asociada a la
% columna que ingresará a la base
var = get(handles.popupmenu_selectvar, 'string');
q = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2));

% recupera el índice de la fila (o componente) del vector columna que
% dejará la base según la razón mínima tal que sea positiva
ratios_aux = ratios;
ratios_aux(ratios <= 0) = Inf; 
[C, p] = min(ratios_aux); 

if C ~= Inf % si existe algún yij que cumple el criterio de factibilidad
    pivote_process(p, q);
else
    errordlg('El conjunto representado por las restricciones no es acotado.', 'El valor objetivo decrece sin limite','modal');
    return;
end
% se actualiza la tabla de la interfaz
All_display = Tableau;
set(handles.table_simplexdisplay, 'data', All_display);

% --- Se ejecuta el método simplex primal
function simplex_dual(ratios, handles)
global Tableau;

% recupera el índice de la variable no básica seleccionada asociada a la
% columna que ingresará a la base
var = get(handles.popupmenu_selectvar, 'string');
p = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2));

% recupera el índice de la fila (o componente) del vector columna que
% dejará la base según la razón mínima tal que sea positiva
ratios_aux = ratios;
ratios_aux(ratios <= 0) = Inf; 
[C, q] = min(ratios_aux); 

if C ~= Inf % si existe algún yij que cumple el criterio de factibilidad
    pivote_process(p, q);
else
    errordlg('El conjunto representado por es vacío.', 'No hay solución','modal');
    return;
end
% se actualiza la tabla de la interfaz
All_display = Tableau;
set(handles.table_simplexdisplay, 'data', All_display);


% --- Función utilitaria: proceso del pivote sobre el elemento Ypq
function pivote_process(p, q)
global Tableau Order_current Dimension;

% Ecuaciones de pivote
Ypq = Tableau(p, q);
fp = Tableau(p, :); % recupera la fila p actual
fp = fp/Ypq;
Tableau(p, :) = fp; % se actualiza la fila p en la tabla del Simplex
% se actualizan las demás filas
for i = 1:Dimension(1)
    if i ~= p
        fi = Tableau(i, :);
        Yiq = Tableau(i, q);
        fi = fi - fp*Yiq;
        Tableau(i, :) = fi; % fila p actualizada
    end
end

% intercambian roles las variables
Order_current(1, Order_current(1,:) == p) = Order_current(1,q);
Order_current(1, q) = p;
    
% --- Executes on button press in pushbutton_modify.
function pushbutton_modify_Callback(hObject,eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem;

% abre ventana para modificar la especificación del problema
handles.setProblem = @setProblem; % se comparte manejador de función
handles.gui_Matrix_problem = Matrix_problem; % se comparte la especificación de problema actual
handles.gui_Problem = Problem('LPApp', handles);
guidata(handles.output, handles);

% --------------------------------------------------------------------
function Simplex_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject, 'Checked'),'off')  
    set(hObject,'Checked','on');
end
set(handles.Simplex_dual,'Checked','off');

% --- Executes on selection change in popupmenu_selectvar.
function popupmenu_selectvar_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to popupmenu_selectvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

calc_ratios(handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_selectvar_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to popupmenu_selectvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider_increment_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to slider_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Tableau Dimension All_display;

% se recupera la variable por incrementar y el incremento
incr = get(hObject, 'Value');
var = get(handles.popupmenu_selectvar, 'string');
J = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2));

% se actualizan los datos desplegados
Aux_Tableau = Tableau;
All_display = get(handles.table_simplexdisplay, 'data');
ratios = All_display(1:(Dimension(1)-1), end);
ratios_aux = ratios;
ratios_aux(ratios <= 0) = Inf; 
[ratio, Y] = min(ratios_aux);  %#ok<NASGU>
% se actualizan las variables básicas según el incremento de la variable no
% básica
Xb = Aux_Tableau(:,end);
Yj = Tableau(:, J);
Aux_Tableau(:,end) = Xb - Yj*incr/100*ratio;
% se actualiza la tabla de la interfaz
All_display = zeros(Dimension(1), Dimension(2)+1);
All_display(:, 1:Dimension(2)) = Aux_Tableau;
All_display(1:(Dimension(1)-1), end) = ratios;
set(handles.table_simplexdisplay, 'data', All_display);

% --- Executes during object creation, after setting all properties.
function slider_increment_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to slider_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pushbutton_multiplesolution.
function pushbutton_multiplesolution_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to pushbutton_multiplesolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_next_Callback(hObject, eventdata, handles)

% --- se ejecuta al seleccionar la opción "Nuevo poblema" del menú "Archivo"
function Newproblem_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Newproblem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ventana de dialogo en donde se solicita las dimensiones del problema
prompt = {'Número de ecuaciones:','Número de variables:'};
dlg_title = 'Dimensiones del problema';
num_lines = 1;
def = {'1','2'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

% se verifica que las dimensiones del problema sean consistentes (m < n)
if ~isempty(answer)
    user_entry1 = str2double(cellstr(answer{1}));
    user_entry2 = str2double(cellstr(answer{2}));

    while (1)
        if (isnan(user_entry1) || user_entry1 < 1) || (isnan(user_entry2) || user_entry2 < 2) || user_entry1 > user_entry2
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            user_entry1 = str2double(cellstr(answer{1}));
            user_entry2 = str2double(cellstr(answer{2}));
        else
            break;
        end
    end
    % se abre la ventana para introducir la especificación del problema
    handles.gui_Matrix_problem = zeros(user_entry1+1, user_entry2+1); % se comparte la especificación del problema
    handles.setProblem = @setProblem; % se comparte el manejador de función    
    handles.gui_Problem = Problem('LPApp', handles);
    guidata(handles.output, handles);
end

% --- actualiza la especificación del problema
function setProblem(Problem, handles)
% Problem especificación del problema
% handles estructura de manejadores y datos de usuario
global Dimension Matrix_problem Tableau Order_current Order_initial All_display;

% se actualizan algunas variables globales previamente descritas
Matrix_problem = Problem;
Dimension = size(Matrix_problem);
Order_initial = handles.Order_initial;
Order_current = Order_initial;
Tableau = Matrix_problem; 

% se convierte la última en términos de Rj =cj - zj
[Y, I] = sort(Order_initial); %#ok<ASGLU>
Tableau_sorted = Matrix_problem(:, I);
c = Tableau_sorted(end,:);
z = Tableau_sorted(end, 1:(Dimension(1)-1))*Tableau_sorted(1:(Dimension(1)-1),:);

R = c - z;
Tableau(end, I) = R;
% se actualiza la tabla de la interfaz
All_display = zeros(Dimension(1), Dimension(2));
All_display(:, 1:Dimension(2)) = Tableau;
set(handles.table_simplexdisplay, 'data', All_display); 
% se rotulan las columnas y las filas de manera
% correspondiente
setTableauTags(handles);
% se calculan las nuevas variables no básicas si las hay
calc_variables(handles);


% ----se rotulan las columnas y las filas de manera
% correspondiente
function setTableauTags(handles)
global Dimension;

colName = cell(Dimension(2)+1, 1);
for i = 1:(Dimension(2)-1)
    colName(i) = cellstr(strcat('X',num2str(i)));
end
colName(Dimension(2)) = cellstr('Yi0'); 
if strcmp(get(handles.Simplex, 'Checked'), 'on')
    colName(Dimension(2)+1) = cellstr('Yi0/Yij');
    colFormat=cell(1,Dimension(2)+1);
    colFormat(1,Dimension(2)+1) = cellstr('rat');
end
set(handles.table_simplexdisplay, 'columnname', colName);

rowName=cell(1,Dimension(1));
for i = 1:(Dimension(1)-1)
    rowName(1,i) = cellstr(strcat('f',num2str(i)));
end
rowName(Dimension(1)) = cellstr('Rj');

if strcmp(get(handles.Simplex_dual, 'Checked'), 'on')
    rowName(Dimension(1)+1) = cellstr('-Rj/Yij');
    colFormat=cell(1,Dimension(2));
end
set(handles.table_simplexdisplay, 'rowname', rowName);

% se especifica el formato de salida de los datos
for i = 1:(Dimension(2))
    colFormat(1,i) = cellstr('rat');
end
set(handles.table_simplexdisplay, 'columnformat', colFormat);


% --- Executes when user attempts to close LPApp.
function LPApp_CloseRequestFcn(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to LPApp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% se cierran todas las ventanas abiertas hijas de la aplicación principal

try
    delete(handles.gui_Problem);
catch %#ok<CTCH>    
end
delete(hObject);

% --------------------------------------------------------------------
function Loadproblem_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Loadproblem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, Y, Z] = uigetfile('*.mat', 'Selecciona un nombre de archivo válido'); %#ok<ASGLU,NASGU>
if filename ~= 0
    S = load(filename);
    handles.setProblem = @setProblem;    
    handles.gui_Matrix_problem = S.Matrix_problem;
    handles.gui_Problem = Problem('LPApp', handles);
    guidata(handles.output, handles);
end

% --------------------------------------------------------------------
function Saveproblem_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to Saveproblem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem; %#ok<NUSED>

uisave('Matrix_problem', 'LProblem');


% --------------------------------------------------------------------
function Simplex_dual_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSL>
% hObject    handle to Simplex_dual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject, 'Checked'),'off')  
    set(hObject,'Checked','on');
end
set(handles.Simplex,'Checked','off');
