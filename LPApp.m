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
    set(handles.Sensibility, 'enable', 'off');
    set(handles.Postoptimality, 'enable', 'off');
    set(handles.slider_increment, 'value', 0);
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['on ';'on ';'on ';'on ';'on ';'off'];   
elseif strcmp(environmentname, 'sol_multiples')
    set(handles.Sensibility, 'enable', 'on');
    set(handles.Postoptimality, 'enable', 'on');
    val_switches = ['on ';'off';'on ';'off';'on ';'on '];
    set(handles.Saveproblem, 'enable', 'on');
elseif strcmp(environmentname, 'end')
    set(handles.Sensibility, 'enable', 'on');
    set(handles.Postoptimality, 'enable', 'on');
    set(handles.slider_increment, 'value', 0);
    val_switches = ['off';'off';'on ';'off';'on '; 'off'];
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.pushbutton_asign, 'Enable', 'off');
    set(handles.pushbutton_asignall, 'Enable', 'off');
    set(handles.pushbutton_watchsolution, 'visible', 'on');    
elseif strcmp(environmentname, 'next_assign')
    set(handles.Sensibility, 'enable', 'off');
    set(handles.Postoptimality, 'enable', 'off');
    setTableauTags(handles, char('Origen', 'Destino', 'Demanda', 'Oferta'));
    set(handles.panel_transportation, 'title', 'Encontrar solución inicial');
    set(handles.pushbutton_asign, 'string', 'Asignar');
    set(handles.pushbutton_asignall, 'string', 'Asignar todo');
    set(handles.popupmenu_selectvar2, 'Enable', 'on');
    set(handles.pushbutton_asign, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on');
    set(handles.pushbutton_watchsolution, 'visible', 'off');  
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['on ';'on ';'on ';'off';'on ';'off'];   
elseif strcmp(environmentname, 'next_calc')    
    setTableauTags(handles, char('Origen', 'Destino', 'V', 'U'));
    set(handles.panel_transportation, 'title', 'Calcular coeficientes relativos');
    set(handles.text_selectvar2, 'string', 'Seleccionar multiplicador inicial');
    set(handles.pushbutton_asign, 'string', 'Calcular');
    set(handles.pushbutton_asignall, 'string', 'Calcular todo');   
    set(handles.popupmenu_selectvar2, 'Enable', 'on');
    set(handles.pushbutton_asign, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on');
    set(handles.pushbutton_watchsolution, 'visible', 'on');
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['on ';'on ';'on ';'off';'on ';'off']; 
elseif strcmp(environmentname, 'next_cicle')
    setTableauTags(handles, char('Origen', 'Destino', 'Demanda', 'Oferta'));
    set(handles.panel_transportation, 'title', 'Calcular ciclo de cambio');
    set(handles.text_selectvar2, 'string', 'Seleccionar variable que entra');
    set(handles.pushbutton_asign, 'string', 'Buscar');
    set(handles.pushbutton_asignall, 'string', 'Buscar todo');
    set(handles.popupmenu_selectvar2, 'Enable', 'on');
    set(handles.pushbutton_asign, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on'); 
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['on ';'on ';'on ';'off';'on ';'off'];   
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
global Dimension Tableau Order_current Solution_initial;
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

if strcmp(get(handles.Simplex_transportation, 'Checked'), 'off')
    % se construye el arreglo de cadenas de las variables
    dim2 = size(IVar);
    Var = cell(dim2(2), 1);
    for i = 1:dim2(2)
        Var(i) = cellstr(strcat('X',num2str(IVar(i))));
    end

    if ~isempty(Var) % en el caso que haya Rj negativos
        if strcmp(get(handles.Simplex, 'Checked'), 'on')
            [Y, p] = min(Rj(IRjneg)); %#ok<ASGLU> % obtiene el índice de variable del Rj más negativo
        elseif strcmp(get(handles.Simplex_dual, 'Checked'), 'on')
            [Y, p] = min(IX0neg); %#ok<ASGLU> % obtiene el índice de variable del Rj más negativo
        end
        set(handles.popupmenu_selectvar, 'string', char(Var)); % despliega las variables no básicas en la interfaz
        set_environment('next', handles); % ajusta el ambiente
        set(handles.popupmenu_selectvar, 'value', p); % se selecciona, por omisión, la variable más negativa
        calc_ratios(handles); % se calculan las razones para aplicar el criterio de factibilidad
    else % verificar si hay Rj iguales a cero
        j = 1;
        Rj = Tableau(end,1:(Dimension(2)-1)); % Recupera los coeficientes de costo reducido (Rj) de la tabla Simplex
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
        if ~isempty(Var) && strcmp(get(handles.Simplex, 'Checked'), 'on') % en el caso que haya Rj = 0        
            set(handles.popupmenu_selectvar, 'string', char(Var));
            set_environment('sol_multiples', handles);
            set(handles.popupmenu_selectvar, 'value', 1);
            calc_ratios(handles);
        else % no hay Soluciones multiples, no hay más opciones
            set(handles.popupmenu_selectvar, 'string', char('',''));
            set_environment('end', handles);
        end
    end
else
    % se construye el arreglo de cadenas de las variables    
    Var = cell(Dimension(2)-1, 1);
    for i = 1:(Dimension(2)-1)
        Var(i) = cellstr(strcat('X1',num2str(i)));
    end
    Solution_initial = 1;
    set_environment('next_assign', handles); % ajusta el ambiente
    set(handles.popupmenu_selectvar2, 'string', char(Var));
    set(handles.popupmenu_selectvar2, 'value', 1);
    calc_nextassignment(handles);
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

% ----------------
function calc_nextassignment(handles)
global Solution Node_current Num_assignation T_Tableau Dimension Matrix_problem T_VarType ...
    Solution_change Minimo Solution_initial Node_NBV;

if Solution_initial == 1
    if all(Node_current==0)
        Num_assignation = 1;
        T_VarType = zeros(Dimension(1)-1, Dimension(2)-1);
        var = get(handles.popupmenu_selectvar2, 'string'); % se recupera la variable no básica seleccionada
        I = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 2)); % se recupera el índice de variable seleccionada 
        J = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 3)); % se recupera el índice de variable seleccionada 
        Node_current = [I, J];
        T_VarType(I, J) = 1;
        Solution(Num_assignation, 1, 1)  = Node_current(1);
        Solution(Num_assignation, 1, 2)  = Node_current(2);
        Limites = [T_Tableau(end, J), T_Tableau(I, end)];
        [Solution(Num_assignation, 1, 3), i]  = min(Limites);
        if (i == 1)        
            T_Tableau(I, end) = T_Tableau(I, end) - T_Tableau(end, J);
            T_Tableau(end, J) =  0;
        else
            T_Tableau(end, J) = T_Tableau(end, J) - T_Tableau(I, end);
            T_Tableau(I, end) =  0;
        end
    else
        Num_assignation = Num_assignation + 1;    
        if T_Tableau(Node_current(1), end) ~= 0
            Solution(Num_assignation, 1, 2)  = mod(Node_current(2),Dimension(2)-1)+1;
            Solution(Num_assignation, 1, 1)  = Node_current(1);
        else
            Solution(Num_assignation, 1, 2)  = Node_current(2);
            Solution(Num_assignation, 1, 1)  = Node_current(1)+1;
        end
        Node_current = [Solution(Num_assignation, 1, 1), Solution(Num_assignation, 1, 2)];
        Limites = [T_Tableau(end, Node_current(2)), T_Tableau(Node_current(1), end)];
        [Solution(Num_assignation, 1, 3), i]  = min(Limites);
        if (i == 1)        
            T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), end) - T_Tableau(end, Node_current(2));
            T_Tableau(end, Node_current(2)) =  0;
        else
            T_Tableau(end, Node_current(2)) = T_Tableau(end, Node_current(2)) - T_Tableau(Node_current(1), end);
            T_Tableau(Node_current(1), end) =  0;
        end        
    end
    T_VarType(Node_current(1), Node_current(2)) = 1;

    All_display = zeros(Dimension(1), Dimension(2));
    for j=1:Num_assignation
        All_display(Solution(j, 1, 1), Solution(j, 1, 2)) = Solution(j, 1, 3);
    end  
    All_display(end, :) = T_Tableau(end, :);
    All_display(:, end) = T_Tableau(:, end);
    set(handles.table_simplexdisplay, 'data', All_display);
    if Num_assignation == Dimension(1)+Dimension(2)-3
        Node_current = [0, 0];
        set_environment('next_calc', handles);
        T_Tableau = zeros(Dimension(1), Dimension(2));
        T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);
        % se construye el arreglo de cadenas de las variables    
        Var = cell(4, 1);    
        Var(1) = cellstr(strcat('U',num2str(Solution(1, 1, 1))));   
        Var(2) = cellstr(strcat('V',num2str(Solution(1, 1, 2))));   
        Var(3) = cellstr(strcat('U',num2str(Solution(Dimension(1)+Dimension(2)-3, 1, 1))));   
        Var(4) = cellstr(strcat('V',num2str(Solution(Dimension(1)+Dimension(2)-3, 1, 2))));       
        set(handles.popupmenu_selectvar2, 'string', char(Var));
        set(handles.popupmenu_selectvar2, 'value', 4);    
    end
else
    var = get(handles.popupmenu_selectvar2, 'string'); % se recupera la variable no básica seleccionada
    I = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 2)); % se recupera variable seleccionada 
    J = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 3)); % se recupera el índice de variable seleccionada
    Node_BV = [I, J];
    if all(Node_current==0)
        Num_assignation = 1;
        T_VarType = zeros(Dimension(1)-1, Dimension(2)-1);
        Node_current = [Solution(Num_assignation, 1, 1), Solution(Num_assignation, 1, 2)];        
        if Node_current(1) == Node_BV(1) && Node_current(2) == Node_BV(2)
           Solution(Num_assignation, 1, 1) = Node_NBV(1);
           Solution(Num_assignation, 1, 2) = Node_NBV(2);
           Solution(Num_assignation, 1, 3) = Minimo;
           T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), end) - Minimo;
           T_Tableau(end, Node_current(2)) = T_Tableau(end, Node_current(2)) - Minimo;
           Node_current = [Solution(Num_assignation, 1, 1), Solution(Num_assignation, 1, 2)];
           T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), end) + Minimo;
           T_Tableau(end, Node_current(2)) = T_Tableau(end, Node_current(2)) + Minimo;
        else
            Solution(Num_assignation, 1, 3) = Solution(Num_assignation, 1, 3) + Solution_change(Num_assignation)*Minimo;
            T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), end) + Solution_change(Num_assignation)*Minimo;
            T_Tableau(end, Node_current(2)) = T_Tableau(end, Node_current(2)) + Solution_change(Num_assignation)*Minimo;
        end        
    else
        Num_assignation = Num_assignation+1;
        Node_current = [Solution(Num_assignation, 1, 1), Solution(Num_assignation, 1, 2)];        
        if Node_current(1) == Node_BV(1) && Node_current(2) == Node_BV(2)
           Solution(Num_assignation, 1, 1) = Node_NBV(1);
           Solution(Num_assignation, 1, 2) = Node_NBV(2);
           Solution(Num_assignation, 1, 3) = Minimo;
           T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), end) - Minimo;
           T_Tableau(end, Node_current(2)) = T_Tableau(end, Node_current(2)) - Minimo;
           Node_current = [Solution(Num_assignation, 1, 1), Solution(Num_assignation, 1, 2)];
           T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), end) + Minimo;
           T_Tableau(end, Node_current(2)) = T_Tableau(end, Node_current(2)) + Minimo;
        else            
            Solution(Num_assignation, 1, 3) = Solution(Num_assignation, 1, 3) + Solution_change(Num_assignation)*Minimo;
            T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), end) + Solution_change(Num_assignation)*Minimo;
            T_Tableau(end, Node_current(2)) = T_Tableau(end, Node_current(2)) + Solution_change(Num_assignation)*Minimo;
        end        
    end
    T_VarType(Node_current(1), Node_current(2)) = 1;
    
    All_display = zeros(Dimension(1), Dimension(2));
    for j=1:Num_assignation
        All_display(Solution(j, 1, 1), Solution(j, 1, 2)) = Solution(j, 1, 3);
    end  
    All_display(end, :) = T_Tableau(end, :);
    All_display(:, end) = T_Tableau(:, end);
    set(handles.table_simplexdisplay, 'data', All_display);
    if Num_assignation == Dimension(1)+Dimension(2)-3
        Node_current = [0, 0];
        set_environment('next_calc', handles);
        T_Tableau = zeros(Dimension(1), Dimension(2));
        T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);
        % se construye el arreglo de cadenas de las variables    
        Var = cell(4, 1);    
        Var(1) = cellstr(strcat('U',num2str(Solution(1, 1, 1))));   
        Var(2) = cellstr(strcat('V',num2str(Solution(1, 1, 2))));   
        Var(3) = cellstr(strcat('U',num2str(Solution(Dimension(1)+Dimension(2)-3, 1, 1))));   
        Var(4) = cellstr(strcat('V',num2str(Solution(Dimension(1)+Dimension(2)-3, 1, 2))));       
        set(handles.popupmenu_selectvar2, 'string', char(Var));
        set(handles.popupmenu_selectvar2, 'value', 4);    
    end
end

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
setrowheaders(handles, char('X', 'X', 'Rj', 'Yi0'));
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
ratios_aux(ratios < 0 | isnan(ratios)) = Inf; 
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
global Tableau Order_current;

% recupera el índice de la variable no básica seleccionada asociada a la
% columna que ingresará a la base
var = get(handles.popupmenu_selectvar, 'string');
p_aux = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2));
p = Order_current(p_aux);

% recupera el índice de la fila (o componente) del vector columna que
% dejará la base según la razón mínima tal que sea positiva
ratios_aux = ratios;
ratios_aux(ratios < 0 | isnan(ratios)) = Inf; 
[C, q] = min(ratios_aux); 

if C ~= Inf % si existe algún yij que cumple el criterio de factibilidad
    pivote_process(p, q);
else
    errordlg('El conjunto representado por las restricciones es vacío.', 'No hay solución','modal');
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
set(handles.Simplex_transportation,'Checked','off');
set(handles.panel_enhancement, 'visible', 'on');
set(handles.panel_transportation, 'visible', 'off');
set(handles.panel_enhancement, 'title', 'Criterio de mejoramiento (Simplex primal)');
set(handles.text_selectvar, 'string', 'Seleccionar variable que entra');

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
global Tableau Dimension Order_initial Order_current Matrix_problem;

% se recupera la variable por incrementar y el incremento
incr = get(hObject, 'Value');
var = get(handles.popupmenu_selectvar, 'string');
J = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2));

% se actualizan los datos desplegados
Aux_Tableau = Tableau;
All_display = get(handles.table_simplexdisplay, 'data');

if strcmp(get(handles.Simplex, 'Checked'), 'on')
    ratios = All_display(1:(Dimension(1)-1), end);
    ratios_aux = ratios;
    ratios_aux(ratios < 0) = Inf; 
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
elseif strcmp(get(handles.Simplex_dual, 'Checked'), 'on')
    ratios = All_display(end, 1:(Dimension(2)-1));
    ratios_aux = ratios;    
    Vector_nonbasic = Order_current(1:Dimension(2)-1) >= Dimension(1);
    ratios_aux(ratios < 0 & Vector_nonbasic) = Inf;
    ratios_aux(ratios == 0 & ~Vector_nonbasic) = Inf; 
    [ratio, Y] = min(ratios_aux);  %#ok<NASGU>
    % se actualizan los coeficientes de costo reducido según la
    % actualización de la solución del dual
    Rj = Aux_Tableau(end,:);
    A = Matrix_problem(1:Dimension(1)-1,1:Dimension(2));
    [Z, I] = sort(Order_initial); %#ok<ASGLU>
    B_inverse = Tableau(:, I(1:Dimension(1)-1));
    U = B_inverse(Order_current(J), 1:(Dimension(1)-1));
    Aux_Tableau(end, :) = Rj + U*A*incr/100*ratio;
    % se actualiza la tabla de la interfaz
    All_display = zeros(Dimension(1)+1, Dimension(2));
    All_display(1:Dimension(1), 1:Dimension(2)) = Aux_Tableau;
    All_display(end, 1:(Dimension(2)-1)) = ratios;
    set(handles.table_simplexdisplay, 'data', All_display);
end


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
if strcmp(get(handles.Simplex_transportation, 'Checked'), 'off')
    prompt = {'Número de ecuaciones:','Número de variables:'};
    dlg_title = 'Dimensiones del problema';
    num_lines = 1;
    def = {'1','2'};
else
    prompt = {'Número de origenes:','Número de destinos:'};
    dlg_title = 'Dimensiones del problema';
    num_lines = 1;
    def = {'2','2'};
end

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
global Dimension Matrix_problem Tableau Order_current Order_initial Node_current Solution T_Tableau;

% se actualizan algunas variables globales previamente descritas
Matrix_problem = Problem;
Dimension = size(Matrix_problem);

if strcmp(get(handles.Simplex_transportation, 'Checked'), 'off')
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
else
    Node_current = [0,0];
    Solution = zeros(Dimension(1)+Dimension(2)-3, 1, 3);
    T_Tableau = Matrix_problem;
end
% se actualiza la tabla de la interfaz
All_display = zeros(Dimension(1), Dimension(2));
if strcmp(get(handles.Simplex_transportation, 'Checked'), 'off')
    % se rotulan las columnas y las filas de manera correspondiente
    setTableauTags(handles, char('X', 'X', 'Rj', 'Yi0'));
    All_display(:, 1:Dimension(2)) = Tableau;
else
    % se rotulan las columnas y las filas de manera correspondiente
    setTableauTags(handles, char('Origen', 'Destino', 'Demanda', 'Oferta'));
    All_display(:, 1:Dimension(2)) = T_Tableau;
end
set(handles.table_simplexdisplay, 'data', All_display); 
% se calculan las nuevas variables no básicas si las hay
calc_variables(handles);


% ----se rotulan las columnas y las filas de manera
% correspondiente
function setTableauTags(handles, headers)
global Dimension;

colName = cell(Dimension(2)+1, 1);
for i = 1:(Dimension(2)-1)
    colName(i) = cellstr(strcat(headers(2,:),num2str(i)));
end
colName(Dimension(2)) = cellstr(headers(4,:)); 
if strcmp(get(handles.Simplex, 'Checked'), 'on')
    colName(Dimension(2)+1) = cellstr('Yi0/Yij');    
end
set(handles.table_simplexdisplay, 'columnname', colName);

setrowheaders(handles, headers);

if strcmp(get(handles.Simplex, 'Checked'), 'on')
    colFormat=cell(1,Dimension(2)+1);    
else
    colFormat=cell(1,Dimension(2));
end
% se especifica el formato de salida de los datos
colFormat(1,:) = cellstr('rat');

set(handles.table_simplexdisplay, 'columnformat', colFormat);

% -------
function setrowheaders(handles, headers)
global Dimension Order_current;

[Z, I] = sort(Order_current); %#ok<ASGLU>
rowName=cell(1,Dimension(1));
for i = 1:(Dimension(1)-1)
    if strcmp(get(handles.Simplex_transportation, 'Checked'), 'off')
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(I(i))));
    else
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end
end
rowName(Dimension(1)) = cellstr(headers(3,:));

if strcmp(get(handles.Simplex_dual, 'Checked'), 'on')
    rowName(Dimension(1)+1) = cellstr('-Rj/Yij');   
end
set(handles.table_simplexdisplay, 'rowname', rowName);

% --- Executes when user attempts to close LPApp.
function LPApp_CloseRequestFcn(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to LPApp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% se cierran todas las ventanas abiertas hijas de la aplicación principal

try
    delete(handles.gui_Problem);
    delete(handles.gui_Postoptimality);
    delete(handles.gui_Sensibility);
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
set(handles.Simplex_transportation,'Checked','off');
set(handles.panel_enhancement, 'visible', 'on');
set(handles.panel_transportation, 'visible', 'off');
set(handles.panel_enhancement, 'title', 'Criterio de mejoramiento (Simplex dual)');
set(handles.text_selectvar, 'string', 'Seleccionar variable que sale');

% ----------
function setProblemAndTableau(Problem, NewTableau, handles)
global Matrix_problem Tableau;

set(handles.table_simplexdisplay, 'data', NewTableau);
Tableau = NewTableau;
Matrix_problem = Problem;
calc_variables(handles);

% --------------------------------------------------------------------
function Postoptimality_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Postoptimality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem Tableau Order_current Order_initial;

handles.gui_tableau = Tableau;
handles.gui_Matrix_problem = Matrix_problem;
handles.Order_current = Order_current;
handles.Order_initial = Order_initial;
handles.setProblemAndTableau = @setProblemAndTableau;
handles.gui_postoptimality = Postoptimality('LPApp', handles);
guidata(handles.output, handles);

% --------------------------------------------------------------------
function Sensibility_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Sensibility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem Tableau Order_current Order_initial;

handles.gui_tableau = Tableau;
handles.gui_Matrix_problem = Matrix_problem;
handles.Order_current = Order_current;
handles.Order_initial = Order_initial;
handles.gui_sensibility = Sensibility('LPApp', handles);
guidata(handles.output, handles);

% --------------------------------------------------------------------
function Simplex_transportation_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_transportation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject, 'Checked'),'off')  
    set(hObject,'Checked','on');
end
set(handles.Simplex,'Checked','off');
set(handles.Simplex_dual,'Checked','off');
set(handles.panel_enhancement, 'visible', 'off');
set(handles.panel_transportation, 'visible', 'on');


% --- Executes on selection change in popupmenu_selectvar2.
function popupmenu_selectvar2_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to popupmenu_selectvar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Node_current Solution T_Tableau Matrix_problem Dimension Num_assignation T_VarType_Aux ...
    Solution_initial T_VarType NewSolution;

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selectvar2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selectvar2

% se calculan las nuevas variables no básicas si las hay
if strcmp(get(handles.pushbutton_asign, 'string'), 'Calcular')
    Node_current = [0, 0];
    set_environment('next_calc', handles);
    T_Tableau = zeros(Dimension(1), Dimension(2));
    T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);
  
    calc_nextmultiplier(handles);    
elseif strcmp(get(handles.pushbutton_asign, 'string'), 'Asignar')
    Node_current = [0,0];
    if Solution_initial == 1
        Solution = zeros(Dimension(1)+Dimension(2)-3, 1, 3);
    else
        NewSolution = Solution;
    end
    T_Tableau = Matrix_problem;

    % se actualiza la tabla de la interfaz
    All_display = zeros(Dimension(1), Dimension(2));
    set(handles.table_simplexdisplay, 'data', All_display); 
    calc_nextassignment(handles);
elseif strcmp(get(handles.pushbutton_asign, 'string'), 'Buscar')
    Num_assignation = 0;
    T_VarType_Aux = T_VarType;
    calc_nextciclevar(handles, [0, 0]);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_selectvar2_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to popupmenu_selectvar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_asign.
function pushbutton_asign_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_asign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Num_assignation;

if strcmp(get(hObject, 'string'), 'Calcular')    
    calc_nextmultiplier(handles);    
elseif strcmp(get(hObject, 'string'), 'Asignar')
    calc_nextassignment(handles);
elseif strcmp(get(hObject, 'string'), 'Buscar')
    calc_nextciclevar(handles, Num_assignation);
end

% --- Executes on button press in pushbutton_asignall.
function pushbutton_asignall_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to pushbutton_asignall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ----------------
function calc_nextmultiplier(handles)
global Solution Node_current Num_assignation T_Tableau Dimension ...
    T_VarType_Aux T_VarType Increment Matrix_problem;

if all(Node_current==0)    
    T_Tableau = zeros(Dimension(1), Dimension(2));
    T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);    
    set_environment('next_calc', handles);    
    var = get(handles.popupmenu_selectvar2, 'string'); % se recupera la variable no básica seleccionada
    V = var(get(handles.popupmenu_selectvar2, 'value'), 1); % se recupera variable seleccionada 
    I = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 2)); % se recupera el índice de variable seleccionada 
    if strcmp(V, 'U')        
        T_Tableau(I, end) = 0; % Asignacion inicial
        if I == 1            
            T_Tableau(end, Solution(1, 1, 2)) = T_Tableau(1, Solution(1, 1, 2))-T_Tableau(1, end);
            T_Tableau(1, Solution(1, 1, 2)) = 0;
            Increment = 1;
            Num_assignation = 1;
            Node_current = [1, Solution(1, 1, 2)];           
        else
            Num_assignation = Dimension(1)+Dimension(2)-3;
            T_Tableau(end, Solution(Num_assignation, 1, 2)) = T_Tableau(Dimension(1)-1, ...
                Solution(Num_assignation, 1, 2))-T_Tableau(Dimension(1)-1, end);
            T_Tableau(Dimension(1)-1, Solution(Num_assignation, 1, 2)) = 0;
            Increment = -1;            
            Node_current = [Dimension(1)-1, Solution(Num_assignation, 1, 2)];           
        end
    elseif strcmp(V, 'V')        
        T_Tableau(end, I) = 0; % Asignacion inicial
        if I == Solution(1, 1, 2)            
            T_Tableau(1, end) = T_Tableau(1, Solution(1, 1, 2))-T_Tableau(end, Solution(1, 1, 2));
            T_Tableau(1, Solution(1, 1, 2)) = 0;
            Increment = 1;
            Num_assignation = 1;
            Node_current = [1, Solution(1, 1, 2)];           
        else
            Num_assignation = Dimension(1)+Dimension(2)-3;
            T_Tableau(Solution(Num_assignation, 1, 1), end) = T_Tableau(Solution(Num_assignation, 1, 1), ...
                Solution(Num_assignation, 1, 2))-T_Tableau(end, Solution(Num_assignation, 1, 2));
            T_Tableau(Solution(Num_assignation, 1, 1), Solution(Num_assignation, 1, 2)) = 0;
            Increment = -1;            
            Node_current = [Dimension(1)-1, Solution(Num_assignation, 1, 2)];           
        end 
    end
else
    Num_assignation = Num_assignation + Increment;  
    if Node_current(1) == Solution(Num_assignation, 1, 1)       
        T_Tableau(end, Solution(Num_assignation, 1, 2)) = T_Tableau(Solution(Num_assignation, 1, 1), ...
            Solution(Num_assignation, 1, 2))-T_Tableau(Solution(Num_assignation, 1, 1), end);
        T_Tableau(Solution(Num_assignation, 1, 1),Solution(Num_assignation, 1, 2)) = 0;  
    else
        T_Tableau(Solution(Num_assignation, 1, 1), end) = T_Tableau(Solution(Num_assignation, 1, 1), ...
            Solution(Num_assignation, 1, 2))-T_Tableau(end, Solution(Num_assignation, 1, 2));
        T_Tableau(Solution(Num_assignation, 1, 1),Solution(Num_assignation, 1, 2)) = 0;
    end
    Node_current = [Solution(Num_assignation, 1, 1), Solution(Num_assignation, 1, 2)];           
end
if (Num_assignation == 1 && Increment == -1) || (Num_assignation == Dimension(1)+Dimension(2)-3 && Increment == 1)
    for i = 1:Dimension(1)-1
        for j = 1:Dimension(2)-1
            if T_VarType(i,j) ~= 1
                T_Tableau(i,j) = T_Tableau(i,j)-T_Tableau(end,j)-T_Tableau(i,end);
            end
        end
    end    
    % se construye el arreglo de cadenas de las variables    
    count = sum(sum(T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) < 0));
    if count ~= 0 % en el caso que haya Rj negativos
        set_environment('next_cicle', handles);
        Var = cell(count, 1);
        minimo = 0; indice = 0;
        i = 1;
        for j = 1:Dimension(1)-1
            for k =1:Dimension(2)-1
                if T_VarType(j, k) == 0
                    if T_Tableau(j, k) < 0
                        if minimo > T_Tableau(j, k)
                            indice = count;
                            minimo = T_Tableau(j, k);
                            Node_current = [j, k];
                        end
                        Var(i) = cellstr(strcat('X',strcat(num2str(j),num2str(k))));
                        i = i + 1; 
                    end
                end
            end
        end

        Num_assignation = 0;
        T_VarType_Aux = T_VarType;
        set(handles.popupmenu_selectvar2, 'string', char(Var));
        set(handles.popupmenu_selectvar2, 'value', indice);  
    else
        set(handles.popupmenu_selectvar2, 'Enable', 'off');
        set(handles.popupmenu_selectvar2, 'value', 1);
        set(handles.popupmenu_selectvar2, 'string', char('', ''));
        set_environment('end', handles);
    end
end
All_display = T_Tableau;
set(handles.table_simplexdisplay, 'data', All_display);


% ------------------------
function calc_nextciclevar(handles, Num_cassignation)
global Node_current Num_assignation T_Tableau T_VarType_Aux Node_NBV Dimension Node_Cicle...
    Empty_dimension Matrix_problem Solution_initial Solution Solution_change Minimo;

if Num_assignation == 0
    Node_Cicle = zeros(Dimension(1)+Dimension(2)-2, 1, 3);    
    Solution_change = zeros(1,Dimension(1)+Dimension(2)-3);
    T_Tableau = zeros(Dimension(1), Dimension(2));
    T_Tableau(end, 1:Dimension(2)) = Matrix_problem(end, 1:Dimension(2));
    T_Tableau(1:Dimension(1), end) = Matrix_problem(1:Dimension(1), end);
    var = get(handles.popupmenu_selectvar2, 'string'); % se recupera la variable no básica seleccionada
    I = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 2)); % se recupera variable seleccionada 
    J = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 3)); % se recupera el índice de variable seleccionada 
    Node_current = [I, J];
    T_Tableau(Node_current(1), Node_current(2)) = 1;
    Node_NBV = Node_current;  
    ind1 = find(T_VarType_Aux(Node_current(1), :));
    dim = size(ind1);
    if ~isempty(ind1)
        if Node_current(2) < ind1(1)
            Linf = 1;
            Lsup = dim(2);
            change = 1;
        else
            Linf = dim(2);
            Lsup = 1;
            change = -1;
        end
        for i=Linf:change:Lsup
            Num_assignation = Num_assignation + 1;                    
            Node_current = [Node_current(1), ind1(i)];
            Node_Cicle(Num_assignation, 1, 1) = Node_current(1); Node_Cicle(Num_assignation, 1, 2) = Node_current(2);
            T_VarType_Aux(Node_current(1), ind1(i)) = 0;
            ind2 = find(T_VarType_Aux(:, ind1(i)),1);
            if ~isempty(ind2)
                T_Tableau(Node_current(1), ind1(i)) = -1;
                Node_Cicle(Num_assignation, 1, 3) = -1;
                for l = 1:Dimension(1)+Dimension(2)-3;
                    if Solution(l, 1, 1) == Node_current(1) && Solution(l, 1, 2) == ind1(i)
                        Solution_change(l) = -1;
                        break;
                    end
                end                
                break;
            else
                T_Tableau(Node_current(1), ind1(i)) = 0;                
            end            
        end
        Empty_dimension = 0;
    else
        ind1 = find(T_VarType_Aux(:, Node_current(2)));
        dim = size(ind1);
        if ~isempty(ind1)
            if Node_current(1) < ind1(1)
                Linf = 1;
                Lsup = dim(2);
                change = 1;
            else
                Linf = dim(2);
                Lsup = 1;
                change = -1;
            end
            for i=Linf:change:Lsup
                Num_assignation = Num_assignation + 1;
                Node_current = [ind1(i), Node_current(2)];
                Node_Cicle(Num_assignation, 1, 1) = Node_current(1); Node_Cicle(Num_assignation, 1, 2) = Node_current(2);
                T_VarType_Aux(ind1(i), Node_current(2)) = 0;
                ind2 = find(T_VarType_Aux(ind1(i), :),1);
                if ~isempty(ind2)
                    T_Tableau(ind1(i), Node_current(2)) = -1;
                    Node_Cicle(Num_assignation, 1, 3) = -1;
                    for l = 1:Dimension(1)+Dimension(2)-3;
                        if Solution(l, 1, 1) == ind1(i) && Solution(l, 1, 2) == Node_current(2)
                            Solution_change(l) = -1;
                            break;
                        end
                    end
                    break;
                else                    
                    T_Tableau(ind1(i), Node_current(2)) = 0;                    
                end            
            end
            Empty_dimension = 1;
        end
    end        
else
    Node_current = [Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)];
    if Empty_dimension == 1
        ind1 = find(T_VarType_Aux(Node_current(1), :));
        dim = size(ind1);
        if ~isempty(ind1)
            if Node_current(2) < ind1(1)
                Linf = 1;
                Lsup = dim(2);  
                change = 1;
            else
                Linf = dim(2);
                Lsup = 1;
                change = -1;
            end
            for i = Linf:change:Lsup                
                Num_assignation = Num_assignation + 1;
                Node_current = [Node_current(1), ind1(i)];
                Node_Cicle(Num_assignation, 1, 1) = Node_current(1); Node_Cicle(Num_assignation, 1, 2) = Node_current(2);
                T_VarType_Aux(Node_current(1), ind1(i)) = 0;
                ind2 = find(T_VarType_Aux(:, ind1(i)),1);
                if ~isempty(ind2)
                    if Node_Cicle(Num_cassignation, 1, 3) == -1
                        T_Tableau(Node_current(1), ind1(i)) = 1;
                        Node_Cicle(Num_assignation, 1, 3) = 1;
                        for l = 1:Dimension(1)+Dimension(2)-3;
                            if Solution(l, 1, 1) == Node_current(1) && Solution(l, 1, 2) == ind1(i)
                                Solution_change(l) = 1;
                                break;
                            end
                        end
                    else
                        T_Tableau(Node_current(1), ind1(i)) = -1;
                        Node_Cicle(Num_assignation, 1, 3) = -1;
                        for l = 1:Dimension(1)+Dimension(2)-3;
                            if Solution(l, 1, 1) == Node_current(1) && Solution(l, 1, 2) == ind1(i)
                                Solution_change(l) = -1;
                                break;
                            end
                        end
                    end                    
                    if (Node_NBV(1) == Node_current(1) || Node_NBV(2) == Node_current(2)) && Num_cassignation > 2
                        Num_assignation = Dimension(1)+Dimension(2)-3;
                    end
                    break;
                else
                    if (Node_NBV(1) == Node_current(1) || Node_NBV(2) == Node_current(2)) && Num_cassignation > 2
                       if Node_Cicle(Num_cassignation, 1, 3) == -1
                            T_Tableau(Node_current(1), ind1(i)) = 1;
                            Node_Cicle(Num_assignation, 1, 3) = 1;
                            for l = 1:Dimension(1)+Dimension(2)-3;
                                if Solution(l, 1, 1) == Node_current(1) && Solution(l, 1, 2) == ind1(i)
                                    Solution_change(l) = 1;
                                    break;
                                end
                            end
                        else
                            T_Tableau(Node_current(1), ind1(i)) = -1;
                            Node_Cicle(Num_assignation, 1, 3) = -1;
                            for l = 1:Dimension(1)+Dimension(2)-3;
                                if Solution(l, 1, 1) == Node_current(1) && Solution(l, 1, 2) == ind1(i)
                                    Solution_change(l) = -1;
                                    break;
                                end
                            end
                       end                       
                       Num_assignation = Dimension(1)+Dimension(2)-3;                                                
                    else
                        T_Tableau(Node_current(1), ind1(i)) = 0;              
                    end                    
                end
            end
            Empty_dimension = 0; 
        else        
            Empty_dimension = 1;
            for k = (Dimension(1)+Dimension(2)-3):-1:1;
                if Node_Cicle(k, 1, 3) == -1 || Node_Cicle(k, 1, 3) == 1
                    Num_cassignation = k;
                    break;
                end
            end
            T_Tableau(Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)) = 0;
            Node_Cicle(Num_cassignation, 1, 3) = 0;
            for l = 1:Dimension(1)+Dimension(2)-3;
                if Solution(l, 1, 1) == Node_Cicle(Num_cassignation, 1, 1) && Solution(l, 1, 2) == Node_Cicle(Num_cassignation, 1, 2)
                    Solution_change(l) = 0;
                    break;
                end
            end  
            calc_nextciclevar(handles,Num_cassignation);
            return;
        end
    else
        ind1 = find(T_VarType_Aux(:, Node_current(2)));
        dim = size(ind1);
        if ~isempty(ind1)
            if Node_current(1) < ind1(1)
                Linf = 1;
                Lsup = dim(2);
                change = 1;
            else
                Linf = dim(2);
                Lsup = 1;
                change = -1;
            end
            for i=Linf:change:Lsup
                Num_assignation = Num_assignation + 1;
                Node_current = [ind1(i), Node_current(2)];
                Node_Cicle(Num_assignation, 1, 1) = Node_current(1); Node_Cicle(Num_assignation, 1, 2) = Node_current(2);
                T_VarType_Aux(ind1(i), Node_current(2)) = 0;
                ind2 = find(T_VarType_Aux(ind1(i), :),1);
                if ~isempty(ind2)                    
                    if Node_Cicle(Num_cassignation, 1, 3) == -1
                        T_Tableau(ind1(i), Node_current(2)) = 1;
                        Node_Cicle(Num_assignation, 1, 3) = 1;
                        for l = 1:Dimension(1)+Dimension(2)-3;
                            if Solution(l, 1, 1) == ind1(i) && Solution(l, 1, 2) == Node_current(2)
                                Solution_change(l) = 1;
                                break;
                            end
                        end
                    else
                        T_Tableau(ind1(i), Node_current(2)) = -1;
                        Node_Cicle(Num_assignation, 1, 3) = -1;
                        for l = 1:Dimension(1)+Dimension(2)-3;
                            if Solution(l, 1, 1) == ind1(i) && Solution(l, 1, 2) == Node_current(2)
                                Solution_change(l) = -1;
                                break;
                            end
                        end
                    end
                    if (Node_NBV(1) == Node_current(1) || Node_NBV(2) == Node_current(2)) && Num_cassignation > 1
                        Num_assignation = Dimension(1)+Dimension(2)-3;
                    end                    
                    break;
                else
                    if (Node_NBV(1) == Node_current(1) || Node_NBV(2) == Node_current(2)) && Num_cassignation > 1
                        if Node_Cicle(Num_cassignation, 1, 3) == -1
                            T_Tableau(ind1(i), Node_current(2)) = 1;
                            Node_Cicle(Num_assignation, 1, 3) = 1;
                            for l = 1:Dimension(1)+Dimension(2)-3;
                                if Solution(l, 1, 1) == ind1(i) && Solution(l, 1, 2) == Node_current(2)
                                    Solution_change(l) = 1;
                                    break;
                                end
                            end
                        else
                            T_Tableau(ind1(i), Node_current(2)) = -1;
                            Node_Cicle(Num_assignation, 1, 3) = -1;
                            for l = 1:Dimension(1)+Dimension(2)-3;
                                if Solution(l, 1, 1) == ind1(i) && Solution(l, 1, 2) == Node_current(2)
                                    Solution_change(l) = -1;
                                    break;
                                end
                            end                            
                       end
                       if Num_assignation == Num_cassignation
                            Num_assignation = Dimension(1)+Dimension(2)-3;
                       end                
                    else
                        T_Tableau(ind1(i), Node_current(2)) = 0;
                    end
                end                                
            end
            Empty_dimension = 1;   
        else
           Empty_dimension = 0;
           for k = (Dimension(1)+Dimension(2)-3):-1:1;
                if Node_Cicle(k, 1, 3) == -1 || Node_Cicle(k, 1, 3) == 1
                    Num_cassignation = k;
                    break;
                end
            end
            T_Tableau(Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)) = 0;
            Node_Cicle(Num_cassignation, 1, 3) = 0;
            for l = 1:Dimension(1)+Dimension(2)-3;
                if Solution(l, 1, 1) == Node_Cicle(Num_cassignation, 1, 1) && Solution(l, 1, 2) == Node_Cicle(Num_cassignation, 1, 2)
                    Solution_change(l) = 0;
                    break;
                end
            end  
            calc_nextciclevar(handles,Num_cassignation);            
            return;
        end
    end
end

if Num_assignation == Dimension(1)+Dimension(2)-3    
    Solution_initial = 0;
    Node_current = [0, 0];
    set_environment('next_assign', handles);
    set(handles.text_selectvar2, 'string', 'Seleccionar variable que sale');
    Solution_Aux = Solution;
    Solution_Aux(Solution_change >= 0, :, 3) =Inf; 
    minimo = min(Solution_Aux(:, :, 3));
    Minimo = minimo;
    ind = find(Solution_Aux(:, :, 3) == minimo);
    dim = size(ind);
    Var = cell(dim(1), 1);
    for i = 1:dim(1)
        Var(i) = cellstr(strcat('X',strcat(num2str(Solution(ind(i), 1, 1)), num2str(Solution(ind(i), 1, 2)))));
    end
    set(handles.popupmenu_selectvar2, 'string', char(Var));
    set(handles.popupmenu_selectvar2, 'value', 1);      
end
All_display = T_Tableau;
set(handles.table_simplexdisplay, 'data', All_display);


% --- Executes on button press in pushbutton_watchsolution.
function pushbutton_watchsolution_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_watchsolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Solution Matrix_problem Dimension ;

ObjectiveValue = 0;
All_display = zeros(Dimension(1),Dimension(2));
for j=1:Dimension(1)+Dimension(2)-3
    All_display(Solution(j, 1, 1), Solution(j, 1, 2)) = Solution(j, 1, 3);
    ObjectiveValue = ObjectiveValue + Solution(j, 1, 3)*Matrix_problem(Solution(j, 1, 1), Solution(j, 1, 2));
end 
All_display(Dimension(1), Dimension(2)) = ObjectiveValue;
set(handles.table_simplexdisplay, 'data', All_display);
