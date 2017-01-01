% Copyright 2011 2012 2016 M.Sc. Porfirio Armando Rodr�guez

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

handles.Method = 0;
handles.latex = '';
handles.latexproblem = '';
handles.latexIIphasesproblem = '';
handles.latexfile = '';

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

% --- Funci�n utilitaria: configura los distintos ambientes de ejecuci�n
function set_environment(environmentname, handles)
% environmentname especifica el tipo de ambiente por configurar
% handles         estructura con manejadores y datos de usuario
global Variables_ind;

if strcmp(environmentname, 'next')
    Variables_ind = 0;
    set(handles.Sensibility, 'enable', 'off');
    set(handles.Postoptimality, 'enable', 'off');
    set(handles.slider_increment, 'value', 0);
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.popupmenu_selectvar3, 'Enable', 'on');
    val_switches = ['on ';'on ';'on ';'on ';'on ';'off'];   
elseif strcmp(environmentname, 'sol_multiples')
    set(handles.Sensibility, 'enable', 'on');
    set(handles.Postoptimality, 'enable', 'on');
    val_switches = ['on ';'off';'on ';'off';'on ';'on '];
    set(handles.popupmenu_selectvar3, 'Enable', 'on');
    set(handles.Saveproblem, 'enable', 'on');
elseif strcmp(environmentname, 'end')
    if handles.Method ~= 3
        set(handles.Sensibility, 'enable', 'on');
        set(handles.Postoptimality, 'enable', 'on');
    end
    set(handles.slider_increment, 'value', 0);
    val_switches = ['off';'off';'on ';'off';'on '; 'off'];
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.Next_value, 'Enable', 'off');
    set(handles.pushbutton_asignall, 'Enable', 'off');
    set(handles.Watch_varval, 'Enable', 'on'); 
    set(handles.popupmenu_selectvar, 'Enable', 'off');
    set(handles.popupmenu_selectvar, 'string', char(' ', ' '));
    set(handles.popupmenu_selectvar2, 'Enable', 'off');
    set(handles.popupmenu_selectvar2, 'string', char(' ', ' '));
    set(handles.popupmenu_selectvar3, 'Enable', 'off');
    set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
elseif strcmp(environmentname, 'next_assign')
    set(handles.Sensibility, 'enable', 'off');
    set(handles.Postoptimality, 'enable', 'off');
    setTableauTags(handles, char('Origen', 'Destino', 'Demanda', 'Oferta'));
    %set(handles.popupmenu_selectvar2, 'enable', 'on');
    %set(handles.popupmenu_selectvar3, 'enable', 'off');
    set(handles.Next_value, 'Label', 'Siguiente valor');
    set(handles.pushbutton_asignall, 'string', 'Asignar todo');
    set(handles.Next_value, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on');
    set(handles.Watch_varval, 'Enable', 'off');  
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['off';'off';'on ';'off';'on ';'off'];   
elseif strcmp(environmentname, 'next_calc')    
    setTableauTags(handles, char('Origen', 'Destino', 'V', 'U'));
    set(handles.Next_value, 'Label', 'Siguiente multiplicador');
    set(handles.pushbutton_asignall, 'string', 'Calcular todo');   
    set(handles.popupmenu_selectvar2, 'Enable', 'on');       
    set(handles.popupmenu_selectvar, 'Enable', 'off');
    set(handles.popupmenu_selectvar, 'string', char(' ', ' '));    
    set(handles.popupmenu_selectvar3, 'Enable', 'off');
    set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
    set(handles.Next_value, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on');
    set(handles.Watch_varval, 'Enable', 'on');
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['off';'off';'on ';'off';'on ';'off']; 
elseif strcmp(environmentname, 'next_cicle')
    setTableauTags(handles, char('Origen', 'Destino', 'Demanda', 'Oferta'));
    set(handles.Next_value, 'Label', 'Siguiente nodo');
    set(handles.pushbutton_asignall, 'string', 'Buscar todo');
    set(handles.popupmenu_selectvar3, 'Enable', 'off');
    set(handles.popupmenu_selectvar, 'Enable', 'on');
    set(handles.Next_value, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on'); 
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['on ';'off';'on ';'off';'on ';'off'];   
end

gui_environment(val_switches, handles);

% --- Funci�n utilitaria: configura el ambiente en cada paso del M�todo
% Simplex
function gui_environment(val_switches, handles)
% val_switches    valores para la propiedad de cada control por configurar
% handles         estructura con manejadores y datos de usuario
set(handles.popupmenu_selectvar, 'Enable', val_switches(1,:));
set(handles.slider_increment, 'Enable', val_switches(2,:));
set(handles.pushbutton_start, 'Enable', val_switches(3,:));
set(handles.pushbutton_next, 'Enable', val_switches(4,:));
set(handles.Modify, 'enable', val_switches(5,:));
set(handles.Multiple_solution, 'Enable', val_switches(6,:));

% --- Executes on button press in pushbutton_start.
function pushbutton_start_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem Order_initial latex;
% Matrix_problem    especificaci�n del problema
% Order_initial     orden inicial de vectores columna de la tabla Simplex, el orden de la
% matrix identidad primero y luego los dem�s vectores

%AGREGADO(27/12/2016)

if handles.istwophases == 1 && handles.isinit_secondphase ~= 1
    handles.latex = handles.latexIIphasesproblem;
    handles.whatphase = 1;
    Matrix_problem = handles.First_Matrix_problem;
    handles.Order_initial = handles.First_Order_initial;    
elseif handles.istwophases == 1 && handles.isinit_secondphase == 1
    handles.isinit_secondphase = 0;
    Matrix_problem = handles.gui_Matrix_problem;
    generar_latexspecification();
    handles.latex = [handles.latex, latex];
    Order_initial = handles.Order_initial;
    generar_latexsimplexbegin()
    handles.latex = [handles.latex, latex];
    guidata(handles.output, handles);
else
    %set(handles.Next_value, 'Enable', 'on');
    handles.latex = handles.latexproblem;
    handles.Order_initial = Order_initial;
end

guidata(handles.output, handles);
%AGREGADO(27/12/2016)

% Configura la tabla inicial del M�todo Simplex
% handles.Order_initial = Order_initial;
setProblem(Matrix_problem, handles)


% --- Funci�n utilitaria: obtiene las variables no b�sicas candidatas a
% ser seleccionadas para convertirse en b�sicas seg�n el criterio de
% mejoramiento de la soluci�n (coeficientes de costo relativo negativos)
function h = calc_variables(handles)
% handles       estructura con manejadores y datos de usuario
global Dimension Tableau Order_current Solution_initial Matrix_problem latex; %#ok<NUSED>
% Dimension     dimensiones de la especificaci�n del problema en t�rminos
% del n�mero de ecuaciones y variables
% Tableau       Tabla del Simplex actual
% Order_current orden actual de vectores columna de la tabla Simplex

if handles.Method == 1
    Rj = Tableau(end,1:(Dimension(2)-1)); % Recupera los coeficientes de costo reducido (Rj) de la tabla Simplex
    [Y, I] = sort(Order_current); %#ok<ASGLU> % obtiene los indices de variables ordenados de menor a mayor
    IRjneg = I(Order_current(Rj < 0));  % obtiene los indices de variables con Rj negativo
    IVar = IRjneg;
elseif handles.Method == 2
    X0 = Tableau(1:(Dimension(1)-1), end); % Recupera los valores de las variables b�sicas de la tabla Simplex
    [Y, I] = sort(Order_current); %#ok<ASGLU> % obtiene los indices de variables ordenados de menor a mayor
    IX0neg = I(X0 < 0);  % obtiene los indices de variables con valores negativos
    IVar =  IX0neg;
end

if handles.Method ~= 3
    % se construye el arreglo de cadenas de las variables
    dim2 = size(IVar);
    Var = cell(dim2(2), 1);
    for i = 1:dim2(2)
        Var(i) = cellstr(strcat('X',num2str(IVar(i))));
    end

    if ~isempty(Var) % en el caso que haya Rj negativos
        if handles.Method == 1
            [Y, p] = min(Rj(IRjneg)); %#ok<ASGLU> % obtiene el �ndice de variable del Rj m�s negativo
        elseif handles.Method == 2
            [Y, p] = min(IX0neg); %#ok<ASGLU> % obtiene el �ndice de variable del Rj m�s negativo
        end
        set(handles.popupmenu_selectvar, 'string', char(Var)); % despliega las variables no b�sicas en la interfaz
        set_environment('next', handles); % ajusta el ambiente
        set(handles.popupmenu_selectvar, 'value', p); % se selecciona, por omisi�n, la variable m�s negativa
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
        
        %MODIFICADO(27/12/2016)            
        maximo = max(find(Order_current < Dimension(1))); %#ok<MXFND>
        if ~isempty(Var) && handles.Method == 1 && ...
                (((handles.istwophases == 1 && handles.whatphase == 1)  && ...
                    maximo > Dimension(2)-handles.maxcanon_vector-2) || ...
                    handles.istwophases == 0 || (handles.istwophases == 1 && ... 
                    handles.whatphase == 2)) % en el caso que haya Rj = 0        
            set(handles.popupmenu_selectvar, 'string', char(Var));
            set_environment('sol_multiples', handles);
            set(handles.popupmenu_selectvar, 'value', 1);
            calc_ratios(handles);
            msgbox('El proceso ha encontrado una soluci�n �ptima. Haga clic en "Solucion m�ltiple" para encontrar otra.','M�todo simplex','modal');
            generar_latexanalisis();
        else % No hay Soluciones m�ltiples, no hay m�s opciones
            if (handles.istwophases == 1 && handles.whatphase == 1) 
                if Tableau(end, end) == 0 && maximo <= Dimension(2)-handles.maxcanon_vector-2
                    msgbox('La primera fase ha terminado. Haga clic en "Iniciar" para la segunda fase.','M�todo simplex de dos fases','modal');
                    generar_latexanalisis('\subsection{An�lisis de resultados}');
                    %Matrix_problem = zeros(size(handles.Orig_Matrix_problem));
                    Matrix_problem = Tableau(:, 1:Dimension(2)-handles.maxcanon_vector-2);
                    Matrix_problem(:, Dimension(2)-handles.maxcanon_vector-1) = Tableau(:,Dimension(2)-handles.maxcanon_vector);
                    Matrix_problem(end, :) = handles.Orig_Matrix_problem(end, :);
                    Order_current = Order_current(1:Dimension(2)-handles.maxcanon_vector-2);                    
                    Order_current(Dimension(2)-handles.maxcanon_vector-1) = Dimension(2)-handles.maxcanon_vector-1;
                    Order_current(Order_current > Dimension(1)-1) = Dimension(1):Dimension(2)-handles.maxcanon_vector-1;
                    handles.gui_Matrix_problem = Matrix_problem;
                    handles.Order_initial = Order_current;
                    handles.whatphase = 2;
                    handles.isinit_secondphase = 1;
                    %guidata(handles.output, handles);
                    %setProblem(Matrix_problem, handles)
                    %return;
                else
                    msgbox('La primera fase ha terminado. El problema original no tiene soluci�n','M�todo simplex de dos fases','modal');
                end
            elseif (handles.istwophases == 1 && handles.whatphase == 2)
                msgbox('La segunda fase ha terminado. Se ha encontrado una soluci�n �ptima','M�todo simplex de dos fases','modal');
                generar_latexanalisis('\subsection{An�lisis de resultados}'); 
            else                               
                generar_latexanalisis('\section{An�lisis de resultados}');              
                msgbox('Se ha encontrado una soluci�n �ptima','M�todo simplex','modal'); 
            end
            set(handles.popupmenu_selectvar, 'string', char(' ',' '));
            set(handles.popupmenu_selectvar3, 'string', char(' ',' '));
            set_environment('end', handles);
        end
        %MODIFICADO(27/12/2016)
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

h = handles;

% --- Funci�n utilitaria: calcula las razones del criterio de factibilidad
function calc_ratios(handles)
% handles   estructura con manejadores y datos de usuario
global Tableau Dimension Order_current;
% All_display condensa los datos que se desplegar�n en la interfaz

var = get(handles.popupmenu_selectvar, 'string'); % se recupera la variable no b�sica seleccionada
dim = size(var(get(handles.popupmenu_selectvar, 'value'), :));
J = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2:dim(2))); % se recupera el �ndice de variable seleccionada
if handles.Method == 1
    % Se calculan las razones para la variable no b�sica seleccionada
    Y0 = Tableau(1:(Dimension(1)-1), end); % se recupera el vector de t�rminos coeficientes
    Yj = Tableau(1:(Dimension(1)-1), J); % se recupera el vector columna asociada a la variable no b�sica
    ratios = Y0./Yj;
    ratios_aux = ratios;
    [Y, I] = sort(Order_current); %#ok<ASGLU> % obtiene los indices de variables ordenados de menor a mayor    
    ratios_aux(ratios < 0 | isnan(ratios) | Yj < 0) = Inf; 
    [C, p] = min(ratios_aux); %#ok<NASGU>
    if C ~= Inf
        ind = find(ratios_aux == C);
        dim = size(ind);
        % se construye el arreglo de cadenas de las variables    
        Var = cell(dim(1), 1);
        for i = 1:dim(1)
            Var(i) = cellstr(strcat('X',num2str(I(ind(i)))));
        end
        set(handles.popupmenu_selectvar3, 'string', char(Var));
        set(handles.popupmenu_selectvar3, 'value', 1);
    else
        set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
        set(handles.popupmenu_selectvar3, 'value', 1);
    end
    
    All_display = zeros(Dimension(1), Dimension(2)+1);
    % se actualiza la tabla de la interfaz
    All_display(:, 1:Dimension(2)) = Tableau;
    All_display(1:(Dimension(1)-1), end) = ratios;
elseif handles.Method == 2
    % Se calculan las razones para la variable b�sica seleccionada
    Fi = Tableau(Order_current(J), 1:(Dimension(2)-1)); % se recupera el vector fila asociada a la variable b�sica
    Rj = Tableau(end, 1:(Dimension(2)-1)); % se recupera el vector de coeficientes reducidos
    ratios = -Rj./Fi;
    ratios_aux = ratios;
    Basic_Var = Order_current < Dimension(1);
    ratios_aux(Basic_Var) = Inf;
    ratios_aux(ratios < 0 | isnan(ratios) | Fi > 0) = Inf; 
    [C, p] = min(ratios_aux); %#ok<NASGU>
    if C ~= Inf
        ind = find(ratios_aux == C);
        dim = size(ind);
        % se construye el arreglo de cadenas de las variables    
        Var = cell(dim(2), 1);
        for i = 1:dim(2)
            Var(i) = cellstr(strcat('X',num2str(ind(i))));
        end
        set(handles.popupmenu_selectvar3, 'string', char(Var));
        set(handles.popupmenu_selectvar3, 'value', 1);
    else
        set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
        set(handles.popupmenu_selectvar3, 'value', 1);
    end
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

if handles.Method == 1
    colFormat=cell(1,Dimension(2)+1);    
else
    colFormat=cell(1,Dimension(2));
end
colFormat(1,:) = cellstr('rat');
set(handles.table_simplexdisplay, 'columnformat', colFormat);
if Solution_initial == 1
    if all(Node_current==0)
        Num_assignation = 1;
        T_VarType = zeros(Dimension(1)-1, Dimension(2)-1);
        var = get(handles.popupmenu_selectvar2, 'string'); % se recupera la variable no b�sica seleccionada
        I = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 2)); % se recupera el �ndice de variable seleccionada 
        J = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 3)); % se recupera el �ndice de variable seleccionada 
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
    ObjectiveValue = 0;
    for j=1:Num_assignation        
        All_display(Solution(j, 1, 1), Solution(j, 1, 2)) = Solution(j, 1, 3);
        ObjectiveValue = ObjectiveValue + Solution(j, 1, 3)*Matrix_problem(Solution(j, 1, 1), Solution(j, 1, 2));
    end  
    All_display(end, :) = T_Tableau(end, :);
    All_display(:, end) = T_Tableau(:, end);   
    All_display(Dimension(1), Dimension(2)) = ObjectiveValue;    
    set(handles.table_simplexdisplay, 'data', All_display);
    
    if Num_assignation == Dimension(1)+Dimension(2)-3
        msgbox('Se ha encontrado una soluci�n inicial.','C�lculo de nueva soluci�n.','modal');
        generar_latexnextsolution('\subsection{Encontrando una soluci�n inicial por sustituci�n hacia atr�s}');
        Node_current = [0, 0];        
        set_environment('next_calc', handles);
        T_Tableau = zeros(Dimension(1), Dimension(2));
        T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);
        % se construye el arreglo de cadenas de las variables    
        Var = cell(2, 1);    
        Var(1) = cellstr(strcat('U',num2str(1)));   
        Var(2) = cellstr(strcat('V',num2str(Dimension(2)-1)));           
        set(handles.popupmenu_selectvar2, 'string', char(Var));
        set(handles.popupmenu_selectvar2, 'value', 2);
    end
    set(handles.popupmenu_selectvar2, 'enable', 'on');
else
    var = get(handles.popupmenu_selectvar3, 'string'); % se recupera la variable no b�sica seleccionada
    I = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 2)); % se recupera variable seleccionada 
    J = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 3)); % se recupera el �ndice de variable seleccionada
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
    ObjectiveValue = 0;
    for j=1:Num_assignation
        All_display(Solution(j, 1, 1), Solution(j, 1, 2)) = Solution(j, 1, 3);
        ObjectiveValue = ObjectiveValue + Solution(j, 1, 3)*Matrix_problem(Solution(j, 1, 1), Solution(j, 1, 2));
    end  
    All_display(end, :) = T_Tableau(end, :);
    All_display(:, end) = T_Tableau(:, end);
    All_display(Dimension(1), Dimension(2)) = ObjectiveValue;
    set(handles.table_simplexdisplay, 'data', All_display);
    if Num_assignation == Dimension(1)+Dimension(2)-3
        msgbox('Se ha encontrado una nueva soluci�n.','C�lculo de nueva soluci�n.','modal');
        generar_latexnextsolution('\subsection{Encontrando una nueva soluci�n por redistribuci�n}');
        Node_current = [0, 0]; 
        set_environment('next_calc', handles);
        T_Tableau = zeros(Dimension(1), Dimension(2));
        T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);
        % se construye el arreglo de cadenas de las variables    
        Var = cell(2, 1);    
        Var(1) = cellstr(strcat('U',num2str(1)));   
        Var(2) = cellstr(strcat('V',num2str(Dimension(2)-1)));           
        set(handles.popupmenu_selectvar2, 'string', char(Var));
        set(handles.popupmenu_selectvar2, 'value', 2);
        %set(handles.popupmenu_selectvar3, 'string', char(' ',' '));
    elseif Num_assignation > 1
        set(handles.popupmenu_selectvar3, 'enable', 'off');
    end
end

% --- Executes on button press in pushbutton_next.
function pushbutton_next_Callback(hObject, eventdata, handles) %#ok<INUSL>
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global latex;

latex = '';
if handles.Method == 1   
    % se ejecuta el m�todo simplex
    simplex_primal(handles);
elseif handles.Method == 2  
    simplex_dual(handles);
end

setrowheaders(handles, char('X', 'X', 'Rj', 'Yi0'));
handles = calc_variables(handles);

handles.latex = [handles.latex, latex];
guidata(handles.output, handles);
    
% --- Se ejecuta el m�todo simplex primal
function simplex_primal(handles)
global Tableau Order_current operations;

% recupera el �ndice de la variable no b�sica seleccionada asociada a la
% columna que ingresar� a la base
var = get(handles.popupmenu_selectvar, 'string');
dim = size(var(get(handles.popupmenu_selectvar, 'value'), :));
q = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2:dim(2)));
% recupera el �ndice de la fila (o componente) del vector columna que
% dejar� la base seg�n la raz�n m�nima tal que sea positiva
var = get(handles.popupmenu_selectvar3, 'string');
dim = size(var(get(handles.popupmenu_selectvar3, 'value'), :));
p_aux = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 2:dim(2)));

if ~isnan(p_aux) % si existe alg�n yij que cumple el criterio de factibilidad 
    p = Order_current(p_aux);
    pivote_process(p, q);
    %AGREGADO(27/12/2016)
    set(handles.listbox_operations, 'string', operations);
    generar_latexsimplexnext()
    %AGREGADO(27/12/201
else
    msgbox('El conjunto representado por las restricciones no es acotado.', 'El valor objetivo decrece sin limite','modal');
    return;
end
% se actualiza la tabla de la interfaz
All_display = Tableau;
set(handles.table_simplexdisplay, 'data', All_display);

% --- Se ejecuta el m�todo simplex primal
function simplex_dual(handles)
global Tableau Order_current operations;

% recupera el �ndice de la variable no b�sica seleccionada asociada a la
% columna que ingresar� a la base
var = get(handles.popupmenu_selectvar, 'string');
dim = size(var(get(handles.popupmenu_selectvar, 'value'), :));
p_aux = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2:dim(2)));
p = Order_current(p_aux);

% recupera el �ndice de la fila (o componente) del vector columna que
% dejar� la base seg�n la raz�n m�nima tal que sea positiva
var = get(handles.popupmenu_selectvar3, 'string');
dim = size(var(get(handles.popupmenu_selectvar3, 'value'), :));
q = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 2:dim(2)));

if ~isnan(q) % si existe alg�n yij que cumple el criterio de factibilidad
    pivote_process(p, q);
    %AGREGADO(27/12/2016)
    set(handles.listbox_operations, 'string', operations);
    generar_latexsimplexnext()
    %AGREGADO(27/12/2016)
else
    msgbox('El conjunto representado por las restricciones es vac�o.', 'No hay soluci�n','modal');
    return;
end

% se actualiza la tabla de la interfaz
All_display = Tableau;
set(handles.table_simplexdisplay, 'data', All_display);


% --- Funci�n utilitaria: proceso del pivote sobre el elemento Ypq
function pivote_process(p, q)
global Tableau Order_current Dimension operations latex;

% Ecuaciones de pivote
Ypq = Tableau(p, q);
fp = Tableau(p, :); % recupera la fila p actual
fp = fp/Ypq;

%AGREGADO(27/12/2016)
[num, den] = numden(sym(1/Ypq, 'r'));
fraction = double([num, den]);
pivote_operation = sprintf('Se efectuaron las siguientes operaciones de pivoteo:\r');
operations = pivote_operation;

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, 'El elemento (', num2str(p), ',', num2str(q), ') de la matriz es el pivote: $\bf{', num2str(Ypq), '}$.']; 
latex = sprintf('%s \r', latex);
latex = [latex, 'Las operaciones de pivoteo son: \\ \\']; 
latex = sprintf('%s \r', latex);

if fraction(2) == 1
    pivote_operation = sprintf('f%d <= (%d) x f%d', p, fraction(1), p);
    latex = [latex, '$f_', num2str(p), '\leftarrow ', '(', num2str(fraction(1)), ')f_', num2str(p), '$\\'];     
else
    pivote_operation = sprintf('f%d <= (%d/%d) x f%d', p, fraction(1), fraction(2), p);
    latex = [latex, '$f_', num2str(p), '\leftarrow ', '(\frac{', num2str(fraction(1)), '}{', num2str(fraction(2)), '})f_', num2str(p), '$\\']; 
end
latex = sprintf('%s \r', latex);
operations = char(operations, pivote_operation);
%AGREGADO(27/12/2016)

Tableau(p, :) = fp; % se actualiza la fila p en la tabla del Simplex
% se actualizan las dem�s filas
for i = 1:Dimension(1)
    if i ~= p
        fi = Tableau(i, :);
        Yiq = Tableau(i, q);
        fi = fi - fp*Yiq;
        
        %AGREGADO(27/12/2016)
        [num, den] = numden(sym(Yiq/Ypq, 'r'));
        fraction = double([num, den]);
        if fraction(2) == 1
            pivote_operation = sprintf('f%d  <= f%d - (%d) x f%d', i, i, fraction(1), p);
            latex = [latex, '$f_', num2str(i), '\leftarrow f_', num2str(i), '-(', num2str(fraction(1)), ')f_', num2str(p), '$\\'];
        else
            pivote_operation = sprintf('f%d  <= f%d - (%d/%d) x f%d', i, i, fraction(1), fraction(2), p);
            latex = [latex, '$f_', num2str(i), '\leftarrow f_', num2str(i), '-(\frac{', num2str(fraction(1)), '}{', num2str(fraction(2)), '})f_', num2str(p), '$\\'];
        end
        latex = sprintf('%s \r', latex);
        operations = char(operations, pivote_operation);
        %AGREGADO(27/12/2016)
        
        Tableau(i, :) = fi; % fila p actualizada
    end
end

% intercambian roles las variables
Order_current(1, Order_current(1,:) == p) = Order_current(1,q);
Order_current(1, q) = p;
    

% --- Executes on selection change in popupmenu_selectvar.
function popupmenu_selectvar_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to popupmenu_selectvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Num_assignation T_VarType_Aux T_VarType;


if handles.Method == 3
    Num_assignation = 0;
    T_VarType_Aux = T_VarType;
    %calc_nextciclevar(handles, [0, 0]);
    calc_nextciclevar(handles, 0);
else
    calc_ratios(handles);
end
    

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

if handles.Method == 1
    ratios = All_display(1:(Dimension(1)-1), end);
    ratios_aux = ratios;
    ratios_aux(ratios < 0) = Inf; 
    [ratio, Y] = min(ratios_aux);  %#ok<NASGU>
    % se actualizan las variables b�sicas seg�n el incremento de la variable no
    % b�sica
    Xb = Aux_Tableau(:,end);
    Yj = Tableau(:, J);
    Aux_Tableau(:,end) = Xb - Yj*incr/100*ratio;
    % se actualiza la tabla de la interfaz
    All_display = zeros(Dimension(1), Dimension(2)+1);
    All_display(:, 1:Dimension(2)) = Aux_Tableau;
    All_display(1:(Dimension(1)-1), end) = ratios;
    set(handles.table_simplexdisplay, 'data', All_display);
elseif handles.Method == 2
    ratios = All_display(end, 1:(Dimension(2)-1));
    ratios_aux = ratios;    
    Vector_nonbasic = Order_current(1:Dimension(2)-1) >= Dimension(1);
    ratios_aux(ratios < 0 & Vector_nonbasic) = Inf;
    ratios_aux(ratios == 0 & ~Vector_nonbasic) = Inf; 
    [ratio, Y] = min(ratios_aux);  %#ok<NASGU>
    % se actualizan los coeficientes de costo reducido seg�n la
    % actualizaci�n de la soluci�n del dual
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


% --- actualiza la especificaci�n del problema
function setProblem(Problem, handles)
% Problem especificaci�n del problema
% handles estructura de manejadores y datos de usuario
global Dimension Matrix_problem Tableau Order_current Order_initial Node_current Solution T_Tableau Num_assignation;

% se actualizan algunas variables globales previamente descritas
Matrix_problem = Problem;
Dimension = size(Matrix_problem);
%AGREGADO(27/12/2016)
set(handles.Next_value, 'Enable', 'on');
set(handles.listbox_operations, 'string', '');
%AGREGADO(27/12/2016)

if handles.Method ~= 3
    Order_initial = handles.Order_initial;
    Order_current = Order_initial;
    Tableau = Matrix_problem; 

    % se convierte la �ltima en t�rminos de Rj =cj - zj
    [Y, I] = sort(Order_initial); %#ok<ASGLU>
    Tableau_sorted = Matrix_problem(:, I);
    c = Tableau_sorted(end,:);
    z = Tableau_sorted(end, 1:(Dimension(1)-1))*Tableau_sorted(1:(Dimension(1)-1),:);

    R = c - z;
    Tableau(end, I) = R;
    %AGREGADO(27/12/2016)
    set(handles.listbox_operations, 'string', ...
        char('La fila de costos unitarios se ha convertido en costos relativos:', 'r <= c - z'));
    %generar_latexsimplexbegin(handles);
    %AGREGADO(27/12/2016)
else
    Node_current = [0,0];
    Num_assignation = 0;
    Solution = zeros(Dimension(1)+Dimension(2)-3, 1, 3);
    T_Tableau = Matrix_problem;
end
% se actualiza la tabla de la interfaz
All_display = zeros(Dimension(1), Dimension(2));
if handles.Method ~= 3
    % se rotulan las columnas y las filas de manera correspondiente
    setTableauTags(handles, char('X', 'X', 'Rj', 'Yi0'));
    All_display(:, 1:Dimension(2)) = Tableau;
else
    % se rotulan las columnas y las filas de manera correspondiente
    setTableauTags(handles, char('Origen', 'Destino', 'Demanda', 'Oferta'));
    All_display(:, 1:Dimension(2)) = T_Tableau;
end
set(handles.table_simplexdisplay, 'data', All_display); 
% se calculan las nuevas variables no b�sicas si las hay
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
if handles.Method == 1
    colName(Dimension(2)+1) = cellstr('Yi0/Yij');    
end
set(handles.table_simplexdisplay, 'columnname', colName);

setrowheaders(handles, headers);

if handles.Method == 1
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
    if handles.Method ~= 3
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(I(i))));
    else
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end
end
rowName(Dimension(1)) = cellstr(headers(3,:));

if handles.Method == 2
    rowName(Dimension(1)+1) = cellstr('-Rj/Yij');   
end
set(handles.table_simplexdisplay, 'rowname', rowName);

% --- Executes when user attempts to close LPApp.
function LPApp_CloseRequestFcn(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to LPApp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% se cierran todas las ventanas abiertas hijas de la aplicaci�n principal

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

[filename, Y, Z] = uigetfile('*.mat', 'Selecciona un nombre de archivo v�lido'); %#ok<NASGU>
if filename ~= 0
    S = load([Y, filename]);
    
    %AGREGADO(27/12/2016)
    handles.istwophases = 0;
    handles.whatphase = 1;
    handles.isinit_secondphase = 0;
    handles.latex = '';
    if (isfield(S, 'Method'))
        if S.Method == 3            
            handles.Method = 3;                       
        elseif S.Method == 1            
            handles.Method = 1;
            set(handles.panel_enhancement, 'title', 'Panel de control (Simplex primal)');                      
        elseif S.Method == 2            
            handles.Method = 2;
            set(handles.panel_enhancement, 'title', 'Panel de control (Simplex dual)');                   
        end
        handles.setProblem = @setProblem;    
        handles.gui_Matrix_problem = S.Matrix_problem;
        handles.gui_Problem = Problem('LPApp', handles);
        guidata(handles.output, handles);
    else
        errordlg('El problema no es compatible con la versi�n actual', 'Problema no compatible');
    end
    %MODIFICADO(27/12/2016)
end

% --------------------------------------------------------------------
function Saveproblem_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Saveproblem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem;

%AGREGADO(27/12/2016)
if handles.Method == 3
    Method = 3; %#ok<NASGU>
elseif handles.Method == 1
    Method = 1; %#ok<NASGU>
    if handles.istwophases == 1
        Matrix_problem = handles.Orig_Matrix_problem;
    end        
elseif handles.Method == 2
    Method = 2; %#ok<NASGU>
end
%AGREGADO(27/12/2016)

uisave(cellstr(char('Matrix_problem', 'Method')), 'LProblem');
%uisave('Matrix_problem', 'LProblem');



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


% --- Executes on button press in pushbutton_asignall.
function pushbutton_asignall_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_asignall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Num_assignation Dimension latex;

if strcmp(get(handles.Next_value, 'Label'), 'Siguiente multiplicador')    
    while (Num_assignation ~= Dimension(1)+Dimension(2)-3)
        calc_nextmultiplier(handles);    
    end
    paso = 2;
    Num_assignation = 0;
elseif strcmp(get(handles.Next_value, 'Label'), 'Siguiente valor')
    while (Num_assignation ~= Dimension(1)+Dimension(2)-3)
        calc_nextassignment(handles);
    end
    paso = 1;
    Num_assignation = 0;
elseif strcmp(get(handles.Next_value, 'Label'), 'Siguiente nodo')
    while (Num_assignation ~= Dimension(1)+Dimension(2)-3)
        calc_nextciclevar(handles, Num_assignation);
    end
    paso = 3;
    Num_assignation = 0;
end
if paso == 1
    handles.latex = [handles.latex, latex];
elseif paso == 2
    handles.latex = [handles.latex, latex];
else
    generar_latexnextcicle();
    handles.latex = [handles.latex, latex];
end
guidata(handles.output, handles);


% ----------------
function calc_nextmultiplier(handles)
global Node_current Num_assignation T_Tableau Dimension ...
    T_VarType_Aux T_VarType By_dimension Matrix_problem latex;

colFormat(1,:) = cellstr('rat');
set(handles.table_simplexdisplay, 'columnformat', colFormat);
if all(Node_current==0)
    T_VarType_Aux = T_VarType;
     Num_assignation = 0;
    T_Tableau = zeros(Dimension(1), Dimension(2));
    T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);
    T_Tableau(end, :) = NaN;
    T_Tableau(:, end) = NaN;
    set_environment('next_calc', handles);    
    var = get(handles.popupmenu_selectvar2, 'string'); % se recupera la variable no b�sica seleccionada
    V = var(get(handles.popupmenu_selectvar2, 'value'), 1); % se recupera variable seleccionada 
    I = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 2)); % se recupera el �ndice de variable seleccionada 
    if strcmp(V, 'U')  
        By_dimension = 0;
        T_Tableau(I, end) = 0; % Asignacion inicial
        ind1 = find(T_VarType_Aux(I, :));
        dim1 = size(ind1);
        for i=1:dim1(2)
            Num_assignation = Num_assignation + 1;                    
            Node_current = [I, ind1(i)];
            T_VarType_Aux(I, ind1(i)) = 0;
            T_Tableau(end, Node_current(2)) = T_Tableau(I, Node_current(2))-T_Tableau(I, end);
            T_Tableau(I, Node_current(2)) = 0;         
        end
    elseif strcmp(V, 'V')
        By_dimension = 1;
        T_Tableau(end, I) = 0; % Asignacion inicial
        ind1 = find(T_VarType_Aux(:, I));
        dim1 = size(ind1);
        for i=1:dim1(1)
            Num_assignation = Num_assignation + 1;                    
            Node_current = [ind1(i), I];
            T_VarType_Aux(ind1(i), I) = 0;
            T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), I)-T_Tableau(end, I);
            T_Tableau(Node_current(1), I) = 0;
        end        
    end
else
    set(handles.popupmenu_selectvar2, 'enable', 'off');
    if By_dimension == 1
        for i = 1:Dimension(1)-1
            if ~isnan(T_Tableau(i, end))
                ind1 = find(T_VarType_Aux(i, :));
                dim1 = size(ind1);
                if ~isempty(ind1)                    
                    for j = 1:dim1(2)
                        Num_assignation = Num_assignation + 1;                    
                        Node_current = [i, ind1(j)];
                        T_VarType_Aux(i, ind1(j)) = 0;
                        T_Tableau(end, Node_current(2)) = T_Tableau(i, Node_current(2))-T_Tableau(i, end);
                        T_Tableau(i, Node_current(2)) = 0;
                    end
                end
            end
        end
        By_dimension = 0;
    else
        for i = 1:Dimension(2)-1
            if ~isnan(T_Tableau(end, i))
                ind1 = find(T_VarType_Aux(:,i));
                dim1 = size(ind1);
                if ~isempty(ind1)
                    for j = 1:dim1(1)
                        Num_assignation = Num_assignation + 1;                    
                        Node_current = [ind1(j), i];
                        T_VarType_Aux(ind1(j), i) = 0;
                        T_Tableau(Node_current(1), end) = T_Tableau(Node_current(1), i)-T_Tableau(end, i);
                        T_Tableau(Node_current(1), i) = 0;
                    end
                end
            end
        end
        By_dimension = 1;
    end         
end
if (Num_assignation == Dimension(1)+Dimension(2)-3)
    for i = 1:Dimension(1)-1
        for j = 1:Dimension(2)-1
            if T_VarType(i,j) ~= 1
                T_Tableau(i,j) = T_Tableau(i,j)-T_Tableau(end,j)-T_Tableau(i,end);
            end
        end
    end
    latex = '';
    generar_latexnextrj();
    % se construye el arreglo de cadenas de las variables    
    count = sum(sum(T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) < 0));
    if count ~= 0 % en el caso que haya Rj negativos
        msgbox('Hay costos reducidos negativos. La soluci�n actual se puede mejorar.','C�lculo de costos reducidos.','modal');
        set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
        set_environment('next_cicle', handles);
        Var = cell(count, 1);
        minimo = 0; indice = 0;
        i = 1;
        for j = 1:Dimension(1)-1
            for k =1:Dimension(2)-1
                if T_VarType(j, k) == 0
                    if T_Tableau(j, k) < 0
                        if minimo > T_Tableau(j, k)
                            indice = i;
                            minimo = T_Tableau(j, k);
                            Node_current = [j, k];
                        end
                        Var(i) = cellstr(strcat('X',strcat(num2str(j),num2str(k))));
                        i = i + 1; 
                    end
                end
            end
        end

        %Num_assignation = 0;
        T_VarType_Aux = T_VarType;
        set(handles.popupmenu_selectvar, 'string', char(Var));
        set(handles.popupmenu_selectvar, 'value', indice);  
    else
        msgbox('La soluci�n actual es �ptima. No hay costos reducidos negativos.','Soluci�n �ptima encontrada.','modal');
        generar_latexnextsolution('\section{An�lisis de Resultados}');
        set_environment('end', handles);
    end
end
All_display = T_Tableau;
set(handles.table_simplexdisplay, 'data', All_display);


% ------------------------
function calc_nextciclevar(handles, Num_cassignation)
global Node_current Num_assignation T_Tableau T_VarType_Aux Node_NBV Dimension Node_Cicle...
    Empty_dimension Matrix_problem Solution_initial Solution Solution_change Minimo;

if handles.Method == 1
    colFormat=cell(1,Dimension(2)+1);    
else
    colFormat=cell(1,Dimension(2));
end
colFormat(1,:) = cellstr('rat');
set(handles.table_simplexdisplay, 'columnformat', colFormat);
if Num_assignation == 0
    Node_Cicle = zeros(Dimension(1)+Dimension(2)-2, 1, 3);    
    Solution_change = zeros(1,Dimension(1)+Dimension(2)-3);
    T_Tableau = zeros(Dimension(1), Dimension(2));
    T_Tableau(end, 1:Dimension(2)) = Matrix_problem(end, 1:Dimension(2));
    T_Tableau(1:Dimension(1), end) = Matrix_problem(1:Dimension(1), end);
    var = get(handles.popupmenu_selectvar, 'string'); % se recupera la variable no b�sica seleccionada
    I = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2)); % se recupera variable seleccionada 
    J = str2double(var(get(handles.popupmenu_selectvar, 'value'), 3)); % se recupera el �ndice de variable seleccionada 
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
                Lsup = dim(1);
                change = 1;
            else
                Linf = dim(1);
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
    set(handles.popupmenu_selectvar, 'enable', 'off');
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
                    if (Node_NBV(1) == Node_current(1) || Node_NBV(2) == Node_current(2)) && Num_cassignation > 1
                        Num_assignation = Dimension(1)+Dimension(2)-3;
                    end
                    Empty_dimension = 0; 
                    break;
                else
                    if (Node_NBV(1) == Node_current(1) || Node_NBV(2) == Node_current(2)) && Num_cassignation > 1
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
                       if Num_assignation-1 == Num_cassignation %%%%
                            Num_assignation = Dimension(1)+Dimension(2)-3;%%%%
                       end%%%%
                       Empty_dimension = 0; 
                       break;%%%%
                    else
                        T_Tableau(Node_current(1), ind1(i)) = 0;              
                        if ~(i == Lsup)
                            Empty_dimension = 0;
                        end
                        %MODIFICADO(27/12/2016)
                        %Node_Cicle(Num_assignation, 1, 3) = 0;
                        %Num_assignation = Num_assignation - 1;                        
                        %MODIFICADO(27/12/2016)
                    end                     
                end
            end  
        else
            Empty_dimension = 0;
            T_Tableau(Node_current(1), Node_current(2)) = 0;
            Node_Cicle(Num_assignation, 1, 3) = 0;
            for l = 1:Dimension(1)+Dimension(2)-3;
                if Solution(l, 1, 1) == Node_current(1) && Solution(l, 1, 2) == Node_current(2)
                    Solution_change(l) = 0;
                    break;
                end
            end
            for k = (Dimension(1)+Dimension(2)-3):-1:1;                
                if Node_Cicle(k, 1, 3) == -1 || Node_Cicle(k, 1, 3) == 1
                    Num_cassignation = k;
                    %MODIFICADO(27/12/2016)
                    Num_assignation = k-1;
                    %MODIFICADO(27/12/2016)
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
            calc_nextciclevar(handles,Num_cassignation-1);
            return;
        end
    else
        ind1 = find(T_VarType_Aux(:, Node_current(2)));
        dim = size(ind1);
        if ~isempty(ind1)
            if Node_current(1) < ind1(1)
                Linf = 1;
                Lsup = dim(1);
                change = 1;
            else
                Linf = dim(1);
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
                    Empty_dimension = 1; 
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
                       if Num_assignation-1 == Num_cassignation
                            Num_assignation = Dimension(1)+Dimension(2)-3;
                       end
                       Empty_dimension = 1; 
                       break;
                    else
                        msgbox('El camino no es un ciclo. Se retornar� al punto anterior.', 'Camino en encrucijada.','modal');
                        T_Tableau(ind1(i), Node_current(2)) = 0;
                        if ~(i == Lsup)
                            Empty_dimension = 1;
                        end
                        %MODIFICADO(27/12/2016)
                        %Node_Cicle(Num_assignation, 1, 3) = 0;
                        %Num_assignation = Num_assignation -1;                        
                        %MODIFICADO(27/12/2016)
                    end
                end                                
            end  
        else
           Empty_dimension = 1;
           T_Tableau(Node_current(1), Node_current(2)) = 0;
           Node_Cicle(Num_assignation, 1, 3) = 0;
           for l = 1:Dimension(1)+Dimension(2)-3;
                if Solution(l, 1, 1) == Node_current(1) && Solution(l, 1, 2) == Node_current(2)
                    Solution_change(l) = 0;
                    break;
                end
           end
           for k = (Dimension(1)+Dimension(2)-3):-1:1;
                if Node_Cicle(k, 1, 3) == -1 || Node_Cicle(k, 1, 3) == 1
                    Num_cassignation = k;
                    %MODIFICADO(27/12/2016)
                    Num_assignation = k - 1;
                    %MODIFICADO(27/12/2016)
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
            calc_nextciclevar(handles,Num_cassignation-1);            
            return;
        end
    end
end

if Num_assignation == Dimension(1)+Dimension(2)-3
    msgbox('Se ha encontrado un ciclo de redistribuci�n.','C�lculo de ciclo de redistribuci�n.','modal');
    Solution_initial = 0;
    Node_current = [0, 0];
    set_environment('next_assign', handles);
    %set(handles.text_selectvar2, 'string', 'Seleccionar variable que sale');
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
    set(handles.Watch_varval, 'Enable', 'on');
    set(handles.popupmenu_selectvar3, 'string', char(Var));
    set(handles.popupmenu_selectvar3, 'value', 1);
    set(handles.popupmenu_selectvar3, 'enable', 'on');
    %set(handles.popupmenu_selectvar, 'enable', 'off');
end
All_display = T_Tableau;
set(handles.table_simplexdisplay, 'data', All_display);


% --- Executes on selection change in popupmenu_selectvar3.
function popupmenu_selectvar3_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to popupmenu_selectvar3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selectvar3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selectvar3
global Node_current Solution T_Tableau Matrix_problem Dimension Num_assignation ...
    Solution_initial NewSolution;

if handles.Method == 3
    Node_current = [0,0];
    Num_assignation = 0;
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
end


% --- Executes during object creation, after setting all properties.
function popupmenu_selectvar3_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to popupmenu_selectvar3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- se ejecuta al seleccionar la opci�n "Nuevo poblema/Simplex primal" del men� "Archivo"
function Simplex_primal_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_primal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

init_method(handles, 1);

% --- se ejecuta al seleccionar la opci�n "Nuevo poblema/Simplex dual" del men� "Archivo"
function Simplex_dual2_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_dual2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

init_method(handles, 2);

% --- se ejecuta al seleccionar la opci�n "Nuevo poblema/Simplex de transporte" del men� "Archivo"
function Simplex_transportation2_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_transportation2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

init_method(handles, 3);

function init_method(handles, method)
% ventana de dialogo en donde se solicita las dimensiones del problema

if method == 1 || method == 2
    prompt = {'N�mero de ecuaciones:','N�mero de variables:'};
    def = {'1','2'};    
elseif method == 3
    prompt = {'N�mero de origenes:','N�mero de destinos:'};
    def = {'2','2'}; 
end

dlg_title = 'Dimensiones del problema';
num_lines = 1;

answer = inputdlg(prompt,dlg_title,num_lines,def);

% se verifica que las dimensiones del problema sean consistentes (m < n)
if ~isempty(answer)    
    if method == 1
        handles.Method = 1;
        set(handles.panel_enhancement, 'title', 'Panel de control (Simplex primal)');
    elseif method == 2
        handles.Method = 2;
        set(handles.panel_enhancement, 'title', 'Panel de control (Simplex dual)');
    elseif method == 3
        handles.Method = 3;
    end
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
        % se abre la ventana para introducir la especificaci�n del problema
        handles.gui_Matrix_problem = zeros(user_entry1+1, user_entry2+1); % se comparte la especificaci�n del problema
        handles.First_Matrix_problem = zeros(user_entry1+1, user_entry2+1);
        handles.istwophases = 0;
        handles.whatphase = 1;
        handles.isinit_secondphase = 0;
        handles.Orig_Matrix_problem = zeros(user_entry1+1, user_entry2+1);
        handles.maxcanon_vector = 0;
        handles.setProblem = @setProblem; % se comparte el manejador de funci�n    
        handles.gui_Problem = Problem('LPApp', handles);
        guidata(handles.output, handles);
    end
    handles.latex = '';
end



% --- Executes on selection change in listbox_operations.
function listbox_operations_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to listbox_operations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_operations contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_operations


% --- Executes during object creation, after setting all properties.
function listbox_operations_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to listbox_operations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Generate_latex_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSL>
% hObject    handle to Generate_latex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

latexfile = '% OJO: El documento est� guardado en codificaci�n de caracteres ISO-8859-1';
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\documentclass[11pt]{article}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '%configuraci�n para idioma espa�ol'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage[latin1]{inputenc}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage{amsfonts}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage{amsmath}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage[spanish,activeacute,es-noindentfirst]{babel}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\title{Archivo generado por Asistente Simplex + 2016}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\author{Profesor M.Sc. Porfirio Armando Rodr�guez}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\date{\today}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\begin{document}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\maketitle'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile,handles.latex];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\end{document}'];

%latex = handles.latex;
handles.latexfile = ['LProblem', num2str(now),'.tex'];
dlmwrite(handles.latexfile, latexfile, '');
msgbox(['Se gener� el siguiente archivo LaTex:', handles.latexfile, '.'],'Archivo LaTex generado', 'modal');

function generar_latexsimplexnext()
global Tableau latex;

dim = size(Tableau);

latex = sprintf('%s \r La tabla que se obtiene es la siguiente: \r', latex);
latex = [latex, '\begin{equation}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\nonumber'];
latex = sprintf('%s \r', latex);
latex = [latex, '\left[ \begin{array}{']; 
for j = 1:dim(2)-1
    latex = [latex,'c'];
end
latex = [latex,'|c'];
latex = sprintf('%s} \r', latex);

for i = 1:dim(1)
    for j = 1:dim(2)
        if j == 1        
            latex = [latex, num2str(Tableau(i, j))]; %#ok<*AGROW>
        else
            latex = [latex,'& ', num2str(Tableau(i, j))];                
        end
    end
    latex = [latex, ' \\'];
    latex = sprintf('%s \r', latex);
end
latex = [latex, '\end{array} \right]']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex);

function generar_latexanalisis(title)
global Tableau latex;

latex = sprintf('%s \r', latex);
latex = [latex, title];
latex = sprintf('%s \r El valor �ptimo es: ', latex);
latex = [latex, '$\bf{', num2str(-Tableau(end, end)),'}$.'];
latex = sprintf('%s \r', latex);


function generar_latexspecification()
global Matrix_problem latex;

dim = size(Matrix_problem);

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, '\section{Segunda Fase del M�todo Simplex}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\subsection{Especificaci�n equivalente del problema original}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{equation}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\nonumber'];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{array}{']; 
for j = 1:dim(2)+1
    latex = [latex,'c']; %#ok<*AGROW>
end
latex = sprintf('%s} \r', latex);
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{minimizar \ '];
inicial = 0;
for j = 1:dim(2) - 1
    if Matrix_problem(dim(1), j) > 0
        if inicial == 1
            latex = [latex, ...
            '+', num2str(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
        else 
            latex = [latex, ...
            num2str(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
        end
        inicial = 1;
    elseif Matrix_problem(dim(1), j) < 0
        latex = [latex, ...
        num2str(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
        inicial = 1;
    end
    if j == dim(2) - 1
        latex = [latex, '} \\'];
    end 
end
latex = sprintf('%s \r', latex);
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{sujeto \ a} \\'];
latex = sprintf('%s \r', latex);
inicial = 0;
for i = 1:dim(1) - 1
    for j = 1:dim(2) - 1
        if Matrix_problem(i, j) > 0
            if inicial == 1
                latex = [latex, ...
                '+', num2str(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
            else
                latex = [latex, ...
                num2str(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
            end
            inicial = 1;
        elseif Matrix_problem(i, j) < 0
            latex = [latex, ...
            num2str(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
            inicial = 1;
        else
            latex = [latex, ' & '];
        end           
    end
    latex = [latex, ' = & ', num2str(Matrix_problem(i, dim(2))), '\\'];
    latex = sprintf('%s \r', latex);
end
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{x_j \ge 0 \ (j = 1, \dots,', num2str(dim(2)-1), ')}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\end{array}']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex);

function generar_latexsimplexbegin()
global Matrix_problem Order_initial latex;

Tableau = Matrix_problem;
dim = size(Tableau);

% se convierte la �ltima en t�rminos de Rj =cj - zj
[Y, I] = sort(Order_initial); %#ok<ASGLU>
Tableau_sorted = Matrix_problem(:, I);
c = Tableau_sorted(end,:);
z = Tableau_sorted(end, 1:(dim(1)-1))*Tableau_sorted(1:(dim(1)-1),:);

R = c - z;
Tableau(end, I) = R;

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, '\subsection{Desarrollo del m�todo simplex}'];
latex = sprintf('%s \r', latex);
latex = sprintf('%s \r La tabla inicial del m�todo simplex es la siguiente: \r', latex);
latex = [latex, '\begin{equation}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\nonumber'];
latex = sprintf('%s \r', latex);
latex = [latex, '\left[ \begin{array}{']; 
for j = 1:dim(2)-1
    latex = [latex,'c']; %#ok<AGROW>
end
latex = [latex,'|c'];
latex = sprintf('%s} \r', latex);

for i = 1:dim(1)
    for j = 1:dim(2)
        if j == 1        
            latex = [latex, num2str(Tableau(i, j))]; %#ok<AGROW>
        else
            latex = [latex,'& ', num2str(Tableau(i, j))];                 %#ok<AGROW>
        end
    end
    latex = [latex, ' \\']; %#ok<AGROW>
    latex = sprintf('%s \r', latex);
end
latex = [latex, '\end{array} \right]']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex);

function generar_latexnextsolution(title)
global Matrix_problem latex Solution;

dim = size(Matrix_problem);
ObjectiveValue = 0;
Matrix_solution = zeros(dim(1),dim(2));
for j=1:dim(1)+dim(2)-3        
    Matrix_solution(Solution(j, 1, 1), Solution(j, 1, 2)) = Solution(j, 1, 3);
    ObjectiveValue = ObjectiveValue + Solution(j, 1, 3)*Matrix_problem(Solution(j, 1, 1), Solution(j, 1, 2));
end
Matrix_solution(dim(1), dim(2)) = ObjectiveValue;

latex = sprintf('%s \r', latex);
latex = [latex, title];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{equation}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\nonumber'];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{tabular}{']; 
for j = 1:2*dim(2)
    latex = [latex,'c']; %#ok<AGROW>
end
latex = sprintf('%s} \r', latex);

for i = 1:dim(1)+1
    if i == 1
        for j = 1:dim(2)
            if j < dim(2)        
                latex = [latex, '&  \multicolumn{2}{c}{Destino ', num2str(j), '} ']; %#ok<AGROW>
            else
                latex = [latex, '&  \multicolumn{1}{c}{Oferta} \\']; %#ok<AGROW>                           
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    elseif i < dim(1)+1
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{Origen ', num2str(i-1),'} ']; %#ok<AGROW>
            elseif j < dim(2)+1
                latex = [latex, '&  \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{$x_{', num2str(i-1),',', num2str(j-1), '}$} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 3:2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
        latex = sprintf('%s \r', latex);
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{} ']; %#ok<AGROW>
            elseif j < dim(2)+1
                latex = [latex, '&  \multicolumn{1}{c}{',num2str(Matrix_solution(i-1, j-1)),'} & \multicolumn{1}{c|}{} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{',num2str(Matrix_problem(i-1, j-1)),'} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    else
        for j = 1:dim(2)
            if j > 1 
                if j < dim(2)
                    latex = [latex, '&  \multicolumn{1}{c}{', num2str(Matrix_problem(i-1, j-1)), '} &  \multicolumn{1}{c|}{} ']; %#ok<AGROW>
                else
                    latex = [latex, '&  \multicolumn{1}{c}{', num2str(Matrix_problem(i-1, j-1)), '} &  \multicolumn{1}{c|}{} '];
                end
            else
                latex = [latex, '\multicolumn{1}{c|}{Demanda} ']; %#ok<AGROW>                           
            end            
        end
        latex = [latex, '&  \multicolumn{1}{c|}{',num2str(Matrix_solution(i-1, j)),'} \\'];
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    end
    latex = sprintf('%s \r', latex);
end

latex = [latex, '\end{tabular}']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex);


function generar_latexnextrj()
global Matrix_problem T_Tableau latex;

dim = size(Matrix_problem);

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, '\subsection{Encontrando los costos reducidos asociados a la soluci�n actual}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{equation}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\nonumber'];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{tabular}{']; 
for j = 1:2*dim(2)
    latex = [latex,'c']; %#ok<AGROW>
end
latex = sprintf('%s} \r', latex);

for i = 1:dim(1)+1
    if i == 1
        for j = 1:dim(2)
            if j < dim(2)        
                latex = [latex, '&  \multicolumn{2}{c}{Destino ', num2str(j), '} ']; %#ok<AGROW>
            else
                latex = [latex, '&  \multicolumn{1}{c}{$u_i$} \\']; %#ok<AGROW>                           
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    elseif i < dim(1)+1
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{Origen ', num2str(i-1),'} ']; %#ok<AGROW>
            elseif j < dim(2)+1
                latex = [latex, '&  \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{$r_{', num2str(i-1),',', num2str(j-1), '}$} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 3:2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
        latex = sprintf('%s \r', latex);
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{} ']; %#ok<AGROW>
            elseif j < dim(2)+1
                latex = [latex, '&  \multicolumn{1}{c}{',num2str(T_Tableau(i-1, j-1)),'} & \multicolumn{1}{c|}{} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{',num2str(T_Tableau(i-1, j-1)),'} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    else
        for j = 1:dim(2)
            if j > 1        
                latex = [latex, '&  \multicolumn{1}{c}{', num2str(T_Tableau(i-1, j-1)), '} &  \multicolumn{1}{c|}{} ']; %#ok<AGROW>
            else
                latex = [latex, '\multicolumn{1}{c|}{$v_j$} ']; %#ok<AGROW>                           
            end            
        end
        latex = [latex, '&  \multicolumn{1}{c|}{} \\'];
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    end
    latex = sprintf('%s \r', latex);
end

latex = [latex, '\end{tabular}']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex);

function generar_latexnextcicle()
global Matrix_problem T_Tableau latex;

dim = size(Matrix_problem);

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, '\subsection{Encontrando un ciclo de redistribuci�n}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{equation}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\nonumber'];
latex = sprintf('%s \r', latex);
latex = [latex, '\begin{tabular}{']; 
for j = 1:2*dim(2)
    latex = [latex,'c']; %#ok<AGROW>
end
latex = sprintf('%s} \r', latex);

for i = 1:dim(1)+1
    if i == 1
        for j = 1:dim(2)
            if j < dim(2)        
                latex = [latex, '&  \multicolumn{2}{c}{Destino ', num2str(j), '} ']; %#ok<AGROW>
            else
                latex = [latex, '&  \multicolumn{1}{c}{Oferta} \\']; %#ok<AGROW>                           
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    elseif i < dim(1)+1
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{Origen ', num2str(i-1),'} ']; %#ok<AGROW>
            elseif j < dim(2)+1
                latex = [latex, '&  \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{$x_{', num2str(i-1),',', num2str(j-1), '}$} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 3:2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
        latex = sprintf('%s \r', latex);
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{} ']; %#ok<AGROW>
            elseif j < dim(2)+1
                latex = [latex, '&  \multicolumn{1}{c}{',num2str(T_Tableau(i-1, j-1)),'} & \multicolumn{1}{c|}{} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{',num2str(T_Tableau(i-1, j-1)),'} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    else
        for j = 1:dim(2)
            if j > 1        
                latex = [latex, '&  \multicolumn{1}{c}{', num2str(Matrix_problem(i-1, j-1)), '} &  \multicolumn{1}{c|}{} ']; %#ok<AGROW>
            else
                latex = [latex, '\multicolumn{1}{c|}{Demanda} ']; %#ok<AGROW>                           
            end            
        end
        latex = [latex, '&  \multicolumn{1}{c|}{} \\'];
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    end
    latex = sprintf('%s \r', latex);
end

latex = [latex, '\end{tabular}']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex);


% --------------------------------------------------------------------
function Modify_Callback(hObject, eventdata, handles) %#ok<DEFNU,,INUSL>
% hObject    handle to Modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem;

% abre ventana para modificar la especificaci�n del problema
handles.setProblem = @setProblem; % se comparte manejador de funci�n
%handles.gui_Matrix_problem = Matrix_problem; % se comparte la especificaci�n de problema actual

%AGREGADO(27/12/2016)
handles.latex = '';
if handles.istwophases == 1
    handles.gui_Matrix_problem = handles.Orig_Matrix_problem;
    handles.isinit_secondphase = 0;
else
    handles.gui_Matrix_problem = Matrix_problem; % se comparte la especificaci�n de problema actual
end
%AGREGADO(27/12/2016)

handles.gui_Problem = Problem('LPApp', handles);
guidata(handles.output, handles);


% --------------------------------------------------------------------
function Close_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD>
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(hObject);


% --------------------------------------------------------------------
function Multiple_solution_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to Multiple_solution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pushbutton_next_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Watch_varval_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Watch_varval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Solution Matrix_problem Dimension;

if handles.Method == 1
    colFormat=cell(1,Dimension(2)+1);    
else
    colFormat=cell(1,Dimension(2));
end
colFormat(1,:) = cellstr('char');
set(handles.table_simplexdisplay, 'columnformat', colFormat);
ObjectiveValue = 0;
All_display = cell(Dimension(1),Dimension(2));
for j=1:Dimension(1)+Dimension(2)-3
    All_display(Solution(j, 1, 1), Solution(j, 1, 2)) = cellstr(strcat('X', strcat(num2str(Solution(j, 1, 1)), num2str(Solution(j, 1, 2))),'|', num2str(Solution(j, 1, 3))));
    ObjectiveValue = ObjectiveValue + Solution(j, 1, 3)*Matrix_problem(Solution(j, 1, 1), Solution(j, 1, 2));
end
All_display(Dimension(1), Dimension(2)) = cellstr(num2str(ObjectiveValue));
    
set(handles.table_simplexdisplay, 'data', All_display);




% --------------------------------------------------------------------
function Next_value_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Next_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Num_assignation Dimension latex;

if strcmp(get(hObject, 'Label'), 'Siguiente multiplicador')  
    paso = 2;
    calc_nextmultiplier(handles);
elseif strcmp(get(hObject, 'Label'), 'Siguiente valor')    
    calc_nextassignment(handles);
    paso = 1;
elseif strcmp(get(hObject, 'Label'), 'Siguiente nodo')    
    calc_nextciclevar(handles, Num_assignation);
    paso = 3;
end

if (Num_assignation == Dimension(1)+Dimension(2)-3)
    Num_assignation = 0;
    latex = '';
    if paso == 1
        handles.latex = [handles.latex, latex];
    elseif paso == 2
        handles.latex = [handles.latex, latex];
    else
        generar_latexnextcicle();
        handles.latex = [handles.latex, latex];
    end
    guidata(handles.output, handles);
end



% --- Executes on selection change in popupmenu_selectvar2.
function popupmenu_selectvar2_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to popupmenu_selectvar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selectvar2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selectvar2
global Node_current Solution T_Tableau Matrix_problem Dimension Num_assignation T_VarType_Aux ...
    Solution_initial T_VarType NewSolution;


% se calculan las nuevas variables no b�sicas si las hay
if strcmp(get(handles.Next_value, 'Label'), 'Siguiente multiplicador')
    Node_current = [0, 0];
    Num_assignation = 0;
    set_environment('next_calc', handles);
    T_Tableau = zeros(Dimension(1), Dimension(2));
    T_Tableau(1:Dimension(1)-1, 1:Dimension(2)-1) = Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1);
    T_VarType_Aux = T_VarType;
    calc_nextmultiplier(handles);    
elseif strcmp(get(handles.Next_value, 'Label'), 'Siguiente valor')
    Node_current = [0,0];
    Num_assignation = 0;
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
elseif strcmp(get(handles.Next_value, 'Label'), 'Siguiente nodo')
    Num_assignation = 0;
    T_VarType_Aux = T_VarType;
    %calc_nextciclevar(handles, [0, 0]);
    calc_nextciclevar(handles, 0);
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
