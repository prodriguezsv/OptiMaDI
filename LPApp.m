% Copyright 2011 2012 2016 2017 M.Sc. Porfirio Armando Rodríguez

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

format long;
handles.Method = 0;
handles.Newmethod = 0;
handles.Submethod = 0;
handles.isonlylecture = 0;
handles.latex = '';
handles.latexproblem = '';
handles.latexIIphasesproblem = '';
handles.latexfile = '';
handles.gui_Matrix_problem = []; % se comparte la especificación del problema
handles.First_Matrix_problem = [];
handles.istwophases = 0;
handles.whatphase = 1;
handles.isinit_secondphase = 0;
handles.Orig_Matrix_problem = [];
handles.maxcanon_vector = 0;
handles.nodelist = cell(1);
handles.nodelist{1} = [];
handles.ispenalties = 0;
handles.ZOptimo = Inf;
handles.XOptimo = [];
handles.Level = 0;
handles.Travel = 1;
handles.Interactive = 1;

% Update handles structure
guidata(hObject, handles);

movegui('center');
% UIWAIT makes LPApp wait for user response (see UIRESUME)
% uiwait(handles.LPApp);

set(handles.table_simplexdisplay, 'data', cell(100,100));

% --- Outputs from this function are returned to the command line.
function varargout = LPApp_OutputFcn(hObject, eventdata, handles) %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

format short;

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Función utilitaria: configura los distintos ambientes de ejecución
function set_environment(environmentname, handles)
% environmentname especifica el tipo de ambiente por configurar
% handles         estructura con manejadores y datos de usuario
global Variables_ind;

if strcmp(environmentname, 'next')
    if handles.Method ~= 3        
        set(handles.Simplex_multiplier, 'Enable', 'on');
    else
        set(handles.Watch_varval, 'Enable', 'off');
    end
    Variables_ind = 0;
    set(handles.Sensibility, 'enable', 'off');
    set(handles.Postoptimality, 'enable', 'off');
    set(handles.slider_increment, 'value', 0);
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.popupmenu_selectvar3, 'Enable', 'on');
    %set(handles.Restriction_nonnegativity, 'Checked', 'off');
    set(handles.Uncut_planes, 'Enable', 'off');   
    set(handles.pushbutton_asignall, 'Enable', 'off');
    set(handles.uipushtool1, 'Enable', 'off');
    set(handles.uitoggletool6, 'enable', 'off');
    val_switches = ['on ';'on ';'on ';'on ';'on ';'off'];   
elseif strcmp(environmentname, 'sol_multiples')
    set(handles.Sensibility, 'enable', 'on');
    set(handles.Postoptimality, 'enable', 'on');
    val_switches = ['on ';'off';'on ';'off';'on ';'on '];
    set(handles.popupmenu_selectvar3, 'Enable', 'on');
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.Simplex_multiplier, 'Enable', 'on');
elseif strcmp(environmentname, 'end')
    if handles.Method ~= 3
        set(handles.Sensibility, 'enable', 'on');
        set(handles.Postoptimality, 'enable', 'on');
        set(handles.Watch_varval, 'Enable', 'off');
        set(handles.Simplex_multiplier, 'Enable', 'on');
    else
        set(handles.Watch_varval, 'Enable', 'on');
        set(handles.Simplex_multiplier, 'Enable', 'off');
    end
    set(handles.slider_increment, 'value', 0);    
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.Next_value, 'Enable', 'off');
    set(handles.pushbutton_asignall, 'Enable', 'off');     
    set(handles.popupmenu_selectvar, 'Enable', 'off');
    set(handles.popupmenu_selectvar, 'string', char(' ', ' '));
    set(handles.popupmenu_selectvar, 'value', 1);
    set(handles.popupmenu_selectvar2, 'Enable', 'off');
    set(handles.popupmenu_selectvar2, 'string', char(' ', ' '));
    set(handles.popupmenu_selectvar2, 'value', 1);
    set(handles.popupmenu_selectvar3, 'Enable', 'off');
    set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
    set(handles.popupmenu_selectvar3, 'value', 1);
    val_switches = ['off';'off';'on ';'off';'on '; 'off'];
elseif strcmp(environmentname, 'next_assign')
    set(handles.Restriction_nonnegativity, 'Enable', 'off');
    set(handles.Uncut_planes, 'Enable', 'off');
    set(handles.Mode_geo3D, 'Enable', 'off');
    set(handles.Sensibility, 'enable', 'off');
    set(handles.Postoptimality, 'enable', 'off');
    setTableauTags(handles, char('O', 'D', 'Rj', 'Ci'));
    %set(handles.popupmenu_selectvar2, 'enable', 'on');
    %set(handles.popupmenu_selectvar3, 'enable', 'off');
    set(handles.Next_value, 'Label', 'Siguiente valor');
    set(handles.pushbutton_asignall, 'string', 'Asignar todo');
    set(handles.Next_value, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on');
    set(handles.Watch_varval, 'Enable', 'off');  
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.uipushtool1, 'Enable', 'off');
    set(handles.uitoggletool6, 'enable', 'off');
    val_switches = ['off';'off';'on ';'off';'on ';'off'];
elseif strcmp(environmentname, 'next_newassign')
    set(handles.Restriction_nonnegativity, 'Enable', 'off');
    set(handles.Uncut_planes, 'Enable', 'off');
    set(handles.Mode_geo3D, 'Enable', 'off');
    set(handles.Sensibility, 'enable', 'off');
    set(handles.Postoptimality, 'enable', 'off');
    setTableauTags(handles, char('O', 'D', 'Rj', 'Ci'));
    set(handles.Next_value, 'Label', 'Siguiente valor');
    set(handles.pushbutton_asignall, 'string', 'Asignar todo');
    set(handles.Next_value, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'off');
    set(handles.Watch_varval, 'Enable', 'off');  
    set(handles.Saveproblem, 'enable', 'on');
    set(handles.uipushtool1, 'Enable', 'off');
    set(handles.uitoggletool6, 'enable', 'off');
    val_switches = ['off';'off';'on ';'on ';'on ';'off'];    
elseif strcmp(environmentname, 'next_calc')    
    setTableauTags(handles, char('O', 'D', 'Vj', 'Ui'));
    set(handles.Next_value, 'Label', 'Siguiente multiplicador');
    set(handles.pushbutton_asignall, 'string', 'Calcular todo');   
    set(handles.popupmenu_selectvar2, 'Enable', 'on');       
    set(handles.popupmenu_selectvar, 'Enable', 'off');
    set(handles.popupmenu_selectvar, 'string', char(' ', ' ')); 
    set(handles.popupmenu_selectvar, 'value', 1);
    set(handles.popupmenu_selectvar3, 'Enable', 'off');
    set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
    set(handles.popupmenu_selectvar3, 'value', 1);
    set(handles.Next_value, 'Enable', 'on');
    set(handles.pushbutton_asignall, 'Enable', 'on');
    set(handles.Watch_varval, 'Enable', 'on');
    set(handles.Saveproblem, 'enable', 'on');
    val_switches = ['off';'off';'on ';'off';'on ';'off']; 
elseif strcmp(environmentname, 'next_cicle')
    setTableauTags(handles, char('O', 'D', 'Rj', 'Ci'));
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

% --- Función utilitaria: configura el ambiente en cada paso del Método
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
% Matrix_problem    especificación del problema
% Order_initial     orden inicial de vectores columna de la tabla Simplex, el orden de la
% matrix identidad primero y luego los demás vectores

%AGREGADO(27/12/2016)

if handles.istwophases == 1 && handles.isinit_secondphase ~= 1
    handles.latex = handles.latexIIphasesproblem;
    handles.whatphase = 1;
    Matrix_problem = handles.First_Matrix_problem;
    handles.Order_initial = handles.First_Order_initial;  
    
    setProblem(Matrix_problem, handles);
elseif handles.istwophases == 1 && handles.isinit_secondphase == 1    
    handles.isinit_secondphase = 0;
    Matrix_problem = handles.gui_Matrix_problem;
    generar_latexspecification();
    handles.latex = [handles.latex, latex];
    Order_initial = handles.Order_initial;
    generar_latexsimplexbegin()
    handles.latex = [handles.latex, latex];
    guidata(handles.output, handles);
    
    handles.isonlylecture = 1; latex = '';
    handles.setProblem = @setProblem; % se comparte el manejador de función    
    handles.gui_Problem = Problem('LPApp', handles); 
    %handles.isinit_secondphase = 0;
    %guidata(handles.output, handles);
elseif handles.Method == 4
    Matrix_problem = handles.First_Matrix_problem;
    handles.Order_initial = handles.First_Order_initial;
    handles.nodelist = cell(1);
    handles.nodelist{1} = [];
    handles.ispenalties = 0;
    handles.ZOptimo = Inf;
    handles.XOptimo = [];
    handles.Level = 0;
    handles.latex = handles.latexproblem;
    set(handles.pushbutton_next, 'String', 'Mejorar solucion');
    set(handles.popupmenu_selectvar2, 'Enable', 'on');
    Var = char('LIFO-a', 'LIFO-d', 'FIFO-a', 'FIFO-d');
    set(handles.popupmenu_selectvar2, 'string', Var);
    set(handles.popupmenu_selectvar2, 'value', 1);
    setProblem(Matrix_problem, handles);
else
    %set(handles.Next_value, 'Enable', 'on');
    handles.latex = handles.latexproblem;
    handles.Order_initial = Order_initial;
    
    setProblem(Matrix_problem, handles);
end

guidata(handles.output, handles);
%AGREGADO(27/12/2016)

% Configura la tabla inicial del Método Simplex
% handles.Order_initial = Order_initial;
% se abre la ventana para introducir la especificación del problema

if strcmp(get(handles.Mode_geo3D, 'Checked'), 'on')
    %set(handles.Restriction_nonnegativity, 'Checked', 'off')
    %set(handles.Restriction_nonnegativity, 'Enable', 'on')
    trace3D(handles);
end


% --- Función utilitaria: obtiene las variables no básicas candidatas a
% ser seleccionadas para convertirse en básicas según el criterio de
% mejoramiento de la solución (coeficientes de costo relativo negativos)
function h = calc_variables(handles)
% handles       estructura con manejadores y datos de usuario
global Dimension Tableau Order_current Solution_initial Matrix_problem latex output_cod output_sol output_val ...
    VarSelect VarSelect2 VarSelect3 Jp Jq; %#ok<NUSED>
% Dimension     dimensiones de la especificación del problema en términos
% del número de ecuaciones y variables
% Tableau       Tabla del Simplex actual
% Order_current orden actual de vectores columna de la tabla Simplex

n = Dimension(2)-1;
m = Dimension(1)-1;

if handles.Method == 1 || handles.Submethod == 1
    %Rj = Tableau(end,1:(Dimension(2)-1)); % Recupera los coeficientes de costo reducido (Rj) de la tabla Simplex
    Rj = Tableau(m+1,1:n); % Recupera los coeficientes de costo reducido (Rj) de la tabla Simplex
    [Y, I] = sort(Order_current); %#ok<ASGLU> % obtiene los indices de variables ordenados de menor a mayor
    IRjneg = I(Order_current(Rj < 0));  % obtiene los indices de variables con Rj negativo
    IVar = IRjneg;
elseif handles.Method == 2 || handles.Submethod == 2
    %X0 = Tableau(1:(Dimension(1)-1), end); % Recupera los valores de las variables básicas de la tabla Simplex
    X0 = Tableau(1:m, end); % Recupera los valores de las variables básicas de la tabla Simplex
    [Y, I] = sort(Order_current); %#ok<ASGLU> % obtiene los indices de variables ordenados de menor a mayor
    IX0neg = I(X0 < 0);  % obtiene los indices de variables con valores negativos
    IVar =  IX0neg;
end

if handles.Method ~= 3
    % se construye el arreglo de cadenas de las variables
    dim2 = size(IVar);
    Var = cell(dim2(2), 1);
    j = 0;
    if handles.Method == 4    
        j = handles.Level;
    end
    for i = 1:dim2(2)
        Var(i) = cellstr(strcat('X',num2str(IVar(i)+j)));
    end
    
    VarSelect = Var;

    if ~isempty(Var) % en el caso que haya Rj negativos
        if handles.Method == 1 || handles.Submethod == 1
            [Y, p] = min(Rj(IRjneg)); %#ok<ASGLU> % obtiene el índice de variable del Rj más negativo
        elseif handles.Method == 2 || handles.Submethod == 2
            [Y, p] = min(IX0neg); %#ok<ASGLU> % obtiene el índice de variable del X0i más negativo
        end
        if handles.Interactive == 1
        set(handles.popupmenu_selectvar, 'string', char(Var)); % despliega las variables no básicas en la interfaz
        set_environment('next', handles); % ajusta el ambiente
        set(handles.popupmenu_selectvar, 'value', p); % se selecciona, por omisión, la variable más negativa
        end
        Jp = p; calc_ratios(handles); % se calculan las razones para aplicar el criterio de factibilidad
        [Z, I] = sort(Order_current); %#ok<ASGLU>   
        bacic_sol = zeros(Dimension(2)-1,1);
        bacic_sol(I(1:Dimension(1)-1)) = Tableau(1:Dimension(1)-1,Dimension(2));
        output_cod = 0; output_sol = bacic_sol; output_val = -Tableau(end,end);
    else % verificar si hay Rj iguales a cero
        j = 1;
        %Rj = Tableau(end,1:(Dimension(2)-1)); % Recupera los coeficientes de costo reducido (Rj) de la tabla Simplex
        Rj = Tableau(end,1:n); % Recupera los coeficientes de costo reducido (Rj) de la tabla Simplex
        %IRjeqaux = I(Rj(I(Dimension(1):(Dimension(2)-1)))==0);
        IRjeqaux = I(Rj(I(m+1:n))==0);
        dim3 = size(IRjeqaux);
        IRjeq = zeros(1, dim3(2));
        for i = m+1:n
            if Rj(I(i)) == 0
                IRjeq(j) = I(i);
                j = j+1;
            end
        end

        dim2 = size(IRjeq);
        Var = cell(dim2(2), 1);
        j = 0;
        if handles.Method == 4    
            j = handles.Level;
        end
        for i = 1:dim2(2)
            Var(i) = cellstr(strcat('X',num2str(IRjeq(i)+j)));
        end
        VarSelect = Var;
        
        %MODIFICADO(27/12/2016)            
        %maximo = max(find(Order_current < Dimension(1))); %#ok<MXFND>
        maximo = max(find(Order_current < m+1)); %#ok<MXFND>
        if ~isempty(Var) && ((handles.Method == 1 || handles.Method == 2)  || ...
                    ((handles.istwophases == 1 || handles.ispenalties == 1) && maximo > n -(m-handles.maxcanon_vector))) % en el caso que haya Rj = 0        
            if handles.Interactive == 1
            if handles.Method == 1 || handles.Submethod == 1
                set(handles.popupmenu_selectvar, 'string', char(Var));
            elseif handles.Method == 2 || handles.Submethod == 2
                set(handles.popupmenu_selectvar, 'string', ['X', num2str(maximo)]);
            end
            set_environment('sol_multiples', handles);
            set(handles.popupmenu_selectvar, 'value', 1);
            end
            Jp = 1; calc_ratios(handles);
            l = 0;
            if (handles.ispenalties == 1)
                l = Dimension(1)-1 - handles.maxcanon_vector;
            end
            [Z, I] = sort(Order_current); %#ok<ASGLU>  
            bacic_sol = zeros(Dimension(2)- 1 - l,1);
            bacic_sol(I(1:Dimension(1)-1)) = Tableau(1:Dimension(1)-1,Dimension(2));
            output_cod = 4; output_sol = bacic_sol; output_val = -Tableau(end,end);
            if handles.Interactive == 1
            if strcmp(get(handles.Mode_geo3D, 'Checked'), 'off')
                msgbox(char('El proceso ha encontrado una solución óptima. Haga clic en "Solucion múltiple" para encontrar otra.'),char('Método simplex'),'modal');
            end
            %latex = '';
            if handles.Method == 4
                generar_latexanalisis('\subsection{Análisis de resultados}');
            else
                generar_latexanalisis('\section{Análisis de resultados}');
            end            
            %handles.latex = [handles.latex, latex];
            %guidata(handles.output, handles);
            end
        else % No hay Soluciones múltiples, no hay más opciones
            if (handles.istwophases == 1 && handles.whatphase == 1) 
                if Tableau(end, end) == 0 && maximo <= n-(m-handles.maxcanon_vector) %%OJO
                    if handles.Interactive == 1
                    if strcmp(get(handles.Mode_geo3D, 'Checked'), 'off')
                        msgbox('La primera fase ha terminado. Haga clic en "Iniciar" para la segunda fase.',char('Método simplex de dos fases'),'modal');
                    end
                    end
                    [Z, I] = sort(Order_current); %#ok<ASGLU> 
                    bacic_sol = zeros(Dimension(2)-1,1);
                    bacic_sol(I(1:Dimension(1)-1)) = Tableau(1:Dimension(1)-1,Dimension(2));
                    output_cod = 0; output_sol = bacic_sol; output_val = -Tableau(end,end); 
                    if handles.Interactive == 1
                    %latex = '';
                    
                    generar_latexanalisis('\subsection{Análisis de resultados}');
                    %handles.latex = [handles.latex, latex];
                    %guidata(handles.output, handles);
                    end
                    %Matrix_problem = zeros(size(handles.Orig_Matrix_problem)); 
                    %mcv = handles.maxcanon_vector;
                    Matrix_problem = Tableau(:, 1:(n+1-(m-handles.maxcanon_vector)));
                    Matrix_problem(:, n+1-(m-handles.maxcanon_vector)) = Tableau(:,n+1);
                    Matrix_problem(end, :) = handles.Orig_Matrix_problem(end, :);
                    Order_current = Order_current(1:(n+1-(m-handles.maxcanon_vector)));
                    Order_current(n+1-(m-handles.maxcanon_vector)) = n+1-(m-handles.maxcanon_vector);
                    Order_current(Order_current > Dimension(1)-1) = m+1:(n+1-(m-handles.maxcanon_vector));
                    handles.gui_Matrix_problem = Matrix_problem;
                    Dimension = size(Matrix_problem);
                    handles.Order_initial = Order_current;
                    handles.whatphase = 2;
                    handles.isinit_secondphase = 1;
                    %guidata(handles.output, handles);
                    %setProblem(Matrix_problem, handles)
                    %return;
                else
                    output_cod = 3; output_sol = []; output_val = [];
                    if handles.Interactive == 1
                    msgbox(char('La primera fase ha terminado. El problema original no tiene solución.'),char('Método simplex de dos fases'),'modal');
                    end
                    generar_latexanalisis('\subsection{Análisis de resultados}');
                end
            elseif (handles.istwophases == 1 && handles.whatphase == 2)
                [Z, I] = sort(Order_current); %#ok<ASGLU> 
                bacic_sol = zeros(Dimension(2)-1,1);
                bacic_sol(I(1:Dimension(1)-1)) = Tableau(1:Dimension(1)-1,Dimension(2));
                output_cod = 1; output_sol = bacic_sol; output_val = -Tableau(end,end);    
                if handles.Interactive == 1
                if strcmp(get(handles.Mode_geo3D, 'Checked'), 'off')
                    msgbox(char('La segunda fase ha terminado. Se ha encontrado una solución óptima'),char('Método simplex de dos fases'),'modal');
                end         
                %latex = '';
                generar_latexanalisis('\subsection{Análisis de resultados}'); 
                %handles.latex = [handles.latex, latex];
                %guidata(handles.output, handles);
                end
            else
                if (handles.ispenalties == 1) && maximo > n-(m-handles.maxcanon_vector) %%OJO                        
                    output_cod = 3; output_sol = []; output_val = [];
                    if handles.Interactive == 1
                    msgbox(char('El problema original no tiene solución.'),char('Método simplex de penalizaciones'),'modal');
                    end
                    generar_latexanalisis('\subsection{Análisis de resultados}');
                else
                    l = 0;
                    if (handles.ispenalties == 1)
                        l = Dimension(1)-1 - handles.maxcanon_vector;
                    end
                    [Z, I] = sort(Order_current); %#ok<ASGLU>  
                    bacic_sol = zeros(Dimension(2) - 1 - l,1);
                    bacic_sol(I(1:Dimension(1)-1)) = Tableau(1:Dimension(1)-1,Dimension(2));
                    output_cod = 1; output_sol = bacic_sol; output_val = -Tableau(end,end);
                    if handles.Interactive == 1
                    %latex = '';
                    if handles.Method == 4
                        generar_latexanalisis('\subsection{Análisis de resultados}');
                    else
                        generar_latexanalisis('\section{Análisis de resultados}');
                    end
                    %handles.latex = [handles.latex, latex];
                    %guidata(handles.output, handles);
                    if strcmp(get(handles.Mode_geo3D, 'Checked'), 'off')
                        msgbox(char('Se ha encontrado una solución óptima'),char('Método simplex'),'modal'); 
                    end
                    end
                end
            end
            if handles.Interactive == 1
            set(handles.popupmenu_selectvar, 'string', char(' ',' '));
            set(handles.popupmenu_selectvar3, 'string', char(' ',' '));
            set_environment('end', handles);
            end
        end
        %MODIFICADO(27/12/2016)
    end
else
    % se construye el arreglo de cadenas de las variables    
    Var = cell(n, 1);
    for i = 1:n
        Var(i) = cellstr(strcat('X1',num2str(i)));
    end
    VarSelect2= Var;
    Solution_initial = 1;
    set_environment('next_assign', handles); % ajusta el ambiente
    set(handles.popupmenu_selectvar2, 'string', char(Var));
    set(handles.popupmenu_selectvar2, 'value', 1);
    calc_nextassignment(handles);
end

h = handles;

% --- Función utilitaria: calcula las razones del criterio de factibilidad
function calc_ratios(handles)
% handles   estructura con manejadores y datos de usuario
global Tableau Dimension Order_current VarSelect VarSelect3 Jp Jq;
% All_display condensa los datos que se desplegarán en la interfaz


n = Dimension(2)-1;
m = Dimension(1)-1;

j = 0;
if handles.Method == 4    
    j = handles.Level;
end
if handles.Interactive == 1
var = get(handles.popupmenu_selectvar, 'string'); % se recupera la variable no básica seleccionada
dim = size(var(get(handles.popupmenu_selectvar, 'value'), :));
J = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2:dim(2))); % se recupera el índice de variable seleccionada
else
    var = VarSelect; % se recupera la variable no básica seleccionada
    dim = size(var(Jp, :));
    J = str2double(var(Jp, 2:dim(2))); % se recupera el índice de variable seleccionada
end
if handles.Method == 1 || handles.Submethod == 1
    % Se calculan las razones para la variable no básica seleccionada
    Y0 = Tableau(1:m, end); % se recupera el vector de términos coeficientes
    Yj = Tableau(1:m, J-j); % se recupera el vector columna asociada a la variable no básica
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
            Var(i) = cellstr(strcat('X',num2str(I(ind(i))+j)));
        end
        VarSelect3 = Var; Jq = 1;
        if handles.Interactive == 1
        set(handles.popupmenu_selectvar3, 'string', char(Var));
        set(handles.popupmenu_selectvar3, 'value', 1);
        end
    else
        if handles.Interactive == 1
        set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
        set(handles.popupmenu_selectvar3, 'value', 1);
        end
    end
    
    if handles.Interactive == 1
    All_display = zeros(m+1, n+1+1);
    % se actualiza la tabla de la interfaz
    All_display(:, 1:n+1) = Tableau;
    All_display(1:m, end) = ratios;
    
    spreadsheet = cell(100,100);
    spreadsheet(1:m+1, 1:n+1+1) = num2cell(All_display);
    end
elseif handles.Method == 2 || handles.Submethod == 2
    % Se calculan las razones para la variable básica seleccionada
    Fi = Tableau(Order_current(J-j), 1:n); % se recupera el vector fila asociada a la variable básica
    Rj = Tableau(end, 1:n); % se recupera el vector de coeficientes reducidos
    ratios = -Rj./Fi;
    ratios_aux = ratios;
    Basic_Var = Order_current < m+1;
    ratios_aux(Basic_Var) = Inf;
    %ratios_aux(ratios < 0 | isnan(ratios) | Fi > 0) = Inf; 
    ratios_aux(ratios < 0 | isnan(ratios)) = Inf; 
    [C, p] = min(ratios_aux); %#ok<NASGU>
    if C ~= Inf
        ind = find(ratios_aux == C);
        dim = size(ind);
        if handles.Interactive == 1
        % se construye el arreglo de cadenas de las variables    
        Var = cell(dim(2), 1);
        for i = 1:dim(2)
            Var(i) = cellstr(strcat('X',num2str(ind(i)+j)));
        end
        VarSelect3 = Var; Jq = 1;
        set(handles.popupmenu_selectvar3, 'string', char(Var));
        set(handles.popupmenu_selectvar3, 'value', 1);
        end
    else
        if handles.Interactive == 1
        set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
        set(handles.popupmenu_selectvar3, 'value', 1);
        end
    end
    if handles.Interactive == 1
    All_display = zeros(m+1+1, n+1);
    % se actualiza la tabla de la interfaz
    All_display(1:m+1, :) = Tableau;
    All_display(end, 1:n) = ratios;
    
    spreadsheet = cell(100,100);
    spreadsheet(1:m+1+1, 1:n+1) = num2cell(All_display);
    end
end


%MODIFICADO 1/01/2017
set(handles.table_simplexdisplay, 'data', spreadsheet);
%set(handles.table_simplexdisplay, 'data', All_display); 
%MODIFICADO 1/01/2017

% ----------------
function calc_nextassignment(handles)
global Solution Node_current Num_assignation T_Tableau Dimension Matrix_problem T_VarType ...
    Solution_change Minimo Solution_initial Node_NBV ;

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
    ObjectiveValue = 0;
    for j=1:Num_assignation        
        All_display(Solution(j, 1, 1), Solution(j, 1, 2)) = Solution(j, 1, 3);
        ObjectiveValue = ObjectiveValue + Solution(j, 1, 3)*Matrix_problem(Solution(j, 1, 1), Solution(j, 1, 2));
    end  
    All_display(end, :) = T_Tableau(end, :);
    All_display(:, end) = T_Tableau(:, end);   
    All_display(Dimension(1), Dimension(2)) = ObjectiveValue;
    spreadsheet = cell(100,100);
    spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
    set(handles.table_simplexdisplay, 'data', spreadsheet);
    %set(handles.table_simplexdisplay, 'data', All_display);
    
    if Num_assignation == Dimension(1)+Dimension(2)-3
        msgbox(char('Se ha encontrado una solución inicial.'),char('Cálculo de nueva solución.'),'modal');
        generar_latexnextsolution('\subsection{Encontrando una solución inicial por sustitución hacia atrás}');
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
    var = get(handles.popupmenu_selectvar3, 'string'); % se recupera la variable no básica seleccionada
    I = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 2)); % se recupera variable seleccionada 
    J = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 3)); % se recupera el índice de variable seleccionada
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
    
    spreadsheet = cell(100,100);
    spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
    set(handles.table_simplexdisplay, 'data', spreadsheet);
    %set(handles.table_simplexdisplay, 'data', All_display);
    if Num_assignation == Dimension(1)+Dimension(2)-3
        msgbox(char('Se ha encontrado una nueva solución.'),char('Cálculo de nueva solución.'),'modal');
        generar_latexnextsolution('\subsection{Encontrando una nueva solución por redistribución}');
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
global latex Num_assignation Dimension;

if handles.Method == 1
    if handles.Interactive == 1
    latex = '';
    end
    % se ejecuta el método simplex
    simplex_primal(handles);
elseif handles.Method == 2
    if handles.Interactive == 1
    latex = '';
    end
    simplex_dual(handles);
elseif handles.Method == 3
    if strcmp(get(handles.Next_value, 'Label'), 'Siguiente valor')
        while (Num_assignation ~= Dimension(1)+Dimension(2)-3)
            calc_nextassignment(handles);
        end

        Num_assignation = 0;
    end
    handles.latex = [handles.latex, latex];
elseif handles.Method == 4
    if handles.Interactive == 1
    latex = '';
    end
    branch_and_bound(handles);
end

if handles.Interactive == 1
if handles.Method ~= 3 && handles.Method ~= 4
    if strcmp(get(handles.Mode_geo3D, 'Checked'), 'on')
        trace3D(handles);
        set(handles.popupmenu_selectvar2, 'Enable', 'on');
        set(handles.pushbutton_asignall, 'Enable', 'on');
    end
end
end
    
% --- Se ejecuta el método simplex primal
function [output_cod, output_sol, output_val, h] = simplex_primal(handles)
global Tableau Order_current operations Dimension latex output_cod output_sol output_val ...
    VarSelect VarSelect3 Jp Jq; %#ok<REDEF>

%output_cod = 0; output_sol = [];  output_val = [];
% recupera el índice de la variable no básica seleccionada asociada a la
% columna que ingresará a la base
j = 0;
if handles.Method == 4
    j = handles.Level;
end
if handles.Interactive == 1
var = get(handles.popupmenu_selectvar, 'string');
dim = size(var(get(handles.popupmenu_selectvar, 'value'), :));
q = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2:dim(2)))-j;
% recupera el índice de la fila (o componente) del vector columna que
% dejará la base según la razón mínima tal que sea positiva
var = get(handles.popupmenu_selectvar3, 'string');
dim = size(var(get(handles.popupmenu_selectvar3, 'value'), :));
p_aux = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 2:dim(2)));
else
var = VarSelect;
dim = size(var(Jp, :));
q = str2double(var(Jp, 2:dim(2)))-j;
% recupera el índice de la fila (o componente) del vector columna que
% dejará la base según la razón mínima tal que sea positiva
var = VarSelect3;
dim = size(var(Jq, :));
p_aux = str2double(var(Jq, 2:dim(2)));    
end

if ~isnan(p_aux) % si existe algún yij que cumple el criterio de factibilidad 
    p = Order_current(p_aux-j);
    pivote_process(p, q, handles);
    if handles.Interactive == 1
    %AGREGADO(27/12/2016)
    set(handles.listbox_operations, 'string', operations);
    generar_latexsimplexnext()
    %AGREGADO(27/12/201
    end
else
    output_cod = 2; output_sol = [];  output_val = []; %#ok<NASGU>
    if handles.Interactive == 1
    operations = char('El conjunto representado por las restricciones no es acotado:', 'El valor objetivo decrece sin límite.');
    if strcmp(get(handles.Mode_geo3D, 'Checked'), 'off')
        msgbox('El conjunto representado por las restricciones no es acotado.', char('El valor objetivo decrece sin límite'),'modal');
    end
    latex = '';
    if handles.Method == 4
        generar_latexanalisis('\subsection{Análisis de resultados}');
    else
        generar_latexanalisis('\section{Análisis de resultados}');
    end
    handles.latex = [handles.latex, latex];
    guidata(handles.output, handles);
    set(handles.listbox_operations, 'string', operations);
    end
    h = handles; return;
end
if handles.Interactive == 1
% se actualiza la tabla de la interfaz
All_display = Tableau;
spreadsheet = cell(100,100);
spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
set(handles.table_simplexdisplay, 'data', spreadsheet);
%set(handles.table_simplexdisplay, 'data', All_display);
setrowheaders(handles, char('X', 'X', 'Rj', 'Yi0'));
end
handles = calc_variables(handles);
if handles.Interactive == 1
handles.latex = [handles.latex, latex];
guidata(handles.output, handles);
end
h = handles;

% --- Se ejecuta el método simplex dual
function [output_cod, output_sol, output_val, h] = simplex_dual(handles)
global Tableau Order_current operations Dimension latex output_cod output_sol output_val ...
    VarSelect VarSelect3 Jp Jq; %#ok<REDEF>

%output_cod = 0; output_sol = [];  output_val = [];
j = 0;
if handles.Method == 4    
    j = handles.Level;
end
if handles.Interactive == 1
% recupera el índice de la variable no básica seleccionada asociada a la
% columna que ingresará a la base
var = get(handles.popupmenu_selectvar, 'string');
dim = size(var(get(handles.popupmenu_selectvar, 'value'), :));
p_aux = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2:dim(2)));
p = Order_current(p_aux-j);

% recupera el índice de la fila (o componente) del vector columna que
% dejará la base según la razón mínima tal que sea positiva
var = get(handles.popupmenu_selectvar3, 'string');
dim = size(var(get(handles.popupmenu_selectvar3, 'value'), :));
q = str2double(var(get(handles.popupmenu_selectvar3, 'value'), 2:dim(2)))-j;
else
% recupera el índice de la variable no básica seleccionada asociada a la
% columna que ingresará a la base
var = VarSelect;
dim = size(var(Jp, :));
p_aux = str2double(var(Jp, 2:dim(2)));
p = Order_current(p_aux-j);

% recupera el índice de la fila (o componente) del vector columna que
% dejará la base según la razón mínima tal que sea positiva
var = VarSelect3;
dim = size(var(Jq, :));
q = str2double(var(Jq, 2:dim(2)))-j;    
end

if ~isnan(q) % si existe algún yij que cumple el criterio de factibilidad
    pivote_process(p, q, handles);
    %AGREGADO(27/12/2016)
    set(handles.listbox_operations, 'string', operations);
    generar_latexsimplexnext()
    %AGREGADO(27/12/2016)
else
    output_cod = 3; output_sol = [];  output_val = []; %#ok<NASGU>
    if handles.Interactive == 1
    operations = char('El conjunto representado por las restricciones es vacío:', 'No hay solución.');
    if strcmp(get(handles.Mode_geo3D, 'Checked'), 'off')
        msgbox(char('El conjunto representado por las restricciones es vacío.'), char('No hay solución'),'modal');
    end
    latex = '';
    if handles.Method == 4
        generar_latexanalisis('\subsection{Análisis de resultados}');
    else
        generar_latexanalisis('\section{Análisis de resultados}');
    end
    handles.latex = [handles.latex, latex];
    guidata(handles.output, handles);
    set(handles.listbox_operations, 'string', operations);
    end
    h = handles; return;
end

if handles.Interactive == 1
% se actualiza la tabla de la interfaz
All_display = Tableau;
spreadsheet = cell(100,100);
spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
set(handles.table_simplexdisplay, 'data', spreadsheet);
%set(handles.table_simplexdisplay, 'data', All_display);
setrowheaders(handles, char('X', 'X', 'Rj', 'Yi0'));
end
handles = calc_variables(handles);
%calc_variables(handles);
if handles.Interactive == 1
handles.latex = [handles.latex, latex];
guidata(handles.output, handles);
end
h = handles;


% --- Función utilitaria: proceso del pivote sobre el elemento Ypq
function pivote_process(p, q, handles)
global Tableau Order_current Dimension operations latex;

%format rat;
% Ecuaciones de pivote
Ypq = Tableau(p, q);
fp = Tableau(p, :); % recupera la fila p actual
fp_t = fp/Ypq;

if handles.Interactive == 1
%AGREGADO(27/12/2016)
%[num, den] = numden(sym(1/Ypq, 'r'));
%fraction = double([num, den]);
pivote_operation = sprintf('Se efectuaron las siguientes operaciones de pivoteo:\r');
operations = pivote_operation;

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, 'El elemento (', num2str(p), ',', num2str(q), ') de la matriz es el pivote: $\bf{', rats(Ypq), '}$.']; 
latex = sprintf('%s \r', latex);
latex = [latex, 'Las operaciones de pivoteo son: \\ \\']; 
latex = sprintf('%s \r', latex);

%if fraction(2) == 1    
    pivote_operation = sprintf('f%d <= (%s) x f%d', p, rats(1/Ypq), p);    
    latex = [latex, '$f_{', num2str(p), '} \leftarrow ', '(', rats(1/Ypq), ')f_{', num2str(p), '}$\\'];     
%else
%    pivote_operation = sprintf('f%d <= (%d/%d) x f%d', p, fraction(1), fraction(2), p);
%    latex = [latex, '$f_', num2str(p), '\leftarrow ', '(\frac{', num2str(fraction(1)), '}{', num2str(fraction(2)), '})f_', num2str(p), '$\\']; 
%end
latex = sprintf('%s \r', latex);
operations = char(operations, pivote_operation);
%AGREGADO(27/12/2016)
end

Tableau(p, :) = fp_t; % se actualiza la fila p en la tabla del Simplex
% se actualizan las demás filas
for i = 1:Dimension(1)
    if i ~= p
        fi = Tableau(i, :);
        Yiq = Tableau(i, q);
        fi = eval(sym(str2num(rats(fi)))) - eval(sym(fp*str2num(rats(Yiq/Ypq)))); %#ok<ST2NM>
        
        if handles.Interactive == 1
        %AGREGADO(27/12/2016)
        %[num, den] = numden(sym(Yiq/Ypq, 'r'));
        %fraction = double([num, den]);
        %if fraction(2) == 1
            pivote_operation = sprintf('f%d  <= f%d - (%s) x f%d', i, i, rats(Yiq/Ypq), p);
            latex = [latex, '$f_{', num2str(i), '} \leftarrow f_{', num2str(i), '}-(', rats(Yiq/Ypq), ')f_{', num2str(p), '}$\\'];
        %else
        %    pivote_operation = sprintf('f%d  <= f%d - (%d/%d) x f%d', i, i, fraction(1), fraction(2), p);
        %    latex = [latex, '$f_', num2str(i), '\leftarrow f_', num2str(i), '-(\frac{', num2str(fraction(1)), '}{', num2str(fraction(2)), '})f_', num2str(p), '$\\'];
        %end
        latex = sprintf('%s \r', latex);
        operations = char(operations, pivote_operation);
        %AGREGADO(27/12/2016)
        end
        
        Tableau(i, :) = fi; % fila p actualizada
    end
end

% intercambian roles las variables
Order_current(1, Order_current(1,:) == p) = Order_current(1,q);
Order_current(1, q) = p;


% --- Se ejecuta el método de ramificación y cota
function branch_and_bound(handles)
global Matrix_problem latex simplex output_cod output_sol output_val; 

Dim = size(handles.Orig_Matrix_problem);

len1 = size(handles.nodelist);
if ~isempty(handles.nodelist)
    if handles.Travel == 1 || handles.Travel == 2
        node = handles.nodelist{len1(2)};
    else
        node = handles.nodelist{1};
    end
    len = length(node);
    handles.Level = len;
    if len >= 0 && len < Dim(2)-1
        if len > 0 && simplex == 0
            aux = handles.Orig_Matrix_problem;
            
            PesoActual = aux(1:Dim(1)-2,1:len)*node;
            subproblem = [aux(1:Dim(1)-1,len+1:Dim(2)-1), ...
                aux(1:Dim(1)-1,end) - [PesoActual;0]; ...
                aux(Dim(1), len+1:Dim(2))];            
            
            %generar_latexspecification();
            %handles.latex = [handles.latex, latex];    
            %generar_latexsimplexbegin()
            %handles.latex = [handles.latex, latex];
            %guidata(handles.output, handles);
            
            handles.gui_Matrix_problem = subproblem;
            Matrix_problem = handles.gui_Matrix_problem;

            handles.isonlylecture = 1;
            handles.setProblem = @setProblem; % se comparte el manejador de función    
            handles.gui_Problem = Problem('LPApp', handles);
                        
            if handles.Travel == 1 || handles.Travel == 2 
                handles.nodelist = handles.nodelist(1:len1(2)-1);
            else
                handles.nodelist = handles.nodelist(2:len1(2));
            end
            simplex = 1;
            set(handles.pushbutton_next, 'String', 'Mejorar solución');
            guidata(handles.output, handles);
        else
            %latex = '';
            if handles.Submethod == 1            
                % se ejecuta el método simplex
                [cod, sol, val, handles] = simplex_primal(handles);
            else
                % se ejecuta el método simplex dual
                [cod, sol, val, handles] = simplex_dual(handles);
            end
            if cod == 0 || cod == 4
                simplex = 1;
            elseif cod == 1
                if ~isempty(node)
                    ZActual = handles.Orig_Matrix_problem(Dim(1)-1,1:len)*node;                    
                    B = ZActual + ceil(val);
                else
                    B = ceil(val);
                    set(handles.popupmenu_selectvar2, 'Enable', 'off');
                end
                s = get(handles.listbox_operations, 'string');
                set(handles.listbox_operations, 'string', ...
                     char(s, ['La cota de la solucion óptima extendida relajada es: ', num2str(B),'.']));
                set(handles.pushbutton_next, 'Enable', 'On');
                if all(ceil(sol) == sol)
                    if B < handles.ZOptimo
                        handles.XOptimo = [node; sol(1:Dim(2)-1-len)];
                        handles.ZOptimo = B;
                        
                        %%Contraer rama                  
                        set(handles.pushbutton_next, 'String', 'Contraer');
                        msgbox(char('La solución extendida ', mat2str(handles.XOptimo'), 'mejora el valor objetivo.'), ...
                            char('Método de ramificación y poda'),'modal'); 
                        latex = sprintf('%s \r', '');                                                
                        latex = [latex, ['La solución extendida ', mat2str(handles.XOptimo'), ' mejora el valor objetivo.']]; 
                        latex = sprintf('%s Nueva cota de valor óptimo: %s. \r', latex, num2str(B));
                        handles.latex = [handles.latex, latex];
                    else
                        %%Podar rama                  
                        set(handles.pushbutton_next, 'String', 'Podar');
                        msgbox(char('La solución extendida ', mat2str([node; sol(1:Dim(2)-1-len)]'), 'no mejora el valor objetivo.'), ...
                            char('Método de ramificación y poda'),'modal'); 
                        latex = sprintf('%s \r', '');                                                
                        latex = [latex, ['La solución extendida ', mat2str([node; sol(1:Dim(2)-1-len)]'), ' no mejora el valor objetivo.']]; 
                        latex = sprintf('%s Valor objetivo alcanzado: %s. \r', latex, num2str(B));
                        handles.latex = [handles.latex, latex];
                    end
                                       
                    if handles.Travel == 1 || handles.Travel == 2 
                        handles.nodelist = handles.nodelist(1:len1(2)-1); 
                    else
                        handles.nodelist = handles.nodelist(2:len1(2));
                    end
                    guidata(handles.output, handles);
                else
                    if B < handles.ZOptimo
                        Bound_node = handles.Orig_Matrix_problem(Dim(1), len+1);
                        for i = 0:Bound_node
                            if handles.Travel == 1 
                                handles.nodelist{len1(2)+i} = [node; Bound_node-i];
                            elseif handles.Travel == 2
                                handles.nodelist{len1(2)+i} = [node; i];
                            elseif handles.Travel == 3
                                handles.nodelist{len1(2)+i+1} = [node; Bound_node-i];
                            else
                                handles.nodelist{len1(2)+i+1} = [node; i];
                            end                            
                        end
                        if handles.Travel == 3 || handles.Travel == 4
                            handles.nodelist = handles.nodelist(2:len1(2)+i+1);
                        end
                        %%Ramificar                        
                        set(handles.pushbutton_next, 'String', 'Ramificar');
                        msgbox(char('Alguna solución extendida de ', mat2str(node'), ' puede mejorar el valor objetivo.'), ...
                            char('Método de ramificación y poda'),'modal');
                        latex = sprintf('%s \r', '');                                                
                        latex = [latex, ['Alguna solución extendida de ', mat2str(node'), ' puede mejorar el valor objetivo.']]; 
                        latex = sprintf('%s Cota alcanzada: %s. \r', latex, num2str(B));
                        handles.latex = [handles.latex, latex];
                        guidata(handles.output, handles);
                    else
                        %%Podar rama                        
                        set(handles.pushbutton_next, 'String', 'Podar');
                        msgbox(char('Cualquier solución extendida de ', mat2str(node'), ' no mejora el valor objetivo.'), ...
                            char('Método de ramificación y poda'),'modal'); 
                        latex = sprintf('%s \r', '');                                                
                        latex = [latex, ['Cualquier solución extendida de ', mat2str(node'), ' no mejora el valor objetivo.']]; 
                        latex = sprintf('%s Cota alcanzada: %s. \r', latex, num2str(B));
                        handles.latex = [handles.latex, latex];
                        if handles.Travel == 1 || handles.Travel == 2
                            handles.nodelist = handles.nodelist(1:len1(2)-1);
                        else
                            handles.nodelist = handles.nodelist(2:len1(2));
                        end
                        guidata(handles.output, handles);
                    end
                end
                simplex = 0;
                handles.ispenalties = 0; %% por si se ocupó el metodo de penalizaciones
                guidata(handles.output, handles);
            elseif cod == 2
                Bound_node = handles.Orig_Matrix_problem(Dim(1), len+1);
                for i = 0:Bound_node
                    if handles.Travel == 1 
                        handles.nodelist{len1(2)+i} = [node; Bound_node-i];
                    elseif handles.Travel == 2
                        handles.nodelist{len1(2)+i} = [node; i];
                    elseif handles.Travel == 3
                        handles.nodelist{len1(2)+i+1} = [node; Bound_node-i];
                    else
                        handles.nodelist{len1(2)+i+1} = [node; i];
                    end                            
                end
                if handles.Travel == 3 || handles.Travel == 4
                    handles.nodelist = handles.nodelist(2:len1(2)+i+1);
                end
                %%Ramificar
                simplex = 0;
                set(handles.pushbutton_next, 'Enable', 'On');
                set(handles.pushbutton_next, 'String', 'Ramificar');
                msgbox(char('Cualquier solución extendida de ', mat2str(node'), ' no es acotada.'), ...
                    char('Método de ramificación y poda'),'modal');
                latex = sprintf('%s \r', '');                                                
                latex = [latex, ['Cualquier solución extendida de ', mat2str(node'), ' no es acotada.']];                
                handles.latex = [handles.latex, latex];
                handles.ispenalties = 0; %% por si se ocupó el metodo de penalizaciones
                guidata(handles.output, handles);
            elseif cod == 3
                %%Podar rama
                simplex = 0;
                set(handles.pushbutton_next, 'Enable', 'On');
                set(handles.pushbutton_next, 'String', 'Podar');
                msgbox(char('Cualquier solución extendida de ', mat2str(node'), ' no es factible.'), ...
                    char('Método de ramificación y poda'),'modal');
                latex = sprintf('%s \r', '');                                                
                latex = [latex, ['Cualquier solución extendida de ', mat2str(node'), ' no es factible.']];                
                handles.latex = [handles.latex, latex];
                if handles.Travel == 1 || handles.Travel == 2
                    handles.nodelist = handles.nodelist(1:len1(2)-1);
                else
                    handles.nodelist = handles.nodelist(2:len1(2));
                end
                handles.ispenalties = 0; %% por si se ocupó el metodo de penalizaciones
                guidata(handles.output, handles);
            end        
        end
    else
        aux = handles.Orig_Matrix_problem;
        PesoActual = aux(1:Dim(1)-2,1:Dim(2)-1)*node;
        if all(PesoActual <= aux(1:Dim(1)-2,Dim(2)))
            ZActual = aux(Dim(1)-1,1:Dim(2)-1)*node;
            if ZActual < handles.ZOptimo
                handles.XOptimo = node;
                handles.ZOptimo = ZActual;
                msgbox(char('La solución ', mat2str(node'), 'mejora el valor objetivo.'), ...
                    char('Método de ramificación y poda'),'modal');
                latex = sprintf('%s \\subsection{Verificación de solución en el nivel %i} \r', '', handles.Level);
                latex = [latex, ['La solución extendida ', mat2str(handles.XOptimo'), ' mejora el valor objetivo.']]; 
                latex = sprintf('%s Nueva cota de valor óptimo: %s. \r', latex, rats(handles.ZOptimo));
                handles.latex = [handles.latex, latex];
                
                set(handles.pushbutton_next, 'Enable', 'On');
                set(handles.pushbutton_next, 'String', 'Contraer');
            else
                msgbox(char('La solución ', mat2str(node'), 'no mejora el valor objetivo.'), ...
                    char('Método de ramificación y poda'),'modal');
                latex = sprintf('%s \\subsection{Verificación de solución en el nivel %i} \r', '', handles.Level);
                latex = [latex, ['La solución extendida ', mat2str(node'), ' no mejora el valor objetivo.']]; 
                latex = sprintf('%s Valor objetivo alcanzado: %s. \r', latex, rats(ZActual));
                handles.latex = [handles.latex, latex];
                set(handles.pushbutton_next, 'Enable', 'On');
                set(handles.pushbutton_next, 'String', 'Podar');
            end
        else
            msgbox(char('La solución ', mat2str(node'), 'no es factible.'), ...
                    char('Método de ramificación y poda'),'modal');
            latex = sprintf('%s \\subsection{Verificación de solución en el nivel %i} \r', '', handles.Level);
            latex = [latex, ['Cualquier solución extendida de ', mat2str(node'), ' no es factible.']];                
            handles.latex = [handles.latex, latex];
            set(handles.pushbutton_next, 'Enable', 'On');
            set(handles.pushbutton_next, 'String', 'Podar');
        end
        
        if handles.Travel == 1 || handles.Travel == 2
            handles.nodelist = handles.nodelist(1:len1(2)-1);
        else
            handles.nodelist = handles.nodelist(2:len1(2));
        end
        %simplex = 1;        
        guidata(handles.output, handles);
        %%Evaluar hojas
    end
else
    handles.Level = 0;
    set(handles.pushbutton_next, 'Enable', 'Off');
    set(handles.pushbutton_next, 'String', 'Mejorar solucion');
    if handles.ZOptimo ~= Inf
        output_cod = 1; output_sol = handles.XOptimo; output_val = handles.ZOptimo;        
        msgbox(char('La solución óptima es ', mat2str(handles.XOptimo'), ' y su valor objetivo es ', ...
            [mat2str(handles.ZOptimo'),'.']), char('Método de ramificación y poda'),'modal');
        s = get(handles.listbox_operations, 'string');
        set(handles.listbox_operations, 'string', ...
             char(s, 'La solución óptima es ', mat2str(handles.XOptimo'), ' y su valor objetivo es: ', ...
            [rats(handles.ZOptimo'),'.']));
    else
        output_cod = 3; output_sol = []; output_val = [];                
        msgbox(char('El problema no tiene soluciones factibles.'), ...
                    char('Método de ramificación y poda'),'modal');
    end
    latex = '';
    generar_latexanalisis('\section{Análisis de resultados}');
    handles.latex = [handles.latex, latex];
    %%Se ha recorrido todo el árbol en LIFO
    guidata(handles.output, handles);
end
    

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
spreadsheet = cell(100,100);
spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(Aux_Tableau);
set(handles.table_simplexdisplay, 'data', spreadsheet);
All_display = Aux_Tableau;

if handles.Method == 1
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
    
    spreadsheet = cell(100,100);
    spreadsheet(1:Dimension(1), 1:Dimension(2)+1) = num2cell(All_display);
    set(handles.table_simplexdisplay, 'data', spreadsheet);
    %set(handles.table_simplexdisplay, 'data', All_display);
elseif handles.Method == 2
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
    
    spreadsheet = cell(100,100);
    spreadsheet(1:Dimension(1)+1, 1:Dimension(2)) = num2cell(All_display);
    set(handles.table_simplexdisplay, 'data', spreadsheet);
    %set(handles.table_simplexdisplay, 'data', All_display);
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


% --- actualiza la especificación del problema
function setProblem(Problem, handles)
% Problem especificación del problema
% handles estructura de manejadores y datos de usuario
global Dimension Matrix_problem Tableau Order_current Order_initial Node_current Solution T_Tableau Num_assignation latex;

% se actualizan algunas variables globales previamente descritas

Matrix_problem = Problem;
Dimension = size(Matrix_problem);
n = Dimension(2)-1;
m = Dimension(1)-1;

if handles.Interactive == 1
%AGREGADO(27/12/2016)
set(handles.listbox_operations, 'string', '');
%AGREGADO(27/12/2016)

if handles.Method == 2 || handles.Submethod == 2
    set(handles.text_selectvar, 'string', 'Variable saliente');
    set(handles.text_selectvar3, 'string', 'Variable entrante');
else
    set(handles.text_selectvar, 'string', 'Variable entrante');
    set(handles.text_selectvar3, 'string', 'Variable saliente');
end
end

if handles.Method ~= 3
    Order_initial = handles.Order_initial;
    Order_current = Order_initial;
    Tableau = Matrix_problem; 

    % se convierte la última en términos de Rj =cj - zj
    [Y, I] = sort(Order_initial); %#ok<ASGLU>
    Tableau_sorted = Matrix_problem(1:m+1, I);
    c = Tableau_sorted(end,:);
    z = Tableau_sorted(end, 1:(Dimension(1)-1))*Tableau_sorted(1:(Dimension(1)-1),:);

    R = eval(sym(str2num(rat(c)))) - eval(sym(str2num(rat(z)))); %#ok<ST2NM>
    Tableau(end, I) = R;
    if handles.Interactive == 1
    %AGREGADO(27/12/2016)
    set(handles.listbox_operations, 'string', ...
        char('La fila de costos unitarios se ha convertido en costos relativos:', 'r <= c - z'));
    %generar_latexsimplexbegin(handles);
    %AGREGADO(27/12/2016)
    end
else
    set(handles.Next_value, 'Enable', 'on');
    Node_current = [0,0];
    Num_assignation = 0;
    Solution = zeros(Dimension(1)+Dimension(2)-3, 1, 3);
    T_Tableau = Matrix_problem;
end
if handles.Interactive == 1
% se actualiza la tabla de la interfaz
All_display = zeros(m+1, n+1);
if handles.Method ~= 3
    % se rotulan las columnas y las filas de manera correspondiente
    setTableauTags(handles, char('X', 'X', 'Rj', 'Yi0'));
    All_display(:, 1:Dimension(2)) = Tableau;
else
    % se rotulan las columnas y las filas de manera correspondiente
    setTableauTags(handles, char('O', 'D', 'Rj', 'Ci'));
    All_display(:, 1:Dimension(2)) = T_Tableau;
end
%AGREGADO 1/01/2017
spreadsheet = cell(100,100);
spreadsheet(1:m+1, 1:Dimension(2)) = num2cell(All_display);
set(handles.table_simplexdisplay, 'data', spreadsheet);

%set(handles.Restriction_nonnegativity, 'Enable', 'off');
%set(handles.Uncut_planes, 'Enable', 'off');
%set(handles.Saveimage, 'Enable', 'off');

set(handles.Generate_latex, 'Enable', 'on');
end

%set(handles.table_simplexdisplay, 'data', All_display); 
%AGREGADO 1/01/2017
% se calculan las nuevas variables no básicas si las hay
calc_variables(handles);

if handles.Interactive == 1
handles.latex = [handles.latex, latex];
guidata(handles.output, handles); 
if Dimension(2)- Dimension(1) == 3 && handles.istwophases == 0
    set(handles.Mode_geo3D, 'Enable', 'on'); 
    %set(handles.pushbutton_asignall, 'Enable', 'on')
else
    set(handles.Mode_geo3D, 'Enable', 'off');
    set(handles.Saveimage, 'Enable', 'off');
end
end


% ----se rotulan las columnas y las filas de manera
% correspondiente
function setTableauTags(handles, headers)
global Dimension;

colName = cell(Dimension(2)+1, 1);
j = 0;
if handles.Method == 4    
    j = handles.Level;
end
for i = (1+j):Dimension(2)-1+j
    colName(i-j) = cellstr(strcat(headers(2,:),num2str(i)));
end
colName(Dimension(2)) = cellstr(headers(4,:)); 
if handles.Method == 1 || handles.Submethod == 1
    colName(Dimension(2)+1) = cellstr('Yi0/Yij');
    %colName(Dimension(2)+2) = cellstr('lambdai');    
end
set(handles.table_simplexdisplay, 'columnname', colName);

setrowheaders(handles, headers);

if handles.Method == 1 || handles.Submethod == 1
    colFormat=cell(1,Dimension(2)+1);    
else
    colFormat=cell(1,Dimension(2)+1);
end
% se especifica el formato de salida de los datos
colFormat(1,:) = cellstr('rat');

set(handles.table_simplexdisplay, 'columnformat', colFormat);

% -------
function setrowheaders(handles, headers)
global Dimension Order_current;

[Z, I] = sort(Order_current); %#ok<ASGLU>
rowName=cell(1,Dimension(1));
j = 0;
if handles.Method == 4    
    j = handles.Level;
end
for i = 1:(Dimension(1)-1)
    if handles.Method ~= 3
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(I(i)+j)));
    else
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end
end
rowName(Dimension(1)) = cellstr(headers(3,:));

if handles.Method == 2 || handles.Submethod == 2
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

[filename, Y, Z] = uigetfile('*.mat', 'Selecciona un nombre de archivo válido', './Problemas'); %#ok<NASGU>
if filename ~= 0
    S = load([Y, filename]);
    
    %AGREGADO(27/12/2016)
    %handles.istwophases = 0;
    %handles.whatphase = 1;
    %handles.isinit_secondphase = 0;
    %handles.latex = '';
    if (isfield(S, 'Method'))
        handles.Newmethod = S.Method;
        handles.setProblem = @setProblem;    
        handles.gui_Matrix_problem = S.Matrix_problem;
        handles.gui_Problem = Problem('LPApp', handles);
        guidata(handles.output, handles);
        %handles.latex = '';
    else
        errordlg(char('El problema no es compatible con la versión actual'), 'Problema no compatible');
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
    if handles.istwophases == 1 || handles.ispenalties == 1
        Matrix_problem = handles.Orig_Matrix_problem;
    end        
elseif handles.Method == 2
    Method = 2; %#ok<NASGU>
    if handles.istwophases == 1 || handles.ispenalties == 1
        Matrix_problem = handles.Orig_Matrix_problem;
    end
elseif handles.Method == 4
    Method = 4; %#ok<NASGU>
    Matrix_problem = handles.Orig_Matrix_problem;
end
%AGREGADO(27/12/2016)

mkdir('Problemas');
uisave(cellstr(char('Matrix_problem', 'Method')), './Problemas/LProblem1');
%uisave('Matrix_problem', 'LProblem');



% ----------
function setProblemAndTableau(Problem, NewTableau, handles)
global Matrix_problem Tableau;

dim = size(NewTableau);
spreadsheet = cell(100,100);
spreadsheet(1:dim(1), 1:dim(2)) = num2cell(NewTableau);
set(handles.table_simplexdisplay, 'data', spreadsheet);
%set(handles.table_simplexdisplay, 'data', NewTableau);

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
global Num_assignation Dimension latex handles_surf handles_norm;

if strcmp(get(handles.Mode_geo3D, 'Checked'),'off')
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
else
    plano = get(handles.popupmenu_selectvar2, 'value');
    if strcmp(get(handles_surf(plano), 'Visible'),'off')
        set(handles_surf(plano), 'Visible','on');
        set(handles_norm(plano), 'Visible','on');
        set(handles.pushbutton_asignall, 'string','Quitar');
    else
        set(handles_surf(plano), 'Visible','off');
        set(handles_norm(plano), 'Visible','off');
        set(handles.pushbutton_asignall, 'string','Poner');
    end
end


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
    var = get(handles.popupmenu_selectvar2, 'string'); % se recupera la variable no básica seleccionada
    V = var(get(handles.popupmenu_selectvar2, 'value'), 1); % se recupera variable seleccionada 
    I = str2double(var(get(handles.popupmenu_selectvar2, 'value'), 2)); % se recupera el índice de variable seleccionada 
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
        %all_nan = 0;
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
                    %all_nan = 0;
                end
            %else
            %    all_nan = 1;
            end
        end
        %if all_nan == 1
        %    [f, c] = find(T_VarType_Aux);
        %    dim = size(f);
        %    if dim(1) == 1
        %        [f, c] = find(T_VarType(:, c));
        %        T_VarType_Aux(f(1), c(1)) = 1;
        %    end
        %end
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
        msgbox(char('Hay costos reducidos negativos. La solución actual se puede mejorar.'),char('Cálculo de costos reducidos.'),'modal');
        set(handles.popupmenu_selectvar3, 'string', char(' ', ' '));
        set(handles.popupmenu_selectvar3, 'value', 1);
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
        msgbox(char('La solución actual es óptima. No hay costos reducidos negativos.'),char('Solución óptima encontrada.'),'modal');
        generar_latexnextsolution('\section{Análisis de Resultados}');
        set_environment('end', handles);
    end
end
All_display = T_Tableau;
spreadsheet = cell(100,100);
spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
set(handles.table_simplexdisplay, 'data', spreadsheet);
%set(handles.table_simplexdisplay, 'data', All_display);


% ------------------------
function calc_nextciclevar(handles, Num_cassignation)
global Node_current Num_assignation T_Tableau T_VarType_Aux T_VarType Node_NBV Dimension Node_Cicle...
    Empty_dimension Matrix_problem Solution_initial Solution Solution_change Minimo;

if handles.Method == 1
    colFormat=cell(1,Dimension(2)+1);    
else
    colFormat=cell(1,Dimension(2));
end
colFormat(1,:) = cellstr('rat');
set(handles.table_simplexdisplay, 'columnformat', colFormat);

if Num_assignation == 0 || Num_cassignation == 0
    Num_assignation = 0;
    if all(T_VarType_Aux == 0)
        T_VarType_Aux = T_VarType;
    end
    
    Node_Cicle = zeros(Dimension(1)+Dimension(2)-2, 1, 3);    
    Solution_change = zeros(1,Dimension(1)+Dimension(2)-3);
    T_Tableau = zeros(Dimension(1), Dimension(2));
    T_Tableau(end, 1:Dimension(2)) = Matrix_problem(end, 1:Dimension(2));
    T_Tableau(1:Dimension(1), end) = Matrix_problem(1:Dimension(1), end);
    var = get(handles.popupmenu_selectvar, 'string'); % se recupera la variable no básica seleccionada
    I = str2double(var(get(handles.popupmenu_selectvar, 'value'), 2)); % se recupera variable seleccionada 
    J = str2double(var(get(handles.popupmenu_selectvar, 'value'), 3)); % se recupera el índice de variable seleccionada 
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
        Empty_dimension = 0; % Se ha seleccionado una variable básica en la fila actual
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
            Empty_dimension = 1; % Se ha seleccionado una variable básica en la columna actual
        end
    end        
else
    set(handles.popupmenu_selectvar, 'enable', 'off');
    Node_current = [Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)];
    if Empty_dimension == 1 % Se ha seleccionado una variable básica en la columna actual
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
                    Empty_dimension = 0; % Se ha seleccionado una variable básica en la fila actual
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
                       if Num_assignation-1 >= Num_cassignation %%%%OJO
                            Num_assignation = Dimension(1)+Dimension(2)-3;%%%%
                       end%%%%
                       Empty_dimension = 0; % Se ha seleccionado una variable básica en la fila actual
                       break;%%%%
                    else
                        T_Tableau(Node_current(1), ind1(i)) = 0;              
                        if ~(i == Lsup)
                            Empty_dimension = 0; % Se ha seleccionado una variable básica en la fila actual
                        else
                            Empty_dimension = 1; % Se ha seleccionado una variable básica en la columna actual
                        end
                        %MODIFICADO(27/12/2016)
                        %Node_Cicle(Num_assignation, 1, 3) = 0; %%%%%OJO
                        %Num_assignation = Num_assignation - 1;                        
                        %MODIFICADO(27/12/2016)
                    end                     
                end
            end  
        else
            Empty_dimension = 0; % Se ha seleccionado una variable básica en la fila actual
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
            Node_current = [Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)];
            ind1 = find(T_VarType_Aux(:, Node_current(2)));
            ind_aux = find(T_Tableau(1:Dimension(1)-1, Node_current(2))); %%%%OJO
            
            if isempty(ind1) || ind1(1) < Node_current(1) || length(ind_aux) > 1
                T_Tableau(Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)) = 0;
                Node_Cicle(Num_cassignation, 1, 3) = 0;
                for l = 1:Dimension(1)+Dimension(2)-3;
                    if Solution(l, 1, 1) == Node_Cicle(Num_cassignation, 1, 1) && Solution(l, 1, 2) == Node_Cicle(Num_cassignation, 1, 2)
                        Solution_change(l) = 0;
                        break;
                    end
                end
                
                calc_nextciclevar(handles,Num_cassignation-1);
            else
                Num_assignation = Num_cassignation;
                calc_nextciclevar(handles,Num_cassignation);
            end            
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
                    Empty_dimension = 1; % Se ha seleccionado una variable básica en la columna actual
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
                       if Num_assignation-1 >= Num_cassignation %%%OJO
                            Num_assignation = Dimension(1)+Dimension(2)-3;
                       end
                       Empty_dimension = 1; % Se ha seleccionado una variable básica en la columna actual
                       break;
                    else
                        %%%msgbox(char('El camino no es un ciclo. Se retornará al punto anterior.'), 'Camino en encrucijada.','modal');
                        T_Tableau(ind1(i), Node_current(2)) = 0;
                        if ~(i == Lsup)
                            Empty_dimension = 1; % Se ha seleccionado una variable básica en la columna actual
                        else
                            Empty_dimension = 0; % Se ha seleccionado una variable básica en la fila actual
                        end
                        %MODIFICADO(27/12/2016)
                        %Node_Cicle(Num_assignation, 1, 3) = 0;
                        %Num_assignation = Num_assignation -1;                        
                        %MODIFICADO(27/12/2016)
                    end
                end                                
            end  
        else
           Empty_dimension = 1; % Se ha seleccionado una variable básica en la columna actual
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
           Node_current = [Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)];
           ind2 = find(T_VarType_Aux(Node_current(1), :));
           ind_aux = find(T_Tableau(Node_current(1), 1:Dimension(2)-1)); %%%%%%OJO
           
           if isempty(ind2) || ind2(1) < Node_current(2) || length(ind_aux) > 1
                T_Tableau(Node_Cicle(Num_cassignation, 1, 1), Node_Cicle(Num_cassignation, 1, 2)) = 0;
                Node_Cicle(Num_cassignation, 1, 3) = 0;
                for l = 1:Dimension(1)+Dimension(2)-3;
                    if Solution(l, 1, 1) == Node_Cicle(Num_cassignation, 1, 1) && Solution(l, 1, 2) == Node_Cicle(Num_cassignation, 1, 2)
                        Solution_change(l) = 0;
                        break;
                    end
                end  
                calc_nextciclevar(handles,Num_cassignation-1);
           else
                Num_assignation = Num_cassignation;
                calc_nextciclevar(handles,Num_cassignation);
           end
           return;
        end
    end
end

if Num_assignation == Dimension(1)+Dimension(2)-3
    msgbox(char('Se ha encontrado un ciclo de redistribución.'),char('Cálculo de ciclo de redistribución.'),'modal');
    Solution_initial = 0;
    Node_current = [0, 0];
    set_environment('next_newassign', handles);
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
spreadsheet = cell(100,100);
spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
set(handles.table_simplexdisplay, 'data', spreadsheet);
%set(handles.table_simplexdisplay, 'data', All_display);


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
    spreadsheet = cell(100,100);
    spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
    set(handles.table_simplexdisplay, 'data', spreadsheet);
    %set(handles.table_simplexdisplay, 'data', All_display); 
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

% --- se ejecuta al seleccionar la opción "Nuevo poblema/Simplex primal" del menú "Archivo"
function Simplex_primal_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_primal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Newmethod = 1;
init_method(handles);

% --- se ejecuta al seleccionar la opción "Nuevo poblema/Simplex dual" del menú "Archivo"
function Simplex_dual2_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_dual2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Newmethod = 2;
init_method(handles);

% --- se ejecuta al seleccionar la opción "Nuevo poblema/Simplex de transporte" del menú "Archivo"
function Simplex_transportation2_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_transportation2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Newmethod = 3;
init_method(handles);

function init_method(handles)
% ventana de dialogo en donde se solicita las dimensiones del problema

if handles.Newmethod == 1 || handles.Newmethod == 2 || handles.Newmethod == 4
    prompt = {char('Número de ecuaciones:'),char('Número de variables:')};
    def = {'1','2'};
elseif handles.Newmethod == 3
    prompt = {char('Número de origenes:'),char('Número de destinos:')};
    def = {'2','2'}; 
end

dlg_title = 'Dimensiones del problema';
num_lines = 1;

answer = inputdlg(prompt,dlg_title,num_lines,def);

% se verifica que las dimensiones del problema sean consistentes (m < n)
if ~isempty(answer)    
    user_entry1 = str2double(cellstr(answer{1}));
    user_entry2 = str2double(cellstr(answer{2}));

    while (1)
        if (isnan(user_entry1) || user_entry1 < 1) || (isnan(user_entry2) || user_entry2 < 2) || (user_entry1 > user_entry2 && ...
                handles.Newmethod ~= 3 && handles.Newmethod ~= 4) || (user_entry1 > 96 || user_entry2 > 96) || ...
                (user_entry1 > user_entry2 + user_entry1 && handles.Newmethod == 4)
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if isempty(answer)
                return;
            end
            user_entry1 = str2double(cellstr(answer{1}));
            user_entry2 = str2double(cellstr(answer{2}));
        else
            break;
        end
    end
    if handles.Newmethod == 1 || handles.Newmethod == 2 || handles.Newmethod == 3
        handles.gui_Matrix_problem = zeros(user_entry1+1, user_entry2+1);
    elseif handles.Newmethod == 4
        handles.gui_Matrix_problem = zeros(user_entry1+2, user_entry2+1);
        handles.gui_Matrix_problem(user_entry1+2, 1:user_entry2) = ones(1, user_entry2);
    end
    % se abre la ventana para introducir la especificación del problema
    handles.setProblem = @setProblem; % se comparte el manejador de función    
    handles.gui_Problem = Problem('LPApp', handles);
    guidata(handles.output, handles);
    %handles.latex = '';
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

latexfile = '% OJO: El documento está guardado en codificación de caracteres ISO-8859-1';
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '% ¡¡Compilar de LaTex->DVI->pdf!!'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\documentclass[letterpaper]{article}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage[textwidth=16cm,textheight=24cm]{geometry}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '%configuración para idioma español'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage[latin1]{inputenc}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage{amsfonts}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage{amsmath}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage[spanish,activeacute,es-noindentfirst]{babel}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage{graphicx}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\usepackage{epstopdf}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\title{Archivo generado por Asistente Simplex + 2017}'];
latexfile = sprintf('%s \r', latexfile);
latexfile = [latexfile, '\author{Profesor M.Sc. Porfirio Armando Rodríguez}'];
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
mkdir('LaTex');
[FileName,PathName,FilterIndex] = uiputfile({'*.tex','Archivo de imagen (*.tex)'},'Guardar informe LaTex',...
          './LaTex/ReporteLaTex1.tex'); %#ok<ASGLU,NASGU>
handles.latexfile = ['./LaTex/', FileName];

if FileName ~= 0
    dlmwrite(handles.latexfile, latexfile, '');
end

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
            %latex = [latex, num2str(Tableau(i, j))]; %#ok<*AGROW>
            latex = [latex, rats(Tableau(i, j))];
        else
            %latex = [latex,'& ', num2str(Tableau(i, j))];
            latex = [latex,'& ', rats(Tableau(i, j))];
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
global latex output_cod output_sol output_val;

latex = sprintf('%s \r', latex);
latex = [latex, title];
if output_cod == 1 || output_cod == 4 || output_cod == 0
    latex = sprintf('%s La solución óptima es: ', latex);
    latex = [latex, '$[\textbf{', rats(output_sol'),'}]$.'];
    latex = sprintf('%s \r El valor óptimo es: ', latex);
    %latex = [latex, '$\bf{', num2str(-Tableau(end, end)),'}$.'];
    latex = [latex, '$\bf{', rats(output_val),'}$. '];
    latex = sprintf('%s \r', latex);
elseif output_cod == 2    
    latex = [latex, 'El conjunto representado por las restricciones no es acotado: {\bf El valor objetivo decrece sin límite.} \\']; 
    latex = sprintf('%s \r', latex);
elseif output_cod == 3
    latex = sprintf('%s El conjunto representado por las restricciones es vacío. No tiene soluciones factibles. \r', latex);
end


function generar_latexspecification()
global Matrix_problem latex;

dim = size(Matrix_problem);

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, '\section{Segunda Fase del Método Simplex}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\subsection{Especificación equivalente del problema original}'];
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
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{\min \ '];
inicial = 0;
for j = 1:dim(2) - 1
    if Matrix_problem(dim(1), j) == 1
        coef = '';
    elseif Matrix_problem(dim(1), j) == -1
        coef = '-';
    else
        coef = rats(Matrix_problem(dim(1), j));
    end
    if Matrix_problem(dim(1), j) > 0
        if inicial == 1
            %latex = [latex, ...
            %'+', num2str(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
            latex = [latex, '+', coef, 'x_{',num2str(j),'}'];
        else 
            %latex = [latex, ...
            %num2str(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
            latex = [latex, coef, 'x_{',num2str(j),'}'];
        end
        inicial = 1;
    elseif Matrix_problem(dim(1), j) < 0
        %latex = [latex, ...
        %num2str(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
        latex = [latex, coef, 'x_{',num2str(j),'}'];
        inicial = 1;
    end
    if j == dim(2) - 1
        latex = [latex, '} \\'];
    end 
end
latex = sprintf('%s \r', latex);
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{\text{sujeto  a}} \\'];
latex = sprintf('%s \r', latex);
inicial = 0;
for i = 1:dim(1) - 1
    for j = 1:dim(2) - 1
        if Matrix_problem(i, j) == 1
            coef = '';
        elseif Matrix_problem(i, j) == -1
            coef = '-';
        else
            coef = rats(Matrix_problem(i, j));
        end
        if Matrix_problem(i, j) > 0
            if inicial == 1
                %latex = [latex, ...
                %'+', num2str(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
                latex = [latex, '+', coef, 'x_{',num2str(j), '} & '];
            else
                %latex = [latex, ...
                %num2str(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
                latex = [latex, coef, 'x_{',num2str(j), '} & '];
            end
            inicial = 1;
        elseif Matrix_problem(i, j) < 0
            %latex = [latex, ...
            %num2str(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
            latex = [latex, coef, 'x_{',num2str(j), '} & '];
            inicial = 1;
        else
            latex = [latex, ' & '];
        end           
    end
    inicial = 0;
    %latex = [latex, ' = & ', num2str(Matrix_problem(i, dim(2))), '\\'];
    latex = [latex, ' = & ', rats(Matrix_problem(i, dim(2))), '\\'];
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

% se convierte la última en términos de Rj =cj - zj
[Y, I] = sort(Order_initial); %#ok<ASGLU>
Tableau_sorted = Matrix_problem(:, I);
c = Tableau_sorted(end,:);
z = Tableau_sorted(end, 1:(dim(1)-1))*Tableau_sorted(1:(dim(1)-1),:);

R = c - z;
Tableau(end, I) = R;

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, '\subsection{Desarrollo del método simplex}'];
latex = sprintf('%s \r', latex);
latex = sprintf('%s \r La tabla inicial del método simplex es la siguiente: \r', latex);
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
            %latex = [latex, num2str(Tableau(i, j))]; %#ok<AGROW>
            latex = [latex, rats(Tableau(i, j))]; %#ok<AGROW>
        else
            %latex = [latex,'& ', num2str(Tableau(i, j))]; %#ok<AGROW>
            latex = [latex,'& ', rats(Tableau(i, j))]; %#ok<AGROW>
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
                latex = [latex, '&  \multicolumn{2}{c}{$D_{', num2str(j), '}$} ']; %#ok<AGROW>
            else
                latex = [latex, '&  \multicolumn{1}{c}{$C_i$} \\']; %#ok<AGROW>                           
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    elseif i < dim(1)+1
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{$O_{', num2str(i-1),'}$} ']; %#ok<AGROW>
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
                latex = [latex, '&  \multicolumn{1}{c}{',rats(Matrix_solution(i-1, j-1)),'} & \multicolumn{1}{c|}{} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{',rats(Matrix_problem(i-1, j-1)),'} \\']; %#ok<AGROW>                
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
                    latex = [latex, '&  \multicolumn{1}{c}{', rats(Matrix_problem(i-1, j-1)), '} &  \multicolumn{1}{c|}{} ']; %#ok<AGROW>
                else
                    latex = [latex, '&  \multicolumn{1}{c}{', rats(Matrix_problem(i-1, j-1)), '} &  \multicolumn{1}{c|}{} '];
                end
            else
                latex = [latex, '\multicolumn{1}{c|}{$R_j$} ']; %#ok<AGROW>                           
            end            
        end
        latex = [latex, '&  \multicolumn{1}{c|}{',rats(Matrix_solution(i-1, j)),'} \\'];
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
latex = [latex, '\subsection{Encontrando los costos reducidos asociados a la solución actual}'];
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
                latex = [latex, '&  \multicolumn{2}{c}{$D_{', num2str(j), '}$} ']; %#ok<AGROW>
            else
                latex = [latex, '&  \multicolumn{1}{c}{$U_i$} \\']; %#ok<AGROW>                           
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    elseif i < dim(1)+1
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{$O_{', num2str(i-1),'}$} ']; %#ok<AGROW>
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
                latex = [latex, '&  \multicolumn{1}{c}{',rats(T_Tableau(i-1, j-1)),'} & \multicolumn{1}{c|}{} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{',rats(T_Tableau(i-1, j-1)),'} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    else
        for j = 1:dim(2)
            if j > 1        
                latex = [latex, '&  \multicolumn{1}{c}{', rats(T_Tableau(i-1, j-1)), '} &  \multicolumn{1}{c|}{} ']; %#ok<AGROW>
            else
                latex = [latex, '\multicolumn{1}{c|}{$V_j$} ']; %#ok<AGROW>                           
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
latex = [latex, '\subsection{Encontrando un ciclo de redistribución}'];
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
                latex = [latex, '&  \multicolumn{2}{c}{$D_{', num2str(j), '}$} ']; %#ok<AGROW>
            else
                latex = [latex, '&  \multicolumn{1}{c}{$C_i$} \\']; %#ok<AGROW>                           
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    elseif i < dim(1)+1
        for j = 1:dim(2)+1
            if j == 1        
                latex = [latex, '\multicolumn{1}{c|}{$O_{', num2str(i-1),'}$} ']; %#ok<AGROW>
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
                latex = [latex, '&  \multicolumn{1}{c}{',rats(T_Tableau(i-1, j-1)),'} & \multicolumn{1}{c|}{} ']; %#ok<AGROW>                
            else
                latex = [latex, '& \multicolumn{1}{c|}{',rats(T_Tableau(i-1, j-1)),'} \\']; %#ok<AGROW>                
            end
        end
        latex = sprintf('%s \r', latex);
        for j = 2:2*dim(2)
            latex = [latex, '\cline{', num2str(j), '-', num2str(j), '} ']; %#ok<AGROW>            
        end
    else
        for j = 1:dim(2)
            if j > 1        
                latex = [latex, '&  \multicolumn{1}{c}{', rats(Matrix_problem(i-1, j-1)), '} &  \multicolumn{1}{c|}{} ']; %#ok<AGROW>
            else
                latex = [latex, '\multicolumn{1}{c|}{$R_j$} ']; %#ok<AGROW>                           
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

% abre ventana para modificar la especificación del problema
handles.setProblem = @setProblem; % se comparte manejador de función
%handles.gui_Matrix_problem = Matrix_problem; % se comparte la especificación de problema actual

%AGREGADO(27/12/2016)
handles.latex = '';
if handles.istwophases == 1
    handles.gui_Matrix_problem = handles.Orig_Matrix_problem;
    %handles.isinit_secondphase = 0;
    %handles.whatphase = 1;
elseif handles.Method == 4 || handles.ispenalties == 1
    %m = handles.Orig_Matrix_problem;
    handles.gui_Matrix_problem = handles.Orig_Matrix_problem;
else
    handles.gui_Matrix_problem = Matrix_problem; % se comparte la especificación de problema actual
end
%AGREGADO(27/12/2016)

handles.Newmethod = handles.Method;
handles.gui_Problem = Problem('LPApp', handles);
guidata(handles.output, handles);



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
 
spreadsheet = cell(100,100);
spreadsheet(1:Dimension(1), 1:Dimension(2)) = All_display;
set(handles.table_simplexdisplay, 'data', spreadsheet);
%set(handles.table_simplexdisplay, 'data', All_display);




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
    Solution_initial T_VarType NewSolution handles_surf;

if strcmp(get(handles.Mode_geo3D, 'Checked'),'off')
    if handles.Method == 4 && handles.Level == 0
        handles.Travel = get(handles.popupmenu_selectvar2, 'value');      
        guidata(handles.output, handles);
    else        
        % se calculan las nuevas variables no básicas si las hay
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
            spreadsheet = cell(100,100);
            spreadsheet(1:Dimension(1), 1:Dimension(2)) = num2cell(All_display);
            set(handles.table_simplexdisplay, 'data', spreadsheet);
            %set(handles.table_simplexdisplay, 'data', All_display); 
            calc_nextassignment(handles);
        elseif strcmp(get(handles.Next_value, 'Label'), 'Siguiente nodo')
            Num_assignation = 0;
            T_VarType_Aux = T_VarType;
            %calc_nextciclevar(handles, [0, 0]);
            calc_nextciclevar(handles, 0);
        end
    end
else
    plano = get(handles.popupmenu_selectvar2, 'value');    
    if strcmp(get(handles_surf(plano), 'visible'), 'on')
        set(handles.pushbutton_asignall, 'string','Quitar');
    else
        set(handles.pushbutton_asignall, 'string','Poner');
    end
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


% --------------------------------------------------------------------
function Mode_geo3D_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Mode_geo3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global handles_surf handles_norm x_minmax y_minmax z_minmax;

if strcmp(get(hObject, 'Checked'), 'on')
    set(handles.table_simplexdisplay, 'Visible', 'on');
    set(handles.axes_simplex3D, 'Visible', 'off');
    set(handles.text_selectvar2, 'string', 'Inicial');
    set(hObject, 'Checked', 'off');
    set(handles.uipushtool1, 'Enable', 'off');
    set(handles.uitoggletool6, 'enable', 'off');
    set(handles.popupmenu_selectvar2, 'Enable', 'off');
    set(handles.popupmenu_selectvar2, 'string', char(' ', ' '));
    set(handles.popupmenu_selectvar2, 'value', 1);
    set(handles.pushbutton_asignall, 'String', 'Asignar todas');
    set(handles.pushbutton_asignall, 'Enable', 'off');
    set(handles.Restriction_nonnegativity, 'Checked', 'off');
    set(handles.Restriction_nonnegativity, 'Enable', 'off');
    set(handles.Uncut_planes, 'Enable', 'off');
    set(handles.Saveimage, 'Enable', 'off');
    hold off;
    handles_surf = [];
    handles_norm = [];
else
    handles_surf = [];
    handles_norm = [];
    cla(handles.axes_simplex3D);  
    x_minmax = zeros(1,2);
    y_minmax = zeros(1,2);
    z_minmax = zeros(1,2);
    
    set(handles.Restriction_nonnegativity, 'Enable', 'on');
    trace3D(handles);
    set(hObject, 'Checked', 'on');
    set(handles.uipushtool1, 'Enable', 'on');
    set(handles.uitoggletool6, 'enable', 'on');
    set(handles.axes_simplex3D, 'Visible', 'on');
    set(handles.Uncut_planes, 'Enable', 'off');    
    set(handles.table_simplexdisplay, 'Visible', 'off');
    set(handles.popupmenu_selectvar2, 'Enable', 'on');    
    set(handles.text_selectvar2, 'string', 'Plano');    
    set(handles.pushbutton_asignall, 'String', 'Quitar');
    set(handles.pushbutton_asignall, 'Enable', 'on');   
    set(handles.Saveimage, 'Enable', 'on');
end

function trace3D(handles)
global Matrix_problem Dimension Tableau handles_surf handles_norm x_minmax y_minmax z_minmax Order_current points_set table h;
%Construimos el sistema de coordenadas en la ventana
view(3); grid on;
xlabel(handles.axes_simplex3D, 'Eje X1');
ylabel(handles.axes_simplex3D, 'Eje X2');
zlabel(handles.axes_simplex3D, 'Eje X3');

% Obtenemos los puntos de intersección con los ejes y graficamos
% los planos

table = zeros(Dimension(1)+3, 4);
table(1:Dimension(1), 1:4) = Matrix_problem(1:Dimension(1),1:4); %%1
table(1:Dimension(1),end) = Matrix_problem(1:Dimension(1),end); %%1
table(Dimension(1)+1,:) = [0 0 1 0];
table(Dimension(1)+2,:) = [0 1 0 0];
table(Dimension(1)+3,:) = [1 0 0 0];

Basic_Order_current_3D = find(Order_current(1:Dimension(1)-1) < Dimension(1));
solution = zeros(1,Dimension(2)-1);
if ~isempty(Basic_Order_current_3D)    
    basic_solution = Tableau(1:Dimension(1)-1, end);
    sorted_basic_solution = basic_solution(Order_current(Order_current < Dimension(1)));    
    solution(Order_current(1:Dimension(2)-1) < Dimension(1)) = sorted_basic_solution;
    costos = Matrix_problem(end,1:Dimension(2)-1);
    Objective_function = costos(Basic_Order_current_3D)*solution(Basic_Order_current_3D)';
    table(Dimension(1),4) = Objective_function;    
else
    table(Dimension(1),4) = 0;
end
%table(end,end) = -Tableau(end,end);

if isempty(handles_surf)    
    handles_surf = zeros(1,Dimension(1)+3);
    handles_norm = zeros(1,Dimension(1));
    points_set = cell(1, Dimension(1)+3);
else
    delete(h);
end

% se construye el arreglo de índices de los planos    
Var = cell(Dimension(1), 1);
for i = 1:(Dimension(1))
    Var(i) = cellstr(['f',num2str(i), ': (', strtrim(rats(table(i,1))), ')X1 + (',strtrim(rats(table(i,2))), ')X2 + (', strtrim(rats(table(i,3))), ')X3 = ', strtrim(rats(table(i,4)))]);
end
set(handles.popupmenu_selectvar2, 'string', char(Var));
set(handles.popupmenu_selectvar2, 'value', 1);


for i=1:Dimension(1)
    if strcmp(get(handles.Mode_geo3D, 'Checked'),'off') || i == Dimension(1)
        isthereintercept = zeros(1,3); 
        z11 = 0;
        if table(i, 3) ~= 0
            isthereintercept(3) = 1;
            x11 = 0; y11 =0;  syms z11; %#ok<NASGU>
            z = solve([num2str(table(i, 1)), '*x11+', ...
                num2str(table(i, 2)), '*y11+', num2str(table(i, 3)), '*z11=', ...
                num2str(table(i, 4))], z11); z11=eval(z);
        end
        y21 = 0;
        if table(i, 2) ~= 0
            isthereintercept(2) = 1;
            x21 = 0; z21 =0;  syms y21; %#ok<NASGU>
            y = solve([num2str(table(i, 1)),'*x21+', ...
                num2str(table(i, 2)), '*y21+', num2str(table(i, 3)), '*z21=', ...
                num2str(table(i, 4))], y21); y21=eval(y);
        end
        x31 = 0;
        if table(i, 1) ~= 0
            isthereintercept(1) = 1;
            y31 = 0; z31 =0; syms x31;
            x = solve([num2str(table(i, 1)),'*x31+', ...
                num2str(table(i, 2)), '*y31+', num2str(table(i, 3)), '*z31=', ...
                num2str(table(i, 4))], x31); x31=eval(x);
        end              
        
        if all(isthereintercept > 0)
            if all([x31 y21 z11] > 0) 
                x1 = 0; y1 = 0; z1 = z11;
                x2 = 0; y2 = y21; z2 = 0;                
                x3 = x31; y3 = 0; z3 = 0;
                if i == Dimension(1)
                    outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                    x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);                    
                else
                    x4 = x31; y4 = 0; z4 = 0;               
                end
                points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 1);
            elseif any([x31 y21 z11] > 0)
                if all([x31 y21] > 0)
                    inner_point1 = 4*[x31 0 0]-3*[0 0 z11];
                    inner_point2 = 4*[0 y21 0]-3*[0 0 z11];
                    x1 = inner_point1(1); y1 = inner_point1(2); z1 = inner_point1(3);
                    x2 = inner_point2(1); y2 = inner_point2(2); z2 = inner_point2(3);
                    x3 = x31; y3 = 0; z3 = 0;
                    x4 = 0; y4 = y21; z4 = 0;
                    points_set{i} = struct('points', [x2 y2 z2; x4 y4 z4; x3 y3 z3; x1 y1 z1], 'tipo', 2);
                elseif all([x31 z11] >0)
                    inner_point1 = 4*[x31 0 0]-3*[0 y21 0];
                    inner_point2 = 4*[0 0 z11]-3*[0 y21 0];
                    x1 = inner_point1(1); y1 = inner_point1(2); z1 = inner_point1(3);
                    x2 = inner_point2(1); y2 = inner_point2(2); z2 = inner_point2(3);
                    x3 = 0; y3 = 0; z3 = z11;
                    x4 = x31; y4 = 0; z4 = 0;
                    points_set{i} = struct('points', [x3 y3 z3; x2 y2 z2; x4 y4 z4; x1 y1 z1], 'tipo', 3);
                elseif all([y21 z11] > 0)
                    inner_point1 = 4*[0 y21 0]-3*[x31 0 0];
                    inner_point2 = 4*[0 0 z11]-3*[x31 0 0];
                    x1 = inner_point1(1); y1 = inner_point1(2); z1 = inner_point1(3);
                    x2 = inner_point2(1); y2 = inner_point2(2); z2 = inner_point2(3);
                    x3 = 0; y3 = 0; z3 = z11;
                    x4 = 0; y4 = y21; z4 = 0;
                    points_set{i} = struct('points', [x3 y3 z3; x4 y4 z4; x1 y1 z1; x2 y2 z2], 'tipo', 4);
                elseif x31 > 0
                    inner_point1 = 4*[x31 0 0]-3*[0 0 z11];
                    inner_point2 = 4*[x31 0 0]-3*[0 y21 0];
                    x1 = inner_point1(1); y1 = inner_point1(2); z1 = inner_point1(3);
                    x2 = inner_point2(1); y2 = inner_point2(2); z2 = inner_point2(3);
                    x3 = x31; y3 = 0; z3 = 0;
                    if i == Dimension(1)
                        outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                        x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);                    
                    else
                        x4 = x31; y4 = 0; z4 = 0;               
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 5);
                elseif y21 > 0
                    inner_point1 = 4*[0 y21 0]-3*[0 0 z11];
                    inner_point2 = 4*[0 y21 0]-3*[x31 0 0];
                    x1 = inner_point1(1); y1 = inner_point1(2); z1 = inner_point1(3);
                    x2 = inner_point2(1); y2 = inner_point2(2); z2 = inner_point2(3);
                    x3 = 0; y3 = y21; z3 = 0;
                    if i == Dimension(1)
                        outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                        x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);                    
                    else
                        x4 = 0; y4 = y21; z4 = 0;               
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 6);
                elseif z11 > 0
                    inner_point1 = 4*[0 0 z11]-3*[0 y21 0];
                    inner_point2 = 4*[0 0 z11]-3*[x31 0 0];
                    x1 = inner_point1(1); y1 = inner_point1(2); z1 = inner_point1(3);
                    x2 = inner_point2(1); y2 = inner_point2(2); z2 = inner_point2(3);
                    x3 = 0; y3 = 0; z3 = z11;
                    if i == Dimension(1)
                        outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                        x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);                    
                    else
                        x4 = 0; y4 = 0; z4 = z11;               
                    end
                    points_set{i} = struct('points', [x3 y3 z3; x2 y2 z2; x1 y1 z1; x4 y4 z4], 'tipo', 7);
                end
                %%%desarrollar
            elseif all([x31 y21 z11] == 0)
                x_minmax(2) = max([x_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
                y_minmax(2) = max([y_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
                z_minmax(2) = max([z_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
                isthereinterceptline = zeros(1,3);
                
                x11 = x_minmax(2); z11 =0; syms y11;
                y = solve([num2str(table(i, 1)), '*x11+', ...
                num2str(table(i, 2)), '*y11+', num2str(table(i, 3)), '*z11=', ...
                num2str(table(i, 4))], y11); y11=eval(y);
                if y11 > 0
                    isthereinterceptline(1) = 1; % Plano XY
                end
                
                x21 = x_minmax(2); y21 = 0; syms z21;
                z = solve([num2str(table(i, 1)), '*x21+', ...
                num2str(table(i, 2)), '*y21+', num2str(table(i, 3)), '*z21=', ...
                num2str(table(i, 4))], z21); z21=eval(z);
                if z21 > 0
                    isthereinterceptline(2) = 1; % Plano XZ
                end
                
                x31 = 0; z31 =z_minmax(2); syms y31;
                y = solve([num2str(table(i, 1)), '*x31+', ...
                num2str(table(i, 2)), '*y31+', num2str(table(i, 3)), '*z31=', ...
                num2str(table(i, 4))], y31); y31=eval(y);
                if y31 > 0
                    isthereinterceptline(3) = 1; % Plano YZ
                end
                
                if isthereinterceptline(1) == 1 && isthereinterceptline(2) == 1
                    x1 = x11; y1 = y11; z1 = z11;
                    x2 = x21; y2 = y21; z2 = z21;
                    x3 = 0; y3 = 0; z3 = 0;
                    if i == Dimension(1)
                        outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                        x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);
                    else
                        x4 = 0; y4 = 0; z4 = 0;
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 8);
                elseif isthereinterceptline(1) == 1 && isthereinterceptline(3) == 1
                    x1 = x11; y1 = y11; z1 = z11;
                    x2 = x31; y2 = y31; z2 = z31;
                    x3 = 0; y3 = 0; z3 = 0;
                    if i == Dimension(1)
                        outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                        x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);
                    else
                        x4 = 0; y4 = 0; z4 = 0;
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 9);
                elseif isthereinterceptline(2) == 1 && isthereinterceptline(3) == 1
                    x1 = x21; y1 = y21; z1 = z21;
                    x2 = x31; y2 = y31; z2 = z31;
                    x3 = 0; y3 = 0; z3 = 0;
                    if i == Dimension(1)
                        outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                        x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);
                    else
                        x4 = 0; y4 = 0; z4 = 0;
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 10);
                else
                    if i == Dimension(1)
                        x1 = x11; y1 = y11; z1 = z11;
                        x2 = x21; y2 = y21; z2 = z21;
                        x3 = x31; y3 = y31; z3 = z31;                        
                        outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                        x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);                    
                    else
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;
                        x3 = 0; y3 = 0; z3 = 0;   
                        x4 = 0; y4 = 0; z4 = 0;
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 11);
                end
                %%% desarrollar
            else % si todos son negativos
                if i == Dimension(1)
                    x1 = x31; y1 = 0; z1 = 0;
                    x2 = 0; y2 = y21; z2 = 0;
                    x3 = 0; y3 = 0; z3 = z11;
                    outer_point = -1*[x1 y1 z1]+1*[x2 y2 z2]+1*[x3 y3 z3];
                    x4 = outer_point(1); y4 = outer_point(2); z4 = outer_point(3);                    
                else
                    x1 = 0; y1 = 0; z1 = 0;
                    x2 = 0; y2 = 0; z2 = 0;
                    x3 = 0; y3 = 0; z3 = 0;   
                    x4 = 0; y4 = 0; z4 = 0;
                end
                points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 12);
            end            
        else
            x_minmax(2) = max([x_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
            y_minmax(2) = max([y_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
            z_minmax(2) = max([z_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
            if isthereintercept(1) == 1 && isthereintercept(2) == 1
                if all([x31 y21] > 0) %% plano paralelo a Z
                    x1 = 0; y1 = y21; z1 = 0;
                    x2 = x31; y2 = 0; z2 = 0;
                    x3 = 0; y3 = y21; z3 = z_minmax(2);
                    x4 = x31; y4 = 0; z4 = z_minmax(2);
                    points_set{i} = struct('points', [x3 y3 z3; x1 y1 z1; x2 y2 z2; x4 y4 z4], 'tipo', 13);
                elseif any([x31 y21] > 0) 
                    if x31 >0 %% plano paralelo a Z con y31 < 0
                        inner_point1 = 4*[x31 0 0]-3*[0 y21 0];                        
                        
                        x1 = x31; y1 = 0; z1 = z_minmax(2);
                        x2 = inner_point1(1); y2 = inner_point1(2); z2 = z_minmax(2);
                        x3 = x31; y3 = 0; z3 = 0;
                        x4 = inner_point1(1); y4 = inner_point1(2); z4 = inner_point1(3);
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 14);
                    else %% plano paralelo a Z con x31 < 0
                        inner_point1 = 2*[0 y21 0]-1*[x31 0 0];                        
                       
                        x1 = 0; y1 = y21; z1 = z_minmax(2);
                        x2 = inner_point1(1); y2 = inner_point1(2); z2 = z_minmax(2);
                        x3 = 0; y3 = y21; z3 = 0;
                        x4 = inner_point1(1); y4 = inner_point1(2); z4 = inner_point1(3);
                        points_set{i} = struct('points', [x1 y1 z1; x3 y3 z3; x4 y4 z4; x2 y2 z2], 'tipo', 15);
                    end
                elseif all([x31 y21] == 0) %%%% Intersecta con eje Z
                    x11 = x_minmax(2); z11 = 0; syms y11; %#ok<NASGU>
                    y = solve([num2str(table(i, 1)), '*x11+', ...
                    num2str(table(i, 2)), '*y11+', num2str(table(i, 3)), '*z11=', ...
                    num2str(table(i, 4))], y11); y11=eval(y);
                
                    if y11 > 0
                        x1 = x11; y1 = y11; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;                            
                        x3 = x11; y3 = y11; z3 = z_minmax(2);
                        x4 = 0; y4 = 0; z4 = z_minmax(2);                       
                    else
                        if i == Dimension(1)                                                    
                            x1 = x11; y1 = y11; z1 = 0;
                            x2 = 0; y2 = 0; z2 = 0;                            
                            x3 = x11; y3 = y11; z3 = z_minmax(2);
                            x4 = 0; y4 = 0; z4 = z_minmax(2);
                        else
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = 0; y2 = 0; z2 = 0;
                            x3 = 0; y3 = 0; z3 = 0;   
                            x4 = 0; y4 = 0; z4 = 0;
                        end
                    end
                    points_set{i} = struct('points', [x4 y4 z4; x3 y3 z3; x2 y2 z2; x1 y1 z1], 'tipo', 16);
                else %%% si son negativos
                    if i == Dimension(1)                        
                        x1 = x31; y1 = 0; z1 = 0;
                        x2 = 0; y2 = y21; z2 = 0;
                        x3 = x31; y3 = 0; z3 = z_minmax(2);   
                        x4 = 0; y4 = y21; z4 = z_minmax(2);
                    else
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;
                        x3 = 0; y3 = 0; z3 = 0;   
                        x4 = 0; y4 = 0; z4 = 0;
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 17);
                end
            elseif isthereintercept(1) == 1 && isthereintercept(3) == 1
                if all([x31 z11] > 0) %% plano paralelo a Y
                    x1 = 0; y1 = 0; z1 = z11;
                    x2 = x31; y2 = 0; z2 = 0;
                    x3 = 0; y3 = y_minmax(2); z3 = z11;
                    x4 = x31; y4 = y_minmax(2); z4 = 0;
                    points_set{i} = struct('points', [x1 y1 z1; x3 y3 z3; x2 y2 z2; x4 y4 z4], 'tipo', 18);
                elseif any([x31 z11] > 0) 
                    if x31 >0 %% plano paralelo a Y con z31 < 0
                        inner_point1 = 4*[x31 0 0]-3*[0 0 z11];  
                        
                        x1 = x31; y1 = y_minmax(2); z1 = 0;
                        x2 = inner_point1(1); y2 = y_minmax(2); z2 = inner_point1(3);
                        x3 = x31; y3 = 0; z3 = 0;
                        x4 = inner_point1(1); y4 = inner_point1(2); z4 = inner_point1(3); 
                        points_set{i} = struct('points', [x4 y4 z4; x1 y1 z1; x3 y3 z3; x2 y2 z2], 'tipo', 19);
                    else %% plano paralelo a Y con x31 < 0
                        inner_point1 = 2*[0 0 z11]-1*[x31 0 0];                        
                        
                        x1 = 0; y1 = y_minmax(2); z1 = z11;
                        x2 = inner_point1(1); y2 = y_minmax(2); z2 = inner_point1(3);
                        x3 = 0; y3 = 0; z3 = z11;
                        x4 = inner_point1(1); y4 = inner_point1(2); z4 = inner_point1(3);
                        points_set{i} = struct('points', [x3 y3 z3; x1 y1 z1; x4 y4 z4; x4 y4 z4], 'tipo', 20);
                    end                    
                elseif all([x31 z11] == 0) %%%% Intersecta con eje Y
                    x21 = x_minmax(2); y21 = 0; syms z21; %#ok<NASGU>
                    z = solve([num2str(table(i, 1)), '*x21+', ...
                    num2str(table(i, 2)), '*y21+', num2str(table(i, 3)), '*z21=', ...
                    num2str(table(i, 4))], z21); z21=eval(z);
                    if z21 > 0                                                
                        x1 = 0; y1 = y_minmax(2); z1 = 0;
                        x2 = x21; y2 = y_minmax(2); z2 = z21;                        
                        x3 = 0; y3 = 0; z3 = 0;
                        x4 = x21; y4 = 0; z4 = z21;
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 21);
                    else
                        if i == Dimension(1)                        
                            x1 = 0; y1 = y_minmax(2); z1 = 0;
                            x2 = x21; y2 = y_minmax(2); z2 = z21;                        
                            x3 = 0; y3 = 0; z3 = 0;
                            x4 = x21; y4 = 0; z4 = z21;
                        else
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = 0; y2 = 0; z2 = 0;
                            x3 = 0; y3 = 0; z3 = 0;   
                            x4 = 0; y4 = 0; z4 = 0;
                        end
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 38);
                    end                    
                else %%%si son negativos
                    if i == Dimension(1)                        
                        x1 = x31; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = z11;
                        x3 = x31; y3 = y_minmax(2); z3 = 0;   
                        x4 = 0; y4 = y_minmax(2); z4 = z11;
                    else
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;
                        x3 = 0; y3 = 0; z3 = 0;   
                        x4 = 0; y4 = 0; z4 = 0;
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 22);
                end
            elseif isthereintercept(2) == 1 && isthereintercept(3) == 1
                if all([y21 z11] > 0) %% plano paralelo a X
                    x1 = 0; y1 = y21; z1 = 0;
                    x2 = x_minmax(2); y2 = y21; z2 = 0;
                    x3 = 0; y3 = 0; z3 = z11;
                    x4 = x_minmax(2); y4 = 0; z4 = z11;
                    %points_set{i} = struct('points', [x3 y3 z3; x1 y1 z1; x2 y2 z2; x4 y4 z4], 'tipo', 23);
                    points_set{i} = struct('points', [x3 y3 z3; x1 y1 z1; x4 y4 z4; x2 y2 z2], 'tipo', 23);
                elseif any([y21 z11] > 0)
                    if y21 >0 %% plano paralelo a X con z31 < 0
                        inner_point1 = 4*[0 y21 0]-3*[0 0 z11];                        
                        
                        x1 = x_minmax(2); y1 = y21; z1 = 0;
                        x2 = 0; y2 = y21; z2 = 0;
                        x3 = x_minmax(2); y3 = inner_point1(2); z3 = inner_point1(3);
                        x4 = inner_point1(1); y4 = inner_point1(2); z4 = inner_point1(3);
                        points_set{i} = struct('points', [x4 y4 z4; x2 y2 z2; x1 y1 z1; x3 y3 z3], 'tipo', 24);
                    else %% plano paralelo a Y con y31 < 0
                        inner_point1 = 2*[0 0 z11]-1*[0 y21 0];                        
                        
                        x1 = x_minmax(2); y1 = 0; z1 = z11;
                        x2 = 0; y2 = 0; z2 = z11;
                        x3 = x_minmax(2); y3 = inner_point1(2); z3 = inner_point1(3);
                        x4 = inner_point1(1); y4 = inner_point1(2); z4 = inner_point1(3);
                        points_set{i} = struct('points', [x2 y2 z2; x4 y4 z4; x1 y1 z1; x3 y3 z3], 'tipo', 25);
                    end                    
                elseif all([y21 z11] == 0) %%%% Con eje X
                    x31 = 0; y31 = y_minmax(2); syms z31; %#ok<NASGU>
                    z = solve([num2str(table(i, 1)), '*x31+', ...
                    num2str(table(i, 2)), '*y31+', num2str(table(i, 3)), '*z31=', ...
                    num2str(table(i, 4))], z31); z31=eval(z);
                    if z31 > 0
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = x_minmax(2); y2 = 0; z2 = 0;
                        x3 = 0; y3 = y31; z3 = z31;
                        x4 = x_minmax(2); y4 = y31; z4 = z31;
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 26);
                    else
                        if i == Dimension(1)                        
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = x_minmax(2); y2 = 0; z2 = 0;
                            x3 = 0; y3 = y31; z3 = z31;
                            x4 = x_minmax(2); y4 = y31; z4 = z31;
                        else
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = 0; y2 = 0; z2 = 0;
                            x3 = 0; y3 = 0; z3 = 0;   
                            x4 = 0; y4 = 0; z4 = 0;
                        end
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 27);
                    end                  
                else %%% si son negativos
                    if i == Dimension(1)                        
                        x1 = 0; y1 = y21; z1 = 0;
                        x2 = x_minmax(2); y2 = y21; z2 = 0;
                        x3 = 0; y3 = 0; z3 = z11;                           
                        x4 = x_minmax(2); y4 = 0; z4 = z11;
                    else
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;
                        x3 = 0; y3 = 0; z3 = 0;   
                        x4 = 0; y4 = 0; z4 = 0;
                    end
                    points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 28);
                end
            else
                x_minmax(2) = max([x_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
                y_minmax(2) = max([y_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
                z_minmax(2) = max([z_minmax(2), 0.50*max(table(1:Dimension(1),end)), 1]);
                if isthereintercept(1) == 1
                    if x31 > 0 %%% Plano paralelo a YZ
                        x1 = x31; y1 = 0; z1 = 0;   
                        x2 = x31; y2 = 0; z2 = z_minmax(2); 
                        x3 = x31; y3 = y_minmax(2); z3 = 0;
                        x4 = x31; y4 = y_minmax(2); z4 = z_minmax(2);
                        points_set{i} = struct('points', [x2 y2 z2; x4 y4 z4; x1 y1 z1; x3 y3 z3], 'tipo', 29); %%%% Referencia Plano XZ
                    elseif x31 == 0 % Plano YZ
                        if i == Dimension(1)                         
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = 0; y2 = y_minmax(2); z2 = 0;
                            x3 = 0; y3 = 0; z3 = z_minmax(2);   
                            x4 = 0; y4 = y_minmax(2); z4 = z_minmax(2);                       
                        else
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = 0; y2 = 0; z2 = 0;
                            x3 = 0; y3 = 0; z3 = 0;   
                            x4 = 0; y4 = 0; z4 = 0;
                        end
                        points_set{i} = struct('points', [x3 y3 z3; x2 y2 z2; x1 y1 z1; x4 y4 z4], 'tipo', 30);
                    else
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;
                        x3 = 0; y3 = 0; z3 = 0;   
                        x4 = 0; y4 = 0; z4 = 0;
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 31);
                    end
                elseif isthereintercept(2) == 1
                    if y21 > 0  %%% Plano paralelo a XZ                                          
                        x1 = 0; y1 = y21; z1 = z_minmax(2);                   
                        x2 = x_minmax(2); y2 = y21; z2 = z_minmax(2);                    
                        x3 = 0; y3 = y21; z3 = 0;
                        x4 = x_minmax(2); y4 = y21; z4 = 0;
                        points_set{i} = struct('points', [x1 y1 z1; x3 y3 z3; x4 y4 z4; x2 y2 z2], 'tipo', 32);
                    elseif y21 == 0 % Plano XZ
                        if i == Dimension(1) 
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = x_minmax(2); y2 = 0; z2 = 0;                        
                            x3 = 0; y3 = 0; z3 = z_minmax(2);   
                            x4 = x_minmax(2); y4 = 0; z4 = z_minmax(2);
                        else
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = 0; y2 = 0; z2 = 0;
                            x3 = 0; y3 = 0; z3 = 0;   
                            x4 = 0; y4 = 0; z4 = 0;
                        end
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 33);
                    else
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;
                        x3 = 0; y3 = 0; z3 = 0;   
                        x4 = 0; y4 = 0; z4 = 0;
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 34);
                    end
                elseif isthereintercept(3) == 1
                    if z11 > 0  %%% Plano paralelo a XY
                        x1 = 0; y1 = y_minmax(2); z1 = z11;
                        x2 = x_minmax(2); y2 = y_minmax(2); z2 = z11;
                        x3 = 0; y3 = 0; z3 = z11;
                        x4 = x_minmax(2); y4 = 0; z4 = z11;
                        points_set{i} = struct('points', [x3 y3 z3; x1 y1 z1; x4 y4 z4; x2 y2 z2], 'tipo', 35);
                    elseif z11 == 0 % Plano XY
                        if i == Dimension(1) 
                            x1 = 0; y1 = 0; z1 = 0; 
                            x2 = x_minmax(2); y2 = 0; z2 = 0;
                            x3 = 0; y3 = y_minmax(2); z3 = 0;                          
                            x4 = x_minmax(2); y4 = y_minmax(2); z4 = 0;                        
                        else
                            x1 = 0; y1 = 0; z1 = 0;
                            x2 = 0; y2 = 0; z2 = 0;
                            x3 = 0; y3 = 0; z3 = 0;   
                            x4 = 0; y4 = 0; z4 = 0;
                        end
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 36);
                    else
                        x1 = 0; y1 = 0; z1 = 0;
                        x2 = 0; y2 = 0; z2 = 0;
                        x3 = 0; y3 = 0; z3 = 0;   
                        x4 = 0; y4 = 0; z4 = 0;
                        points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 37);
                    end
                end
            end        
        end
        
        if i == Dimension(1)
            if handles_surf(i) ~= 0
                set(handles_surf(i), 'Visible', 'off');
                delete(handles_surf(i))
                set(handles_norm(i), 'Visible', 'off');
                delete(handles_norm(i))
            end
            if all(isthereintercept == 0)
                x1 = 0; y1 = 0; z1 = 0;
                x2 = 0; y2 = 0; z2 = 0;
                x3 = 0; y3 = 0; z3 = 0;   
                x4 = 0; y4 = 0; z4 = 0;
                points_set{i} = struct('points', [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4], 'tipo', 0);
            end
        end
        
        %p = [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4];
        p = points_set{i}.points;
        if ~(all(p(1,:) == p(2,:)) || all( p(1,:) == p(3,:)) || all(p(2,:) == p(3,:)))
            %points_set_nonred = points_set;
            K = convex_hull(p, points_set{i}.points);
            K_dim = size(K);
            if i == Dimension(1)
                c = 'y';
            else
                c = 'g';
            end
            if K_dim(1) == 5
                handles_surf(i) = patch('xdata',[p(K(1),1) p(K(2),1) p(K(3),1) p(K(4),1) p(K(5),1)], 'ydata',[p(K(1),2) p(K(2),2) p(K(3),2) p(K(4),2) p(K(5),2)], ...
                    'zdata', [p(K(1),3) p(K(2),3) p(K(3),3) p(K(4),3) p(K(5),3)], 'FaceColor', c); hold on; 
            else
                handles_surf(i) = patch('xdata',[p(K(1),1) p(K(2),1) p(K(3),1) p(K(4),1)], 'ydata',[p(K(1),2) p(K(2),2) p(K(3),2) p(K(4),2)], ...
                    'zdata', [p(K(1),3) p(K(2),3) p(K(3),3) p(K(4),3)], 'FaceColor', c); hold on; 
            end
        end        
        
        %mean_point = 0.3*[x1 y1 z1]+0.3*[x2 y2 z2]+0.4*[x3 y3 z3];
        mean_point = 0.3*p(1,:)+0.3*p(2,:)+0.4*p(3,:);
        handles_norm(i) = quiver3(handles.axes_simplex3D, mean_point(1),mean_point(2),mean_point(3), table(i, 1), table(i, 2), table(i, 3), 'Color', 'red');
        
        x_minmax(1) = min([x_minmax(1), x1, x2, x3, x4]);
        x_minmax(2) = max([x_minmax(2), x1, x2, x3, x4]);
        y_minmax(1) = min([y_minmax(1), y1, y2, y3, y4]);
        y_minmax(2) = max([y_minmax(2), y1, y2, y3, y4]);
        z_minmax(1) = min([z_minmax(1), z1, z2, z3, z4]);
        z_minmax(2) = max([z_minmax(2), z1, z2, z3, z4]);                
    end                          
end

points_set{Dimension(1)+1} = struct('points', [x_minmax(2) y_minmax(2) 0; 0 y_minmax(2) 0; x_minmax(2) 0 0; 0 0 0], 'tipo', 36); %XY
points_set{Dimension(1)+2} = struct('points', [0 0 z_minmax(2);  x_minmax(2) 0 z_minmax(2); x_minmax(2) 0 0; 0 0 0], 'tipo',33); %XZ
points_set{Dimension(1)+3} = struct('points', [0 0 z_minmax(2); 0 y_minmax(2) 0; 0 y_minmax(2) z_minmax(2); 0 0 0], 'tipo',30); %YZ
h = plot3(handles.axes_simplex3D, solution(1), solution(2), solution(3), 'rs');

if strcmp(get(handles.Mode_geo3D, 'Checked'), 'off')
    handles_surf(Dimension(1)+1) = patch('xdata',[0 x_minmax(2) x_minmax(2) 0],'ydata', [0 0 y_minmax(2) y_minmax(2)], ...
        'zdata', [0 0 0 0], 'FaceColor', 'b', 'visible', 'off', 'facealpha', .5); hold on; % XY
    handles_surf(Dimension(1)+2) = patch('xdata', [0 x_minmax(2) x_minmax(2) 0], 'ydata',[0 0 0 0], ... 
        'zdata', [0 0 z_minmax(2) z_minmax(2)], 'FaceColor', 'b', 'visible', 'off', 'facealpha', .5); hold on; % XZ
    handles_surf(Dimension(1)+3) = patch('xdata', [0 0 0 0],'ydata', [0 y_minmax(2) y_minmax(2) 0], ...
        'zdata', [0 0 z_minmax(2) z_minmax(2)], 'FaceColor', 'b', 'visible', 'off', 'facealpha', .5); hold on; %YZ
else
    if strcmp(get(handles.Restriction_nonnegativity, 'Checked'), 'off')   
        visibility = 'off';
    else
        visibility = 'on';
    end

    set(handles_surf(Dimension(1)+1), 'visible', visibility);
    set(handles_surf(Dimension(1)+2), 'visible', visibility);
    set(handles_surf(Dimension(1)+3), 'visible', visibility);
end

if strcmp(get(handles.uipushtool1, 'Enable'), 'on')
    set(handles.uipushtool1, 'Enable', 'on');
    set(handles.uitoggletool6, 'enable', 'on');
    set(handles.Uncut_planes, 'Enable', 'off');    
else
    set(handles.uipushtool1, 'Enable', 'off');
    set(handles.uitoggletool6, 'enable', 'off');
    set(handles.Uncut_planes, 'Enable', 'on');
end



% --------------------------------------------------------------------
function Control_panel_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Control_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.Control_panel, 'Checked'), 'off')
    set(handles.panel_enhancement, 'Visible', 'on')
    set(handles.Control_panel, 'Checked', 'on')
else
    set(handles.panel_enhancement, 'Visible', 'off')
    set(handles.Control_panel, 'Checked', 'off')
end


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles) %#ok<DEFNU,INUSL>
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global points_set Dimension table handles_surf handles_norm x_minmax y_minmax z_minmax minxyz points_set_nonred m_i m_j;
%e = 0.00005; %%%% Cambiar si se cambia precisión

points_set_nonred = points_set;
points_set_convex = points_set_nonred;
table_copy = table;

minxyz = [x_minmax(2) y_minmax(2) z_minmax(2)];
if strcmp(get(handles.uitoggletool6, 'state'), 'on')
    indexset1 = [randperm(Dimension(1)-1), Dimension(1)+1, Dimension(1)+2, Dimension(1)+3];
    indexset2 = randperm(Dimension(1)-1);
    indexset2_copy = indexset2;
else    
    indexset1 = [1:Dimension(1)-1, Dimension(1)+1, Dimension(1)+2, Dimension(1)+3];
    indexset2 = 1:Dimension(1)-1; 
    indexset2_copy = indexset2;
    %indexset1 = [1, 2, Dimension(1)+1, Dimension(1)+2, Dimension(1)+3];
    %indexset2 = [2, 1];  
    %indexset2_copy = indexset2;
end
    s = get(handles.listbox_operations, 'string');
    set(handles.listbox_operations, 'string', ...
        char(s, mat2str(indexset1), mat2str(indexset2)));
for i = indexset1
    if i ~= Dimension(1)
        all_points = [];
        if i > Dimension(1)
            indexset2 = indexset2_copy;
        else
            indexset2 = setdiff(indexset2,indexset1(1:find(indexset1==i)));
        end
        pasar = 1;
        for j = indexset2            
            if j ~= i && ~all(table_copy(i,:) == table_copy(j,:))
                syms x1 y1 z1 x2 y2 z2 x3 y3 z3;
                points_actual_i = points_set_convex{i}.points;
                points_actual_j = points_set_convex{j}.points;
                d_i = table_copy(i, 1:3); % vector director de i
                d_j = table_copy(j, 1:3); % vector director de j
                
                % Plano XY 
                
                if points_set_convex{i}.tipo == 1 || points_set_convex{i}.tipo == 2 || points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 13 || ...
                   points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                   points_set_convex{i}.tipo == 29 || points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 4 ||points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || ... 
                   points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 18 || points_set_convex{i}.tipo == 36 || points_set_convex{i}.tipo == 23 || points_set_convex{i}.tipo == 21 || ...       
                   points_set_convex{i}.tipo == 26
                    if points_set_convex{j}.tipo == 1 || points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 13 || ...
                       points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                       points_set_convex{j}.tipo == 29 || points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 4 ||points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || ... 
                       points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 18 || points_set_convex{j}.tipo == 36 || points_set_convex{j}.tipo == 23 || points_set_convex{j}.tipo == 21 || ...   
                       points_set_convex{j}.tipo == 26
                                           
                        N = cross([table_copy(i, 1), table_copy(i, 2), table_copy(i, 3)], [table_copy(j, 1), table_copy(j, 2), table_copy(j, 3)]);
                        if  any(N(:) ~= 0)
                            sol = solve([num2str(table_copy(i, 1)), '*x1+', num2str(table_copy(i, 2)), '*y1 =', num2str(table_copy(i, 4))], ...
                            [num2str(table_copy(j, 1)), '*x1+', num2str(table_copy(j, 2)), '*y1 =', num2str(table_copy(j, 4))], x1, y1);
                        else
                            sol = [];
                        end
                        
                        if (points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29 || ... %%%% Orden X-Y-Z-W
                                points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 18 || points_set_convex{i}.tipo == 16)
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([2 4 3 1],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([2 4 3 1],:);    
                        end
                        if (points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 29 || ...
                                points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 18 || points_set_convex{j}.tipo == 16)
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([2 4 3 1],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([2 4 3 1],:);    
                        end
                        if points_set_convex{i}.tipo == 26
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([4 2 1 3],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([4 2 1 3],:);    
                        end
                        if points_set_convex{j}.tipo == 26
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([4 2 1 3],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([4 2 1 3],:);    
                        end
                        if points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 21
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([2 1 3 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([2 1 3 4],:);    
                        end
                        if points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 21
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([2 1 3 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([2 1 3 4],:);    
                        end                        
                        
                        %%%% El vector director ya está considerado al eliminar casos redundantes cuando encontramos
                        %%%% interceptos y los casos acá desarrollados.
                        %%%% Agregar casos degenerados que pasan por el
                        %%%% origen (Caso 8, 9 y 16)                    
                        %m_i = (points_actual_i(2,2)-points_actual_i(3,2))/(points_actual_i(2,1)-points_actual_i(3,1));
                        %if isnan(m_i)
                            m_i = (points_set_nonred{i}.points(2,2)-points_set_nonred{i}.points(3,2))/(points_set_nonred{i}.points(2,1)-points_set_nonred{i}.points(3,1));
                        %end
                        %m_j = (points_actual_j(2,2)-points_actual_j(3,2))/(points_actual_j(2,1)-points_actual_j(3,1));
                        %if isnan(m_j)
                            m_j = (points_set_nonred{j}.points(2,2)-points_set_nonred{j}.points(3,2))/(points_set_nonred{j}.points(2,1)-points_set_nonred{j}.points(3,1));
                        %end
                        %if ((m_i < 0 || m_i == Inf)&& (m_j < 0 || m_j == Inf)) || ...
                        if (m_i < 0 && m_j < 0) || ...
                                (points_set_convex{i}.tipo == 36 || points_set_convex{i}.tipo == 33 || points_set_convex{i}.tipo == 30)
                                %((m_i == 0 || m_i == Inf) && (points_set_convex{i}.tipo == 36))
                            if ~isempty(sol)
                                sol = [sol.x1 sol.y1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2));
                                    if x > 0 && y > 0                                        
                                        %if points_set_nonred{i}.points(3, 1) <= points_set_nonred{j}.points(3, 1) %%% m_i <= m_j
                                        if m_i <= m_j
                                            if points_set_convex{i}.points(3, 3) == 0 && points_set_convex{j}.points(2, 3) == 0 %%%%%
                                                if points_actual_i(2, 1) <= x && points_actual_i(2, 2) >= y%%%11/01  Agregar tipo 
                                                    %%%%%%%
                                                    %points_actual_j_before = points_set_convex{j}.points;
                                                    %[points_set_convex{i}.points, points_set_convex{j}.points] = ZParalelPlanesInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                    %    x, y, table_copy, i, j, 2, 3, 4, 1);
                                                    %%%%%%%

                                                    points_set_convex{i}.points(2, :) = [x y 0];

                                                    %if points_set_convex{j}.points(3, 1) >= x && points_set_convex{j}.points(3, 2) <= y%%%11/01  Agregar tipo 
                                                    %    if all(points_actual_j_before(3, :) == points_actual_j_before(4, :))
                                                    %         points_set_convex{j}.points(4, :) = [x y 0]; %%%7/01 
                                                    %    end
                                                    %    points_set_convex{j}.points(3, :) = [x y 0]; 
                                                    %end
                                                %else
                                                end
                                                if points_actual_j(3, 1) >= x && points_actual_j(3, 2) <= y%%%11/01  Agregar tipo 
                                                    %%%%%%%
                                                    %points_actual_j_before = points_set_convex{j}.points;
                                                    %[points_set_convex{j}.points, points_set_convex{i}.points] = ZParalelPlanesInterceptXY(handles, points_set_convex{j}.points, points_set_convex{i}.points, ...
                                                    %    x, y, table_copy, i, j, 3, 2, 4, 1);
                                                    %%%%%%%

                                                    if all(points_actual_j(3, :) == points_actual_j(4, :))
                                                         points_set_convex{j}.points(4, :) = [x y 0]; %%%7/01 
                                                    end
                                                    points_set_convex{j}.points(3, :) = [x y 0];                                                                                                                                                                        
                                                end                                                
                                                %%%%%
                                                %%%Desarrollar planos tipo 2
                                                %%%%%
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 2 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 13) || ...
                                                    points_set_convex{i}.tipo == 13 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 13)
                                                    if points_actual_i(2, 1) <= x && points_actual_i(2, 2) >= y || ...
                                                        points_actual_j(3, 1) >= x && points_actual_j(3, 2) <= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(2, 1) <= x && points_actual_i(2, 2) >= y
                                                                        points_set_convex{i}.points(1, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(3, 1) >= x && points_actual_j(3, 2) <= y
                                                                        if points_set_convex{i}.tipo == 13 && points_set_convex{j}.tipo == 13 || ...
                                                                            points_set_convex{i}.tipo == 2 && points_set_convex{j}.tipo == 2
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            end
                                        else
                                            if points_set_convex{i}.points(2, 3) == 0 && points_set_convex{j}.points(3, 3) == 0 %%%%
                                                if points_actual_i(3, 2) <= y && points_actual_i(3, 1) >= x %%%%11/01
                                                    %%%%%%
                                                    %points_actual_i_before = points_set_convex{i}.points;
                                                    %[points_set_convex{i}.points, points_set_convex{j}.points] = ZParalelPlanesInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                    %    x, y, table_copy, i, j, 3, 2, 4, 1);
                                                    %%%%%%% 

                                                    if all(points_actual_i(3, :) == points_actual_i(4, :))
                                                         points_set_convex{i}.points(4, :) = [x y 0]; %%%7/01 
                                                    end
                                                    points_set_convex{i}.points(3, :) = [x y 0];

                                                    %if points_set_convex{j}.points(2, 2) >= y && points_set_convex{j}.points(2, 1) <= x %%%%11/01
                                                    %    points_set_convex{j}.points(2, :) = [x y 0];
                                                    %end
                                                %else
                                                end
                                                if points_actual_j(2, 2) >= y && points_actual_j(2, 1) <= x %%%%11/01
                                                    %%%%%%
                                                    %[points_set_convex{j}.points, points_set_convex{i}.points] = ZParalelPlanesInterceptXY(handles, points_set_convex{j}.points, points_set_convex{i}.points, ...
                                                    %    x, y, table_copy, i, j, 2, 3, 4, 1);
                                                    %%%%%%% 

                                                    points_set_convex{j}.points(2, :) = [x y 0];                                                        
                                                end
                                                %%%%%
                                                %%%Desarrollar planos tipo 2
                                                %%%%%
                                                if points_set_convex{i}.tipo == 2 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 13) || ...
                                                    points_set_convex{i}.tipo == 13 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 13)
                                                    if points_actual_i(3, 2) <= y && points_actual_i(3, 1) >= x || ...
                                                        points_actual_j(2, 2) >= y && points_actual_j(2, 1) <= x
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(2, 2) >= y && points_actual_j(2, 1) <= x
                                                                        points_set_convex{j}.points(1, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(3, 2) <= y && points_actual_i(3, 1) >= x
                                                                        if points_set_convex{i}.tipo == 13 && points_set_convex{j}.tipo == 13 || ...
                                                                            points_set_convex{i}.tipo == 2 && points_set_convex{j}.tipo == 2
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end                                                                    
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            end
                                        end                                    
                                    else
                                        %%%%%%
                                        caso = 1;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            if points_set_nonred{i}.points(3, 1) <= points_set_nonred{j}.points(3, 1)
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 3, 4, 1);
                                            else
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 2, 4, 1);
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXY(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%                                    
                                end
                            else
                                if points_set_nonred{i}.points(3, 1) <= points_set_nonred{j}.points(3, 1)
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 3, 4, 1);
                                else
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 2, 4, 1);
                                end
                            end                        
                        elseif ((m_i < 0 ) && ((m_j > 0 || m_j == Inf) && (points_set_convex{j}.tipo ==3 || ...
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 || ...
                                (d_j(1) >= 0 && d_j(2) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 16)))) || ...
                                ((m_j < 0 ) && ((m_i > 0 || m_i == Inf) && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(2) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 16))))
                            %%%%%%%
                            if ~isempty(sol)
                                sol = [sol.x1 sol.y1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y =eval(sol(2));
                                    if (x >= 0 && y > 0) || (x > 0 && y >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 0 1]); C = dot(N, [0 0 1]);

                                        if m_i < 0
                                            %if points_set_convex{i}.points(2, 3) == 0 && points_set_convex{j}.points(3, 3) == 0 %%%%%
                                                if points_actual_i(3, 1) >= x  && points_actual_i(3, 2) <= y
                                                    if points_set_convex{i}.points(3, 3) == 0                                                     
                                                        if all(points_set_convex{i}.points(3, :) == points_set_convex{i}.points(4, :))
                                                            points_set_convex{i}.points(4, :) = [x y 0];
                                                        end
                                                        points_set_convex{i}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x y points_set_convex{i}.points(4, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_j(2, 1) >= x  && points_actual_j(2, 2) >= y
                                                    if points_set_convex{j}.points(2, 3) == 0 
                                                        points_set_convex{j}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(1, :) = [x y points_set_convex{j}.points(1, 3)];
                                                        %points_set_convex{j}.points(1, :) = [x y 0];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 2 && (points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 8 || ...
                                                        points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 29) || ...
                                                        points_set_convex{i}.tipo == 13 && (points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 8 || ...
                                                        points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 29)
                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 2) <= y || ...
                                                        points_actual_j(2, 1) >= x  && points_actual_j(2, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 2) <= y
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(2, 1) >= x  && points_actual_j(2, 2) >= y
                                                                        if points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 29
                                                                            points_set_convex{j}.points(1, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end                                                                                                
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if points_set_convex{i}.points(3, 3) == 0 && points_set_convex{j}.points(2, 3) == 0 %%%%
                                                if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) <= y
                                                    if points_set_convex{j}.points(3, 3) == 0
                                                        if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                            points_set_convex{j}.points(4, :) = [x y 0];
                                                        end
                                                        points_set_convex{j}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x y points_set_convex{j}.points(4, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_i(2, 1) >= x  && points_actual_i(2, 2) >= y
                                                    if points_set_convex{i}.points(2, 3) == 0 
                                                        points_set_convex{i}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(1, :) = [x y points_set_convex{i}.points(1, 3)];
                                                        %points_set_convex{i}.points(1, :) = [x y 0];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 2 && (points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 8 || ...
                                                        points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 29) || ...
                                                        points_set_convex{j}.tipo == 13 && (points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 8 || ...
                                                        points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 29)
                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) <= y || ...
                                                       points_actual_i(2, 1) >= x  && points_actual_i(2, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) <= y
                                                                        points_set_convex{j}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(2, 1) >= x  && points_actual_i(2, 2) >= y
                                                                        if points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 29
                                                                            points_set_convex{i}.points(1, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,y]==0)
                                            %%%%%%
                                        %    caso = 2;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                        %                i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 2;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1  
                                            if m_i < 0
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXY(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < 0
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 2, 4, 1);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 3, 4, 1);
                                    %%%%%%%%%%%%%
                                end
                            end                            
                        elseif (m_i < 0 && (m_j >= 0 && (points_set_convex{j}.tipo ==4 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 23 ||...
                                points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(1) <= 0 && d_j(2) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 16)))) ...
                                || ...
                               (m_j < 0 && (m_i >= 0 && (points_set_convex{i}.tipo ==4 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 23 ||...
                                points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(1) <= 0 && d_i(2) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 16))))
                            %%%%%%%
                            if ~isempty(sol)
                                sol = [sol.x1 sol.y1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y =eval(sol(2));
                                    if (x >= 0 && y > 0) || (x > 0 && y >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 0 1]); C = dot(N, [0 0 1]);
                                        
                                        if m_i < 0 
                                            %if points_set_convex{i}.points(3, 3) == 0 && points_set_convex{j}.points(2, 3) == 0 %%%%
                                                if points_actual_i(2, 1) <= x  && points_actual_i(2, 2) >= y
                                                    if points_set_convex{i}.points(2, 3) == 0 
                                                        points_set_convex{i}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(1, :) = [x y points_set_convex{i}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_set_convex{j}.tipo == 16
                                                    points_actual_j = points_actual_j([4 3 2 1], :);
                                                    points_set_convex{j}.points = points_set_convex{j}.points([4 3 2 1], :);
                                                end
                                                
                                                if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) >= y
                                                    if points_set_convex{j}.points(3, 3) == 0
                                                        if all(points_actual_j(3, :) == points_actual_j(4, :))
                                                             points_set_convex{j}.points(4, :) = [x y 0]; %%%7/01 
                                                        end
                                                        points_set_convex{j}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x y points_set_convex{j}.points(4, 3)];
                                                        %points_set_convex{j}.points(4, :) = [x y 0];
                                                    %else                                                        
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end                                                                                        

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 2 && (points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 8 || ...
                                                          points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 15 ||  points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 32) || ...
                                                        points_set_convex{i}.tipo == 13 && (points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 8 || ...
                                                          points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 15 ||  points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 32)
                                                    if points_actual_i(2, 1) <= x  && points_actual_i(2, 2) >= y || ...
                                                        points_actual_j(3, 1) >= x  && points_actual_j(3, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(2, 1) <= x  && points_actual_i(2, 2) >= y
                                                                        points_set_convex{i}.points(1, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) >= y
                                                                        if points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 32
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end

                                                if points_set_convex{j}.tipo == 16
                                                    points_actual_j = points_actual_j([1 2 3 4], :); 
                                                    points_set_convex{j}.points = points_set_convex{j}.points([1 2 3 4], :);
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if points_set_convex{i}.points(2, 3) == 0 && points_set_convex{j}.points(3, 3) == 0 %%%%%
                                                if points_actual_j(2, 1) <= x  && points_actual_j(2, 2) >= y
                                                    if points_set_convex{j}.points(2, 3) == 0 
                                                        points_set_convex{j}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(1, :) = [x y points_set_convex{j}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_set_convex{i}.tipo == 16
                                                    points_actual_i = points_actual_i([4 3 2 1], :);  
                                                    points_set_convex{i}.points = points_set_convex{i}.points([4 3 2 1], :);
                                                end
                                                if points_actual_i(3, 1) >= x  && points_actual_i(3, 2) >= y
                                                    if points_set_convex{i}.points(3, 3) == 0 
                                                        if all(points_actual_i(3, :) == points_actual_i(4, :))
                                                             points_set_convex{i}.points(4, :) = [x y 0]; %%%7/01 
                                                        end
                                                        points_set_convex{i}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x y points_set_convex{i}.points(4, 3)];
                                                        %points_set_convex{i}.points(4, :) = [x y 0];
                                                    %else                                                        
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 2 && (points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 8 || ...
                                                          points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 15 ||  points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 32) || ...
                                                        points_set_convex{j}.tipo == 13 && (points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 8 || ...
                                                          points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 15 ||  points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 32)
                                                    if points_actual_j(2, 1) <= x  && points_actual_j(2, 2) >= y || ...
                                                       points_actual_i(3, 1) >= x  && points_actual_i(3, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(2, 1) <= x  && points_actual_j(2, 2) >= y
                                                                        points_set_convex{j}.points(1, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 2) >= y
                                                                        if points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 32
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end

                                                if points_set_convex{i}.tipo == 16
                                                    points_actual_i = points_actual_i([1 2 3 4], :);
                                                    points_set_convex{i}.points = points_set_convex{i}.points([1 2 3 4], :);
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,y]==0)
                                            %%%%%%
                                        %    caso = 3;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                        %                i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 3;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1 
                                            if m_i < 0 
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXY(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < 0 
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 3, 4, 1);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 2, 4, 1);
                                    %%%%%%%%%%%%%
                                end
                            end                            
                        elseif (((m_i >= 0 && (points_set_convex{i}.tipo ==4 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 23 ||...
                                points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(1) <= 0 && d_i(2) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || ...
                                points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) && ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==3 || ... %%%% m_j > 0 
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 || ...
                                (d_j(1) >= 0 && d_j(2) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || ...
                                points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))))) ...
                                ||...
                               (((m_j >= 0 && (points_set_convex{j}.tipo ==4 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 23 ||...
                                points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(1) <= 0 && d_j(2) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || ...
                                points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) && ...
                                (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(2) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || ...
                                points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))))
                            if ~isempty(sol)
                                sol = [sol.x1 sol.y1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2));
                                    if (x > 0 && y > 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 0 1]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 0 1]);

                                        if m_i < m_j
                                            %if points_set_convex{i}.points(2, 3) == 0 && points_set_convex{j}.points(3, 3) == 0 %%%%
                                                if points_actual_i(3, 1) >= x && points_actual_i(3, 2) >= y
                                                    if all(points_actual_i(3, :) == points_actual_i(4, :))
                                                        points_set_convex{i}.points(4, :) = [x y 0]; %%%7/01 
                                                    end
                                                    points_set_convex{i}.points(3, :) = [x y 0];
                                                    
                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x y points_set_convex{i}.points(4, 3)];                                            
                                                    %else                                                    
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                    %end
                                                end                                                

                                                if points_actual_j(2, 1) >= x && points_actual_j(2, 2) >= y
                                                    points_set_convex{j}.points(2, :) = [x y 0];
                                                    
                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(1, :) = [x y points_set_convex{j}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 5 && (points_set_convex{i}.tipo == 32 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 6) || ...
                                                    points_set_convex{j}.tipo == 14 && (points_set_convex{i}.tipo == 32 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 6) || ...
                                                    points_set_convex{j}.tipo == 29 && (points_set_convex{i}.tipo == 32 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 6)
                                                    if points_actual_j(2, 1) >= x && points_actual_j(2, 2) >= y || ...
                                                       points_actual_i(3, 1) >= x && points_actual_i(3, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(2, 1) >= x && points_actual_j(2, 2) >= y
                                                                        points_set_convex{j}.points(1, :) = eval(sol); %%%4
                                                                    end
                                                                    if points_actual_i(3, 1) >= x && points_actual_i(3, 2) >= y
                                                                        if points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if points_set_convex{i}.points(3, 3) == 0 && points_set_convex{j}.points(2, 3) == 0 %%%%%
                                                if points_actual_i(2, 1) >= x && points_actual_i(2, 2) >= y
                                                    points_set_convex{i}.points(2, :) = [x y 0];
                                                    
                                                %    if abs(C) < e
                                                %        points_set_convex{i}.points(1, :) = [x y points_set_convex{i}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                %    end
                                                end                                                

                                                if points_actual_j(3, 1) >= x && points_actual_j(3, 2) >= y
                                                    if all(points_actual_j(3, :) == points_actual_j(4, :))
                                                        points_set_convex{j}.points(4, :) = [x y 0]; %%%7/01 
                                                    end
                                                    points_set_convex{j}.points(3, :) = [x y 0];
                                                    
                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x y points_set_convex{j}.points(4, 3)];
                                                    %else                                                    
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                    %end
                                                end 
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 5 && (points_set_convex{j}.tipo == 32 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 6) || ...
                                                    points_set_convex{i}.tipo == 14 && (points_set_convex{j}.tipo == 32 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 6) || ...
                                                    points_set_convex{i}.tipo == 29 && (points_set_convex{j}.tipo == 32 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 6)
                                                    if points_actual_i(2, 1) >= x && points_actual_i(2, 2) >= y || ...
                                                       points_actual_j(3, 1) >= x && points_actual_j(3, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(2, 1) >= x && points_actual_i(2, 2) >= y
                                                                        points_set_convex{i}.points(1, :) = eval(sol); %%4
                                                                    end
                                                                    if points_actual_j(3, 1) >= x && points_actual_j(3, 2) >= y
                                                                        if points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 32
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                    else
                                        %%%%%%
                                        caso = 4;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1 
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXY(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 2, 4, 1);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 3, 4, 1);
                                    %%%%%%%%%%%%%
                                end
                            end                                                                        
                        elseif ((m_i >= 0 && (points_set_convex{i}.tipo ==4 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 23 ||...
                                points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(1) <= 0 && d_i(2) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) &&...
                                ((m_j >= 0 && (points_set_convex{j}.tipo ==4 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 23 ||...
                                points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(1) <= 0 && d_j(2) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))))...
                                ||...
                                ((m_j >= 0 && (points_set_convex{j}.tipo ==4 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 23 ||...
                                points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(1) <= 0 && d_j(2) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) && ...
                                ((m_i >= 0 && (points_set_convex{i}.tipo ==4 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 23 ||...
                                points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(1) <= 0 && d_i(2) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) %%%% Agregar tipo a ambos planos
                            if ~isempty(sol)
                                sol = [sol.x1 sol.y1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2));
                                    if (x >= 0 && y > 0) || (x > 0 && y >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 0 1]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 0 1]);
                                        
                                        if m_i < m_j
                                            %if points_set_convex{i}.points(3, 3) == 0 && points_set_convex{j}.points(2, 3) == 0 %%%%%
                                                if points_actual_i(2, 1) <= x  && points_actual_i(2, 2) <= y
                                                    if points_set_convex{i}.points(2, 3) == 0 
                                                        points_set_convex{i}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end
                                                    
                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(1, :) = [x y points_set_convex{i}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) >= y
                                                    if points_set_convex{j}.points(3, 3) == 0 
                                                        if all(points_actual_j(3, :) == points_actual_j(4, :))
                                                            points_set_convex{j}.points(4, :) = [x y 0]; %%%7/01 
                                                        end
                                                        points_set_convex{j}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x y points_set_convex{j}.points(4, 3)];
                                                    %else                                                        
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 6 && (points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32) || ...
                                                    points_set_convex{j}.tipo == 15 && (points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32) || ...
                                                    points_set_convex{j}.tipo == 32 && (points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32)
                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) >= y || ...
                                                       points_actual_i(2, 1) <= x  && points_actual_i(2, 2) <= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 2) >= y
                                                                        points_set_convex{j}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(2, 1) <= x  && points_actual_i(2, 2) <= y
                                                                        if points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32
                                                                            points_set_convex{i}.points(1, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if points_set_convex{i}.points(2, 3) == 0 && points_set_convex{j}.points(3, 3) == 0 %%%%
                                                if points_actual_j(2, 1) <= x  && points_actual_j(2, 2) <= y
                                                    if points_set_convex{j}.points(2, 3) == 0 
                                                        points_set_convex{j}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(1, :) = [x y points_set_convex{j}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_i(3, 1) >= x  && points_actual_i(3, 2) >= y
                                                    if points_set_convex{i}.points(3, 3) == 0 
                                                        if all(points_actual_i(3, :) == points_actual_i(4, :))
                                                            points_set_convex{i}.points(4, :) = [x y 0]; %%%7/01 
                                                        end
                                                        points_set_convex{i}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x y points_set_convex{i}.points(4, 3)];
                                                    %else                                                        
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 6 && (points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32) || ...
                                                    points_set_convex{j}.tipo == 15 && (points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32) || ...
                                                    points_set_convex{j}.tipo == 32 && (points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32)
                                                    if points_actual_j(2, 1) <= x  && points_actual_j(2, 2) <= y || ...
                                                       points_actual_i(3, 1) >= x  && points_actual_i(3, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(2, 1) <= x  && points_actual_j(2, 2) <= y
                                                                        points_set_convex{j}.points(1, :) = eval(sol); %%4
                                                                    end
                                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 2) >= y
                                                                        if points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 32
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,y]==0)
                                            %%%%%%
                                        %    caso = 5;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                        %                i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 5;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXY(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 3, 4, 1);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 2, 4, 1);
                                    %%%%%%%%%%%%%
                                end
                            end
                        elseif (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(2) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) && ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==3 || ...
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 || ...
                                (d_j(1) >= 0 && d_j(2) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))))...
                                || ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==3 || ...
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 || ...
                                (d_j(1) >= 0 && d_j(2) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) && ...
                                (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(2) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))))
                            if ~isempty(sol)
                                sol = [sol.x1 sol.y1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2));
                                    if (x > 0 && y > 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 0 1]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 0 1]);
                                        
                                        if m_i < m_j
                                            %if points_set_convex{i}.points(3, 3) == 0 && points_set_convex{j}.points(2, 3) == 0 %%%%%
                                                if points_actual_i(2, 1) >= x  && points_actual_i(2, 2) >= y
                                                    if points_set_convex{i}.points(2, 3) == 0 
                                                        points_set_convex{i}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(1, :) = [x y points_set_convex{i}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_j(3, 1) <= x  && points_actual_j(3, 2) <= y
                                                    if points_set_convex{j}.points(3, 3) == 0
                                                        if all(points_actual_j(3, :) == points_actual_j(4, :))
                                                            points_set_convex{j}.points(4, :) = [x y 0]; %%%7/01 
                                                        end
                                                        points_set_convex{j}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x y points_set_convex{j}.points(4, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 5 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29) || ...
                                                    points_set_convex{j}.tipo == 14 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29) || ...
                                                    points_set_convex{j}.tipo == 29 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29)
                                                    if points_actual_j(3, 1) <= x  && points_actual_j(3, 2) <= y || ...
                                                       points_actual_i(2, 1) >= x  && points_actual_i(2, 2) >= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(3, 1) <= x  && points_actual_j(3, 2) <= y
                                                                        points_set_convex{j}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(2, 1) >= x  && points_actual_i(2, 2) >= y
                                                                        if points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29
                                                                            points_set_convex{i}.points(1, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if points_set_convex{i}.points(2, 3) == 0 && points_set_convex{j}.points(3, 3) == 0 %%%%
                                                if points_actual_j(2, 1) >= x  && points_actual_j(2, 2) >= y
                                                    if points_set_convex{j}.points(2, 3) == 0 
                                                        points_set_convex{j}.points(2, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(1, :) = [x y points_set_convex{j}.points(1, 3)];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_i(3, 1) <= x  && points_actual_i(3, 2) <= y
                                                    if points_set_convex{i}.points(3, 3) == 0 
                                                        if all(points_actual_i(3, :) == points_actual_i(4, :))
                                                            points_set_convex{i}.points(4, :) = [x y 0]; %%%7/01 
                                                        end
                                                        points_set_convex{i}.points(3, :) = [x y 0];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x y 0];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x y points_set_convex{i}.points(4, 3)];
                                                    %else                                                        
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 5 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29) || ...
                                                    points_set_convex{j}.tipo == 14 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29) || ...
                                                    points_set_convex{j}.tipo == 29 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29)
                                                    if points_actual_j(2, 1) >= x  && points_actual_j(2, 2) >= y || ...
                                                       points_actual_i(3, 1) <= x  && points_actual_i(3, 2) <= y
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(3)-max([points_actual_i(:,3);points_actual_j(:,3)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(2, 1) >= x  && points_actual_j(2, 2) >= y
                                                                        points_set_convex{j}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(3, 1) <= x  && points_actual_i(3, 2) <= y
                                                                        if points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,y]==0)
                                            %%%%%%
                                        %    caso = 6;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                        %                i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 6;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_convex{i}, points_set_convex{j},  x, y, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 3, 4, 1);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 2, 4, 1);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXY(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 3, 4, 1);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXY(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 2, 4, 1);
                                    %%%%%%%%%%%%%
                                end
                            end
                        end
                        if (points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29 || ...
                                points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 18 || points_set_convex{i}.tipo == 16)
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([4 1 3 2],:); 
                            points_actual_i(1:4, :) = points_actual_i([4 1 3 2],:); 
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([4 1 3 2],:);    
                        end
                        if (points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 29 || ...
                                points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 18 || points_set_convex{j}.tipo == 16)
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([4 1 3 2],:);
                            points_actual_j(1:4, :) = points_actual_j([4 1 3 2],:);  
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([4 1 3 2],:);    
                        end  
                        if points_set_convex{i}.tipo == 26
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([3 2 4 1],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([3 2 4 1],:);    
                        end
                        if points_set_convex{j}.tipo == 26
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([3 2 4 1],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([3 2 4 1],:);    
                        end
                        if points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 21
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([2 1 3 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([2 1 3 4],:);    
                        end
                        if points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 21
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([2 1 3 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([2 1 3 4],:);    
                        end
                    end
                end

                %Plano XZ
                 if points_set_convex{i}.tipo == 1 || points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 14 || ...
                         points_set_convex{i}.tipo == 18 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 29 || points_set_convex{i}.tipo == 35 || ...
                         points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 2 || points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || ...
                         points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 13 || points_set_convex{i}.tipo == 33 || points_set_convex{i}.tipo == 23 || ... 
                         points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 16
                     if points_set_convex{j}.tipo == 1 || points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 14 || ...
                             points_set_convex{j}.tipo == 18 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 29 || points_set_convex{j}.tipo == 35 || ...
                             points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || ...
                             points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 13 || points_set_convex{j}.tipo == 33 || points_set_convex{j}.tipo == 23 || ... 
                             points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 16
                         
                        N = cross([table_copy(i, 1), table_copy(i, 2), table_copy(i, 3)], [table_copy(j, 1), table_copy(j, 2), table_copy(j, 3)]);
                        if  any(N(:) ~= 0)
                            sol = solve([num2str(table_copy(i, 1)), '*x2+', num2str(table_copy(i, 3)), '*z2 =', num2str(table_copy(i, 4))], ...
                            [num2str(table_copy(j, 1)), '*x2+', num2str(table_copy(j, 3)), '*z2 =', num2str(table_copy(j, 4))], x2, z2);
                        else
                            sol = [];
                        end
                        
                        if points_set_convex{i}.tipo == 4 %%%% Orden X-Z-Y-W
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([1 2 4 3],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([1 2 4 3],:);    
                        end
                        if points_set_convex{j}.tipo == 4
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([1 2 4 3],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([1 2 4 3],:);    
                        end
                        if points_set_convex{i}.tipo == 2 %|| points_set_convex{i}.tipo == 16 %%%% Orden X-Z-Y-W
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([4 2 3 1],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([4 2 3 1],:);    
                        end
                        if points_set_convex{j}.tipo == 2 %|| points_set_convex{i}.tipo == 16
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([4 2 3 1],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([4 2 3 1],:);    
                        end
                        if points_set_convex{i}.tipo == 10 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([3 2 1 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([3 2 1 4],:);    
                        end
                        if points_set_convex{j}.tipo == 10
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([3 2 1 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([3 2 1 4],:);    
                        end  
                        if points_set_convex{i}.tipo == 8 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([2 1 3 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([2 1 3 4],:);    
                        end
                        if points_set_convex{j}.tipo == 8
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([2 1 3 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([2 1 3 4],:);    
                        end
                        if points_set_convex{i}.tipo == 26 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([2 3 1 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([2 3 1 4],:);    
                        end
                        if points_set_convex{j}.tipo == 26
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([2 3 1 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([2 3 1 4],:);    
                        end
                        
                        %%%% El vector director ya está considerado al eliminar casos redundantes cuando encontramos
                        %%%% interceptos y los casos acá desarrollados.
                        %%%% Agregar casos degenerados que pasan por el
                        %%%% origen (Caso 8, 10 y 21) ya está                      
                        %m_i = (points_actual_i(1,3)-points_actual_i(3,3))/(points_actual_i(1,1)-points_actual_i(3,1));
                        %if isnan(m_i)
                            m_i = (points_set_nonred{i}.points(1,3)-points_set_nonred{i}.points(3,3))/(points_set_nonred{i}.points(1,1)-points_set_nonred{i}.points(3,1));
                        %end
                        %m_j = (points_actual_j(1,3)-points_actual_j(3,3))/(points_actual_j(1,1)-points_actual_j(3,1));
                        %if isnan(m_j)
                            m_j = (points_set_nonred{j}.points(1,3)-points_set_nonred{j}.points(3,3))/(points_set_nonred{j}.points(1,1)-points_set_nonred{j}.points(3,1));
                        %end
                        if (m_i < 0 && m_j < 0 ) || ...
                                (points_set_convex{i}.tipo == 36 || points_set_convex{i}.tipo == 33 || points_set_convex{i}.tipo == 30)
                                %((m_i == 0 || m_i == Inf) && (points_set_convex{i}.tipo == 33))
                            if ~isempty(sol)
                                sol = [sol.x2 sol.z2];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); z=eval(sol(2));
                                    if x > 0 && z > 0                                    
                                        %if points_set_nonred{i}.points(3, 1) <= points_set_nonred{j}.points(3, 1)
                                        if m_i <= m_j
                                            if (points_set_convex{i}.points(3, 2) == 0 || (points_set_convex{i}.points(3, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                               (points_set_convex{j}.points(1, 2) == 0 || (points_set_convex{j}.points(1, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_i(1, 1) <= x && points_actual_i(1, 3) >= z %%%%11/01
                                                    if points_set_convex{i}.points(1, 2) == 0 
                                                        %%%%%%
                                                        %points_actual_j_before = points_set_convex{j}.points;
                                                        %[points_set_convex{i}.points, points_set_convex{j}.points] = YParalelPlanesInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                        %    x, z, table_copy, i, j, 1, 3, 4, 2);
                                                        %%%%%%%
                                                        points_set_convex{i}.points(1, :) = [x 0 z];

                                                        %if points_set_convex{j}.points(3, 2) == 0
                                                            %if points_actual_j(3, 1) >= x && points_actual_j(3, 3) <= z %%%%11/01
                                                        %        if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                        %             points_set_convex{j}.points(4, :) = [x 0 z]; %%%7/01 
                                                        %        end
                                                        %        points_set_convex{j}.points(3, :) = [x 0 z];
                                                            %end
                                                        %else
                                                            %if points_actual_j(4, 1) >= x && points_actual_j(4, 3) <= z %%%%11/01
                                                        %        points_set_convex{j}.points(4, :) = [x 0 z];
                                                            %end
                                                        %end
                                                    else %if points_set_convex{i}.points(4, 2) == 0 
                                                        %if points_actual_i(4, 1) <= x && points_actual_i(4, 3) >= z %% Plano i
                                                            %%%%%%
                                                         %   points_actual_j_before = points_set_convex{j}.points;
                                                         %   [points_set_convex{i}.points, points_set_convex{j}.points] = YParalelPlanesInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                         %       x, z, table_copy, i, j, 4, 3, 1, 2);
                                                            %%%%%%%                                                                                                                
                                                            dim = size(points_set_convex{i}.points);
                                                            if dim(1) == 4
                                                                points_set_convex{i}.points(4, :) = [x 0 z];
                                                            else
                                                                points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                            end

                                                         %   if points_set_convex{j}.points(3, 2) == 0 
                                                                %if points_actual_j(3, 1) >= x && points_actual_j(3, 3) <= z %%%%11/01
                                                         %           if all(points_actual_j_before(3, :) == points_actual_j_before(4, :))
                                                         %                points_set_convex{j}.points(4, :) = [x 0 z]; %%%7/01 
                                                         %           end
                                                         %           points_set_convex{j}.points(3, :) = [x 0 z];
                                                                %end
                                                         %   else
                                                                %if points_actual_j(4, 1) >= x && points_actual_j(4, 3) <= z %%%%11/01                                                            
                                                         %           points_set_convex{j}.points(4, :) = [x 0 z];
                                                                %end
                                                         %   end
                                                        %else
                                                            %%%%%
                                                        %end
                                                    end
                                                end
                                                if points_actual_j(3, 1) >= x && points_actual_j(3, 3) <= z%%%11/01  Agregar tipo 
                                                    %points_actual_j_before = points_set_convex{j}.points;
                                                    %[points_set_convex{j}.points, points_set_convex{i}.points] = YParalelPlanesInterceptXZ(handles, points_set_convex{j}.points, points_set_convex{i}.points, ...
                                                    %    x, z, table_copy, i, j, 3, 1, 4, 2);
                                                    %%%%%%%

                                                    if points_set_convex{j}.points(3, 2) == 0
                                                        %if points_actual_j(3, 1) >= x && points_actual_j(3, 3) <= z %%%%11/01
                                                            if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                                 points_set_convex{j}.points(4, :) = [x 0 z]; %%%7/01 
                                                            end
                                                            points_set_convex{j}.points(3, :) = [x 0 z];
                                                        %end
                                                    else
                                                        %if points_actual_j(4, 1) >= x && points_actual_j(4, 3) <= z %%%%11/01
                                                            
                                                            dim = size(points_set_convex{j}.points);
                                                            if dim(1) == 4
                                                                points_set_convex{j}.points(4, :) = [x 0 z];
                                                            else
                                                                points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                            end
                                                        %end
                                                    end
                                                end
                                                %%%%%
                                                %%%Desarrollar planos tipo 3
                                                %%%%%
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 3 && (points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 18) || ...
                                                    points_set_convex{i}.tipo == 18 && (points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 18)
                                                    if points_actual_i(1, 1) <= x && points_actual_i(1, 3) >= z || ...
                                                        points_actual_j(3, 1) >= x && points_actual_j(3, 3) <= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 1) <= x && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(2, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(3, 1) >= x && points_actual_j(3, 3) <= z_aux
                                                                        if points_set_convex{i}.tipo == 18 && points_set_convex{j}.tipo == 18 || ...
                                                                            points_set_convex{i}.tipo == 3 && points_set_convex{j}.tipo == 3
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            end
                                        else
                                            if (points_set_convex{i}.points(1, 2) == 0 || (points_set_convex{i}.points(1, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                               (points_set_convex{j}.points(3, 2) == 0 || (points_set_convex{j}.points(3, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_j(1, 3) >= z && points_actual_j(1, 1) <= x %%%11/01
                                                    if points_set_convex{j}.points(1, 2) == 0
                                                        %points_actual_i_before = points_set_convex{i}.points;
                                                        %%%%%%
                                                        %[points_set_convex{i}.points, points_set_convex{j}.points] = YParalelPlanesInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                        %    x, z, table_copy, i, j, 3, 1, 4, 2);
                                                        %%%%%%% 
                                                                                                                
                                                        points_set_convex{j}.points(1, :) = [x 0 z];

                                                        %if points_set_convex{j}.points(1, 2) == 0
                                                            %if points_actual_j(1, 3) >= z && points_actual_j(1, 1) <= x %%%11/01
                                                        %        points_set_convex{j}.points(1, :) = [x 0 z];
                                                            %end
                                                        %else
                                                            %if points_actual_j(4, 3) >= z && points_actual_j(4, 1) <= x %%%11/01
                                                        %        points_set_convex{j}.points(4, :) = [x 0 z];
                                                            %end
                                                        %end
                                                    else %if points_set_convex{i}.points(3, 2) == 0
                                                        %if points_actual_i(4, 1) <= x && points_actual_i(4, 3) >= z %% Plano i
                                                            %%%%%%
                                                            %[points_set_convex{i}.points, points_set_convex{j}.points] = YParalelPlanesInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                            %    x, z, table_copy, i, j, 4, 1, 3, 2);
                                                            %%%%%%%                                                                                                                
                                                            dim = size(points_set_convex{j}.points);
                                                            if dim(1) == 4
                                                                points_set_convex{j}.points(4, :) = [x 0 z]; %%%7/01 
                                                            else
                                                                points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                            end
                                                            %if points_set_convex{j}.points(1, 2) == 0
                                                                %if points_actual_j(1, 1) >= x && points_actual_j(1, 3) <= z %%%%11/01
                                                            %        points_set_convex{j}.points(1, :) = [x 0 z];
                                                                %end
                                                            %else
                                                                %if points_actual_j(4, 1) >= x && points_actual_j(4, 3) <= z %%%%11/01
                                                            %        points_set_convex{j}.points(4, :) = [x 0 z];
                                                                %end
                                                            %end
                                                        %else
                                                            %%%%%%%
                                                        %end     
                                                    end
                                                end
                                                if points_actual_i(3, 1) >= x && points_actual_i(3, 3) <= z%%%11/01  Agregar tipo 
                                                    %%%%%%
                                                    %[points_set_convex{j}.points, points_set_convex{i}.points] = YParalelPlanesInterceptXZ(handles, points_set_convex{j}.points, points_set_convex{i}.points, ...
                                                    %    x, z, table_copy, i, j, 1, 3, 4, 2);
                                                    %%%%%%% 

                                                    if points_set_convex{i}.points(3, 2) == 0
                                                        if all(points_set_convex{i}.points(3, :) == points_set_convex{i}.points(4, :))
                                                            points_set_convex{i}.points(4, :) = [x 0 z];
                                                        end
                                                        %if points_actual_j(1, 3) >= z && points_actual_j(1, 1) <= x %%%11/01
                                                        points_set_convex{i}.points(3, :) = [x 0 z];
                                                        %end
                                                    else
                                                        %if points_actual_j(4, 3) >= z && points_actual_j(4, 1) <= x %%%11/01
                                                        dim = size(points_set_convex{i}.points);
                                                        if dim(1) == 4
                                                            points_set_convex{i}.points(4, :) = [x 0 z];
                                                        else
                                                            points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                        end                                                            
                                                        %end
                                                    end
                                                end
                                                %%%%%
                                                %%%Desarrollar planos tipo 3
                                                %%%%%
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 3 && (points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 18) || ...
                                                    points_set_convex{j}.tipo == 18 && (points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 18)
                                                    if points_actual_j(1, 1) <= x && points_actual_j(1, 3) >= z || ...
                                                        points_actual_i(3, 1) >= x && points_actual_i(3, 3) <= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(1, 1) <= x && points_actual_j(1, 3) >= z_aux
                                                                        points_set_convex{j}.points(2, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(3, 1) >= x && points_actual_i(3, 3) <= z_aux
                                                                        if points_set_convex{i}.tipo == 18 && points_set_convex{j}.tipo == 18 || ...
                                                                            points_set_convex{i}.tipo == 3 && points_set_convex{j}.tipo == 3
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            end
                                        end                                    
                                    else
                                        %%%%%%
                                        caso = 1;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            if points_set_nonred{i}.points(3, 1) <= points_set_nonred{j}.points(3, 1)
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 3, 4, 2);
                                            else
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 1, 4, 2);
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%                                    
                                end
                            else
                                if points_set_nonred{i}.points(3, 1) <= points_set_nonred{j}.points(3, 1)
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 3, 4, 2);
                                else
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 1, 4, 2);
                                end
                            end
                        elseif ((m_i < 0 ) && ((m_j > 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ... %%|| m_i == Inf
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 || ...
                                (d_j(1) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 21)))) ...
                                || ...
                                ((m_j < 0 ) && ((m_i > 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ... %%|| m_j == Inf
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 21))))
                            %%%%%%%
                            if ~isempty(sol)
                                sol = [sol.x2 sol.z2];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); z=eval(sol(2));
                                    if (x >= 0 && z > 0) || (x > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 1 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 1 0]);
                                        
                                        if m_i < 0
                                            %if (points_set_convex{i}.points(1, 2) == 0 || (points_set_convex{i}.points(1, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(3, 2) == 0 || (points_set_convex{j}.points(3, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_i(3, 1) >= x  && points_actual_i(3, 3) <= z
                                                    if points_set_convex{i}.points(3, 2) == 0
                                                        if all(points_set_convex{i}.points(3, :) == points_set_convex{i}.points(4, :))
                                                            points_set_convex{i}.points(4, :) = [x 0 z];
                                                        end
                                                        points_set_convex{i}.points(3, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x points_set_convex{i}.points(4, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_j(1, 1) >= x  && points_actual_j(1, 3) >= z
                                                    if points_set_convex{j}.points(1, 2) == 0 
                                                        points_set_convex{j}.points(1, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(2, :) = [x points_set_convex{j}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 3 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 8 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 29) || ...
                                                        points_set_convex{i}.tipo == 18 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 8 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 29)
                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 3) <= z || ...
                                                        points_actual_j(1, 1) >= x  && points_actual_j(1, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z + 1; 
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 3) <= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(1, 1) >= x  && points_actual_j(1, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 29
                                                                            points_set_convex{j}.points(2, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(3, 2) == 0 || (points_set_convex{i}.points(3, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(1, 2) == 0 || (points_set_convex{j}.points(1, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_j(3, 1) >= x  && points_actual_j(3, 3) <= z
                                                    if points_set_convex{j}.points(3, 2) == 0
                                                        if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                            points_set_convex{j}.points(4, :) = [x 0 z];
                                                        end
                                                        points_set_convex{j}.points(3, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x points_set_convex{j}.points(4, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_i(1, 1) >= x  && points_actual_i(1, 3) >= z
                                                    if points_set_convex{i}.points(1, 2) == 0 
                                                        points_set_convex{i}.points(1, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(2, :) = [x points_set_convex{i}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 3 && (points_set_convex{i}.tipo == 2 || points_set_convex{i}.tipo == 8 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 29) || ...
                                                        points_set_convex{j}.tipo == 18 && (points_set_convex{i}.tipo == 2 || points_set_convex{i}.tipo == 8 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 29)
                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 3) <= z || ...
                                                       points_actual_i(1, 1) >= x  && points_actual_i(1, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 3) <= z_aux
                                                                        points_set_convex{j}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(1, 1) >= x  && points_actual_i(1, 3) >= z_aux
                                                                        if points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 29
                                                                            points_set_convex{i}.points(2, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %   [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,z]==0)
                                            %%%%%%
                                        %    caso = 2;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 2;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1 
                                            if m_i < 0
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%% 
                                end
                            else
                                if m_i < 0
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 1, 4, 2);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 3, 4, 2);
                                    %%%%%%%%%%%%%
                                end
                            end                            
                        elseif (m_i < 0 && (m_j >= 0 && (points_set_convex{j}.tipo ==4 || ... %%m_i >= 0
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 23 || ...
                                points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(1) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 21))))...
                                || ...
                               (m_j < 0 && (m_i >= 0 && (points_set_convex{i}.tipo ==4 || ... %%m_j >= 0
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 23 || ...
                                points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(1) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 21))))
                            if ~isempty(sol)
                                sol = [sol.x2 sol.z2];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); z=eval(sol(2));
                                    if (x >= 0 && z > 0) || (x > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 1 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 1 0]);
                                        
                                        if m_i < 0
                                            %if (points_set_convex{i}.points(3, 2) == 0 || (points_set_convex{i}.points(3, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(1, 2) == 0 || (points_set_convex{j}.points(1, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_i(1, 1) <= x  && points_actual_i(1, 3) >= z
                                                    if points_set_convex{i}.points(1, 2) == 0 
                                                        points_set_convex{i}.points(1, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(2, :) = [x points_set_convex{i}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_j(3, 1) >= x  && points_actual_j(3, 3) >= z
                                                    if points_set_convex{j}.points(3, 2) == 0
                                                        if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                            points_set_convex{j}.points(4, :) = [x 0 z];
                                                        end
                                                        points_set_convex{j}.points(3, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x points_set_convex{j}.points(4, 2) z];
                                                    %else                                                        
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 3 && (points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 8 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 35) || ...
                                                        points_set_convex{i}.tipo == 18 && (points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 8 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 35)
                                                    if points_actual_i(1, 1) <= x  && points_actual_i(1, 3) >= z || ...
                                                        points_actual_j(3, 1) >= x  && points_actual_j(3, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 1) <= x  && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(2, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(1, 2) == 0 || (points_set_convex{i}.points(1, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(3, 2) == 0 || (points_set_convex{j}.points(3, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_j(1, 1) <= x  && points_actual_j(1, 3) >= z
                                                    if points_set_convex{j}.points(1, 2) == 0 
                                                        points_set_convex{j}.points(1, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(2, :) = [x points_set_convex{j}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_i(3, 1) >= x  && points_actual_i(3, 3) >= z
                                                    if points_set_convex{i}.points(3, 2) == 0 
                                                        if all(points_set_convex{i}.points(3, :) == points_set_convex{i}.points(4, :))
                                                            points_set_convex{i}.points(4, :) = [x 0 z];
                                                        end
                                                        points_set_convex{i}.points(3, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x points_set_convex{i}.points(4, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 3 && (points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 8 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 35) || ...
                                                        points_set_convex{j}.tipo == 18 && (points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 8 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 35)
                                                    if points_actual_j(1, 1) <= x  && points_actual_j(1, 3) >= z || ...
                                                       points_actual_i(3, 1) >= x  && points_actual_i(3, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(1, 1) <= x  && points_actual_j(1, 3) >= z_aux
                                                                        points_set_convex{j}.points(2, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 3) >= z_aux
                                                                        if points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 35
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,z]==0)
                                            %%%%%%
                                        %    caso = 3;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 3;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1 
                                            if m_i < 0
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%% 
                                end
                            else
                                if m_i < 0
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 3, 4, 2);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 1, 4, 2);
                                    %%%%%%%%%%%%%
                                end
                            end
                        elseif ((m_i >= 0 && (points_set_convex{i}.tipo ==4 || ...
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 23 || ...
                                points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(1) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) && ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ...
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 ||...
                                (d_j(1) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) ...
                                ||...
                                ((m_j >= 0 && (points_set_convex{j}.tipo ==4 || ...
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 23 || ...
                                points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(1) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) && ...
                                (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ...
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))))
                            if ~isempty(sol)
                                sol = [sol.x2 sol.z2];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); z=eval(sol(2));
                                    if (x > 0 && z > 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 1 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 1 0]);
                                        
                                        if m_i < m_j
                                            %if (points_set_convex{i}.points(1, 2) == 0 || (points_set_convex{i}.points(1, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(3, 2) == 0 || (points_set_convex{j}.points(3, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_i(3, 1) >= x && points_actual_i(3, 3) >= z
                                                    if all(points_set_convex{i}.points(3, :) == points_set_convex{i}.points(4, :))
                                                        points_set_convex{i}.points(4, :) = [x 0 z];
                                                    end
                                                    points_set_convex{i}.points(3, :) = [x 0 z];
                                                    
                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x points_set_convex{i}.points(4, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                    %end
                                                end                                                

                                                if points_actual_j(1, 1) >= x && points_actual_j(1, 3) >= z
                                                    points_set_convex{j}.points(1, :) = [x 0 z];
                                                    
                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(2, :) = [x points_set_convex{j}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 7 && (points_set_convex{i}.tipo == 35 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 5) || ...
                                                    points_set_convex{j}.tipo == 19 && (points_set_convex{i}.tipo == 35 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 5) || ...
                                                    points_set_convex{j}.tipo == 29 && (points_set_convex{i}.tipo == 35 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 5)
                                                    if points_actual_j(1, 1) >= x && points_actual_j(1, 3) >= z || ...
                                                       points_actual_i(3, 1) >= x && points_actual_i(3, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(1, 1) >= x && points_actual_j(1, 3) >= z_aux
                                                                        points_set_convex{j}.points(2, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(3, 1) >= x && points_actual_i(3, 3) >= z_aux
                                                                        if points_set_convex{i}.tipo == 35
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %   [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(3, 2) == 0 || (points_set_convex{i}.points(3, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(1, 2) == 0 || (points_set_convex{j}.points(1, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_i(1, 1) >= x && points_actual_i(1, 3) >= z
                                                    points_set_convex{i}.points(1, :) = [x 0 z];
                                                    
                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(2, :) = [x points_set_convex{i}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(s, handles.listbox_operations, 'string', ...
                                                    %        char('La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                    %end
                                                end                                                

                                                if points_actual_j(3, 1) >= x && points_actual_j(3, 3) >= z
                                                    if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                        points_set_convex{j}.points(4, :) = [x 0 z];
                                                    end
                                                    points_set_convex{j}.points(3, :) = [x 0 z];
                                                    
                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x points_set_convex{j}.points(4, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 7 && (points_set_convex{j}.tipo == 35 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 5) || ...
                                                    points_set_convex{i}.tipo == 19 && (points_set_convex{j}.tipo == 35 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 5) || ...
                                                    points_set_convex{i}.tipo == 29 && (points_set_convex{j}.tipo == 35 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 5)
                                                    if points_actual_i(1, 1) >= x && points_actual_i(1, 3) >= z || ...
                                                       points_actual_j(3, 1) >= x && points_actual_j(3, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 1) >= x && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(2, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(3, 1) >= x && points_actual_j(3, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                    else
                                        %%%%%%
                                        caso = 4;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                                    i, j,caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 1, 4, 2);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 3, 4, 2);
                                    %%%%%%%%%%%%%
                                end
                            end                                                                    
                        elseif ((m_i >= 0 && (points_set_convex{i}.tipo ==4 || ...
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 23 || ...
                                points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(1) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) && ...
                                ((m_j >= 0 && (points_set_convex{j}.tipo ==4 || ...
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 23 || ...
                                points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(1) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))))...
                                ||...
                               ((m_j >= 0 && (points_set_convex{j}.tipo ==4 || ...
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 23 || ...
                                points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(1) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) && ...
                                ((m_i >= 0 && (points_set_convex{i}.tipo ==4 || ...
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 23 || ...
                                points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(1) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))))
                            if ~isempty(sol)
                                sol = [sol.x2 sol.z2];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); z=eval(sol(2));
                                    if (x >= 0 && z > 0) || (x > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 1 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 1 0]);
                                        
                                        if m_i < m_j
                                            %if (points_set_convex{i}.points(3, 2) == 0 || (points_set_convex{i}.points(3, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(1, 2) == 0 || (points_set_convex{j}.points(1, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_i(1, 1) <= x  && points_actual_i(1, 3) <= z
                                                    if points_set_convex{i}.points(1, 2) == 0 
                                                        points_set_convex{i}.points(1, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(2, :) = [x points_set_convex{i}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_j(3, 1) >= x  && points_actual_j(3, 3) >= z
                                                    if points_set_convex{j}.points(3, 2) == 0
                                                        if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                            points_set_convex{j}.points(4, :) = [x 0 z];
                                                        end
                                                        points_set_convex{j}.points(3, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [x points_set_convex{j}.points(4, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 7 && (points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 35) || ...
                                                    points_set_convex{j}.tipo == 20 && (points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 35) || ...
                                                    points_set_convex{j}.tipo == 35 && (points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 35)
                                                    if points_actual_i(1, 1) <= x  && points_actual_i(1, 3) <= z || ...
                                                       points_actual_j(3, 1) >= x  && points_actual_j(3, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 1) <= x  && points_actual_i(1, 3) <= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(3, 1) >= x  && points_actual_j(3, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(1, 2) == 0 || (points_set_convex{i}.points(1, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.points(3, 2) == 0 || (points_set_convex{j}.points(3, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                                if points_actual_j(1, 1) <= x  && points_actual_j(1, 3) <= z
                                                    if points_set_convex{j}.points(1, 2) == 0 
                                                        points_set_convex{j}.points(1, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(2, :) = [x points_set_convex{j}.points(2, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end

                                                if points_actual_i(3, 1) >= x  && points_actual_i(3, 3) >= z
                                                    if points_set_convex{i}.points(3, 2) == 0
                                                        if all(points_set_convex{i}.points(3, :) == points_set_convex{i}.points(4, :))
                                                            points_set_convex{i}.points(4, :) = [x 0 z];
                                                        end                                                        
                                                        points_set_convex{i}.points(3, :) = [x 0 z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [x points_set_convex{i}.points(4, 2) z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.'));
                                                    %end
                                                end
                                                
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 7 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 35) || ...
                                                    points_set_convex{i}.tipo == 20 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 35) || ...
                                                    points_set_convex{i}.tipo == 35 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 35)
                                                    if points_actual_j(1, 1) <= x  && points_actual_j(1, 3) <= z || ...
                                                       points_actual_i(3, 1) >= x  && points_actual_i(3, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(1, 1) <= x  && points_actual_j(1, 3) <= z_aux
                                                                        points_set_convex{j}.points(2, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(3, 1) >= x  && points_actual_i(3, 3) >= z_aux
                                                                        if points_set_convex{i}.tipo == 35
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,z]==0)
                                            %%%%%%
                                        %    caso = 5;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 5;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1 
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 3, 4, 2);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 1, 4, 2);
                                    %%%%%%%%%%%%%
                                end
                            end
                        elseif (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ...
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) && ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ...
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 || ...
                                (d_j(1) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))))...
                                || ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ...
                                points_set_convex{j}.tipo == 5 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 14 || points_set_convex{j}.tipo == 19 || points_set_convex{j}.tipo == 29 || ...
                                (d_j(1) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 8 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) && ...
                                (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ...
                                points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 14 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29 || ...
                                (d_i(1) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 8 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))))
                            if ~isempty(sol)
                                sol = [sol.x2 sol.z2];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); z=eval(sol(2));
                                    if (x >= 0 && z > 0) || (x > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [0 1 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [0 1 0]);
                                        
                                        if m_i < m_j
                                            %if (points_set_convex{i}.points(3, 2) == 0 || (points_set_convex{i}.points(3, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.point
                                            %   s(1, 2) == 0 || (points_set_convex{j}.points(1, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                            if points_actual_i(1, 1) >= x  && points_actual_i(1, 3) >= z
                                                if points_set_convex{i}.points(1, 2) == 0 
                                                    points_set_convex{i}.points(1, :) = [x 0 z];
                                                else
                                                    dim = size(points_set_convex{i}.points);
                                                    points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                end

                                                %if abs(C) < e
                                                %    points_set_convex{i}.points(2, :) = [x points_set_convex{i}.points(2, 2) z];
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                %end
                                            end

                                            if points_actual_j(3, 1) <= x  && points_actual_j(3, 3) <= z
                                                if points_set_convex{j}.points(3, 2) == 0 
                                                    if all(points_set_convex{j}.points(3, :) == points_set_convex{j}.points(4, :))
                                                        points_set_convex{j}.points(4, :) = [x y 0];
                                                    end
                                                    points_set_convex{j}.points(3, :) = [x 0 z];
                                                else
                                                    dim = size(points_set_convex{j}.points);
                                                    points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                end

                                                %if abs(C2) < e
                                                %    points_set_convex{j}.points(4, :) = [x points_set_convex{j}.points(4, 2) z];
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                %end
                                            end
                                                
                                            %%%%%% Modificar (Revisar)
                                            if points_set_convex{j}.tipo == 5 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29) || ...
                                                points_set_convex{j}.tipo == 19 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29) || ...
                                                points_set_convex{j}.tipo == 29 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29)
                                                if points_actual_i(1, 1) >= x  && points_actual_i(1, 3) >= z || ...
                                                   points_actual_j(3, 1) <= x  && points_actual_j(3, 3) <= z
                                                    z_aux = z;
                                                    syms x4 y4 z4 z;
                                                    sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                    if ~isempty(sol)
                                                        sol = [sol.x4 sol.y4 sol.z4];
                                                        if ~isempty(symvar(sol))
                                                            %points_set_convex{j}.points(4, :) = [0 y z];
                                                            %z = z + 1;
                                                            z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                            z = z_sol; p_new = eval(eval(sol));
                                                            if all(p_new(:) >= 0)
                                                                if points_actual_i(1, 1) >= x  && points_actual_i(1, 3) >= z_aux
                                                                    points_set_convex{i}.points(2, :) = eval(sol); %%3
                                                                end
                                                                if points_actual_j(3, 1) <= x  && points_actual_j(3, 3) <= z_aux
                                                                    if points_set_convex{j}.tipo == 29
                                                                        points_set_convex{j}.points(4, :) = eval(sol); %%2
                                                                    else
                                                                        dim = size(points_set_convex{j}.points);
                                                                        points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(1, 2) == 0 || (points_set_convex{i}.points(1, 2) > 0 && points_set_convex{i}.points(4, 2) == 0 )) && ... %%%% OJO
                                            %   (points_set_convex{j}.point
                                            %   s(3, 2) == 0 || (points_set_convex{j}.points(3, 2) > 0 && points_set_convex{j}.points(4, 2) == 0 )) %%%% OJO
                                            if points_actual_j(1, 1) >= x  && points_actual_j(1, 3) >= z
                                                if points_set_convex{j}.points(1, 2) == 0 
                                                    points_set_convex{j}.points(1, :) = [x 0 z];
                                                else
                                                    dim = size(points_set_convex{j}.points);
                                                    points_set_convex{j}.points(dim(1)+1, :) = [x 0 z];
                                                end

                                                %if abs(C2) < e
                                                %    points_set_convex{j}.points(2, :) = [x points_set_convex{j}.points(2, 2) z];
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                %end
                                            end

                                            if points_actual_i(3, 1) <= x  && points_actual_i(3, 3) <= z
                                                if points_set_convex{i}.points(3, 2) == 0
                                                    if all(points_set_convex{i}.points(3, :) == points_set_convex{i}.points(4, :))
                                                        points_set_convex{i}.points(4, :) = [x y 0];
                                                    end
                                                    points_set_convex{i}.points(3, :) = [x 0 z];
                                                else
                                                    dim = size(points_set_convex{i}.points);
                                                    points_set_convex{i}.points(dim(1)+1, :) = [x 0 z];
                                                end

                                                %if abs(C) < e
                                                %    points_set_convex{i}.points(4, :) = [x points_set_convex{i}.points(4, 2) z];
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.'));
                                                %end
                                            end
                                                
                                            %%%%%% Modificar (Revisar)
                                            if points_set_convex{j}.tipo == 5 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29) || ...
                                                points_set_convex{j}.tipo == 19 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29) || ...
                                                points_set_convex{j}.tipo == 29 && (points_set_convex{i}.tipo == 5 || points_set_convex{i}.tipo == 19 || points_set_convex{i}.tipo == 29)                                                
                                                if points_actual_j(1, 1) >= x  && points_actual_j(1, 3) >= z || ...
                                                   points_actual_i(3, 1) <= x  && points_actual_i(3, 3) <= z
                                                    z_aux = z;
                                                    syms x4 y4 z4 z;
                                                    sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                    if ~isempty(sol)
                                                        sol = [sol.x4 sol.y4 sol.z4];
                                                        if ~isempty(symvar(sol))
                                                            %points_set_convex{j}.points(4, :) = [0 y z];
                                                            %z = z + 1;
                                                            z_sol = solve(sol(2)-max([points_actual_i(:,2);points_actual_j(:,2)]),z);
                                                            z = z_sol; p_new = eval(eval(sol));
                                                            if all(p_new(:) >= 0)
                                                                if points_actual_j(1, 1) >= x  && points_actual_j(1, 3) >= z_aux
                                                                    points_set_convex{j}.points(2, :) = eval(sol); %%3
                                                                end
                                                                if points_actual_i(3, 1) <= x  && points_actual_i(3, 3) <= z_aux
                                                                    if points_set_convex{i}.tipo == 29
                                                                        points_set_convex{i}.points(4, :) = eval(sol);%%2
                                                                    else
                                                                        dim = size(points_set_convex{i}.points);
                                                                        points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([x,z]==0)
                                            %%%%%%
                                        %    caso = 6;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 6;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_convex{i}, points_set_convex{j},  x, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1 
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 3, 4, 2);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 3, 1, 4, 2);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptXZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 3, 4, 2);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptXZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 3, 1, 4, 2);
                                    %%%%%%%%%%%%%
                                end
                            end
                        end
                        if points_set_convex{i}.tipo == 4
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([1 2 4 3],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([1 2 4 3],:);    
                        end
                        if points_set_convex{j}.tipo == 4
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([1 2 4 3],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([1 2 4 3],:);    
                        end
                        if points_set_convex{i}.tipo == 2 %|| points_set_convex{i}.tipo == 16 %%%% Orden X-Z-Y-W
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([4 2 3 1],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([4 2 3 1],:);    
                        end
                        if points_set_convex{j}.tipo == 2 %|| points_set_convex{j}.tipo == 16
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([4 2 3 1],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([4 2 3 1],:);    
                        end
                        if points_set_convex{i}.tipo == 10 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([3 2 1 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([3 2 1 4],:);    
                        end
                        if points_set_convex{j}.tipo == 10
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([3 2 1 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([3 2 1 4],:);    
                        end
                        if points_set_convex{i}.tipo == 8 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([2 1 3 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([2 1 3 4],:);    
                        end
                        if points_set_convex{j}.tipo == 8
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([2 1 3 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([2 1 3 4],:);    
                        end
                        if points_set_convex{i}.tipo == 26 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([3 1 2 4],:); 
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([3 1 2 4],:);    
                        end
                        if points_set_convex{j}.tipo == 26
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([3 1 2 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([3 1 2 4],:);    %[3 1 2 4]
                        end
                     end                     
                 end

                 %Plano YZ
                 
                  if points_set_convex{i}.tipo == 1 || points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 15 || ...
                     points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 23 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 32 || ...
                     points_set_convex{i}.tipo == 35 || points_set_convex{i}.tipo == 2 || points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || ...
                     points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 18 || points_set_convex{i}.tipo == 30 || points_set_convex{i}.tipo == 13 || ... %points_set_convex{i}.tipo == 8 || ...
                     points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21
                      if points_set_convex{j}.tipo == 1 || points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 15 || ...
                        points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 23 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 32 || ...
                        points_set_convex{j}.tipo == 35 || points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || ...
                        points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 18 || points_set_convex{j}.tipo == 30 || points_set_convex{j}.tipo == 13 || ... %points_set_convex{j}.tipo == 8 || ...
                        points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21
                        N = cross([table_copy(i, 1), table_copy(i, 2), table_copy(i, 3)], [table_copy(j, 1), table_copy(j, 2), table_copy(j, 3)]);
                        if  any(N(:) ~= 0)
                            sol = solve([num2str(table_copy(i, 2)), '*y3+', num2str(table_copy(i, 3)), '*z3 =', num2str(table_copy(i, 4))], ...
                                [num2str(table_copy(j, 2)), '*y3+', num2str(table_copy(j, 3)), '*z3 =', num2str(table_copy(j, 4))], y3, z3);
                        else
                            sol = [];
                        end
                        if points_set_convex{i}.tipo == 10 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([3 2 1 4],:);%[2 1 3 4]
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([3 2 1 4],:);    
                        end
                        if points_set_convex{j}.tipo == 10 
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([3 2 1 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([3 2 1 4],:);    
                        end
                        if points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([1 3 2 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([1 3 2 4],:);    
                        end
                        if points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([1 3 2 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([1 3 2 4],:);    
                        end
                        if points_set_convex{i}.tipo == 9
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([2 3 1 4],:);
                            points_actual_i = points_set_convex{i}.points;
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([2 3 1 4],:);    
                        end
                        if points_set_convex{j}.tipo == 9
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([2 3 1 4],:);
                            points_actual_j = points_set_convex{j}.points;
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([2 3 1 4],:);    
                        end

                        %%%% El vector director ya está considerado al eliminar casos redundantes cuando encontramos
                        %%%% interceptos y los casos acá desarrollados.
                        %%%% Agregar casos degenerados que pasan por el
                        %%%% origen (Caso 9, 10 y 26) ya está                        
                        %m_i = (points_actual_i(1,3)-points_actual_i(2,3))/(points_actual_i(1,2)-points_actual_i(2,2));
                        %if isnan(m_i)
                            m_i = (points_set_nonred{i}.points(1,3)-points_set_nonred{i}.points(2,3))/(points_set_nonred{i}.points(1,2)-points_set_nonred{i}.points(2,2));
                        %end
                        %m_j = (points_actual_j(1,3)-points_actual_j(2,3))/(points_actual_j(1,2)-points_actual_j(2,2));
                        %if isnan(m_j)
                            m_j = (points_set_nonred{j}.points(1,3)-points_set_nonred{j}.points(2,3))/(points_set_nonred{j}.points(1,2)-points_set_nonred{j}.points(2,2));
                        %end
                        if (m_i < 0 && m_j < 0) || ...
                                (points_set_convex{i}.tipo == 36 || points_set_convex{i}.tipo == 33 || points_set_convex{i}.tipo == 30)
                                %((m_i == 0 || m_i == Inf) && (points_set_convex{i}.tipo == 30))
                            if ~isempty(sol)
                                sol = [sol.y3 sol.z3];
                                if isempty(symvar(sol))
                                    y=eval(sol(1)); z=eval(sol(2));
                                    if y > 0 && z > 0                                    
                                        %if points_set_nonred{i}.points(2, 2) <= points_set_nonred{j}.points(2, 2)
                                        if m_i <= m_j
                                            if (points_set_convex{i}.points(2, 1) == 0 || (points_set_convex{i}.points(2, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                                (points_set_convex{j}.points(1, 1) == 0 || (points_set_convex{j}.points(1, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                %if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z %%%11/01
                                                if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z
                                                    if points_set_convex{i}.points(1, 1) == 0 
                                                        %%%%%%
                                                        %[points_set_convex{i}.points, points_set_convex{j}.points] = XParalelPlanesInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                        %    y, z, table_copy, i, j, 1, 2, 4, 3);
                                                        %%%%%%%                                                    
                                                        points_set_convex{i}.points(1, :) = [0 y z];

                                                        %if points_set_convex{j}.points(2, 1) == 0
                                                            %if points_actual_j(2, 2) <= y  && points_actual_j(2, 3) >= z %%%11/01
                                                        %        points_set_convex{j}.points(2, :) = [0 y z];
                                                            %end
                                                        %else
                                                            %if points_actual_j(4, 2) <= y  && points_actual_j(4, 3) >= z %%%11/01
                                                        %        points_set_convex{j}.points(4, :) = [0 y z];
                                                            %end
                                                        %end
                                                    else %if points_set_convex{i}.points(4, 1) == 0 
                                                        %if points_actual_i(4, 2) <= y  && points_actual_i(4, 3) >= z %%%11/01
                                                            %%%%%%
                                                            %[points_set_convex{i}.points, points_set_convex{j}.points] = XParalelPlanesInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                            %    y, z, table_copy, i, j, 4, 2, 1, 3);
                                                            %%%%%%%                                                                                                                
                                                            dim = size(points_set_convex{i}.points);
                                                            if dim(1) == 4
                                                                points_set_convex{i}.points(4, :) = [0 y z];
                                                            else
                                                                points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                            end

                                                            %if points_set_convex{j}.points(2, 1) == 0
                                                                %if points_actual_j(2, 2) <= y  && points_actual_j(2, 3) >= z %%%11/01
                                                            %        points_set_convex{j}.points(2, :) = [0 y z];
                                                                %end
                                                            %else
                                                                %if points_actual_j(4, 2) <= y  && points_actual_j(4, 3) >= z %%%11/01
                                                            %        points_set_convex{j}.points(4, :) = [0 y z];
                                                                %end
                                                            %end
                                                        %else
                                                            %%%%
                                                        %end 
                                                    end
                                                end
                                                
                                                if points_actual_j(2, 2) >= y && points_actual_j(2, 3) <= z
                                                    %if points_actual_j(2, 2) >= y && points_actual_j(2, 3) <= z%%%11/01  Agregar tipo
                                                    %%%%%%
                                                    %[points_set_convex{j}.points, points_set_convex{i}.points] = XParalelPlanesInterceptYZ(handles, points_set_convex{j}.points, points_set_convex{i}.points, ...
                                                    %    y, z, table_copy, i, j, 2, 1, 4, 3);
                                                    %%%%%%%                                                                                                               

                                                    if points_set_convex{j}.points(2, 1) == 0
                                                        %if points_actual_j(2, 2) <= y  && points_actual_j(2, 3) >= z %%%11/01
                                                            points_set_convex{j}.points(2, :) = [0 y z];
                                                        %end
                                                    else
                                                        %if points_actual_j(4, 2) <= y  && points_actual_j(4, 3) >= z %%%11/01
                                                            
                                                            dim = size(points_set_convex{j}.points);
                                                            if dim(1) == 4
                                                                points_set_convex{j}.points(4, :) = [0 y z];
                                                            else
                                                                points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                            end
                                                        %end
                                                    end
                                                end
                                                %%%%%
                                                %%%Desarrollar planos tipo 4
                                                %%%%%
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{i}.tipo == 4 && (points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 23) || ...
                                                    points_set_convex{i}.tipo == 23 && (points_set_convex{j}.tipo == 4 || points_set_convex{j}.tipo == 23)
                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z || ...
                                                        points_actual_j(2, 2) >= y && points_actual_j(2, 3) <= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(3, :) = eval(sol); % 4
                                                                    end
                                                                    if points_actual_j(2, 2) >= y && points_actual_j(2, 3) <= z_aux
                                                                        if points_set_convex{i}.tipo == 23 && points_set_convex{j}.tipo == 23 || ...
                                                                                points_set_convex{i}.tipo == 4 && points_set_convex{j}.tipo == 4
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            end
                                        else
                                            if (points_set_convex{i}.points(1, 1) == 0 || (points_set_convex{i}.points(1, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                               (points_set_convex{j}.points(2, 1) == 0 || (points_set_convex{j}.points(2, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                %if points_actual_i(2, 3) <= z && points_actual_i(2, 2) >= y%%%11/01
                                                if points_actual_i(2, 3) <= z && points_actual_i(2, 2) >= y%%%11/01
                                                    if points_set_convex{i}.points(2, 1) == 0
                                                        %%%%%%
                                                        %[points_set_convex{i}.points, points_set_convex{j}.points] = XParalelPlanesInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                        %    y, z, table_copy, i, j, 2, 1, 4, 3);
                                                        %%%%%%%                                                                                                        
                                                        points_set_convex{i}.points(2, :) = [0 y z];

                                                        %if points_set_convex{j}.points(1, 1) == 0
                                                            %if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z %%%11/01
                                                        %        points_set_convex{j}.points(1, :) = [0 y z];
                                                            %end
                                                        %else
                                                            %if points_actual_j(4, 2) <= y  && points_actual_j(4, 3) >= z %%%11/01
                                                        %        points_set_convex{j}.points(4, :) = [0 y z];
                                                            %end
                                                        %end
                                                    else %if points_set_convex{i}.points(4, 1) == 0
                                                        %if points_actual_i(4, 3) <= z && points_actual_i(4, 2) >= y%%%11/01
                                                            %%%%%%
                                                            %[points_set_convex{i}.points, points_set_convex{j}.points] = XParalelPlanesInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, ...
                                                            %    y, z, table_copy, i, j, 4, 1, 2, 3);
                                                            %%%%%%%                                                                                                                
                                                            dim = size(points_set_convex{i}.points);
                                                            if dim(1) == 4
                                                                points_set_convex{i}.points(4, :) = [0 y z];
                                                            else
                                                                points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                            end
                                                            %if points_set_convex{j}.points(1, 1) == 0
                                                                %if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z %%%11/01
                                                            %        points_set_convex{j}.points(1, :) = [0 y z];
                                                                %end
                                                            %else
                                                                %if points_actual_j(4, 2) <= y  && points_actual_j(4, 3) >= z %%%11/01
                                                            %        points_set_convex{j}.points(4, :) = [0 y z];
                                                                %end
                                                            %end
                                                        %else
                                                            %%%%%%
                                                        %end  
                                                    end
                                                end
                                                
                                                if points_actual_j(1, 2) <= y && points_actual_j(1, 3) >= z
                                                    %if points_actual_j(1, 2) >= y && points_actual_j(1, 3) <= z%%%11/01  Agregar tipo
                                                    %%%%%%
                                                    %[points_set_convex{j}.points, points_set_convex{i}.points] = XParalelPlanesInterceptYZ(points_set_convex{j}.points, points_set_convex{i}.points, ...
                                                    %    y, z, table_copy, i, j, 1, 2, 4, 3);
                                                    %%%%%%%

                                                    if points_set_convex{j}.points(1, 1) == 0
                                                        %if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z %%%11/01
                                                            points_set_convex{j}.points(1, :) = [0 y z];
                                                        %end
                                                    else
                                                        %if points_actual_j(4, 2) <= y  && points_actual_j(4, 3) >= z %%%11/01                                                            
                                                            dim = size(points_set_convex{j}.points);
                                                            if dim(1) == 4
                                                                points_set_convex{j}.points(4, :) = [0 y z];
                                                            else
                                                                points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                            end
                                                        %end
                                                    end
                                                end
                                                %%%%%
                                                %%%Desarrollar planos tipo 4
                                                %%%%%
                                                %%%%%% Modificar (Revisar)
                                                if points_set_convex{j}.tipo == 4 && (points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 23) || ...
                                                    points_set_convex{j}.tipo == 23 && (points_set_convex{i}.tipo == 4 || points_set_convex{i}.tipo == 23)
                                                    if points_actual_j(1, 2) <= y && points_actual_j(1, 3) >= z || ...
                                                        points_actual_i(2, 2) >= y && points_actual_i(2, 3) <= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %y = y + 1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(1, 2) <= y && points_actual_j(1, 3) >= z_aux
                                                                        points_set_convex{j}.points(3, :) = eval(sol); %
                                                                    end
                                                                    if points_actual_i(2, 2) >= y && points_actual_i(2, 3) <= z_aux
                                                                        if points_set_convex{i}.tipo == 23 && points_set_convex{j}.tipo == 23 || ...
                                                                                points_set_convex{i}.tipo == 4 && points_set_convex{j}.tipo == 4
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            end
                                        end                                    
                                    else
                                        %%%%%%
                                        caso = 1;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%% 
                                        if IsPlaneCrossIntercept == 1
                                            if points_set_nonred{i}.points(2, 2) <= points_set_nonred{j}.points(2, 2)
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 2, 4, 3);
                                            else
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 1, 4, 3);  
                                            end
                                        end
                                    end                 
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptYZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%                                    
                                end
                            else
                                if points_set_nonred{i}.points(2, 2) <= points_set_nonred{j}.points(2, 2)
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 2, 4, 3);
                                else
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 1, 4, 3);  
                                end
                            end
                        elseif ((m_i < 0 ) && ((m_j > 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(2) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 26)))) ...
                                || ...
                                ((m_j < 0 ) && ((m_i > 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(2) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 26))))
                            %%%%%%%
                            if ~isempty(sol)
                                sol = [sol.y3 sol.z3];
                                if isempty(symvar(sol))
                                    y=eval(sol(1)); z=eval(sol(2));
                                    if (y >= 0 && z > 0) || (y > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [1 0 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [1 0 0]);
                                        
                                        if m_i < 0
                                            %if (points_set_convex{i}.points(1, 1) == 0 || (points_set_convex{i}.points(1, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %   (points_set_convex{j}.points(2, 1) == 0 || (points_set_convex{j}.points(2, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) <= z
                                                    if points_set_convex{i}.points(2, 1) == 0 
                                                        points_set_convex{i}.points(2, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                if points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z
                                                    if points_set_convex{j}.points(1, 1) == 0 
                                                        points_set_convex{j}.points(1, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C2) < e                                                        
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                %%%%%% Modificar
                                                if points_set_convex{i}.tipo == 4 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 9 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32) || ...
                                                        points_set_convex{i}.tipo == 23 && (points_set_convex{j}.tipo == 2 || points_set_convex{j}.tipo == 9 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32)
                                                    if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) <= z || ...
                                                        points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) <= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(2, 1) == 0 || (points_set_convex{i}.points(2, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %    (points_set_convex{j}.points(1, 1) == 0 || (points_set_convex{j}.points(1, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) <= z
                                                    if points_set_convex{j}.points(2, 1) == 0 
                                                        points_set_convex{j}.points(2, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                if points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z
                                                    if points_set_convex{i}.points(1, 1) == 0 
                                                        points_set_convex{i}.points(1, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 2 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                %%%%%% Modificar
                                                if points_set_convex{j}.tipo == 4 && (points_set_convex{i}.tipo == 2 || points_set_convex{i}.tipo == 9 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 32) || ...
                                                        points_set_convex{j}.tipo == 23 && (points_set_convex{i}.tipo == 2 || points_set_convex{i}.tipo == 9 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 32)
                                                    if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) <= z || ...
                                                       points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z         
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) <= z_aux
                                                                        points_set_convex{j}.points(4, :) = eval(sol);
                                                                    end

                                                                    if points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z_aux
                                                                        if points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 32
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                             %   [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                             %       i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([y,z]==0)
                                            %%%%%%
                                        %    caso = 2;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end                                                                                
                                    else
                                        %%%%%%
                                        caso = 2;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            %%%%%%%
                                            if m_i < 0
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptYZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                %%%%%%
                                if m_i < 0
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 1, 4, 3);  
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 2, 4, 3);
                                    %%%%%%%%%%%%%
                                end
                            end 
                        elseif (m_i < 0 && (m_j >= 0 && (points_set_convex{j}.tipo ==3 || ...
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 20 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(2) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 26)))) ...
                                ||...
                               (m_j < 0 && (m_i >= 0 && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 20 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(2) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 26))))
                            %%%%%                            
                            if ~isempty(sol)
                                sol = [sol.y3 sol.z3];
                                if isempty(symvar(sol))
                                    y=eval(sol(1)); z=eval(sol(2));
                                    if (y >= 0 && z > 0) || (y > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [1 0 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [1 0 0]);
                                        if m_i < 0
                                            %if (points_set_convex{i}.points(2, 1) == 0 || (points_set_convex{i}.points(2, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %    (points_set_convex{j}.points(1, 1) == 0 || (points_set_convex{j}.points(1, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z
                                                    if points_set_convex{i}.points(1, 1) == 0 
                                                        points_set_convex{i}.points(1, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z
                                                    if points_set_convex{j}.points(2, 1) == 0 
                                                        points_set_convex{j}.points(2, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                %%%%%%
                                                if points_set_convex{i}.tipo == 4 && (points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 9 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 35) || ...
                                                        points_set_convex{i}.tipo == 23 && (points_set_convex{j}.tipo == 3 || points_set_convex{j}.tipo == 9 || ...
                                                        points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 35)
                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z || ...
                                                        points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z+1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(1, 1) == 0 || (points_set_convex{i}.points(1, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                             %  (points_set_convex{j}.points(2, 1) == 0 || (points_set_convex{j}.points(2, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z
                                                    if points_set_convex{j}.points(1, 1) == 0 
                                                        points_set_convex{j}.points(1, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z
                                                    if points_set_convex{i}.points(2, 1) == 0 
                                                        points_set_convex{i}.points(2, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 3 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                %%%%%%
                                                if points_set_convex{j}.tipo == 4 && (points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 9 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 35) || ...
                                                        points_set_convex{j}.tipo == 23 && (points_set_convex{i}.tipo == 3 || points_set_convex{i}.tipo == 9 || ...
                                                        points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 35)
                                                    if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z    || ...
                                                       points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z  
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{j}.points(4, :) = [0 y z];
                                                                %z = z + 1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z_aux
                                                                        points_set_convex{j}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z_aux
                                                                        if points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 26 || points_set_convex{i}.tipo == 35
                                                                            points_set_convex{i}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{i}.points);
                                                                            points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %       i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([y,z]==0)
                                            %%%%%%
                                        %    caso = 3;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 3;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            %%%%%%%
                                            if m_i < 0
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptYZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < 0
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 2, 4, 3);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 1, 4, 3);  
                                    %%%%%%%%%%%%%
                                end
                            end                             
                        elseif ((m_i >= 0 && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 20 ||  points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(2) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))) && ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(2) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))))) ...
                                ||...
                               ((m_j >= 0 && (points_set_convex{j}.tipo ==3 || ...
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 20 ||  points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(2) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))) && ...
                               (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(2) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))))
                            if ~isempty(sol)
                                sol = [sol.y3 sol.z3];
                                if isempty(symvar(sol))
                                    y=eval(sol(1)); z=eval(sol(2));
                                    if (y > 0 && z > 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [1 0 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [1 0 0]);
                                        
                                        if m_i < m_j
                                            %if (points_set_convex{i}.points(1, 1) == 0 || (points_set_convex{i}.points(1, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %   (points_set_convex{j}.points(2, 1) == 0 || (points_set_convex{j}.points(2, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z
                                                    points_set_convex{i}.points(2, :) = [0 y z];
                                                    
                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end
                                                
                                                if points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z
                                                    points_set_convex{j}.points(1, :) = [0 y z];
                                                    
                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end
                                                
                                                if points_set_convex{i}.tipo == 7 && (points_set_convex{j}.tipo == 35 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 6) || ...
                                                    points_set_convex{i}.tipo == 24 && (points_set_convex{j}.tipo == 35 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 6) || ...
                                                    points_set_convex{i}.tipo == 32 && (points_set_convex{j}.tipo == 35 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 6)
                                                    if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z || ...
                                                        points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z+1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(2, 1) == 0 || (points_set_convex{i}.points(2, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %    (points_set_convex{j}.points(1, 1) == 0 || (points_set_convex{j}.points(1, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z
                                                    points_set_convex{i}.points(1, :) = [0 y z];
                                                    
                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end                                                

                                                if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z
                                                    points_set_convex{j}.points(2, :) = [0 y z];
                                                    
                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end
                                                
                                                if points_set_convex{j}.tipo == 7 && (points_set_convex{i}.tipo == 35 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 6) || ...
                                                    points_set_convex{j}.tipo == 24 && (points_set_convex{i}.tipo == 35 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 6) || ...
                                                    points_set_convex{j}.tipo == 32 && (points_set_convex{i}.tipo == 35 || points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 6)
                                                    if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z || ...
                                                        points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z+1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                    else
                                        %%%%%%
                                        caso = 4;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            %%%%%%
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptYZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 1, 4, 3);  
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 2, 4, 3);
                                    %%%%%%%%%%%%%
                                end
                            end                                                                     
                        elseif ((m_i >= 0 && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 20 ||  points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(2) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))) && ...
                                (m_j >= 0 && (points_set_convex{j}.tipo ==3 || ...
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 20 ||  points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(2) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) ...
                                || ...
                               ((m_j >= 0 && (points_set_convex{j}.tipo ==3 || ...
                                points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 18 || ...
                                points_set_convex{j}.tipo == 20 ||  points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35 || ...
                                (d_j(2) <= 0 && d_j(3) >= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))) && ...
                                (m_i >= 0 && (points_set_convex{i}.tipo ==3 || ...
                                points_set_convex{i}.tipo == 7 || points_set_convex{i}.tipo == 18 || ...
                                points_set_convex{i}.tipo == 20 ||  points_set_convex{i}.tipo == 25 || points_set_convex{i}.tipo == 35 || ...
                                (d_i(2) <= 0 && d_i(3) >= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26)))) %%%% Agregar tipo a ambos planos
                            %%%%%%
                            if ~isempty(sol)
                                sol = [sol.y3 sol.z3];
                                if isempty(symvar(sol))
                                    y=eval(sol(1)); z=eval(sol(2));
                                    if (y >= 0 && z > 0) || (y > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [1 0 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [1 0 0]);
                                        
                                        if m_i < m_j
                                            %if (points_set_convex{i}.points(2, 1) == 0 || (points_set_convex{i}.points(2, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %    (points_set_convex{j}.points(1, 1) == 0 || (points_set_convex{j}.points(1, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) <= z
                                                    if points_set_convex{i}.points(1, 1) == 0 
                                                        points_set_convex{i}.points(1, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z]; 
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z
                                                    if points_set_convex{j}.points(2, 1) == 0 
                                                        points_set_convex{j}.points(2, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end
                                                
                                                if points_set_convex{i}.tipo == 7 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35) || ...
                                                    points_set_convex{i}.tipo == 25 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35) || ...
                                                    points_set_convex{i}.tipo == 35 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35)
                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) <= z || ...
                                                        points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z+1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(2, 2) >= y  && points_actual_j(2, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(1, 1) == 0 || (points_set_convex{i}.points(1, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %   (points_set_convex{j}.points(2, 1) == 0 || (points_set_convex{j}.points(2, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                                if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) <= z
                                                    if points_set_convex{j}.points(1, 1) == 0 
                                                        points_set_convex{j}.points(1, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{j}.points);
                                                        points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C2) < e
                                                    %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z]; 
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end

                                                if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z
                                                    if points_set_convex{i}.points(2, 1) == 0 
                                                        points_set_convex{i}.points(2, :) = [0 y z];
                                                    else
                                                        dim = size(points_set_convex{i}.points);
                                                        points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                    end

                                                    %if abs(C) < e
                                                    %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z]; 
                                                    %else
                                                    %    s = get(handles.listbox_operations, 'string');
                                                    %    set(handles.listbox_operations, 'string', ...
                                                    %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 5 de tipo de plano en desarrollo.')); 
                                                    %end
                                                end
                                                
                                                if points_set_convex{i}.tipo == 7 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35) || ...
                                                    points_set_convex{i}.tipo == 25 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35) || ...
                                                    points_set_convex{i}.tipo == 35 && (points_set_convex{j}.tipo == 7 || points_set_convex{j}.tipo == 25 || points_set_convex{j}.tipo == 35)
                                                    if points_actual_i(2, 2) >= y  && points_actual_i(2, 3) >= z || ...
                                                        points_actual_j(1, 2) <= y  && points_actual_j(1, 3) <= z
                                                        z_aux = z;
                                                        syms x4 y4 z4 z;
                                                        sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                    [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                        if ~isempty(sol)
                                                            sol = [sol.x4 sol.y4 sol.z4];
                                                            if ~isempty(symvar(sol))
                                                                %points_set_convex{i}.points(4, :) = [0 y z];
                                                                %z = z+1;
                                                                z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                                z = z_sol; p_new = eval(eval(sol));
                                                                if all(p_new(:) >= 0)
                                                                    if points_actual_i(1, 2) <= y  && points_actual_i(1, 3) >= z_aux
                                                                        points_set_convex{i}.points(4, :) = eval(sol);
                                                                    end
                                                                    if points_actual_j(1, 2) <= y  && points_actual_j(1, 3) >= z_aux
                                                                        if points_set_convex{j}.tipo == 35
                                                                            points_set_convex{j}.points(4, :) = eval(sol);
                                                                        else
                                                                            dim = size(points_set_convex{j}.points);
                                                                            points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            %else
                                            %    %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([y,z]==0)
                                            %%%%%%
                                        %    caso = 5;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 5;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            %%%%%%%
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptYZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 2, 4, 3);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 1, 4, 3);  
                                    %%%%%%%%%%%%%
                                end
                            end
                        elseif (((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(2) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))) && ...
                                ((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(2) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26)))) ...
                                || ...
                                (((m_j >= 0 || m_j == Inf) && (points_set_convex{j}.tipo ==2 || ...
                                points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 13 || ...
                                points_set_convex{j}.tipo == 15 || points_set_convex{j}.tipo == 24 || points_set_convex{j}.tipo == 32 || ...
                                (d_j(2) >= 0 && d_j(3) <= 0) && (points_set_convex{j}.tipo == 9 || points_set_convex{j}.tipo == 10 || points_set_convex{j}.tipo == 16  || points_set_convex{j}.tipo == 21 || points_set_convex{j}.tipo == 26))) && ...
                                ((m_i >= 0 || m_i == Inf) && (points_set_convex{i}.tipo ==2 || ...
                                points_set_convex{i}.tipo == 6 || points_set_convex{i}.tipo == 13 || ...
                                points_set_convex{i}.tipo == 15 || points_set_convex{i}.tipo == 24 || points_set_convex{i}.tipo == 32 || ...
                                (d_i(2) >= 0 && d_i(3) <= 0) && (points_set_convex{i}.tipo == 9 || points_set_convex{i}.tipo == 10 || points_set_convex{i}.tipo == 16  || points_set_convex{i}.tipo == 21 || points_set_convex{i}.tipo == 26))))
                            if ~isempty(sol)
                                sol = [sol.y3 sol.z3];
                                if isempty(symvar(sol))
                                    y=eval(sol(1)); z=eval(sol(2));
                                    if (y >= 0 && z > 0) || (y > 0 && z >= 0)
                                        %N = cross(points_actual_i(1,:)-points_actual_i(2,:), points_actual_i(1,:)-points_actual_i(3,:));
                                        %C = dot(N, [1 0 0]);
                                        %N2 = cross(points_actual_j(1,:)-points_actual_j(2,:), points_actual_j(1,:)-points_actual_j(3,:));
                                        %C2 = dot(N2, [1 0 0]);
                                        
                                        if m_i < m_j
                                            %if (points_set_convex{i}.points(2, 1) == 0 || (points_set_convex{i}.points(2, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %    (points_set_convex{j}.poin
                                            %    ts(1, 1) == 0 || (points_set_convex{j}.points(1, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                            if points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z
                                                if points_set_convex{i}.points(1, 1) == 0 
                                                    points_set_convex{i}.points(1, :) = [0 y z];
                                                else
                                                    dim = size(points_set_convex{i}.points);
                                                    points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                end

                                                %if abs(C) < e
                                                %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z]; 
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.')); 
                                                %end
                                            end

                                            if points_actual_j(2, 2) <= y  && points_actual_j(2, 3) <= z
                                                if points_set_convex{j}.points(2, 1) == 0 
                                                    points_set_convex{j}.points(2, :) = [0 y z];
                                                else
                                                    dim = size(points_set_convex{j}.points);
                                                    points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                end

                                                %if abs(C2) < e
                                                %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z];
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.')); 
                                                %end
                                            end
                                                
                                            if points_set_convex{i}.tipo == 6 && (points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32) || ...
                                                points_set_convex{i}.tipo == 26 && (points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32) || ...
                                                points_set_convex{i}.tipo == 32 && (points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32)                                                
                                                if points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z || ...
                                                    points_actual_j(2, 2) <= y  && points_actual_j(2, 3) <= z
                                                    z_aux = z;
                                                    syms x4 y4 z4 z;
                                                    sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                    if ~isempty(sol)
                                                        sol = [sol.x4 sol.y4 sol.z4];
                                                        if ~isempty(symvar(sol))
                                                            %points_set_convex{i}.points(4, :) = [0 y z];
                                                            %z = z+1;
                                                            z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                            z = z_sol; p_new = eval(eval(sol));
                                                            if all(p_new(:) >= 0)
                                                                if points_actual_i(1, 2) >= y  && points_actual_i(1, 3) >= z_aux
                                                                    points_set_convex{i}.points(4, :) = eval(sol); %%3
                                                                end
                                                                if points_actual_j(2, 2) <= y  && points_actual_j(2, 3) <= z_aux
                                                                    if points_set_convex{j}.tipo == 32
                                                                        points_set_convex{j}.points(4, :) = eval(sol); %%3
                                                                    else
                                                                        dim = size(points_set_convex{j}.points);
                                                                        points_set_convex{j}.points(dim(1)+1, :) = eval(sol);
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            %end
                                        else
                                            %if (points_set_convex{i}.points(1, 1) == 0 || (points_set_convex{i}.points(1, 1) > 0 && points_set_convex{i}.points(4, 1) == 0 )) && ...%%%% OJO
                                            %   (points_set_convex{j}.point
                                            %   s(2, 1) == 0 || (points_set_convex{j}.points(2, 1) > 0 && points_set_convex{j}.points(4, 1) == 0 )) %%%% OJO
                                            if points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z
                                                if points_set_convex{j}.points(1, 1) == 0 
                                                    points_set_convex{j}.points(1, :) = [0 y z];
                                                else
                                                    dim = size(points_set_convex{j}.points);
                                                    points_set_convex{j}.points(dim(1)+1, :) = [0 y z];
                                                end

                                                %if abs(C2) < e
                                                %    points_set_convex{j}.points(4, :) = [points_set_convex{j}.points(4, 1) y z]; 
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.')); 
                                                %end
                                            end

                                            if points_actual_i(2, 2) <= y  && points_actual_i(2, 3) <= z
                                                if points_set_convex{i}.points(2, 1) == 0 
                                                    points_set_convex{i}.points(2, :) = [0 y z];
                                                else
                                                    dim = size(points_set_convex{i}.points);
                                                    points_set_convex{i}.points(dim(1)+1, :) = [0 y z];
                                                end

                                                %if abs(C) < e
                                                %    points_set_convex{i}.points(4, :) = [points_set_convex{i}.points(4, 1) y z]; 
                                                %else
                                                %    s = get(handles.listbox_operations, 'string');
                                                %    set(handles.listbox_operations, 'string', ...
                                                %        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 6 de tipo de plano en desarrollo.')); 
                                                %end
                                            end
                                                
                                            if points_set_convex{i}.tipo == 6 && (points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32) || ...
                                                points_set_convex{i}.tipo == 26 && (points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32) || ...
                                                points_set_convex{i}.tipo == 32 && (points_set_convex{j}.tipo == 6 || points_set_convex{j}.tipo == 26 || points_set_convex{j}.tipo == 32)
                                                if points_actual_i(2, 2) <= y  && points_actual_i(2, 3) <= z || ...
                                                    points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z
                                                    z_aux = z;
                                                    syms x4 y4 z4 z;
                                                    sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                                                                [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
                                                    if ~isempty(sol)
                                                        sol = [sol.x4 sol.y4 sol.z4];
                                                        if ~isempty(symvar(sol))
                                                            %points_set_convex{i}.points(4, :) = [0 y z];
                                                            %z = z+1;
                                                            z_sol = solve(sol(1)-max([points_actual_i(:,1);points_actual_j(:,1)]),z);
                                                            z = z_sol; p_new = eval(eval(sol));
                                                            if all(p_new(:) >= 0)
                                                                if points_actual_i(2, 2) <= y  && points_actual_i(2, 3) <= z_aux
                                                                    if points_set_convex{i}.tipo == 32
                                                                        points_set_convex{i}.points(4, :) = eval(sol); %%3                                                                    
                                                                    else
                                                                        dim = size(points_set_convex{i}.points);
                                                                        points_set_convex{i}.points(dim(1)+1, :) = eval(sol);
                                                                    end                                                                    
                                                                end
                                                                if points_actual_j(1, 2) >= y  && points_actual_j(1, 3) >= z_aux
                                                                    points_set_convex{j}.points(4, :) = eval(sol); %%3
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                            %else
                                                %%%%%%%%%%%%%
                                            %    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                            %        i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            %end
                                        end
                                        %if any([y,z]==0)
                                            %%%%%%
                                        %    caso = 6;
                                        %    [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                        %            i, j, caso); %#ok<NASGU>
                                            %%%%%%%
                                        %end
                                    else
                                        %%%%%%
                                        caso = 6;
                                        [points_set_convex{i}, points_set_convex{j}, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_convex{i}, points_set_convex{j},  y, z, table_copy, ...
                                                    i, j, caso);
                                        %%%%%%%
                                        if IsPlaneCrossIntercept == 1
                                            %%%%%%
                                            if m_i < m_j
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 1, 2, 4, 3);
                                                %%%%%%%%%%%%%
                                            else
                                                %%%%%%%%%%%%%
                                                [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                                    i, j, 2, 1, 4, 3);  
                                                %%%%%%%%%%%%%
                                            end
                                        end
                                    end
                                else
                                    %%%%%%
                                    [points_set_convex{i}, points_set_convex{j}, table_copy, pasar] = LineInterceptYZ(points_set_convex{i}, points_set_convex{j}, table_copy, ...
                                                    i, j, pasar);
                                    %%%%%%%
                                end
                            else
                                if m_i < m_j
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 1, 2, 4, 3);
                                    %%%%%%%%%%%%%
                                else
                                    %%%%%%%%%%%%%
                                    [points_set_convex{i}.points, points_set_convex{j}.points, table_copy] = PlaneCrossInterceptYZ(handles, points_set_convex{i}.points, points_set_convex{j}.points, table_copy, ...
                                        i, j, 2, 1, 4, 3);  
                                    %%%%%%%%%%%%%
                                end
                            end                            
                        end
                        if points_set_convex{i}.tipo == 10 %%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([3 2 1 4],:);                            
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([3 2 1 4],:);    
                        end
                        if points_set_convex{j}.tipo == 10
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([3 2 1 4],:);                            
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([3 2 1 4],:);    
                        end 
                        if points_set_convex{i}.tipo == 16 || points_set_convex{i}.tipo == 21%%%% Orden
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([1 3 2 4],:);                            
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([1 3 2 4],:);    
                        end
                        if points_set_convex{j}.tipo == 16 || points_set_convex{j}.tipo == 21
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([1 3 2 4],:);                            
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([1 3 2 4],:);    
                        end
                        if points_set_convex{i}.tipo == 9
                            points_set_convex{i}.points(1:4, :) = points_set_convex{i}.points([3 1 2 4],:);                            
                            points_set_nonred{i}.points(1:4, :) = points_set_nonred{i}.points([3 1 2 4],:);    
                        end
                        if points_set_convex{j}.tipo == 9
                            points_set_convex{j}.points(1:4, :) = points_set_convex{j}.points([3 1 2 4],:);
                            points_set_nonred{j}.points(1:4, :) = points_set_nonred{j}.points([3 1 2 4],:);    
                        end
                      end                   
                  end
                  if i < Dimension(1)
                      if any(str2num(rats(points_set_convex{i}.points*table_copy(j, 1:3)')) > table_copy(j, 4)) && ...
                          any(str2num(rats(points_set_convex{i}.points*table_copy(j, 1:3)')) < table_copy(j, 4)) && ...
                            all(str2num(rats(points_set_convex{j}.points*table_copy(i, 1:3)')) <= table_copy(i, 4)) %#ok<ST2NM> %%%OJO con la precisión
                                
                                nonintersect = str2num(rats(points_set_convex{i}.points*table_copy(j, 1:3)')) > table_copy(j, 4); %#ok<ST2NM> %& ...
                                    %all(points_set_convex{i}.points' > 0)';  
                                points_set_convex_i_aux = points_set_convex{i}.points; 
                                
                                %points_set_convex{i}.points = points_set_convex{i}.points(~nonintersect', :);
                                i_first= find(~nonintersect',1);
                                if ~isempty(i_first) 
                                    copia_points = [];
                                    for k = 1:sum(nonintersect)
                                        copia_points = [copia_points;points_set_convex{i}.points(i_first, :)];
                                    end
                                    points_set_convex{i}.points(nonintersect', :) = copia_points;
                                else
                                    points_set_convex{i}.points = [];
                                end
                                points_i = unique(points_set_convex{i}.points, 'rows');    
                                dim_aux = size(points_i);
                                if dim_aux(1) > 0 && dim_aux(1) < 3
                                    points_set_convex{i}.points = [points_i; ...
                                        points_set_convex{j}.points(str2num(rats(points_set_convex{j}.points*table_copy(i, 1:3)')) == table_copy(i, 4),:)]; %#ok<ST2NM>
                                    points_i = unique(points_set_convex{i}.points, 'rows');
                                    dim_aux = size(points_i);
                                    if dim_aux(1) == 3
                                        points_set_convex{i}.points = [points_i;points_i(3,:)];
                                    elseif dim_aux(1) < 3
                                        points_set_convex{i}.points = points_set_convex_i_aux;
                                    end
                                elseif isempty(points_set_convex{i}.points)
                                    points_set_convex{i}.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_convex{i}.points;
                                    set(handles_surf(i), 'Visible', 'off');
                                    delete(handles_surf(i));
                                    set(handles_norm(i), 'Visible', 'off');
                                    delete(handles_norm(i));
                                    handles_surf(i) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
                                    handles_norm(i) = quiver3(0, 0, 0, table_copy(i, 1), table_copy(i, 2), table_copy(i, 3), 'Color', 'red');
                                end

                      elseif any(str2num(rats(points_set_convex{j}.points*table_copy(i, 1:3)')) > table_copy(i, 4)) && ...
                          any(str2num(rats(points_set_convex{j}.points*table_copy(i, 1:3)')) < table_copy(i, 4)) && ...
                              all(str2num(rats(points_set_convex{i}.points*table_copy(j, 1:3)')) <= table_copy(j, 4)) %#ok<ST2NM> %%%OJO con la precisión
                            
                                nonintersect = str2num(rats(points_set_convex{j}.points*table_copy(i, 1:3)')) > table_copy(i, 4); %#ok<ST2NM> %& ...
                                        %all(points_set_convex{j}.points' > 0)';  
                                points_set_convex_j_aux = points_set_convex{j}.points; 
                                
                                %points_set_convex{j}.points = points_set_convex{j}.points(~nonintersect', :);
                                j_first= find(~nonintersect',1);
                                if ~isempty(j_first)
                                    copia_points = [];
                                    for k = 1:sum(nonintersect)
                                        copia_points = [copia_points;points_set_convex{j}.points(j_first, :)];
                                    end
                                    points_set_convex{j}.points(nonintersect', :) = copia_points;                                    
                                else
                                    points_set_convex{j}.points = [];
                                end
                                points_j = unique(points_set_convex{j}.points, 'rows');
                                dim_aux = size(points_j);
                                if dim_aux(1) > 0 && dim_aux(1) < 3
                                    points_set_convex{j}.points = [points_j; ...
                                        points_set_convex{i}.points(str2num(rats(points_set_convex{i}.points*table_copy(j, 1:3)')) == table_copy(j, 4),:)]; %#ok<ST2NM>
                                    points_j = unique(points_set_convex{j}.points, 'rows');
                                    dim_aux = size(points_j);                                    
                                    if dim_aux(1) == 3
                                        points_set_convex{j}.points = [points_j;points_j(3,:)];
                                    elseif dim_aux(1) < 3
                                        points_set_convex{j}.points = points_set_convex_j_aux;
                                    end
                                elseif isempty(points_set_convex{j}.points)
                                    points_set_convex{j}.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_convex{j}.points;
                                    set(handles_surf(j), 'Visible', 'off');
                                    delete(handles_surf(j));
                                    set(handles_norm(j), 'Visible', 'off');
                                    delete(handles_norm(j));
                                    handles_surf(j) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
                                    handles_norm(j) = quiver3(0, 0, 0, table_copy(j, 1), table_copy(j, 2), table_copy(j, 3), 'Color', 'red');
                                end
                      elseif all(str2num(rats(points_set_convex{j}.points*table_copy(i, 1:3)')) >= table_copy(i, 4)) && ...
                          all(str2num(rats(points_set_convex{i}.points*table_copy(j, 1:3)')) <= table_copy(j, 4)) %#ok<ST2NM>
                            points_set_convex{j}.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_convex{j}.points;
                            set(handles_surf(j), 'Visible', 'off');
                            delete(handles_surf(j));
                            set(handles_norm(j), 'Visible', 'off');
                            delete(handles_norm(j));
                            handles_surf(j) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
                            handles_norm(j) = quiver3(0, 0, 0, table_copy(j, 1), table_copy(j, 2), table_copy(j, 3), 'Color', 'red');
                      elseif all(str2num(rats(points_set_convex{i}.points*table_copy(j, 1:3)')) >= table_copy(j, 4)) && ...
                          all(str2num(rats(points_set_convex{j}.points*table_copy(i, 1:3)')) <= table_copy(i, 4)) %#ok<ST2NM>
                            points_set_convex{i}.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_convex{i}.points;
                            set(handles_surf(i), 'Visible', 'off');
                            delete(handles_surf(i));
                            set(handles_norm(i), 'Visible', 'off');
                            delete(handles_norm(i));
                            handles_surf(i) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
                            handles_norm(i) = quiver3(0, 0, 0, table_copy(i, 1), table_copy(i, 2), table_copy(i, 3), 'Color', 'red');
                      end
                  end
            end
            
            p1 = points_set_convex{i}.points;
            if i > Dimension(1)
                if ~all(p1(1, :) == p1(2, :)) && ~all(p1(1, :)==p1(3, :)) && ~all(p1(2, :) == p1(3, :))
                    if ~all(points_set_convex{i}.points(:) == points_set_nonred{i}.points(:))
                        points_i = unique(points_set_convex{i}.points, 'rows');
                        dim_i = size(points_i);
                        points_all = unique(all_points, 'rows');
                        dim_all = size(points_all);
                        if dim_i(1) >= 3 && dim_all(1) >= 3
                            d_i = cross(points_i(1, :) - points_i(2, :), points_i(1, :) - points_i(3, :));
                            d_all = cross(points_all(1, :) - points_all(2, :), points_all(1, :) - points_all(3, :));
                            if points_set_convex{i}.tipo == 36
                                if d_i(3) > 0
                                    d_i = -d_i;
                                end
                                if d_all(3) > 0
                                    d_all = -d_all;
                                end
                            elseif points_set_convex{i}.tipo == 33
                                if d_i(2) > 0
                                    d_i = -d_i;
                                end
                                if d_all(2) > 0
                                    d_all = -d_all;
                                end
                            elseif points_set_convex{i}.tipo == 30
                                if d_i(1) > 0
                                    d_i = -d_i;
                                end
                                if d_all(1) > 0
                                    d_all = -d_all;
                                end
                            end
                            if all(str2num(rats(points_i*d_all')) >= str2num(rats(points_all(1,:)*d_all'))) && ...
                                    any(str2num(rats(points_i*d_all')) > str2num(rats(points_all(1,:)*d_all'))) %#ok<ST2NM>
                                all_points = [points_set_convex{i}.points];
                            elseif ~(all(str2num(rats(points_all*d_i')) >= str2num(rats(points_i(1,:)*d_i'))) && ...
                                    any(str2num(rats(points_all*d_i')) > str2num(rats(points_i(1,:)*d_i')))) %#ok<ST2NM>
                                all_points = [all_points;points_set_convex{i}.points];   
                            end
                        elseif isempty(all_points)
                            all_points = points_set_convex{i}.points;   
                        end
                    end
                end
            end
        end
        p = points_set_convex{i}.points;
        if i > Dimension(1)
            if ~isempty(all_points)
                p = all_points;
                pasar = 1;
            end
        end
        if (all(p(1, :) == p(2, :)) || all(p(1, :)==p(3, :)) || all(p(2, :) == p(3, :))) && all(p(3, :) == p(4, :))
            pasar = 0;
        end

        if  any(p(:)~=0) && pasar                               
            mean_point = 0.3*p(1,:) + 0.3*p(2,:) + 0.4*p(3,:);
            if i <= Dimension(1)
                K = convex_hull(p, points_set_nonred{i}.points);
            else
                K = convex_hull(p, p);
            end
            K_dim = size(K);
            if i < Dimension(1)
                set(handles_surf(i), 'Visible', 'off');
                delete(handles_surf(i));
                set(handles_norm(i), 'Visible', 'off');
                delete(handles_norm(i));

                handles_surf(i) = patch('xdata',p(K(1:K_dim(1)-1),1), 'ydata',p(K(1:K_dim(1)-1),2), ...
                        'zdata', p(K(1:K_dim(1)-1),3), 'FaceColor', 'g'); hold on;
                handles_norm(i) = quiver3(handles.axes_simplex3D, mean_point(1),mean_point(2),mean_point(3), table_copy(i, 1), table_copy(i, 2), table_copy(i, 3), 'Color', 'red');
            else
                set(handles_surf(i), 'Visible', 'off');
                delete(handles_surf(i));

                if strcmp(get(handles.Restriction_nonnegativity, 'Checked'), 'off')
                    visibility = 'off';
                else
                    visibility = 'on';
                end
                handles_surf(i) = patch('xdata',p(K(1:K_dim(1)-1),1), 'ydata',p(K(1:K_dim(1)-1),2), ...
                        'zdata', p(K(1:K_dim(1)-1),3), 'FaceColor', 'b', 'visible', visibility, 'facealpha', .5); hold on;                        
            end
        end
    end
end

set(handles.Uncut_planes, 'Enable', 'on');
set(handles.uipushtool1, 'Enable', 'off');
set(handles.uitoggletool6, 'enable', 'off');


function [points_set_i_con,points_set_j_con, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXY(points_set_i_con, points_set_j_con, x, y, table_copy, ...
    i, j, caso)
global points_set_nonred m_i m_j handles_surf handles_norm;
e = 0.00005; %%%% Cambiar si se cambia precisión
IsPlaneCrossIntercept = 0;

d_i = table_copy(i, 1:3); % vector director de i
d_j = table_copy(j, 1:3); % vector director de j
if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))
    
    %points_set_i_con.points = points_set_j_con.points;
    
    %table_copy(i, :) = table_copy(j, :);
    %points_set_nonred{i} = points_set_nonred{j};
    %points_set_i_con.tipo =points_set_j_con.tipo;
    points_set_i_con.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_i_con.points;
    set(handles_surf(i), 'Visible', 'off');
    delete(handles_surf(i));
    set(handles_norm(i), 'Visible', 'off');
    delete(handles_norm(i));
    handles_surf(i) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
    handles_norm(i) = quiver3(0, 0, 0, table_copy(i, 1), table_copy(i, 2), table_copy(i, 3), 'Color', 'red');
elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))    
    %points_set_j_con.points = points_set_i_con.points;
    
    %table_copy(j, :) = table_copy(i, :);
    %points_set_nonred{j} = points_set_nonred{i};
    %points_set_j_con.tipo =points_set_i_con.tipo; 
    points_set_j_con.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_j_con.points;
    set(handles_surf(j), 'Visible', 'off');
    delete(handles_surf(j));
    set(handles_norm(j), 'Visible', 'off');
    delete(handles_norm(j));
    handles_surf(j) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
    handles_norm(j) = quiver3(0, 0, 0, table_copy(j, 1), table_copy(j, 2), table_copy(j, 3), 'Color', 'red');
else    
    N = cross(points_set_i_con.points(1, :)-points_set_i_con.points(2, :), points_set_i_con.points(1, :)-points_set_i_con.points(3, :));
    if any(N ~= 0) 
        C = dot(N, [0 0 1]);
    else
        C = 1;
    end
    N2 = cross(points_set_j_con.points(1, :)-points_set_j_con.points(2, :), points_set_j_con.points(1, :)-points_set_j_con.points(3, :));    
    if any(N2 ~= 0) 
        C2 = dot(N2, [0 0 1]);
    else
        C2 = 1;
    end

    if (abs(C) < e) && (abs(C2) < e)
        if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
                all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))
            points_set_i_con.points = points_set_j_con.points;
            
            table_copy(i, :) = table_copy(j, :);
            points_set_nonred{i} = points_set_nonred{j};
            points_set_i_con.tipo =points_set_j_con.tipo;
        elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
                all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))
            points_set_j_con.points = points_set_i_con.points;
            
            table_copy(j, :) = table_copy(i, :);
            points_set_nonred{j} = points_set_nonred{i};
            points_set_j_con.tipo =points_set_i_con.tipo;
        end
    elseif ((abs(C) < e) || (abs(C2) < e)) && x == 0 && y ~= 0
            points_set_i_con.points(2, :) = [x y 0];
            points_set_i_con.points(4, :) = points_set_i_con.points(3, :);

            points_set_j_con.points(2, :) = [x y 0];
            points_set_j_con.points(4, :) = [x y 0];                                            
    elseif ((abs(C) < e) || (abs(C2) < e)) && y == 0 && x ~= 0
            points_set_i_con.points(3, :) = [x y 0];                                                
            points_set_i_con.points(4, :) = [x y 0];

            points_set_j_con.points(3, :) = [x y 0];                                                
            points_set_j_con.points(4, :) = points_set_i_con.points(3, :);
    elseif (caso == 4 || caso == 5 || caso == 6) && x == 0 && y == 0 %&& (points_set_nonred{i}.tipo == 10 || points_set_nonred{j}.tipo == 10)           
            syms x4 y4 z4 z;
            sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                        [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
            if ~isempty(sol)
                sol = [sol.x4 sol.y4 sol.z4];
                if ~isempty(symvar(sol))
                    z_sol = solve(sol(3)-max([points_set_i_con.points(:,3);points_set_j_con.points(:,3)]),z); %#ok<NODEF>
                    z = z_sol; p_new = eval(eval(sol)); %#ok<NASGU>
                    if all(p_new(:) >= 0)
                        if m_i < m_j
                            if (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 9 || points_set_i_con.tipo == 26) && ...
                                (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 9 || points_set_j_con.tipo == 21)
                                if (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 9) && points_set_i_con.tipo == 26  %%Con eje X                                         
                                    if (caso == 4 && d_j(2) >=0 && d_i(2) <= 0)
                                        points_set_j_con.points(1, :) = eval(sol);
                                        points_set_i_con.points(4, :) = eval(sol); %%3    
                                    elseif (caso == 6 && d_j(2) <=0 && d_i(2) <= 0)
                                        points_set_j_con.points(1, :) = eval(sol);
                                        points_set_i_con.points(2, :) = eval(sol); %%3
                                    end
                                elseif (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 9) && points_set_j_con.tipo == 21 %%Con eje Y
                                    if (caso == 4 && d_i(2) <=0 && d_j(2) >= 0)
                                        points_set_i_con.points(2, :) = eval(sol); 
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(2) <=0 && d_j(2) <= 0)
                                        points_set_i_con.points(1, :) = eval(sol); 
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    end
                                end
                            end
                        else
                            if (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 9 || points_set_j_con.tipo == 26) && ...
                                (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 9 || points_set_i_con.tipo == 21) || ...
                                (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 9) && (points_set_j_con.tipo == 21)
                                if (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 9) && points_set_j_con.tipo == 26 %%Con eje X
                                    if (caso == 4 && d_i(2) >=0 && d_j(2) <= 0)
                                        points_set_i_con.points(1, :) = eval(sol); 
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(2) <=0 && d_j(2) <= 0)
                                        points_set_i_con.points(1, :) = eval(sol); 
                                        points_set_j_con.points(2, :) = eval(sol); %%3
                                    end
                                elseif (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 9) && points_set_i_con.tipo == 21 %%Con eje Y
                                    if (caso == 4 && d_j(2) <=0 && d_i(2) >= 0)
                                        points_set_j_con.points(2, :) = eval(sol); 
                                        points_set_i_con.points(4, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_j(2) <=0 && d_i(2) <= 0)
                                        points_set_j_con.points(1, :) = eval(sol); 
                                        points_set_i_con.points(4, :) = eval(sol); %%3
                                    end
                                elseif (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 9) && (points_set_j_con.tipo == 21)
                                    if (caso == 4 && d_i(2) <=0 && d_j(2) >= 0)
                                        points_set_i_con.points(2, :) = eval(sol); 
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(2) <=0 && d_j(2) <= 0)
                                        points_set_i_con.points(1, :) = eval(sol); 
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    end
                                end
                            end
                        end
                    end
                end
            end
    else
        IsPlaneCrossIntercept = 1;
    end
end

function [points_set_i_con,points_set_j_con, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptXZ(points_set_i_con, points_set_j_con, x, z, table_copy, ...
    i, j, caso)
global points_set_nonred m_i m_j handles_surf handles_norm;
e = 0.00005; %%%% Cambiar si se cambia precisión
IsPlaneCrossIntercept = 0;

d_i = table_copy(i, 1:3); % vector director de i
d_j = table_copy(j, 1:3); % vector director de j
%m_i = (points_set_i_con.points(1,3)-points_set_i_con.points(3,3))/(points_set_i_con.points(1,1)-points_set_i_con.points(3,1));
%m_j = (points_set_j_con.points(1,3)-points_set_j_con.points(3,3))/(points_set_j_con.points(1,1)-points_set_j_con.points(3,1));
if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))    
    %points_set_i_con.points = points_set_j_con.points;    
    %table_copy(i, :) = table_copy(j, :);
    %points_set_nonred{i} = points_set_nonred{j};
    %points_set_i_con.tipo = points_set_j_con.tipo;
    
    points_set_i_con.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p =  points_set_i_con.points;
    set(handles_surf(i), 'Visible', 'off');
    delete(handles_surf(i));
    set(handles_norm(i), 'Visible', 'off');
    delete(handles_norm(i));
    handles_surf(i) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
    handles_norm(i) = quiver3(0, 0, 0, table_copy(i, 1), table_copy(i, 2), table_copy(i, 3), 'Color', 'red');
elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))    
    %points_set_j_con.points = points_set_i_con.points;
    
    %table_copy(j, :) = table_copy(i, :);
    %points_set_nonred{j} = points_set_nonred{i};
    %points_set_j_con.tipo = points_set_i_con.tipo;
    
    points_set_j_con.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p =  points_set_j_con.points;
    set(handles_surf(j), 'Visible', 'off');
    delete(handles_surf(j));
    set(handles_norm(j), 'Visible', 'off');
    delete(handles_norm(j));
    handles_surf(j) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
    handles_norm(j) = quiver3(0, 0, 0, table_copy(j, 1), table_copy(j, 2), table_copy(j, 3), 'Color', 'red');
else
    N = cross(points_set_i_con.points(1, :)-points_set_i_con.points(2, :), points_set_i_con.points(1, :)-points_set_i_con.points(3, :));
    if any(N ~= 0) 
        C = dot(N, [0 1 0]);
    else
        C = 1;
    end
    N2 = cross(points_set_j_con.points(1, :)-points_set_j_con.points(2, :), points_set_j_con.points(1, :)-points_set_j_con.points(3, :));
    if any(N2 ~= 0) 
        C2 = dot(N2, [0 1 0]);
    else
        C2 = 1;
    end

    if (abs(C) < e) && (abs(C2) < e)
        if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
                all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))
            points_set_i_con.points = points_set_j_con.points;
             
            table_copy(i, :) = table_copy(j, :);
            points_set_nonred{i} = points_set_nonred{j};
            points_set_i_con.tipo = points_set_j_con.tipo;
        elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
                all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))
            points_set_j_con.points = points_set_i_con.points;
            
            table_copy(j, :) = table_copy(i, :);
            points_set_nonred{j} = points_set_nonred{i};
            points_set_j_con.tipo = points_set_i_con.tipo;
        end 
    elseif ((abs(C) < e) || (abs(C2) < e)) && x == 0 && z ~= 0
            points_set_i_con.points(1, :) = [x 0 z];
            points_set_i_con.points(4, :) = points_set_i_con.points(3, :);


            points_set_j_con.points(1, :) = [x 0 z];                                                    
            points_set_j_con.points(4, :) = [x 0 z];
    elseif ((abs(C) < e) || (abs(C2) < e)) && z == 0 && x ~= 0
            points_set_i_con.points(3, :) = [x 0 z];                                                
            points_set_i_con.points(4, :) = [x 0 z];                                                

            points_set_j_con.points(3, :) = [x 0 z];                                                
            points_set_j_con.points(4, :) = points_set_j_con.points(3, :);
    elseif (caso == 4 || caso == 5 || caso == 6) && x == 0 && z == 0 %&& (points_set_nonred{i}.tipo == 10 || points_set_nonred{j}.tipo == 10)        
            syms x4 y4 z4 z;
            sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                        [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
            if ~isempty(sol)
                sol = [sol.x4 sol.y4 sol.z4];
                if ~isempty(symvar(sol))
                    z_sol = solve(sol(2)-max([points_set_i_con.points(:,2);points_set_j_con.points(:,2)]),z);
                    z = z_sol; p_new = eval(eval(sol)); %#ok<NASGU>
                    if all(p_new(:) >= 0)
                        if m_i <= m_j
                            if (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 10 || points_set_i_con.tipo == 26) && ...
                                (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 10 || points_set_j_con.tipo == 16)
                                if (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 10) && points_set_j_con.tipo == 16  %%Con eje Z                                       
                                    if (caso == 4 && d_j(3) <=0 && d_i(3) >= 0)
                                        points_set_i_con.points(3, :) = eval(sol);
                                        points_set_j_con.points(1, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_j(3) <=0 && d_i(3) <= 0)
                                        points_set_i_con.points(3, :) = eval(sol);
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    end
                                elseif (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 10) && points_set_i_con.tipo == 26 %%Con eje X
                                    if (caso == 4 && d_i(3) <=0 && d_j(3) >= 0)
                                        points_set_j_con.points(2, :) = eval(sol); 
                                        points_set_i_con.points(2, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(3) <=0 && d_j(3) <= 0)
                                        points_set_j_con.points(2, :) = eval(sol); 
                                        points_set_i_con.points(1, :) = eval(sol); %%3
                                    end
                                end
                            end
                        else
                            if (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 10 || points_set_j_con.tipo == 26) && ...
                                (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 10 || points_set_i_con.tipo == 16) 
                                if (points_set_j_con.tipo == 8 || points_set_j_con.tipo == 10) && points_set_i_con.tipo == 16 %%Con eje Z
                                    if (caso == 4 && d_i(3) <=0 && d_j(3) >= 0)
                                        points_set_j_con.points(3, :) = eval(sol); 
                                        points_set_i_con.points(1, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(3) <=0 && d_j(3) <= 0)
                                        points_set_j_con.points(3, :) = eval(sol); 
                                        points_set_i_con.points(4, :) = eval(sol); %%3
                                    end
                                elseif (points_set_i_con.tipo == 8 || points_set_i_con.tipo == 10) && points_set_j_con.tipo == 26 %%Con eje X
                                    if (caso == 4 && d_i(3) >=0 && d_j(3) <= 0)
                                        points_set_i_con.points(2, :) = eval(sol); 
                                        points_set_j_con.points(2, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(3) <=0 && d_j(3) <= 0)
                                        points_set_i_con.points(2, :) = eval(sol); 
                                        points_set_j_con.points(1, :) = eval(sol); %%3
                                    end
                                end
                            end
                        end
                    end
                end
            end
    else
        IsPlaneCrossIntercept = 1;
    end
end

function [points_set_i_con,points_set_j_con, table_copy, IsPlaneCrossIntercept] = NonPositiveInterceptYZ(points_set_i_con, points_set_j_con, y, z, table_copy, ...
    i, j, caso)
global points_set_nonred m_i m_j handles_surf handles_norm;
e = 0.00005; %%%% Cambiar si se cambia precisión
IsPlaneCrossIntercept = 0;

d_i = table_copy(i, 1:3); % vector director de i
d_j = table_copy(j, 1:3); % vector director de j
%m_i = (points_set_i_con.points(1,3)-points_set_i_con.points(2,3))/(points_set_i_con.points(1,2)-points_set_i_con.points(2,2));
%m_j = (points_set_j_con.points(1,3)-points_set_j_con.points(2,3))/(points_set_j_con.points(1,2)-points_set_j_con.points(2,2));
if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))    
    %points_set_i_con.points = points_set_j_con.points;
    
    %table_copy(i, :) = table_copy(j, :);
    %points_set_nonred{i} = points_set_nonred{j};
    %points_set_i_con.tipo = points_set_j_con.tipo;
    points_set_i_con.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_i_con.points;
    set(handles_surf(i), 'Visible', 'off');
    delete(handles_surf(i));
    set(handles_norm(i), 'Visible', 'off');
    delete(handles_norm(i));
    handles_surf(i) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
    handles_norm(i) = quiver3(0, 0, 0, table_copy(i, 1), table_copy(i, 2), table_copy(i, 3), 'Color', 'red');
elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))    
    %points_set_j_con.points = points_set_i_con.points;
    
    %table_copy(j, :) = table_copy(i, :);
    %points_set_nonred{j} = points_set_nonred{i};
    %points_set_j_con.tipo = points_set_i_con.tipo;
    points_set_j_con.points = [0 0 0;0 0 0;0 0 0;0 0 0]; p = points_set_j_con.points;
    set(handles_surf(j), 'Visible', 'off');
    delete(handles_surf(j));
    set(handles_norm(j), 'Visible', 'off');
    delete(handles_norm(j));
    handles_surf(j) = patch('xdata',p(:,1), 'ydata',p(:,2), 'zdata', p(:,3), 'FaceColor', 'g'); hold on;
    handles_norm(j) = quiver3(0, 0, 0, table_copy(j, 1), table_copy(j, 2), table_copy(j, 3), 'Color', 'red');
else
    N = cross(points_set_i_con.points(1, :)-points_set_i_con.points(2, :), points_set_i_con.points(1, :)-points_set_i_con.points(3, :));    
    if any(N ~= 0) 
        C = dot(N, [1 0 0]);
    else
        C = 1;
    end
    N2 = cross(points_set_j_con.points(1, :)-points_set_j_con.points(2, :), points_set_j_con.points(1, :)-points_set_j_con.points(3, :));    
    if any(N2 ~= 0) 
        C2 = dot(N2, [1 0 0]);
    else
        C2 = 1;
    end

    if (abs(C) < e) && (abs(C2) < e)
        if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
                all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))
            points_set_i_con.points = points_set_j_con.points;
            
            table_copy(i, :) = table_copy(j, :);
            points_set_nonred{i} = points_set_nonred{j};
            points_set_i_con.tipo = points_set_j_con.tipo;
        elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
                all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))
            points_set_j_con.points = points_set_i_con.points;
            
            table_copy(j, :) = table_copy(i, :);
            points_set_nonred{j} = points_set_nonred{i};
            points_set_j_con.tipo = points_set_i_con.tipo;
        end
    elseif ((abs(C) < e) ||  (abs(C2) < e)) && y == 0 && z ~= 0
            points_set_i_con.points(1, :) = [0 y z];
            points_set_i_con.points(4, :) = points_set_i_con.points(3, :);

            points_set_j_con.points(1, :) = [0 y z];                                                
            points_set_j_con.points(4, :) = points_set_j_con.points(3, :);
    elseif ((abs(C) < e) ||  (abs(C2) < e)) && z == 0 && y ~= 0
            points_set_i_con.points(2, :) = [0 y z];                                                
            points_set_i_con.points(4, :) = points_set_i_con.points(3, :);

            points_set_j_con.points(2, :) = [0 y z];                                                
            points_set_j_con.points(4, :) = points_set_j_con.points(3, :);
    elseif (caso == 4 || caso == 5 || caso == 6) && y == 0 && z == 0 %&& (points_set_nonred{i}.tipo == 10 || points_set_nonred{j}.tipo == 10)                 
            syms x4 y4 z4 z;
            sol = solve([num2str(table_copy(i, 1)), '*x4+', num2str(table_copy(i, 2)), '*y4+', num2str(table_copy(i, 3)), '*z4 =', num2str(table_copy(i, 4))], ...
                        [num2str(table_copy(j, 1)), '*x4+', num2str(table_copy(j, 2)), '*y4+', num2str(table_copy(j, 3)), '*z4 =', num2str(table_copy(j, 4))], x4, y4, z4);
            if ~isempty(sol)
                sol = [sol.x4 sol.y4 sol.z4];
                if ~isempty(symvar(sol))
                    z_sol = solve(sol(1)-max([points_set_i_con.points(:,1);points_set_j_con.points(:,1)]),z);
                    z = z_sol; p_new = eval(eval(sol)); %#ok<NASGU>
                    if all(p_new(:) >= 0)
                        if m_i <= m_j
                            if (points_set_i_con.tipo == 9 || points_set_i_con.tipo == 10 || points_set_i_con.tipo == 21) && ...
                                (points_set_j_con.tipo == 9 || points_set_j_con.tipo == 10 || points_set_j_con.tipo == 16)                                                        
                                if (points_set_i_con.tipo == 9 || points_set_i_con.tipo == 10) && points_set_j_con.tipo == 16  %%Con eje Z                                         
                                    if (caso == 4 && d_j(3) <=0 && d_i(3) >= 0)
                                        points_set_i_con.points(2, :) = eval(sol);
                                        points_set_j_con.points(1, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_j(3) <=0 && d_i(3) <= 0)
                                        points_set_i_con.points(2, :) = eval(sol);
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    end
                                elseif (points_set_j_con.tipo == 9 || points_set_j_con.tipo == 10) && points_set_i_con.tipo == 21  %%Con eje Y
                                    if (caso == 4 && d_i(3) >=0 && d_j(3) <= 0)
                                        points_set_j_con.points(1, :) = eval(sol);
                                        points_set_i_con.points(1, :) = eval(sol); %%3
                                    elseif (caso == 5 && d_i(3) >=0 && d_j(3) >= 0)
                                        points_set_j_con.points(1, :) = eval(sol);
                                        points_set_i_con.points(4, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(3) <=0 && d_j(3) <= 0)
                                        points_set_j_con.points(3, :) = eval(sol);
                                        points_set_i_con.points(1, :) = eval(sol); %%3
                                    end
                                end
                            end
                        else
                            if (points_set_j_con.tipo == 9 || points_set_j_con.tipo == 10 || points_set_j_con.tipo == 21) && ...
                                (points_set_i_con.tipo == 9 || points_set_i_con.tipo == 10 || points_set_i_con.tipo == 16)   
                                if (points_set_j_con.tipo == 9 || points_set_j_con.tipo == 10) && points_set_i_con.tipo == 16 %%Con eje Z                                          
                                    if (caso == 4 && d_i(3) <=0 && d_j(3) >= 0)
                                        points_set_j_con.points(2, :) = eval(sol);
                                        points_set_i_con.points(1, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_i(3) <=0 && d_j(3) <= 0)
                                        points_set_j_con.points(2, :) = eval(sol);
                                        points_set_i_con.points(4, :) = eval(sol); %%3
                                    end
                                elseif (points_set_i_con.tipo == 9 || points_set_i_con.tipo == 10) && points_set_j_con.tipo == 21  %%Con eje Y                                         
                                    if (caso == 4 && d_j(3) >=0 && d_i(3) <= 0)
                                        points_set_i_con.points(1, :) = eval(sol);
                                        points_set_j_con.points(1, :) = eval(sol); %%3
                                    elseif (caso == 5 && d_j(3) >=0 && d_i(3) >= 0)
                                        points_set_i_con.points(1, :) = eval(sol);
                                        points_set_j_con.points(4, :) = eval(sol); %%3
                                    elseif (caso == 6 && d_j(3) <=0 && d_i(3) <= 0)
                                        points_set_i_con.points(3, :) = eval(sol);
                                        points_set_j_con.points(1, :) = eval(sol); %%3
                                    end
                                end
                            end
                        end
                    end
                end
            end  
    else
        IsPlaneCrossIntercept = 1;
    end
end


function [points_set_i_con,points_set_j_con, table_copy, pasar] =  LineInterceptXY(points_set_i_con,points_set_j_con, table_copy, ...
     i, j, pasar)
global points_set_nonred;
e = 0.00005; %%%% Cambiar si se cambia precisión

if points_set_nonred{i}.points(3, 1) >= points_set_nonred{j}.points(3, 1) && ...
        points_set_nonred{i}.points(2, 2) >= points_set_nonred{j}.points(2, 2)
    if points_set_i_con.tipo == 36     %%% Plano XY
       ind = SearchPoints(points_set_j_con.points, 1);
       points = unique(points_set_j_con.points(ind,:), 'rows');
       dim = size(points);
       if dim(1) == 2
           points_set_i_con.points(1, :) = points(1,:);
           points_set_i_con.points(2, :) = points(2,:);
           points_set_i_con.points(3, :) = [0 0 0];
           points_set_i_con.points(4, :) = [0 0 0];
           
           N = cross(points_set_j_con.points(1, :)-points_set_j_con.points(2, :), points_set_j_con.points(1, :)-points_set_j_con.points(3, :));
           if N(3)~= 0
                C = dot(N, [0 0 1]);
           else
               C = 1;
           end
           if N(2)~= 0
               C2 = dot(N, [0 1 0]);
           else
               C2 = 1;
           end
           if N(1)~= 0
                C3 = dot(N, [1 0 0]);
           else
               C3 = 1;
           end
           if abs(C) < e && any(points(:,3) ~= 0)
               points_set_i_con.points(4, :) = [0 0 max([points(1,3) points(2,3)])];
           elseif abs(C2) < e && any(points(:,2) ~= 0)
               points_set_i_con.points(4, :) = [0 max([points(1,2) points(2,2)]) 0];
           elseif abs(C3) < e && any(points(:,1) ~= 0)
               points_set_i_con.points(4, :) = [max([points(1,1) points(2,1)]) 0 0];
           end
       else
           pasar = 0;
       end       
    elseif points_set_i_con.tipo == 13
       points_set_i_con.points(3, :) = points_set_j_con.points(3, :);
       if all(points_set_j_con.points(4, :) == points_set_j_con.points(3, :))
            points_set_i_con.points(4, :) = points_set_j_con.points(3, :);
       else
           points_set_i_con.points(4, :) = points_set_j_con.points(4, :);
       end
       points_set_i_con.points(2, :) = points_set_j_con.points(2, :);
       
       if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
               all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))
            points_set_i_con.points = points_set_j_con.points;
            
            table_copy(i, :) = table_copy(j, :);
            points_set_nonred{i} = points_set_nonred{j};
            points_set_i_con.tipo = points_set_j_con.tipo;
       end
    end
elseif points_set_nonred{i}.points(3, 1) < points_set_nonred{j}.points(3, 1) && ...
        points_set_nonred{i}.points(2, 2) < points_set_nonred{j}.points(2, 2)
    if points_set_j_con.tipo == 36     %%% Plano XY
       ind = SearchPoints(points_set_i_con.points, 1);
       points = unique(points_set_i_con.points(ind,:), 'rows');
       dim = size(points);
       if dim(1) == 2
           points_set_j_con.points(1, :) = points(1,:);
           points_set_j_con.points(2, :) = points(2,:);
           points_set_j_con.points(3, :) = [0 0 0];
           points_set_j_con.points(4, :) = [0 0 0];
       else
           pasar = 0;
       end       
    elseif points_set_j_con.tipo == 13
       points_set_j_con.points(3, :) = points_set_i_con.points(3, :);
       if all(points_set_i_con.points(4, :) == points_set_i_con.points(3, :))
            points_set_j_con.points(4, :) = points_set_i_con.points(3, :);
       else
           points_set_j_con.points(4, :) = points_set_i_con.points(4, :);
       end
       points_set_j_con.points(2, :) = points_set_i_con.points(2, :);
       
       if all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
               all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))
            points_set_j_con.points = points_set_i_con.points;
            
            table_copy(j, :) = table_copy(i, :);
            points_set_nonred{j} = points_set_nonred{i};
            points_set_j_con.tipo = points_set_i_con.tipo;
       end
    end
end

function [points_set_i_con,points_set_j_con, table_copy, pasar] =  LineInterceptXZ(points_set_i_con,points_set_j_con, table_copy, ...
     i, j, pasar)
global points_set_nonred;
e = 0.00005; %%%% Cambiar si se cambia precisión

if points_set_nonred{i}.points(1, 3) < points_set_nonred{j}.points(1, 3) && ...
            points_set_nonred{i}.points(3, 1) < points_set_nonred{j}.points(3, 1)
  if points_set_j_con.tipo == 33
       ind = SearchPoints(points_set_i_con.points, 2);
       points = unique(points_set_i_con.points(ind,:), 'rows');
       dim = size(points);
       if dim(1) == 2
           points_set_j_con.points(1, :) = points(1,:);
           points_set_j_con.points(2, :) = points(2,:);
           points_set_j_con.points(3, :) = [0 0 0];
           points_set_j_con.points(4, :) = [0 0 0];
       else
           pasar = 0;
       end
  elseif points_set_j_con.tipo == 18
      points_set_j_con.points(3, :) = points_set_i_con.points(3, :);
      if all(points_set_i_con.points(4, :) == points_set_i_con.points(3, :))
        points_set_j_con.points(4, :) = points_set_i_con.points(3, :);
      else
          points_set_j_con.points(4, :) = points_set_i_con.points(4, :);
      end
      
      if all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
          all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))
        points_set_j_con.points = points_set_i_con.points;
                
        table_copy(j, :) = table_copy(i, :);
        points_set_nonred{j} = points_set_nonred{i};
        points_set_j_con.tipo = points_set_i_con.tipo;
      end
      points_set_j_con.points(1, :) = points_set_i_con.points(1, :);
  end
elseif points_set_nonred{i}.points(1, 3) >= points_set_nonred{j}.points(1, 3) && ...
            points_set_nonred{i}.points(3, 1) >= points_set_nonred{j}.points(3, 1)
  if points_set_i_con.tipo == 33
       ind = SearchPoints(points_set_j_con.points, 2);
       points = unique(points_set_j_con.points(ind,:), 'rows');
       dim = size(points);
       if dim(1) == 2
           points_set_i_con.points(1, :) = points(1,:);
           points_set_i_con.points(2, :) = points(2,:);
           points_set_i_con.points(3, :) = [0 0 0];
           points_set_i_con.points(4, :) = [0 0 0];
           
           N = cross(points_set_j_con.points(1, :)-points_set_j_con.points(2, :), points_set_j_con.points(1, :)-points_set_j_con.points(3, :));
           if N(3)~= 0
                C = dot(N, [0 0 1]);
           else
               C = 1;
           end
           if N(2)~= 0
               C2 = dot(N, [0 1 0]);
           else
               C2 = 1;
           end
           if N(1)~= 0
                C3 = dot(N, [1 0 0]);
           else
               C3 = 1;
           end
               
           if abs(C) < e && any(points(:,3) ~= 0)
               points_set_i_con.points(4, :) = [0 0 max([points(1,3) points(2,3)])];
           elseif abs(C2) < e && any(points(:,2) ~= 0)
               points_set_i_con.points(4, :) = [0 max([points(1,2) points(2,2)]) 0];
           elseif abs(C3) < e && any(points(:,1) ~= 0)
               points_set_i_con.points(4, :) = [max([points(1,1) points(2,1)]) 0 0];
           end
       else
           pasar = 0;
       end      
  elseif points_set_i_con.tipo == 18
      points_set_i_con.points(3, :) = points_set_j_con.points(3, :);
      if all(points_set_j_con.points(4, :) == points_set_j_con.points(3, :))
        points_set_i_con.points(4, :) = points_set_j_con.points(3, :);
      else
          points_set_i_con.points(4, :) = points_set_j_con.points(4, :);
      end
      
      if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
              all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))
        points_set_i_con.points = points_set_j_con.points;
        
        table_copy(i, :) = table_copy(j, :);
        points_set_nonred{i} = points_set_nonred{j};
        points_set_i_con.tipo = points_set_j_con.tipo;
      end
      points_set_i_con.points(1, :) = points_set_j_con.points(1, :);
  end
end
 
function [points_set_i_con,points_set_j_con, table_copy, pasar] =  LineInterceptYZ(points_set_i_con,points_set_j_con, table_copy, ...
     i, j, pasar)
global points_set_nonred;
e = 0.00005; %%%% Cambiar si se cambia precisión

if points_set_nonred{i}.points(2, 2) >= points_set_nonred{j}.points(2, 2) && ...
            points_set_nonred{i}.points(1, 3) >= points_set_nonred{j}.points(1, 3)
   if points_set_i_con.tipo == 30
       ind = SearchPoints(points_set_j_con.points, 3);
       points = unique(points_set_j_con.points(ind,:), 'rows');
       dim = size(points);
       if dim(1) == 2
           points_set_i_con.points(1, :) = points(1,:);
           points_set_i_con.points(2, :) = points(2,:);
           points_set_i_con.points(3, :) = [0 0 0];
           points_set_i_con.points(4, :) = [0 0 0];
           
           N = cross(points_set_j_con.points(1, :)-points_set_j_con.points(2, :), points_set_j_con.points(1, :)-points_set_j_con.points(3, :));
           if N(3)~= 0
                C = dot(N, [0 0 1]);
           else
               C = 1;
           end
           if N(2)~= 0
               C2 = dot(N, [0 1 0]);
           else
               C2 = 1;
           end
           if N(1)~= 0
                C3 = dot(N, [1 0 0]);
           else
               C3 = 1;
           end
           if abs(C) < e
               points_set_i_con.points(4, :) = [0 0 max([points(1,3) points(2,3)])];
           elseif abs(C2) < e
               points_set_i_con.points(4, :) = [0 max([points(1,2) points(2,2)]) 0];
           elseif abs(C3) < e
               points_set_i_con.points(4, :) = [max([points(1,1) points(2,1)]) 0 0];
           end
       else
           pasar = 0;
       end                        
   elseif points_set_i_con.tipo == 23
       points_set_i_con.points(2, :) = points_set_j_con.points(2, :);
       points_set_i_con.points(1, :) = points_set_j_con.points(1, :);
       
       if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
               all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))
           points_set_i_con.points = points_set_j_con.points;
           
           table_copy(i, :) = table_copy(j, :);
           points_set_nonred{i} = points_set_nonred{j};
           points_set_i_con.tipo = points_set_j_con.tipo;
       end
   end
elseif points_set_nonred{i}.points(2, 2) < points_set_nonred{j}.points(2, 2) && ...
            points_set_nonred{i}.points(1, 3) < points_set_nonred{j}.points(1, 3)
   if points_set_j_con.tipo == 30
       ind = SearchPoints(points_set_i_con.points, 3);
       points = unique(points_set_i_con.points(ind,:), 'rows');
       dim = size(points);
       if dim(1) == 2
           points_set_j_con.points(1, :) = points(1,:);
           points_set_j_con.points(2, :) = points(2,:);           
           points_set_j_con.points(3, :) = [0 0 0];
           points_set_j_con.points(4, :) = [0 0 0];
       else
           pasar = 0;
       end                                    
   elseif points_set_j_con.tipo == 23
       points_set_j_con.points(2, :) = points_set_i_con.points(2, :);
       points_set_j_con.points(1, :) = points_set_i_con.points(1, :);
       
       if all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
               all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))
           points_set_j_con.points = points_set_i_con.points;
           
           table_copy(j, :) = table_copy(i, :);
           points_set_nonred{j} = points_set_nonred{i};
           points_set_j_con.tipo = points_set_i_con.tipo;
       end
   end
end

function ind =  SearchPoints(points_sets, plano)

if plano == 1
    ind = find(points_sets(:,3) == 0);
    points = unique(points_sets(ind,:), 'rows'); 
    dim = size(points);
    if dim(1) == 1
        ind2 = find(points_sets(:,3) ~= 0);
        [V, i] = min(points_sets(ind2,3)); %#ok<ASGLU>
        ind = [ind; ind2(i)];
    end   
elseif plano == 2
    ind = find(points_sets(:,2) == 0);
    points = unique(points_sets(ind,:), 'rows');
    dim = size(points);
    if dim(1) == 1
        ind2 = find(points_sets(:,2) ~= 0);
        [V, i] = min(points_sets(ind2,2)); %#ok<ASGLU>
        ind = [ind; ind2(i)];
    end
elseif plano == 3
    ind = find(points_sets(:,1) == 0);
    points = unique(points_sets(ind,:), 'rows');
    dim = size(points);
    if dim(1) == 1
        ind2 = find(points_sets(:,1) ~= 0);
        [V, i] = min(points_sets(ind2,1)); %#ok<ASGLU>
        ind = [ind; ind2(i)];
    end
end

function [points_set_i,points_set_j, table_copy] = PlaneCrossInterceptXY(handles, points_set_i, points_set_j, table_copy, i, j, ...
    index_1, index_2, index_3, index_4)
global points_set_nonred;
syms x1 y1 z1;
e = 0.00005; %%%% Cambiar si se cambia precisión

N2 = cross(points_set_j(1, :)-points_set_j(2, :), points_set_j(1, :)-points_set_j(3, :));
C2 = dot(N2, [0 0 1]);
N = cross(points_set_i(1, :)-points_set_i(2, :), points_set_i(1, :)-points_set_i(3, :));
C = dot(N, [0 0 1]);

if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))    
    points_set_i = points_set_j;
    
    table_copy(i, :) = table_copy(j, :);
    points_set_nonred{i} = points_set_nonred{j};
    %points_set_i_con.tipo = points_set_j_con.tipo;
elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))    
    points_set_j = points_set_i;
    
    table_copy(j, :) = table_copy(i, :);
    points_set_nonred{j} = points_set_nonred{i};
    %points_set_j_con.tipo = points_set_i_con.tipo;
elseif (abs(C) < e) && (abs(C2) < e)
    %%%% Desarrollar
    s = get(handles.listbox_operations, 'string');
    set(handles.listbox_operations, 'string', ...
         char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de planos paralelos a Z en desarrollo.'));
elseif (abs(C2) < e && points_set_nonred{i}.tipo == 1)
    %if points_set_i(index_2, index_coor) < points_set_j(index_2, index_coor)                                                            
    if (points_set_i(3, 3) > 0 || (points_set_i(3, 3) == 0 && points_set_i(4, 3) > 0)) ||...
            (points_set_i(2, 3) > 0 || (points_set_i(2, 3) == 0 && points_set_i(4, 3) > 0))
        if index_1 == 2
            if all(points_set_i(3, :) == points_set_i(4, :))                
                [points_set_i, points_set_j] = calcIntersect2(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);            
            else                                    
                [points_set_i, points_set_j] = calcIntersect(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            end
        elseif index_1 == 3
            if all(points_set_i(3, :) == points_set_i(4, :))
                [points_set_i, points_set_j] = calcIntersect2(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            else
                [points_set_i, points_set_j] = calcIntersect(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            end
        end                    
    elseif (points_set_j(2, 3) > 0 || (points_set_j(2, 3) == 0 && points_set_j(4, 3) > 0)) ||...
            (points_set_j(3, 3) > 0 || (points_set_j(3, 3) == 0 && points_set_j(4, 3) > 0))
        [points_set_i, points_set_j] = calcIntersect3(points_set_i, points_set_j, table_copy, i, j);        
    elseif points_set_i(2, 3) == 0 && points_set_i(3, 3) == 0
        [points_set_i, points_set_j] = calcIntersect3(points_set_i, points_set_j, table_copy, i, j);                
    end
elseif (abs(C) < e && points_set_nonred{j}.tipo == 1)   
    if (points_set_j(3, 3) > 0 || (points_set_j(3, 3) == 0 && points_set_j(4, 3) > 0)) ||...
            (points_set_j(2, 3) > 0 || (points_set_j(2, 3) == 0 && points_set_j(4, 3) > 0))
        if index_1 == 2           
            if all(points_set_j(3, :) == points_set_j(4, :))
                [points_set_j, points_set_i] = calcIntersect2(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            else
                [points_set_j, points_set_i] = calcIntersect(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            end
        elseif index_1 == 3
            if all(points_set_j(3, :) == points_set_j(4, :))                
                [points_set_j, points_set_i] = calcIntersect2(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            else
                [points_set_j, points_set_i] = calcIntersect(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            end
        end                    
    elseif (points_set_i(3, 3) > 0 || (points_set_i(3, 3) == 0 && points_set_i(4, 3) > 0)) ||...
            (points_set_i(2, 3) > 0 || (points_set_i(2, 3) == 0 && points_set_i(4, 3) > 0))
        [points_set_j, points_set_i] = calcIntersect3(points_set_j, points_set_i, table_copy, j, i);
                
    elseif points_set_j(2, 3) == 0 && points_set_j(3, 3) == 0       
        [points_set_j, points_set_i] = calcIntersect3(points_set_j, points_set_i, table_copy, j, i);
    end
elseif ((abs(C) < e && points_set_nonred{j}.tipo ~= 1) || (abs(C2) < e && points_set_nonred{i}.tipo ~= 1))
    %%%%% Desarrollar
    s = get(handles.listbox_operations, 'string');
    set(handles.listbox_operations, 'string', ...
            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 7 de plano paralelo a Z y de tipo distinto a 1 con limites interiores en desarrollo.'));
end


%%%%%%
function [points_set_i, points_set_j] = calcIntersect(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4)
global points_set_nonred;
%%%%%%%%

if points_set_i(4, 1) == 0
    index_copy = 2;
    frontier = points_set_i(4,:);
    dim = size(points_set_i);
    if dim(1) == 4
        frontier = [frontier; points_set_i(2,:); points_set_i(3,:)];
    else
        for k = 5:dim(1)
            frontier = [frontier; points_set_i(k,:)];
        end
        frontier = [frontier; points_set_i(2,:); points_set_i(3,:)];
    end
elseif points_set_i(4, 2) == 0
    index_copy = 3;
    frontier = points_set_i(4,:);
    dim = size(points_set_i);
    if dim(1) == 4
        frontier = [frontier; points_set_i(3,:); points_set_i(2,:)];
    else
        for k = 5:dim(1)
            frontier = [frontier; points_set_i(k,:)];
        end
        frontier = [frontier; points_set_i(3,:); points_set_i(2,:)];
    end
else
    index_copy = 4;
    frontier = points_set_i(2,:);
    dim = size(points_set_i);
    if dim(1) == 4
        frontier = [frontier; points_set_i(4,:); points_set_i(3,:)];
    else
        for k = 5:dim(1)
            frontier = [frontier; points_set_i(k,:)];
        end
        frontier = [frontier; points_set_i(4,:); points_set_i(3,:)];
    end
end

outer_points = frontier*table_copy(j, 1:3)' > table_copy(j, 4);
if ~all(outer_points(:)==0) && ~all(outer_points(:)==1)
    index_first = find(outer_points(:), 1, 'first');
    index_last = find(outer_points(:), 1, 'last');
    if index_first > 1
        sol = findIntercept(frontier(index_first-1, :), frontier(index_first, :), table_copy, i, j);
        if ~isempty(sol)
            sol = [sol.x1 sol.y1 sol.z1];
            if isempty(symvar(sol))                                                
                x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                if all([x y z] > 0)                                    
                    if index_first == dim(1)-1
                        frontier_aux = frontier;
                        frontier = frontier(1:index_first-1,:);
                        frontier(index_first,:) = [x y z];
                        frontier = [frontier; frontier_aux(index_first:dim(1)-1,:)];
                    else
                        frontier(index_first, :) = [x y z];
                    end                   
                    m_j = (points_set_nonred{j}.points(index_2, 2)-points_set_nonred{j}.points(index_1, 2))/...
                            (points_set_nonred{j}.points(index_2, 1)-points_set_nonred{j}.points(index_1, 1));
                    m_i = (frontier(index_first-1, 2)-frontier(index_first, 2))/...
                        (frontier(index_first-1, 1)-frontier(index_first, 1));
                    
                    if index_copy == 2 || index_copy == 4
                        if (m_i > 0) && m_j < 0
                            index_1 = 2; index_2 = 3;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = 3; index_2 = 2;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = 3; index_2 = 2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = 3; index_2 = 2;
                        end
                    elseif index_copy == 3
                        if (m_i > 0) && m_j < 0
                            index_1 = 3; index_2 = 2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = 2; index_2 = 3;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = 2; index_2 = 3;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = 2; index_2 = 3;
                        end
                    end
                    if index_1 == 2
                        if (points_set_j(index_2, 1) <= x && points_set_j(index_2, 2) >= y) && ...
                            (points_set_j(index_2, 1) <= x <= points_set_j(index_1, 1) && points_set_j(index_2, 2) <= y <= points_set_j(index_1, 2) || ...
                            points_set_j(index_1, 1) <= x <= points_set_j(index_2, 1) && points_set_j(index_1, 2) <= y <= points_set_j(index_2, 2))
                            if points_set_j(index_4, 1) == points_set_j(index_2, 1) && points_set_j(index_4, 2) == points_set_j(index_2, 2)
                                points_set_j(index_4, :) = [x y z];
                            else
                                points_set_j(index_3, :) = [x y z];
                            end
                            points_set_j(index_2, :) = [x y 0];
                        else
                            if points_set_j(index_4, 1) == points_set_j(index_2, 1) && points_set_j(index_4, 2) == points_set_j(index_2, 2)
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                            end 
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)                                                                                                 
                                        frontier(index_first, :) = [x y z];
                                    end
                                end
                            end
                        end
                        
                        if index_first >= dim(1)-2 
                            if points_set_j(index_4, 1) == points_set_j(index_2, 1) && points_set_j(index_4, 2) == points_set_j(index_2, 2)
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_3,:), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_4,:), table_copy, i, j);
                            end        
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] >= 0)
                                        if index_first == dim(1)-1
                                            %if frontier(dim(1), 1) >= x && frontier(dim(1), 2) >= y
                                                frontier(dim(1), :) = [x y z];
                                            %end
                                        else
                                            %if frontier(dim(1)-1, 1) >= x && frontier(dim(1)-1, 2) >= y
                                                frontier(dim(1)-1, :) = [x y z];
                                            %end
                                        end
                                    end
                                end
                            end
                        end
                    elseif index_1 == 3
                        if (points_set_j(index_2, 1) <= x && points_set_j(index_2, 2) >= y && ...                            
                            (points_set_j(index_2, 1) <= x <= points_set_j(index_1, 1) && points_set_j(index_2, 2) <= y <= points_set_j(index_1, 2) || ...
                            points_set_j(index_1, 1) <= x <= points_set_j(index_2, 1) && points_set_j(index_1, 2) <= y <= points_set_j(index_2, 2)))
                        %if (points_set_j(index_1, 1) <= x && points_set_j(index_1, 2) >= y)
                            if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
                                points_set_j(index_3, :) = [x y z];
                            else
                                points_set_j(index_4, :) = [x y z];
                            end
                            points_set_j(index_2, :) = [x y 0];
                            %minxyz(3) = min([minxyz(3) z]);
                        else
                            if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
                                %frontier(index_first, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) z_1];
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                            else
                                %frontier(index_first, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) z_1];
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                            end  
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)
                                        frontier(index_first,:) = [x y z];
                                    end
                                end
                            end
                        end                        

                        if index_first >= dim(1)-2 
                            if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_4,:), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_3,:), table_copy, i, j);
                            end
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] >= 0)
                                        if index_first == dim(1)-1
                                            %if frontier(dim(1), 1) >= x && frontier(dim(1), 2) >= y
                                                frontier(dim(1),:) = [x y z];
                                            %end
                                        else
                                            %if frontier(dim(1)-1, 1) >= x && frontier(dim(1)-1, 2) >= y
                                                frontier(dim(1)-1, :) = [x y z];
                                            %end
                                        end
                                    end
                                end
                            end
                        end
                    end                                        
                end
            end
        end
    end
    if index_last < dim(1)-1
        sol = findIntercept(frontier(index_last+1, :), frontier(index_last, :), table_copy, i, j);
        if ~isempty(sol)
            sol = [sol.x1 sol.y1 sol.z1];
            if isempty(symvar(sol))                                                
                x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                if all([x y z] > 0)                           
                    if index_last == index_first
                        frontier_aux = frontier;
                        frontier = frontier(1:index_last,:);
                        frontier(index_last+1,:) = [x y z];
                        frontier = [frontier; frontier_aux(index_last+1:dim(1)-1,:)];
                    else
                        frontier(index_last, :) = [x y z];
                    end
                    
                    m_j = (points_set_nonred{j}.points(index_2, 2)-points_set_nonred{j}.points(index_1, 2))/...
                            (points_set_nonred{j}.points(index_2, 1)-points_set_nonred{j}.points(index_1, 1));
                    m_i = (frontier(index_last+1, 2)-frontier(index_last, 2))/...
                        (frontier(index_last+1, 1)-frontier(index_last, 1));
                    
                    if index_copy == 2 || index_copy == 4
                        if (m_i > 0) && m_j < 0
                            index_1 = 2; index_2 = 3;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = 3; index_2 = 2;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = 3; index_2 = 2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = 3; index_2 = 2;
                        end
                    elseif index_copy == 3
                        if (m_i > 0) && m_j < 0
                            index_1 = 3; index_2 = 2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = 2; index_2 = 3;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = 2; index_2 = 3;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = 2; index_2 = 3;
                        end
                    end
                    if index_1 == 2
                        if (points_set_j(index_1, 1) <= x && points_set_j(index_1, 2) >= y && ...                            
                            (points_set_j(index_2, 1) <= x <= points_set_j(index_1, 1) && points_set_j(index_2, 2) <= y <= points_set_j(index_1, 2) || ...
                            points_set_j(index_1, 1) <= x <= points_set_j(index_2, 1) && points_set_j(index_1, 2) <= y <= points_set_j(index_2, 2)))
                        %if (points_set_j(index_2, 1) <= x && points_set_j(index_2, 2) >= y)
                            if points_set_j(index_4, 1) == points_set_j(index_2, 1) && points_set_j(index_4, 2) == points_set_j(index_2, 2)
                                points_set_j(index_3, :) = [x y z];
                            else
                                points_set_j(index_4, :) = [x y z];
                            end
                            points_set_j(index_1, :) = [x y 0];
                            %minxyz(3) = min([minxyz(3) z]);
                        else
                            if points_set_j(index_4, 1) == points_set_j(index_2, 1) && points_set_j(index_4, 2) == points_set_j(index_2, 2)
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_3, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_4, :), table_copy, i, j);
                            end                                        
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)
                                        if index_last == index_first
                                            frontier(index_last+1,:) = [x y z];
                                        else
                                            frontier(index_last, :) = [x y z];
                                        end
                                    end
                                end
                            end
                        end
                        if index_last == 1 
                            if points_set_j(index_4, 1) == points_set_j(index_2, 1) && points_set_j(index_4, 2) == points_set_j(index_2, 2)
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                            end
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] >= 0)
                                        %if frontier(index_last, 1) >= x && frontier(index_last, 2) >= y
                                            frontier(index_last, :) = [x y z];
                                        %end
                                    end
                                end
                            end
                        end
                    elseif index_1 == 3
                        
                        if ((points_set_j(index_1, 1) <= x && points_set_j(index_1, 2) >= y) &&...
                            (points_set_j(index_2, 1) <= x <= points_set_j(index_1, 1) && points_set_j(index_2, 2) <= y <= points_set_j(index_1, 2) || ...
                            points_set_j(index_1, 1) <= x <= points_set_j(index_2, 1) && points_set_j(index_1, 2) <= y <= points_set_j(index_2, 2)))
                        %if (points_set_j(index_2, 1) <= x && points_set_j(index_2, 2) >= y)
                            if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
                                points_set_j(index_4, :) = [x y z];
                            else
                                points_set_j(index_3, :) = [x y z];
                            end
                            points_set_j(index_1, :) = [x y 0];
                            %minxyz(3) = min([minxyz(3) z]);
                        else
                            if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_4, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_3, :), table_copy, i, j);
                            end                                        
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)
                                        %if index_last == index_first
                                            frontier(index_last+1,:) = [x y z];
                                        %else
                                            frontier(index_last, :) = [x y z];
                                        %end
                                    end
                                end
                            end
                        end
                    end
                    if index_last == 1 
                        if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
                            sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                        else
                            sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                        end
                        if ~(isempty(sol))
                            sol = [sol.x1 sol.y1 sol.z1];
                            if isempty(symvar(sol))
                                x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                if all([x y z] >= 0)
                                    %if frontier(index_last, 1) >= x && frontier(index_last, 2) >= y
                                        frontier(index_last, :) = [x y z];
                                    %end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if index_last - index_first > 1
        frontier = [frontier(1:index_first,:);frontier(index_last:dim(1)-1,:)];
    end
    dim = size(frontier);
    points_set_i = points_set_i(1, :);
    if index_copy == 3
        points_set_i = [points_set_i;frontier(dim(1),:);frontier(dim(1)-1,:);frontier(1,:)];
    elseif index_copy == 2
        points_set_i = [points_set_i;frontier(dim(1)-1,:);frontier(dim(1),:);frontier(1,:)];
    else
        points_set_i = [points_set_i;frontier(1,:);frontier(dim(1),:);frontier(dim(1)-1,:)];
    end
    for k = 2:dim(1)-2
        points_set_i(3+k, :) = frontier(k,:);
    end
elseif all(outer_points(:)==1)
    points_set_i = points_set_i(1:4, :);
    if index_copy == 3
        points_set_i(3,:) = points_set_i(4, :);
        points_set_i(4,:) = points_set_i(3, :);
    elseif index_copy == 2
        points_set_i(2,:) = points_set_i(4, :);
        points_set_i(4,:) = points_set_i(3, :);
    else                        
        points_set_i(4,:) = points_set_i(3, :);
    end
    [points_set_i, points_set_j] = calcIntersect3(points_set_i, points_set_j, table_copy, i, j);
    %%%%%
end



%%%%%%%%
function [points_set_i, points_set_j] = calcIntersect2(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4)
%%%%%%%
sol = findIntercept(points_set_i(index_1, :), points_set_i(index_2, :), table_copy, i, j);
if ~isempty(sol) 
    sol = [sol.x1 sol.y1 sol.z1];
    if isempty(symvar(sol))                                                
        x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
        if all([x y z] > 0)
            if index_1 == 2
                if points_set_i(index_1, 1) <= x && points_set_i(index_1, 2) >= y                                
                    points_set_i(4, :) = [x y z];                           
                end
            else
                if points_set_i(index_1, 1) >= x && points_set_i(index_1, 2) <= y                                
                    points_set_i(4, :) = [x y z];                           
                end
            end
            m_j = (points_set_j(index_1, 2)-points_set_j(index_2, 2))/...
                (points_set_j(index_1, 1)-points_set_j(index_2, 1));
            m_i = (points_set_i(index_1, 2)-points_set_i(index_2, 2))/...
                (points_set_i(index_1, 1)-points_set_i(index_2, 1));          
            
            outer_points = [points_set_i(index_1, :); points_set_i(index_2, :)]*table_copy(j, 1:3)' > table_copy(j, 4);
            index_first = find(outer_points(:)==0, 1, 'first');
            if (m_i < 0 && (m_j < 0 || abs(m_j) == Inf) && index_first == 1)
                index_first = 3; index_second = 2;
            elseif (m_i < 0 && (m_j < 0 || abs(m_j) == Inf) && index_first == 2)
                index_first = 2; index_second = 3;
            elseif (m_i < 0 && m_j > 0)
                index_first = 2; index_second = 3;            
            end
            
            if (points_set_j(index_first, 1) <= x && points_set_j(index_first, 2) >= y) || ...
                    (points_set_j(index_first, 1) >= x && points_set_j(index_first, 2) <= y) && ...
                    (points_set_j(index_2, 1) <= x <= points_set_j(index_1, 1) && points_set_j(index_2, 2) <= y <= points_set_j(index_1, 2) || ...
                            points_set_j(index_1, 1) <= x <= points_set_j(index_2, 1) && points_set_j(index_1, 2) <= y <= points_set_j(index_2, 2))
            %if points_set_j(index_1, 1) <= x && points_set_j(index_1, 2) >= y
                if points_set_j(index_3, 1) == points_set_j(index_second, 1) && points_set_j(index_3, 2) == points_set_j(index_second, 2)
                %if points_set_j(index_3, 1) == points_set_j(index_2, 1) && points_set_j(index_3, 2) == points_set_j(index_2, 2)
                    points_set_j(index_4, :) = [x y z];
                else
                    points_set_j(index_3, :) = [x y z];                  
                end
                points_set_j(index_first, :) = [x y 0];
                %minxyz(3) = min([minxyz(3) z]);                                
            else
                if points_set_j(index_3, 1) == points_set_j(index_second, 1) && points_set_j(index_3, 2) == points_set_j(index_second, 2)
                %if points_set_j(index_3, 1) == points_set_j(index_2, 1) && points_set_j(index_3, 2) == points_set_j(index_2, 2)
                    %points_set_i(4, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) z_1];
                    sol = findIntercept(points_set_j(index_first, :), points_set_j(index_4, :), table_copy, i, j);
                else
                    %points_set_i(4, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) z_1];
                    sol = findIntercept(points_set_j(index_first, :), points_set_j(index_3, :), table_copy, i, j);
                end                                
                if ~(isempty(sol))
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))
                        x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                        if all([x y z] > 0)
                            points_set_i(4, :) = [x y z];
                        end
                    end
                end
            end
            if points_set_j(index_3, 1) == points_set_j(index_second, 1) && points_set_j(index_3, 2) == points_set_j(index_second, 2)        
                sol = findIntercept(points_set_j(index_3, :), points_set_j(index_second, :), table_copy, i, j);
            else
                sol = findIntercept(points_set_j(index_4, :), points_set_j(index_second, :), table_copy, i, j);                    
            end
            if ~(isempty(sol))
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))
                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                    if all([x y z] > 0)
                        points_set_i(index_2, :) = [x y z];
                    end
                end
            end
        end        
    end        
end

%%%%%
function [points_set_i, points_set_j] = calcIntersect3(points_set_i, points_set_j, table_copy, i, j)
global points_set_nonred;

if points_set_j(2,1) > 0
    sol = findIntercept(points_set_j(2, :), points_set_j(1, :), table_copy, i, j);
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] > 0)
                points_set_i(4, :) = [x y z];
                points_set_j(1, :) = [x y z];
            end
        end
    end
    sol = findIntercept(points_set_nonred{j}.points(2,:), points_set_nonred{j}.points(1,:), table_copy, i, j);
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                points_set_i(2, :) = [x y z];                
            end
        end
    end
else
    sol = findIntercept(points_set_j(2, :), points_set_j(1, :), table_copy, i, j);
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                points_set_i(2, :) = [x y z];
                points_set_j(1, :) = [x y z];
            end
        end
    end
end

if points_set_j(3,2) > 0
    sol = findIntercept(points_set_j(3, :), points_set_j(4, :), table_copy, i, j);
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] > 0)
                if all(points_set_i(4, :) == points_set_i(3, :))
                    points_set_i(4, :) = [x y z];
                else
                    points_set_i(5, :) = [x y z];
                end
                points_set_j(4, :) = [x y z];
            end
        end
    end
    sol = findIntercept(points_set_nonred{j}.points(3,:), points_set_nonred{j}.points(4,:), table_copy, i, j);
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                if all(points_set_i(3, :) == points_set_i(4, :))
                    points_set_i(4, :) = [x y z];
                end
                points_set_i(3, :) = [x y z];                
            end
        end
    end
else
    sol = findIntercept(points_set_j(3, :), points_set_j(4, :), table_copy, i, j);
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                 if all(points_set_i(3, :) == points_set_i(4, :))
                    points_set_i(4, :) = [x y z];
                end
                points_set_i(3, :) = [x y z];
                points_set_j(4, :) = [x y z];
            end
        end
    end
end


%%%%%%%
function sol = findIntercept(p1, p2, table_copy, i, j)
syms x1 y1 z1;

d = p1 - p2;

sol = solve([num2str(table_copy(i, 1)), '*x1+', num2str(table_copy(i, 2)), '*y1+', num2str(table_copy(i, 3)), '*z1 =', num2str(table_copy(i, 4))], ...
    [num2str(table_copy(j, 1)), '*x1+', num2str(table_copy(j, 2)), '*y1+', num2str(table_copy(j, 3)), '*z1 =', num2str(table_copy(j, 4))], ...
    [num2str(d(2)), '*(x1-', num2str(p1(1)), ')-',num2str(d(1)), '*(y1-', num2str(p1(2)), ') = 0'], x1, y1, z1);
if ~isempty(sol)                     
    sol1 = [sol.x1 sol.y1 sol.z1];
else
    sol1 = x1;                 
end
if isempty(sol) || ~isempty(symvar(sol1)) 
    sol = solve([num2str(table_copy(i, 1)), '*x1+', num2str(table_copy(i, 2)), '*y1+', num2str(table_copy(i, 3)), '*z1 =', num2str(table_copy(i, 4))], ...
        [num2str(table_copy(j, 1)), '*x1+', num2str(table_copy(j, 2)), '*y1+', num2str(table_copy(j, 3)), '*z1 =', num2str(table_copy(j, 4))], ...
        [num2str(d(3)), '*(x1-', num2str(p1(1)), ')-',num2str(d(1)), '*(z1-', num2str(p1(3)), ') = 0'], x1, y1, z1);
end
if ~isempty(sol)                     
    sol1 = [sol.x1 sol.y1 sol.z1];
else
    sol1 = x1;                 
end
if isempty(sol) || ~isempty(symvar(sol1)) 
    sol = solve([num2str(table_copy(i, 1)), '*x1+', num2str(table_copy(i, 2)), '*y1+', num2str(table_copy(i, 3)), '*z1 =', num2str(table_copy(i, 4))], ...
        [num2str(table_copy(j, 1)), '*x1+', num2str(table_copy(j, 2)), '*y1+', num2str(table_copy(j, 3)), '*z1 =', num2str(table_copy(j, 4))], ...
        [num2str(d(3)), '*(y1-', num2str(p1(2)), ')-',num2str(d(2)), '*(z1-', num2str(p1(3)), ') = 0'], x1, y1, z1);
end


function [points_set_i,points_set_j, table_copy] = PlaneCrossInterceptXZ(handles, points_set_i, points_set_j, table_copy, i, j, ...
    index_1, index_2, index_3, index_4)
global points_set_nonred;
syms x1 y1 z1;
e = 0.00005; %%%% Cambiar si se cambia precisión

N2 = cross(points_set_j(1, :)-points_set_j(2, :), points_set_j(1, :)-points_set_j(3, :));
C2 = dot(N2, [0 1 0]);
N = cross(points_set_i(1, :)-points_set_i(2, :), points_set_i(1, :)-points_set_i(3, :));
C = dot(N, [0 1 0]);

if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))    
    points_set_i = points_set_j;
    
    table_copy(i, :) = table_copy(j, :);
    points_set_nonred{i} = points_set_nonred{j};
    %points_set_i_con.tipo = points_set_j_con.tipo;
elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))    
    points_set_j = points_set_i;
    
    table_copy(j, :) = table_copy(i, :);
    points_set_nonred{j} = points_set_nonred{i};
    %points_set_j_con.tipo = points_set_i_con.tipo;
elseif (abs(C) < e) && (abs(C2) < e)
    %%%% Desarrollar
    s = get(handles.listbox_operations, 'string');
    set(handles.listbox_operations, 'string', ...
         char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de planos paralelos a Z en desarrollo.'));
elseif (abs(C2) < e && points_set_nonred{i}.tipo == 1)
    %if points_set_i(index_2, index_coor) < points_set_j(index_2, index_coor)                                                            
    if (points_set_i(3, 2) > 0 || (points_set_i(3, 2) == 0 && points_set_i(4, 2) > 0)) ||...
            (points_set_i(1, 2) > 0 || (points_set_i(1, 2) == 0 && points_set_i(4, 2) > 0))
        if index_1 == 1
            if all(points_set_i(3, :) == points_set_i(4, :))                
                [points_set_i, points_set_j] = calcIntersect2_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);            
            else                                    
                [points_set_i, points_set_j] = calcIntersect_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            end
        elseif index_1 == 3
            if all(points_set_i(3, :) == points_set_i(4, :))
                [points_set_i, points_set_j] = calcIntersect2_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            else
                [points_set_i, points_set_j] = calcIntersect_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            end
        end                    
    elseif (points_set_j(1, 2) > 0 || (points_set_j(1, 2) == 0 && points_set_j(4, 2) > 0)) ||...
            (points_set_j(3, 2) > 0 || (points_set_j(3, 2) == 0 && points_set_j(4, 2) > 0))
        [points_set_i, points_set_j] = calcIntersect3_XZ_YZ(points_set_i, points_set_j, table_copy, i, j);        
    elseif points_set_i(1, 2) == 0 && points_set_i(3, 2) == 0
        [points_set_i, points_set_j] = calcIntersect3_XZ_YZ(points_set_i, points_set_j, table_copy, i, j);                
    end
elseif (abs(C) < e && points_set_nonred{j}.tipo == 1)   
    if (points_set_j(3, 2) > 0 || (points_set_j(3, 2) == 0 && points_set_j(4, 2) > 0)) ||...
            (points_set_j(1, 2) > 0 || (points_set_j(1, 2) == 0 && points_set_j(4, 2) > 0))
        if index_1 == 1           
            if all(points_set_j(3, :) == points_set_j(4, :))
                [points_set_j, points_set_i] = calcIntersect2_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            else
                [points_set_j, points_set_i] = calcIntersect_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            end
        elseif index_1 == 3
            if all(points_set_j(3, :) == points_set_j(4, :))                
                [points_set_j, points_set_i] = calcIntersect2_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            else
                [points_set_j, points_set_i] = calcIntersect_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            end
        end                    
    elseif (points_set_i(3, 2) > 0 || (points_set_i(3, 2) == 0 && points_set_i(4, 2) > 0)) ||...
            (points_set_i(1, 2) > 0 || (points_set_i(1, 2) == 0 && points_set_i(4, 2) > 0))
        [points_set_j, points_set_i] = calcIntersect3_XZ_YZ(points_set_j, points_set_i, table_copy, j, i);
                
    elseif points_set_j(1, 2) == 0 && points_set_j(3, 2) == 0       
        [points_set_j, points_set_i] = calcIntersect3_XZ__YZ(points_set_j, points_set_i, table_copy, j, i);
    end
elseif ((abs(C) < e && points_set_nonred{j}.tipo ~= 1) || (abs(C2) < e && points_set_nonred{i}.tipo ~= 1))
    %%%%% Desarrollar
    s = get(handles.listbox_operations, 'string');
    set(handles.listbox_operations, 'string', ...
            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 7 de plano paralelo a Z y de tipo distinto a 1 con limites interiores en desarrollo.'));
end


%%%%%%
function [points_set_i, points_set_j] = calcIntersect_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4)
global points_set_nonred;
%%%%%%%%

if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
    index_var1 = 1; index_var2 = 3;
    index_comp1 = 1; index_comp2 = 3;
elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
    index_var1 = 1; index_var2 = 2;
    index_comp1 = 2; index_comp2 = 3;
end
if points_set_i(index_3, index_comp1) == 0
    index_copy = index_var1;
    frontier = points_set_i(index_3,:);
    dim = size(points_set_i);
    if dim(1) == 4
        frontier = [frontier; points_set_i(index_var1,:); points_set_i(index_var2,:)];
    else
        for k = 5:dim(1)
            frontier = [frontier; points_set_i(k,:)];
        end
        frontier = [frontier; points_set_i(index_var1,:); points_set_i(index_var2,:)];
    end
elseif points_set_i(index_3, index_comp2) == 0
    index_copy = index_var2;
    frontier = points_set_i(index_3,:);
    dim = size(points_set_i);
    if dim(1) == 4
        frontier = [frontier; points_set_i(index_var2,:); points_set_i(index_var1,:)];
    else
        for k = 5:dim(1)
            frontier = [frontier; points_set_i(k,:)];
        end
        frontier = [frontier; points_set_i(index_var2,:); points_set_i(index_var1,:)];
    end
else
    index_copy = 4;
    frontier = points_set_i(index_var1,:);
    dim = size(points_set_i);
    if dim(1) == 4
        frontier = [frontier; points_set_i(index_3,:); points_set_i(index_var2,:)];
    else
        for k = 5:dim(1)
            frontier = [frontier; points_set_i(k,:)];
        end
        frontier = [frontier; points_set_i(index_3,:); points_set_i(index_var2,:)];
    end
end

outer_points = frontier*table_copy(j, 1:3)' > table_copy(j, 4);
if ~all(outer_points(:)==0) && ~all(outer_points(:)==1)
    index_first = find(outer_points(:), 1, 'first');
    index_last = find(outer_points(:), 1, 'last');
    if index_first > 1
        sol = findIntercept(frontier(index_first-1, :), frontier(index_first, :), table_copy, i, j);
        if ~isempty(sol)
            sol = [sol.x1 sol.y1 sol.z1];
            if isempty(symvar(sol))                                                
                x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                if all([x y z] > 0)                                    
                    if index_first == dim(1)-1
                        frontier_aux = frontier;
                        frontier = frontier(1:index_first-1,:);
                        frontier(index_first,:) = [x y z];
                        frontier = [frontier; frontier_aux(index_first:dim(1)-1,:)];
                    else
                        frontier(index_first, :) = [x y z];
                    end                   
                    m_j = (points_set_nonred{j}.points(index_2, index_comp2)-points_set_nonred{j}.points(index_1, index_comp2))/...
                            (points_set_nonred{j}.points(index_2, index_comp1)-points_set_nonred{j}.points(index_1, index_comp1));
                    m_i = (frontier(index_first-1, index_comp2)-frontier(index_first, index_comp2))/...
                        (frontier(index_first-1, index_comp1)-frontier(index_first, index_comp1));
                    
                    if index_copy == index_var1 || index_copy == index_3 %%%% 4
                        if (m_i > 0) && m_j < 0
                            index_1 = index_var1; index_2 = index_var2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = index_var2; index_2 = index_var1;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = index_var2; index_2 = index_var1;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = index_var2; index_2 = index_var1;
                        end
                    elseif index_copy == index_var2
                        if (m_i > 0) && m_j < 0
                            index_1 = index_var2; index_2 = index_var1;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = index_var1; index_2 = index_var2;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = index_var1; index_2 = index_var2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = index_var1; index_2 = index_var2;
                        end
                    end
                    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                        var_first = x;                
                    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                        var_first = y;
                    end
                    if index_1 == index_var1
                        if (points_set_j(index_2, index_comp1) <= var_first && points_set_j(index_2, index_comp2) >= z) && ...
                            (points_set_j(index_2, index_comp1) <= var_first <= points_set_j(index_1, index_comp1) && points_set_j(index_2, index_comp2) <= z <= points_set_j(index_1, index_comp2) || ...
                            points_set_j(index_1, index_comp1) <= var_first <= points_set_j(index_2, index_comp1) && points_set_j(index_1, index_comp2) <= z <= points_set_j(index_2, index_comp2))
                            if points_set_j(index_4, index_comp1) == points_set_j(index_2, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_2, index_comp2)
                                points_set_j(index_4, :) = [x y z];
                            else
                                points_set_j(index_3, :) = [x y z];
                            end
                            if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                                points_set_j(index_2, :) = [x 0 z];
                            elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                                points_set_j(index_2, :) = [0 y z];
                            end
                        else
                            if points_set_j(index_4, index_comp1) == points_set_j(index_2, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_2, index_comp2)
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                            end 
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)                                                                                                 
                                        frontier(index_first, :) = [x y z];
                                    end
                                end
                            end
                        end
                        
                        if index_first >= dim(1)-2 
                            if points_set_j(index_4, index_comp1) == points_set_j(index_2, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_2, index_comp2)
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_3,:), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_4,:), table_copy, i, j);
                            end        
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] >= 0)
                                        if index_first == dim(1)-1
                                            %if frontier(dim(1), 1) >= x && frontier(dim(1), 2) >= y
                                                frontier(dim(1), :) = [x y z];
                                            %end
                                        else
                                            %if frontier(dim(1)-1, 1) >= x && frontier(dim(1)-1, 2) >= y
                                                frontier(dim(1)-1, :) = [x y z];
                                            %end
                                        end
                                    end
                                end
                            end
                        end
                    elseif index_1 == index_var2
                        if (points_set_j(index_2, index_comp1) <= var_first && points_set_j(index_2, index_comp2) >= z && ...                            
                            (points_set_j(index_2, index_comp1) <= var_first <= points_set_j(index_1, index_comp1) && points_set_j(index_2, index_comp2) <= z <= points_set_j(index_1, index_comp2) || ...
                            points_set_j(index_1, index_comp1) <= var_first <= points_set_j(index_2, index_comp1) && points_set_j(index_1, index_comp2) <= z <= points_set_j(index_2, index_comp2)))
                        %if (points_set_j(index_1, 1) <= x && points_set_j(index_1, 2) >= y)
                            if points_set_j(index_4, index_comp1) == points_set_j(index_1, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_1, index_comp2)
                                points_set_j(index_3, :) = [x y z];
                            else
                                points_set_j(index_4, :) = [x y z];
                            end
                            if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                                points_set_j(index_2, :) = [x 0 z];
                            elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                                points_set_j(index_2, :) = [0 y z];
                            end
                            %minxyz(3) = min([minxyz(3) z]);
                        else
                            if points_set_j(index_4, index_comp1) == points_set_j(index_1, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_1, index_comp2)
                                %frontier(index_first, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) z_1];
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                            else
                                %frontier(index_first, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) z_1];
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                            end  
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)
                                        frontier(index_first,:) = [x y z];
                                    end
                                end
                            end
                        end                        

                        if index_first >= dim(1)-2 
                            if points_set_j(index_4, index_comp1) == points_set_j(index_1, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_1, index_comp2)
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_4,:), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_nonred{j}.points(index_1,:), points_set_nonred{j}.points(index_3,:), table_copy, i, j);
                            end
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] >= 0)
                                        if index_first == dim(1)-1
                                            %if frontier(dim(1), 1) >= x && frontier(dim(1), 2) >= y
                                                frontier(dim(1),:) = [x y z];
                                            %end
                                        else
                                            %if frontier(dim(1)-1, 1) >= x && frontier(dim(1)-1, 2) >= y
                                                frontier(dim(1)-1, :) = [x y z];
                                            %end
                                        end
                                    end
                                end
                            end
                        end
                    end                                        
                end
            end
        end
    end
    if index_last < dim(1)-1
        sol = findIntercept(frontier(index_last+1, :), frontier(index_last, :), table_copy, i, j);
        if ~isempty(sol)
            sol = [sol.x1 sol.y1 sol.z1];
            if isempty(symvar(sol))                                                
                x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                if all([x y z] > 0)                           
                    if index_last == index_first
                        frontier_aux = frontier;
                        frontier = frontier(1:index_last,:);
                        frontier(index_last+1,:) = [x y z];
                        frontier = [frontier; frontier_aux(index_last+1:dim(1)-1,:)];
                    else
                        frontier(index_last, :) = [x y z];
                    end
                    
                    m_j = (points_set_nonred{j}.points(index_2, index_comp2)-points_set_nonred{j}.points(index_1, index_comp2))/...
                            (points_set_nonred{j}.points(index_2, index_comp1)-points_set_nonred{j}.points(index_1, index_comp1));
                    m_i = (frontier(index_last+1, index_comp2)-frontier(index_last, index_comp2))/...
                        (frontier(index_last+1, index_comp1)-frontier(index_last, index_comp1));
                    
                    if index_copy == index_var1 || index_copy == index_3 %%%% 4
                        if (m_i > 0) && m_j < 0
                            index_1 = index_var1; index_2 = index_var2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = index_var2; index_2 = index_var1;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = index_var2; index_2 = index_var1;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = index_var2; index_2 = index_var1;
                        end
                    elseif index_copy == index_var2
                        if (m_i > 0) && m_j < 0
                            index_1 = index_var2; index_2 = index_var1;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j < 0
                            index_1 = index_var1; index_2 = index_var2;
                        elseif (m_i > 0) && m_j > 0
                            index_1 = index_var1; index_2 = index_var2;
                        elseif (m_i < 0 || abs(m_i) == Inf) && m_j > 0
                            index_1 = index_var1; index_2 = index_var2;
                        end
                    end
                    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                        var_first = x;                
                    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                        var_first = y;
                    end
                    if index_1 == index_var1
                        if (points_set_j(index_1, index_comp1) <= var_first && points_set_j(index_1, index_comp2) >= z && ...                            
                            (points_set_j(index_2, index_comp1) <= var_first <= points_set_j(index_1, index_comp1) && points_set_j(index_2, index_comp2) <= z <= points_set_j(index_1, index_comp2) || ...
                            points_set_j(index_1, index_comp1) <= var_first <= points_set_j(index_2, index_comp1) && points_set_j(index_1, index_comp2) <= z <= points_set_j(index_2, index_comp2)))
                        %if (points_set_j(index_2, 1) <= x && points_set_j(index_2, 2) >= y)
                            if points_set_j(index_4, index_comp1) == points_set_j(index_2, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_2, index_comp2)
                                points_set_j(index_3, :) = [x y z];
                            else
                                points_set_j(index_4, :) = [x y z];
                            end
                            if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                                points_set_j(index_1, :) = [x 0 z];
                            elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                                points_set_j(index_1, :) = [0 y z];
                            end
                            %minxyz(3) = min([minxyz(3) z]);
                        else
                            if points_set_j(index_4, index_comp1) == points_set_j(index_2, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_2, index_comp2)
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_3, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_4, :), table_copy, i, j);
                            end                                        
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)
                                        if index_last == index_first
                                            frontier(index_last+1,:) = [x y z];
                                        else
                                            frontier(index_last, :) = [x y z];
                                        end
                                    end
                                end
                            end
                        end
                        if index_last == 1 
                            if points_set_j(index_4, index_comp1) == points_set_j(index_2, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_2, index_comp2)
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                            end
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] >= 0)
                                        %if frontier(index_last, 1) >= x && frontier(index_last, 2) >= y
                                            frontier(index_last, :) = [x y z];
                                        %end
                                    end
                                end
                            end
                        end
                    elseif index_1 == index_var2
                        
                        if ((points_set_j(index_1, index_comp1) <= var_first && points_set_j(index_1, index_comp2) >= z) &&...
                            (points_set_j(index_2, index_comp1) <= var_first <= points_set_j(index_1, index_comp1) && points_set_j(index_2, index_comp2) <= z <= points_set_j(index_1, index_comp2) || ...
                            points_set_j(index_1, index_comp1) <= var_first <= points_set_j(index_2, index_comp1) && points_set_j(index_1, index_comp2) <= z <= points_set_j(index_2, index_comp2)))
                        %if (points_set_j(index_2, 1) <= x && points_set_j(index_2, 2) >= y)
                            if points_set_j(index_4, index_comp1) == points_set_j(index_1, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_1, index_comp2)
                                points_set_j(index_4, :) = [x y z];
                            else
                                points_set_j(index_3, :) = [x y z];
                            end
                            if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                                points_set_j(index_1, :) = [x 0 z];
                            elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                                points_set_j(index_1, :) = [0 y z];
                            end                            
                            %minxyz(3) = min([minxyz(3) z]);
                        else
                            if points_set_j(index_4, index_comp1) == points_set_j(index_1, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_1, index_comp2)
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_4, :), table_copy, i, j);
                            else
                                sol = findIntercept(points_set_j(index_1, :), points_set_j(index_3, :), table_copy, i, j);
                            end                                        
                            if ~(isempty(sol))
                                sol = [sol.x1 sol.y1 sol.z1];
                                if isempty(symvar(sol))
                                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                    if all([x y z] > 0)
                                        %if index_last == index_first
                                            frontier(index_last+1,:) = [x y z];
                                        %else
                                            frontier(index_last, :) = [x y z];
                                        %end
                                    end
                                end
                            end
                        end
                    end
                    if index_last == 1 
                        if points_set_j(index_4, index_comp1) == points_set_j(index_1, index_comp1) && points_set_j(index_4, index_comp2) == points_set_j(index_1, index_comp2)
                            sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
                        else
                            sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
                        end
                        if ~(isempty(sol))
                            sol = [sol.x1 sol.y1 sol.z1];
                            if isempty(symvar(sol))
                                x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                                if all([x y z] >= 0)
                                    %if frontier(index_last, 1) >= x && frontier(index_last, 2) >= y
                                        frontier(index_last, :) = [x y z];
                                    %end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if index_last - index_first > 1
        frontier = [frontier(1:index_first,:);frontier(index_last:dim(1)-1,:)];
    end
    dim = size(frontier);
    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
        points_set_i = points_set_i(2, :);
        if index_copy == index_var2
            points_set_i = [frontier(dim(1),:);points_set_i;frontier(dim(1)-1,:);frontier(1,:)];
        elseif index_copy == index_var1
            points_set_i = [frontier(dim(1)-1,:);points_set_i;frontier(dim(1),:);frontier(1,:)];
        else
            points_set_i = [frontier(dim(1),:);points_set_i;frontier(1,:);frontier(dim(1)-1,:)];
        end       
    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
        points_set_i = points_set_i(3, :);
        if index_copy == index_var2
            points_set_i = [frontier(dim(1),:);frontier(dim(1)-1,:);points_set_i;frontier(1,:)];
        elseif index_copy == index_var1
            points_set_i = [frontier(dim(1)-1,:);frontier(dim(1),:);points_set_i;frontier(1,:)];
        else
            points_set_i = [frontier(dim(1),:);frontier(1,:);points_set_i;frontier(dim(1)-1,:)];
        end
    end
    for k = 2:dim(1)-2
        points_set_i(3+k, :) = frontier(k,:);
    end
elseif all(outer_points(:)==1)
    points_set_i = points_set_i(1:4, :);
    if index_copy == index_var2
        points_set_i(index_var2,:) = points_set_i(4, :);
        points_set_i(4,:) = points_set_i(3, :);
    elseif index_copy == index_var1
        points_set_i(index_var1,:) = points_set_i(4, :);
        points_set_i(4,:) = points_set_i(3, :);
    else                        
        points_set_i(4,:) = points_set_i(3, :);
    end
    [points_set_i, points_set_j] = calcIntersect3_XZ_YZ(points_set_i, points_set_j, table_copy, i, j);
    %%%%%
end


%%%%%%%%
function [points_set_i, points_set_j] = calcIntersect2_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4)
%%%%%
if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
    index_var1 = 1; index_var2 = 3;
    index_comp1 = 1; index_comp2 = 3;
elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
    index_var1 = 1; index_var2 = 2;
    index_comp1 = 2; index_comp2 = 3;
end

sol = findIntercept(points_set_i(index_1, :), points_set_i(index_2, :), table_copy, i, j);
if ~isempty(sol) 
    sol = [sol.x1 sol.y1 sol.z1];
    if isempty(symvar(sol))                                                
        x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
        if all([x y z] > 0)
            if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                var_first = x;                
            elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                var_first = y;
            end
            if index_1 == index_var1
                if points_set_i(index_1, index_comp1) <= var_first && points_set_i(index_1, index_comp2) >= z 
                    points_set_i(4, :) = [x y z];        
                end
            else                
                if points_set_i(index_1, index_comp1) >= var_first && points_set_i(index_1, index_comp2) <= z 
                    points_set_i(4, :) = [x y z];                           
                end                
            end
            m_j = (points_set_j(index_1, index_comp2)-points_set_j(index_2, index_comp2))/...
                (points_set_j(index_1, index_comp1)-points_set_j(index_2, index_comp1));
            m_i = (points_set_i(index_1, index_comp2)-points_set_i(index_2, index_comp2))/...
                (points_set_i(index_1, index_comp1)-points_set_i(index_2, index_comp1));          
            
            outer_points = [points_set_i(index_1, :); points_set_i(index_2, :)]*table_copy(j, 1:3)' > table_copy(j, 4);
            index_first = find(outer_points(:)==0, 1, 'first');
            if (m_i < 0 && (m_j < 0 || abs(m_j) == Inf) && index_first == 1)
                index_first = index_var2; index_second = index_var1;
            elseif (m_i < 0 && (m_j < 0 || abs(m_j) == Inf) && index_first == 2)
                index_first = index_var1; index_second = index_var2;
            elseif (m_i < 0 && m_j > 0)
                index_first = index_var1; index_second = index_var2;            
            end
            
            if (points_set_j(index_first, index_comp1) <= var_first && points_set_j(index_first, index_comp2) >= z) || ...
                    (points_set_j(index_first, index_comp1) >= var_first && points_set_j(index_first, index_comp2) <= z) && ...
                    (points_set_j(index_2, index_comp1) <= var_first <= points_set_j(index_1, index_comp1) && points_set_j(index_2, index_comp2) <= z <= points_set_j(index_1, index_comp2) || ...
                            points_set_j(index_1, index_comp1) <= var_first <= points_set_j(index_2, index_comp1) && points_set_j(index_1, index_comp2) <= z <= points_set_j(index_2, index_comp2))
            %if points_set_j(index_1, 1) <= x && points_set_j(index_1, 2) >= y
                if points_set_j(index_3, index_comp1) == points_set_j(index_second, index_comp1) && points_set_j(index_3, index_comp2) == points_set_j(index_second, index_comp2)
                %if points_set_j(index_3, 1) == points_set_j(index_2, 1) && points_set_j(index_3, 2) == points_set_j(index_2, 2)
                    points_set_j(index_4, :) = [x y z];
                else
                    points_set_j(index_3, :) = [x y z];                  
                end
                if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                    points_set_j(index_first, :) = [x 0 z];
                elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                    points_set_j(index_first, :) = [0 y z];
                end                
                %minxyz(3) = min([minxyz(3) z]);                                
            else
                if points_set_j(index_3, index_comp1) == points_set_j(index_second, index_comp1) && points_set_j(index_3, index_comp2) == points_set_j(index_second, index_comp2)
                %if points_set_j(index_3, 1) == points_set_j(index_2, 1) && points_set_j(index_3, 2) == points_set_j(index_2, 2)
                    %points_set_i(4, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) z_1];
                    sol = findIntercept(points_set_j(index_first, :), points_set_j(index_4, :), table_copy, i, j);
                else
                    %points_set_i(4, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) z_1];
                    sol = findIntercept(points_set_j(index_first, :), points_set_j(index_3, :), table_copy, i, j);
                end                                
                if ~(isempty(sol))
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))
                        x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                        if all([x y z] > 0)
                            points_set_i(4, :) = [x y z];
                        end
                    end
                end
            end
            if points_set_j(index_3, index_comp1) == points_set_j(index_second, index_comp1) && points_set_j(index_3, index_comp2) == points_set_j(index_second, index_comp2)        
                sol = findIntercept(points_set_j(index_3, :), points_set_j(index_second, :), table_copy, i, j);
            else
                sol = findIntercept(points_set_j(index_4, :), points_set_j(index_second, :), table_copy, i, j);                    
            end
            if ~(isempty(sol))
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))
                    x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
                    if all([x y z] > 0)
                        points_set_i(index_2, :) = [x y z];
                    end
                end
            end
        end        
    end        
end

%%%%%
function [points_set_i, points_set_j] = calcIntersect3_XZ_YZ(points_set_i, points_set_j, table_copy, i, j)
global points_set_nonred;

if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
    index_var1 = 1; index_var2 = 3;
    index_comp1 = 1; index_comp2 = 3;
elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
    index_var1 = 1; index_var2 = 2;
    index_comp1 = 2; index_comp2 = 3;
end

if points_set_j(index_var1,index_comp1) > 0
    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
        sol = findIntercept(points_set_j(index_var1, :), points_set_j(index_4, :), table_copy, i, j);    
    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
        sol = findIntercept(points_set_j(index_var1, :), points_set_j(index_3, :), table_copy, i, j);    
    end
    
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] > 0)
                points_set_i(4, :) = [x y z]; %%%%index_3 = 4
                if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                    points_set_j(index_4, :) = [x y z];
                elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                    points_set_j(index_3, :) = [x y z];
                end                
            end
        end
    end
    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
        sol = findIntercept(points_set_nonred{j}.points(index_var1,:), points_set_nonred{j}.points(index_4,:), table_copy, i, j);
    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
        sol = findIntercept(points_set_nonred{j}.points(index_var1,:), points_set_nonred{j}.points(index_3,:), table_copy, i, j);
    end
    
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                points_set_i(index_var1, :) = [x y z];
            end
        end
    end
else
    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
        sol = findIntercept(points_set_j(index_var1, :), points_set_j(index_4, :), table_copy, i, j);    
    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
        sol = findIntercept(points_set_j(index_var1, :), points_set_j(index_3, :), table_copy, i, j);    
    end
    
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                points_set_i(index_var1, :) = [x y z];
                if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                    points_set_j(index_4, :) = [x y z];
                elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                    points_set_j(index_3, :) = [x y z];
                end                
            end
        end
    end
end

if points_set_j(index_var2,index_comp2) > 0
    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
        sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
        sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
    end    

    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] > 0)
                if all(points_set_i(4, :) == points_set_i(3, :))
                    points_set_i(4, :) = [x y z];
                else
                    points_set_i(5, :) = [x y z];
                end
                if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                    points_set_j(index_3, :) = [x y z];                    
                elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                    points_set_j(index_4, :) = [x y z];
                end 
            end
        end
    end
    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
        sol = findIntercept(points_set_nonred{j}.points(index_var2,:), points_set_nonred{j}.points(index_3,:), table_copy, i, j);
    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
        sol = findIntercept(points_set_nonred{j}.points(index_var2,:), points_set_nonred{j}.points(index_4,:), table_copy, i, j);
    end
    
    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                if index_2 == 3
                    if all(points_set_i(3, :) == points_set_i(4, :))
                        points_set_i(4, :) = [x y z];
                    end
                end
                points_set_i(index_2, :) = [x y z];                
            end
        end
    end
else
    if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
        sol = findIntercept(points_set_j(index_2, :), points_set_j(index_3, :), table_copy, i, j);
    elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
        sol = findIntercept(points_set_j(index_2, :), points_set_j(index_4, :), table_copy, i, j);
    end

    if ~(isempty(sol))
        sol = [sol.x1 sol.y1 sol.z1];
        if isempty(symvar(sol))
            x=eval(sol(1)); y=eval(sol(2)); z = eval(sol(3));
            if all([x y z] >= 0)
                if index_2 == 3
                    if all(points_set_i(3, :) == points_set_i(4, :))
                        points_set_i(4, :) = [x y z];
                    end
                end
                points_set_i(index_var2, :) = [x y z];
                if (index_1 == 1 && index_2 == 3) || (index_1 == 3 && index_2 == 1)
                    points_set_j(index_3, :) = [x y z];                    
                elseif (index_1 == 1 && index_2 == 2) || (index_1 == 2 && index_2 == 1)
                    points_set_j(index_4, :) = [x y z];
                end 
            end
        end
    end
end


%%%%%%
function [points_set_i,points_set_j, table_copy] = PlaneCrossInterceptYZ(handles, points_set_i, points_set_j, table_copy, i, j, ...
    index_1, index_2, index_3, index_4)
global points_set_nonred;
syms x1 y1 z1;
e = 0.00005; %%%% Cambiar si se cambia precisión

N2 = cross(points_set_j(1, :)-points_set_j(2, :), points_set_j(1, :)-points_set_j(3, :));
C2 = dot(N2, [1 0 0]);
N = cross(points_set_i(1, :)-points_set_i(2, :), points_set_i(1, :)-points_set_i(3, :));
C = dot(N, [1 0 0]);

if all(points_set_nonred{i}.points*table_copy(j, 1:3)' >= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' <= table_copy(i, 4))    
    points_set_i = points_set_j;
    
    table_copy(i, :) = table_copy(j, :);
    points_set_nonred{i} = points_set_nonred{j};
    %points_set_i_con.tipo = points_set_j_con.tipo;
elseif all(points_set_nonred{i}.points*table_copy(j, 1:3)' <= table_copy(j, 4)) && ...
        all(points_set_nonred{j}.points*table_copy(i, 1:3)' >= table_copy(i, 4))    
    points_set_j = points_set_i;
    
    table_copy(j, :) = table_copy(i, :);
    points_set_nonred{j} = points_set_nonred{i};
    %points_set_j_con.tipo = points_set_i_con.tipo;
elseif (abs(C) < e) && (abs(C2) < e)
    %%%% Desarrollar
    s = get(handles.listbox_operations, 'string');
    set(handles.listbox_operations, 'string', ...
         char(s, 'La reconstrucción puede ser incompleta.', 'Caso 4 de planos paralelos a Z en desarrollo.'));
elseif (abs(C2) < e && points_set_nonred{i}.tipo == 1)
    %if points_set_i(index_2, index_coor) < points_set_j(index_2, index_coor)                                                            
    if (points_set_i(2, 1) > 0 || (points_set_i(2, 1) == 0 && points_set_i(4, 1) > 0)) ||...
            (points_set_i(1, 1) > 0 || (points_set_i(1, 1) == 0 && points_set_i(4, 1) > 0))
        if index_1 == 1
            if all(points_set_i(3, :) == points_set_i(4, :))                
                [points_set_i, points_set_j] = calcIntersect2_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);            
            else                                    
                [points_set_i, points_set_j] = calcIntersect_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            end
        elseif index_1 == 2
            if all(points_set_i(3, :) == points_set_i(4, :))
                [points_set_i, points_set_j] = calcIntersect2_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            else
                [points_set_i, points_set_j] = calcIntersect_XZ_YZ(points_set_i, points_set_j, table_copy, i, j, index_1, index_2, index_3, index_4);
            end
        end                    
    elseif (points_set_j(1, 1) > 0 || (points_set_j(1, 1) == 0 && points_set_j(4, 1) > 0)) ||...
            (points_set_j(1, 1) > 0 || (points_set_j(2, 1) == 0 && points_set_j(4, 1) > 0))
        [points_set_i, points_set_j] = calcIntersect3_XZ_YZ(points_set_i, points_set_j, table_copy, i, j);        
    elseif points_set_i(1, 1) == 0 && points_set_i(2, 1) == 0
        [points_set_i, points_set_j] = calcIntersect3_XZ_YZ(points_set_i, points_set_j, table_copy, i, j);                
    end
elseif (abs(C) < e && points_set_nonred{j}.tipo == 1)   
    if (points_set_j(2, 1) > 0 || (points_set_j(2, 1) == 0 && points_set_j(4, 1) > 0)) ||...
            (points_set_j(1, 1) > 0 || (points_set_j(1, 1) == 0 && points_set_j(4, 1) > 0))
        if index_1 == 1           
            if all(points_set_j(3, :) == points_set_j(4, :))
                [points_set_j, points_set_i] = calcIntersect2_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            else
                [points_set_j, points_set_i] = calcIntersect_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            end
        elseif index_1 == 2
            if all(points_set_j(3, :) == points_set_j(4, :))                
                [points_set_j, points_set_i] = calcIntersect2_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            else
                [points_set_j, points_set_i] = calcIntersect_XZ_YZ(points_set_j, points_set_i, table_copy, j, i, index_1, index_2, index_3, index_4);
            end
        end                    
    elseif (points_set_i(2, 1) > 0 || (points_set_i(2, 1) == 0 && points_set_i(4, 1) > 0)) ||...
            (points_set_i(1, 1) > 0 || (points_set_i(1, 1) == 0 && points_set_i(4, 1) > 0))
        [points_set_j, points_set_i] = calcIntersect3_XZ_YZ(points_set_j, points_set_i, table_copy, j, i);
                
    elseif points_set_j(1, 1) == 0 && points_set_j(2, 1) == 0       
        [points_set_j, points_set_i] = calcIntersect3_XZ_YZ(points_set_j, points_set_i, table_copy, j, i);
    end
elseif ((abs(C) < e && points_set_nonred{j}.tipo ~= 1) || (abs(C2) < e && points_set_nonred{i}.tipo ~= 1))
    %%%%% Desarrollar
    s = get(handles.listbox_operations, 'string');
    set(handles.listbox_operations, 'string', ...
            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 7 de plano paralelo a Z y de tipo distinto a 1 con limites interiores en desarrollo.'));
end


%%%%%
function [points_set_i, points_set_j] = ZParalelPlanesInterceptXY(handles, points_set_i, points_set_j, ...
    x, y, table_copy, i, j, index_1, index_2, index_3, index_4)
%global minxyz;
e = 0.00005; %%%% Cambiar si se cambia precisión

syms x1 y1 z1;

N2 = cross(points_set_j(1,:)-points_set_j(2,:), points_set_j(1,:)-points_set_j(3,:));
N = cross(points_set_i(1,:)-points_set_i(2,:), points_set_i(1,:)-points_set_i(3,:));   
C2 = dot(N2, [0 0 1]); C = dot(N, [0 0 1]);

if (abs(C2) < e) && (abs(C) < e)         
    if points_set_i(index_3,1) == points_set_i(index_1,1) && ...
            points_set_i(index_3,2) == points_set_i(index_1,2)
        if points_set_i(index_4, 3) == points_set_i(index_3, 3) %%%% Prueba 20/01  22/01
            %minxyz(3) = min([points_set_j(4, 3) points_set_i(index_4, 3)]);
            if points_set_j(index_4, 3) == points_set_j(index_3, 3)                
                points_set_i(index_3, :) = [x y points_set_i(index_4, 3)];
                points_set_i(index_4, :) = [points_set_i(index_4, 1) points_set_i(index_4, 2) points_set_i(index_4, 3)];

                if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x                    
                    points_set_j(index_4, :) = [x y points_set_j(index_4, 3)];
                    points_set_j(index_3, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) points_set_j(index_3, 3)];
                end
            else                    
                %x1 = x; y1 = y; syms z1;
                sol = findIntercept(points_set_j(index_3,:), points_set_j(index_4,:), table_copy, i, j);
                if ~isempty(sol)
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))                                                
                        x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                        if all([x_1 y_1 z_1] > 0)
                            if points_set_j(3, 2) <= y_1 && points_set_j(3, 1) >= x_1
                                points_set_j(index_4, :) = [x_1 y_1 z_1];
                            end
                            
                            points_set_i(index_3, :) = [x_1 y_1 z_1];
                            points_set_i(index_4, :) = [points_set_i(index_4, 1) points_set_i(index_4, 2) z_1];
                        end
                    else
                        if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                           points_set_j(index_4, :) = [x y points_set_j(index_4, 3)];
                        end
                        
                        points_set_i(index_3, :) = [x y points_set_i(index_3, 3)];
                        
                        %%%%% Desarrollar
                        s = get(handles.listbox_operations, 'string');
                        set(handles.listbox_operations, 'string', ...
                            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                    end
                else
                    if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                        points_set_j(index_4, :) = [x y points_set_j(index_4, 3)];
                    end
                    
                    points_set_i(index_3, :) = [x y points_set_i(index_3, 3)];
                    
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                         char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            end
        else
            %x1 = x; y1 = y; syms z1; 
            sol = findIntercept(points_set_i(index_3,:), points_set_i(index_4,:), table_copy, i, j);
            if ~isempty(sol)
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))                                                
                    x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                    if all([x_1 y_1 z_1] > 0)                                                    
                        points_set_i(index_3, :) = [x_1 y_1 z_1];

                        if points_set_j(3, 2) <= y_1 && points_set_j(3, 1) >= x_1
                            %minxyz(3) = min([minxyz(3) z]);
                            if points_set_j(index_4, 3) == points_set_j(index_3, 3)
                                points_set_j(index_4, :) = [x_1 y_1 z_1];
                                points_set_j(index_3, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) z_1]; %%%
                            else
                                points_set_j(index_4, :) = [x_1 y_1 z_1];                                
                            end
                        end
                    end
                else
                    points_set_i(index_3, :) = [x y points_set_i(index_3, 3)];
                    %minxyz(3) = min([minxyz(3) z]);
                    if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                        if points_set_j(index_4, 3) == points_set_j(index_3, 3)
                            points_set_j(index_4, :) = [x y points_set_j(index_4, 3)];                            
                        else
                            points_set_j(index_4, :) = [x y points_set_j(index_4, 3)];                            
                        end
                    end
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            else
                points_set_i(index_3, :) = [x y points_set_i(index_3, 3)];
                %minxyz(3) = min([minxyz(3) z]);
                if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                    if points_set_j(index_4, 3) == points_set_j(index_3, 3)
                        points_set_j(index_4, :) = [x y points_set_j(index_4, 3)];                        
                    else
                        points_set_j(index_4, :) = [x y points_set_j(index_4, 3)];                        
                    end
                end
                s = get(handles.listbox_operations, 'string');
                set(handles.listbox_operations, 'string', ...
                    char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
            end
        end
    else
        if points_set_i(index_4, 3) == points_set_i(index_3, 3) %%%% Prueba 20/01  22/01
            %minxyz(3) = min([points_set_j(4, 3) points_set_i(index_4, 3)]);
            if points_set_j(index_4, 3) == points_set_j(index_3, 3)                
                points_set_i(index_4, :) = [x y points_set_i(index_4, 3)];
                points_set_i(index_3, :) = [points_set_i(index_3, 1) points_set_i(index_3, 2) points_set_i(index_3, 3)];

                if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x                    
                    points_set_j(index_3, :) = [x y points_set_j(index_3, 3)];
                    points_set_j(index_4, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) points_set_j(index_4, 3)];
                end
            else                    
                %x1 = x; y1 = y; syms z1;
                sol = findIntercept(points_set_j(index_3,:), points_set_j(index_4,:), table_copy, i, j);
                if ~isempty(sol)
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))                                                
                        x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                        if all([x_1 y_1 z_1] > 0) 
                            if points_set_j(3, 2) <= y_1 && points_set_j(3, 1) >= x_1
                                points_set_j(index_3, :) = [x_1 y_1 z_1];
                            end
                            
                            points_set_i(index_4, :) = [x_1 y_1 z_1];
                            points_set_i(index_3, :) = [points_set_i(index_3, 1) points_set_i(index_3, 2) z_1];
                        end
                    else
                        if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                            points_set_j(index_3, :) = [x y points_set_j(index_3, 3)];
                        end
                        
                        points_set_i(index_4, :) = [x y points_set_i(index_4, 3)];
                        
                        s = get(handles.listbox_operations, 'string');
                        set(handles.listbox_operations, 'string', ...
                            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                    end
                else
                    if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                        points_set_j(index_3, :) = [x y points_set_j(index_3, 3)];
                    end
                    
                    points_set_i(index_4, :) = [x y points_set_i(index_4, 3)];
                    
                    s = get(handles.listbox_operations, 'string');
                    set(s, handles.listbox_operations, 'string', ...
                            char('La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            end
        else
            %x1 = x; y1 = y; syms z1;
            sol = findIntercept(points_set_i(index_3,:), points_set_i(index_4,:), table_copy, i, j);
            if ~isempty(sol)
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))                                                
                    x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                    if all([x_1 y_1 z_1] > 0)                                                    
                        points_set_i(index_4, :) = [x_1 y_1 z_1];
                        %minxyz(3) = min([minxyz(3) z]);
                        if points_set_j(3, 2) <= y_1 && points_set_j(3, 1) >= x_1
                            if points_set_j(index_4, 3) == points_set_j(index_3, 3)
                                points_set_j(index_3, :) = [x_1 y_1 z_1];
                                points_set_j(index_4, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) z_1];
                            else
                                points_set_j(index_3, :) = [x_1 y_1 z_1];                                
                            end
                        end
                    end
                else
                    points_set_i(index_4, :) = [x y points_set_i(index_4, 3)];
                    %minxyz(3) = min([minxyz(3) z]);
                    if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                        if points_set_j(index_4, 3) == points_set_j(index_3, 3)
                            points_set_j(index_3, :) = [x y points_set_j(index_3, 3)];                            
                        else
                            points_set_j(index_3, :) = [x y points_set_j(index_3, 3)];                            
                        end
                    end
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                         char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            else
                points_set_i(index_4, :) = [x y points_set_i(index_4, 3)];
                %minxyz(3) = min([minxyz(3) z]);
                if points_set_j(3, 2) <= y && points_set_j(3, 1) >= x
                    if points_set_j(index_4, 3) == points_set_j(index_3, 3)
                        points_set_j(index_3, :) = [x y points_set_j(index_3, 3)];                        
                    else
                        points_set_j(index_3, :) = [x y points_set_j(index_3, 3)];                        
                    end
                end
                s = get(handles.listbox_operations, 'string');
                set(handles.listbox_operations, 'string', ...
                     char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
            end
        end
    end
elseif (abs(C) < e)
    C2_Y = dot(N2, [0 1 0]); C2_X = dot(N2, [1 0 0]); C2_Z = dot(N2, [0 0 1]);

    if (abs(C2_Z) < e)
        %%%% Desarrollar
        s = get(handles.listbox_operations, 'string');
        set(handles.listbox_operations, 'string', ...
           char(s, 'La reconstrucción puede ser incompleta.', 'Caso 15 de colapso de plano paralelo a Z.'));
    %    if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
    %        if points_set_j(index_3, 1) == points_set_j(index_1, 1) && points_set_j(index_3, 2) == points_set_j(index_1, 2)
    %            points_set_j(index_3, :) = [x y 0];
    %        end
    %    else
    %        if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
    %            points_set_j(index_4, :) = [x y 0];
    %        end
    %    end
    end
    if (index_1 == 2 && index_2 == 3) && (abs(C2_X) < e)
        if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
            if x < points_set_j(index_3, 1)
                points_set_j(index_3, :) = [x points_set_j(index_3, 2) points_set_j(index_3, 3)];
            end                
        else
            if x < points_set_j(index_4, 1)            
                points_set_j(index_4, :) = [x points_set_j(index_4, 2) points_set_j(index_4, 3)];
            end
        end
    elseif (index_1 == 3 && index_2 == 2)  && (abs(C2_Y) < e)
        if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
            if y < points_set_j(index_3, 2)                
                points_set_j(index_3, :) = [points_set_j(index_3, 1) y points_set_j(index_3, 3)];
            end
        else
            if y < points_set_j(index_4, 2)
                points_set_j(index_4, :) = [points_set_j(index_4, 1) y points_set_j(index_4, 3)];
            end
        end
    end
elseif (abs(C2) < e)
    C_Y = dot(N, [0 1 0]); C_X = dot(N, [1 0 0]); C_Z = dot(N, [0 0 1]);
    
    if (abs(C_Z) < e)
        %%%% Desarrollar
        s = get(handles.listbox_operations, 'string');
        set(handles.listbox_operations, 'string', ...
               char(s, 'La reconstrucción puede ser incompleta.', 'Caso 15 de colapso de plano paralelo a Z.'));
    %    if points_set_i(index_3, 1) == points_set_i(index_2, 1) && points_set_i(index_3, 2) == points_set_i(index_2, 2)
    %        if points_set_i(index_4, 1) == points_set_i(index_1, 1) && points_set_i(index_4, 2) == points_set_i(index_1, 2)
    %            points_set_i(index_4, :) = [x y 0];
    %        end
    %    else
    %        if points_set_i(index_3, 1) == points_set_i(index_1, 1) && points_set_i(index_3, 2) == points_set_i(index_1, 2)
    %            points_set_i(index_3, :) = [x y 0];
    %        end
    %    end
    end
    if (index_1 == 2 && index_2 == 3) && (abs(C_Y) < e)
        if points_set_i(index_3, 1) == points_set_i(index_2, 1) && points_set_i(index_3, 2) == points_set_i(index_2, 2)
            if y < points_set_i(index_4, 2)
                points_set_i(index_4, :) = [points_set_i(index_4, 1) y points_set_i(index_4, 3)];
            end                
        else
            if y < points_set_i(index_3, 2)            
                points_set_i(index_3, :) = [points_set_i(index_3, 1) y points_set_i(index_3, 3)];
            end
        end
    elseif (index_1 == 3 && index_2 == 2)  && (abs(C_X) < e)              
        if points_set_i(index_3, 1) == points_set_i(index_2, 1) && points_set_i(index_3, 2) == points_set_i(index_2, 2)
            if x < points_set_i(index_4, 1)                
                points_set_i(index_4, :) = [x points_set_i(index_4, 2) points_set_i(index_4, 3)];
            end
        else
            if x < points_set_i(index_3, 1)
                points_set_i(index_3, :) = [x points_set_i(index_3, 2) points_set_i(index_3, 3)];
            end
        end  
    end        
end


function [points_set_i,points_set_j] = YParalelPlanesInterceptXZ(handles, points_set_i, points_set_j, ...
    x, z, table_copy, i, j, index_1, index_2, index_3, index_4)
%global minxyz;
e = 0.00005; %%%% Cambiar si se cambia precisión

syms x1 y1 z1;

N2 = cross(points_set_j(1,:)-points_set_j(2,:), points_set_j(1,:)-points_set_j(3,:));
N = cross(points_set_i(1,:)-points_set_i(2,:), points_set_i(1,:)-points_set_i(3,:));
C2 = dot(N2, [0 1 0]); C = dot(N, [0 1 0]);

if (abs(C2) < e) && (abs(C) < e)
    if points_set_i(index_3,1) == points_set_i(index_1,1) && ...
            points_set_i(index_3,3) == points_set_i(index_1,3)
        if points_set_i(index_4, 2) == points_set_i(index_3, 2) %%%% Prueba 20/01  22/01               
            if points_set_j(index_4, 2) == points_set_j(index_3, 2)                
                points_set_i(index_3, :) = [x points_set_i(index_3, 2), z];
                points_set_i(index_4, :) = [points_set_i(index_4, 1) points_set_i(index_4, 2) points_set_i(index_4, 3)];
                
                %if points_set_j(3, 3) <= z && points_set_j(3, 1) >= x
                %%%%Revisar esta condición y el resto de la función
                    points_set_j(index_4, :) = [x points_set_j(index_4, 2) z];
                    points_set_j(index_3, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) points_set_j(index_3, 3)];
                %end
            else                   
                %x1 = x; z1 = z; syms y1;
                sol = findIntercept(points_set_j(index_3,:), points_set_j(index_4,:), table_copy, i, j);
                if ~isempty(sol)
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))                                                
                        x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                        if all([x_1 y_1 z_1] > 0)                                                    
                            points_set_j(index_4, :) = [x_1 y_1 z_1];
                            %minxyz(3) = min([minxyz(3) z]);                            
                            points_set_i(index_3, :) = [x_1 y_1 z_1];
                            points_set_i(index_4, :) = [points_set_i(index_4, 1) y_1 points_set_i(index_4, 3)];
                        end
                    else
                        points_set_j(index_4, :) = [x points_set_j(index_4, 2) z];

                        points_set_i(index_3, :) = [x points_set_i(index_3, 2) z];
                        s = get(handles.listbox_operations, 'string');
                        set(handles.listbox_operations, 'string', ...
                            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                    end
                else
                    points_set_j(index_4, :) = [x points_set_j(index_4, 2) z];

                    points_set_i(index_3, :) = [x points_set_i(index_3, 2) z];
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            end
        else
            %x1 = x; z1 = z; syms y1;
            sol = findIntercept(points_set_i(index_3,:), points_set_i(index_4,:), table_copy, i, j);
            if ~isempty(sol)
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))                                                
                    x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                    if all([x_1 y_1 z_1] > 0)                                                    
                        points_set_i(index_3, :) = [x_1 y_1 z_1];
                        %minxyz(3) = min([minxyz(3) z]);
                        if points_set_j(index_4, 2) == points_set_j(index_3, 2)
                            points_set_j(index_4, :) = [x_1 y_1 z_1];
                            points_set_j(index_3, :) = [points_set_j(index_3, 1) y_1 points_set_j(index_3, 3)];
                        else
                            points_set_j(index_4, :) = [x_1 y_1 z_1];
                        end
                    end
                else
                    points_set_i(index_3, :) = [x points_set_i(index_3, 2) z];
                    %minxyz(3) = min([minxyz(3) z]);
                    if points_set_j(index_4, 2) == points_set_j(index_3, 2)
                        points_set_j(index_4, :) = [x points_set_j(index_4, 2) z];                            
                    else
                        points_set_j(index_4, :) = [x points_set_j(index_4, 2) z];
                    end
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            else
                points_set_i(index_3, :) = [x points_set_i(index_3, 2) z];
                %minxyz(3) = min([minxyz(3) z]);
                if points_set_j(index_4, 2) == points_set_j(index_3, 2)
                    points_set_j(index_4, :) = [x points_set_j(index_4, 2) z];                            
                else
                    points_set_j(index_4, :) = [x points_set_j(index_4, 2) z];
                end
                s = get(handles.listbox_operations, 'string');
                set(handles.listbox_operations, 'string', ...
                     char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
            end
        end
    else
        if points_set_i(index_4, 2) == points_set_i(index_3, 2) %%%% Prueba 20/01  22/01               
            if points_set_j(index_4, 2) == points_set_j(index_3, 2)                
                points_set_i(index_4, :) = [x points_set_i(index_4, 2) z];
                points_set_i(index_3, :) = [points_set_i(index_3, 1) points_set_i(index_3, 2) points_set_i(index_3, 3)];

                points_set_j(index_3, :) = [x points_set_j(index_3, 2) z];
                points_set_j(index_4, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) points_set_j(index_4, 3)];
            else
                %x1 = x; z1 = z; syms y1;
                sol = findIntercept(points_set_j(index_3,:), points_set_j(index_4,:), table_copy, i, j);
                if ~isempty(sol)
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))                                                
                        x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                        if all([x_1 y_1 z_1] > 0)                                                    
                            points_set_j(index_3, :) = [x_1 y_1 z_1];
                            %minxyz(3) = min([minxyz(3) z]);
                            points_set_i(index_4, :) = [x_1 y_1 z_1];
                            points_set_i(index_3, :) = [points_set_i(index_3, 1) y_1 points_set_i(index_3, 3)];
                        end
                    else
                        points_set_j(index_3, :) = [x points_set_j(index_3, 2) z];

                        points_set_i(index_4, :) = [x points_set_i(index_4, 2) z];
                        s = get(handles.listbox_operations, 'string');
                        set(handles.listbox_operations, 'string', ...
                            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                    end
                else
                    points_set_j(index_3, :) = [x points_set_j(index_3, 2) z];

                    points_set_i(index_4, :) = [x points_set_i(index_4, 2) z];
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            end
        else
            %x1 = x; z1 = z; syms y1;
            sol = findIntercept(points_set_i(index_3,:), points_set_i(index_4,:), table_copy, i, j);
            if ~isempty(sol)
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))                                                
                    x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                    if all([x_1 y_1 z_1] > 0)                                                    
                        points_set_i(index_4, :) = [x_1 y_1 z_1];
                        %minxyz(3) = min([minxyz(3) z]);
                        if points_set_j(index_4, 2) == points_set_j(index_3, 2)
                            points_set_j(index_3, :) = [x_1 y_1 z_1];
                            points_set_j(index_4, :) = [points_set_j(index_4, 1) y_1 points_set_j(index_4, 3)];
                        else
                            points_set_j(index_3, :) = [x_1 y_1 z_1];
                        end
                    end
                else
                    points_set_i(index_4, :) = [x points_set_i(index_4, 2) z];
                    %minxyz(3) = min([minxyz(3) z]);
                    if points_set_j(index_4, 2) == points_set_j(index_3, 2)
                        points_set_j(index_3, :) = [x points_set_j(index_3, 2) z];                        
                    else
                        points_set_j(index_3, :) = [x points_set_j(index_3, 2) z];
                    end
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            else
                points_set_i(index_4, :) = [x points_set_i(index_4, 2) z];
                %minxyz(3) = min([minxyz(3) z]);
                if points_set_j(index_4, 2) == points_set_j(index_3, 2)
                    points_set_j(index_3, :) = [x points_set_j(index_3, 2) z];                    
                else
                    points_set_j(index_3, :) = [x points_set_j(index_3, 2) z];
                end
                s = get(handles.listbox_operations, 'string');
                set(handles.listbox_operations, 'string', ...
                    char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
            end
        end
    end
elseif (abs(C) < e)
    C2_Z = dot(N2, [0 0 1]); C2_X = dot(N2, [1 0 0]); C2_Y = dot(N2, [0 1 0]);


    if (abs(C2_Y) < e)
        %%%% Desarrollar
        s = get(handles.listbox_operations, 'string');
        set(handles.listbox_operations, 'string', ...
            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 16 de colapso de plano paralelo a Y.'));
    %    if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 2) == points_set_j(index_1, 2)
    %        if points_set_j(index_3, 1) == points_set_j(index_1, 1) && points_set_j(index_3, 3) == points_set_j(index_1, 3)
    %            points_set_j(index_3, :) = [x 0 z];
    %        end
    %    else
    %        if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
    %            points_set_j(index_4, :) = [x 0 z];
    %        end
    %    end
    end
    if ((index_1 == 1 && index_2 == 3)||(index_1 == 4 && index_2 == 3)) && (abs(C2_X) < e)
        if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
            %minxyz(3) = min([minxyz(3) min([z points_set_i(index_3, 3)])]);
            if x < points_set_j(index_3, 1)                             
                points_set_j(index_3, :) = [x points_set_j(index_3, 2) points_set_j(index_3, 3)];
            end
        else
            if x < points_set_j(index_4, 1)
                points_set_j(index_4, :) = [x points_set_j(index_4, 2) points_set_j(index_4, 3)];
            end        
        end
    elseif ((index_1 == 3 && index_2 == 1)||(index_1 == 4 && index_2 == 1))  && (abs(C2_Z) < e)
        if points_set_j(index_4, 1) == points_set_j(index_1, 1) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
            %minxyz(1) = min([minxyz(1) min([x  points_set_i(index_3, 1)])]);
            if z < points_set_j(index_3, 3)
                points_set_j(index_3, :) = [points_set_j(index_3, 1)  points_set_j(index_3, 2) z];
            end        
        else
            if z < points_set_j(index_4, 3)
                points_set_j(index_4, :) = [points_set_j(index_4, 1)  points_set_j(index_4, 2) z];
            end        
        end  
    end
elseif (abs(C2) < e)
    C_Z = dot(N, [0 0 1]); C_X = dot(N, [1 0 0]); C_Y = dot(N, [0 1 0]);
    
    if (abs(C_Y) < e)
        %%%% Desarrollar
        s = get(handles.listbox_operations, 'string');
        set(handles.listbox_operations, 'string', ...
             char(s, 'La reconstrucción puede ser incompleta.', 'Caso 16 de colapso de plano paralelo a Y.'));
    %    if points_set_i(index_3, 1) == points_set_i(index_2, 1) && points_set_i(index_3, 3) == points_set_i(index_2, 3)
    %        if points_set_i(index_4, 1) == points_set_i(index_1, 1) && points_set_i(index_4, 3) == points_set_i(index_1, 3)
    %            points_set_i(index_4, :) = [x 0 z];
    %        end
    %    else
    %        if points_set_i(index_3, 1) == points_set_i(index_1, 1) && points_set_i(index_3, 3) == points_set_i(index_1, 3)
    %            points_set_i(index_3, :) = [x 0 z];
    %        end
    %    end
    end
    if ((index_1 == 1 && index_2 == 3)||(index_1 == 4 && index_2 == 3)) && (abs(C_Z) < e)
        if points_set_i(index_3, 1) == points_set_i(index_2, 1) && points_set_i(index_3, 3) == points_set_i(index_2, 3)
            %minxyz(3) = min([minxyz(3) min([z points_set_i(index_3, 3)])]);
            if z < points_set_i(index_4, 3)                             
                points_set_i(index_4, :) = [points_set_i(index_4, 1) points_set_i(index_4, 2) z];
            end        
        else
            if z < points_set_i(index_3, 3)
                points_set_i(index_3, :) = [points_set_i(index_3, 1) points_set_i(index_3, 2) z];
            end        
        end
    elseif ((index_1 == 3 && index_2 == 1)||(index_1 == 4 && index_2 == 1))  && (abs(C_X) < e)
        if points_set_i(index_3, 1) == points_set_i(index_2, 1) && points_set_i(index_3, 3) == points_set_i(index_2, 3)
            %minxyz(1) = min([minxyz(1) min([x  points_set_i(index_3, 1)])]);
            if x < points_set_i(index_4, 1)            
                points_set_i(index_4, :) = [x  points_set_i(index_4, 2)  points_set_i(index_4, 3)];
            end        
        else
            if x < points_set_i(index_3, 1)            
                points_set_i(index_3, :) = [x  points_set_i(index_3, 2)  points_set_i(index_3, 3)];
            end        
        end  
    end
end


function [points_set_i,points_set_j] = XParalelPlanesInterceptYZ(handles, points_set_i, points_set_j, ...
    y, z, table_copy, i, j, index_1, index_2, index_3, index_4)
%global minxyz;
e = 0.00005; %%%% Cambiar si se cambia precisión

%syms x1 y1 z1;
N2 = cross(points_set_j(1,:)-points_set_j(2,:), points_set_j(1,:)-points_set_j(3,:));
N = cross(points_set_i(1,:)-points_set_i(2,:), points_set_i(1,:)-points_set_i(3,:));
C2 = dot(N2, [1 0 0]); C = dot(N, [1 0 0]);

if (abs(C2) < e) && (abs(C) < e)
    if points_set_i(index_3,2) == points_set_i(index_1,2) && ...
            points_set_i(index_3,3) == points_set_i(index_1,3)
        if points_set_i(index_4, 1) == points_set_i(index_3, 1) %%%% Prueba 20/01  22/01
            if points_set_j(index_4, 1) == points_set_j(index_3, 1) %%%% Prueba 20/01  22/01                
                points_set_i(index_3, :) = [points_set_i(index_3, 1) y z];
                points_set_i(index_4, :) = [points_set_i(index_4, 1) points_set_i(index_4, 2)  points_set_i(index_4, 3)];
                
                %if points_set_j(2, 3) <= z && points_set_j(2, 2) >= y
                %%%%Revisar esta condición y el resto de la función
                    points_set_j(index_4, :) = [points_set_j(index_4, 1) y z];
                    points_set_j(index_3, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) points_set_j(index_3, 3)];
                %end
            else
                %y1 = y; z1 = z; syms x1;
                sol = findIntercept(points_set_j(index_3,:), points_set_j(index_4,:), table_copy, i, j);
                if ~isempty(sol)
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))                                                
                        x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                        if all([x_1 y_1 z_1] > 0)                                                    
                            points_set_j(index_4, :) = [x_1 y_1 z_1];
                            %minxyz(3) = min([minxyz(3) z]);
                            points_set_i(index_3, :) = [x_1 y_1 z_1];
                            points_set_i(index_4, :) = [x_1 points_set_i(index_4, 2)  points_set_i(index_4, 3)];
                        end
                    else
                        points_set_j(index_4, :) = [points_set_j(index_4, 1) y z];
                        %minxyz(3) = min([minxyz(3) z]);
                        points_set_i(index_3, :) = [points_set_i(index_3, 1) y z];
                        s = get(handles.listbox_operations, 'string');
                        set(handles.listbox_operations, 'string', ...
                            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                    end
                else                        
                    points_set_j(index_4, :) = [points_set_j(index_4, 1) y z];
                    %minxyz(3) = min([minxyz(3) z]);
                    points_set_i(index_3, :) = [points_set_i(index_3, 1) y z];
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            end
        else                
            %y1 = y; z1 = z; syms x1;
            sol = findIntercept(points_set_i(index_3,:), points_set_i(index_4,:), table_copy, i, j);
            if ~isempty(sol)
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))                                                
                    x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                    if all([x_1 y_1 z_1] > 0)                                                    
                        points_set_i(index_3, :) = [x_1 y_1 z_1];
                        %minxyz(3) = min([minxyz(3) z]);
                        if points_set_j(index_4, 1) == points_set_j(index_3, 1)
                            points_set_j(index_4, :) = [x_1 y_1 z_1];
                            points_set_j(index_3, :) = [x_1 points_set_j(index_3, 2) points_set_j(index_3, 3)];
                        else
                            points_set_j(index_4, :) = [x_1 y_1 z_1];
                        end
                    end
                else
                    points_set_i(index_3, :) = [points_set_i(index_3, 1) y z];
                    %minxyz(3) = min([minxyz(3) z]);
                    if points_set_j(index_4, 1) == points_set_j(index_3, 1)
                        points_set_j(index_4, :) = [points_set_j(index_4, 1) y z];                            
                    else
                        points_set_j(index_4, :) = [points_set_j(index_4, 1) y z];
                    end
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            else
                points_set_i(index_3, :) = [points_set_i(index_3, 1) y z];
                %minxyz(3) = min([minxyz(3) z]);
                if points_set_j(index_4, 1) == points_set_j(index_3, 1)
                    points_set_j(index_4, :) = [points_set_j(index_4, 1) y z];                            
                else
                    points_set_j(index_4, :) = [points_set_j(index_4, 1) y z];
                end
                s = get(handles.listbox_operations, 'string');
                set(handles.listbox_operations, 'string', ...
                    char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
            end
        end
    else
        if points_set_i(index_4, 1) == points_set_i(index_3, 1) %%%% Prueba 20/01  22/01
            if points_set_j(index_4, 1) == points_set_j(index_3, 1) %%%% Prueba 20/01  22/01                
                points_set_i(index_4, :) = [points_set_i(index_4, 1) y z];
                points_set_i(index_3, :) = [points_set_i(index_3, 1) points_set_i(index_3, 2) points_set_i(index_3, 3)];

                points_set_j(index_3, :) = [points_set_j(index_3, 1) y z];
                points_set_j(index_4, :) = [points_set_j(index_4, 1) points_set_j(index_4, 2) points_set_j(index_4, 3)];
            else                   
                %y1 = y; z1 = z; syms x1;
                sol = findIntercept(points_set_j(index_3,:), points_set_j(index_4,:), table_copy, i, j);
                if ~isempty(sol)
                    sol = [sol.x1 sol.y1 sol.z1];
                    if isempty(symvar(sol))                                                
                        x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                        if all([x_1 y_1 z_1] > 0)                                                    
                            points_set_j(index_3, :) = [x_1 y_1 z_1];
                            %minxyz(3) = min([minxyz(3) z]);
                            points_set_i(index_4, :) = [x_1 y_1 z_1];
                            points_set_i(index_3, :) = [x_1 points_set_i(index_3, 2) points_set_i(index_3, 3)];
                        end
                    else
                        points_set_j(index_3, :) = [points_set_j(index_3, 1) y z];
                        %minxyz(3) = min([minxyz(3) z]);
                        points_set_i(index_4, :) = [points_set_i(index_4, 1) y z];
                        s = get(handles.listbox_operations, 'string');
                        set(handles.listbox_operations, 'string', ...
                            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                    end
                else
                    points_set_j(index_3, :) = [points_set_j(index_3, 1) y z];
                    %minxyz(3) = min([minxyz(3) z]);
                    points_set_i(index_4, :) = [points_set_i(index_4, 1) y z];
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            end
        else
            %y1 = y; z1 = z; syms x1;
            sol = findIntercept(points_set_i(index_3,:), points_set_i(index_4,:), table_copy, i, j);
            if ~isempty(sol)
                sol = [sol.x1 sol.y1 sol.z1];
                if isempty(symvar(sol))                                                
                    x_1=eval(sol(1)); y_1=eval(sol(2)); z_1 = eval(sol(3));
                    if all([x_1 y_1 z_1] > 0)                                                    
                        points_set_i(index_4, :) = [x_1 y_1 z_1];
                        %minxyz(3) = min([minxyz(3) z]);
                        if points_set_j(index_4, 1) == points_set_j(index_3, 1)
                            points_set_j(index_3, :) = [x_1 y_1 z_1];
                            points_set_j(index_4, :) = [x_1 points_set_j(index_4, 2) points_set_j(index_4, 3)];
                        else
                            points_set_j(index_3, :) = [x_1 y_1 z_1];
                        end
                    end
                else
                    points_set_i(index_4, :) = [points_set_i(index_4, 1) y z];
                    %minxyz(3) = min([minxyz(3) z]);
                    if points_set_j(index_4, 1) == points_set_j(index_3, 1)
                        points_set_j(index_3, :) = [points_set_j(index_3, 1) y z];                            
                    else
                        points_set_j(index_3, :) = [points_set_j(index_3, 1) y z];
                    end
                    s = get(handles.listbox_operations, 'string');
                    set(handles.listbox_operations, 'string', ...
                        char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
                end
            else
                points_set_i(index_4, :) = [points_set_i(index_4, 1) y z];
                %minxyz(3) = min([minxyz(3) z]);
                if points_set_j(index_4, 1) == points_set_j(index_3, 1)
                    points_set_j(index_3, :) = [points_set_j(index_3, 1) y z];                            
                else
                    points_set_j(index_3, :) = [points_set_j(index_3, 1) y z];
                end
                s = get(handles.listbox_operations, 'string');
                set(handles.listbox_operations, 'string', ...
                    char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de falta de precisión en cálculo de intersecciones.'));
            end
        end
    end                                                        
elseif (abs(C) < e)
    C2_Z = dot(N2, [0 0 1]); C2_Y = dot(N2, [0 1 0]); C2_X = dot(N2, [1 0 0]);

    if (abs(C2_X) < e)
        %%%% Desarrollar
        s = get(handles.listbox_operations, 'string');
        set(handles.listbox_operations, 'string', ...
            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de colapso de plano paralelo a X.'));
    %    if points_set_j(index_4, 2) == points_set_j(index_1, 2) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
    %        if points_set_j(index_3, 2) == points_set_j(index_1, 2) && points_set_j(index_3, 3) == points_set_j(index_1, 3)
    %            points_set_j(index_3, :) = [0 y z];
    %        end
    %    else
    %        if points_set_j(index_4, 2) == points_set_j(index_1, 2) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
    %            points_set_j(index_4, :) = [0 y z];
    %        end
    %    end
    end
    if ((index_1 == 1 && index_2 == 2)||(index_1 == 4 && index_2 == 2)) && (abs(C2_Y) < e)
        if points_set_j(index_4, 2) == points_set_j(index_1, 2) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
            if y < points_set_j(index_3, 2)              
                points_set_j(index_3, :) = [points_set_j(index_3, 1) y points_set_j(index_3, 3)];
            end        
        else
            if y < points_set_j(index_4, 2)
                points_set_j(index_4, :) = [points_set_j(index_4, 1) y points_set_j(index_4, 3)];
            end
        end
    elseif ((index_1 == 2 && index_2 == 1)||(index_1 == 4 && index_2 == 1)) && (abs(C2_Z) < e)
        if points_set_j(index_4, 2) == points_set_j(index_1, 2) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
            if z < points_set_j(index_3, 3)
                points_set_j(index_3, :) = [points_set_j(index_3, 1) points_set_j(index_3, 2) z];
            end        
        else
            if z < points_set_j(index_4, 3)            
                points_set_j(index_4, :) = [points_set_j(index_4, 1) y points_set_j(index_4, 2) z];
            end
        end
    end
elseif (abs(C2) < e)
    C_Z = dot(N, [0 0 1]); C_Y = dot(N, [0 1 0]); C_X = dot(N, [1 0 0]);
    
    if (abs(C_X) < e)
        %%%% Desarrollar
        s = get(handles.listbox_operations, 'string');
        set(handles.listbox_operations, 'string', ...
            char(s, 'La reconstrucción puede ser incompleta.', 'Caso 14 de colapso de plano paralelo a X.'));
    %    if points_set_i(index_3, 2) == points_set_i(index_2, 2) && points_set_i(index_3, 3) == points_set_i(index_2, 3)
    %        if points_set_j(index_4, 2) == points_set_j(index_1, 2) && points_set_j(index_4, 3) == points_set_j(index_1, 3)
    %            points_set_i(index_4, :) = [0 y z];
    %        end
    %    else
    %        if points_set_j(index_3, 2) == points_set_j(index_1, 2) && points_set_j(index_3, 3) == points_set_j(index_1, 3)
    %            points_set_i(index_3, :) = [0 y z];
    %        end
    %    end
    end
    if ((index_1 == 1 && index_2 == 2)||(index_1 == 4 && index_2 == 2)) && (abs(C_Z) < e)
        if points_set_i(index_3, 2) == points_set_i(index_2, 2) && points_set_i(index_3, 3) == points_set_i(index_2, 3)
            if z < points_set_i(index_4, 3)              
                points_set_i(index_4, :) = [points_set_i(index_4, 1) points_set_i(index_4, 2) z];
            end        
        else
            if z < points_set_i(index_3, 3)            
                points_set_i(index_3, :) = [points_set_i(index_3, 1) points_set_i(index_3, 2) z];
            end
        end
    elseif ((index_1 == 2 && index_2 == 1)||(index_1 == 4 && index_2 == 1)) && (abs(C_Y) < e)
        if points_set_i(index_3, 2) == points_set_i(index_2, 2) && points_set_i(index_3, 3) == points_set_i(index_2, 3)
            if y < points_set_i(index_4, 2)
                points_set_i(index_4, :) = [points_set_i(index_4, 1) y points_set_i(index_4, 3)];
            end        
        else
            if y < points_set_i(index_3, 2)            
                points_set_i(index_3, :) = [points_set_i(index_3, 1) y points_set_i(index_3, 3)];
            end
        end
    end
end


%%%
%%% Calcula el orden de los nodos para trazar un polígono convexo
function ind = convex_hull(points, ps_i)
%global points_set_nonred;


e = 0.00005; %%%% Cambiar si se cambia precisión
%ps_i = points_set_nonred{i}.points;
N3 = cross(ps_i(1,:)-ps_i(2,:), ps_i(1,:)-ps_i(3,:));
if any(N3 ~= 0)
    C3 = dot(N3, [0 0 1]);
else
    C3 = 1;
end
N2 = cross(ps_i(1,:)-ps_i(2,:), ps_i(1,:)-ps_i(3,:));
if any(N2 ~= 0)
    C2 = dot(N2, [0 1 0]);
else
    C2 = 1;
end
N1 = cross(ps_i(1,:)-ps_i(2,:), ps_i(1,:)-ps_i(3,:));
if any(N1 ~= 0)
    C1 = dot(N1, [1 0 0]);
else
    C1 = 1;
end

%if (abs(C) < e && ~(abs(points(1, 1) - points(3, 1) < e) && abs(points(3, 1) - points(2, 1)) < e) && ...
%        ~(abs(points(1, 2) - points(3, 2) < e) && abs(points(3, 2) - points(2, 2)) < e))    
if abs(C3) < e
    %points_copy = sortrows(points(1:3,:), 3); %%%%ojo%%%
    %if ~(abs(points_copy(1, 1) - points_copy(2, 1)) < e  || abs(points_copy(1, 2) - points_copy(2, 2)) < e)
    %    punto1 = points_copy(1,:);
    %    punto2 = points_copy(2,:);
    %elseif ~(abs(points_copy(2, 3) - points_copy(3, 3)) < e)
    %    punto1 = points_copy(1,:);
    %    punto2 = points_copy(3,:);
    %else
    %    punto1 = points_copy(2,:);
    %    punto2 = points_copy(3,:);
    %end

    %if (abs(points_copy(1,3) - points_copy(2,3)) < e && abs(points_copy(2,3) - points_copy(3,3)) < e)
    %        theta = 0;
    %else
        theta = pi/2;
    %end

    %u = punto1 - punto2;    
    u = ps_i(3,:) - ps_i(2,:);
    R_u = [cos(theta)+u(1)^2*(1-cos(theta))           u(1)*u(2)*(1-cos(theta))-u(3)*sin(theta)    u(1)*u(3)*(1-cos(theta))+u(2)*sin(theta); ...
           u(2)*u(1)*(1-cos(theta))+u(3)*sin(theta)   cos(theta)+u(2)^2*(1-cos(theta))            u(2)*u(3)*(1-cos(theta))-u(1)*sin(theta); ...
           u(3)*u(1)*(1-cos(theta))-u(2)*sin(theta)   u(3)*u(2)*(1-cos(theta))+u(1)*sin(theta)    cos(theta)+u(3)^2*(1-cos(theta))];
    new_points = R_u*points';
    points = new_points';
%end
    
else
    %if ~(abs(ps_i(1, 3) - ps_i(2, 3)) < e || abs(ps_i(2, 3) - ps_i(3, 3)) < e || abs(ps_i(1, 3) - ps_i(3, 3)) < e)
    if ~(abs(C1) < e && abs(C2) < e) 
        syms z1 z2;
        %%%Encontrar el ángulo de rotación con respecto al eje X
        coef1 = sqrt(ps_i(1, 2)^2 + ps_i(1, 3)^2);
        if ps_i(1, 2) == 0
            angulo1 = 'NaN';    
        else
            angulo1 = atan(ps_i(1, 3)/ps_i(1, 2));
        end

        coef2 = sqrt(ps_i(2, 2)^2 + ps_i(2, 3)^2);
        if ps_i(2, 2) == 0
            angulo2 = 'NaN';
        else
            angulo2 = atan(ps_i(2, 3)/ps_i(2, 2));
        end
        if ~(abs(C2) < e && abs(C3) < e)
            theta = 0;
            if coef1 ~= 0 && coef2 ~= 0
                if (abs(C1) < e && abs(C3) < e)
                    theta = pi/2;        
                elseif ~(strcmp(angulo1, 'NaN')||strcmp(angulo2, 'NaN'))
                    sol = solve(['asin(z1/', num2str(coef1), ')-', num2str(angulo1), '-asin(z1/', num2str(coef2), ')+', num2str(angulo2)], z1);
                    if ~isempty(sol)
                        z = eval(sol);
                        theta = asin(min(abs(z(1)))/coef1) - angulo1;
                    end
                elseif strcmp(angulo1, 'NaN') && strcmp(angulo2, 'NaN')
                    sol = solve(['acos(z1/', num2str(ps_i(1, 3)), ')+', 'acos(z1/', num2str(ps_i(2, 3)), ')'], z1);
                    if ~isempty(sol)
                        z = eval(sol);
                        theta = asin(min(abs(z(1)))/coef2) - angulo2;            
                    end
                elseif strcmp(angulo1, 'NaN')
                    sol = solve(['acos(z1/', num2str(ps_i(1, 3)), ')', '-asin(z1/', num2str(coef2), ')+', num2str(angulo2)], z1);
                    if ~isempty(sol)
                        z = eval(sol);
                        theta = asin(min(abs(z(1)))/coef2) - angulo2;            
                    end
                elseif strcmp(angulo2, 'NaN')
                    sol = solve(['acos(z1/', num2str(ps_i(2, 3)), ')', '-asin(z1/', num2str(coef1), ')+', num2str(angulo1)], z1);
                    if ~isempty(sol)
                        z = eval(sol);
                        theta = asin(min(abs(z(1)))/coef1) - angulo1;            
                    end 
                end
            elseif any([coef1 coef2] ~= 0)
                theta = pi/2;
            end

            R_x = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            new_points = R_x*points';
            points = new_points'; 

            %syms x1 x2 x3;
            %sol = solve([num2str(points(1, 1)), '*x1+', num2str(points(1, 2)), '*x2+', num2str(points(1, 3)), '*x3=', num2str(points(1, 3))], ...
            %    [num2str(points(2, 1)), '*x1+', num2str(points(2, 2)), '*x2+', num2str(points(2, 3)), '*x3=', num2str(points(1, 3))], ...
            %    [num2str(points(3, 1)), '*x1+', num2str(points(3, 2)), '*x2+', num2str(points(3, 3)), '*x3=', num2str(points(1, 3))], x1, x2, x3);
            %sol1 = [sol.x1 sol.x2 sol.x3];
            %z = eval(sol1);
            %syms psi1 phi1 thet1;

            %sol = solve(['cos(psi1)*sin(thet1)*cos(phi1)+sin(phi1)*sin(phi1)=', num2str(z(1))], ...
            %            ['cos(psi1)*sin(thet1)*sin(phi1)-sin(psi1)*cos(phi1)=', num2str(z(2))], ...
            %            ['cos(psi1)*sin(thet1)=', num2str(z(3))], psi1, phi1, thet1);
            %sol2 = [sol.psi1 sol.phi1 sol.thet1];
            %z = eval(sol2);            
        else        
            theta = pi/2;      

            R_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
            new_points = R_y*points';
            points = new_points';
        end    
    else 
        if (abs(C2) < e && abs(C3) < e)
            theta = pi/2;       

            R_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
            new_points = R_y*points';
            points = new_points';               
        elseif (abs(C1) < e && abs(C3) < e) 
            theta = pi/2;

            R_x = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            new_points = R_x*points';
            points = new_points'; 
        else
            theta = pi/2;

            R_x = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            new_points = R_x*points';
            points = new_points';
        end
    end
end
ind = convhull(points(:,1), points(:,2));


% --------------------------------------------------------------------
function Restriction_nonnegativity_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD>
% hObject    handle to Restriction_nonnegativity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Dimension handles_surf;

if strcmp(get(hObject, 'Checked'), 'off')
    set(handles_surf(Dimension(1)+1), 'visible', 'on');
    set(handles_surf(Dimension(1)+2), 'visible', 'on');
    set(handles_surf(Dimension(1)+3), 'visible', 'on');
    set(hObject, 'Checked', 'on');
else
    set(handles_surf(Dimension(1)+1), 'visible', 'off');
    set(handles_surf(Dimension(1)+2), 'visible', 'off');
    set(handles_surf(Dimension(1)+3), 'visible', 'off');
    set(hObject, 'Checked', 'off');
end

% --------------------------------------------------------------------
function Uncut_planes_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Uncut_planes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global handles_surf handles_norm;

handles_surf = [];
handles_norm = [];
cla(handles.axes_simplex3D);
set(handles.listbox_operations, 'string', ...
    char('Descripción 3D del problema.', 'Planos sin corte.'));

set(handles.Mode_geo3D, 'Checked','off'); %Truco para correr función trace3D

trace3D(handles);    

set(handles.Mode_geo3D, 'Checked','on');

%set(handles.Uncut_planes, 'Enable','off');
set(handles.uipushtool1, 'Enable', 'on');
set(handles.uitoggletool6, 'enable', 'on');
set(hObject, 'Enable','off');
set(handles.pushbutton_asignall, 'Enable', 'on');


% --------------------------------------------------------------------
function Saveimage_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Saveimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global latex;

mkdir('Recursos');
[FileName,PathName,FilterIndex] = uiputfile({'*.eps','Archivo de imagen (*.eps)'},'Guardar imagen',...
          './Recursos/NuevoPoliedro1.eps'); %#ok<ASGLU,NASGU>
if FileName ~= 0
    pos1 = get(handles.panel_enhancement, 'Position');
    pos1_aux = pos1;
    pos2 = get(handles.LPApp, 'Position');
    pos1(1) = pos2(3);
    pos1(2) = pos2(4);
    set(handles.panel_enhancement, 'Position', pos1);
    saveas(gcf, ['./Recursos/', FileName], 'eps');
    
    suffix = num2str(now);
    latex = '';
    latex = sprintf('%s \r La representación gráfica 3D ', latex);
    latex = [latex, '(\ref{fig:poli',suffix, '}) es el problema equivalente en forma de desigualdades.'];
    latex = sprintf('%s \r', latex);
    latex = [latex, '\begin{figure}[htbp]'];
    latex = sprintf('%s \r', latex);
    latex = [latex, '\centering'];
    latex = sprintf('%s \r', latex);
    latex = [latex, '\includegraphics [width = 0.8\textwidth] {', ['../Recursos/', FileName], '}']; 
    latex = sprintf('%s \r', latex);
    latex = [latex, '\caption{Representación gráfica 3D (X1, X2, X3).} \label{fig:poli', suffix, '}']; 
    latex = sprintf('%s \r', latex);
    latex = [latex, '\end{figure}'];
    latex = sprintf('%s \r', latex);
    handles.latex = [handles.latex, latex];
    guidata(handles.output, handles);
    set(handles.panel_enhancement, 'Position', pos1_aux);
end


% --------------------------------------------------------------------
function Simplex_multiplier_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Simplex_multiplier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Order_initial Tableau Order_current Matrix_problem;

%Order_initial = LPApphandle.Order_initial(1, 1:(dim(2)-1));
Dimension = size (Tableau);
%Basic_vector_initial = Order_initial < Dimension(1);
%bk = LPApphandle.gui_Matrix_problem(1:(dim(1)-1), end);
%yk = LPApphandle.gui_tableau(1:(dim(1)-1), end);
%BInv = Tableau(1:(Dimension(1)-1), Basic_vector_initial==1);
[Y, I] = sort(Order_initial); %#ok<ASGLU>
Tableau_sorted = Tableau(:, I);
[Y, J] = sort(Order_current); %#ok<ASGLU>
Matrix_sorted = Matrix_problem(:, J);
c = Matrix_sorted(end,:);
BInv = Tableau_sorted(1:(Dimension(1)-1), 1:(Dimension(1)-1));

All_display = zeros(Dimension(1), Dimension(2)+2);
% se actualiza la tabla de la interfaz
sdisplay = get(handles.table_simplexdisplay, 'data');
T = cell2mat(sdisplay(1:Dimension(1), 1:Dimension(2)+1));
dim = size(T);

All_display(:, 1:dim(2)) = cell2mat(sdisplay(1:dim(1), 1:dim(2)));
lambda = c(1:(Dimension(1)-1))*BInv;
%All_display(1:(Dimension(1)-1), end) = lambda(Order_initial(Order_initial < Dimension(1)));
All_display(1:(Dimension(1)-1), end) = lambda;

spreadsheet = cell(100,100);
spreadsheet(1:Dimension(1), 1:Dimension(2)+2) = num2cell(All_display);
set(handles.table_simplexdisplay, 'data', spreadsheet);

colName = cell(Dimension(2)+2, 1);
for i = 1:(Dimension(2)-1)
    colName(i) = cellstr(strcat('X',num2str(i)));
end
colName(Dimension(2)) = cellstr('Yi0'); 
if handles.Method == 1 || handles.Submethod == 1
    colName(Dimension(2)+1) = cellstr('Yi0/Yij');
    colName(Dimension(2)+2) = cellstr('lambdai');
else
    colName(Dimension(2)+1) = cellstr('lambdai');
end
set(handles.table_simplexdisplay, 'columnname', colName);

%setrowheaders(handles, headers);

if handles.Method == 1 || handles.Submethod == 1
    colFormat=cell(1,Dimension(2)+2);    
else
    colFormat=cell(1,Dimension(2)+1);
end
% se especifica el formato de salida de los datos
colFormat(1,:) = cellstr('rat');

set(handles.table_simplexdisplay, 'columnformat', colFormat);


% --------------------------------------------------------------------
function Branch_and_bound_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Branch_and_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Newmethod = 4;
init_method(handles);
