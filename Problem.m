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

if handles.LPApphandle.Newmethod ~= 3
    setheading(handles, char('f', 'X', 'Yi0', 'Z', 'lambdai', 'Ui'));
else
    setheading(handles, char('O', 'D', 'Ci', 'Rj'));
end
if (handles.LPApphandle.istwophases == 1 || handles.LPApphandle.Newmethod == 4) && handles.LPApphandle.isonlylecture == 1
    set(handles.pushbutton_cancel, 'Enable', 'off');
    set(hObject, 'WindowStyle', 'modal');
    if handles.LPApphandle.Newmethod == 4
        set(handles.uitoggletool1, 'Enable', 'off');
        set(handles.uitoggletool1, 'State', 'on');
        set(handles.uipushtool1, 'TooltipString', 'Eliminar fila');
        set(handles.uipushtool2, 'Enable', 'off');
    end
else
    set(handles.pushbutton_cancel, 'Enable', 'on');
    set(hObject, 'WindowStyle', 'normal');
end
% UIWAIT makes problem wait for user response (see UIRESUME)
% uiwait(handles.Problem);

function setheading(handles, headers)
global Dimension;

% se obtienen las dimensiones del problema
Dimension = size(handles.LPApphandle.gui_Matrix_problem);

%All_display = zeros(Dimension(1), Dimension(2));

% se rotulan las filas y columnas de las tablas
colName = cell(Dimension(2), 1);
if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
        handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 4
    j = 0;
    if handles.LPApphandle.Newmethod == 4 && handles.LPApphandle.isonlylecture == 1
        Dim = size(handles.LPApphandle.Orig_Matrix_problem);
        j = (Dim(2)-1)- (Dimension(2) -1);
    end
    for i = (1+j):Dimension(2)-1+j
        colName(i-j) = cellstr(strcat(headers(2,:),num2str(i)));
    end
    
    colName(Dimension(2)) = cellstr(headers(3,:));
elseif handles.LPApphandle.Newmethod == 5
    for i = 1:Dimension(2)-2
        colName(i) = cellstr(strcat(headers(2,:),num2str(i)));
    end
    
    colName(Dimension(2)-1) = cellstr(headers(3,:));
    colName(Dimension(2)) = cellstr(headers(5,:));
end   

set(handles.table_problem, 'columnname', colName);

rowName=cell(1,Dimension(1));
if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
        handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 5
    for i = 1:Dimension(1)-1
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end
    
    rowName(Dimension(1)) = cellstr(headers(4,:));
elseif handles.LPApphandle.Newmethod == 4
    for i = 1:Dimension(1)-2
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end
    
    rowName(Dimension(1)-1) = cellstr(headers(4,:));
    rowName(Dimension(1)) = cellstr(headers(6,:));
end

set(handles.table_problem, 'rowname', rowName);

% se establece el formato de las columnas
colFormat=cell(1,Dimension(2));
for i = 1:Dimension(2)
    colFormat(1,i) = cellstr('rat');
end
set(handles.table_problem, 'columnformat', colFormat);

if handles.LPApphandle.istwophases ~= 1
    % se configuran las columnas como editables
    colEdit = ones(1,Dimension(2));
    set(handles.table_problem, 'columneditable', (colEdit == 1));
end
% se despliega la tabla para el ingreso de datos
All_display(1:Dimension(1), 1:Dimension(2)) = handles.LPApphandle.gui_Matrix_problem;
set(handles.table_problem, 'data', All_display);

% --- Outputs from this function are returned to the command line.
function varargout = Problem_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles, 'output')
    varargout{1} = handles.output;
else
    close(hObject);
end

% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles) %#ok<INUSL>
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem Order_initial latex Dimension;
% Matrix_problem    guarda la especificación del problema
% Orden_initial     registra el orden de las columnas según la matriz
% identidad

% recupera la especificación dada por el usuario
All_display = get(handles.table_problem, 'data');
Matrix_problem = All_display;

Dimension = size(Matrix_problem);

if ~(handles.LPApphandle.istwophases == 1 && handles.LPApphandle.isonlylecture == 1) %%%%OJO
    if handles.LPApphandle.isonlylecture ~= 1
        handles.LPApphandle.latexproblem = '';
        handles.LPApphandle.latexIIphasesproblem = '';
        handles.LPApphandle.latexfile = '';
        handles.LPApphandle.gui_Matrix_problem = []; % se comparte la especificación del problema   
        handles.LPApphandle.First_Matrix_problem = [];
        handles.LPApphandle.istwophases = 0;    
        handles.LPApphandle.whatphase = 1;
        handles.LPApphandle.isinit_secondphase = 0;
        handles.LPApphandle.Orig_Matrix_problem = [];
        handles.LPApphandle.maxcanon_vector = 0;
        handles.LPApphandle.nodelist = cell(1);
        handles.LPApphandle.nodelist{1} = [];
        handles.LPApphandle.ispenalties = 0;
        handles.LPApphandle.ZOptimo = Inf;
        handles.LPApphandle.XOptimo = [];
        handles.LPApphandle.Level = 0;
        handles.LPApphandle.Travel = 1;
        handles.Interactive = 1;
        set(handles.LPApphandle.pushbutton_next, 'String', 'Mejorar solucion');
    end
    % se verifica que la especificación este en forma correcta para el método
    % correspondiente
    resp = iscorrectform(handles);
    if ~resp
        return;
    elseif resp == 2 %% Dual
        handles.LPApphandle.Submethod = 2;
        
        if handles.LPApphandle.isonlylecture ~= 1
            handles.LPApphandle.Orig_Matrix_problem = Matrix_problem;
            set(handles.LPApphandle.text_selectvar2, 'string', 'Recorrido');
            set(handles.LPApphandle.pushbutton_next, 'String', 'Mejorar solucion');
            % se construye el arreglo de índices de recorridos
            Var = char('LIFO-a', 'LIFO-d', 'FIFO-a', 'FIFO-d');
            set(handles.LPApphandle.popupmenu_selectvar2, 'string', Var);
            set(handles.LPApphandle.popupmenu_selectvar2, 'value', 1);
            set(handles.LPApphandle.popupmenu_selectvar2, 'Enable', 'on');
        end
        
        aux =  zeros(Dimension(1)-1, Dimension(2)+Dimension(1)-2);
        aux(1:Dimension(1)-2,:) = [Matrix_problem(1:Dimension(1)-2, 1:Dimension(2)-1), eye(Dimension(1)-2), ...
            Matrix_problem(1:Dimension(1)-2, Dimension(2))];
        aux(Dimension(1)-1,1:Dimension(2)-1) = Matrix_problem(Dimension(1)-1, ...
            1:Dimension(2)-1);
        Matrix_problem = aux;        
        Dimension = size(Matrix_problem);
        Order_initial = zeros(1, Dimension(2));
        
        iscanonicform(handles);
        if handles.LPApphandle.isonlylecture ~= 1
            handles.LPApphandle.First_Matrix_problem = Matrix_problem;
            handles.LPApphandle.First_Order_initial = Order_initial;
        end
        guidata(handles.LPApphandle.output, handles.LPApphandle);
    elseif resp == 3 %% Primal
        handles.LPApphandle.Submethod = 1;

        if handles.LPApphandle.isonlylecture ~= 1
            handles.LPApphandle.Orig_Matrix_problem = Matrix_problem;
            set(handles.LPApphandle.text_selectvar2, 'string', 'Recorrido');
            set(handles.LPApphandle.pushbutton_next, 'String', 'Mejorar solucion');
            % se construye el arreglo de índices de recorridos
            Var = char('LIFO-a', 'LIFO-d', 'FIFO-a', 'FIFO-d');
            set(handles.LPApphandle.popupmenu_selectvar2, 'string', Var);
            set(handles.LPApphandle.popupmenu_selectvar2, 'value', 1);
            set(handles.LPApphandle.popupmenu_selectvar2, 'Enable', 'on');
        end
        
        aux =  zeros(Dimension(1)-1, Dimension(2)+Dimension(1)-2);
        aux(1:Dimension(1)-2,:) = [Matrix_problem(1:Dimension(1)-2, 1:Dimension(2)-1), eye(Dimension(1)-2), ...
            Matrix_problem(1:Dimension(1)-2, Dimension(2))];
        aux(Dimension(1)-1,1:Dimension(2)-1) = Matrix_problem(Dimension(1)-1, ...
            1:Dimension(2)-1);
        Matrix_problem = aux;
        Dimension = size(Matrix_problem);
        Order_initial = zeros(1, Dimension(2));
        
        iscanonicform(handles);
        if handles.LPApphandle.isonlylecture ~= 1
            handles.LPApphandle.First_Matrix_problem = Matrix_problem;
            handles.LPApphandle.First_Order_initial = Order_initial;
        end
        guidata(handles.LPApphandle.output, handles.LPApphandle);
    elseif resp == 4 %% Primal penalizaciones
        handles.LPApphandle.Submethod = 1;
        
        if handles.LPApphandle.isonlylecture ~= 1
            handles.LPApphandle.Orig_Matrix_problem = Matrix_problem;
            set(handles.LPApphandle.text_selectvar2, 'string', 'Recorrido');
            set(handles.LPApphandle.pushbutton_next, 'String', 'Mejorar solucion');
            % se construye el arreglo de índices de recorridos
            Var = char('LIFO-a', 'LIFO-d', 'FIFO-a', 'FIFO-d');
            set(handles.LPApphandle.popupmenu_selectvar2, 'string', Var);
            set(handles.LPApphandle.popupmenu_selectvar2, 'value', 1);
            set(handles.LPApphandle.popupmenu_selectvar2, 'Enable', 'on');
        end
        
        aux =  zeros(Dimension(1)-1, Dimension(2)+Dimension(1)-2);
        aux(1:Dimension(1)-2,:) = [Matrix_problem(1:Dimension(1)-2, 1:Dimension(2)-1), eye(Dimension(1)-2), ...
            Matrix_problem(1:Dimension(1)-2, Dimension(2))];
        aux(Dimension(1)-1,1:Dimension(2)-1) = Matrix_problem(Dimension(1)-1, ...
            1:Dimension(2)-1);
        Matrix_problem = aux;
        Dimension = size(Matrix_problem);
        Order_initial = zeros(1, Dimension(2));
        
        % Convertir a forma estándar
        X0 = Matrix_problem(1:(Dimension(1)-1), end);
        for i = 1:Dimension(1)-1
            if X0(i) < 0
                Matrix_problem(i, :)= -Matrix_problem(i, :);
            end
        end
        
        iscanonicform(handles);
        
        identidad = eye(Dimension(1)-1);
        
        handles.LPApphandle.ispenalties = 1;
        %handles.LPApphandle.Orig_Matrix_problem = Matrix_problem;
        daux = size(Order_initial(Order_initial>0));
        handles.LPApphandle.maxcanon_vector = daux(2);
        aux = Matrix_problem(1:Dimension(1)-1, Dimension(2));
        %Matrix_problem(Dimension(1), :) = zeros(1, Dimension(2));
        j=0;
        maxi = max(abs(Matrix_problem(Dimension(1), :)));
        for i = (1:Dimension(1)-1)
            ind = find(Order_initial==i, 1);
            if isempty(ind)
                Matrix_problem(1:Dimension(1)-1, Dimension(2) + j) = identidad(:, i);
                Matrix_problem(Dimension(1), Dimension(2) + j) = maxi(1)*100;
                Order_initial(Dimension(2) + j) = i;                
                j=j+1;
            end
        end
        Order_initial(Order_initial == 0) = Dimension(1):Dimension(2) + j - 1;
        Order_initial(Dimension(2) + j) = Dimension(2) + j;
        Matrix_problem(1:Dimension(1)-1, Dimension(2) + j) = aux;
        if handles.LPApphandle.isonlylecture ~= 1
            handles.LPApphandle.First_Matrix_problem = Matrix_problem;
            handles.LPApphandle.First_Order_initial = Order_initial;
        end
        guidata(handles.LPApphandle.output, handles.LPApphandle);
    elseif resp == 5   %% Dos fases     
        identidad = eye(Dimension(1)-1);
        
        handles.LPApphandle.istwophases = 1;
        handles.LPApphandle.whatphase = 1;
        handles.LPApphandle.isinit_secondphase = 0;
        handles.LPApphandle.Orig_Matrix_problem = Matrix_problem;
        daux = size(Order_initial(Order_initial>0));
        handles.LPApphandle.maxcanon_vector = daux(2);
        aux = Matrix_problem(1:Dimension(1)-1, Dimension(2));
        Matrix_problem(Dimension(1), :) = zeros(1, Dimension(2));
        j=0;       
        %for i = (Order+1:Dim(1)-1)
        for i = (1:Dimension(1)-1)
            ind = find(Order_initial==i, 1);
            if isempty(ind)
                Matrix_problem(1:Dimension(1)-1, Dimension(2) + j) = identidad(:, i);
                Matrix_problem(Dimension(1), Dimension(2) + j) = 1;
                Order_initial(Dimension(2) + j) = i;                
                j=j+1;
            end
        end
        Order_initial(Order_initial == 0) = Dimension(1):Dimension(2) + j - 1;
        Order_initial(Dimension(2) + j) = Dimension(2) + j;
        Matrix_problem(1:Dimension(1)-1, Dimension(2) + j) = aux;
        %Dimension = size(Matrix_problem);
        handles.LPApphandle.First_Matrix_problem = Matrix_problem;
        handles.LPApphandle.First_Order_initial = Order_initial;
        guidata(handles.LPApphandle.output, handles.LPApphandle);
    elseif resp == 6    %% Penalizaciones            
        identidad = eye(Dimension(1)-1);
        
        handles.LPApphandle.ispenalties = 1;
        handles.LPApphandle.Orig_Matrix_problem = Matrix_problem;
        daux = size(Order_initial(Order_initial>0));
        handles.LPApphandle.maxcanon_vector = daux(2);
        aux = Matrix_problem(1:Dimension(1)-1, Dimension(2));
        %Matrix_problem(Dimension(1), :) = zeros(1, Dimension(2));
        j=0;       
        maxi = max(Matrix_problem(Dimension(1), :));
        for i = (1:Dimension(1)-1)
            ind = find(Order_initial==i, 1);
            if isempty(ind)
                Matrix_problem(1:Dimension(1)-1, Dimension(2) + j) = identidad(:, i);
                Matrix_problem(Dimension(1), Dimension(2) + j) = maxi(1)*100;
                Order_initial(Dimension(2) + j) = i;                
                j=j+1;
            end
        end
        Order_initial(Order_initial == 0) = Dimension(1):Dimension(2) + j - 1;
        Order_initial(Dimension(2) + j) = Dimension(2) + j;
        Matrix_problem(1:Dimension(1)-1, Dimension(2) + j) = aux;        
        guidata(handles.LPApphandle.output, handles.LPApphandle);
    end

    handles.LPApphandle.Method = handles.LPApphandle.Newmethod;
    
    set(handles.LPApphandle.Restriction_nonnegativity, 'checked', 'off');
    set(handles.LPApphandle.Restriction_nonnegativity, 'Enable', 'off');
    set(handles.LPApphandle.Uncut_planes, 'Enable', 'off');

    if handles.LPApphandle.Method ~= 3
        latex = '';
        if handles.LPApphandle.isonlylecture ~= 1
            generar_latexspecification('\section{Especificación del problema}', handles);
            handles.LPApphandle.latex = latex;   
            %l = handles.LPApphandle.latex;
            guidata(handles.LPApphandle.output, handles.LPApphandle);
            %l = handles.LPApphandle.latex;
             %AGREGADO(27/12/2016)
            if handles.LPApphandle.istwophases == 1
                latex = '\section{Primera fase del Método Simplex}';
                generar_latexspecification('\subsection{Especificación del problema asociado}', handles);
                handles.LPApphandle.latex = [handles.LPApphandle.latex, latex];
                generar_latexsimplexbegin('\subsection{Desarrollo del método simplex}');
                handles.LPApphandle.latex = [handles.LPApphandle.latex, latex];
                handles.LPApphandle.latexIIphasesproblem = handles.LPApphandle.latex;
                guidata(handles.LPApphandle.output, handles.LPApphandle);        
            elseif handles.LPApphandle.Method == 4
                latex = '';
                generar_latexspecification(['\section{Desarrollo del Método de ramificación y poda}', sprintf('%s \r', ''), ...
                    '\subsection{Especificación del subproblema en el nivel ', num2str(handles.LPApphandle.Level) ,'}'], handles);   
                handles.LPApphandle.latex = [handles.LPApphandle.latex, latex]; 
                generar_latexsimplexbegin('\subsection{Desarrollo del Método Simplex}');                
                handles.LPApphandle.latex = [handles.LPApphandle.latex, latex]; 
                handles.LPApphandle.latexproblem = handles.LPApphandle.latex;
                %l = handles.LPApphandle.latex;                
                guidata(handles.LPApphandle.output, handles.LPApphandle); 
            else
                generar_latexsimplexbegin('\section{Desarrollo del Método Simplex}');
                handles.LPApphandle.latex = [handles.LPApphandle.latex, latex];  
                %l = handles.LPApphandle.latex;
                handles.LPApphandle.latexproblem = handles.LPApphandle.latex;
                guidata(handles.LPApphandle.output, handles.LPApphandle); 
            end
        elseif handles.LPApphandle.Method == 4
            generar_latexspecification(['\subsection{Especificación del subproblema en el nivel ', ...
                num2str(handles.LPApphandle.Level) ,'}'], handles);
            handles.LPApphandle.latex = [handles.LPApphandle.latex, latex];
            generar_latexsimplexbegin('\subsection{Desarrollo del Método Simplex}');
            handles.LPApphandle.latex = [handles.LPApphandle.latex, latex]; 
            %l = handles.LPApphandle.latex;
            guidata(handles.LPApphandle.output, handles.LPApphandle);
        end
    %AGREGADO(27/12/2016)
    else
        generar_latextransportespecification()
        handles.LPApphandle.latex = latex;
        handles.LPApphandle.latexproblem = latex;
        guidata(handles.LPApphandle.output, handles.LPApphandle);
    end

    handles.LPApphandle.Order_initial = Order_initial; % se comparte el orden inicial
    latex = '';
    handles.LPApphandle.isonlylecture = 0;
    guidata(handles.LPApphandle.output, handles.LPApphandle);
else
    handles.LPApphandle.isonlylecture = 0;
    guidata(handles.LPApphandle.output, handles.LPApphandle);
end
%guidata(handles.LPApphandle.output, handles.LPApphandle); 
handles.LPApphandle.setProblem(Matrix_problem, handles.LPApphandle); % se establece el problema actual
if strcmp(get(handles.LPApphandle.Mode_geo3D, 'Checked'), 'on')
    set(handles.LPApphandle.Mode_geo3D, 'Checked', 'off');
    set(handles.LPApphandle.axes_simplex3D, 'visible', 'off');
    set(handles.LPApphandle.table_simplexdisplay, 'visible', 'on');
end
delete(handles.Problem);

% ---Funcion utilitaria: Verifica si se ha ingresado la formulación del
% problema de manera correcta
function iscorrect = iscorrectform(handles)
global Matrix_problem Order_initial Dimension;

Dimension = size(Matrix_problem);
Order_initial = zeros(1, Dimension(2));

if handles.LPApphandle.Newmethod == 3
    oferta = Matrix_problem(1:(Dimension(1)-1), end);
    demanda = Matrix_problem(end, 1:(Dimension(2)-1));
    if sum(oferta)~= sum(demanda) && ~isdegeneratedsolution(handles)
        if sum(oferta) - sum(demanda) > 0
            resp = questdlg(char('El problema no está balanceado. ¿Desea agregar un destino?'),'Problema no balanceado', ...
                'De acuerdo','No','Cancelar','Cancelar');
            if strcmp(resp, 'De acuerdo')
                uipushtool2_ClickedCallback(handles.uipushtool2, handles, handles);
                All_display = get(handles.table_problem, 'data');
                % se obtienen las dimensiones del problema
                Dimension = size(All_display);
                All_display(end, Dimension(2)-1) = sum(oferta) - sum(demanda);
                set(handles.table_problem, 'data', All_display);
            end
        else
            resp = questdlg(char('El problema no está balanceado. ¿Desea agregar un origen?'),'Problema no balanceado', ...
                'De acuerdo','No','Cancelar','Cancelar');
            if strcmp(resp, 'De acuerdo')
                uipushtool1_ClickedCallback(handles.uipushtool1, handles, handles);
                All_display = get(handles.table_problem, 'data');
                % se obtienen las dimensiones del problema
                Dimension = size(All_display);
                All_display(Dimension(1)- 1, end) = sum(demanda) - sum(oferta);
                set(handles.table_problem, 'data', All_display);
            end
        end
                    
        iscorrect = 0;
        return;
    end
    if all(Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1) == 0)
        errordlg(char('El problema no tiene una función objetivo definida.'),char('Método de transporte'),'modal');
        iscorrect = 0;
        return;
    end
    if ~isfactiblesolution(handles)
        iscorrect = 0;
        return;
    end  
    if isdegeneratedsolution(handles)        
        errordlg(char('La oferta (demanda) de algún origen (destino) es cero.'),char('Origen o destino inútil'),'modal');
        iscorrect = 0;
        return;    
    end
elseif handles.LPApphandle.Newmethod == 1 
    if rank(full(Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1))) ~= Dimension(1)-1
       errordlg('El problema tiene restricciones redundantes.',char('Método Simplex'),'modal');
       iscorrect = 0;
       return;
    end
    
    resp = iscanonicform(handles);
    if resp == 0
        iscorrect = 0;
        return;
    elseif resp == 5
        iscorrect = 5;
        
        resp = isfactiblesolution(handles);
        if ~resp
            iscorrect = 0;
            return;    
        end
        if isdegeneratedsolution(handles)        
            msgbox(char('La solución inicial es degenerada.'),char('Método simplex'),'modal');
        end
        return;
    elseif resp == 6
        iscorrect = 6;
        
        resp = isfactiblesolution(handles);
        if ~resp
            iscorrect = 0;
            return;    
        end
        if isdegeneratedsolution(handles)        
            msgbox(char('La solución inicial es degenerada.'),char('Método simplex'),'modal');
        end
        return; 
    end
    
    if isdegeneratedsolution(handles)        
        msgbox(char('La solución inicial es degenerada.'),char('Método simplex'),'modal');
    end
elseif handles.LPApphandle.Newmethod == 2
    if rank(full(Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1))) ~= Dimension(1)-1
       errordlg('El problema tiene restricciones redundantes.',char('Método Simplex'),'modal');
       iscorrect = 0;
       return;
    end            
    
    resp = iscanonicform(handles);
    if resp == 0
        iscorrect = 0;
        return;
    elseif resp == 5
        iscorrect = 5;
        
        resp = isfactiblesolution(handles);
        if ~resp
            iscorrect = 0;
            return;    
        end
        if isdegeneratedsolution(handles)        
            msgbox(char('La solución inicial es degenerada.'),char('Método simplex'),'modal');
        end
        return;
    elseif resp == 6
        iscorrect = 6;
        
        resp = isfactiblesolution(handles);
        if ~resp
            iscorrect = 0;
            return;    
        end
        if isdegeneratedsolution(handles)        
            msgbox(char('La solución inicial es degenerada.'),char('Método simplex'),'modal');
        end
        return;
    end
    
    if isdegeneratedsolution(handles)        
        msgbox(char('La solución inicial es degenerada.'),char('Método simplex'),'modal');
    end
    
elseif handles.LPApphandle.Newmethod == 4
    if handles.LPApphandle.isonlylecture ~= 1
        %standarproblem = [Matrix_problem(1:Dimension(1)-2, 1:Dimension(2)-1), eye(1:Dimension(1)-2)];
        if rank(full(Matrix_problem(1:Dimension(1)-2, 1:Dimension(2)-1))) ~= Dimension(1)-2
           errordlg('El problema tiene restricciones redundantes.',char('Método de ramificación y poda'),'modal');
           iscorrect = 0;
           return;
        end
    end
    resp = isfactiblesolution(handles);
    if ~resp
        iscorrect = 0;
        return;
    elseif resp == 2
        iscorrect = 2;
        
        if isdegeneratedsolution(handles)        
            msgbox(char('La solución inicial es degenerada.'),char('Método de ramificación y poda'),'modal');
        end
        return;        
    elseif resp == 3
        iscorrect = 3;
        
        if isdegeneratedsolution(handles)        
            msgbox(char('La solución inicial es degenerada.'),char('Método de ramificación y poda'),'modal');
        end
        return;
    elseif resp == 4
        iscorrect = 4;
        if isdegeneratedsolution(handles)        
            msgbox(char('La solución inicial es degenerada.'),char('Método de ramificación y poda'),'modal');
        end
        return;
    end     
end

iscorrect = 1;

% -------
function response = isfactiblesolution(handles)
global Matrix_problem Dimension Order_initial;

% se verifica que la solución inicial sea factible
if handles.LPApphandle.Newmethod == 1    
    X0 = Matrix_problem(1:(Dimension(1)-1), end);
    if any(X0 < 0)
        errordlg(char('La solución inicial no es primal factible.'),char('Método Simplex Primal no aplicable'),'modal');
        response = 0;
        return;
    end
elseif handles.LPApphandle.Newmethod == 2
    % se convierte la última en términos de Rj =cj - zj
    [Y, I] = sort(Order_initial); %#ok<ASGLU>
    Tableau_sorted = Matrix_problem(:, I);
    c = Tableau_sorted(end,:);
    z = Tableau_sorted(end, 1:(Dimension(1)-1))*Tableau_sorted(1:(Dimension(1)-1),:);

    %Rj = c - z; 
    Rj = eval(sym(str2num(rat(c)))) - eval(sym(str2num(rat(z)))); %#ok<ST2NM>
    if any(Rj < 0)
        response = 0;
        errordlg(char('La solución inicial no es dual factible.'),char('Método Simplex Dual no aplicable'),'modal');
        return;
    end
elseif handles.LPApphandle.Newmethod == 3
    oferta = Matrix_problem(1:(Dimension(1)-1), end);
    demanda = Matrix_problem(end, 1:(Dimension(2)-1));
    if any(oferta < 0) || any(demanda < 0)        
        response = 0;
        errordlg(char('La solución inicial no será factible.'),char('Método Simplex no aplicable'),'modal');
        return;
    end
elseif handles.LPApphandle.Newmethod == 4
    X0 = Matrix_problem(1:(Dimension(1)-2), end);
    % se convierte la última en términos de Rj =cj - zj
    c_aux =  [Matrix_problem(Dimension(1)-1, 1:Dimension(2)-1), zeros(1,Dimension(1)-2)];
    M_aux =  [Matrix_problem(1:Dimension(1)-2, 1:Dimension(2)-1), eye(Dimension(1)-2)];
    Order_aux = [Dimension(1)-1:Dimension(2)-1+Dimension(1)-2, 1:Dimension(1)-2];
    
    [Y, I] = sort(Order_aux); %#ok<ASGLU>
    Tableau_sorted = M_aux(:, I);
    c = c_aux(:,I);
    z = c(1, 1:Dimension(1)-2)*Tableau_sorted(1:(Dimension(1)-2),1:Dimension(2)-1+Dimension(1)-2);
    %Rj = c - z;
    Rj = eval(sym(str2num(rat(c)))) - eval(sym(str2num(rat(z)))); %#ok<ST2NM>
    
    Ui = Matrix_problem(end, 1:Dimension(2)-1);
    if all(Ui >= 0)
        if any(X0 < 0) && all(Rj >= 0)
            response = 2;

            msgbox(char('La solución inicial es dual factible. Se aplicará el método simplex dual en los nodos.'), ...
                char('Método de ramificación y cota'),'modal');        
            return;
        elseif all(X0 >= 0)
            response = 3;

            msgbox(char('La solución inicial es primal factible. Se aplicará el método simplex primal en los nodos.'), ...
                char('Método de ramificación y cota'),'modal');                
            return;
        elseif any(X0 < 0)
            resp = questdlg(char('La solución inicial no es primal ni dual factible.',...
                ' ¿Desea transformar el problema y aplicar el método de penalizaciones?'), ...
                char('Método de ramificación y cota'), 'De acuerdo','No','Cancelar','Cancelar');
            if ~isempty(resp)
                if strcmp(resp, 'De acuerdo')
                    response = 4;                     
                elseif strcmp(resp, 'No') || strcmp(resp, 'Cancelar')
                    response = 0;
                end
            else
                response = 0;
            end
            %MODIFICADO(27/12/2016)
            return;
        end
    else
        response = 0;
        
        msgbox(char('El problema no está en forma estándar.'), ...
                char('Método de ramificación y cota'),'modal');
        return;
    end
end
response = 1;

% ---------

function response = isdegeneratedsolution(handles)
global Matrix_problem Dimension;

% se verifica que la solución no sea degenerada
if handles.LPApphandle.Newmethod ~= 3
    if handles.LPApphandle.Newmethod == 2 || handles.LPApphandle.Newmethod == 1
        X0 = Matrix_problem(1:(Dimension(1)-1), end);
    elseif handles.LPApphandle.Newmethod == 4
        X0 = Matrix_problem(1:(Dimension(1)-2), end);
    elseif handles.LPApphandle.Newmethod == 5
        X0 = Matrix_problem(1:(Dimension(1)-1), Dimension(2)-1);
    end
    
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
function response = iscanonicform(handles)
global Matrix_problem Order_initial Dimension;


n = Dimension(2)-1;
m = Dimension(1)-1;

% se construye la matriz identidad
identidad = eye(m);
% se verifica que la matriz identidad este inmersa en la especifición del
% problema
response = zeros(1,m);
for i =1:m
    for j = 1:n            
        if identidad(:, i) == Matrix_problem(1:m, j)
            Order_initial(j) = i;
            response(i) = 1;
            break;        
        end
    end

end
if ~all(response(:)==1)
    if handles.LPApphandle.Newmethod ~= 4
        %MODIFICADO(27/12/2016)
        resp = questdlg(char('La matriz no está en forma canónica. ¿Cuál método desea aplicar?'),char('Método simplex'), ...
            'Dos fases','Penalizaciones','Cancelar','Cancelar');
        if ~isempty(resp)
            if strcmp(resp, 'Dos fases')
                response = 5;            
            elseif strcmp(resp, 'Penalizaciones')
                response = 6;
            elseif strcmp(resp, 'Cancelar')
                response = 0;
            end
        else
            response = 0;
        end
        %MODIFICADO(27/12/2016)
    end
    return;
end
Order_initial(Order_initial == 0) = m+1:n+1; % se etiquetan las demás columnas
response = 1;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.Problem);

function generar_latexspecification(title, handles)
global Matrix_problem latex;

dim = size(Matrix_problem);
All_display = get(handles.table_problem, 'data');
% se obtienen las dimensiones del problema
Dim_sub = size(All_display);
k = 0; l=0;
integral = ''; sep_sistema = ' = & ';
if handles.LPApphandle.Newmethod == 4 && handles.LPApphandle.isonlylecture == 1
    Dim_aux = size(handles.LPApphandle.Orig_Matrix_problem);
    k = (Dim_aux(2) - 1)- (Dim_sub(2) - 1);
elseif handles.LPApphandle.Newmethod == 4 && handles.LPApphandle.isonlylecture ~= 1 && ...
        isempty(strfind(title, 'nivel'))
    if handles.LPApphandle.ispenalties == 1
        l = dim(1)-1+(dim(1)-1-handles.LPApphandle.maxcanon_vector);
    else
        l = dim(1)-1;
    end
    integral = 'x_j \in \mathbb{Z}';
    sep_sistema = '\le &';
elseif (handles.LPApphandle.istwophases == 1 && isempty(strfind(title, 'asociado')) ...
        || (handles.LPApphandle.ispenalties == 1 && handles.LPApphandle.Newmethod ~=4))
    l = dim(1) -1 - handles.LPApphandle.maxcanon_vector;
end

latex = sprintf('%s \r', latex);
latex = [latex, title];
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
for j = 1:dim(2) - 1- l
    if (handles.LPApphandle.istwophases == 1 && isempty(strfind(title, 'asociado')) ...
            || (handles.LPApphandle.ispenalties == 1 && handles.LPApphandle.Newmethod ~=4))
        if handles.LPApphandle.Orig_Matrix_problem(dim(1), j) == 1
            coef = '';
        elseif handles.LPApphandle.Orig_Matrix_problem(dim(1), j) == -1
            coef = '-';
        else
            coef = rats(handles.LPApphandle.Orig_Matrix_problem(dim(1), j));
        end
        if handles.LPApphandle.Orig_Matrix_problem(dim(1), j) > 0
            if inicial == 1
                latex = [latex, '+', coef, 'x_{',num2str(j+k),'}'];
            else 
                latex = [latex, coef, 'x_{',num2str(j+k),'}'];
            end
            inicial = 1;
        elseif handles.LPApphandle.Orig_Matrix_problem(dim(1), j) < 0
            latex = [latex, coef, 'x_{',num2str(j+k),'}'];
            inicial = 1;
        end
    else
        if Matrix_problem(dim(1), j) == 1
            coef = '';
        elseif Matrix_problem(dim(1), j) == -1
            coef = '-';
        else
            coef = rats(Matrix_problem(dim(1), j));
        end

        if Matrix_problem(dim(1), j) > 0
            if inicial == 1
                latex = [latex, '+', coef, 'x_{',num2str(j+k),'}'];
            else 
                latex = [latex, coef, 'x_{',num2str(j+k),'}'];
            end
            inicial = 1;
        elseif Matrix_problem(dim(1), j) < 0
            latex = [latex, coef, 'x_{',num2str(j+k),'}'];
            inicial = 1;
        end
    end
    if j == dim(2) - 1-l
        latex = [latex, '} \\'];
    end 
end
latex = sprintf('%s \r', latex);
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{\text{sujeto a}} \\'];
latex = sprintf('%s \r', latex);
inicial = 0;
for i = 1:dim(1) - 1
    for j = 1:dim(2) - 1 - l
        if handles.LPApphandle.Newmethod == 4 && handles.LPApphandle.isonlylecture ~= 1 && ...
            isempty(strfind(title, 'nivel')) || (handles.LPApphandle.istwophases == 1 && ...
            isempty(strfind(title, 'asociado')) || ...
            (handles.LPApphandle.ispenalties == 1 && handles.LPApphandle.Newmethod ~=4))
            if handles.LPApphandle.Orig_Matrix_problem(i, j) == 1
                coef = '';
            elseif handles.LPApphandle.Orig_Matrix_problem(i, j) == -1
                coef = '-';
            else
                coef = rats(handles.LPApphandle.Orig_Matrix_problem(i, j));
            end
            if handles.LPApphandle.Orig_Matrix_problem(i, j) > 0
                if inicial == 1
                    latex = [latex, '+', coef, 'x_{',num2str(j+k), '} & '];
                else
                    latex = [latex, coef, 'x_{',num2str(j+k), '} & '];
                end
                inicial = 1;
            elseif handles.LPApphandle.Orig_Matrix_problem(i, j) < 0
                latex = [latex, coef, 'x_{',num2str(j+k), '} & '];
                inicial = 1;
            else
                latex = [latex, ' & '];
            end
        else 
            if Matrix_problem(i, j) == 1
                coef = '';
            elseif Matrix_problem(i, j) == -1
                coef = '-';
            else
                coef = rats(Matrix_problem(i, j));
            end
            if Matrix_problem(i, j) > 0
                if inicial == 1
                    latex = [latex, '+', coef, 'x_{',num2str(j+k), '} & '];
                else
                    latex = [latex, coef, 'x_{',num2str(j+k), '} & '];
                end
                inicial = 1;
            elseif Matrix_problem(i, j) < 0
                latex = [latex, coef, 'x_{',num2str(j+k), '} & '];
                inicial = 1;
            else
                latex = [latex, ' & '];
            end
        end                         
    end
    inicial = 0;
    if handles.LPApphandle.Newmethod == 4 && handles.LPApphandle.isonlylecture ~= 1 && ...
            isempty(strfind(title, 'nivel'))
        latex = [latex, sep_sistema, rats(handles.LPApphandle.Orig_Matrix_problem(i, end)), '\\'];
    else
        latex = [latex, sep_sistema, rats(Matrix_problem(i, dim(2))), '\\'];
    end
    latex = sprintf('%s \r', latex);
end

latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{x_j \ge 0,' , integral, ' \ (j = ', num2str(1+k), ',\dots,', num2str(dim(2)-1+k-l), ')}'];
latex = sprintf('%s \r', latex);
latex = [latex, '\end{array}']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex); 

    
function generar_latexsimplexbegin(title)
global Matrix_problem Order_initial latex;

Tableau = Matrix_problem;
dim = size(Tableau);

% se convierte la última en términos de Rj =cj - zj
[Y, I] = sort(Order_initial); %#ok<ASGLU>
Tableau_sorted = Matrix_problem(:, I);
c = Tableau_sorted(end,:);
z = Tableau_sorted(end, 1:(dim(1)-1))*Tableau_sorted(1:(dim(1)-1),:);

%R = c - z;
R = eval(sym(str2num(rat(c)))) - eval(sym(str2num(rat(z)))); %#ok<ST2NM>
Tableau(end, I) = R;

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, title];
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
            latex = [latex, rats(Tableau(i, j))]; %#ok<AGROW>
        else
            latex = [latex,'& ', rats(Tableau(i, j))];                 %#ok<AGROW>
        end
    end
    latex = [latex, ' \\']; %#ok<AGROW>
    latex = sprintf('%s \r', latex);
end
latex = [latex, '\end{array} \right]']; 
latex = sprintf('%s \r', latex);
latex = [latex, '\end{equation}'];
latex = sprintf('%s \r', latex);


function generar_latextransportespecification()
global Matrix_problem latex;

dim = size(Matrix_problem);

latex = '';
latex = sprintf('%s \r', latex);
latex = [latex, '\section{Especificación del Problema de Transporte}'];
latex = sprintf('%s \r', latex);
latex = [latex, 'La siguiente tabla representa los costos unitarios asociados a cada variable: \\'];
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
                latex = [latex, '&  \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{$c_{', num2str(i-1),',', num2str(j-1), '}$} ']; %#ok<AGROW>                
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
                latex = [latex, '&  \multicolumn{1}{c}{',rats(Matrix_problem(i-1, j-1)),'} & \multicolumn{1}{c|}{} ']; %#ok<AGROW>                
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
latex = [latex, '\section{Desarrollo del Método de Transporte}'];
latex = sprintf('%s \r', latex);



% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles) %#ok<INUSL>
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global Dimension;

All_display = get(handles.table_problem, 'data');
% se obtienen las dimensiones del problema
Dimension = size(All_display);

if handles.LPApphandle.Newmethod ~= 3
    headers = char('f', 'X', 'Yi0', 'Z', 'lambdai', 'Ui');        
else
    headers = char('O', 'D', 'Oferta', 'Demanda');
end
  
if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
                handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 4
    n = Dimension(2)-1;
elseif handles.LPApphandle.Newmethod == 5
    n = Dimension(2)-2;
end
if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
                handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 5
    m = Dimension(1)-1;
elseif handles.LPApphandle.Newmethod == 4
    m = Dimension(1)-2;
end

if strcmp(get(handles.uitoggletool1, 'State'), 'off');
    if  ~((m+1 > n && (handles.LPApphandle.Newmethod ~= 3 || handles.LPApphandle.Newmethod ~= 4)) || m+1 > 96 || ...
            (m+1 > n + m + 1 && handles.LPApphandle.Newmethod == 4))    
        All_display_aux = All_display;
        All_display(1:Dimension(1)+1,1:Dimension(2)) = zeros(Dimension(1)+1,Dimension(2));
        All_display([1:m, m+2:Dimension(1)+1],1:Dimension(2)) = ...
            All_display_aux(1:Dimension(1),1:Dimension(2));        
    end
else
    % se construye el arreglo de cadenas de las ecuaciones
    if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
        handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 5
        Var = cell(Dimension(1)-1, 1);
        for i = 1:Dimension(1)-1
            Var(i) = cellstr(strcat('f',num2str(i)));
        end
    elseif handles.LPApphandle.Newmethod == 4
        Var = cell(Dimension(1)-2, 1);
        for i = 1:Dimension(1)-2
            Var(i) = cellstr(strcat('f',num2str(i)));
        end
    end
    [s,v] = listdlg('PromptString','Selecciona una ecuación:', 'SelectionMode','single',...
                'ListString', Var, 'InitialValue', m);
     if v ~=0
         if ~((m-1 > n && (handles.LPApphandle.Newmethod ~= 3 || handles.LPApphandle.Newmethod ~= 4)) || m-1 < 1  || ...
                 (m-1 > n + m - 1 && handles.LPApphandle.Newmethod == 4))            
            if s > 1
                All_display = All_display([1:s-1, s+1:Dimension(1)],1:Dimension(2));
            elseif s == 1
                All_display = All_display(2:Dimension(1), 1:Dimension(2));
            end
         end
     end
end

Dimension = size(All_display);
    
rowName=cell(1,Dimension(1));
if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
        handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 5
    for i = 1:Dimension(1)-1
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end

    rowName(Dimension(1)) = cellstr(headers(4,:));
elseif handles.LPApphandle.Newmethod == 4
    for i = 1:Dimension(1)-2
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end

    rowName(Dimension(1)-1) = cellstr(headers(4,:));
    rowName(Dimension(1)) = cellstr(headers(6,:));
    handles.LPApphandle.gui_Matrix_problem(Dimension(1), 1:Dimension(2)-1) = ones(1, Dimension(2)-1);
end

set(handles.table_problem, 'rowname', rowName);

% se despliega la tabla para el ingreso de datos

set(handles.table_problem, 'data', All_display);




% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles) %#ok<INUSL>
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global Dimension;

All_display = get(handles.table_problem, 'data');
% se obtienen las dimensiones del problema
Dimension = size(All_display);

if handles.LPApphandle.Newmethod ~= 3
    headers = char('f', 'X', 'Yi0', 'Z', 'lambdai', 'Ui');        
else
    headers = char('O', 'D', 'Oferta', 'Demanda');        
end

if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
                handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 4
    n = Dimension(2)-1;
elseif handles.LPApphandle.Newmethod == 5
    n = Dimension(2)-2;
end
if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
                handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 5
    m = Dimension(1)-1;
elseif handles.LPApphandle.Newmethod == 4
    m = Dimension(1)-2;
end

if strcmp(get(handles.uitoggletool1, 'State'), 'off');
    if  ~((m > n+1 && (handles.LPApphandle.Newmethod ~= 3 || handles.LPApphandle.Newmethod ~= 4)) || ...
            n+1 > 96 || (m > n+1 + m && handles.LPApphandle.Newmethod == 4))
        All_display_aux = All_display;
        All_display(1:Dimension(1),1:Dimension(2)+1) = zeros(Dimension(1),Dimension(2)+1);
        All_display(1:Dimension(1), [1:n, n+2:Dimension(2)+1]) = ...
            All_display_aux(1:Dimension(1),1:Dimension(2));        
    end
else
    % se construye el arreglo de cadenas de las variables
    if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
        handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 4
        Var = cell(Dimension(2)-1, 1);
        j = 0;
        if handles.LPApphandle.Newmethod == 4 && handles.LPApphandle.isonlylecture == 1
            Dim = size(handles.LPApphandle.Orig_Matrix_problem);
            j = (Dim(2)-1)- (Dimension(2) -1);
        end
        
        for i = 1+j:Dimension(2)-1+j
            Var(i-j) = cellstr(strcat('X',num2str(i)));
        end
    elseif handles.LPApphandle.Newmethod == 5
        Var = cell(Dimension(2)-2, 1);
        for i = 1:Dimension(2)-2
            Var(i) = cellstr(strcat('X',num2str(i)));
        end
    end
    [s,v] = listdlg('PromptString','Selecciona una variable:', 'SelectionMode','single',...
                'ListString', Var, 'InitialValue', n);
     if v ~=0
         if ~((m > n-1 && (handles.LPApphandle.Newmethod ~= 3 || handles.LPApphandle.Newmethod ~= 4)) || n-1 < 2 || ...
                 (m > n-1 + m && handles.LPApphandle.Newmethod == 4))
            if s > 1
                All_display = All_display(1:Dimension(1), [1:s-1, s+1:Dimension(2)]);
            elseif s == 1
                All_display = All_display(1:Dimension(1), 2:Dimension(2));
            end
         end
     end
end

%set(handles.table_problem, 'data', All_display);
Dimension = size(All_display);

colName = cell(Dimension(2), 1);
if handles.LPApphandle.Newmethod == 1 || handles.LPApphandle.Newmethod == 2 || ...
        handles.LPApphandle.Newmethod == 3 || handles.LPApphandle.Newmethod == 4
    j = 0;
    if handles.LPApphandle.Newmethod == 4 && handles.LPApphandle.isonlylecture == 1
        Dim = size(handles.LPApphandle.Orig_Matrix_problem);
        j = (Dim(2)-1)- (Dimension(2) -1);
    end
    for i = (1+j):Dimension(2)-1+j
        colName(i-j) = cellstr(strcat(headers(2,:),num2str(i)));
    end    

    colName(Dimension(2)) = cellstr(headers(3,:));
elseif handles.LPApphandle.Newmethod == 5
    for i = 1:Dimension(2)-2
        colName(i) = cellstr(strcat(headers(2,:),num2str(i)));
    end

    colName(Dimension(2)-1) = cellstr(headers(3,:));
    colName(Dimension(2)) = cellstr(headers(5,:));
end   

set(handles.table_problem, 'columnname', colName);

% se establece el formato de las columnas
colFormat=cell(1,Dimension(2));
for i = 1:Dimension(2)
    colFormat(1,i) = cellstr('rat');
end
set(handles.table_problem, 'columnformat', colFormat);

if handles.LPApphandle.istwophases ~= 1
    % se configuran las columnas como editables
    colEdit = ones(1,Dimension(2));
    set(handles.table_problem, 'columneditable', (colEdit == 1));
end
%se despliega la tabla para el ingreso de datos
set(handles.table_problem, 'data', All_display);



% --- Executes when user attempts to close Problem.
function Problem_CloseRequestFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to Problem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if (handles.LPApphandle.Newmethod == 4 || handles.LPApphandle.istwophases == 1) && ...
        handles.LPApphandle.isonlylecture == 1
    pushbutton_ok_Callback(hObject, eventdata, handles);
else
    delete(hObject);
end


% --------------------------------------------------------------------
function uitoggletool1_OffCallback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipushtool1, 'TooltipString', 'Agregar fila');
set(handles.uipushtool2, 'TooltipString', 'Agregar columna');
set(handles.uitoggletool1, 'TooltipString', 'Elimnar fila o columna');


% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipushtool1, 'TooltipString', 'Eliminar fila');
set(handles.uipushtool2, 'TooltipString', 'Eliminar columna');
set(handles.uitoggletool1, 'TooltipString', 'Agregar fila o columna');
