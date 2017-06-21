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
    setheading(handles, char('f', 'X', 'Yi0', 'Z', 'lambdai'));
else
    setheading(handles, char('Origen', 'Destino', 'Oferta', 'Demanda'));
end
% UIWAIT makes problem wait for user response (see UIRESUME)
% uiwait(handles.Problem);

function setheading(handles, headers)
global Dimension;

% se obtienen las dimensiones del problema
Dimension = size(handles.LPApphandle.gui_Matrix_problem);
Dim = size (headers);
if Dim(1) == 5
    col_plus = 1;
    All_display = zeros(Dimension(1), Dimension(2));    
else
    col_plus = 0;
end

% se rotulan las filas y columnas de las tablas
colName = cell(Dimension(2) + col_plus, 1);
for i = 1:Dimension(2)-1
    colName(i) = cellstr(strcat(headers(2,:),num2str(i)));
end
colName(Dimension(2)) = cellstr(headers(3,:));
if Dim(1) == 5
    colName(Dimension(2)+1) = cellstr(headers(5,:));
end

set(handles.table_problem, 'columnname', colName);

rowName=cell(1,Dimension(1));
for i = 1:Dimension(1)-1
    rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
end
rowName(Dimension(1)) = cellstr(headers(4,:));
set(handles.table_problem, 'rowname', rowName);
% se establece el formato de las columnas
colFormat=cell(1,Dimension(2)+col_plus);
for i = 1:(Dimension(2)+col_plus)
    colFormat(1,i) = cellstr('rat');
end
set(handles.table_problem, 'columnformat', colFormat);
% se configuran las columnas como editables
colEdit = ones(1,Dimension(2));
set(handles.table_problem, 'columneditable', (colEdit == 1));
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
function pushbutton_ok_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Matrix_problem Order_initial latex;
% Matrix_problem    guarda la especificación del problema
% Orden_initial     registra el orden de las columnas según la matriz
% identidad

% recupera la especificación dada por el usuario
All_display = get(handles.table_problem, 'data');
Matrix_problem = All_display;
Dim = size(Matrix_problem);
ident = eye(Dim(1)-1);

% se verifica que la especificación este en forma correcta para el método
% correspondiente
if ~iscorrectform(handles)
    return;
end

handles.LPApphandle.Method = handles.LPApphandle.Newmethod;

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

set(handles.LPApphandle.Restriction_nonnegativity, 'checked', 'off');
set(handles.LPApphandle.Restriction_nonnegativity, 'Enable', 'off');
set(handles.LPApphandle.Uncut_planes, 'Enable', 'off');

if handles.LPApphandle.Method ~= 3
    latex = '';
    generar_latexspecification('\section{Especificación del problema}');
    handles.LPApphandle.latex = latex;    
    guidata(handles.LPApphandle.output, handles.LPApphandle);
else
    generar_latextransportespecification()
    handles.LPApphandle.latex = latex;
    handles.LPApphandle.latexproblem = latex;
    guidata(handles.LPApphandle.output, handles.LPApphandle);
end

%AGREGADO(27/12/2016)
if handles.LPApphandle.Method ~= 3 && handles.LPApphandle.Method ~= 2
    resp = iscanonicform();
    if resp == 0
        handles.LPApphandle.istwophases = 1;
        handles.LPApphandle.whatphase = 1;
        handles.LPApphandle.isinit_secondphase = 0;
        handles.LPApphandle.Orig_Matrix_problem = Matrix_problem;
        daux = size(Order_initial(Order_initial>0));
        handles.LPApphandle.maxcanon_vector = daux(2);
        aux = Matrix_problem(1:Dim(1)-1, Dim(2));
        Matrix_problem(Dim(1), :) = zeros(1, Dim(2));
        j=0;       
        %for i = (Order+1:Dim(1)-1)
        for i = (1:Dim(1)-1)
            ind = find(Order_initial==i, 1);
            if isempty(ind)
                Matrix_problem(1:Dim(1)-1, Dim(2) + j) = ident(:, i);
                Matrix_problem(Dim(1), Dim(2) + j) = 1;
                Order_initial(Dim(2) + j) = i;                
                j=j+1;
            end
        end
        Order_initial(Order_initial == 0) = Dim(1):Dim(2) + j - 1;
        Order_initial(Dim(2) + j) = Dim(2) + j;
        Matrix_problem(1:Dim(1)-1, Dim(2) + j) = aux;
        %Dimension = size(Matrix_problem);
        handles.LPApphandle.First_Matrix_problem = Matrix_problem;
        handles.LPApphandle.First_Order_initial = Order_initial;
        
        latex = '\section{Primera fase del Método Simplex}';
        generar_latexspecification('\subsection{Especificación del problema asociado}');
        handles.LPApphandle.latex = [handles.LPApphandle.latex, latex];
        generar_latexsimplexbegin('\subsection{Desarrollo del método simplex}');
        handles.LPApphandle.latex = [handles.LPApphandle.latex, latex];
        handles.LPApphandle.latexIIphasesproblem = handles.LPApphandle.latex;
        guidata(handles.LPApphandle.output, handles.LPApphandle);        
    else        
        if resp == 2
            return;
        end
        generar_latexsimplexbegin('\section{Desarrollo del Método Simplex}');
        handles.LPApphandle.latex = [handles.LPApphandle.latex, latex];    
        handles.LPApphandle.latexproblem = latex;
        guidata(handles.LPApphandle.output, handles.LPApphandle); 
    end
end
if handles.LPApphandle.Method ~= 3
    if isdegeneratedsolution(handles)        
        errordlg(char('La solución inicial es degenerada.'),char('Método simplex'),'modal');
    end
end
%AGREGADO(27/12/2016)

handles.LPApphandle.Order_initial = Order_initial; % se comparte el orden inicial
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
        
elseif rank(full(Matrix_problem(1:Dimension(1)-1, 1:Dimension(2)-1))) ~= Dimension(1)-1
   errordlg('El problema tiene restricciones redundantes.',char('Método Simplex'),'modal');
   iscorrect = 0;
   return;
end

if handles.LPApphandle.Newmethod == 2 
    resp = iscanonicform();
    if resp == 0 || resp == 2
        iscorrect = 0;
        return;
    end
end

if ~isfactiblesolution(handles)
    iscorrect = 0;
    return;
end

if isdegeneratedsolution(handles)
    if handles.LPApphandle.Newmethod == 3
        errordlg(char('La oferta (demanda) de algún origen (destino) es cero.'),char('Origen o destino inútil'),'modal');
        iscorrect = 0;
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

    Rj = c - z;    
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
end
response = 1;

% ---------

function response = isdegeneratedsolution(handles)
global Matrix_problem Dimension;

% se verifica que la solución no sea degenerada
if handles.LPApphandle.Newmethod ~= 3
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
response = zeros(1,Dimension(1)-1);
for i =1:(Dimension(1)-1)
    for j = 1:(Dimension(2)-1)            
        if identidad(:, i) == Matrix_problem(1:Dimension(1)-1, j)
            Order_initial(j) = i;
            response(i) = 1;
            break;        
        end
    end

end
if ~all(response(:)==1)
    %MODIFICADO(27/12/2016)
    resp = questdlg(char('La matriz no está en forma canónica. ¿Desea aplicar el método de dos fases?'),char('Método simplex de dos fases'), ...
        'De acuerdo','No','Cancelar','Cancelar');
    if ~isempty(resp)
        if strcmp(resp, 'De acuerdo')
            response = 0;                
        elseif strcmp(resp, 'No') || strcmp(resp, 'Cancelar')
            response = 2;
        end
    else
        response = 2;
    end
    %MODIFICADO(27/12/2016)
    return;
end
Order_initial(Order_initial == 0) = Dimension(1):Dimension(2); % se etiquetan las demás columnas
response = 1;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.Problem);

function generar_latexspecification(title)
global Matrix_problem latex;

dim = size(Matrix_problem);

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
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{minimizar \ '];
inicial = 0;
for j = 1:dim(2) - 1
    if Matrix_problem(dim(1), j) > 0
        if inicial == 1
            latex = [latex, ...
            '+', rats(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
        else 
            latex = [latex, ...
            rats(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
        end
        inicial = 1;
    elseif Matrix_problem(dim(1), j) < 0
        latex = [latex, ...
        rats(Matrix_problem(dim(1), j)), 'x_',num2str(j)];
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
                '+', rats(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
            else
                latex = [latex, ...
                rats(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
            end
            inicial = 1;
        elseif Matrix_problem(i, j) < 0
            latex = [latex, ...
            rats(Matrix_problem(i, j)), 'x_',num2str(j), ' & '];
            inicial = 1;
        else
            latex = [latex, ' & '];
        end           
    end
    inicial = 0;
    latex = [latex, ' = & ', rats(Matrix_problem(i, dim(2))), '\\'];
    latex = sprintf('%s \r', latex);
end
latex = [latex, '\multicolumn{',num2str(dim(2)+1),'}{l}{x_j \ge 0 \ (j = 1, \dots,', num2str(dim(2)-1), ')}'];
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

R = c - z;
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

if  ~((Dimension(1) > Dimension(2) - 1 && handles.LPApphandle.Newmethod ~= 3) || (Dimension(1) > 96 || Dimension(2) - 1 > 96))
    if handles.LPApphandle.Newmethod ~= 3
        headers = char('f', 'X', 'Yi0', 'Z', 'lambdai');
        %col_plus = 1;
    else
        headers = char('Origen', 'Destino', 'Oferta', 'Demanda');
        %col_plus = 0;
    end
    
    All_display_aux = All_display;
    All_display(1:Dimension(1)+1,1:Dimension(2)) = zeros(Dimension(1)+1,Dimension(2));
    All_display([1:(Dimension(1)-1), (Dimension(1)+1)],1:Dimension(2)) = ...
        All_display_aux(1:Dimension(1),1:Dimension(2));
    
    set(handles.table_problem, 'data', All_display);
    

    rowName=cell(1,Dimension(1)+1);
    for i = 1:Dimension(1)
        rowName(1,i) = cellstr(strcat(headers(1,:),num2str(i)));
    end
    rowName(Dimension(1)+1) = cellstr(headers(4,:));
    set(handles.table_problem, 'rowname', rowName);
    
    % se despliega la tabla para el ingreso de datos

    set(handles.table_problem, 'data', All_display);
end




% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles) %#ok<INUSL>
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global Dimension;

All_display = get(handles.table_problem, 'data');
% se obtienen las dimensiones del problema
Dimension = size(All_display);

if  ~((Dimension(1) - 1 > Dimension(2) && handles.LPApphandle.Newmethod ~= 3) || (Dimension(1) > 96 || Dimension(2) - 1 > 96))
    if handles.LPApphandle.Newmethod ~= 3
        headers = char('f', 'X', 'Yi0', 'Z', 'lambdai');
        col_plus = 1;
    else
        headers = char('Origen', 'Destino', 'Oferta', 'Demanda');
        col_plus = 0;
    end
    
    All_display_aux = All_display;
    All_display(1:Dimension(1),1:Dimension(2)+1) = zeros(Dimension(1),Dimension(2)+1);
    All_display(1:Dimension(1), [1:(Dimension(2)-1), (Dimension(2)+1)]) = ...
        All_display_aux(1:Dimension(1),1:Dimension(2));
    
    set(handles.table_problem, 'data', All_display);
    % se rotulan las filas y columnas de las tablas
    colName = cell(Dimension(2)+1, 1);
    for i = 1:Dimension(2)
        colName(i) = cellstr(strcat(headers(2,:),num2str(i)));
    end
    colName(Dimension(2)+1) = cellstr(headers(3,:));
    if col_plus == 1
        colName(Dimension(2)+2) = cellstr(headers(5,:));
    end
    set(handles.table_problem, 'columnname', colName);

    
    % se establece el formato de las columnas
    colFormat=cell(1,Dimension(2)+1);
    for i = 1:(Dimension(2)+1)
        colFormat(1,i) = cellstr('rat');
    end
    set(handles.table_problem, 'columnformat', colFormat);
    % se configuran las columnas como editables
    if col_plus == 1
        colEdit = ones(1,Dimension(2)+1); %OJO
    else
        colEdit = ones(1,Dimension(2)+1);
    end
    set(handles.table_problem, 'columneditable', (colEdit == 1));
    % se despliega la tabla para el ingreso de datos

    set(handles.table_problem, 'data', All_display);
end
