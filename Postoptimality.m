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

function varargout = Postoptimality(varargin)
% POSTOPTIMO M-file for Postoptimo.fig
%      POSTOPTIMO, by itself, creates a new POSTOPTIMO or raises the existing
%      singleton*.
%
%      H = POSTOPTIMO returns the handle to a new POSTOPTIMO or the handle to
%      the existing singleton*.
%
%      POSTOPTIMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTOPTIMO.M with the given input arguments.
%
%      POSTOPTIMO('Property','Value',...) creates a new POSTOPTIMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Postoptimo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Postoptimo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Postoptimo

% Last Modified by GUIDE v2.5 30-Apr-2012 22:15:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Postoptimality_OpeningFcn, ...
                   'gui_OutputFcn',  @Postoptimality_OutputFcn, ...
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


% --- Executes just before Postoptimo is made visible.
function Postoptimality_OpeningFcn(hObject, eventdata, handles, varargin)  %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Postoptimo (see VARARGIN)

% Choose default command line output for Postoptimality
handles.output = hObject;

movegui(hObject, 'center');

LPApphandle = find(strcmp(varargin, 'LPApp'));
if ~isempty(LPApphandle)
   handles.LPApphandle = varargin{LPApphandle+1};
end
% Update handles structure
guidata(hObject, handles);

dim = size(handles.LPApphandle.gui_Matrix_problem);
colName = cell(dim(2), 1);
for i = 1:(dim(2)-1)
    colName(i) = cellstr(strcat('X',num2str(i)));
end
colName(dim(2)) = cellstr('Yi0');
set(handles.table_tableau, 'columnname', colName);
set(handles.table_problem, 'columnname', colName);

colFormat=cell(1,dim(2));
for i = 1:dim(2)
    colFormat(1,i) = cellstr('rat');
end
set(handles.table_tableau, 'columnformat', colFormat);
set(handles.table_problem, 'columnformat', colFormat);

rowName=cell(1,dim(1));
for i = 1:(dim(1)-1)
    rowName(1,i) = cellstr(strcat('f',num2str(i)));
end
[Z, I] = sort(handles.LPApphandle.Order_current); %#ok<ASGLU>
rowName=cell(1,dim(1));
for i = 1:(dim(1)-1)
    rowName(1,i) = cellstr(strcat('X',num2str(I(i))));
end
rowName(dim(1)) = cellstr('Rj');

set(handles.table_tableau, 'rowname', rowName);
rowName(dim(1)) = cellstr('Z');
set(handles.table_problem, 'rowname', rowName);

colEdit = ones(1,dim(2));
set(handles.table_problem, 'columneditable', (colEdit == 1));

set(handles.table_problem, 'data', handles.LPApphandle.gui_Matrix_problem);
set(handles.table_tableau, 'data', handles.LPApphandle.gui_tableau);


% UIWAIT makes Postoptimo wait for user response (see UIRESUME)
% uiwait(handles.postoptimality);


% --- Outputs from this function are returned to the command line.
function varargout = Postoptimality_OutputFcn(hObject, eventdata, handles)   %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_loadproblem.
function pushbutton_loadproblem_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_loadproblem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dim = size(handles.LPApphandle.gui_Matrix_problem);
All_display_tableau = get(handles.table_tableau, 'data');
b = All_display_tableau(1:(dim(1)-1), end);
Rj = All_display_tableau(end, 1:(dim(2)-1));
if any(b < 0) && all(Rj >= 0)
    if strcmp(get(handles.LPApphandle.Simplex, 'Checked'), 'on')
        msgbox('Seleccione el método Simplex dual', 'Cambie de método', 'help');
        return;
    end
elseif all(b > 0) && any(Rj < 0)
    if strcmp(get(handles.LPApphandle.Simplex_dual, 'Cheched'), 'on')
        msgbox('Seleccione el método Simplex primal', 'Cambie de método', 'help');
        return;
    end
end

handles.LPApphandle.setProblemAndTableau(get(handles.table_problem, 'data'), get(handles.table_tableau, 'data'), handles.LPApphandle);
delete(handles.postoptimality);


% --- Executes on button press in pushbutton_restart.
function pushbutton_restart_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to pushbutton_restart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA

set(handles.table_problem, 'data', handles.LPApphandle.gui_Matrix_problem)
set(handles.table_tableau, 'data', handles.LPApphandle.gui_tableau)

% --- Executes when entered data in editable cell(s) in table_problem.
function table_problem_CellEditCallback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to table_problem (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
dim = size(handles.LPApphandle.gui_Matrix_problem);
Order_current = handles.LPApphandle.Order_current;
Basic_vector = Order_current < dim(1);
Order_initial = handles.LPApphandle.Order_initial;
Basic_vector_initial = Order_initial < dim(1);
All_display_problem = get(handles.table_problem, 'data');
All_display_tableau = get(handles.table_tableau, 'data');
cellindices = eventdata.Indices;
row = cellindices(1);
column = cellindices(2);
if row < dim(1) && column < dim(2)
    errordlg('Característica aún no implementada', 'No implementado');
    All_display_problem(row, column) = eventdata.PreviousData;
    set(handles.table_problem, 'data', All_display_problem);
    return;
end

if row == dim(1) && column ~= dim(2)
    if Basic_vector(column)
        NewRj = zeros(1,dim(2)-1);
        for i = 1:(dim(2)-1)
            if ~Basic_vector(i)
                NewRj(i) = All_display_tableau(row, i) - (eventdata.PreviousData - eventdata.NewData)*All_display_tableau(column, i);
            end
        end
        All_display_tableau(end, 1:(dim(2)-1)) = NewRj;
        All_display_tableau(end, end) = All_display_problem(end, Basic_vector==1)*All_display_tableau(1:(dim(1)-1), end);
    else
        All_display_tableau(end,column) = All_display_tableau(end,column)-(eventdata.NewData - eventdata.PreviousData);        
    end
elseif column == dim(2) && row ~= dim(1)
    bprime = All_display_problem(1:(dim(1)-1), end);
    b = All_display_problem(1:(dim(1)-1), end);
    b(row) = eventdata.PreviousData;
    All_display_tableau(1:(dim(1)-1), end) = All_display_tableau(1:(dim(1)-1), end)+All_display_tableau(1:(dim(1)-1), Basic_vector_initial==1)*(bprime-b);
    All_display_tableau(end, end) = All_display_problem(end, Basic_vector==1)*All_display_tableau(1:(dim(1)-1), end);
end
set(handles.table_tableau, 'data', All_display_tableau);
