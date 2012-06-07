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

function varargout = Sensibility(varargin)
% SENSIBILITY M-file for Sensibility.fig
%      SENSIBILITY, by itself, creates a new SENSIBILITY or raises the existing
%      singleton*.
%
%      H = SENSIBILITY returns the handle to a new SENSIBILITY or the handle to
%      the existing singleton*.
%
%      SENSIBILITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SENSIBILITY.M with the given input arguments.
%
%      SENSIBILITY('Property','Value',...) creates a new SENSIBILITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sensibility_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sensibility_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sensibility

% Last Modified by GUIDE v2.5 30-Apr-2012 10:42:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Sensibility_OpeningFcn, ...
                   'gui_OutputFcn',  @Sensibility_OutputFcn, ...
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


% --- Executes just before Sensibility is made visible.
function Sensibility_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sensibility (see VARARGIN)

% Choose default command line output for Sensibility
handles.output = hObject;

movegui(hObject, 'center');

LPApphandle = find(strcmp(varargin, 'LPApp'));
if ~isempty(LPApphandle)
   handles.LPApphandle = varargin{LPApphandle+1};
end
% Update handles structure
guidata(hObject, handles);

dim = size(handles.LPApphandle.gui_Matrix_problem);
colName = cellstr(char('Valor actual', 'Cota inferior', 'Cota superior'));
set(handles.table_invariableInterval, 'columnname', colName);

colFormat=cellstr(char('rat', 'rat', 'rat'));
set(handles.table_invariableInterval, 'columnformat', colFormat');

rowName=cell(1,dim(2)+dim(1)-2);
for i = 1:(dim(2)-1)
    rowName(1,i) = cellstr(strcat('c',num2str(i)));
end

for i = dim(2):(dim(2)+dim(1)-2)
    rowName(1,i) = cellstr(strcat('b',num2str(i-(dim(2)-1))));
end

set(handles.table_invariableInterval, 'rowname', rowName);

colEdit = zeros(1,3);
set(handles.table_invariableInterval, 'columneditable', (colEdit == 1));

setInvaribleInterval(handles.LPApphandle, handles)

% UIWAIT makes Sensibility wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function setInvaribleInterval(LPApphandle, handles)
dim = size(LPApphandle.gui_Matrix_problem);

Order_current = LPApphandle.Order_current(1, 1:(dim(2)-1));
Basic_vector = Order_current < dim(1);
ck = LPApphandle.gui_Matrix_problem(end, 1:(dim(2)-1));
rj = LPApphandle.gui_tableau(end, Basic_vector==0);
alphakj = LPApphandle.gui_tableau(1:(dim(1)-1), Basic_vector==0);
minck = zeros(dim(2)-1, 1); 
maxck = zeros(dim(2)-1, 1); 
for i = 1:(dim(2)-1)
    if i < dim(1)
        extreme = ck(i) - rj./alphakj(i,:);
        if ~isempty(extreme(alphakj(i,:) > 0))
            minck(i) = max(extreme(alphakj(i,:) > 0));
        else
            minck(i) = Inf;
        end
        if ~isempty(extreme(alphakj(i,:) < 0))
            maxck(i) = min(extreme(alphakj(i,:) < 0));
        else
            maxck(i) = Inf;
        end
    else
        minck(i) = ck(i)-rj(i-(dim(1)-1));
        maxck(i) = Inf;
    end
end
Order_initial = LPApphandle.Order_initial(1, 1:(dim(2)-1));
Basic_vector_initial = Order_initial < dim(1);
bk = LPApphandle.gui_Matrix_problem(1:(dim(1)-1), end);
yk = LPApphandle.gui_tableau(1:(dim(1)-1), end);
BInv = LPApphandle.gui_tableau(1:(dim(1)-1), Basic_vector_initial==1);

minbk = zeros(dim(1)-1, 1); 
maxbk = zeros(dim(1)-1, 1); 
for i = 1:(dim(1)-1)
    extreme = bk(i) - yk./BInv(:,i);
    if ~isempty(extreme(BInv(:,i) > 0))
        minbk(i) = max(extreme(BInv(:,i) > 0));
    else
        minbk(i) = Inf;
    end
    if ~isempty(extreme(BInv(:,i) < 0))
        maxbk(i) = min(extreme(BInv(:,i) < 0));
    else
        maxbk(i) = Inf;
    end
end

All_display = zeros(dim(2)+dim(1)-2, 3);
All_display(1:(dim(2)-1), 1) = ck;
All_display(1:(dim(2)-1), 2) = minck;
All_display(1:(dim(2)-1), 3) = maxck;
All_display(dim(2):(dim(2)+dim(1)-2), 1) = bk;
All_display(dim(2):(dim(2)+dim(1)-2), 2) = minbk;
All_display(dim(2):(dim(2)+dim(1)-2), 3) = maxbk;
set(handles.table_invariableInterval, 'data', All_display);

% --- Outputs from this function are returned to the command line.
function varargout = Sensibility_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
