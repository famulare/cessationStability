function [S,E] = xls2struct(fn,varargin)
% XLS2STRUCT Reads xls-file and converts into a struct array
%
% [S,E] = XLS2STRUCT(FILE) reads the Excel-file FILE and converts the columns
% into separate fields in the structure S. The first line of the file
% must contain unique headers for all of the columns, the rest of the file
% is data. 
%
% Column names can contain spaces, these are converted to underscores in 
% order to obtain proper variable names.
%
% S is the plane representation of the table. E is a second structure with
% an element-wise representation. I.e. notation S.Time(i) versus E(i).Time.
%
% XLS2STRUCT(FILE,SHEET) reads SHEET instead of the default sheet.
%
% Example:
%     City    Time   Temp
%    Dallas   12      98
%    Tulsa    13      99
%    Boise    14      97
%
% converts to
% S.City = {'Dallas','Tulsa','Boise'}';
% S.Time = [12 13 14]';
% S.Temp = [98 99 97]';
%
% E(1)
% ans = 
%    City: 'Dallas'
%    Time: 12
%    Temp: 98
%
% See also: XLSREAD.

% CHANGES
% aha, greatly simplified using Matlab 7 features, 21-jan-05
% aha, minor bug fix and added ' ' to '_' conversion, 25-nov-2004 
% aha, original version, 20-oct-2003


if ~exist(fn,'file')
  error(['Could not open file: ' fn])
end
%if nargin<2, varargin = ''; end
  
[t1,t2,raw] = xlsread(fn,varargin{:});

h = raw(1,:);
h = strrep(h,' ','_'); h = strrep(h,'"','_');
t = raw(2:end,:);

% PLANE storage
S = struct(h{1},[]);
for i = 1:length(h)
  
  c = t(:,i);
  inum = cellfun('isclass',c,'double');
    
  if all(inum)
    % all numeric
    S.(h{i}) = cell2mat(c);
  else
    % all text or mixed
    S.(h{i}) = c;
  end
end

% ELEMENTWISE storage
E = cell2struct(t,h,2);








