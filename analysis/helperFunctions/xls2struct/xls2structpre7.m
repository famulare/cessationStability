function [S,E] = xls2struct(fn,sn)
% XLS2STRUCT Reads xls-file and converts into a struct array
%
% [S,E] = XLS2STRUCT(FILE) reads the Excel-file FILE and converts the columns
% into separate fields in the structure S. The first line of the file
% must contain unique headers for all of the columns. Names may contain
% spaces but no special characters that might be interpreted by Matlab. The
% rest of the file is data. 
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

% aha, ver 1.0, 20-oct-2003
% aha, 5-jan-04, blanks in header names allowed

if ~exist(fn,'file')
  error(['Could not open file: ' fn])
end
if exist('sn','var') & ~isempty(sn)
  [t1,t2] = xlsread(fn,sn);
else
  [t1,t2] = xlsread(fn);
end

if (size(t1,1)==size(t2,1)) 
  t1 = t1(2:end,:);
end
h = t2(1,:);
t2 = t2(2:end,:);

% ugly workaround for leading and trailing text columns
if size(t1,2)<size(t2,2)
  i = 1;
  % this loop does the leading blanks: for each column that contains only
  % text add one column of NaN's 
  while 1
    if  any(cellfun('isempty',t2(:,i)))
      break
    end
    t1 = [NaN*ones(size(t1,1),1) t1];
    i = i+1;
  end
  % the rest must be trailing columns
  t1 = [t1, NaN*ones(size(t1,1),size(t2,2)-size(t1,2))];
end

nr = size(t1,1);


% PLANE storage
S = struct(h{1},[]);
for i = 1:length(h)
  
  % replace blanks in header name with underscore
  h{i} = strrep(h{i},' ','_');
  
  d1 = t1(:,i);
  d2 = t2(:,i);
  
  i1 = isnan(d1);
  i2 = cellfun('isempty',d2);
  
  if all(i2)
    % all numeric
    S.(h{i}) = d1;
  elseif all(i1)
    % all text
    S.(h{i}) = d2;
  else
    % mixed, conv numbers to cells
    d2(~i1) = mat2cell(d1(~i1),ones(sum(~i1),1),1);
    S.(h{i}) = d2;
  end
end

% ELEMENTWISE storage
i1 = isnan(t1); 
[mm,nn] = size(t1(~i1));
if isempty(t2)
  t2 = cell(size(t1));
end
t2(~i1) = mat2cell(t1(~i1),ones(mm,1),ones(nn,1));
E = cell2struct(t2,h,2);








