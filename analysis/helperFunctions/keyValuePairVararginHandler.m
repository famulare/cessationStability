function options = keyValuePairVararginHandler(options,inputs)
% inputs
%   options stuct: key -value pairs
%   inputs = varargin of calling function
%

if numel(inputs)==1 && isstruct(inputs{1})
    f=fieldnames(inputs{1});
    d=struct2cell(inputs{1});
    clear inputs
    inputs={};
    for k=1:length(f)
        inputs{end+1}=f{k};
        inputs{end+1}=d(k);
    end
end
    

for pair = reshape(inputs,2,[])
    inpName = pair{1}; 
    if any(strcmp(inpName, fieldnames(options)))
        options.(inpName) = pair{2};
    else  error('%s is not a recognized parameter name.\n', inpName);
    end;
end;

end