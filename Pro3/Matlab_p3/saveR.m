function saveR(filename, varargin)

if(nargin < 1), error('Requires at least one input arguments.'); end

if(nargin < 2),
    vars = evalin('caller', 'who');
else
    vars = varargin;
end

fid = fopen(filename,'wt');

for var_index = 1:length(vars),
    var_name = vars{var_index};
    var_namestr = ['"' var_name '" <-'];
    
    var_value = evalin('caller', vars{var_index});
    
    var_size = size(var_value);
    var_sizestr = mat2str(var_size(:));
    var_sizestr = strrep(var_sizestr(2:end-1),';',', ');
    
    if(isstruct(var_value))
        warning('SAVER:structure','R cannot handle structures. File "%s" will be written but will not contain variable "%s".',filename,var_name);
        continue;
    elseif(iscell(var_value)),
        var_valuestr = sprintf('"%s", ',var_value{:});
        var_valuestr = ['structure(c(' var_valuestr(1:end-2) '), .Dim = c(' var_sizestr '))'];
    elseif(isscalar(var_value)),
        var_valuestr = num2str(var_value);
        var_valuestr = strrep(var_valuestr,'NaN','NA');
    else
        var_valuestr = mat2str(var_value(:));
        var_valuestr = strrep(var_valuestr,'NaN','NA');
        var_valuestr = ['structure(c(' strrep(var_valuestr(2:end-1),';',', ') '), .Dim = c(' var_sizestr '))'];         
    end
      
    fprintf(fid, '%s\n%s\n', var_namestr, var_valuestr);
    
end

fclose(fid);

end