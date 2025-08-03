function vessels = findFolders(folder, name, except)
if nargin < 3
    except = '';
end

vessels = {};
if (isa(name,'char') || size(name, 2) == 1)
    name = char(name);
    
    dirinfo = dir(folder);
    curlength = size(vessels,2);
    fname = fullfile(folder, name);
    fnamee = fullfile(folder, except);
    if size(dir(fname),1)>0
        dirs = dir(fname);
        dirse = dir(fnamee);
        if(size(dirse,1)==0 || size(except,1) == 0)
            for i = 1: size(dir(fname),1)
                if (exist(fullfile(folder, dirs(i).name), 'file'))>0
                    vessels{curlength+1} = fullfile(folder, dirs(i).name);
                    curlength = curlength + 1;
                end
            end
        end
    end
    for K = 1 : length(dirinfo)
        thisdir = dirinfo(K).name;
        if (~(strcmp(thisdir,'.') || strcmp(thisdir, '..')))
            curfolder = fullfile(folder,thisdir);
            %       if (exist([curfolder, '\vessel.tiff'], 'file'))
            %         vessels{curlength+1} = [curfolder, '\vessel.tiff'];
            %         curlength =curlength + 1;
            %       end
            vessels = [vessels, findFolders(curfolder, name, except)];
        end
    end
else
    for i = 1:size(name,2)
        vessels = [vessels, findFolders(folder, name{i}, except)];
    end
        
    % if (~exist('vessels'))
    %     vessels = 0;
    % end
end
end