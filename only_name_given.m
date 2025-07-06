function [Xname,idx] = only_name_given(all_name,idx)
% -----------------------------------------
number_of_names = length(all_name);

new_idx = ones(number_of_names,1);
new_idx(idx) = 0;

for i=1:number_of_names
    if new_idx(i) == 0
        name{i} = all_name{i};
    else
        name{i} = [];
    end
end

Xname = name;
for i = number_of_names:-1:2
    if ~isempty(name{i})
        temp = strsplit(name{i},'_');
        if ~strcmp(temp(1),Xname)
            Xname(i) = temp(1);
        else
            Xname{i} = [];
        end
    else
        Xname{i} = [];
    end
end
Xname = Xname';