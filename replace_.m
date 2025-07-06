function new_id = replace_(id)

% replace '_' to '-'

n = length(id);
for i = 1 : n
   str = id{i}; 
   idx = find(str == '_');
   str(idx) = '-';
   new_id{i} = str;
end
new_id = new_id';