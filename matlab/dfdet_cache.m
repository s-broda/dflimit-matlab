function v = dfdet_cache(op,name,key,val)
%DFDET_CACHE  The caches of the general engine: prepared designs, outer paths, refined lines.
%   v = DFDET_CACHE('get',name,key) is the value stored under key in the cache name ([] if none);
%   DFDET_CACHE('put',name,key,val) stores it (each cache keeps at most 64 entries, the oldest dropped beyond);
%   DFDET_CACHE('clear') drops every entry (the values do not change; the next call builds them again).
persistent C O
if isempty(C), C = struct(); O = struct(); end
v = [];
switch op
    case 'get'
        if isfield(C,name) && isKey(C.(name),key), v = C.(name)(key); end
    case 'put'
        if ~isfield(C,name), C.(name) = containers.Map('KeyType','char','ValueType','any'); O.(name) = {}; end
        m = C.(name);
        if ~isKey(m,key), O.(name){end+1} = key; end
        m(key) = val; %#ok<NASGU>
        while numel(O.(name)) > 64
            remove(C.(name),O.(name){1}); O.(name)(1) = [];
        end
        v = val;
    case 'clear'
        C = struct(); O = struct();
end
end
