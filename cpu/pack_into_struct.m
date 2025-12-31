function S = pack_into_struct(varargin)

S = struct();

n = numel(varargin);
for ii=1:n
    S.(inputname(ii)) = varargin{ii};
end


end