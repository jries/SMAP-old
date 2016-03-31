function out=pluginnames(varargin)
plugs=plugin;
for k=1:length(varargin)
    plugs=plugs.(varargin{k});
end
out=fieldnames(plugs);
end