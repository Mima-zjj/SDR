function val = get_opt(opts, name, defaultVal)

if isfield(opts, name)
    val = opts.(name);
else
    val = defaultVal;
end