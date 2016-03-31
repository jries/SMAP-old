function destination=copystructReduce(source,ind)

    if ~isempty(source)
    fn=fieldnames(source);
    for k=1:length(fn)
        if length(ind)==length(source.(fn{k}))
        destination.(fn{k})=source.(fn{k})(ind);
        end
    end
    end
