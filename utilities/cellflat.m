%cellflat - flatten nested cells
%
% data = cellflat(data)

function data = cellflat(data)
  try
      data = cellfun(@cellflat,data,'un',0);
      if any(cellfun(@iscell,data))
          data = [data{:}];
      end
  catch
      % a non-cell node, so simply return node data as-is
  end
end
