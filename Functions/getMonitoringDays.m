function dd = getMonitoringDays(start, finish)
%
% INPUT
%
% start : first monitoring date
% finish : last monitoring date
%
% OUTPUT
%
% dd : vector of monitoring dates between start and finish

    dd(1) = start;
    add = start;

    while ~ismonday(add)  % add days until i get to monday as requested
        add = add + caldays(1);
    end

    idx = 2;
    while add<finish
        dd(idx) = add;
        idx = idx + 1;
        add = add + caldays(7);
    end

    dd = [dd finish];

end