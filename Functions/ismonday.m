function check = ismonday(date)
%
% INPUT
%
% date : datetime
%
% OUTPUT
% 
% check : true if date is monday
    
    check = (weekday(date) == 2);
end