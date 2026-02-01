function [patname, chan, clnr] = unpack_Uname(unit_name);

H = strfind(unit_name,'_');;
patname = unit_name(1:H(1)-1);
chan = str2num(unit_name(H(1)+1:H(2)-1));
clnr = str2num(unit_name(H(2)+1:end));

