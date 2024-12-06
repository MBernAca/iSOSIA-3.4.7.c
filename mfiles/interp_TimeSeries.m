function [interp_values] = interp_TimeSeries(Xvalue2interp,Yvalue2interp,Xinterp_array)

% First check that times are unique
[Unique_x,ia,ic] = unique(Xvalue2interp);
unique_y = Yvalue2interp(ia);

%interpolate 
interp_values = interp1(Unique_x,unique_y,Xinterp_array);

end