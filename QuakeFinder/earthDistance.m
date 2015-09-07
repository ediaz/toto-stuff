function dist = earthDistance(loc1, loc2)
% Computes the distance, in km, along the surface of the Earth between the two
% locations, each specified as a [lat lon] vector.
%  Uses a round Earth approximation.
%
% Updated to accept station ID as string for second input argument

earthRadius = 6371.01;

deltaLon = abs( loc1(2) - loc2(2) );
if( deltaLon > 180 ) 
	deltaLon = 360 - deltaLon;
end
%if( loc1(2) - loc2(2) < 0 )
%	deltaLon = deltaLon * -1;
%end

% Vincenty formula
dist = earthRadius * atan2( sqrt( ( cosd(loc1(1)) * sind(deltaLon) )^2 + ( cosd(loc2(1)) * sind(loc1(1)) - sind(loc2(1)) * cosd(loc1(1)) * cosd(deltaLon) )^2 ), sind(loc2(1)) * sind(loc1(1)) + cosd(loc2(1)) * cosd(loc1(1)) * cosd(deltaLon) );

%[dist,az] = distance(loc1(1),loc1(2),loc2(1),loc2(2))
%dist = earthRadius * dist;

