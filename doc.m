%% regular_3D_polygon_incircle
%
% Function to compute and display the incircle of a regular polygon in 3D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
%% Syntax
%
% regular_3D_polygon_incircle(P);
%
% regular_3D_polygon_incircle(P, nb_samples);
%
% regular_3D_polygon_incircle(P, nb_samples, option_display);
%
% [C, I, r] = regular_3D_polygon_incircle(P, nb_samples, option_display);
%
%
%% Description
%
%
% regular_3D_polygon_incircle(P) computes and displays the incircle of polygon P.
%
% regular_3D_polygon_incircle(P, nb_samples) uses nb_samples to draw the
% circle.
%
% regular_3D_polygon_incircle(P, nb_samples, option_display) displays the
% circle when option_display is set either to logical true or real numeric 1, and
% doesn't when it is set to logical false or real numeric 0.
%
% [C, I, r] = regular_3D_polygon_incircle(P, nb_samples, option_display) stores
% the results in [C, I, r] vector.
%
%
%% See also
%
% | <https://fr.mathworks.com/matlabcentral/fileexchange/119778-triangle-incircle-3d-2d?s_tid=prof_contriblnk triangle incircle> |
%
%
%% Input arguments
%
%        [ |  |  | ]
% - P = [Px Py Pz] : real matrix double. size(P,1) > 2. 2 <= size(P,2) <= 3. P vertex coordinates. 
%        [ |  |  | ]
%
% - nb_samples : integer scalar double. The number of samples to draw the
%                incircle. nb_samples >= 3.
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
%% Output arguments
%
%        [- Cx -]
% - C = [- Cy -]: real matrix double. The incircle coordinates. size(R) = [size(A,1), nb_samples].
%        [- Cz -]
%
%        [Ix]
% - I = [Iy] : real column vector double. 2 <= numel(I) <= 3. The incircle centre.
%        [Iz]
%
% - r : real scalar double. the incircle radius.
%
%
%% Example #1
% Square/rhombus in 3D space
nb_samples = 60;
option_display = true;
V1 = [1 1 0];
V2 = [0 0 sqrt(2)];
V3 = [-1 -1 0];
V4 = [0 0 -sqrt(2)];
P = cat(1,V1,V2,V3,V4);
regular_3D_polygon_incircle(P, nb_samples, option_display);
view(-99, 20);

%% Example #2
% Regular hexagon in 3D space
V1 = [-sqrt(2)/3 sqrt(6)/3 -1/3];
V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3];
V5 = [0 0 1];
V2 = (2/3)*(V1+V3) - (1/3)*V5;
V4 = (2/3)*(V3+V5) - (1/3)*V1;
V6 = (2/3)*(V1+V5) - (1/3)*V3;
P = cat(1,V1,V2,V3,V4,V5,V6);
regular_3D_polygon_incircle(P, nb_samples, option_display);