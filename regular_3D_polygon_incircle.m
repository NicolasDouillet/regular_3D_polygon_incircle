function [C, I, r] = regular_3D_polygon_incircle(P, nb_samples, option_display)
%% Function to compute and display the incircle of a regular polygon in 3D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
% Syntax
%
% regular_3D_polygon_incircle(P);
% regular_3D_polygon_incircle(P, nb_samples);
% regular_3D_polygon_incircle(P, nb_samples, option_display);
% [C, I, r] = regular_3D_polygon_incircle(P, nb_samples, option_display);
%
%
% Description
%
%
% regular_3D_polygon_incircle(P) computes and displays the incircle of polygon P.
% regular_3D_polygon_incircle(P, nb_samples) uses nb_samples to draw the
% circle.
% regular_3D_polygon_incircle(P, nb_samples, option_display) displays the
% circle when option_display is set either to logical true or real numeric 1, and
% doesn't when it is set to logical false or real numeric 0.
% [C, I, r] = regular_3D_polygon_incircle(P, nb_samples, option_display) stores
% the results in [C, I, r] vector.
%
% See also INCENTER CIRCUMCENTER
%
%
% Input arguments
%
%       [ | |  |]
% - P = [Py Py Pz] : real matrix double. 2 <= size(A,2) <= 3. P vertex coordinates. 
%       [ | |  |]
%
% - nb_samples : integer scalar double. The number of samples to draw the
%                incircle. nb_samples >= 3.
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
% Output arguments
%
%       [- Cx -]
% - C = [- Cy -]: real matrix double. The incircle coordinates. size(R) = [size(A,1), nb_samples].
%       [- Cz -]
%
%       [Ix]
% - I = [Iy] : real column vector double. 2 <= numel(I) <= 3. The incircle centre.
%       [Iz]
%
% - r : real scalar double. the incircle radius.
%
%
% Example #1
% Square/rhombus in 3D space
% nb_samples = 60;
% option_display = true;
% V1 = [1 1 0];
% V2 = [0 0 sqrt(2)];
% V3 = [-1 -1 0];
% V4 = [0 0 -sqrt(2)];
% P = cat(1,V1,V2,V3,V4);
% regular_3D_polygon_incircle(P, nb_samples, option_display);
% view(-99, 20);
%
% Example #2
% Regular hexagon in 3D space
% V1 = [-sqrt(2)/3 sqrt(6)/3 -1/3];
% V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3];
% V5 = [0 0 1];
% V2 = (2/3)*(V1+V3) - (1/3)*V5;
% V4 = (2/3)*(V3+V5) - (1/3)*V1;
% V6 = (2/3)*(V1+V5) - (1/3)*V3;
% P = cat(1,V1,V2,V3,V4,V5,V6);
% regular_3D_polygon_incircle(P, nb_samples, option_display);


%% Input parsing
assert(nargin > 0, 'Not enought input arguments.');
assert(nargin < 6, 'Too many input arguments.');

if nargin < 3
    
    option_display = true;
    
    if nargin < 2
       
        nb_samples = 60;
        
    end
    
end


dimension = size(P,2);
assert(dimension > 1 && dimension < 4,'Input point set must be of dimension 2 or 3.');


%% Body
N = size(P,1); % number of vertices / edges of the polygon

if dimension < 3 % one padding in 2D case    
    
    P = cat(2,P,ones(size(P,1),1));    
    
end

% Polygon centre
I = mean(P,1);

% Polygon circumcircle radius
R = sqrt(sum((P(1,:)-I).^2,2));

% Polygon incircle radius
r = R*cos(pi/N);

% Compute 2D circle
theta = linspace(0,2*pi,nb_samples);

Cx = r*cos(theta);
Cy = r*sin(theta);
Cz = zeros(1,numel(theta));


% Rotate circle
if dimension > 2
    
    
    % One normal vector to the polygon plane
    % ABC triangle director and normal vectors computation     
    n = cross(P(2,:)-P(1,:),P(3,:)-P(1,:));
    
    % Vector u to rotate around
    k = [0 0 1]';
    u = cross(k,n)/norm(cross(k,n));
    
    % Angle between k and u
    alpha = atan2(norm(cross(k,n)),dot(k,n));    
    
    % 3D rotation matrix around u vector
    Rm = @(delta)[u(1)^2+cos(delta)*(1-u(1)^2) (1-cos(delta))*u(1)*u(2)-u(3)*sin(delta) (1-cos(delta))*u(1)*u(3)+u(2)*sin(delta);
                  (1-cos(delta))*u(1)*u(2)+u(3)*sin(delta) u(2)^2+cos(delta)*(1-u(2)^2) (1-cos(delta))*u(2)*u(3)-u(1)*sin(delta);
                  (1-cos(delta))*u(1)*u(3)-u(2)*sin(delta) (1-cos(delta))*u(2)*u(3)+u(1)*sin(delta) u(3)^2+cos(delta)*(1-u(3)^2)];
                  
    C = (Rm(alpha) * cat(1,Cx,Cy,Cz))' + I;
    
else % if dimension == 2
    
    size(cat(1,Cx,Cy,Cz)')
    size(I)
    
    
    C = cat(1,Cx,Cy,Cz)' + I;
    
    % one simplifications in 2D case
    C = C(:,1:2); 
    I = I(1:2);
    
end


%% Display
if option_display
    
    figure    
        
    if dimension > 2
        
        line([P(:,1); P(1,1)],[P(:,2); P(1,2)],[P(:,3); P(1,3)],'Color',[1 0 0],'Linewidth',2), hold on;
        line(C(:,1),C(:,2),C(:,3),'Color',[0 0 1],'Linewidth',2), hold on;    
        view(3);
    
    else % if dimension == 2
        
        line([P(:,1); P(1,1)],[P(:,2); P(1,2)],'Color',[1 0 0],'Linewidth',2), hold on;
        line(C(:,1),C(:,2),'Color',[0 0 1],'Linewidth',2), hold on; 
        view(2);
        
    end
    
    axis equal, axis tight;    
    
end


end % regular_3D_polygon_incircle