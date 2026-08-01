% Visualization code for Horn designs. Includes sliders for realtime 
% geometry manipulation. 
% Code by Mostafa Sedky
% August 1, 2026

function pMUT_Horn_GUI
clc; close all;

%% ===================== PARAMETERS =====================
p.pmut_diameter = 0.7e-3;
p.pmut_radius   = p.pmut_diameter/2;

p.r_base = p.pmut_radius * 2.4;
p.r_top  = 0.0018;
p.height = 0.0031;

p.p = 4.93;
p.q = 6.52;

p.A_theta = 0;
p.m_theta = 0;

p.A_x = 0;
p.A_y = 0;
p.A_z = 0;
p.m_z = 0;

% Off-axis + tilt parameters
p.dx_top    = 0;
p.dy_top    = 0;
p.shift_pow = 1.5;

p.exit_tilt_x = 0;   % radians
p.exit_tilt_y = 0;   % radians
p.tilt_pow   = 1.2; % smooth transition base→exit


p.n_theta = 200;
p.n_z     = 150;

pad = 1.5;   % 50% extra space
%% ===================== FIGURE =====================
fig = figure('Name','pMUT Acoustic Horn Designer',...
    'Color','w',...
    'NumberTitle','off',...
    'Position',[100 100 1200 700]);

ax = axes(fig,'Position',[0.05 0.1 0.55 0.8]);
hold(ax,'on');

%% HARD AXIS + CAMERA LOCK 
set(ax, ...
    'DataAspectRatio',[1 1 1], ...
    'DataAspectRatioMode','manual', ...
    'PlotBoxAspectRatioMode','manual', ...
    'CameraViewAngleMode','manual');

axis(ax,'vis3d');
grid(ax,'on');
view(ax,3);

xlabel(ax,'X (mm)');
ylabel(ax,'Y (mm)');
zlabel(ax,'Z (mm)');

camlight(ax,'headlight');
lighting(ax,'gouraud');
material(ax,'dull');

%% ===================== UI PANEL =====================
panel = uipanel(fig,'Title','Horn Parameters',...
    'FontSize',10,...
    'Position',[0.65 0.1 0.32 0.8]);

ctrl = struct();
editBoxes = struct();  % Store edit box handles
y = 0.92; dy = 0.065;

addSlider('r_base',p.pmut_radius,p.pmut_radius * 3,p.r_base)
addSlider('r_top',0.0002,0.004,p.r_top)
addSlider('height',0.001,0.0031,p.height)
addSlider('p',0.1,10,p.p)
addSlider('q',0.1,10,p.q)
addSlider('dx_top',-0.002,0.002, p.dx_top)
addSlider('dy_top',-0.002,0.002, p.dy_top)
addSlider('exit_tilt_x',-0.5,0.5, p.exit_tilt_x)
addSlider('exit_tilt_y',-0.5,0.5, p.exit_tilt_y)
addSlider('tilt_pow',0,2, p.tilt_pow)
addSlider('shift_pow',0,2, p.shift_pow)



uicontrol(panel,'Style','pushbutton',...
    'String','Export STL',...
    'FontSize',11,...
    'Units','normalized',...
    'Position',[0.2 0.02 0.6 0.07],...
    'Callback',@exportSTL);

%% ===================== INITIAL DRAW =====================
[X,Y,Z] = generateHorn(p);
surfHandle = surf(ax,X,Y,Z,...
    'FaceColor',[0.2 0.6 0.9],...
    'EdgeColor','none',...
    'FaceAlpha',0.35);

% Freeze limits AFTER first render
axis(ax,'tight');
xlim0 = xlim(ax)*pad;
ylim0 = ylim(ax)*pad;
zlim0 = zlim(ax)*pad;

%% ===================== CALLBACKS =====================
    function addSlider(name,minv,maxv,val)
        uicontrol(panel,'Style','text',...
            'String',name,...
            'Units','normalized',...
            'Position',[0.05 y 0.35 0.04],...
            'HorizontalAlignment','left');

        ctrl.(name) = uicontrol(panel,'Style','slider',...
            'Min',minv,'Max',maxv,'Value',val,...
            'Units','normalized',...
            'Position',[0.42 y 0.45 0.04],...
            'Callback',@(src,evt) updateSurface(name));

        editBoxes.(name) = uicontrol(panel,'Style','edit',...
            'String',num2str(val,'%.4g'),...
            'Units','normalized',...
            'Position',[0.89 y 0.1 0.04],...
            'Callback',@(src,~) editSlider(name,src));

        y = y - dy;
    end

    function editSlider(name,src)
        v = str2double(src.String);
        if ~isnan(v)
            ctrl.(name).Value = v;
            updateSurface(name);
        end
    end

    function updateSurface(changedName)
        % Update all parameters from sliders
        f = fieldnames(ctrl);
        for k = 1:numel(f)
            p.(f{k}) = ctrl.(f{k}).Value;
        end
        
        % Update the edit box text for the changed slider
        if nargin > 0 && ~isempty(changedName)
            editBoxes.(changedName).String = num2str(p.(changedName),'%.4g');
        end

        [X,Y,Z] = generateHorn(p);
        set(surfHandle,'XData',X,'YData',Y,'ZData',Z);

        % Re-apply locked limits (IMPORTANT)
        set(ax,'XLim',xlim0,'YLim',ylim0,'ZLim',zlim0);
        drawnow;
    end

    function exportSTL(~,~)
        [X,Y,Z] = generateHorn(p);

        scaleF = 1000; % m → mm
        X = X*scaleF; Y = Y*scaleF; Z = Z*scaleF;

        [F,V] = surf2patch(X,Y,Z,'triangles');
        TR = triangulation(F,V);

        [file,path] = uiputfile('*.stl','Save STL');
        if isequal(file,0), return; end

        stlwrite(TR,fullfile(path,file));
        msgbox('STL exported (mm units)','Success');
    end
end

%% ===================== GEOMETRY =====================
function [X,Y,Z] = generateHorn(p)

z = linspace(0,p.height,p.n_z);
theta = linspace(0,2*pi,p.n_theta);

[X,Y,Z] = deal(zeros(p.n_z,p.n_theta));

% Base and exit normals
n_base = [0;0;1];

nx = sin(p.exit_tilt_y);
ny = -sin(p.exit_tilt_x);
nz = cos(p.exit_tilt_x)*cos(p.exit_tilt_y);
n_exit = [nx; ny; nz];
n_exit = n_exit / norm(n_exit);

% Radius profile
r = p.r_base*(1-z/p.height).^(1/p.p) + ...
    p.r_top*(z/p.height).^(1/p.q);

% Centerline lateral shift
s = (z/p.height).^p.shift_pow;
cx = p.dx_top*s;
cy = p.dy_top*s;

Zmax = -inf;

for i = 1:length(z)

    t = (z(i)/p.height)^p.tilt_pow;

    % Interpolated normal
    n = (1-t)*n_base + t*n_exit;
    n = n / norm(n);

    % Local orthonormal frame
    ref = [1;0;0];
    if abs(dot(ref,n)) > 0.9
        ref = [0;1;0];
    end
    e1 = cross(n,ref); e1 = e1/norm(e1);
    e2 = cross(n,e1);

    for j = 1:length(theta)
        P = [cx(i); cy(i); z(i)] + ...
            r(i)*(cos(theta(j))*e1 + sin(theta(j))*e2);

        X(i,j) = P(1);
        Y(i,j) = P(2);
        Z(i,j) = P(3);

        Zmax = max(Zmax, Z(i,j));
    end
end

% HARD HEIGHT CONSTRAINT (GLOBAL)
Z = Z * (p.height / Zmax);

end