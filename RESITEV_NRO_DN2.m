clear; clc; close all;

%% ---------------------------------------------------------
%  Paths to input files (ABSOLUTE PATHS)
%% ---------------------------------------------------------
dataPath = "";

file_nodes = dataPath + "vozlisca_temperature_dn2_14.txt";
file_cells = dataPath + "celice_dn2_14.txt";
file_size  = dataPath + "velikost_Nx_Ny.txt";

%% ---------------------------------------------------------
%  Read grid size Nx, Ny
%% ---------------------------------------------------------
fid = fopen(file_size, 'r');
Nx = fscanf(fid, "Nx = %d");
Ny = fscanf(fid, "Ny = %d");
fclose(fid);

fprintf("Nx = %d, Ny = %d\n", Nx, Ny);

%% ---------------------------------------------------------
%  Read node data
%% ---------------------------------------------------------
fid = fopen(file_nodes, 'r');

header = fgetl(fid); % "x,y,T"

Nx_f = sscanf(fgetl(fid), "st. koordinat v x-smeri: %d");
Ny_f = sscanf(fgetl(fid), "st. koordinat v y-smeri: %d");
Nv   = sscanf(fgetl(fid), "st. vseh vozlisc: %d");

fprintf("Prebran Nx_f=%d, Ny_f=%d, Nv=%d\n", Nx_f, Ny_f, Nv);

raw = textscan(fid, "%f%f%f", Nv, "Delimiter", ",");
fclose(fid);

x = raw{1};
y = raw{2};
T = raw{3};

fprintf("Prebranih vozlišč: %d\n", length(x));

%% ---------------------------------------------------------
%  Read cell definitions
%% ---------------------------------------------------------
fid = fopen(file_cells, 'r');
fgetl(fid); % "pt1,pt2,pt3,pt4"

Nc = sscanf(fgetl(fid), "st. celic: %d");

C = textscan(fid, "%d%d%d%d", Nc, "Delimiter", ",");
cell_ids = [C{1}, C{2}, C{3}, C{4}];
fclose(fid);

fprintf("Prebranih celic: %d\n", Nc);

%% ---------------------------------------------------------
%  Interpolation point
%% ---------------------------------------------------------
xp = 0.403;
yp = 0.503;

%% ---------------------------------------------------------
%  Method 1: scatteredInterpolant
%% ---------------------------------------------------------
tic;
F1 = scatteredInterpolant(x, y, T, 'linear');
T1 = F1(xp, yp);
time1 = toc;

%% ---------------------------------------------------------
%  Method 2: griddedInterpolant (correct NDGRID)
%% ---------------------------------------------------------
tic;

% Convert linear vectors to NDGRID matrices
Xv = unique(x);
Yv = unique(y);

[Yg, Xg] = ndgrid(Yv, Xv);     % NDGRID is correct for MATLAB

Tgrid = reshape(T, Ny, Nx);    % Ny rows, Nx columns (as given)

F2 = griddedInterpolant(Yg, Xg, Tgrid, "linear");
T2 = F2(yp, xp);

time2 = toc;

%% ---------------------------------------------------------
%  Method 3: Manual bilinear interpolation
%% ---------------------------------------------------------
tic;
T3 = NaN;

for c = 1:Nc
    ids = cell_ids(c, :);

    bl = ids(1); br = ids(2); tr = ids(3); tl = ids(4);

    xmin = x(bl); xmax = x(br);
    ymin = y(bl); ymax = y(tl);

    if xp>=xmin && xp<=xmax && yp>=ymin && yp<=ymax

        T11 = T(bl);
        T21 = T(br);
        T22 = T(tr);
        T12 = T(tl);

        K1 = (xmax-xp)/(xmax-xmin)*T11 + (xp-xmin)/(xmax-xmin)*T21;
        K2 = (xmax-xp)/(xmax-xmin)*T12 + (xp-xmin)/(xmax-xmin)*T22;

        T3 = (ymax-yp)/(ymax-ymin)*K1 + (yp-ymin)/(ymax-ymin)*K2;
        break;
    end
end

time3 = toc;

%% ---------------------------------------------------------
%  Find maximum temperature and its coordinates
%% ---------------------------------------------------------
[Tmax, idx_max] = max(T);
xmaxT = x(idx_max);
ymaxT = y(idx_max);

%% ---------------------------------------------------------
%  Output results
%% ---------------------------------------------------------
fprintf("\n=== REZULTATI ===\n");
fprintf("Method 1 (scatteredInterpolant): T = %.6f, time = %.6f s\n", T1, time1);
fprintf("Method 2 (griddedInterpolant):   T = %.6f, time = %.6f s\n", T2, time2);
fprintf("Method 3 (manual bilinear):      T = %.6f, time = %.6f s\n", T3, time3);

fprintf("\nNajvečja temperatura: %.6f °C\n", Tmax);
fprintf("Koordinate: (%.6f, %.6f)\n", xmaxT, ymaxT);