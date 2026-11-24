clc; clear; close all; repoRoot = fileparts(mfilename('fullpath')); RGB = orderedcolors("gem");

%% Negative electrode invasive-measurements (Graphite, LixC6)

% theoretical capacity
n_atoms     = 6;                        % Number of atoms per formula unit [-]
m_atomic    = 12.01115;                 % Atomic mass [g/mol]
M           = sum(n_atoms.*m_atomic);   % Molecular mass [g/mol]
n           = 1;                        % number of reactive electrons per formula [-]
F           = 96485.3329;               % Faraday constant[sA/mol]
q_th        = n*F/(3.6*M);              % Theoretical capacity of material [mAh/g]

% weight
m_singleCoatingCC       = 0.05303;                                      % x + y [g]
m_doubleSideCoatingCC   = 0.08860;                                      % 2x + y [g]
m_coating               = m_doubleSideCoatingCC - m_singleCoatingCC;    % x = 2x+y - (x+y) [g]
m_CC                    = m_singleCoatingCC - m_coating;                % y [g]

% practical capacity
QsampleNeg              = m_coating*q_th;                               % Sample capacity [mAh]

% current collector density validation (Cu "copper")
L_singleCoatingCC       = 0.09588;                                      % x + y [mm]
L_doubleSideCoatingCC   = 0.18400;                                      % 2x + y [mm]
L_coating               = L_doubleSideCoatingCC - L_singleCoatingCC;    % x = 2x+y - (x+y) [mm] 
L_CC                    = L_singleCoatingCC - L_coating;                % y [mm]
A                       = pi*(18/2/10)^2;                               % Sample area (D = 18mm) [cm^2]
rhoCC                   = m_CC/(A*L_CC/10);                             % Measured density [g/cm^3]
rhoCC_lit               = 8.96;                                         % Literature density [g/cm^3];
relError                = abs(rhoCC-rhoCC_lit)./rhoCC_lit*100           % Relative error [%]

Q_r = QsampleNeg

Q_r_negative_mAh = Q_r; % Retain the result without writing a file.

%% Positive electrode invasive-measurements (NMC811, LiNi0.8Mn0.1Co0.1O2)

% theoretical capacity
n_atoms     = [0.8 0.1 0.1 2];                 % Number of atoms per formula unit [-]
m_atomic    = [58.71 54.938 58.9332 15.9994];  % Atomic mass [g/mol]
M           = sum(n_atoms.*m_atomic);          % Molecular mass [g/mol]
n           = 1;                               % number of reactive electrons per formula [-]
F           = 96485.3329;                      % Faraday constant[sA/mol]
q_th        = n*F/(3.6*M);                     % Theoretical capacity of material

% mass calculation
m_singleCoatingCC       = 0.06853;                                      % x + y [g]
m_doubleSideCoatingCC   = 0.12315;                                      % 2x + y [g]
m_coating               = m_doubleSideCoatingCC - m_singleCoatingCC;    % x = 2x+y - (x+y) [g]
m_CC                    = m_singleCoatingCC - m_coating;                % y [g]

% practical capacity
QsamplePos              = m_coating*q_th;                               % Sample capacity [mAh]

% current collector density diagnostic (assumed aluminum)
L_singleCoatingCC       = 0.08400;                                      % x+y [mm]
L_doubleSideCoatingCC   = 0.14650;                                      % 2x + y [mm]
L_coating               = L_doubleSideCoatingCC - L_singleCoatingCC;    % x = 2x+y - (x+y) [mm] 
L_CC                    = L_singleCoatingCC - L_coating;                % y [mm]
A                       = pi*(18/2/10)^2;                               % Sample area (D = 18mm) [cm^2]
rhoCC                   = m_CC/(A*L_CC/10);                             % Measured density [g/cm^3]
rhoCC_lit               = 2.70;                                         % Literature density [g/cm^3];
relError                = abs(rhoCC-rhoCC_lit)./rhoCC_lit*100           % Relative error [%]

Q_r = QsamplePos

Q_r_positive_mAh = Q_r; % Retain the result without writing a file.
