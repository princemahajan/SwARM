function [ GMST, ERS ] = JD2GMST( JD_UTC )
%JD2GMST Converts given Julian Date to Greenwich Meast Sidereal Time
%(radians) and Earth Rotational Speed (radians per sec).
% Author: Bharat Mahajan (https://github.com/princemahajan)

% Earth Orientation Parameters Data
% Copied from https://datacenter.iers.org/eop/-/somos/5Rgv/latest/214 on
% 01/01/2000 and 05/22/2016 0h UTC
% FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2F12.6)
%**********************************************************************************
%  
%      Date      MJD      x          y        UT1-UTC       LOD         dX        dY        x Err     y Err   UT1-UTC Err  LOD Err     dX Err       dY Err  
%                         "          "           s           s          "         "           "          "          s         s            "           "
%     (0h UTC)
% 2000   1   1  51544   0.043242   0.377915   0.3555094   0.0009762  -0.000108  -0.000103   0.000141   0.000158  0.0000111  0.0000244    0.000023    0.000047
% 2016   5  22  57530   0.071177   0.493709  -0.1708580   0.0014622   0.000159   0.000193   0.000042   0.000037  0.0000078  0.0000180    0.000038    0.000045
%

UT1_UTC = 0.3555094; % seconds

% UT1
JD_UT1 = JD_UTC + UT1_UTC/86400; % days

% Julian centuries elapsed since J2000.0 (See Vallado topic: Sidereal Time)
J2000_EPOCH = 2451545; % Jan 1, 2000 12:00:00 TT
Tut1 = (JD_UT1 - J2000_EPOCH)/36525; % centuries

% Greenwich Mean Sidereal Time in seconds
GMST = 67310.54841 + (876600*3600 + 8640184.812866)*Tut1 + 0.093104*Tut1^2 - 6.2e-6*Tut1^3;

% convert to radians
GMST = mod(mod(GMST, 86400)*2*pi/86400,2*pi);

% Earth Rotational speed in sidereal day per solar day (See Vallado topic: solar time and universal time)
ERS = 1.002737909350795 + 5.9006e-11*Tut1 - 5.9e-15*Tut1^2;

% ERS in radians per sec
ERS = ERS*2*pi/86400;

end

