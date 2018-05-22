function [] = errorInfo(out,maxg,eta,varargin)

% errorInfo(out,maxg,eta,varargin)
%
% procedure to write error information during iteration either on the
% screen or file
%
% out = 1 -> open a file called errorfile and attend error information,
%   else -> write error information on the screen
% maxg = maximal local error indicator
% eta = global error indicator
% varargin = arbitrary (e.g. 1 if error occurs due to maximal level)
%
% Bettina Schieche, 27th of July 2010

if out
    fid = fopen('errorfile','a'); % open file
else
    fid = 1; % channel for screen
end

if nargin == 5
    fprintf(fid,'Maximal level reached!!\n'); % special error information
end
    
fprintf(fid,'Global error indicator: %5e\n',eta);
fprintf(fid,'Largest local error indicator: %5e\n',maxg);  

if out
    fclose(fid);
end