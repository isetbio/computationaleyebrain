function name = lsScarlet(webdir, pattern)
% list files on scarlet server
%
%   function lsScarlet(webdir, [pattern])
%
% Input:
%   webdir:  web directory string
%   pattern: file type pattern string
%   
% Output:
%   name: struct array containing file names and links
%
% Notes:
%   This function is using regular expression to parse the html string.
%   Thus, it highly depend on the php code on the server. If the server
%   code get changed, this might not work. In other words, this is just a
%   quick and dirty solution for getting the list of files on the scarlet
%   web server
%
% Example:
%   webdir = 'http://scarlet.stanford.edu/validation/SCIEN/ISETBIO/VESA';
%   name = lsScarlet(webdir, '.bmp');
%
% (HJ) ISETBIO TEAM, 2015

%% Check inputs
if notDefined('webdir'), error('Web url required'); end
if notDefined('pattern'), pattern = []; end

%% Read and parse html string
p    = '<a[^>]*href="(?<link>[^"]*)">(?<name>[^<]*)</a>';
str  = webread(webdir);
name = regexp(str, p, 'names');

%% Filter by user input pattern
if ~isempty(pattern)
    indx = arrayfun(@(x) ~isempty(strfind(x.name, pattern)), name);
    name = name(indx);
end

end