% An example of how to download an image from the VESA data set provided by
% David Hoffman
%

%
% We want a way to do a dir on the webdir
%

webdir = 'http://scarlet.stanford.edu/validation/SCIEN/ISETBIO/VESA/img';

urlwrite(fullfile(webdir,'t_1024x1024Cr_x10.bmp'),'test.bmp');
img = imread('test.bmp'); imshow(img)

