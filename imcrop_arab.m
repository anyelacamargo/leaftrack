% This program segments Arabidopsis rosettes
% RGB image is converted into HSV channel
% Five regions of interest are selected
% Segmentation is done by ROI
% an output file is produced
% Anyela Camargo, Feb 2016.

%Function main
function main()
  
 rootname =  'M:\anyela\repo\leaftrack';
 outputfile = 'M:\anyela\repo\result.txt';
 cv = 1;
 if cv == 0 
    cv = getCoordinates();
    save('coordinatesAT12Big.mat', 'cv' );
 end
 load('coordinatesAT12Big.mat');
 rosetteArea(rootname, outputfile, cv)
  

%Segment image
% function rosetteArea(rootname, outputfile)
%   rootname   = folder where images are stored
%   outputfile = file where results will be written
function rosetteArea(rootname, outputfile, cv)  
    vl = {};
    % ROI by pot in six pots tray. Two rows and three columns. 
    % Middle pot first row is empty
    
    rd = dir(rootname);
    
    fileID = fopen(outputfile,'w');
    fprintf(fileID,'%s, %s, %s, %s, %s\n', 'Filename', 'tray', 'loc', 'area', 'Excen');
    %b = [44,46,48,52,58,60,67,70,71,73,75,77,79,82,83,85,86,87];
    %b = [74, 87];
    for i=4:length(rd)-1
        z = 1;
        name0 = rd(i).name;
        %sd1 = dir(strcat(rootname, name0));
        fname = strcat(rootname, '\', name0);
        I= imread(fname);
        %Filename
        strcat(rootname, name0)
        %Number of pots in tray
        C = zeros(size(I(:,:,1)));
        
        for x=1:length(cv)
            [BWHO, a] = processImage(I,cv{x}, x,  3, 0.52, 0.99); 
            C = C + BWHO;
            for j=1:length(a)
                savedata(fileID, i, x, a(j).Area, a(j).Eccentricity, fname)
            end
            close all;               
        end
        %imshow(C);
    end
    fclose(fileID)

    
 function[BWHO,a] = processImage(I, loc, pos, channel, min, max) 
    HSV = rgb2hsv(I);
    ISEG = roiEllipse(HSV(:,:,channel), loc);
    BWH = roicolor(ISEG, min,max);
    BW2 = imfill(BWH,'holes');
    BWHO = bwareaopen(BW2, 6000); % 4500 was before
    %imshow(BWHO);
    a = filterDataAll(BWHO);
    %imshow(X);
    % Plot output images
    %figure, imshow(I), hold on
    %himage = imshow(X);
    %set(himage, 'AlphaData', 0.5);
    %title(sprintf('N = %s', pos));
            
    
% Save data in file
% function savedata(filename, loc, cc, fname)
%   fname   = File name where segmented image is
%   loc      = Location of pot in Tray
%   c        = Area of segmented area
%   filename = Output file where results are stored
function savedata(outputfile, tray, loc, desc1, desc2, fname)
  fprintf(outputfile,'%s,%2d,%2d, %12.0f, %2.4f\n', fname, tray, loc, desc1, desc2);
    
    
%Filter out missclassified pixels
% function[X,a]
%   X =  Binary matrix, 1 is selected area  
%   a = area of selected object
function[X,a] = filterData(BWHO)
    b = [];
    cc = regionprops(BWHO, 'Area', 'Eccentricity');
    [a,v] = sort([cc.Area]);
    for i=length(v):-1:(length(v)-1)
        i, cc(v(i)).Eccentricity;
        if(cc(v(i)).Eccentricity < 0.80)
            b = v(i);
            break;
        end
    end
    CC = bwconncomp(BWHO);
    L = labelmatrix(CC);
    %[a,b] = max([cc.Area]);
    X = zeros(size(L));
    o = find(L == b);
    a = cc(b).Area;
    X(o) = 1;
   
function[cc] = filterDataAll(BWHO)
    %b = [];
    cc = regionprops(BWHO, 'Area', 'Eccentricity');
    
 
 %Select target area
 % I    = Image 
 % loc  = Target location in Image I
 % ISEG = Selected image with binary background
 
function[ISEG] = roiEllipse(I, loc)
    h_im = imshow(I);
    e = imellipse(gca,loc);
    BW = createMask(e,h_im);
    ISEG = I.*BW;
  
    
function[ISEG] = roiEllipseManueal(I, loc)
    cv = {};
    v1 = [510.588418430885 495.937110834371 419.905354919053 432.707347447073];
    v2 = [1499.13211678832 497.241605839416 414.201459854015 414.201459854015];
    v3 = [505.648905109489 1202.58467153285 414.201459854015 414.201459854015];
    v4 = [994.886861313869 1193.5802919708 414.201459854015 414.201459854015];
    v5 = [1508.13649635036 1196.58175182482 414.201459854015 414.201459854015];
    cv = [cv, v1, v2, v3, v4, v5];
    h_im = imshow(I);
    e = imellipse(gca,cv{1});
    BW = createMask(e,h_im);
    ISEG = I.*BW;
    

% Select ellipse from Image I
function selectROI(I)
    imshow(I);
    h = imellipse;
    % Z:\lemnatec_data\by_experiment\AT7\AT7-02211\2014-01-31\00_VIS_tv_000-0-0-0.png;


 % Get coordinates of pots
function[cv] = getCoordinates()
    [FileName,PathName,FilterIndex] = uigetfile();
    fname = strcat(PathName, FileName);
    I = imread(fname);
    cv = {};
    R=1;
    while R~=0
        imshow(I);
        h = imellipse
        wait(h);
        cv = [cv, getPosition(h)];
        R = input('Another image?')
    end
    save( 'filename', 'cv' );
    
    
    I = imread('K:\Phendata\arabidopsis\AT12\DSC_0022.JPG');
    c = 1;
    HSV = rgb2hsv(I);
    BWH = roicolor(HSV(:,:,c), 0.19,0.25);
    imshow(BWH);
    ISEG = roiEllipse(BWH, cv{5});
    BW2 = imfill(ISEG,'holes');
    cc = regionprops(BW2, 'Area', 'Eccentricity')
    
    
    
    
    