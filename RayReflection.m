%Daniel Adea
%20499515
%Ray Casting
%This script models ray casting in a customizeable 2-d environment 

%clear workspace
clear all; clc; close all;

%Set Bounding Variables

maxBounces = 4;
minIntensity = 0.2;
c = 100;
castStart = 0;
castEnd = 5*pi/4;
castingRadius = castEnd - castStart;   
numRays = 50;
stepSize = castingRadius/numRays;
syms x y

%Set wall locations 
wallSetting = 1;

switch wallSetting
    case 1
        numWalls = 5;
        numCurves = 0;
    %create structure to hold wall data
    walls(1:numWalls) = struct('P1', [], 'P2', [], 'phi', [], 'lambda', [], 'ID', []);
        
    walls(1).P1 = [0, 0];
    walls(1).P2 = [4, 0]; 
    walls(2).P1 = [4, 0];
    walls(2).P2 = [4, 6];
    walls(3).P1 = [4, 6];
    walls(3).P2 = [0, 6];
    walls(4).P1 = [0, 6];
    walls(4).P2 = [0, 0];
    walls(5).P1 = [1, 3.5];
    walls(5).P2 = [4, 3.5];
    
    %set wall reflection coefficients
    walls(1).lambda = 1;
    walls(2).lambda = 1;
    walls(3).lambda = 1;
    walls(4).lambda = 1;
    walls(5).lambda = 1;
    
    case 5
        numWalls = 5;
        numCurves = 0;
    %create structure to hold wall data
    walls(1:numWalls) = struct('P1', [], 'P2', [], 'phi', [], 'lambda', [], 'ID', []);
        
    walls(1).P1 = [0, 0];
    walls(1).P2 = [3, 3]; 
    walls(2).P1 = [3, 3];
    walls(2).P2 = [4, 6];
    walls(3).P1 = [4, 6];
    walls(3).P2 = [-1, 6];
    walls(4).P1 = [-1, 6];
    walls(4).P2 = [0, 0];
    walls(5).P1 = [1, 3.5];
    walls(5).P2 = [4, 3.5];
    
    %set wall reflection coefficients
    walls(1).lambda = 1;
    walls(2).lambda = 1;
    walls(3).lambda = 1;
    walls(4).lambda = 1;
    walls(5).lambda = 1;
    case 4
        numWalls = 8;
        numCurves = 0;
    %create structure to hold wall data
    walls(1:numWalls) = struct('P1', [], 'P2', [], 'phi', [], 'lambda', [], 'ID', []);
        
    walls(1).P1 = [0, 1];
    walls(1).P2 = [4, 0]; 
    walls(2).P1 = [4, 0];
    walls(2).P2 = [4, 6];
    walls(3).P1 = [4, 6];
    walls(3).P2 = [0, 8];
    walls(4).P1 = [0, 8];
    walls(4).P2 = [0, 5];
    walls(5).P1 = [2, 5];
    walls(5).P2 = [4, 3.5];
    walls(6).P1 = [0.1, 7];
    walls(6).P2 = [0.1, 5];
    walls(7).P1 = [0.1, 3];
    walls(7).P2 = [0.1, 1];
    walls(8).P1 = [0, 1];
    walls(8).P2 = [0, 3];
    
    %set wall reflection coefficients
    walls(1).lambda = 1;
    walls(2).lambda = 1;
    walls(3).lambda = 1;
    walls(4).lambda = 1;
    walls(5).lambda = 0.4;
    walls(6).lambda = 0.4;
    walls(7).lambda = 0.6;
    walls(8).lambda = 1; 
    case 3
        numWalls = 5;
        numCurves = 0;
    %create structure to hold wall data
    walls(1:numWalls) = struct('P1', [], 'P2', [], 'phi', [], 'lambda', [], 'ID', []);
        
    walls(1).P1 = [0, 0];
    walls(1).P2 = [4, 0]; 
    walls(2).P1 = [4, 0];
    walls(2).P2 = [4, 6];
    walls(3).P1 = [4, 6];
    walls(3).P2 = [3.5, 7];
    walls(4).P1 = [0, 6];
    walls(4).P2 = [0, 0];
    walls(5).P1 = [1, 3.5];
    walls(5).P2 = [4, 3.5];
    
    %set wall reflection coefficients
    walls(1).lambda = 1;
    walls(2).lambda = 1;
    walls(3).lambda = 1;
    walls(4).lambda = 1;
    walls(5).lambda = 1;    
    case 2    
    numWalls = 5;
    numCurves = 1;
    
    %create structure to hold wall data
    walls(1:numWalls) = struct('P1', [], 'P2', [], 'phi', [], 'lambda', [], 'ID', []);
    %phi is wall angle, lambda is reflection coefficient
    
    curves(1:numCurves) = struct('x1', [], 'x2', [], 'fx', [], 'lambda', [], 'ID', []);

    walls(1).P1 = [0, 0];
    walls(1).P2 = [4, 0]; 
    walls(2).P1 = [4, 0];
    walls(2).P2 = [4, 6];
    walls(3).P1 = [4, 6];
    walls(3).P2 = [4, 8];
    walls(4).P1 = [0, 6];
    walls(4).P2 = [0, 0];
    walls(5).P1 = [0, 6];
    walls(5).P2 = [0, 8];

    curves(1).x1 = 0;
    curves(1).x2 = 4;
    curves(1).fx = @(xx) cos(xx)+7;

    curves(2).x1 = 0;
    curves(2).x2 = 2;
    curves(2).fx = @(xx) cos(xx)+3;

    %set wall reflection coefficients
    walls(1).lambda = 1;
    walls(2).lambda = 1;
    walls(3).lambda = 1;
    walls(4).lambda = 1;
    walls(5).lambda = 1;
    curves(1).lambda = 1;
    curves(2).lambda = 1;
end

%initialize array to be overridden by intersection variables each bounce
intersections = zeros(numWalls, 3);

%axes for curves

%   Find axes for walls
%initialize variables to hold bounds
maxX = walls(1).P1(1);
maxY = walls(1).P1(2);
minX = maxX;
minY = maxY;



%check for smallest/largest values for bounds
for k = 1:numWalls
    xWallPts = [walls(k).P1(1),walls(k).P2(1)];
    yWallPts = [walls(k).P1(2),walls(k).P2(2)];
    for p =1:2
        if xWallPts(p)>maxX
            maxX = xWallPts(p);
        end
        if yWallPts(p)>maxY
            maxY = yWallPts(p);
        end
        if xWallPts(p)<minX
            minX = xWallPts(p);
        end
        if yWallPts(p)<minY
            minY = yWallPts(p);
        end
    end
end

%set up for plotting
plot(0,0, 'MarkerSize',  0.1);
axis([minX-1, maxX+1, minY-1, maxY+1]);
hold on;

%Plot walls and assign the wall angle   
for k = 1:numWalls
    walls(k).phi = atan( (walls(k).P2(2) - walls(k).P1(2)) / (walls(k).P2(1) - walls(k).P1(1)) );
    
    xWallPts = [walls(k).P1(1),walls(k).P2(1)];
    yWallPts = [walls(k).P1(2),walls(k).P2(2)];
    plot(xWallPts, yWallPts, 'b.-', 'Linewidth', 3);
end

%plot curves
for k = 1:numCurves
    fplot(curves(k).fx,[curves(k).x1, curves(k).x2], 'b','Linewidth', 3);
end

%set origin
whichOrig = input('Enter 1 if you would like to choose the origin, 2 if not:');
if whichOrig == 1
    origin = ginput(1);
else    
    origin = [1, 3];
end

%Create an array for the bounces of one ray
currBounces = zeros(maxBounces+2, 3);

%iterate through every angle 
for theta = castStart:stepSize:castEnd
%theta = 0;

%reset variables for a new ray

%whichBounce is the bounce we're on, and we're looking for the point it is
%bouncing to aka at 0 we're at origin looking for where it hits the first
%time
whichBounce = 0;
currWall = 0;
bounceAlive = true;
P3 = origin;
I = 1;
outWindow = false;
    %iterate through the bounces of one ray
    while whichBounce <= maxBounces+1 && bounceAlive 
        if I < minIntensity
            bounceAlive = false;
        end
        
        if outWindow
            currBounces(whichBounce+1,1) = P4(1);
            currBounces(whichBounce+1,2) = P4(2);
            currBounces(whichBounce+1,3) = I;
            plot(currBounces(whichBounce:whichBounce+1,1), currBounces(whichBounce:whichBounce+1,2), 'Color', [0 1-currBounces(whichBounce,3) 0]);
            drawnow;
            %pause(1);
            bounceAlive = false;
        else
        
        
        %reset variables for a new bounce
        P4 = P3 + 100*[cos(theta), sin(theta)];
        tMin = 1;
        wallHit = 0;
        curveHit = 0;
        %check each wall to see which one is bounced on
        for i = 1:numWalls
            intersections(i, 1:3) = intersection(walls(i).P1(1), walls(i).P1(2), walls(i).P2(1), walls(i).P2(2), P3(1), P3(2), P4(1), P4(2));
            if i~=currWall && intersections(i,1) == 1
                if tMin >= intersections(i,3)
                    tMin = intersections(i,3);
                    wallHit = i;
                end
            end
        end
        
        %check each curve to see if one is bounced on 
        for n = 1:numCurves
            m = (P4(2) - P3(2)) / (P4(1) - P3(1));
            fX = @(xx) m*(xx - P3(1)) + P3(2);
            yFun = @(xx) fX(xx) - curves(n).fx(xx);
%            x = fzero(yFun, [P3(1) P3(2)]);
            try

            if P4(1) < P3(1)
                [sol, err, flag]= fsolve(yFun, P4(1));
            else
                [sol, err, flag]= fsolve(yFun, P3(1));    
            end

                if flag> 0        

                x = sol;
                y = fX(x);
                else
                    [sol, err, flag]= fsolve(yFun, 0);
                        x= sol;
                        y = fX(x);
                end
             catch ME
                 x= 0 ;
                  y = 0;
            end
            
             if (x < curves(1).x1 && x > curves(1).x2) || (x > curves(1).x1 && x < curves(1).x2)
                tDist = sqrt((x-P3(1))^2 + (y-P3(2))^2);
                totDist = sqrt((P4(1)-P3(1))^2 + (P4(2)-P3(2))^2);
                t = tDist/totDist;
                if tMin>t && t>1*10^-6 
                    tMin = t;
                    wallHit = 0;
                    curveHit = n; 
                end
            end               
        end
       
        %get values for bounce; then plot
            currBounces(whichBounce+1,1) = P3(1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           currBounces(whichBounce+1,1) = P3(1);
            currBounces(whichBounce+1,2) = P3(2);
            currBounces(whichBounce+1,3) = I;
            if whichBounce ~= 0
                plot(currBounces(whichBounce:whichBounce+1,1), currBounces(whichBounce:whichBounce+1,2), 'Color', [0 1-currBounces(whichBounce,3) 0]);
                drawnow;
                %pause(1);
            end 
            whichBounce = whichBounce +1;
        %if hit wall
        if wallHit ~= 0

            %set up for next bounce
            theta = -theta + 2*walls(wallHit).phi + rand(1)/60;
            I = I*walls(wallHit).lambda;
            P3 = walls(wallHit).P1 + intersections(wallHit, 2) * ( walls(wallHit).P2 - walls(wallHit).P1 );
            currWall = wallHit;
        %if hit curve
        elseif curveHit ~=0
            
            %set up for next bounce
            h = 1*10^-6;
            curvePhi = (curves(curveHit).fx(x+h) - curves(curveHit).fx(x-h))/(2*h);    
            theta = -theta + 2*curvePhi;
            I = I*curves(curveHit).lambda;

            P3 = [x y];    
            currWall = 0;
        %if out window
        else
            outWindow = true; 


        end  
        end
    end
%plot each ray    
%plot(currBounces(1:whichBounce,1), currBounces(1:whichBounce,2));
%drawnow;
%pause(2)
end          
hold off
        
        