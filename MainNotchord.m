close all;

origHex = [[1, 0];[1/2,sqrt(3)/2];[-1/2, sqrt(3)/2];[-1 0];[-1/2 -sqrt(3)/2];[1/2 -sqrt(3)/2]];

Vertices = [];
VertexIDs = {};

rowsize = 1;
columnsize = 30;

IDTypes = [];

%Make row of hexagons
for i = 1:columnsize
    rightTranslation = (i - 1)*3;
    newVertices = origHex + [rightTranslation 0];
    Vertices = [Vertices;newVertices];
    VertexIDs = [VertexIDs, {(1:6)+(i-1)*6}];
end

%Fill in gaps between gaps using a shifted row
for i = 1:length(Vertices)
    Vertices = [Vertices;Vertices(i,:)+[1.5, -sqrt(3)/2]];
    if (mod(i,6) == 0)
       VertexIDs = [VertexIDs {VertexIDs{end}+6}];
    end
end

%Replicate row pair for many rows
for i = 1:rowsize
   Vertices = [Vertices;[Vertices(1:360,:) + [0, -sqrt(3)]*i]];
   for j = 1:columnsize*2
      VertexIDs = [VertexIDs {VertexIDs{end}+6}];
   end
end

%Join vertices sharing a vertex into a single vertex
for i = 1:length(VertexIDs)
    currIDs = VertexIDs{i};
    for j = 1:length(currIDs)
        currVert = Vertices(currIDs(j),:);
        pos = find(abs(sum(Vertices - currVert,2)) < 10e-8);
        if (~isempty(pos))
            VertexIDs{i}(j) = pos(1);
        end
    end
end

v = VideoWriter('D:\Vertex Model\vm_notochord.avi');
open(v);

lambda = 30;
a0 = 3*sqrt(3)/2;
c0 = 2*sqrt(pi*a0);
beta  = 5;
rng(1)

gamma = 1;
nu = 1;
deltat = 0.000005;
gammarandom = 3;
hi = 1;

T1min = 0.1;

AllVels = {};

h = VisualizeVertexModel(VertexIDs,Vertices);

as = zeros(length(VertexIDs),1);
for i = 1:length(as)
    if (rand(1) < 0.85)
        as(i) = a0/6;
    else
        as(i) = a0;
    end
end

frames = [];
for step = 1:1000
    disp(string(step));
    oldPositions = Vertices;
    forces = zeros(size(Vertices,1),size(Vertices,2));
    areas = zeros(length(VertexIDs),1);
    for i = 1:length(VertexIDs)
        areas(i) = area(polyshape(Vertices(VertexIDs{i},1),Vertices(VertexIDs{i},2)));
    end
    
    perims = zeros(length(VertexIDs),1);
    for i = 1:length(VertexIDs)
        perims(i) = perimeter(polyshape(Vertices(VertexIDs{i},1),Vertices(VertexIDs{i},2)));
    end
    
    for i = 1:length(VertexIDs)
        
        currIDs = VertexIDs{i};
        for j = 1:length(currIDs)
            k = currIDs(j);
            currVert = Vertices(k,:);
            if (j == length(currIDs))
                nextInd = 1;
            else
                nextInd = j+1;
            end
            if (j == 1)
                prevInd = length(currIDs);
            else
                prevInd = j-1;
            end
            
            nextVert = Vertices(currIDs(nextInd),:);
            prevVert = Vertices(currIDs(prevInd),:);
            
            deltaiAl = 1/2*[nextVert(2)-prevVert(2),prevVert(1)-nextVert(1)];
            
            a = as(i);
            
            AreaTerm = 2*lambda*deltaiAl*(areas(i)-a);   
            
            deltaidl1 = 1/norm(prevVert - currVert) * (currVert - prevVert);
            deltaidl2 = 1/norm(nextVert - currVert) * (currVert - nextVert);
            
            PerimeterTerm = 2*beta*(perims(i) - 2*sqrt(pi*a))*(deltaidl1+deltaidl2);

            ContractileTerm = gamma*(deltaiAl)*areas(i);
            
            RandomTerm = gammarandom * randn(1,2);
                        
            forces(k,:) = forces(k,:) - AreaTerm - PerimeterTerm - RandomTerm;
            
        end
    end
        
    Vertices = Vertices + forces*deltat/nu;  
    
    %Implement T1 Vertex Transitions
    prevT1s = [0 0];

    for i = 1:length(VertexIDs)
        currIDs = VertexIDs{i};
        for j = 1:length(currIDs)
            if (j == length(currIDs))
                next = 1;
            else
                next = j + 1;
            end
            dist = norm(Vertices(currIDs(j),:) - Vertices(currIDs(next),:));
            
            if (dist < T1min && sum([VertexIDs{:}] == currIDs(j)) >= 2 && sum([VertexIDs{:}] == currIDs(next)) >= 2)
                if (any(sum([currIDs(j) currIDs(next)] == prevT1s, 2) == 2) || any(sum([currIDs(next) currIDs(j)] == prevT1s, 2) == 2))
                    continue;
                end
                                
                cells2points = [];
                cells1point = [];
                for m = 1:length(VertexIDs)
                    tempIDs = VertexIDs{m};
                    if (sum(ismember([currIDs(j) currIDs(next)], tempIDs)) == 2)
                        cells2points = [cells2points m];
                    elseif (sum(ismember([currIDs(j) currIDs(next)], tempIDs)) == 1)
                        cells1point = [cells1point m];
                    end
                end
                
                p1x = Vertices(currIDs(j),1);
                p1y = Vertices(currIDs(j),2);
                
                p2x = Vertices(currIDs(next),1);
                p2y = Vertices(currIDs(next),2);
                
                mix = mean([p1x p2x]);
                miy = mean([p1y p2y]);
                
                slope = (p2y-p1y) / (p2x-p1x);
                slope = -1/slope;
                
                ydist = abs(p1y - p2y);
                xvalues = [mix - ydist/2 * hi, mix + ydist/2 * hi];
                
                yvalues = slope * (xvalues - mix) + miy;
                
                Vertices = [Vertices;[xvalues(1) yvalues(1)]];
                Vertices = [Vertices;[xvalues(2) yvalues(2)]];
                
                title(string(step))
                for m = 1:length(cells2points)
                    
                    tempIDs = VertexIDs{cells2points(m)};
                    
                    jindex = find(tempIDs == currIDs(j));
                    nextindex = find(tempIDs == currIDs(next));
                                        
                    tempIDs1 = tempIDs;
                    tempIDs1(jindex) = length(Vertices) - 1;
                    tempIDs1(nextindex) = [];

                    tempIDs2 = tempIDs;
                    tempIDs2(jindex) = length(Vertices);
                    tempIDs2(nextindex) = [];
                    
                    [in, on] = inpolygon(Vertices(length(Vertices)-1,1),Vertices(length(Vertices)-1,2),Vertices(tempIDs,1),Vertices(tempIDs,2));
                    if (in || on)
                        tempIDs = tempIDs1;
                    else
                        tempIDs = tempIDs2;
                    end
%                     if (perimeter(polyshape(Vertices(tempIDs1,1),Vertices(tempIDs1,2))) < perimeter(polyshape(Vertices(tempIDs2,1),Vertices(tempIDs2,2))))
%                         tempIDs = tempIDs1;
%                     else
%                         tempIDs = tempIDs2;
%                     end
            
                    VertexIDs{cells2points(m)} = tempIDs;
                end
                
                for m = 1:length(cells1point)
                    
                    tempIDs = VertexIDs{cells1point(m)};
                    
                    pointindex = find((tempIDs == currIDs(j)) + (tempIDs == currIDs(next)));
                    
                    if (pointindex == 1)
                        prevInd = length(tempIDs);
                    else
                        prevInd = pointindex - 1;
                    end
                    
                    if (pointindex == length(tempIDs))
                        nextInd = 1;
                    else
                        nextInd = pointindex + 1;
                    end
                    
                    for n = 1:length(cells2points)
                        if (ismember(tempIDs(prevInd), VertexIDs{cells2points(n)}))
                            if (ismember(length(Vertices) - 1, VertexIDs{cells2points(n)}))
                                tempIDs = [tempIDs(1:pointindex - 1) length(Vertices) - 1 length(Vertices) tempIDs(pointindex + 1:end)];
                                break
                            elseif (ismember(length(Vertices), VertexIDs{cells2points(n)}))
                                tempIDs = [tempIDs(1:pointindex - 1) length(Vertices) length(Vertices) - 1 tempIDs(pointindex + 1:end)];
                                break
                            end
                        end
                        
                        if (ismember(tempIDs(nextInd), VertexIDs{cells2points(n)}))
                            if (ismember(length(Vertices) - 1, VertexIDs{cells2points(n)}))
                                tempIDs = [tempIDs(1:pointindex - 1) length(Vertices) length(Vertices) - 1 tempIDs(pointindex + 1:end)];
                                break
                            elseif (ismember(length(Vertices), VertexIDs{cells2points(n)}))
                                tempIDs = [tempIDs(1:pointindex - 1) length(Vertices) - 1 length(Vertices) tempIDs(pointindex + 1:end)];
                                break
                            end
                        end
                    end
                    VertexIDs{cells1point(m)} = tempIDs;
                    
                end
                
                prevT1s = [prevT1s;[length(Vertices)-1 length(Vertices)]];
                prevT1s = [prevT1s;[length(Vertices) length(Vertices)-1]];
                
                break;                        
            end
        end
    end
    
    h = VisualizeVertexModel(VertexIDs,Vertices);
    frame = getframe(h);
    writeVideo(v, frame);
end

close(v);



