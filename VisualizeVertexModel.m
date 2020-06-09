function h = VisualizeVertexModel(VertexIDs, Vertices, VertexTypes)
    figure(1)
    close(1)
    h = figure(1);
    hold on;
    
    for i = 1:length(VertexIDs)
        currIDs = VertexIDs{i};
        if (VertexTypes(i) == 1)
            fill(Vertices(currIDs,1)',Vertices(currIDs,2)','c');
        end
        if (VertexTypes(i) == 2)
            fill(Vertices(currIDs,1)',Vertices(currIDs,2)','b');
        end
        if (VertexTypes(i) == 3)
            fill(Vertices(currIDs,1)',Vertices(currIDs,2)','b');
        end
        if (VertexTypes(i) == 4)
            fill(Vertices(currIDs,1)',Vertices(currIDs,2)','g');
        end
    end
   
end