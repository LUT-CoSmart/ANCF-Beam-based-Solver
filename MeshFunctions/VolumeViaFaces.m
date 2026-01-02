function V = VolumeViaFaces(vertices, faces)

    Vsum = 0.0;
    
    c = mean(vertices, 1);
    v = vertices - c; % change to the "center"
    
    % Gauss-Ostrogradsky's theorem over the closed surface
    for k = 1:size(faces,1)
        i0 = faces(k,1); i1 = faces(k,2); i2 = faces(k,3);
        p0 = v(i0,:);    p1 = v(i1,:);    p2 = v(i2,:);
        Vsum = Vsum + dot(p0, cross(p1, p2)) / 6.0;
    end
    
    V = abs(Vsum);  % absolute volume
end