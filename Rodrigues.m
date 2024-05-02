function R=Rodrigues(w)
    if norm(w)==0, R=eye(3,3);
        else R=eye(3,3)+sin(norm(w))*Skew(w)/norm(w)+(1-cos(norm(w)))*Skew(w)*Skew(w)/(norm(w)^2);
    end
end
