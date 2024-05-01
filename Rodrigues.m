function R = Rodrigues(w)
    if norm(w)==0, R=eye(3,3);
        else R=eye(3,3)+sin(norm(w))*hat(w)/norm(w)+(1-cos(norm(w)))*hat(w)*hat(w)/(norm(w)^2);
    end
end