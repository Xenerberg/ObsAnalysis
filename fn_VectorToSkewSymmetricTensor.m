function SkewSymmetricTensor = fn_VectorToSkewSymmetricTensor(v)
    %SkewSymmetricTensor = zeros(3,3);
    SkewSymmetricTensor(1,1) = sym(0);
    SkewSymmetricTensor(2,2) = sym(0);
    SkewSymmetricTensor(3,3) = sym(0);
    SkewSymmetricTensor(1,2) = -v(3);
    SkewSymmetricTensor(1,3) = v(2);
    SkewSymmetricTensor(2,3) = -v(1);
    SkewSymmetricTensor(2,1) = v(3);
    SkewSymmetricTensor(3,1) = -v(2);
    SkewSymmetricTensor(3,2) = v(1);
end