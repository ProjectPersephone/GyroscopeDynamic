function X_cross = CROSS( X )
% Cross of a Matrix

X_cross = zeros(3,3);
X_cross = [  0    -X(3)   X(2);
           X(3)     0    -X(1);
           -X(2)   X(1)   0
];

end