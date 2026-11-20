function d=det3(M)
%DET3  Determinant of a 3 x 3 matrix, written out (shared by IIPAR and DFX_CORE).
d=M{1,1}.*(M{2,2}.*M{3,3}-M{2,3}.*M{3,2})-M{1,2}.*(M{2,1}.*M{3,3}-M{2,3}.*M{3,1})+M{1,3}.*(M{2,1}.*M{3,2}-M{2,2}.*M{3,1});
end
