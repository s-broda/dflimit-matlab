function A=adj3n(M)
%ADJ3N  Adjugate of a 3 x 3 matrix with the sign convention of TR_LIMDIST (shared by IIPAR and DFX_CORE).
A=[M(2,2)*M(3,3)-M(2,3)*M(3,2), -(M(1,2)*M(3,3)-M(1,3)*M(3,2)), M(1,2)*M(2,3)-M(1,3)*M(2,2);
  -(M(2,1)*M(3,3)-M(2,3)*M(3,1)), M(1,1)*M(3,3)-M(1,3)*M(3,1), -(M(1,1)*M(2,3)-M(1,3)*M(2,1));
   M(2,1)*M(3,2)-M(2,2)*M(3,1), -(M(1,1)*M(3,2)-M(1,2)*M(3,1)), M(1,1)*M(2,2)-M(1,2)*M(2,1)];
end
