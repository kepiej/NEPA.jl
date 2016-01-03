# DEA is a special case of the DDF where gx or gy is zero depending on the orientation
DEA = (X,Y,input,RStype::RS) -> input ? DDF(X,Y,X,zeros(size(Y)),RStype) : DDF(X,Y,zeros(size(X)),Y,RStype)

# Define shorthand functions depending on the returns to scale
DEA_CRS = (X,Y,input) -> DEA(X,Y,input,CRS())
DEA_VRS = (X,Y,input) -> DEA(X,Y,input,VRS())
DEA_NIRS = (X,Y,input) -> DEA(X,Y,input,NIRS())
DEA_NDRS = (X,Y,input) -> DEA(X,Y,input,NDRS())
