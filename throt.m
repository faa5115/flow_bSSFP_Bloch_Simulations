function Rth=throt(phi,theta) 
 Rth=zrot(theta)*xrot(phi)*zrot(-theta);