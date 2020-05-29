notes_update.txt

July 22, 2019:
The mcxyz.c was updated to fix a error:

line 388: replace
	uz = sqrt(1 - ux*ux + uy*uy);
with 
	uz = sqrt(1 - ux*ux - uy*uy);

thanks to Anh Phong Tran, Northeastern University.
