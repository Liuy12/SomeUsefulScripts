rglsave <- function(filedir, play=T, moviesave=T){
  require('rgl')
  rgl.snapshot(paste(filedir,'png',sep='.'),'png')
  rgl.postscript(paste(filedir,'pdf',sep='.'),'pdf')
  M <- par3d("userMatrix") 
  if(play)
    play3d( par3dinterp( userMatrix=list(M,rotate3d(M, angle=pi,x=0, y=0,z= 1))), duration=10)
  if(moviesave)
    movie3d(par3dinterp( userMatrix=list(M,rotate3d(M, angle=pi,x=0, y=0,z= 1))),duration=5,movie=filedir)
}
