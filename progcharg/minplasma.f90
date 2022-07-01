!En fonction de omplasma (variant de ompmin à ompmax), cherche le minimum de dispersion
!de la branche de plasmons (voir figure 2 dans [3])
PROGRAM minplasma
 USE dspecBCS
 IMPLICIT NONE
 REAL(QP) :: omplasma
 INTEGER nn
 CHARACTER(len=90) fichier
 REAL(QP) ompmin,ompmax,domp,tolx
 REAL(QP) zm,qm
 INTEGER iomp,iompdep,nomp
 LOGICAL nouveauf

 open(10,file='minplasma.inp')
  read(10,*)ompmin
  read(10,*)ompmax
  read(10,*)nomp
  read(10,*)qm
  read(10,*)zm
  read(10,*)fichier
  read(10,*)nouveauf
 close(10)

 blaBCS=.TRUE.
 blaBCS=.FALSE.
 blaMin=.TRUE.
 blaZ  =.FALSE.
 blaZ  =.FALSE.

 EPS=1.0e-12_qp

 if(nomp==0)then
  domp=0.0
 else
  domp=(ompmax-ompmin)/nomp
 endif

 write(6,*)'ompmin,ompmax,domp=',ompmin,ompmax,domp
 write(6,*)'fichier=',trim(fichier)

 if(nouveauf)then
  call system("rm "//trim(fichier)//".dat")
  call system("rm "//trim(fichier)//".info")
  open(10,file=trim(fichier)//".info")
   write(10,*)"!omplasma,qm,zm"
   write(10,*)  omplasma,qm,zm
  close(10)
 endif

 do iomp=0,nomp
  omplasma=ompmin+domp*iomp
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)"iomp=",iomp
  write(6,*)
  call mindisp(omplasma,qm,zm)
  write(6,*)
  write(6,*)"omplasma,qm,zm=",omplasma,qm,zm
 
  open(12,file=trim(fichier)//".dat",POSITION="APPEND")
   write(12,*)omplasma,qm,zm
  close(12)

 enddo
END PROGRAM minplasma
