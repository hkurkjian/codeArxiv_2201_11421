! Pour un ensemble de points (omp,q) 
! omp variant dans (ompmin,ompmax)
! q   variant dans (qmin,qmax)
! cherche une solution de det M=0
! si spr=.FALSE., la racine est cherchée dans le prolongement analytique par la fenêtre "sec"
PROGRAM spectreBCS
 USE dspecBCS
 IMPLICIT NONE
 REAL(QP) omplasma,ompmin,ompmax,domp,om,rez,imz,resr(1:4)
 COMPLEX(QPC) z,resc(1:4)
 REAL(QP) q,qmin,qmax,dq,tol(1:2)
 INTEGER iq,iqdep,nq,iomp,iompdep,nomp,sec

 LOGICAL nouveauf,spr,higgs
 INTEGER nn
 CHARACTER(len=90) fichier

 open(10,file='spectreBCS.inp')
  read(10,*)ompmin  !Paramètre de la grille de points (omp,q)
  read(10,*)ompmax
  read(10,*)iompdep
  read(10,*)nomp
  read(10,*)qmin
  read(10,*)qmax
  read(10,*)iqdep
  read(10,*)nq
  read(10,*)rez     !Supposition initiale sur la racine en ompmin,qmin (initialise la recherche par la méthode de Newton)
  read(10,*)imz
  read(10,*)fichier !Nom du fichier de sortie
  read(10,*)spr     !Si .TRUE. on cherche un racine réelle entre 0 et 2. Sinon on prolonge analytiquement par "sec"
  read(10,*)sec
  read(10,*)nouveauf!.TRUE. pour écraser trim(fichier)."dat" s'il existe
  read(10,*)higgs   !.TRUE. pour résoudre M22=0 plutôt que M11*M33-M13**2=0
 close(10)
 
 z=cmplx(rez,imz,QPC)
 om=rez

!Nivo de bavardage
 blaBCS=.TRUE.
 blaBCS=.FALSE.
 blaMin=.TRUE.
 blaZ  =.FALSE.
 blaZ  =.TRUE.

 EPS=1.0e-10_qp

 if(nq==0)then
  dq=0.0
 else
  dq=(qmax-qmin)/nq
 endif
 q=qmin

 if(nomp==0)then
  domp=0.0
 else
  domp=(ompmax-ompmin)/nomp
 endif
 omplasma=ompmin

 write(6,*)'ompmin,ompmax,domp=',ompmin,ompmax,domp
 write(6,*)'qmin,qmax,dq=',qmin,qmax,dq
 write(6,*)'fichier=',trim(fichier)

 if(nouveauf)then
  call system("rm "//trim(fichier)//".dat")
  call system("rm "//trim(fichier)//".info")
  open(10,file=trim(fichier)//".info")
   write(10,*)"!omplasma,qmin,qmax,nq"
   write(10,*)  ompmin,ompmax,nomp,qmin,qmax,nq
  close(10)
 endif

!Boucle sur omplasma
 do iomp=iompdep,nomp
 omplasma=ompmin+domp*iomp
 write(6,*)'--------------------'
 write(6,*)
 write(6,*)"iomp,omplasma=",iomp,omplasma
 write(6,*)
!Boucle sur q
 do iq=iqdep,nq
  q=qmin+dq*iq
  !Precision de la recherche de racines
  tol(1)=4.0_qp*q**2.0_qp*1.e-17_qp
  tol(2)=1.e-20_qp
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)"iq,q,omang(q)=",iq,q,omang(q)
  write(6,*)
  if(spr)then !Racine réelle par ZeroBCS
   call ZeroBCS (q,omplasma,om,tol,resr)
   z=cmplx(om,0.0,QPC)
   resc=cmplx(resr,0.0,QPC)
  else !Racine complexe par ZeroBCSc
   call ZeroBCSc(q,omplasma,sec,z,tol,resc,higgs)
  endif
  write(6,*)
  write(6,*)"q,z=",q,z
 
  open(12,file=trim(fichier)//".dat",POSITION="APPEND")
   write(12,*)q,real(z),imag(z),real(resc),imag(resc)
  close(12)

 enddo  
 enddo  
END PROGRAM spectreBCS
