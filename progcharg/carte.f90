! Crée une carte de M pour
! z=re+i*im avec re variant de remin à remax 
!                im variant de immin à immax 
! si im<0 on utilise le prolongement analytique de M à travers la fenêtre "sec"
!
! L'acquisition de points est parallélisée. La commande:
! export OMP_NUM_THREADS=n
! indique au terminal d'utiliser n threads
PROGRAM carte
USE OMP_LIB
USE dspecBCS
IMPLICIT NONE

REAL(QP) :: im,re,q,omp
COMPLEX(QPC) z,det,M(1:4)
REAL(QP) dimm,dre
REAL(QP) immin,immax,remin,remax
INTEGER nim,nre,nimdeb,iim,ire,nredeb,ic,sec
CHARACTER(len=90) suffixe
LOGICAL nvofich

!k 0 to 3
!zk 1 to 6
open(10,file='carte.inp')
 read(10,*)omp
 read(10,*)q
 read(10,*)sec      ! Secteur pour le prolongement analytique
 read(10,*)immin    ! Paramètres de la carte
 read(10,*)immax
 read(10,*)nim,nimdeb
 read(10,*)remin
 read(10,*)remax
 read(10,*)nre,nredeb
 read(10,*)EPS      !dspecBCS/EPS
 read(10,*)blaBCS
 read(10,*)suffixe !nom du fichier de sortie (sans le ".dat")
 read(10,*)nvofich !TRUE to overwrite output files
close(10)


write(6,*)'--------------------'
write(6,*)
write(6,*)'Programme carte'
write(6,*)
write(6,*)'omp,q,omang(q)=',omp,q,omang(q)
write(6,*)'suffixe=',suffixe
write(6,*)'immin=',immin
write(6,*)'immax=',immax
write(6,*)'nim,nimdeb=',nim,nimdeb
write(6,*)'remin=',remin
write(6,*)'remax=',remax
write(6,*)'nre,nredeb=',nre,nredeb
write(6,*)'precision:',EPS


if(nim==0)then
 dimm=0.0
else
 dimm=(immax-immin)/nim
endif

if(nre==0)then
 dre=0.0
else
 dre=(remax-remin)/nre
endif

if(nvofich) call system("rm det"     //trim(suffixe)//".dat")


iim=nimdeb
ire=nredeb
do ire=nredeb,nre
 re=remin+dre*ire
!Parallélisation avec OpenMP
!$OMP  PARALLEL DO &
!$OMP& PRIVATE(im,det,M) &
!$OMP& SHARED(re,ire,immin,remin,dimm,dre) SCHEDULE(DYNAMIC)
 do iim=nimdeb,nim
!$OMP CRITICAL
  im=immin+dimm*iim
  write(6,*)"im,re=",im,re
!$OMP END CRITICAL

 z=cmplx(re,im,kind=qpc)
 call matpro  (z,q,omp,sec,M,det,.FALSE.)
! write(6,*)"Re avt prolong:",real(M)
! write(6,*)"Im avt prolong:",imag(M)
! write(6,*)"Nozieres:   "   ,imag(2*PI*iiq*rhopro(z,q,sec))

!$OMP CRITICAL
 open(23,file="det"     //trim(suffixe)//".dat",POSITION="APPEND")
  write(6,*)"rez,imz,1/abs(det)**2=",re,im,1/abs(det)**2
  write(23,*)omp,q,re,im,1/abs(det)**2
 close(23)
!$OMP END CRITICAL
 enddo
!$OMP END PARALLEL DO
!saut de ligne entre les registres de même re
 open(23,file="det"     //trim(suffixe)//".dat",POSITION="APPEND")
  write(6,*)
  write(23,*)
 close(23)
enddo

END PROGRAM carte
