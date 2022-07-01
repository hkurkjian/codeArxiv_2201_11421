! Prend un suite de valeur de M et des fonctions de réponse en om+i0^+ pour omega variant de remin à remax
! Le résultat est écrit dans un fichier trim(suffixe)//".dat"
PROGRAM coupe
USE dspecBCS
IMPLICIT NONE

REAL(QP) :: re,q,omp
COMPLEX(QPC) det,M(1:4)
REAL(QP) detr,Mr(1:4)
REAL(QP) dre
REAL(QP) remin,remax
INTEGER nre,ire,nredeb
CHARACTER(len=90) suffixe
LOGICAL nvofich

!k 0 to 3
!zk 1 to 6
open(10,file='coupe.inp')
 read(10,*)omp
 read(10,*)q
 read(10,*)remin    !Paramètres de la coupe
 read(10,*)remax
 read(10,*)nre,nredeb
 read(10,*)EPS      !dspecBCS/EPS
 read(10,*)blaBCS
 read(10,*)suffixe  !nom du fichier de sortie (sans le ".dat")
 read(10,*)nvofich  !TRUE to overwrite output files
close(10)


write(6,*)'--------------------'
write(6,*)
write(6,*)'Programme coupe'
write(6,*)
write(6,*)'omp,q=',omp,q
write(6,*)'suffixe=',suffixe
write(6,*)'remin=',remin
write(6,*)'remax=',remax
write(6,*)'nre,nredeb=',nre,nredeb
write(6,*)'precision:',EPS


if(nre==0)then
 dre=0.0
else
 dre=(remax-remin)/nre
endif

if(nvofich) call system("rm "//trim(suffixe)//".dat")


ire=nredeb
do ire=nredeb,nre
 re=remin+dre*ire
 write(6,*)"re=",re

 call mat_r  (re,q,omp,Mr,detr)
 M=Mr-iiq*PI*rho(re,q)
 det=M(1)*M(3)-M(4)**2

 open(23,file=trim(suffixe)//".dat",POSITION="APPEND")
  write(6,*)"rez,1/abs(det)**2=",re,1/abs(det)**2
  write(23,*)omp,q,re,1/abs(det)**2,real(M(3)/det),imag(M(3)/det),real(1/M(2)),imag(1/M(2))
 close(23)
enddo

END PROGRAM coupe
