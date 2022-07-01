!Unités: z -> z/Delta, om -> om/Delta, q -> \bar{q}=k_F*q/2m*Delta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE dspecBCS
 USE nrtype
 USE nrutil, ONLY : ecrit
 USE modsim
 USE recettes
 IMPLICIT NONE
 LOGICAL blaBCS,blaZ,blaMin   !verbose level
 REAL(QP) EPS !integral precision
 REAL(QP), PARAMETER, DIMENSION(1:4) :: s=(/1.0_qp,1.0_qp,1.0_qp,-1.0_qp/) !signe dans l’intégrale sur omega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE mindisp(omplasma,q,zm) !Extrait le minimum zm<omplasma atteint en q
!de la branche de plasmons lorsque alpha<0
!Utilise la méthode de Newton
!q doit être initialisé (q=0.1 est en général un bon choix) 
 USE recettes, ONLY : rtsafe
 IMPLICIT NONE
 REAL(QP), INTENT(IN) :: omplasma
 REAL(QP), INTENT(INOUT) :: q
 REAL(QP), INTENT(INOUT) :: zm
 
 REAL(QP), DIMENSION(1:1) :: qdep
 REAL(QP) tolx

 if(blaMin)then
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)'minder'
  write(6,*)'omplasma=',omplasma
  write(6,*)
  write(6,*)'--------------------'
 endif
 
 tolx=1.e-11_qp
 qdep=(/q/)
 call mnewt(20,qdep,tolx,1.e-12_qp,dzqdq) !On cherche une racine de la dérivée dz_q/dq
 
 q=qdep(1)

 if(blaMin)then
  write(6,*)'q,omplasma,dzm=',q,omplasma,zm-omplasma
 endif
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dzqdq(qm,dzq,ddzq)
!Calcule dz_q/dq et d^2z_q/dq^2 pour la méthode de Newton
  IMPLICIT NONE
  REAL(QP), DIMENSION(:), INTENT(IN) :: qm
  REAL(QP), DIMENSION(:), INTENT(OUT) :: dzq
  REAL(QP), DIMENSION(:,:), INTENT(OUT) :: ddzq

  REAL(QP) zqP,zqM,dq,tol(1:2),resbid(1:4)
  
  if(blaZ)then
   write(6,*)'----------------------'
   write(6,*)'q=',qm(1)
  endif
  
  dq=qm(1)*1.0e-5

  zqP=zm
  zqM=zm
  tol(1)=4.0_qp*q**2.0_qp*1.e-15_qp
  tol(2)=q*1.e-12_qp

  call ZeroBCS(qm(1)   ,omplasma,zm ,tol,resbid)
  write(6,*)'q,zq=',qm(1),zm
  call ZeroBCS(qm(1)+dq,omplasma,zqP,tol,resbid)
  write(6,*)'qP,zP=',qm(1)+dq,zqP
  call ZeroBCS(qm(1)-dq,omplasma,zqM,tol,resbid)
  write(6,*)'qM,zM=',qm(1)-dq,zqM
  
  dzq =(zqP-zqM)/(2*dq)
  ddzq=(zqP+zqM-2*zm)/dq**2
  if(blaMin)then
   write(6,*)'q,zq,dzq,ddzq=',qm(1),zm,dzq,ddzq
   write(6,*)
  endif
  END SUBROUTINE dzqdq
 END SUBROUTINE mindisp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE ZeroBCS(q,omplasma,om,tol,residu)
! En fonction de q et omplasma, trouve la racine **réelle** om<2*Delta à l'équation M11*M33-M13^2=0
! Utilise la méthode de newton, om doit être initialisé (ce qui peut être non trivial)
! tol=(/tolx,tolf/) erreur tolérée sur om et f(om)=M11*M33-M13^2 respectivement
! residu: le residu de chi=M^-1
! residu(1)=lim_{z->om} (z-om)*chi_{11}
! residu(2)=lim_{z->om} (z-om)*chi_{33}
! residu(3)=lim_{z->om} (z-om)*chi_{13}
 USE recettes, ONLY : mnewt
 IMPLICIT NONE
 REAL(QP), INTENT(IN) :: q,omplasma,tol(1:2)
 REAL(QP), INTENT(INOUT) :: om
 REAL(QP), INTENT(OUT) :: residu(1:4)

 REAL(QP) :: M(1:4),dM(1:4),ddet
 REAL(QP), DIMENSION(1:1) :: omdep
 REAL(QP) tolx
 
 
 if(blaZ)then
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)'ZeroBCS'
  write(6,*)
  
  write(6,*)'--------------------'
  write(6,*)'q=',q
  write(6,*)'om=',om
 endif
 
 omdep=(/om/)
 call mnewt(20,omdep,tol(1),tol(2),fdet)

 residu(1)= M(3)/ddet
 residu(2)= 1/M(2)
 residu(3)= M(1)/ddet
 residu(4)=-M(4)/ddet
 
 om=omdep(1)
 if(blaZ)then
  write(6,*)'--------------------'
  write(6,*)'omfinal=',om
  write(6,*)'M=',M
  write(6,*)'residu=',residu
 endif
 
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fdet(oom,det,Jdet)
! Calcule det=M11*M33-M13^2 et sa dérivée par rapport à omega pour la méthode de Newton
! Partage les valeurs de M et ddet avec la sousroutine ZeroBCS parente
  IMPLICIT NONE
  REAL(QP), DIMENSION(:), INTENT(IN) :: oom
  REAL(QP), DIMENSION(:), INTENT(OUT) :: det
  REAL(QP), DIMENSION(:,:), INTENT(OUT) :: Jdet
  
  if(blaZ)then
   write(6,*)'----------------------'
   write(6,*)'oom=',oom
  endif
  
  call mat_r  (oom(1),q,omplasma,M,det(1))
  call dmat_r (oom(1),q,omplasma,dM)
  
  ddet=dM(1)*M(3)+dM(3)*M(1)-2*dM(4)*M(4)
  Jdet=ddet
  if(blaZ)then
   write(6,*)'om,det,Jdet=',oom(1),det,Jdet
  endif
  END SUBROUTINE fdet
 END SUBROUTINE ZeroBCS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE mat_r(z,q,omplasma,M,det)
! En fonction de z,q et omplasma, calcule la matrice M des fluctuations densité-Delta (et det, son déterminant). 
! par intégration sur l'énergie omega.
! Routine spécialisée pour z réel, qui renvoie Re(M) en z+i*0^+
! M(1) -> M11
! M(2) -> M22
! M(3) -> M33
! M(4) -> M13
! det=M(1)*M(3)-M(4)*2
 REAL(QP), INTENT(IN) :: z,q,omplasma
 REAL(QP), INTENT(OUT) :: M(1:4),det
 
 REAL(QP) bmax,ommax,rhoz(1:4)
 REAL(QP) I1(1:4),Iinf(1:4),Icomp(1:4)

 REAL(QP) res(1:1,1:3),bornes(1:4)
 INTEGER  choix(1:3)

 ommax=max(2*z,2*omang(q)) !On intègre analytiquement le contre-terme de om=0 à ommax
 bmax=1.0e7_qp             !Borne sup de l'intégration (pour om>bmax on intègre intmat analytiquement)
 
!Densité spectrale en om=z, utilisée pour compenser la divergence du terme en 1/(z-om)
 rhoz=rho(z,q)
 

!Découpage de l'intégrale
 bornes=(/2.0_qp,   omang(q),ommax,  bmax/)
 choix   =         (/msql,   msql,   minf/)    !changements de variable sur l'intervalle [bornes(i),bornes(i+1)]
 res(1,:)=         (/-1.0_qp,-1.0_qp,-1.0_qp/) !res(i)=+1 => le dénominateur s'annule dans l'intervalle [bornes(i),bornes(i+1)].

!Compensation des segments résonants
 if((z<omang(q)).AND.(2.0_qp<z))then 
  res(1,1)=1.0_qp
  Icomp=-rhoz*log((omang(q)-z)/(z-2.0_qp))   !=int_2^omang     dom*rho(z)/(z-om), intégrale du terme qui régularise la résonance
 elseif(omang(q)<z)then 
  res(1,2)=1.0_qp
  Icomp=-rhoz*log((ommax-z)/(z-omang(q)))    !=int_omang^ommax dom*rho(z)/(z-om)
 else
  Icomp(:)=0.0_qp
 endif
 I1=decoupe(intmat,bornes,4,res(1:1,:),choix,EPS,blaBCS)


!Contribution des om>bmax
 Iinf(1)=  -z**2*q/(4*bmax**2)-intcterme(ommax)
 Iinf(2)=  -z**2*q/(4*bmax**2)-intcterme(ommax)
 Iinf(3)=  -q     /bmax**2
 Iinf(4)=  -z*q   /(2.0_qp*bmax**2)
 if(blaBCS) write(6,*) "Iinf=",Iinf
 if(blaBCS) write(6,*) "Icomp=",Icomp
 
! Finalisation
 M=I1+Icomp+Iinf
 M(3)=M(3)-2*q**2/(3*omplasma**2) !terme 1/g_33=1/2V(q) dans M33

 if(blaBCS) write(6,*)"q,z,M=",q,z,M
 if(blaBCS) write(6,*)

 det=M(1)*M(3)-M(4)**2
 
!Integrande
 CONTAINS
  FUNCTION intmat(om,arg,m)
   INTEGER,  INTENT(IN) :: m
   REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg
   REAL(QP), DIMENSION(size(om),m)       ::  intmat
   REAL(QP) rhom(1:4),ome,cterm(1:4),reson
   INTEGER is
 
   reson=arg(1)
   do is=1,size(om)
    ome=om(is)
    rhom=rho(ome,q)
 
  ! Contre-terme pour la convergence UV
    cterm(:)=0.0_qp
    if(ome>ommax) cterm(1:2)=-0.5_qp/sqrt(ome**2-4)
 
    if(reson<0.0_qp)then 
  ! Pas de résonance
      intmat(is,:)= rhom           *(1/(z-ome)-s/(z+ome))-cterm
!      write(6,*)"ome,intmat(is,1)*om**2=",ome,intmat(is,1)*ome**2
!      write(6,*)"ome,rhom(1)*(1/(z-ome)-s/(z+ome))*om,cterm(1)=",ome,rhom(1)*(1/(z-ome)-s/(z+ome))*ome,cterm(1)
!      write(6,*)"ome,rhom(1)*(1/(z-ome)-s/(z+ome))*om,cterm(1)=",ome,rhom(1)*(1/(z-ome)-s(1)/(z+ome))*ome,cterm(1)
!      write(6,*)"rhom(1)=",ome,rhom(1)
    else
  ! On compense la résonance en retranchant rhoz
      intmat(is,:)=(rhom-rhoz)/(z-ome)-s*rhom/(z+ome)
    endif
 
   enddo
  END FUNCTION intmat
 END SUBROUTINE mat_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE dmat_r(z,q,omplasma,dM)
! Calcule dM/domega, utilisé par Zero pour la méthode de Newton (plutôt qu'une différence finie)
 REAL(QP), INTENT(IN) :: z,q,omplasma
 REAL(QP), INTENT(OUT) :: dM(1:4)
 
 REAL(QP) bmax,ommax
 REAL(QP) I1(1:4),Iinf(1:4)
 REAL(QP) bornes(1:4),res(1:1,1:3)
 INTEGER  choix(1:3)

 if(z>2.0_qp) stop "Erreur: z>2 dans dmat_r"
 ommax=2*omang(q)
 bmax=1.0e7_qp
 
 if(blaBCS) write(6,*) "dmat_r"
 if(blaBCS) write(6,*) 

!Decoupage de l'intégrale
 bornes=(/2.0_qp,   omang(q),ommax,  bmax/)
 choix   =         (/msql,   msql,   minf/)    !changements de variable sur l'intervalle [bornes(i),bornes(i+1)]
 res(1,:)=         (/-1.0_qp,-1.0_qp,-1.0_qp/) !A priori jamais de résonance pour dM/dz (on se restreint à z<2).
!Integrale
 I1=decoupe(intdmat,bornes,4,res(1:1,:),choix,EPS,blaBCS)
 
!Contribution des om>bmax
 Iinf(1)=  -z**2*q  /(2*bmax**2)
 Iinf(2)=  -z**2*q  /(2*bmax**2)
 Iinf(3)=  -z*q     /bmax**4
 Iinf(4)=  -q       /(2*bmax**2)
 
!Finalisation
 dM=I1+Iinf
 if(blaBCS) write(6,*)"q,z,dM=",q,z,dM
 if(blaBCS) write(6,*)
 
!Integrande
 CONTAINS
  FUNCTION intdmat(om,arg,m)
   INTEGER,  INTENT(IN) :: m!m=size(I1)=4 ici
   REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg
   REAL(QP), DIMENSION(size(om),m)       ::  intdmat
   REAL(QP) rhom(1:4),ome
   INTEGER is
 
   do is=1,size(om)
    ome=om(is)
    rhom=rho(ome,q)
 
    intdmat(is,:)= -rhom           *(1/(z-ome)**2-s/(z+ome)**2)
 
   enddo
  END FUNCTION intdmat
 END SUBROUTINE dmat_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE ZeroBCSc(q,omplasma,sec,z,tol,residu,higgs)
! Cherche la racine complexe de M11*M33-M13^2=0 (ou de M22=0 si higgs=0)
! dans le prolongement analytique.
! sec=1 prolongement par la fenêtre II
! sec=2 prolongement par la fenêtre III
! Utilise la méthode de Newton globalement convergente à 2D (cf. [2] section 9.7) : z doit quand même être bien initialisé
! residu renvoi le résidu de M^{-1} en z
! tol(1)=tolérance sur x, tol(2)=tolérance sur det
 USE recettes, ONLY : mnewt,newt
 IMPLICIT NONE
 REAL(QP),     INTENT(IN) :: q,omplasma,tol(1:2)
 COMPLEX(QPC), INTENT(INOUT) :: z
 COMPLEX(QPC), INTENT(OUT) :: residu(1:4)
 INTEGER     , INTENT(IN) :: sec
 LOGICAL,      INTENT(IN) :: higgs

 COMPLEX(QPC) :: M(1:4),det,ddet
 REAL(QP), DIMENSION(1:2) :: zdep
 LOGICAL verif
 
 REAL(QP), DIMENSION(1:2) :: detvec
 REAL(QP), DIMENSION(1:2,1:2) :: Jdet
 
 if(blaZ)then
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)'ZeroBCSc'
  write(6,*)
  
  write(6,*)'--------------------'
  write(6,*)'q=',q
  write(6,*)'z=',real(z),imag(z)
  write(6,*)'omplasma=',omplasma
 endif
 
 zdep=(/real(z),imag(z)/) !version vectorielle de z

!Méthode de Newton globalement convergente
 call newt (zdep,tol(2),verif,detpro2)

 z=cmplx(zdep(1),zdep(2),kind=qpc)
!Valeur de det et M en la valeur finale de z
 call detpro1 (zdep,detvec,Jdet)

!Résidu
 residu(1)= M(3)/ddet
 residu(2)= 1.0_qp/M(2)
 residu(3)= M(1)/ddet
 residu(4)=-M(4)/ddet

 if(blaZ)then
  write(6,*)'----------------------'
  write(6,*)'zfinal=',z
  write(6,*)'re M=',real(M)
  write(6,*)'im M=',imag(M)
  write(6,*)'re  residu=',real(residu)
  write(6,*)'im  residu=',imag(residu)
  write(6,*)'abs residu=', abs(residu)
  write(6,*)'convergence?=',verif
 endif
 
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fonction complexe (ou plutôt 2D) dont newt doit trouver un racine
  FUNCTION detpro2(zz)
  IMPLICIT NONE
  REAL(QP), DIMENSION(:), INTENT(IN) :: zz
  REAL(QP), DIMENSION(size(zz)) :: detpro2

  COMPLEX(QPC) Mb(1:4)
  COMPLEX(QPC) z,dz,det,detP,detM
  REAL(QP) drez,dimz

  if(blaZ)then
   write(6,*)'----------------------'
   write(6,*)'zz=',zz
  endif
 

  z  =cmplx(zz(1)     ,zz(2)     ,kind=qpc)
  
  call matpro  (z   ,q,omplasma,sec,M ,det, higgs)
  if(blaZ) write(6,*)'z,det=',z,det
  detpro2=(/real(det),imag(det)/)
  END FUNCTION detpro2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE detpro1(zz,detvec,Jdet)
! Fonction 2D pour la méthode de Newton bête
! Utilisée par ZeroBCSc pour calculer le residu
! Partage M,det et ddet avec la sousroutine parente
  IMPLICIT NONE
  REAL(QP), DIMENSION(:), INTENT(IN) :: zz
  REAL(QP), DIMENSION(:), INTENT(OUT) :: detvec
  REAL(QP), DIMENSION(:,:), INTENT(OUT) :: Jdet
  
  COMPLEX(QPC) Mb(1:4)
  COMPLEX(QPC) z,dz,det,detP,detM
  REAL(QP) drez,dimz

  if(blaZ)then
   write(6,*)'----------------------'
   write(6,*)'zz=',zz
  endif
 

  drez=zz(1)*1.e-7_qp
  dimz=zz(2)*1.e-7_qp
  dz=drez+iiq*dimz
  
  z  =cmplx(zz(1)     ,zz(2)     ,QPC)
  
  call matpro  (z   ,q,omplasma,sec,M ,det  ,higgs)
  if(blaZ) write(6,*)'z,det=',z,det
  call matpro  (z+dz,q,omplasma,sec,Mb,detP ,higgs)
  if(blaZ) write(6,*)'z+dz,detP=',z+dz,detP
  call matpro  (z-dz,q,omplasma,sec,Mb,detM ,higgs)
  if(blaZ) write(6,*)'z-dz,detM=',z-dz,detM

  detvec(1)=real(det)
  detvec(2)=imag(det)
  
  ddet=(detP-detM)/(2.0_qp*dz)
  
  Jdet(1,1)=real(ddet)
  Jdet(2,1)=imag(ddet)
  Jdet(1,2)=-Jdet(2,1)
  Jdet(2,2)=Jdet(1,1)
  
  if(blaZ) write(6,*)'ddet=',real(ddet),imag(ddet)
  
  END SUBROUTINE detpro1
 END SUBROUTINE ZeroBCSc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE matpro(z,q,omplasma,sec,M,det,higgs)
! Combine mat_c et rhopro pour former la matrice M prolongée analytiquement
! sec=1 prolongement par la fenêtre II
! sec=2 prolongement par la fenêtre III
! renvoie det=M(1)*M(3)-M(4)**2 ou M(2) selon que higgs est faux ou vrai
 COMPLEX(QPC), INTENT(IN) :: z
 REAL(QP)    , INTENT(IN) :: q,omplasma
 INTEGER     , INTENT(IN) :: sec
 COMPLEX(QPC), INTENT(OUT) :: M(1:4),det
 LOGICAL,      INTENT(IN) :: higgs
 
 call mat_c(z,q,omplasma,M,det)
 if(imag(z)<0.0_qp) M=M-2*PI*iiq*rhopro(z,q,sec)
 det=M(1)*M(3)-M(4)**2
 if(higgs) det=M(2)
 if(blaBCS) write(6,*)
 if(blaBCS) write(6,*) "re(Mpro)=",real(M)
 if(blaBCS) write(6,*) "im(Mpro)=",imag(M)
 if(blaBCS) write(6,*) "detpro=",real(det),imag(det)
 if(blaBCS) write(6,*)
 END SUBROUTINE matpro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE mat_c(z,q,omplasma,M,det)
! Calcule la matrice M pour z complexe.
! Formellement identique à mat_r sauf que z,M et det deviennent complexes
 COMPLEX(QPC), INTENT(IN) :: z
 REAL(QP)    , INTENT(IN) :: q,omplasma
 COMPLEX(QPC), INTENT(OUT) :: M(1:4),det
 
 REAL(QP) bmax,ommax,rhoz(1:4),ommil,ommil2,ommil3
 COMPLEX(QPC) I1(1:4),I1b(1:4),I2(1:4),I2b(1:4),I2c(1:4),I3(1:4),Iinf(1:4),Icomp(1:4)

 REAL(QP) bornes(1:7),res(1:1,1:6)
 INTEGER  choix(1:6)

 ommax=max(2*real(z),2*omang(q))
 bmax=1.0e7_qp
 
!Densité spectrale en omega=real(z) pour compenser une quasidivergence
 rhoz=rho(real(z),q)
 
 I1   (:)=0.0_qp
 Icomp(:)=0.0_qp
 Iinf (:)=0.0_qp

!Dans le cas où z est complexe, a priori pas besoin de compenser l'intégrale (voir exception plus bas)
 res(1,:)=       (/-1.0_qp,-1.0_qp, -1.0_qp,-1.0_qp,-1.0_qp,-1.0_qp/) !res(i)=+1 => le dénominateur s'annule dans l'intervalle [bornes(i),bornes(i+1)].

!On subdivise les intervalles [2,omang] et [omang,ommax] pour isoler les singularités aux bords
!et pour gérer le cas où z s'approche de l'axe réel de sorte que l'intégrale est quasidivergente en real(z)
 if((2.0_qp<real(z)).AND.(real(z)<omang(q)))then
  ommil=real(z)
  if(abs(imag(z))<0.05_qp) res(1:1,1:2)=1.0_qp !Dans ce cas, om compense l'intégrale de 2 à omang
  if(abs(imag(z))<0.05_qp) Icomp=-rhoz*log((z-omang(q))/(z-2.0_qp)) 
 else
  ommil=(omang(q)+2.0_qp)/2.0_qp
 endif
!
 if((omang(q)<real(z)).AND.(real(z)<ommax))then
  ommil3=real(z)
 else
  ommil3=(omang(q)+ommax)/2.0_qp
 endif
 ommil2=(ommil3+omang(q))/2.0_qp

!Découpage de l'intégrale
 bornes=(/2.0_qp,  ommil,  omang(q),ommil2, ommil3, ommax,  bmax/)
 choix   =       (/msqu,   msql,    msql,   msqu,   msql,   rinf/)    !changements de variable sur l'intervalle [bornes(i),bornes(i+1)]
 I1=decoupe(intmatc,bornes,4,res(1:1,:),choix,EPS,blaBCS)

! Contributions des om>bmax
 Iinf(1)=  -z**2*q/(4*bmax**2)-intcterme(ommax)
 Iinf(2)=  -z**2*q/(4*bmax**2)-intcterme(ommax)
 Iinf(3)=  -q     /bmax**2
 Iinf(4)=  -z*q   /(2.0_qp*bmax**2)

! Finalisation
 if(blaBCS) write(6,*) "re(Iinf)=",real(Iinf)
 if(blaBCS) write(6,*) "im(Iinf)=",imag(Iinf)
 if(blaBCS) write(6,*) "re(Icomp)=",real(Icomp)
 if(blaBCS) write(6,*) "im(Icomp)=",imag(Icomp)
 
 M=I1+Icomp+Iinf
 M(3)=M(3)-2*q**2/(3*omplasma**2)

 if(blaBCS) write(6,*)"q,z=",q,real(z),imag(z)
 if(blaBCS) write(6,*)"re M=",real(M)
 if(blaBCS) write(6,*)"im M=",imag(M)
 if(blaBCS) write(6,*)

 det=M(1)*M(3)-M(4)**2
 
 CONTAINS
  FUNCTION intmatc(om,arg,m)
 ! energy integrand
  INTEGER,  INTENT(IN) :: m
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg
  COMPLEX(QPC), DIMENSION(size(om),m) ::  intmatc
  REAL(QP) rhom(1:4),ome,cterm(1:4),reson
  INTEGER is
 
  reson=arg(1)
  do is=1,size(om)
   ome=om(is)
   rhom=rho(ome,q)
   cterm(:)=0.0_qp
   if(ome>ommax) cterm(1:2)=-0.5_qp/sqrt(ome**2-4)
 
   if(reson<0.0_qp)then 
     intmatc(is,:)= rhom           *(1/(z-ome)-s/(z+ome))-cterm
   else !On compense l'intégrale
     intmatc(is,:)=(rhom-rhoz)/(z-ome)-s*rhom/(z+ome)
   endif
 
  enddo
  END FUNCTION intmatc
 END SUBROUTINE mat_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE dmat_c(z,q,omplasma,dM)
! dM/dz dans le cas complexe
! Peut servir pour ZeroBCSc mais n'est pas utilisée dans la version actuelle
 REAL(QP)    , INTENT(IN) :: q,omplasma
 COMPLEX(QPC), INTENT(IN) :: z
 COMPLEX(QPC), INTENT(OUT) :: dM(1:4)
 
 REAL(QP) bmax,ommax
 COMPLEX(QPC) I1(1:4),Iinf(1:4)

 REAL(QP) bornes(1:4),res(1:1,1:3)
 INTEGER  choix(1:3)

 ommax=2*omang(q)
 bmax=1.0e7_qp
 
!Decoupage de l'intégrale
 bornes=(/2.0_qp,   omang(q),ommax,  bmax/)
 choix   =         (/msql,   msql,   minf/)    !changements de variable sur l'intervalle [bornes(i),bornes(i+1)]
 res(1,:)=         (/-1.0_qp,-1.0_qp,-1.0_qp/) !res(i)=+1 => le dénominateur s'annule dans l'intervalle [bornes(i),bornes(i+1)].
 I1=decoupe(intdmatc,bornes,4,res(1:1,:),choix,EPS,blaBCS)

 Iinf(1)=  -z**2*q/(2*bmax**2)
 Iinf(2)=  -z**2*q/(2*bmax**2)
 Iinf(3)=  -z*q   /bmax**4
 Iinf(4)=  -q     /(2*bmax**2)
 
 dM=I1+Iinf
 if(blaBCS) write(6,*)"q,z,dM=",q,z,dM
 if(blaBCS) write(6,*)
 
 CONTAINS
  FUNCTION intdmatc(om,arg,m)
  INTEGER,      INTENT(IN) :: m
  REAL(QP),     INTENT(IN), DIMENSION(:)  ::  om,arg
  COMPLEX(QPC), DIMENSION(size(om),m)       ::  intdmatc
  REAL(QP)      rhom(1:4),ome
  INTEGER is
 
  do is=1,size(om)
   ome=om(is)
   rhom=rho(ome,q)
 
   intdmatc(is,:)= -rhom           *(1/(z-ome)**2-s/(z+ome)**2)
 
  enddo
  END FUNCTION intdmatc
 END SUBROUTINE dmat_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION rho1(om,q) 
! Densité spectrale rho_ij dans le secteur II (2<om<omang)
! Utilise les fonctions elle et ellf de recettes.f90 (convention Gradshteyn)
! Noter que la dépendence en q est triviale
 REAL(QP), INTENT(IN) :: om,q
 REAL(QP) rho1(1:4)
 
 REAL(QP) elliptE,elliptF
 
 elliptE=elle(PI/2.0_qp,mm(om))
 elliptF=ellf(PI/2.0_qp,mm(om))
 
 rho1(1)=om*elliptE/8.0_qp/q
 rho1(2)=rho1(1)-elliptF/(2.0_qp*om)/q
 rho1(3)=rho1(1)
 rho1(4)=elliptF/4.0_qp/q
 END FUNCTION rho1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION rho2(om,q) 
! Densité spectrale rho_ij dans le secteur III (om>omang)
 REAL(QP), INTENT(IN) :: om,q
 REAL(QP) rho2(1:4)
 
 REAL(QP) elliptE,elliptF
 
 elliptE=elle(theta(om,q),mm(om))
 elliptF=ellf(theta(om,q),mm(om))
 
 rho2(1)=om*elliptE/8.0_qp/q
 rho2(2)=rho2(1)-elliptF/(2.0_qp*om)/q
 rho2(3)=rho2(1)-sqrt((om**2-4*q**2-4)/(om**2-4*q**2))/4
 rho2(4)=elliptF/4.0_qp/q
 END FUNCTION rho2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION rho(om,q) 
! Densité spectrale rho_ij
! Combine rho1 et rho2
 REAL(QP), INTENT(IN) :: om,q
 REAL(QP) rho(1:4)
 
 if(om<2.0_qp)then
  rho(:)=0.0_qp
 elseif(om<omang(q))then
  rho=rho1(om,q)
 else
  rho=rho2(om,q)
 endif
 
 END FUNCTION rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION rhopro(z,q,sec) 
! Prolongement analytique de la fonction rho à z complexe
! Utilise la méthode des bornes d'intégration complexe
 REAL(QP),     INTENT(IN) :: q
 COMPLEX(QPC), INTENT(IN) :: z
 INTEGER, INTENT(IN) :: sec
 COMPLEX(QPC) rhopro(1:4)

 REAL(QP) :: tr1,tr2
 COMPLEX(QPC) :: u1,u2


 if(sec==1)then 

! On intègre pour u parcourant le segment [0,sqrt(z**2-4)/2q] du plan complexe
! On fait le changement de variable u=t*sqrt(z**2-4)/2q -> l'intégrale pour t dans [0,1] ne présente pas de difficulté particulière

  rhopro=qromovcq(intrhopro1,0.0_qp,1.0_qp,4,(/bidon/),midsquvcq,EPS)
  rhopro=rhopro*sqrt(z**2-4)/(2*q)

 elseif(sec==2)then
  if(.TRUE.)then
!  if((imag(z)<0.2_qp).AND.(real(z)<omang(q)).AND.(2.0_qp<real(z)))then

! L'intégrale est a priori sur [0,1] mais son intégrande devient complexe.w
! Elle peut présenter des difficultés quand imag(z) est petit et pas dans l'intervalle [omang,+infini[ (voir la condition "if" inutilisée plus haut)
! Notamment au voisinage de
! real(tr1) avec tr1=sqrt(z**2-4)/(2*q)
! real(tr2) avec tr2=z/(2*q)
! On contourne la difficulté en déformant [0,1] en un contour complexe 0 -> u1=-i/2 -> u2=1-i/2 -> 1

   u1=cmplx(0.0_qp,-0.5_qp,kind=qpc)
   u2=cmplx(1.0_qp,-0.5_qp,kind=qpc)
   rhopro=       qromovcq(intrhopro2,0.0_qp,1.0_qp   ,4,(/1.5_qp/),midsquvcq,EPS)*u1
   rhopro=rhopro+qromovcq(intrhopro2,0.0_qp,1.0_qp   ,4,(/2.5_qp/),midsquvcq,EPS)*(u2-u1)
   rhopro=rhopro+qromovcq(intrhopro2,0.0_qp,1.0_qp   ,4,(/3.5_qp/),midsquvcq,EPS)*(1.0_qp-u2)

  else
! Au cas où on voudrait intégrer directement sur [0,1]
   rhopro=       qromovcq(intrhopro2,0.0_qp,1.0_qp,4,(/0.5_qp/),midsquvcq,EPS)
  endif

 else
  stop "Mauvais secteur dans rhopro"
 endif
 if(blaBCS) write(6,*)"re rhopro=",real(rhopro)
 if(blaBCS) write(6,*)"im rhopro=",imag(rhopro)

 CONTAINS
! Intégrande secteur II
  FUNCTION intrhopro1(tt,arg,m)
  INTEGER,      INTENT(IN) :: m
  REAL(QP),     INTENT(IN), DIMENSION(:)  ::  tt,arg
  COMPLEX(QPC), DIMENSION(size(tt),m)       ::  intrhopro1

  COMPLEX(QPC)  r1(size(tt)),r2(size(tt))
  INTEGER is
 
!  do is=1,size(tt)
!   t=tt(is)
 
! Voir eqs. 134 et 139-142 dans BCSandRPAwithCoulomb.tex
   r1=sqrt(z**2-tt**2*(z**2-4))
   r2=sqrt((z**2-4)*(1-tt**2))

   intrhopro1(:,1)= r1/(4*r2)     
   intrhopro1(:,2)= r2/(4*r1)     
   intrhopro1(:,3)= z**2/r2/r1**3 
   intrhopro1(:,4)= z/(2*r1*r2)   
 
!  enddo
  END FUNCTION intrhopro1
! Intégrande secteur III
  FUNCTION intrhopro2(tt,arg,m)
  INTEGER,      INTENT(IN) :: m
  REAL(QP),     INTENT(IN), DIMENSION(:)  ::  tt,arg
  COMPLEX(QPC), DIMENSION(size(tt),m)       ::  intrhopro2

  COMPLEX(QPC), DIMENSION(size(tt)) :: u,r1,r2
!  COMPLEX(QPC), DIMENSION(m)       :: cterme(1:4)
  INTEGER is
 
!  do is=1,size(tt)
!   t=tt(is)

! Changements de variable dans les différents segments du contour complexe
   if(floor(arg(1))==0)then
    u=tt*cmplx(1.0_qp,0.0_qp,kind=qpc)
   elseif(floor(arg(1))==1)then
    u=tt*u1
   elseif(floor(arg(1))==2)then
    u=u1+(u2-u1)*tt
   elseif(floor(arg(1))==3)then
    u=u2+(1.0_qp-u2)*tt
   endif
 
   r1=sqrt(z**2-4*q**2*u**2)
   r2=sqrt(z**2-4-4*q**2*u**2)

   intrhopro2(:,1)= r1/(4*r2)
   intrhopro2(:,2)= r2/(4*r1)
   intrhopro2(:,3)= z**2/r2/r1**3
   intrhopro2(:,4)= z/(2*r1*r2)
 
!  enddo
  END FUNCTION intrhopro2
 END FUNCTION rhopro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION mm(om)      !argument module des intégrales elliptiques. (Attention convention Gradshteyn mm_Mathematica=mm**2)
 REAL(QP), INTENT(IN) :: om
 REAL(QP) mm
 mm=sqrt(om**2-4)/om
END FUNCTION mm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION theta(om,q) !argument angulaire des intégrales elliptiques quand om>omang
 REAL(QP), INTENT(IN) :: om,q
 REAL(QP) theta
 theta=asin(2*q/sqrt(om**2-4))
 END FUNCTION theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION omang(q) !Point anguleux om2 (le seul dans la limite BCS)
 REAL(QP), INTENT(IN) :: q
 REAL(QP) omang
 omang=2*sqrt(1+q**2)
 END FUNCTION omang
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION intcterme(bsup) !=int_2^{bsup} (-0.5)/sqrt(om**2-4)
 REAL(QP), INTENT(IN) :: bsup
 REAL(QP) intcterme
 intcterme= -0.5_qp*acotanh(bsup/sqrt(-4+bsup**2))
 END FUNCTION intcterme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION acotanh(x)
 REAL(QP), INTENT(IN) :: x
 REAL(QP) acotanh
 acotanh= log((x+1.0_qp)/(x-1.0_qp))/2.0_qp
 END FUNCTION acotanh
END MODULE dspecBCS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
