

Program pop
IMPLICIT NONE

REAL(kind=8),ALLOCATABLE :: eff_Q(:,:)    ! QTL effect of 1001 generation
REAL(kind=8),ALLOCATABLE :: pheno(:,:)      ! phenotype of the training set
REAL(kind=8),ALLOCATABLE :: eff_q_l(:,:)    ! QTL effect from gamma distribution
REAL(kind=8),ALLOCATABLE :: adde(:,:)       ! add effect 
REAL(kind=8),ALLOCATABLE :: coe(:,:,:)      !coe--coefficient of the qtl genotype used to calculating the add effect
REAL(kind=8),ALLOCATABLE :: adde1002(:,:),pheno1002(:,:)    ! breeding value of the 1002 generation
REAL(kind=8),ALLOCATABLE :: ebvm(:,:),ebvs(:,:),ebvd(:,:)  
REAL(kind=8),ALLOCATABLE :: deltaF(:,:), gg(:,:), cor(:,:), rate_F(:,:), mean_rate_F(:),mean_gg(:),sd_gg(:), ad(:,:), sd_rate_F(:)
REAL(kind=8),ALLOCATABLE :: r_h2(:,:), Var_G(:,:),Var_P(:,:), Var_E(:,:), m_rh2(:),m_VG(:),m_VP(:),m_VE(:),sd_VG(:),sd_rh2(:), m_F(:), sd_F(:),sd_VE(:)

INTEGER,ALLOCATABLE :: m(:,:,:)    ! Ne of 1~1000 
INTEGER,ALLOCATABLE :: m_qs(:,:,:)  ! sire:geno and effect of QTL and SNP 1~1000 generation 
INTEGER,ALLOCATABLE :: m_qsswap(:,:,:)    ! swap variable 
INTEGER,ALLOCATABLE :: m_qd(:,:,:)        ! dam :geno and effect of QTL and SNP 1~1000 generation 
INTEGER,ALLOCATABLE :: m_qdswap(:,:,:)    ! swap variable 
INTEGER,ALLOCATABLE :: m1001(:,:,:)
INTEGER,ALLOCATABLE :: m1001_c(:,:,:)     ! combine the 1,2 type into 0, 1, 2
INTEGER,ALLOCATABLE :: qtl(:,:,:)         ! QTL genotype of the following generation 
INTEGER,ALLOCATABLE :: s1001(:,:,:), d1001(:,:,:)     ! to product 1001
INTEGER,ALLOCATABLE :: s1001swap(:,:,:), d1001swap(:,:,:)     
INTEGER,ALLOCATABLE :: T1(:,:,:),n_eQTL(:,:),T2(:,:,:)
INTEGER,ALLOCATABLE :: s1002(:,:,:), d1002(:,:,:),m1002(:,:,:),m1002_c(:,:,:)
INTEGER,ALLOCATABLE :: s1003(:,:,:), d1003(:,:,:),m1003(:,:,:),m1003_c(:,:,:), h(:,:) 

REAL(kind=8),ALLOCATABLE :: s(:,:),mean(:,:)
INTEGER :: Ne, nGen, nQTL, n1001, nChr, nChr1, nr, nShuf,nShuf1, nu
INTEGER :: a, b, c, d, T(1,1,1), ab(1)                            !using to shuffle the parents
INTEGER :: ind, ind1
REAL(kind=8) :: cr, h2                        !combination ratio
INTEGER :: ir, iGen, iChr, iID, iQTL, iOff, i, j, k, iu
REAL(kind=8) :: VG, VGi, VE, lambda, mean_sd
REAL(kind=8) :: mr_q !mutation rate of QTL and SNP
REAL(kind=8) :: u, v, w, z, w1, wn(1000,1), wg
REAL(kind=8) :: r_eQTL
integer::x, xx
real(kind=8) :: e
integer:: sirem(20000,1),damm(20000,1),unitn, m_neQTL, m_h
real(kind=8):: Fm(20000,1),sigm(20000,1)
real(kind=8):: ps,pd
real(kind=8):: mean_cor, sd_cor, gamma_a, gamma_b, noa_m 

gamma_a = 0.4
gamma_b = 1.66
noa_m = 0.0
nr = 50                            !nmuber of reiteration
Ne = 100                              !Ne of 1~1000 generation
nChr = 10                             !number of Chromosome
n1001 = 1000                          !number of training set
nGen = 1000                           !generation of the ideal population
nQTL = 100                            !number of QTL on each chromosome
VG =  1.0                             !environmental variance 
h2 = 0.3
VE = (1.0-h2)*VG/h2                              !genetic variance 
VGi = VG/(nChr*nQTL)                  !mean locus variance 
lambda = VE/VGi 
mr_q = 2.5E-5                         !mutation rate of QTL
nShuf1 = 200
nShuf = 2000
nu = 20
mean_sd = 0.0
sirem = 0
damm = 0
Fm = 0.0
sigm = 0.0
ps = 0.05
pd = 0.2
cr = 0.5*(1-exp(-2/(REAL(nQTL-1))))
!CALL hmf(1.0/REAL((nm+1)*nQTL+nm+1),cr)      !calculating the combination rate using Haldane's mapping function
CALL random_seed()
ALLOCATE (m(Ne,2*nChr,nQTL),m_qs(Ne/2,2*nChr,nQTL),m_qd(Ne/2,2*nChr,nQTL),m_qsswap(1,2*nChr,nQTL),m_qdswap(1,2*nChr,nQTL),eff_q_l(nChr,nQTL),&
         m1001(n1001,2*nChr,nQTL),m1001_c(n1001,nChr,nQTL), eff_Q(nChr,nQTL),qtl(n1001,nChr,nQTL), T1(1,2*nChr,nQTL),&
         s1001(n1001/2,2*nChr,nQTL),d1001(n1001/2,2*nChr,nQTL), adde(n1001,1),coe(n1001,nChr,nQTL), pheno(n1001,1),&
         s(nr,1),mean(nr,1), n_eQTL(nr,1), m1002(n1001,2*nChr,nQTL),s1002(int(ps*n1001/2),2*nChr,nQTL),d1002(int(pd*n1001/2),2*nChr,nQTL),&
         adde1002(n1001,1),pheno1002(n1001,1),m1002_c(n1001,nChr,nQTL),r_h2(nr,nu),Var_G(nr,nu),Var_P(nr,nu),Var_E(nr,nu),m_rh2(nu),&
         m_VG(nu),m_VP(nu),m_VE(nu),sd_rh2(nu),sd_VG(nu),sd_VE(nu))

ALLOCATE(ebvm(n1001,2),ebvs(int(ps*n1001/2),2),ebvd(int(pd*n1001/2),2), deltaF(nr,nu-1), gg(nr,nu-1), cor(nr,1),rate_F(nr,nu-1),mean_rate_F(nu-1),&
         mean_gg(nu-1),sd_gg(nu-1),ad(nu-1,1),sd_rate_F(nu-1), m_F(nu-1), sd_F(nu-1), h(nr,n1001))
s = 0.0
mean = 0.0
DO ir = 1, nr    !loop of reiteration
eff_Q(:,:) = 0.0
 qtl(:,:,:) = 0
 coe(:,:,:) = 0.0
m = 1     !genotype of qtl and markers are initialized to be 0


DO iGen = 2, nGen    !loop of generation


!***********************************************************mutation occured here**************************


    DO iID = 1, Ne
        DO iChr = 1, 2*nChr
            DO i = 1, nQTL
                
                CALL random_number(u)
                IF(u<= mr_q) m(iID,iChr,i) = 3 - m(iID,iChr,i)    !mutation of the QTLs


            END DO    !i
        END DO    !iChr       
    END DO    !iID


!*********************************************************************shuffle the sires and dams to ensure the random mating 


DO i = 1, nShuf1 
!************************************ the sires
    CALL random_number(u) 
    a = int((Ne*u)/2)+1 
    CALL random_number(u) 
    b = int((Ne*u)/2)+1 

    T1(1,:,:) = m(a,:,:)
    m(a,:,:) = m(b,:,:)    
    m(b,:,:) = T1(1,:,:)
END DO    !i

!************************************ the dams
    DO i = 1, nShuf1
    CALL random_number(u) 
    a = int((Ne*u)/2)+1+Ne/2 
    CALL random_number(u) 
    b = int((Ne*u)/2)+1+Ne/2 
    
    T1(1,:,:) = m(a,:,:)      
    m(a,:,:) = m(b,:,:)    
    m(b,:,:) = T1(1,:,:)


END DO    !i

!*******************************************************adjust the index of the parents**************

!the first half are sires and the other are dams
DO iID = 1, Ne/2
    m_qs(iID,:,:) = m(iID,:,:)
    m_qd(iID,:,:) = m(iID+Ne/2,:,:)
END DO    !iID



!*********************************************************** mating to product  offsprings*********************

DO iID = 1, Ne/2        
   CALL random_number(u)
     a = int((Ne*u)/2) + 1 
   CALL random_number(u)
     b = int((Ne*u)/2) + 1
  DO iOff = 1, 2   !offsprings
    
     m_qsswap(1,:,:) = m_qs(a,:,:)
     m_qdswap(1,:,:) = m_qd(b,:,:)

     IF(iOff == 1) THEN
         ind = iID
     ELSE
         ind = iID+Ne/2
     END IF
     DO iChr = 1, nChr    !Chromosome
!*************************************************************
         
         DO iQTL = 1, nQTL-1   
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----sire
             DO i = iQTL+1, nQTL
                T(1,1,1) = m_qsswap(1,iChr*2-1,i)   
                m_qsswap(1,iChr*2-1,i) = m_qsswap(1,iChr*2,i)
                m_qsswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m(ind,2*iChr-1,:) = m_qsswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---sire
    ELSE
        m(ind,2*iChr-1,:) = m_qsswap(1,2*iChr,:)
    END IF
!*************************************************************

         DO iQTL = 1, nQTL-1
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----dam
             DO i = iQTL+1, nQTL
                T(1,1,1) = m_qdswap(1,iChr*2-1,i)   
                m_qdswap(1,iChr*2-1,i) = m_qdswap(1,iChr*2,i)
                m_qdswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m(ind,2*iChr,:) = m_qdswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---dam
    ELSE
        m(ind,2*iChr,:) = m_qdswap(1,2*iChr,:)
    END IF

    
    END DO    !iChr

  END DO    !iOff 
END DO    !iID


END DO    !iGen


!#############################################################################################################################



!1000-->1001



!############################################################################################################################## 

!***************************************************************shuffle the sires and dams to ensure the random mating 


DO i = 1, nShuf1 
!************************************ the sires
    CALL random_number(u)     
    a = int((Ne*u)/2)+1 
    CALL random_number(u) 
    b = int((Ne*u)/2)+1 
    
    T1(1,:,:) = m(a,:,:)
    m(a,:,:) = m(b,:,:)    
    m(b,:,:) = T1(1,:,:)
END DO    !i
!************************************ the dams
DO i = 1, nShuf1
	CALL random_number(u) 
    a = int((Ne*u)/2)+1+Ne/2
    CALL random_number(u) 
    b = int((Ne*u)/2)+1+Ne/2 
    
    T1(1,:,:) = m(a,:,:)
    m(a,:,:) = m(b,:,:)    
    m(b,:,:) = T1(1,:,:)


END DO    !i

!the first half are sires and the other are dams
DO iID = 1, Ne/2
    m_qs(iID,:,:) = m(iID,:,:)
    m_qd(iID,:,:) = m(iID+Ne/2,:,:)
END DO    !iID




!*********************************************************** mating to product  offsprings*********************

DO iID = 1, Ne/2        
    CALL random_number(u)
     a = int((Ne*u)/2) + 1 
   CALL random_number(u)
     b = int((Ne*u)/2) + 1

  DO iOff = 1, 20    !offspring
     m_qsswap(1,:,:) = m_qs(a,:,:)   
     m_qdswap(1,:,:) = m_qd(b,:,:)
     IF(iOff <= 10 ) THEN
         ind = 10*(iID-1)+iOff
     ELSE
         ind = n1001/2+10*(iID-2)+iOff
     END IF
     DO iChr = 1, nChr    !Chromosome

     !write(99,*) m_qsswap(1,:,:)
!*************************************************************
         
         DO iQTL = 1, nQTL-1
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----sire
             DO i = iQTL+1, nQTL
                T(1,1,1) = m_qsswap(1,iChr*2-1,i)   
                m_qsswap(1,iChr*2-1,i) = m_qsswap(1,iChr*2,i)
                m_qsswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m1001(ind,2*iChr-1,:) = m_qsswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---sire
    ELSE
        m1001(ind,2*iChr-1,:) = m_qsswap(1,2*iChr,:)
    END IF
!*************************************************************

         DO iQTL = 1, nQTL-1
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----dam
             DO i = iQTL+1, nQTL
                T(1,1,1) = m_qdswap(1,iChr*2-1,i)   
                m_qdswap(1,iChr*2-1,i) = m_qdswap(1,iChr*2,i)
                m_qdswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m1001(ind,2*iChr,:) = m_qdswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---dam
    ELSE
        m1001(ind,2*iChr,:) = m_qdswap(1,2*iChr,:)
    END IF

    
    END DO    !iChr

  END DO    !iOff 
END DO    !iID




!****************************************************************** finding the effective QTLs and SNPs 

   !open(unit=100,file='effq.txt')
   eff_Q = 0.0
   coe = 0.0 
    DO iChr = 1, nChr
            DO i = 1, nQTL
                x = sum(m1001(:,2*iChr-1,i))+sum(m1001(:,2*iChr,i))
                IF  (x > 2000)   THEN  ! qtl mutation occurred     
                CALL gamma(gamma_a,gamma_b,wg)                   !finding QTL          
                CALL random_number(u)
                    IF (u < 0.5) THEN
                       eff_Q(iChr,i) = wg
                    ELSE
                       eff_Q(iChr,i) = -1*wg
                    END IF
                ELSE
                    eff_Q(iChr,i) = 0.0
                END IF
                    
                 !   write(100,*) 'iChr=',iChr,'i=',i,eff_Q(iChr,i)
                 !   write(100,*)
            END DO     !i
        END DO    !iChr

!********************************************************calculating the QTL genotype, add effect and phenotype of each id

n_eQTL(ir,1) = 0
  DO iChr = 1, nChr
     DO i = 1, nQTL
         IF(eff_Q(iChr,i)/=0.0)  n_eQTL(ir,1) = n_eQTL(ir,1) +1
     END DO    !i
  END DO    !iChr 

! open(unit=111,file='adde.txt')
 
 adde(:,:) = 0.0 
   DO iID = 1, n1001
   x = 0
       DO iChr = 1, nChr
           DO i = 1, nQTL
               xx = m1001(iID,2*iChr-1,i) + m1001(iID,2*iChr,i)    !genotype of QTL
               select case(xx)
               case(2)
                   m1001_c(iID,iChr,i) = 0
                   coe(iID,iChr,i) = -1*eff_Q(iChr,i) 
               case(3)
                   m1001_c(iID,iChr,i) = 1
                   coe(iID,iChr,i) = 0.0
                   x = x + 1
               case(4)
                   m1001_c(iID,iChr,i) = 2
                   coe(iID,iChr,i) = eff_Q(iChr,i)
               end select
               !write(222,*) 'iID',iID,'iChr',iChr, 'i', i/2, 'qtl', qtl(iID,iChr,i/2), coe(iID,iChr,i/2)
               !write(222,*) 
           END DO    !i
       END DO    !iChr
 h(ir,iID) = x
   END DO    !iID
   
   !*************************** phenotype

   !open(unit=200,file='pheno1001.txt')
   DO iID = 1, n1001
       adde(iID,1) = sum(coe(iID,:,:))    !the add effect
   END DO    !iID
       CALL NOA(n1001,noa_m,VE,wn)       
       pheno =adde  + wn

      !write(200,*) 'iID=',iID,pheno(iID,1)
      !write(200,*) 

   w = 0.0 
   z = 0.0 
   w1 = 0.0 
   DO iID = 1, n1001
       w = w + (adde(iID,1)-sum(adde(:,1))/(REAL(n1001)))**2
       z = z + (pheno(iID,1)-sum(pheno(:,1))/(REAL(n1001)))**2 
       w1 = w1 + (wn(iID,1)-sum(wn(:,1))/(REAL(n1001)))**2
   END DO    !iID

   Var_G(ir,1) = w/(REAL(n1001-1))
   Var_P(ir,1) = z/(REAL(n1001-1))
   Var_E(ir,1) = w1/(REAL(n1001-1))
   r_h2(ir,1) = Var_G(ir,1)/Var_P(ir,1)


!***************************************************************shuffle the sires and dams to ensure the random mating to 1002 

DO i = 1, nShuf 
!************************************ the sires
    CALL random_number(u)     
    a = int(n1001*u/2)+1 
    CALL random_number(u) 
    b = int(n1001*u/2)+1 

    T1(1,:,:) = m1001(a,:,:)
    m1001(a,:,:) = m1001(b,:,:)    
    m1001(b,:,:) = T1(1,:,:)
END DO    !i
!************************************ the dams
DO i = 1, nShuf    
	CALL random_number(u) 
    a = int(n1001*u/2)+1+n1001/2
    CALL random_number(u) 
    b = int(n1001*u/2)+1+n1001/2 
    
    T1(1,:,:) = m1001(a,:,:)
    m1001(a,:,:) = m1001(b,:,:)    
    m1001(b,:,:) = T1(1,:,:)


END DO    !i
!the first half are sires and the other are dams
DO iID = 1, n1001/2
    s1001(iID,:,:) = m1001(iID,:,:)
    d1001(iID,:,:) = m1001(iID+n1001/2,:,:)
END DO    !iID



!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!caculating the pedigree of generation 1001
DO iID = 1, n1001
    sirem(iID,1) = 0
    damm(iID,1) = 0 
    Fm(iID,1) = 0.0
    sigm(iID,1) = 0.0
END DO    !iID

!DO i = 1, n1001
 !  Fm(i+1000,1) = 0.0
 !  sigm(i+1000,1) = real(15)/real(100)
!END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!writing the pedigree
open (unit=2,file='pedigree.txt')
DO i = 1, n1001
    write(2,'(3I7,2F10.6)') i, sirem(i,1), damm(i,1), Fm(i,1), sigm(i,1)
END DO    !i
!##########################################################################################################################



!1001-->1002






!#########################################################################################################################


!*********************************************************** mating to product  offsprings********

DO iID = 1, n1001/2        
100    CALL random_number(u)

     a = int(n1001/2*u) + 1 
   CALL random_number(u)
     b = int(n1001/2*u) + 1
 
     !IF (sirem(a,1)==sirem(b,1) .and. damm(a,1)==damm(b,1) ) goto 100
  DO iOff = 1, 2    !offspring
     m_qsswap(1,:,:) = s1001(a,:,:)   
     m_qdswap(1,:,:) = d1001(b,:,:)
     IF(iOff == 1) THEN
         ind = iID
     ELSE
         ind = iID+n1001/2
     END IF
     DO iChr = 1, nChr    !Chromosome

     !write(99,*) m_qsswap(1,:,:)
!*************************************************************
         
         DO iQTL = 1, nQTL-1
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----sire
             DO i = iQTL+1,nQTL
                T(1,1,1) = m_qsswap(1,iChr*2-1,i)   
                m_qsswap(1,iChr*2-1,i) = m_qsswap(1,iChr*2,i)
                m_qsswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m1002(ind,2*iChr-1,:) = m_qsswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---sire
    ELSE
        m1002(ind,2*iChr-1,:) = m_qsswap(1,2*iChr,:)
    END IF
    
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !pedigree
    sirem(ind+1000,1) = a
    damm(ind+1000,1) = b+n1001/2

!*************************************************************

         DO iQTL = 1, nQTL-1
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----dam
             DO i = iQTL+1, nQTL
                T(1,1,1) = m_qdswap(1,iChr*2-1,i)   
                m_qdswap(1,iChr*2-1,i) = m_qdswap(1,iChr*2,i)
                m_qdswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m1002(ind,2*iChr,:) = m_qdswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---dam
    ELSE
        m1002(ind,2*iChr,:) = m_qdswap(1,2*iChr,:)
    END IF

    
    END DO    !iChr

  END DO    !iOff 
END DO    !iID


! print*, 'generatioin', 1001


!******************************************************************calculating the add effect and phenotypic value of 1002

coe = 0.0
adde1002(:,:) = 0.0 
   DO iID = 1, n1001
       DO iChr = 1, nChr
           DO i = 1, nQTL
               xx = m1002(iID,2*iChr-1,i) + m1002(iID,2*iChr,i)    !genotype of QTL
               select case(xx)
               case(2)
                   m1002_c(iID,iChr,i) = 0
                   coe(iID,iChr,i) = -1*eff_Q(iChr,i) 
               case(3)
                   m1002_c(iID,iChr,i) = 1
                   coe(iID,iChr,i) = 0.0
               case(4)
                   m1002_c(iID,iChr,i) = 2
                   coe(iID,iChr,i) = eff_Q(iChr,i)
               end select
               !write(222,*) 'iID',iID,'iChr',iChr, 'i', i/2, 'qtl', qtl(iID,iChr,i/2), coe(iID,iChr,i/2)
               !write(222,*) 
           END DO    !i
       END DO    !iChr
       !write(111,*) 'iID=', iID, adde(iID,1)
   END DO    !iID
   !*************************** phenotype

   !open(unit=200,file='pheno1001.txt')
   DO iID = 1, n1001
       adde1002(iID,1) = sum(coe(iID,:,:))    !the add effect
   END DO    !iID
       CALL NOA(n1001,noa_m,VE,wn)       
       pheno1002 =adde1002  + wn


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!calculating the r_h2 
w = 0.0
z = 0.0
w1 = 0.0
DO iID = 1, n1001
   w = w + (adde1002(iID,1)-sum(adde1002(:,1))/(real(n1001)))**2
   z = z + (pheno1002(iID,1)-sum(pheno1002(:,1))/(real(n1001)))**2
   w1 = w1 + (wn(iID,1)-sum(wn(:,1))/(REAL(n1001)))**2   
END DO    !iID
   Var_G(ir,2) = w/(REAL(n1001-1))
   Var_P(ir,2) = z/(REAL(n1001-1))
   Var_E(ir,2) = w1/(REAL(n1001-1))
   r_h2(ir,2) = Var_G(ir,2)/Var_P(ir,2)



!^^^^^^^^^^^^^^^^^^^^^^^^
!writing the pedigree
DO i = 1, n1001
    write(2,'(3I7,F10.6,F10.6)') i+1000, sirem(i+1000,1), damm(i+1000,1)!, Fm(i+1000,1), sigm(i+1000,1)
END DO    !i

   close(2)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   !calculating inb_co
   call inbreed(1001,2000,deltaF(ir,1))



  

!@@@@@@@@@@@@@@@@@@@@@@@@@@ genetic gain
gg(ir,1) = (sum(adde1002(:,1))-sum(adde(:,1)))/(real(n1001))
ad(1,1) = sum(adde1002(:,1))/(real(n1001))


!##########################################################################################################################


!1003-->1020

!#########################################################################################################################

DO iu = 1, nu-2

DO i = 1, n1001
ebvm(i,1) = i
END DO    !i
ebvm(:,2) = pheno1002(:,1) 
call sort(ebvm,n1001)

!DO j = 1, n1001
 !   write(100,*) ebvm(j,1),ebvm(j,2) 
!END DO

!obtain the ebvs and ebvd
a = 0
b = 0
DO i = 1, n1001
    IF(int(ebvm(i,1))<=n1001/2) THEN
        a = a + 1
        IF (a>int(ps*n1001/2)) exit
		 ebvs(a,:) = ebvm(i,:)
    END IF
 END DO   !i
 DO i = 1, n1001
     IF(int(ebvm(i,1))>n1001/2) THEN
        b = b + 1
        IF(b>int(pd*n1001/2)) exit
		 ebvd(b,:) = ebvm(i,:)
     END IF
END DO    !i

!write(10,*) 'sires'

!DO i = 1, int(ps*n1001/2)
 !   write(100,*) ebvs(i,1),ebvs(i,2) 
!END DO

!write(10,*) 'dams'
!DO i = 1, int(ps*n1001/2)
!    write(100,*) ebvd(i,1),ebvd(i,2) 
!END DO

!obtain the s1002 and d1002
DO i = 1, int(ps*n1001/2)
    s1002(i,:,:) = m1002(int(ebvs(i,1)),:,:)
END DO

DO i = 1, int(pd*n1001/2)
    d1002(i,:,:) = m1002(int(ebvd(i,1)),:,:)
END DO
!*********************************************************** mating to product  offsprings********

!print *, 444444
DO iID = 1, n1001/2        
200  CALL random_number(u)

     a = int((ps*n1001/2-1)*u)+1 
     CALL random_number(u)
     b = int((pd*n1001/2-1)*u)+1 
 
     IF(sirem(a+iu*1000,1)==sirem(b+iu*1000,1) .and. damm(a+iu*1000,1)==damm(b+iu*1000,1))  goto 200
     DO iOff = 1, 2
     m_qsswap(1,:,:) = s1002(a,:,:)   
     m_qdswap(1,:,:) = d1002(b,:,:)

     IF(iOff == 1) THEN
         ind = iID
     ELSE
         ind = iID+n1001/2
     END IF 

     DO iChr = 1, nChr    !Chromosome

!*************************************************************
         
         DO iQTL = 1, nQTL-1
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----sire
             DO i = iQTL+1, nQTL
                T(1,1,1) = m_qsswap(1,iChr*2-1,i)   
                m_qsswap(1,iChr*2-1,i) = m_qsswap(1,iChr*2,i)
                m_qsswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m1002(ind,2*iChr-1,:) = m_qsswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---sire
    ELSE
        m1002(ind,2*iChr-1,:) = m_qsswap(1,2*iChr,:)
    END IF
    
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !pedigree
    sirem(ind+(iu+1)*1000,1) = int(ebvs(a,1))+iu*1000
    damm(ind+(iu+1)*1000,1) = int(ebvd(b,1))+iu*1000
    !*************************************************************

         DO iQTL = 1, nQTL-1
             CALL random_number(u)
             IF (u<=cr) THEN         !deciding whether the chromosome exchange----dam
             DO i = iQTL+1, nQTL
                T(1,1,1) = m_qdswap(1,iChr*2-1,i)   
                m_qdswap(1,iChr*2-1,i) = m_qdswap(1,iChr*2,i)
                m_qdswap(1,iChr*2,i) = T(1,1,1)
             END DO    !i
             END IF
         END DO    !iQTL


!*************************************************************

    CALL random_number(u)
    IF(u < 0.5) THEN
        m1002(ind,2*iChr,:) = m_qdswap(1,2*iChr-1,:)    !deciding which chromosome the offspring will inherit---dam
    ELSE
        m1002(ind,2*iChr,:) = m_qdswap(1,2*iChr,:)
    END IF

    END DO    !iChr
    END DO    !iOff
 END DO    !iID


!print *, 99999
!******************************************************************calculating the add effect and phenotypic value of 1002

coe = 0.0 
adde1002(:,:) = 0.0 
   DO iID = 1, n1001
       DO iChr = 1, nChr
           DO i = 1, nQTL
               xx = m1002(iID,2*iChr-1,i) + m1002(iID,2*iChr,i)    !genotype of QTL
               select case(xx)
               case(2)
                   m1002_c(iID,iChr,i) = 0
                   coe(iID,iChr,i) = -1*eff_Q(iChr,i) 
               case(3)
                   m1002_c(iID,iChr,i) = 1
                   coe(iID,iChr,i) = 0.0
               case(4)
                   m1002_c(iID,iChr,i) = 2
                   coe(iID,iChr,i) = eff_Q(iChr,i)
               end select 

           END DO    !i
       END DO    !iChr
       !write(111,*) 'iID=', iID, adde(iID,1)
   END DO    !iID
   !*************************** phenotype

   DO iID = 1, n1001
       adde1002(iID,1) = sum(coe(iID,:,:))    !the add effect
   END DO    !iID
       CALL NOA(n1001,noa_m,VE,wn)       
       pheno1002 =adde1002  + wn

  !print *, 1212121212 




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!calculating the r_h2 
 w = 0.0 
   z = 0.0 
   w1 = 0.0 
   DO iID = 1, n1001
       w = w + (adde1002(iID,1)-sum(adde1002(:,1))/(REAL(n1001)))**2 
       z = z + (pheno1002(iID,1)-sum(pheno1002(:,1))/(REAL(n1001)))**2 
       w1 = w1 + (wn(iID,1)-sum(wn(:,1))/(REAL(n1001)))**2
   END DO    !iID

  Var_G(ir,iu+2) = w/(REAL(n1001-1))
   Var_P(ir,iu+2) = z/(REAL(n1001-1))
   Var_E(ir,iu+2) = w1/(REAL(n1001-1))
   r_h2(ir,iu+2) = Var_G(ir,iu+2)/Var_P(ir,iu+2)
 



!^^^^^^^^^^^^^^^^^^^^^^^^
!writing the pedigree
open(unit=2,file='pedigree.txt',position='append')
DO i = 1, n1001
    write(2,'(3I7,2F10.6)') i+(iu+1)*1000, sirem(i+(iu+1)*1000,1), damm(i+(iu+1)*1000,1) 
END DO    !i

   close(2)
   !print *, 23232323
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   !calculating inb_co
   call inbreed((iu+1)*1000+1,(iu+2)*1000,deltaF(ir,iu+1))


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ genetic gain  
  ad(iu+1,1) = sum(adde1002(:,1))/(real(n1001))
!print *, 34343434

  print*, 'iu=', iu
  END DO    !iu

  DO i = 2, nu-1
      gg(ir,i) = ad(i,1) - ad(i-1,1)
  END DO    !i

 print *, 'ir=', ir
END DO    !end of reiteration 
!print *, 55555
!-------------------------------------------n_eQTL
m_neQTL = sum(n_eQTL(:,1))/real(nr)
write(999,*) 'n_eQTL=', m_neQTL

!-------------------------------------------h
m_h = sum(h(:,:))/(nr*n1001)
DO ir = 1, nr
write(9,*) 'h=', sum(h(ir,:))/n1001
END DO
write(9,*) 'm_h=', m_h

  !--------------------------Genetic Variance
   DO iu = 1, nu 
   w = 0.0
   z = 0.0
   w1 = 0.0
       m_VG(iu) = sum(Var_G(:,iu))/(real(nr))
       m_VE(iu) = sum(Var_E(:,iu))/(real(nr))
     !  m_VP(nu) = sum(Var_P(:,nu))/real(nu)
       m_rh2(iu) = sum(r_h2(:,iu))/(real(nr))
       DO ir = 1, nr
           w = w + (Var_G(ir,iu)-m_VG(iu))**2
           z = z + (r_h2(ir,iu)-m_rh2(iu))**2
           w1 = w1 + (Var_E(ir,iu)-m_VE(iu))**2
       END DO    !ir
       sd_VG(iu) = w/(real(nr-1))
       sd_rh2(iu) = z/(real(nr-1))
       sd_VE(iu) = w1/(real(nr-1))
   END DO 
   
   open(unit=711,file='VG.txt')
       DO iu = 1, nu
           write(711,*) 'iu=', iu 
           write(711,*) 'VG=', m_VG(iu), 'sd=', sd_VG(iu) 
           write(711,*) 'rh2=', m_rh2(iu), 'sd=', sd_rh2(iu) 
           write(711,*) 'VE=', m_VE(iu), 'sd=', sd_VE(iu)
       END DO    !iu
   close(711) 
   !=============================================================================   
   !writing rate_F

   DO ir = 1, nr
       DO i = 1, nu-1
       IF(i==1) THEN
           rate_F(ir,i) = deltaF(ir,i)
       ELSE
           rate_F(ir,i) = (deltaF(ir,i) - deltaF(ir,i-1))/(1-deltaF(ir,i-1))
       END IF 
       END DO    !i
   END DO    !ir
   
   !-------------------- m_F and V_F
   open(unit=93,file='f.txt')
   DO iu = 1, nu-1
       w = 0.0
       m_F(iu) = sum(deltaF(:,iu))/(real(nr))
       DO ir = 1, nr
           w = w + (deltaF(ir,iu)-m_F(iu))**2
       END DO    !ir
       sd_F(iu) = sqrt(w/(real(nr-1)))
       write(93,*) 'iu=', iu, 'f=', m_F(iu), 'sd=', sd_F(iu)
   END DO    !iu
   close(93)

   !--------------------------------
   open(unit=88,file='rate_F.txt') 
   !rate_F(1:2) = 0.0
   !sd_rate_F(1:2) = 0.0
!print *, 66666
   write(88,*) '!rate_F for artificial selection, h2=', h2
   DO i = 1, nu-1
       e = 0.0 
       mean_rate_F(i) = sum(rate_F(:,i))/(real(nr))
       DO j = 1, nr
           e = e + (rate_F(j,i)-mean_rate_F(i))**2
       END DO    !j
       sd_rate_F(i) = sqrt(e/(real(nr-1)))
       
       write(88,*) 'i=', i, 'rate=', mean_rate_F(i), 'sd=', sd_rate_F(i)
   END DO    !i
      
   close(88)

   !print *, 77777
   !============================================================================
   !writing gg
   open(unit=99,file='gg.txt')
   write(99,*) '!gg for artificial selection, h2=', h2
   DO i = 1, nu-1
       e = 0.0
       mean_gg(i) = sum(gg(:,i))/(real(nr))
       DO j = 1, nr
           e = e + (gg(j,i)-mean_gg(i))**2
       END DO    !j
       sd_gg(i) = sqrt(e/(real(nr-1)))
       write(99,*) 'iG=', i, 'mean_gg=', mean_gg(i), 'sd=', sd_gg(i)
   END DO    !i
   close(99)
   

END Program pop




!#####################################################################################################

!subroutine to calculate the inbreeding coefficent
!######################################################################################################
        SUBROUTINE inbreed(N1,N,mean)
      		!PROGRAM FOR CALCULATING INBREEDING COEFFICIENT

		IMPLICIT NONE
		REAL*8, ALLOCATABLE:: F(:),SIGMA(:)		! CONTAINS INBREEDING COEFFICIENTS
		INTEGER,DIMENSION(:),ALLOCATABLE::ANIM,SIRE,DAM
										! SIRE=IDENTIFICATION NO OF THE SIRE
										! DAM=IDENTIFICATION NO OF THE DAM
		INTEGER :: X(2100000),uu
		REAL*8 	:: T(2100000),D(2100000),FF,s, mean
		INTEGER :: INDEX(2100000)          ! INDICATES WHICH ELEMENTS IN T VECTOR WILL 
		INTEGER :: AA,SS,DD,I,J,K,JJ,KK,G,N,N1,AA1,SS1,DD1,KKK,SSS,DDD
		LOGICAL CC
		           
		OPEN(2,FILE="pedigree.txt")

		ALLOCATE (F(N),ANIM(N),SIRE(N),DAM(N),SIGMA(N))

	
		DO I=1,N1-1
                      ! print *, I
			READ(2,*) ANIM(I),SIRE(I),DAM(I),F(I),SIGMA(I)
		END DO

		DO I=N1,N
			READ(2,*) ANIM(I),SIRE(I),DAM(I)
			F(I)=0.00
			SIGMA(I)=0.00
		END DO
	        !print *, 111111	
		SSS=SIRE(N1-1)
		DDD=DAM(N1-1)
DO_I:	DO I=N1,N
			AA=I
			SS=SIRE(AA)
			DD=DAM(AA)
			CALL BBB(N,AA,SS,DD,SIGMA,F)             ! 求每个个体的孟德尔系数
OUTER:		IF(SS == SSS .AND. DD == DDD) THEN
				F(I)=F(I-1)		!FULL SIBS HAVE THE SAME INBREEDING COEFFICIENTS
			ELSE OUTER
				SSS=SS                  !作为中间代换判断SSS与SS是否相等
				DDD=DD
				K=1
				G=1
				
				CALL AAA(N,K,G,AA,SS,DD,X,T,D,F)
INNER:			IF (SS==0 .OR. DD==0) THEN
					F(K)=0.00       !SS或DD值为0时F(K)值为0
				ELSE INNER
					INDEX(K)=1
					CC=.TRUE.
					JJ=1
DO_CC:				DO WHILE (CC)           !当SS与DD都不为0时
						G=G+1
						CC=.FALSE.
						KK=0
						KKK=K
DO_J:					DO J=1,JJ                  !将个体的所有亲属的T值和D值计算出来
							AA=X(KKK-JJ+J)
							IF(K == 1)AA=X(K)
							SS=SIRE(AA)
							DD=DAM(AA)
							IF(SS /= 0) THEN
								CC=.TRUE.
								AA1=SS
								SS1=SIRE(AA1)
								DD1=DAM(AA1)
								K=K+1
								CALL AAA(N,K,G,AA1,SS1,DD1,X,T,D,F)
								INDEX(K)=1
								KK=KK+1
							ENDIF
							IF(DD /= 0) THEN
								CC=.TRUE.
								AA1=DD
								SS1=SIRE(AA1)
								DD1=DAM(AA1)
								K=K+1
								CALL AAA(N,K,G,AA1,SS1,DD1,X,T,D,F)
								INDEX(K)=1
								KK=KK+1
							ENDIF
						ENDDO DO_J
						JJ=KK
					ENDDO DO_CC

	! ACCUMULATE VALUES IN T VECTOR FOR ANIMALS WITH SAME ID
DO_J2:	DO J=1,K
			IF(INDEX(J) == 1) THEN       !J的SS和DD都存在时
				AA=X(J)
DO_JJ:			DO JJ=J+1,K
					IF(AA == X(JJ)) THEN    !意思是J的祖先与J的祖先有为相同的个体
						INDEX(JJ)=0
						T(J)=T(J)+T(JJ)  !如果是X中有个体相同的，就将T值加起来，赋值到最前面的T值中
					ENDIF
			    ENDDO DO_JJ
			 ENDIF
		ENDDO DO_J2

	!CALCULATE INBREEDING COEEFICIENTS
		  FF=0.D0
		  DO J=1,K
			 FF=FF+T(J)*T(J)*D(J)*INDEX(J)
		  ENDDO
		  F(I)=FF-1
		ENDIF INNER

		ENDIF OUTER

		ENDDO DO_I
			REWIND (2)
		DO I=1,N
			WRITE(2,'(3I7,F10.6,F10.6)') ANIM(I),SIRE(I),DAM(I),F(I),SIGMA(I)
		ENDDO
               
               s = 0.0
               DO i = N1,N
                  s = s + F(i) 
               END DO     
               mean = s/real(N-N1+1)
        CLOSE(2)
  	
		
		DEALLOCATE (F,ANIM,SIRE,DAM,SIGMA)
		END SUBROUTINE inbreed

		SUBROUTINE AAA(N,K,G,AA,SS,DD,X,T,D,F)
			INTEGER :: N,K,G,AA,SS,DD,X(2100000)
			REAL*8 :: T(2100000),D(2100000),F(N)

			X(K)=AA
			T(K)=.5**(G-1)              !计算T，个体与X(K)差N代，T值就为0.5的N次方
			IF(SS/=0 .AND. DD/=0) THEN
				D(K)=.5-.25*(F(SS)+F(DD))   !计算D值	
			ELSEIF(SS/=0 .AND. DD==0) THEN
				D(K)=.75-.25*F(SS)
			ELSEIF(SS==0 .AND. DD/=0) THEN
				D(K)=.75-.25*F(DD)
			ELSE
				D(K)=1.00
			ENDIF
		END SUBROUTINE AAA

        SUBROUTINE BBB(N,AA,SS,DD,SIGMA,F)  ! 求孟德尔抽样系数子程序
		INTEGER :: N,AA,SS,DD
		REAL*8 :: SIGMA(N),F(N)
		IF(SS/=0 .AND. DD/=0) THEN
		SIGMA(AA)=0.25*(1-F(SS))*0.3+0.25*(1-F(DD))*0.3
        ELSEIF(SS/=0 .AND. DD==0) THEN
        SIGMA(AA)=0.25*(1-F(SS))*0.3+0.25*0.3
        ELSEIF(SS==0 .AND. DD/=0) THEN
        SIGMA(AA)=0.25*0.3+0.25*(1-F(DD))*0.3
        ELSE
		SIGMA(AA)=0.15
		END IF
		END SUBROUTINE BBB






 !########################################################################################################



 !###########################################################################################################

SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
          
    CALL SYSTEM_CLOCK(COUNT=clock)
          
    seed = clock  +  37  *  (/ (i  -  1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
          
    DEALLOCATE(seed)
END SUBROUTINE

!SUBROUTINE uniform distribution UNI(lb,ub)
SUBROUTINE UNI(n,lb,ub,y)
    IMPLICIT NONE 
	REAL(kind=8) :: lb, ub, u
	INTEGER :: n, i
	REAL(kind=8):: y(n)
	

	
	DO i = 1, n
	call random_number(u)
	    y(i) = u*(ub-lb) + lb 
	END DO



END SUBROUTINE


!SUBROUTINE normal distribution N(mu,sigma_sq) by rejection sampling
SUBROUTINE NOA(n,mu,sigma_sq,y)
   IMPLICIT NONE								
   INTEGER :: i, n
   REAL(kind=8)::mu,sigma_sq,y(n),y1(n),y2(n),u1,u2
   REAL(kind=8),parameter::pi=3.1415926

   DO i = 1, n
1  CALL random_number(u1)
   CALL random_number(u2)
   y1(i)=20*u1 - 10    ! a = 4
   y2(i)=u2/sqrt(2*pi)
   IF(y2(i) > exp(-(y1(i)**2)/2)/sqrt(2*pi)) goto 1

   y(i)=mu + sqrt(sigma_sq)*y1(i)
   END DO

   RETURN
END SUBROUTINE


!SUROUTINE Poisson distribution with lambda 
SUBROUTINE DPOI(n,lambda,y)
    IMPLICIT NONE
    INTEGER(kind=4) :: x, i, n, y(n)
	REAL(kind=8) :: m, lambda, u

	
    DO i = 1, n
	x = -1
	m = exp(lambda)
    
10	CALL random_number(u)
	x = x+1
	m = m*u
	IF(m < 1) then
	    y(i) = x 
	Else
	    goto 10
        END IF
        
        END DO


END SUBROUTINE


!subroutine of Haldane's mapping function
SUBROUTINE hmf(M,r)
IMPLICIT NONE 
REAL(kind=8) :: M, r        !morgan and combination rate 

r = 0.5*(1-exp(-2*M))
RETURN 

END SUBROUTINE


!SUBROUTINE selection_sort for 2 dimensional array,given the m as the target column, and the largest one is the last one
SUBROUTINE selection_sort(a,n,m)
   IMPLICIT NONE
   INTEGER:: n,m
   INTEGER::i,j  !loop counter
   REAL(kind=8)::mi(1,2),a(n,2)
   REAL(kind=8)::temp(1,2) ! swapping variable

   
   DO i=1,n
        mi(1,:)=a(i,:)  !let a(i,m) be the smallest one
        DO j=i + 1, n
            IF(mi(1,m) > a(j,m)) THEN   !find that a(i,m) is not the smallest one
                temp(1,:)=a(j,:)
                a(j,:)=a(i,:)
                a(i,:)=temp(1,:)
                mi(1,:)=a(i,:)
            END IF
        END DO
   END DO
   RETURN
END SUBROUTINE

SUBROUTINE gamma(alpha,beta,gm)
IMPLICIT NONE


INTEGER ::  i
REAL(kind=8) :: alpha, beta
REAL(kind=8) :: u1, u2, x, y, b, y1, gm
REAL(kind=8), PARAMETER :: e = 2.718281828 

    b = (e+alpha)/e
   

10  CALL random_number(u1)
    CALL random_number(u2)

    y = b*u1
    

    IF (y<=1) THEN
        x = y**(1/alpha)
        IF(u2<=exp((-1)*x)) THEN
            y1 = x
        ELSE
            goto 10
        END IF
    ELSE 
        x = (-1)*log((b-y)/alpha) 
        IF(u2<= x**(alpha-1)) THEN
            y1 = x
        ELSE
            goto 10
        END IF
    END IF

    gm = y1/beta
    
   
    

END SUBROUTINE 
 


subroutine sort(a,n)
  implicit none
  integer :: n
  real(kind=8)::a(n,2)
  integer i,j  ! 循环计数器
  real(kind=8):: ma(1,2)  ! 找出每一轮中的最大值
  real(kind=8):: temp(1,2) ! 交换数据时使用
  do i=1,n-1
    ma(1,:)=a(i,:)     ! 暂时令a(i)是最大值
    do j=i+1,n
      if ( ma(1,2) < a(j,2) ) then   ! 发现a(i)不是最大
        temp(1,:)=a(j,:)        ! 把a(i)\a(j)交换
        a(j,:)=a(i,:)
        a(i,:)=temp(1,:)
        ma(1,:)=a(i,:)
      end if
 end do
  end do                             
  return
end subroutine                
