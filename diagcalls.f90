!-----------------------------------------------------------------------------!
           SUBROUTINE CALL_DSYEV_NO_EIGVECS(H,ndim,eigvals)
!-----------------------------------------------------------------------------!
            implicit none

        integer      (kind = 4)            :: ndim

        real         (kind = 8)            :: H(ndim,ndim) 
        real         (kind = 8)            :: eigvals(ndim)
      
        character    (len  = 1),parameter  :: JOBZ="N"
        character    (len  = 1),parameter  :: UPLO="L"
        integer      (kind = 4)            :: LDA
        integer      (kind = 4)            :: LDWORK
        integer      (kind = 4)            :: INFO
        real         (kind = 8),allocatable:: WORK(:)

      LDA = ndim
      LDWORK=3*ndim-1
      allocate(WORK(LDWORK))
      CALL DSYEV(JOBZ,UPLO,ndim,H,LDA,EIGVALS,WORK,LDWORK,info) 
      deallocate(WORK)

                    END
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
           SUBROUTINE CALL_DSYEV_WT_EIGVECS(H,ndim,eigvals,eigvecs)
!-----------------------------------------------------------------------------!
            implicit none

        integer      (kind = 4)            :: ndim

        real         (kind = 8)            :: H(ndim,ndim) 
        real         (kind = 8)            :: eigvals(ndim)
        real         (kind = 8)            :: eigvecs(ndim,ndim)
      
        character    (len  = 1),parameter  :: JOBZ="V"
        character    (len  = 1),parameter  :: UPLO="L"
        integer      (kind = 4)            :: LDA
        integer      (kind = 4)            :: LDWORK
        integer      (kind = 4)            :: INFO
        real         (kind = 8),allocatable:: WORK(:)

      LDA = ndim
      LDWORK=3*ndim-1
      allocate(WORK(LDWORK))
      CALL DSYEV(JOBZ,UPLO,ndim,H,LDA,EIGVALS,WORK,LDWORK,info) 
      deallocate(WORK)
      eigvecs=H
                    END
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
           SUBROUTINE CALL_DSYEVX_WT_EIGVECS(H,ndim,eigvals,eigvecs,il,iu)
!-----------------------------------------------------------------------------!
            implicit none

        integer      (kind = 4)            :: ndim

        real         (kind = 8)            :: H(ndim,ndim) 
        real         (kind = 8)            :: eigvals(ndim)
        real         (kind = 8)            :: eigvecs(ndim,ndim)
      
        character    (len  = 1),parameter  :: JOBZ="V"
        character    (len  = 1),parameter  :: UPLO="L"
        character    (len  = 1),parameter  :: RANGEM="I"
        integer      (kind = 4)            :: LDA
        integer      (kind = 4)            :: LDWORK
        integer      (kind = 4)            :: INFO
        real         (kind = 8),allocatable:: WORK(:)
        real         (kind = 8)            :: vl,vu
        integer      (kind = 4)            :: il,iu
        integer      (kind = 4)            :: m
        real         (kind = 8)            :: abstol=1.0e-6
        real         (kind = 8)            :: Z(ndim,iu-il+1)
        integer      (kind = 4)            :: iwork(5*ndim),ifail(ndim),ldz

      m= il+iu+1
      LDA = ndim
      ldz=ndim
      LDWORK=8*ndim
      allocate(WORK(LDWORK))
      CALL DSYEVX(JOBZ,RANGEM,UPLO,ndim,H,LDA,VL,VU,IL,IU,ABSTOL,M,EIGVALS,Z&
                 &,LDZ,WORK,LDWORK,IWORK,IFAIL,info) 
      deallocate(WORK)
      eigvecs(1:ndim,1:iu-il+1)=Z
                    END
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
           SUBROUTINE CALL_ZGEEV_NO_EIGVECS(H,ndim,eigvals)
!-----------------------------------------------------------------------------!
            implicit none

        integer      (kind = 4)            :: ndim

        COMPLEX      (kind = 8)            :: H(ndim,ndim) 
        COMPLEX      (kind = 8)            :: eigvals(ndim)
        COMPLEX      (kind = 8)            :: dum(ndim,ndim)
      
        character    (len  = 1),parameter  :: JOBVL="N"
        character    (len  = 1),parameter  :: JOBVR="N"
        integer      (kind = 4)            :: LDA,LDVL,LDVR,LDWORK
        integer      (kind = 4)            :: INFO
        COMPLEX         (kind = 8),allocatable:: WORK(:),WORK2(:)

      LDA = ndim
      LDVL= ndim
      LDVR= ndim
      LDWORK=2*ndim
      allocate(WORK(LDWORK),WORK2(LDWORK))
       CALL ZGEEV(JOBVL,JOBVR,ndim,H,LDA,eigvals,dum,LDVL,dum,LDVR,&
            & WORK, LDWORK, WORK2, INFO)
      deallocate(WORK,WORK2)

                    END
!-----------------------------------------------------------------------------!
           SUBROUTINE CALL_ZGEEV_WT_EIGVECS(H,ndim,eigvals,vecsL,vecsR)
!-----------------------------------------------------------------------------!
            implicit none

        integer      (kind = 4)            :: ndim

        COMPLEX      (kind = 8)            :: H(ndim,ndim) 
        COMPLEX      (kind = 8)            :: eigvals(ndim)
        COMPLEX      (kind = 8)            :: vecsl(ndim,ndim)
        COMPLEX      (kind = 8)            :: vecsr(ndim,ndim)
      
        character    (len  = 1),parameter  :: JOBVL="V"
        character    (len  = 1),parameter  :: JOBVR="V"
        integer      (kind = 4)            :: LDA,LDVL,LDVR,LDWORK
        integer      (kind = 4)            :: INFO
        COMPLEX         (kind = 8),allocatable:: WORK(:),WORK2(:)

      LDA = ndim
      LDVL= ndim
      LDVR= ndim
      LDWORK=2*ndim
      allocate(WORK(LDWORK),WORK2(LDWORK))
       CALL ZGEEV(JOBVL,JOBVR,ndim,H,LDA,eigvals,vecsl,LDVL,vecsr,LDVR,&
            & WORK, LDWORK, WORK2, INFO)
      deallocate(WORK,WORK2)

                    END
!-----------------------------------------------------------------------------!
