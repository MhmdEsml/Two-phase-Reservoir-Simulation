        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 22 16:50:53 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INV__genmod
          INTERFACE 
            FUNCTION INV(A) RESULT(AINV)
              REAL(KIND=4), INTENT(IN) :: A(:,:)
              REAL(KIND=4) :: AINV(SIZE(A,1),SIZE(A,2))
            END FUNCTION INV
          END INTERFACE 
        END MODULE INV__genmod
