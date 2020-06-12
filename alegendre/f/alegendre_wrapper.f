!TODO: write docs, using doxygen if possible
!
!Fortran code
!This is a wrapper of some functions from Bremer's code, designed to enable C binding
!

USE extern/ALegendreEval/alegendreeval
USE ISO_C_BINDING




SUBROUTINE alegendre_eval_init_wrapper (dsize)
	IMPLICIT NONE
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: dsize
	
	CALL alegendre_eval_init (dsize)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_eval_wrapper (dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq)
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: dnu
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: dmu
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: t
	
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: alpha
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: alphader
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: vallogp
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: vallogq
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: valp
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: valq
	
	CALL alegendre_eval (dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_nroots_wrapper (dnu, dmu, nproots, nqroots)
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: dnu
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: dmu
	
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: nproots
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: nqroots
	
	CALL alegendre_nroots (dnu, dmu, nproots, nqroots)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_proot_wrapper (dnu, dmu, j, t)
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: dnu
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: dmu
	DOUBLE PRECISION (C_DOUBLE), VALUE 			:: j
	
	DOUBLE PRECISION (C_DOUBLE), INTENT(OUT) 	:: t
	
	CALL alegendre_proot (dnu, dmu, j, t)
RETURN
END SUBROUTINE
