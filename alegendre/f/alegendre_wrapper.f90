!TODO: write docs, using doxygen if possible
!
!Fortran code
!This is a wrapper of some functions from Bremer's code, designed to enable C binding
!




SUBROUTINE alegendre_eval_init_wrapper(dsize) BIND(C)
	USE ISO_C_BINDING
	USE alegendreeval
	
	REAL (C_DOUBLE), INTENT(OUT) 	:: dsize
	
	CALL alegendre_eval_init (dsize)
	RETURN
END SUBROUTINE 


SUBROUTINE alegendre_eval_wrapper (dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq) BIND(C)
	USE ISO_C_BINDING
	USE alegendreeval

	REAL (C_DOUBLE), VALUE 			:: dnu
	REAL (C_DOUBLE), VALUE 			:: dmu
	REAL (C_DOUBLE), VALUE 			:: t
	
	REAL (C_DOUBLE), INTENT(OUT) 	:: alpha
	REAL (C_DOUBLE), INTENT(OUT) 	:: alphader
	REAL (C_DOUBLE), INTENT(OUT) 	:: vallogp
	REAL (C_DOUBLE), INTENT(OUT) 	:: vallogq
	REAL (C_DOUBLE), INTENT(OUT) 	:: valp
	REAL (C_DOUBLE), INTENT(OUT) 	:: valq
	
	CALL alegendre_eval (dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_nroots_wrapper (dnu, dmu, nproots, nqroots) BIND(C)
	USE ISO_C_BINDING
	USE alegendreeval
	
	REAL (C_DOUBLE), VALUE 			:: dnu
	REAL (C_DOUBLE), VALUE 			:: dmu
	
	INTEGER (C_LONG), INTENT(OUT) 	:: nproots
	INTEGER (C_LONG), INTENT(OUT) 	:: nqroots
	
	CALL alegendre_nroots (dnu, dmu, nproots, nqroots)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_proot_wrapper (dnu, dmu, j, t) BIND(C)
	USE ISO_C_BINDING
	USE alegendreeval
	
	REAL (C_DOUBLE), VALUE 			:: dnu
	REAL (C_DOUBLE), VALUE 			:: dmu
	INTEGER (C_LONG), VALUE 		:: j
	
	REAL (C_DOUBLE), INTENT(OUT) 	:: t
	
	CALL alegendre_proot (dnu, dmu, j, t)
RETURN
END SUBROUTINE