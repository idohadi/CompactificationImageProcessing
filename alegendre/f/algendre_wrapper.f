!TODO: write docs, using doxygen if possible
!
!Fortran code
!This is a wrapper of some functions from Bremer's code, designed to enable C binding
!

use extern/ALegendreEval/alegendreeval


SUBROUTINE alegendre_eval_init_wrapper (dsize)
	!TODO
	
	CALL alegendre_eval_init (dsize)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_eval_wrapper (dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq)
	!TODO
	
	CALL alegendre_eval (dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_nroots_wrapper (dnu, dmu, nproots, nqroots)
	!TODO
	
	CALL alegendre_nroots (dnu, dmu, nproots, nqroots)
RETURN
END SUBROUTINE


SUBROUTINE alegendre_proot_wrapper (dnu, dmu, j, t)
	!TODO
	
	CALL alegendre_proot (dnu, dmu, j, t)
RETURN
END SUBROUTINE
