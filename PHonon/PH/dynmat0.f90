!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dynmat0_new
  !-----------------------------------------------------------------------
  !! This routine computes the part of the dynamical matrix which
  !! does not depend upon the change of the Bloch wavefunctions.
  !! It is a driver which calls the routines \(\texttt{dynmat_##}\) and
  !! \(\texttt{d2ionq}\) to compute the electronic part and
  !! the ionic part respectively.
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat,ntyp => nsp, ityp, zv, tau
  USE cell_base,     ONLY : alat, omega, at, bg
  USE gvect,         ONLY : g, gg, ngm, gcutm
  USE symm_base,     ONLY : irt, s, invs
  USE control_flags, ONLY : modenum, llondon, lxdm, ldftd3
  USE ph_restart,    ONLY : ph_writefile
  USE control_ph,    ONLY : current_iq
  USE control_lr,    ONLY : rec_code_read
  USE qpoint,        ONLY : xq
  USE modes,         ONLY : u, nmodes
  USE partial,       ONLY : done_irr, comp_irr
  USE dynmat,        ONLY : dyn, dyn00, dyn_rec
  USE lr_symm_base,  ONLY : minus_q, irotmq, nsymq, rtau
  USE ldaU,          ONLY : lda_plus_u
  USE io_global,     ONLY : stdout

  implicit none

  integer :: nu_i, nu_j, na_icart, nb_jcart, ierr
  integer :: isym, na, nb, icart, jcart, mu, nu, i
  ! counters

  complex(DP) :: wrk, dynwrk (3 * nat, 3 * nat)
  ! auxiliary space

  complex(DP), allocatable :: dyn_prev(:,:), dyn_now(:,:)
  complex(DP), allocatable :: dyn_before_sym(:,:), dyn_after_sym(:,:)
  complex(DP), allocatable :: dyn_after_pattern(:,:)
  complex(DP) :: block(3,3)

  real(DP) :: norm_prev, norm_now, norm_diff, norm_rel
  real(DP) :: norm_herm
  real(DP) :: norm_sym_change, norm_pattern_change
  real(DP) :: block_norm

  IF ( .NOT. comp_irr(0) .or. done_irr(0) ) RETURN
  IF (rec_code_read > -30 ) RETURN

  call start_clock ('dynmat0')

  allocate(dyn_prev(size(dyn,1), size(dyn,2)))
  allocate(dyn_now (size(dyn,1), size(dyn,2)))

  WRITE(stdout,*) ' '
  WRITE(stdout,*) '================ DEBUG dynmat0_new ================'
  WRITE(stdout,*) 'DEBUG dynmat0: current_iq = ', current_iq
  WRITE(stdout,*) 'DEBUG dynmat0: nat        = ', nat
  WRITE(stdout,*) 'DEBUG dynmat0: nmodes     = ', nmodes
  WRITE(stdout,*) 'DEBUG dynmat0: nsymq      = ', nsymq
  WRITE(stdout,*) 'DEBUG dynmat0: modenum    = ', modenum
  WRITE(stdout,*) 'DEBUG dynmat0: minus_q    = ', minus_q
  WRITE(stdout,*) 'DEBUG dynmat0: irotmq     = ', irotmq
  WRITE(stdout,*) 'DEBUG dynmat0: llondon    = ', llondon
  WRITE(stdout,*) 'DEBUG dynmat0: lxdm       = ', lxdm
  WRITE(stdout,*) 'DEBUG dynmat0: ldftd3     = ', ldftd3
  WRITE(stdout,*) 'DEBUG dynmat0: lda_plus_u = ', lda_plus_u
  WRITE(stdout,'(A,3ES18.8)') 'DEBUG dynmat0: xq = ', xq(1), xq(2), xq(3)

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: small-group operations'
  DO isym = 1, nsymq
     WRITE(stdout,*) 'DEBUG dynmat0: ---- isym = ', isym, ' ----'
     WRITE(stdout,*) 'DEBUG dynmat0: s(:,:,isym)'
     DO i = 1, 3
        WRITE(stdout,'(3I8)') s(i,1,isym), s(i,2,isym), s(i,3,isym)
     ENDDO

     WRITE(stdout,*) 'DEBUG dynmat0: irt(isym,na)'
     DO na = 1, nat
        WRITE(stdout,'(A,I5,A,I5)') 'DEBUG dynmat0: atom ', na, ' -> ', &
             irt(isym,na)
     ENDDO

     WRITE(stdout,*) 'DEBUG dynmat0: rtau(:,isym,na)'
     DO na = 1, nat
        WRITE(stdout,'(A,I5,A,3ES18.8)') 'DEBUG dynmat0: atom ', na, &
             ' rtau = ', rtau(1,isym,na), rtau(2,isym,na), &
             rtau(3,isym,na)
     ENDDO
  ENDDO

  !
  ! Start from dyn00
  !
  call zcopy (9 * nat * nat, dyn00, 1, dyn, 1)

  dyn_prev(:,:) = dyn(:,:)

  norm_now = sqrt(sum(abs(dyn(:,:))**2))
  norm_herm = sqrt(sum(abs(dyn(:,:) - conjg(transpose(dyn(:,:))))**2))

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: after zcopy dyn00 -> dyn'
  WRITE(stdout,*) 'DEBUG dynmat0: ||dyn|| = ', norm_now
  WRITE(stdout,*) 'DEBUG dynmat0: hermiticity = ', norm_herm

  !
  ! first electronic contribution arising from the term  <psi|d2v|psi>
  !
  call dynmat_us()

  dyn_now(:,:) = dyn(:,:)
  norm_prev = sqrt(sum(abs(dyn_prev(:,:))**2))
  norm_now  = sqrt(sum(abs(dyn_now(:,:))**2))
  norm_diff = sqrt(sum(abs(dyn_now(:,:) - dyn_prev(:,:))**2))
  norm_herm = sqrt(sum(abs(dyn_now(:,:) - conjg(transpose(dyn_now(:,:))))**2))
  norm_rel  = 0.d0
  IF (norm_prev > 0.d0) norm_rel = norm_diff / norm_prev

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: after dynmat_us'
  WRITE(stdout,*) 'DEBUG dynmat0: ||dyn|| = ', norm_now
  WRITE(stdout,*) 'DEBUG dynmat0: ||delta dynmat_us|| = ', norm_diff
  WRITE(stdout,*) 'DEBUG dynmat0: relative delta dynmat_us = ', norm_rel
  WRITE(stdout,*) 'DEBUG dynmat0: hermiticity = ', norm_herm

  dyn_prev(:,:) = dyn_now(:,:)

  !
  !   Here the ionic contribution
  !
  call d2ionq (nat, ntyp, ityp, zv, tau, alat, omega, xq, at, bg, g, &
       gg, ngm, gcutm, nmodes, u, dyn)

  dyn_now(:,:) = dyn(:,:)
  norm_prev = sqrt(sum(abs(dyn_prev(:,:))**2))
  norm_now  = sqrt(sum(abs(dyn_now(:,:))**2))
  norm_diff = sqrt(sum(abs(dyn_now(:,:) - dyn_prev(:,:))**2))
  norm_herm = sqrt(sum(abs(dyn_now(:,:) - conjg(transpose(dyn_now(:,:))))**2))
  norm_rel  = 0.d0
  IF (norm_prev > 0.d0) norm_rel = norm_diff / norm_prev

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: after d2ionq'
  WRITE(stdout,*) 'DEBUG dynmat0: ||dyn|| = ', norm_now
  WRITE(stdout,*) 'DEBUG dynmat0: ||delta d2ionq|| = ', norm_diff
  WRITE(stdout,*) 'DEBUG dynmat0: relative delta d2ionq = ', norm_rel
  WRITE(stdout,*) 'DEBUG dynmat0: hermiticity = ', norm_herm

  dyn_prev(:,:) = dyn_now(:,:)

  !
  ! Contribution from the dispersion correction
  !
  IF (llondon .OR. lxdm) THEN
     CALL d2ionq_disp(alat,nat,ityp,at,bg,tau,xq,dynwrk)
     CALL rotate_pattern_add (nat,u,dyn,dynwrk)

     dyn_now(:,:) = dyn(:,:)
     norm_prev = sqrt(sum(abs(dyn_prev(:,:))**2))
     norm_now  = sqrt(sum(abs(dyn_now(:,:))**2))
     norm_diff = sqrt(sum(abs(dyn_now(:,:) - dyn_prev(:,:))**2))
     norm_herm = sqrt(sum(abs(dyn_now(:,:) - conjg(transpose(dyn_now(:,:))))**2))
     norm_rel  = 0.d0
     IF (norm_prev > 0.d0) norm_rel = norm_diff / norm_prev

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: after d2ionq_disp + rotate_pattern_add'
     WRITE(stdout,*) 'DEBUG dynmat0: ||dyn|| = ', norm_now
     WRITE(stdout,*) 'DEBUG dynmat0: ||delta dispersion|| = ', norm_diff
     WRITE(stdout,*) 'DEBUG dynmat0: relative delta dispersion = ', norm_rel
     WRITE(stdout,*) 'DEBUG dynmat0: hermiticity = ', norm_herm

     dyn_prev(:,:) = dyn_now(:,:)

  ELSEIF(ldftd3) THEN
     CALL d2ionq_dispd3(alat,nat,at,xq,dynwrk)
     CALL rotate_pattern_add(nat,u,dyn,dynwrk)

     dyn_now(:,:) = dyn(:,:)
     norm_prev = sqrt(sum(abs(dyn_prev(:,:))**2))
     norm_now  = sqrt(sum(abs(dyn_now(:,:))**2))
     norm_diff = sqrt(sum(abs(dyn_now(:,:) - dyn_prev(:,:))**2))
     norm_herm = sqrt(sum(abs(dyn_now(:,:) - conjg(transpose(dyn_now(:,:))))**2))
     norm_rel  = 0.d0
     IF (norm_prev > 0.d0) norm_rel = norm_diff / norm_prev

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: after d2ionq_dispd3 + rotate_pattern_add'
     WRITE(stdout,*) 'DEBUG dynmat0: ||dyn|| = ', norm_now
     WRITE(stdout,*) 'DEBUG dynmat0: ||delta DFT-D3|| = ', norm_diff
     WRITE(stdout,*) 'DEBUG dynmat0: relative delta DFT-D3 = ', norm_rel
     WRITE(stdout,*) 'DEBUG dynmat0: hermiticity = ', norm_herm

     dyn_prev(:,:) = dyn_now(:,:)

  ELSE
     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: no dispersion correction in dynmat0_new'
  ENDIF

  !
  !   Add non-linear core-correction (NLCC) contribution (if any)
  !
  call dynmatcc()

  dyn_now(:,:) = dyn(:,:)
  norm_prev = sqrt(sum(abs(dyn_prev(:,:))**2))
  norm_now  = sqrt(sum(abs(dyn_now(:,:))**2))
  norm_diff = sqrt(sum(abs(dyn_now(:,:) - dyn_prev(:,:))**2))
  norm_herm = sqrt(sum(abs(dyn_now(:,:) - conjg(transpose(dyn_now(:,:))))**2))
  norm_rel  = 0.d0
  IF (norm_prev > 0.d0) norm_rel = norm_diff / norm_prev

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: after dynmatcc'
  WRITE(stdout,*) 'DEBUG dynmat0: ||dyn|| = ', norm_now
  WRITE(stdout,*) 'DEBUG dynmat0: ||delta dynmatcc|| = ', norm_diff
  WRITE(stdout,*) 'DEBUG dynmat0: relative delta dynmatcc = ', norm_rel
  WRITE(stdout,*) 'DEBUG dynmat0: hermiticity = ', norm_herm

  dyn_prev(:,:) = dyn_now(:,:)

  !
  ! DFPT+U: calculate the bare Hubbard dynamical matrix
  !
  IF (lda_plus_u) THEN
     CALL dynmat_hub_bare()

     dyn_now(:,:) = dyn(:,:)
     norm_prev = sqrt(sum(abs(dyn_prev(:,:))**2))
     norm_now  = sqrt(sum(abs(dyn_now(:,:))**2))
     norm_diff = sqrt(sum(abs(dyn_now(:,:) - dyn_prev(:,:))**2))
     norm_herm = sqrt(sum(abs(dyn_now(:,:) - conjg(transpose(dyn_now(:,:))))**2))
     norm_rel  = 0.d0
     IF (norm_prev > 0.d0) norm_rel = norm_diff / norm_prev

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: after dynmat_hub_bare'
     WRITE(stdout,*) 'DEBUG dynmat0: ||dyn|| = ', norm_now
     WRITE(stdout,*) 'DEBUG dynmat0: ||delta Hubbard bare|| = ', norm_diff
     WRITE(stdout,*) 'DEBUG dynmat0: relative delta Hubbard bare = ', norm_rel
     WRITE(stdout,*) 'DEBUG dynmat0: hermiticity = ', norm_herm

     dyn_prev(:,:) = dyn_now(:,:)
  ENDIF

  !
  !   Symmetrizes the dynamical matrix w.r.t. the small group of q and of
  !   mode. This is done here, because this part of the dynmical matrix is
  !   saved with recover and in the other runs the symmetry group might change
  !
  IF (modenum .ne. 0) THEN

     allocate(dyn_before_sym(size(dyn,1), size(dyn,2)))
     allocate(dyn_after_sym (size(dyn,1), size(dyn,2)))
     allocate(dyn_after_pattern(size(dyn,1), size(dyn,2)))

     dyn_before_sym(:,:) = dyn(:,:)

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: modenum != 0, before symdyn_munu_new'
     WRITE(stdout,*) 'DEBUG dynmat0: ||dyn_before_sym|| = ', &
          sqrt(sum(abs(dyn_before_sym(:,:))**2))

     call symdyn_munu_new (dyn, u, xq, s, invs, rtau, irt, at, bg, &
          nsymq, nat, irotmq, minus_q)

     dyn_after_sym(:,:) = dyn(:,:)

     norm_sym_change = sqrt(sum(abs(dyn_after_sym(:,:) - &
          dyn_before_sym(:,:))**2))

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: after symdyn_munu_new'
     WRITE(stdout,*) 'DEBUG dynmat0: ||dyn_after_sym|| = ', &
          sqrt(sum(abs(dyn_after_sym(:,:))**2))
     WRITE(stdout,*) 'DEBUG dynmat0: ||after_sym - before_sym|| = ', &
          norm_sym_change
     IF (sqrt(sum(abs(dyn_before_sym(:,:))**2)) > 0.d0) THEN
        WRITE(stdout,*) 'DEBUG dynmat0: relative symdyn change = ', &
             norm_sym_change / sqrt(sum(abs(dyn_before_sym(:,:))**2))
     ENDIF

     !
     ! rotate again in the pattern basis
     !
     call zcopy (9 * nat * nat, dyn, 1, dynwrk, 1)

     dyn=(0.d0, 0.d0)

     CALL rotate_pattern_add(nat, u, dyn, dynwrk)

     dyn_after_pattern(:,:) = dyn(:,:)

     norm_pattern_change = sqrt(sum(abs(dyn_after_pattern(:,:) - &
          dyn_after_sym(:,:))**2))

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: after rotate_pattern_add back to pattern basis'
     WRITE(stdout,*) 'DEBUG dynmat0: ||dyn_after_pattern|| = ', &
          sqrt(sum(abs(dyn_after_pattern(:,:))**2))
     WRITE(stdout,*) 'DEBUG dynmat0: ||after_pattern - after_sym|| = ', &
          norm_pattern_change
     IF (sqrt(sum(abs(dyn_after_sym(:,:))**2)) > 0.d0) THEN
        WRITE(stdout,*) 'DEBUG dynmat0: relative pattern-basis change = ', &
             norm_pattern_change / sqrt(sum(abs(dyn_after_sym(:,:))**2))
     ENDIF

     deallocate(dyn_before_sym)
     deallocate(dyn_after_sym)
     deallocate(dyn_after_pattern)

  ELSE
     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG dynmat0: modenum == 0'
     WRITE(stdout,*) 'DEBUG dynmat0: symdyn_munu_new is NOT called in dynmat0_new'
     WRITE(stdout,*) 'DEBUG dynmat0: this part is saved to dyn_rec as-is'
  ENDIF

  !
  ! Dump final index blocks before saving dyn_rec.
  ! Warning: at this stage dyn is the matrix in the internal pattern/mode
  ! representation used in PH, not necessarily a pure Cartesian 3x3
  ! force-constant block.
  !
  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: final dyn before dyn_rec save'
  WRITE(stdout,*) 'DEBUG dynmat0: ||dyn_final|| = ', &
       sqrt(sum(abs(dyn(:,:))**2))
  WRITE(stdout,*) 'DEBUG dynmat0: hermiticity final = ', &
       sqrt(sum(abs(dyn(:,:) - conjg(transpose(dyn(:,:))))**2))

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: final 3x3 index blocks of dyn'
  WRITE(stdout,*) 'DEBUG dynmat0: caution: these are index blocks, not necessarily Cartesian blocks'

  DO na = 1, nat
     DO nb = 1, nat

        DO icart = 1, 3
           DO jcart = 1, 3
              mu = 3 * (na - 1) + icart
              nu = 3 * (nb - 1) + jcart
              block(icart,jcart) = dyn(mu,nu)
           ENDDO
        ENDDO

        block_norm = sqrt(sum(abs(block(:,:))**2))

        WRITE(stdout,*) ' '
        WRITE(stdout,'(A,I5,A,I5,A,ES18.8)') &
             'DEBUG dynmat0: block na = ', na, ' nb = ', nb, &
             ' norm = ', block_norm

        WRITE(stdout,*) 'DEBUG dynmat0: block rows: Re(1) Im(1) Re(2) Im(2) Re(3) Im(3)'
        DO icart = 1, 3
           WRITE(stdout,'(6ES18.8)') &
                real(block(icart,1)), aimag(block(icart,1)), &
                real(block(icart,2)), aimag(block(icart,2)), &
                real(block(icart,3)), aimag(block(icart,3))
        ENDDO

     ENDDO
  ENDDO

  !
  !      call tra_write_matrix('dynmat0 dyn',dyn,u,nat)
  !
  dyn_rec(:,:)=dyn(:,:)

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG dynmat0: dyn_rec has been updated from dyn'
  WRITE(stdout,*) 'DEBUG dynmat0: ||dyn_rec|| = ', &
       sqrt(sum(abs(dyn_rec(:,:))**2))
  WRITE(stdout,*) '================ END DEBUG dynmat0_new ================'
  WRITE(stdout,*) ' '

  done_irr(0) = .TRUE.
  CALL ph_writefile('data_dyn',current_iq,0,ierr)

  deallocate(dyn_prev)
  deallocate(dyn_now)

  call stop_clock ('dynmat0')
  return
end subroutine dynmat0_new