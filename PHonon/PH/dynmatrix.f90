!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dynmatrix_new(iq_)
  !-----------------------------------------------------------------------
  !! This routine is a driver which computes the symmetrized dynamical
  !! matrix at q (and in the star of q) and diagonalizes it.  
  !! It writes the result on a iudyn file and writes the eigenvalues on
  !! output.
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : FPI, BOHR_RADIUS_ANGS
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm, amass, zv
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : modenum
  USE cell_base,     ONLY : at, bg, celldm, ibrav, omega
  USE symm_base,     ONLY : s, sr, irt, nsym, invs, t_rev
  USE dynmat,        ONLY : dyn, w2
  USE noncollin_module, ONLY : nspin_mag
  USE modes,         ONLY : u, nmodes, npert, nirr, num_rap_mode
  USE gamma_gamma,   ONLY : nasr, asr, equiv_atoms, has_equivalent, &
                            n_diff_sites
  USE efield_mod,    ONLY : epsilon, zstareu, zstarue0, zstarue
  USE disp,          ONLY : omega_disp
  USE control_ph,    ONLY : epsil, zue, search_sym, ldisp, &
                            done_zue, always_run, ldiag, done_epsil, done_zeu, xmldyn, &
                            current_iq, qplot
  USE ph_restart,    ONLY : ph_writefile
  USE partial,       ONLY : all_comp, comp_irr, done_irr, nat_todo_input
  USE units_ph,      ONLY : iudyn
  USE noncollin_module, ONLY : m_loc, nspin_mag
  USE output,        ONLY : fildyn
  USE io_dyn_mat,    ONLY : write_dyn_mat_header
  USE ramanm,        ONLY : lraman, ramtns

  USE lr_symm_base,  ONLY : minus_q, irotmq, nsymq, rtau
  USE qpoint,        ONLY : xq
  USE control_lr,    ONLY : lgamma, lgamma_gamma, where_rec, rec_code

  implicit none
  INTEGER, INTENT(IN) :: iq_
  !
  ! ... local variables
  !
  integer :: nq, isq (48), imq, na, nb, nt, imode0, jmode0, irr, jrr, &
       ipert, jpert, mu, nu, i, j, nqq
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)

  real(DP) :: sxq (3, 48), work(3)
  ! list of vectors in the star of q
  real(DP), allocatable :: zstar(:,:,:)
  integer :: icart, jcart, ierr
  integer :: isym, iqstar
  logical :: ldiag_loc

  complex(DP), allocatable :: dyn_before(:,:), dyn_after(:,:), dyn_check(:,:)
  complex(DP) :: block(3,3)

  real(DP) :: norm_before, norm_after, norm_diff
  real(DP) :: norm_idem
  real(DP) :: norm_herm_before, norm_herm_after
  real(DP) :: block_norm
  !
  call start_clock('dynmatrix')
  ldiag_loc=ldiag.OR.(nat_todo_input > 0).OR.all_comp
  !
  !     set all noncomputed elements to zero
  !
  if (.not.lgamma_gamma) then
     imode0 = 0
     do irr = 1, nirr
        jmode0 = 0
        do jrr = 1, nirr
           if (.NOT.done_irr (irr).and..NOT.done_irr (jrr)) then
              do ipert = 1, npert (irr)
                 mu = imode0 + ipert
                 do jpert = 1, npert (jrr)
                    nu = jmode0 + jpert
                    dyn (mu, nu) = CMPLX(0.d0, 0.d0,kind=DP)
                 enddo
              enddo
           elseif (.NOT.done_irr (irr) .AND. done_irr (jrr) ) then
              do ipert = 1, npert (irr)
                 mu = imode0 + ipert
                 do jpert = 1, npert (jrr)
                    nu = jmode0 + jpert
                    dyn (mu, nu) = CONJG(dyn (nu, mu) )
                 enddo
              enddo
           endif
           jmode0 = jmode0 + npert (jrr)
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  else
     do irr = 1, nirr
        if (.NOT.comp_irr(irr)) then
           do nu=1,3*nat
              dyn(irr,nu)=(0.d0,0.d0)
           enddo
        endif
     enddo
  endif

    !
  !   Symmetrizes the dynamical matrix w.r.t. the small group of q
  !
  WRITE(stdout,*) ' '
  WRITE(stdout,*) '================ DEBUG dynmatrix_new ================'
  WRITE(stdout,*) 'DEBUG: current_iq = ', current_iq
  WRITE(stdout,*) 'DEBUG: nat        = ', nat
  WRITE(stdout,*) 'DEBUG: nmodes     = ', nmodes
  WRITE(stdout,*) 'DEBUG: nsym       = ', nsym
  WRITE(stdout,*) 'DEBUG: nsymq      = ', nsymq
  WRITE(stdout,*) 'DEBUG: lgamma     = ', lgamma
  WRITE(stdout,*) 'DEBUG: lgamma_gamma = ', lgamma_gamma
  WRITE(stdout,*) 'DEBUG: minus_q    = ', minus_q
  WRITE(stdout,*) 'DEBUG: irotmq     = ', irotmq
  WRITE(stdout,'(A,3ES18.8)') 'DEBUG: xq = ', xq(1), xq(2), xq(3)

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG: small-group operations used by QE'
  DO isym = 1, nsymq
     WRITE(stdout,*) 'DEBUG: ---- isym = ', isym, ' ----'
     WRITE(stdout,*) 'DEBUG: s(:,:,isym)'
     DO i = 1, 3
        WRITE(stdout,'(3I8)') s(i,1,isym), s(i,2,isym), s(i,3,isym)
     ENDDO

     WRITE(stdout,*) 'DEBUG: irt(isym,na)'
     DO na = 1, nat
        WRITE(stdout,'(A,I5,A,I5)') 'DEBUG: atom ', na, ' -> ', irt(isym,na)
     ENDDO

     WRITE(stdout,*) 'DEBUG: rtau(:,isym,na)'
     DO na = 1, nat
        WRITE(stdout,'(A,I5,A,3ES18.8)') 'DEBUG: atom ', na, ' rtau = ', &
             rtau(1,isym,na), rtau(2,isym,na), rtau(3,isym,na)
     ENDDO
  ENDDO

  IF (lgamma_gamma) THEN

     allocate(dyn_before(size(dyn,1), size(dyn,2)))
     allocate(dyn_after (size(dyn,1), size(dyn,2)))

     dyn_before(:,:) = dyn(:,:)

     norm_before = sqrt(sum(abs(dyn_before(:,:))**2))
     norm_herm_before = sqrt(sum(abs(dyn_before(:,:) - &
          conjg(transpose(dyn_before(:,:))))**2))

     CALL generate_dynamical_matrix (nat, nsym, s, invs, irt, at, bg, &
                       n_diff_sites, equiv_atoms, has_equivalent, dyn)
     IF (asr) CALL set_asr_c(nat,nasr,dyn)

     dyn_after(:,:) = dyn(:,:)

     norm_after = sqrt(sum(abs(dyn_after(:,:))**2))
     norm_diff  = sqrt(sum(abs(dyn_after(:,:) - dyn_before(:,:))**2))
     norm_herm_after = sqrt(sum(abs(dyn_after(:,:) - &
          conjg(transpose(dyn_after(:,:))))**2))

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG: gamma-gamma symmetrization'
     WRITE(stdout,*) 'DEBUG: ||dyn_before|| = ', norm_before
     WRITE(stdout,*) 'DEBUG: ||dyn_after || = ', norm_after
     WRITE(stdout,*) 'DEBUG: ||after-before|| = ', norm_diff
     IF (norm_before > 0.d0) THEN
        WRITE(stdout,*) 'DEBUG: relative change = ', norm_diff / norm_before
     ENDIF
     WRITE(stdout,*) 'DEBUG: hermiticity before = ', norm_herm_before
     WRITE(stdout,*) 'DEBUG: hermiticity after  = ', norm_herm_after

     deallocate(dyn_before)
     deallocate(dyn_after)

  ELSE

     allocate(dyn_before(size(dyn,1), size(dyn,2)))
     allocate(dyn_after (size(dyn,1), size(dyn,2)))
     allocate(dyn_check (size(dyn,1), size(dyn,2)))

     dyn_before(:,:) = dyn(:,:)

     norm_before = sqrt(sum(abs(dyn_before(:,:))**2))
     norm_herm_before = sqrt(sum(abs(dyn_before(:,:) - &
          conjg(transpose(dyn_before(:,:))))**2))

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG: before symdyn_munu_new'
     WRITE(stdout,*) 'DEBUG: ||dyn_before|| = ', norm_before
     WRITE(stdout,*) 'DEBUG: hermiticity before = ', norm_herm_before

     CALL symdyn_munu_new (dyn, u, xq, s, invs, rtau, irt, at, bg, &
          nsymq, nat, irotmq, minus_q)

     dyn_after(:,:) = dyn(:,:)

     norm_after = sqrt(sum(abs(dyn_after(:,:))**2))
     norm_diff  = sqrt(sum(abs(dyn_after(:,:) - dyn_before(:,:))**2))
     norm_herm_after = sqrt(sum(abs(dyn_after(:,:) - &
          conjg(transpose(dyn_after(:,:))))**2))

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG: after symdyn_munu_new'
     WRITE(stdout,*) 'DEBUG: ||dyn_after|| = ', norm_after
     WRITE(stdout,*) 'DEBUG: ||after-before|| = ', norm_diff
     IF (norm_before > 0.d0) THEN
        WRITE(stdout,*) 'DEBUG: relative symdyn change = ', norm_diff / norm_before
     ENDIF
     WRITE(stdout,*) 'DEBUG: hermiticity after = ', norm_herm_after

     !
     ! Idempotency check: symdyn(symdyn(D)) - symdyn(D)
     !
     dyn_check(:,:) = dyn_after(:,:)

     CALL symdyn_munu_new (dyn_check, u, xq, s, invs, rtau, irt, at, bg, &
          nsymq, nat, irotmq, minus_q)

     norm_idem = sqrt(sum(abs(dyn_check(:,:) - dyn_after(:,:))**2))

     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG: idempotency check'
     WRITE(stdout,*) 'DEBUG: ||symdyn(after)-after|| = ', norm_idem
     IF (norm_after > 0.d0) THEN
        WRITE(stdout,*) 'DEBUG: relative idempotency error = ', &
             norm_idem / norm_after
     ENDIF

     !
     ! 3x3 atom-pair blocks after symmetrization.
     ! dyn is a 2D matrix: dyn(mu,nu), where
     ! mu = 3*(na-1) + icart
     ! nu = 3*(nb-1) + jcart
     !
     WRITE(stdout,*) ' '
     WRITE(stdout,*) 'DEBUG: 3x3 atom-pair blocks after symdyn_munu_new'

     DO na = 1, nat
        DO nb = 1, nat

           DO icart = 1, 3
              DO jcart = 1, 3
                 mu = 3 * (na - 1) + icart
                 nu = 3 * (nb - 1) + jcart
                 block(icart,jcart) = dyn_after(mu,nu)
              ENDDO
           ENDDO

           block_norm = sqrt(sum(abs(block(:,:))**2))

           WRITE(stdout,*) ' '
           WRITE(stdout,'(A,I5,A,I5,A,ES18.8)') &
                'DEBUG: block na = ', na, ' nb = ', nb, &
                ' norm = ', block_norm

           WRITE(stdout,*) 'DEBUG: block rows: Re(1) Im(1) Re(2) Im(2) Re(3) Im(3)'
           DO icart = 1, 3
              WRITE(stdout,'(6ES18.8)') &
                   real(block(icart,1)), aimag(block(icart,1)), &
                   real(block(icart,2)), aimag(block(icart,2)), &
                   real(block(icart,3)), aimag(block(icart,3))
           ENDDO

        ENDDO
     ENDDO

     deallocate(dyn_before)
     deallocate(dyn_after)
     deallocate(dyn_check)

  ENDIF
  !
  !  if only one mode is computed write the dynamical matrix and stop
  !
  if (modenum .ne. 0) then
     WRITE( stdout, '(/,5x,"Dynamical matrix:")')
     do nu = 1, 3 * nat
        WRITE( stdout, '(5x,2i5,2f10.6)') modenum, nu, dyn (modenum, nu)
     enddo
     call stop_ph (.true.)
  endif

  IF ( .NOT. ldiag_loc ) THEN
     DO irr=0,nirr
        IF (.NOT.done_irr(irr)) THEN
           IF (.not.ldisp.AND..NOT.always_run) THEN
              WRITE(stdout, '(/,5x,"Stopping because representation", &
                                 & i5, " is not done")') irr
              CALL close_phq(.TRUE.)
              CALL stop_smoothly_ph(.TRUE.)
           ELSE
              WRITE(stdout, '(/5x,"Not diagonalizing because representation", &
                                 & i5, " is not done")') irr
           END IF
           RETURN
        ENDIF
     ENDDO
     ldiag_loc=.TRUE.
  ENDIF
    !
  !   Generates the star of q
  !
  call star_q1(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .TRUE., t_rev )

  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG: star of q'
  WRITE(stdout,*) 'DEBUG: nq  = ', nq
  WRITE(stdout,*) 'DEBUG: imq = ', imq

  DO iqstar = 1, nq
     WRITE(stdout,'(A,I5,A,3ES18.8)') 'DEBUG: sxq(:,', iqstar, ') = ', &
          sxq(1,iqstar), sxq(2,iqstar), sxq(3,iqstar)
  ENDDO

  WRITE(stdout,*) 'DEBUG: isq for all space-group operations'
  DO isym = 1, nsym
     WRITE(stdout,'(A,I5,A,I5)') 'DEBUG: isym = ', isym, ' isq = ', isq(isym)
  ENDDO
  !
  ! write on file information on the system
  !
  IF (xmldyn) THEN
     nqq=nq
     IF (imq==0) nqq=2*nq
     IF (lgamma.AND.done_epsil.AND.done_zeu) THEN
        CALL write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
             nqq, epsilon, zstareu, lraman, ramtns)
     ELSE
        CALL write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass, tau,ityp,m_loc,nqq)
     ENDIF
  ELSE
     CALL write_old_dyn_mat_head(iudyn)
  ENDIF
  !
  !   Rotates and writes on iudyn the dynamical matrices of the star of q
  !
  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'DEBUG: calling q2qstar_ph'
  WRITE(stdout,*) 'DEBUG: dyn passed to q2qstar_ph is the symmetrized dyn at xq'
  WRITE(stdout,'(A,3ES18.8)') 'DEBUG: xq passed to q2qstar_ph = ', &
       xq(1), xq(2), xq(3)

  call q2qstar_ph (dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
       nq, sxq, isq, imq, iudyn)

  WRITE(stdout,*) 'DEBUG: returned from q2qstar_ph'
  WRITE(stdout,*) '====================================================='
  WRITE(stdout,*) ' '

  !
  !   Writes (if the case) results for quantities involving electric field
  !
  if (epsil) call write_epsilon_and_zeu (zstareu, epsilon, nat, iudyn)
  IF (zue.AND..NOT.done_zue) THEN
     IF (lgamma_gamma) THEN
        ALLOCATE(zstar(3,3,nat))
        zstar(:,:,:) = 0.d0
        DO jcart = 1, 3
           DO mu = 1, 3 * nat
              na = (mu - 1) / 3 + 1
              icart = mu - 3 * (na - 1)
              zstar(jcart, icart, na) = zstarue0 (mu, jcart)
           ENDDO
           DO na=1,nat
              work(:)=0.0_DP
              DO icart=1,3
                 work(icart)=zstar(jcart,1,na)*at(1,icart)+ &
                             zstar(jcart,2,na)*at(2,icart)+ &
                             zstar(jcart,3,na)*at(3,icart)
              ENDDO
              zstar(jcart,:,na)=work(:)
           ENDDO
        ENDDO
        CALL generate_effective_charges_c ( nat, nsym, s, invs, irt, at, bg, &
           n_diff_sites, equiv_atoms, has_equivalent, asr, nasr, zv, ityp, &
           ntyp, atm, zstar )
        DO na=1,nat
           do icart=1,3
              zstarue(:,na,icart)=zstar(:,icart,na)
           ENDDO
        ENDDO
        done_zue=.TRUE.
        CALL summarize_zue()
        DEALLOCATE(zstar)
     ELSE
        CALL sym_and_write_zue
     ENDIF
  ELSEIF (lgamma) THEN
     IF (done_zue) CALL summarize_zue()
  ENDIF

  if (lraman) call write_ramtns (iudyn, ramtns)
  !
  !   Diagonalizes the dynamical matrix at q
  !
  IF (ldiag_loc) THEN
     call dyndia (xq, nmodes, nat, ntyp, ityp, amass, iudyn, dyn, w2)
     IF (search_sym) THEN
         CALL find_mode_sym_new (dyn, w2, tau, nat, nsymq, s, sr, irt, xq, &
              rtau, amass, ntyp, ityp, 1, lgamma_gamma, .FALSE., &
              num_rap_mode, ierr)
         CALL print_mode_sym(w2, num_rap_mode, lgamma)
     ENDIF
     IF (qplot) omega_disp(:,current_iq)=w2(:)
  END IF
!
! Here we save the dynamical matrix and the effective charges dP/du on
! the recover file. If a recover file with this very high recover code
! is found only the final result is rewritten on output.
!
  rec_code=30
  where_rec='dynmatrix.'
  CALL ph_writefile('status_ph',current_iq,0,ierr)

  call stop_clock('dynmatrix')
  return
end subroutine dynmatrix_new
