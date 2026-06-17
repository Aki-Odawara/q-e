!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine symdynph_gq_new( xq, phi, s, invs, rtau, irt, nsymq, &
                            nat, irotmq, minus_q )
  !-----------------------------------------------------------------------
  !! This routine receives as input an unsymmetrized dynamical
  !! matrix expressed on the crystal axes and imposes the symmetry
  !! of the small group of q. Furthermore it imposes also the symmetry
  !! q -> -q+G if present.  
  !! February 2020: Update (A. Urru) to include the symmetry operations 
  !! that require the time reversal operator (meaning that TS is a 
  !! symmetry of the crystal). For more information please see: 
  !! Phys. Rev. B 100, 045115 (2019).
  !
  USE kinds, only : DP
  USE constants, ONLY: tpi
  USE symm_base, ONLY : t_rev
  !
  implicit none
  !
  integer :: nat
  !! input: the number of atoms
  integer :: s(3,3,48)
  !! input: the symmetry matrices
  integer :: irt(48,nat)
  !! input: the rotated of each vector
  integer :: invs(48)
  !! input: the inverse of each matrix
  integer :: nsymq
  !! input: the order of the small group
  integer :: irotmq
  !! input: the rotation sending q ->-q+G
  real(DP) :: xq(3)
  !! input: the q point
  real(DP) :: rtau(3,48,nat)
  !! input: the R associated at each t
  logical :: minus_q
  !! input: true if a symmetry q->-q+G
  complex(DP) :: phi(3,3,nat,nat)
  !! inp/out: the matrix to symmetrize
  !
  ! ... local variables
  !
  integer :: isymq, sna, snb, irot, na, nb, ipol, jpol, lpol, kpol, &
           iflb (nat, nat)
  integer :: max_isym, max_na, max_nb
  ! counters, indices, work space

  real(DP) :: arg
  real(DP) :: norm_phi, norm_res, rel_res
  real(DP) :: norm_pair, norm_res_pair, rel_pair
  real(DP) :: max_rel_res
  ! the argument of the phase

  complex(DP) :: phip (3, 3, nat, nat), work (3, 3), fase, faseq (48)
  complex(DP), allocatable :: phig(:,:,:,:)
  ! work space, phase factors
  !
  !    We start by imposing hermiticity
  !
  write(*, *) "### symdynph_gq is entered"
  write(*, *) "Initial Dynamical Matrix="
  write(*, '(6F10.5)') phi
  do na = 1, nat
     do nb = 1, nat
        do ipol = 1, 3
           do jpol = 1, 3
              phi (ipol, jpol, na, nb) = 0.5d0 * (phi (ipol, jpol, na, nb) &
                   + CONJG(phi (jpol, ipol, nb, na) ) )
              phi (jpol, ipol, nb, na) = CONJG(phi (ipol, jpol, na, nb) )
           enddo
        enddo
     enddo
  enddo
  write(*, *) "Imposing hermiticity have done !"
  write(*, *) "Dynamical Matrix after imposing hermiticity ="
  write(*, '(6F10.5)') phi
    !
  ! DEBUG: residual before small-group symmetrization
  ! This checks T_g[phi] - phi using the same convention as this routine.
  !
  if (nsymq > 1) then

     allocate(phig(3,3,nat,nat))

     norm_phi = sqrt(sum(abs(phi(:,:,:,:))**2))

     write(*, *) " "
     write(*, *) "### DEBUG symdynph_gq: residual BEFORE small-group symmetrization"
     write(*, *) "DEBUG: ||phi|| after hermiticity = ", norm_phi

     max_rel_res = 0.d0
     max_isym = 0
     max_na = 0
     max_nb = 0

     do isymq = 1, nsymq

        irot = isymq
        phig(:,:,:,:) = (0.d0, 0.d0)

        do na = 1, nat
           do nb = 1, nat

              sna = irt(irot, na)
              snb = irt(irot, nb)

              arg = 0.d0
              do ipol = 1, 3
                 arg = arg + xq(ipol) * &
                      (rtau(ipol, irot, na) - rtau(ipol, irot, nb))
              enddo
              arg = arg * tpi
              fase = CMPLX(cos(arg), sin(arg), kind=DP)

              do ipol = 1, 3
                 do jpol = 1, 3
                    do kpol = 1, 3
                       do lpol = 1, 3
                          if (t_rev(isymq) == 1) then
                             phig(ipol,jpol,na,nb) = phig(ipol,jpol,na,nb) + &
                                  s(ipol,kpol,irot) * s(jpol,lpol,irot) * &
                                  CONJG(phi(kpol,lpol,sna,snb) * fase)
                          else
                             phig(ipol,jpol,na,nb) = phig(ipol,jpol,na,nb) + &
                                  s(ipol,kpol,irot) * s(jpol,lpol,irot) * &
                                  phi(kpol,lpol,sna,snb) * fase
                          endif
                       enddo
                    enddo
                 enddo
              enddo

           enddo
        enddo

        norm_res = sqrt(sum(abs(phig(:,:,:,:) - phi(:,:,:,:))**2))
        rel_res = 0.d0
        if (norm_phi > 0.d0) rel_res = norm_res / norm_phi

        write(*,*) "DEBUG BEFORE: isymq = ", isymq, &
                   " ||Tg(phi)-phi|| = ", norm_res, &
                   " relative = ", rel_res

        if (rel_res > max_rel_res) then
           max_rel_res = rel_res
           max_isym = isymq
        endif

        do na = 1, nat
           do nb = 1, nat
              norm_pair = sqrt(sum(abs(phi(:,:,na,nb))**2))
              norm_res_pair = sqrt(sum(abs(phig(:,:,na,nb) - phi(:,:,na,nb))**2))
              rel_pair = 0.d0
              if (norm_pair > 0.d0) rel_pair = norm_res_pair / norm_pair

              write(*,*) "DEBUG BEFORE block: isymq=", isymq, &
                         " na=", na, " nb=", nb, &
                         " ||block||=", norm_pair, &
                         " ||res||=", norm_res_pair, &
                         " rel=", rel_pair
           enddo
        enddo

     enddo

     write(*,*) "DEBUG BEFORE: max relative residual = ", max_rel_res, &
                " at isymq = ", max_isym

     deallocate(phig)

  endif

  !
  !    If no other symmetry is present we quit here
  !
  if ( (nsymq == 1) .and. (.not.minus_q) ) then
      write(*, *) "no symmetry"
      return
  endif
  !
  !    Then we impose the symmetry q -> -q+G if present
  !
  if (minus_q) then
     do na = 1, nat
        do nb = 1, nat
           do ipol = 1, 3
              do jpol = 1, 3
                 work(:,:) = (0.d0, 0.d0)
                 sna = irt (irotmq, na)
                 snb = irt (irotmq, nb)
                 arg = 0.d0
                 do kpol = 1, 3
                    arg = arg + (xq (kpol) * (rtau (kpol, irotmq, na) - &
                                              rtau (kpol, irotmq, nb) ) )
                 enddo
                 arg = arg * tpi
                 fase = CMPLX(cos (arg), sin (arg) ,kind=DP)
                 do kpol = 1, 3
                    do lpol = 1, 3
                       work (ipol, jpol) = work (ipol, jpol) + &
                            s (ipol, kpol, irotmq) * s (jpol, lpol, irotmq) &
                            * phi (kpol, lpol, sna, snb) * fase
                    enddo
                 enddo
                 phip (ipol, jpol, na, nb) = (phi (ipol, jpol, na, nb) + &
                      CONJG( work (ipol, jpol) ) ) * 0.5d0
              enddo
           enddo
        enddo
     enddo
     phi = phip
     write(*, *) "Impose the time-reversal symmetry!"
     write(*, *) "Dynamical Matrix after imposing time-reversal symmetry ="
     write(*, '(6F10.5)') phi
  endif
  !
  !    Here we symmetrize with respect to the small group of q
  !
  if (nsymq == 1) return

  iflb (:, :) = 0
  do na = 1, nat
     do nb = 1, nat
        if (iflb (na, nb) == 0) then
           work(:,:) = (0.d0, 0.d0)
           do isymq = 1, nsymq
              irot = isymq
              sna = irt (irot, na)
              snb = irt (irot, nb)
              write(*, *) "rotational index:", irot, "before:", na, "after:", sna
              write(*, *) "rotational index:", irot, "before:", nb, "after:", snb
              arg = 0.d0
              do ipol = 1, 3
                 arg = arg + (xq (ipol) * (rtau(ipol, irot, na) - &
                                           rtau(ipol, irot, nb) ) )
                 write(*, *) "cell return vector of na=", rtau(ipol, irot, na)
                 write(*, *) "cell return vector of nb=", rtau(ipol, irot, nb)
                 write(*, *) "difference btw unit cell=", rtau(ipol, irot, na) - &
                                                          rtau(ipol, irot, nb)
                 write(*, *) "component of wavevec.=", xq(ipol)
              enddo
              write(*, *) "arg=", arg
              arg = arg * tpi
              write(*, *) "atoms number=", na, "and", nb
              write(*, *) "argument=", arg
              faseq (isymq) = CMPLX(cos (arg), sin (arg) ,kind=DP)
              write(*, *) "fase factor=", faseq(isymq)
              do ipol = 1, 3
                 do jpol = 1, 3
                    do kpol = 1, 3
                       do lpol = 1, 3
                          IF (t_rev(isymq)==1) THEN
                             work (ipol, jpol) = work (ipol, jpol) + &
                                  s (ipol, kpol, irot) * s (jpol, lpol, irot) &
                           * CONJG(phi (kpol, lpol, sna, snb) * faseq (isymq))
                          ELSE
                             work (ipol, jpol) = work (ipol, jpol) + &
                                  s (ipol, kpol, irot) * s (jpol, lpol, irot) &
                                 * phi (kpol, lpol, sna, snb) * faseq (isymq)
                          ENDIF
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           do isymq = 1, nsymq
              irot = isymq
              sna = irt (irot, na)
              snb = irt (irot, nb)
              do ipol = 1, 3
                 do jpol = 1, 3
                    phi (ipol, jpol, sna, snb) = (0.d0, 0.d0)
                    do kpol = 1, 3
                       do lpol = 1, 3
                          IF (t_rev(isymq)==1) THEN
                             phi(ipol,jpol,sna,snb)=phi(ipol,jpol,sna,snb) &
                             + s(ipol,kpol,invs(irot))*s(jpol,lpol,invs(irot))&
                               * CONJG(work (kpol, lpol)*faseq (isymq))
                          ELSE
                             phi(ipol,jpol,sna,snb)=phi(ipol,jpol,sna,snb) &
                             + s(ipol,kpol,invs(irot))*s(jpol,lpol,invs(irot))&
                               * work (kpol, lpol) * CONJG(faseq (isymq) )
                          ENDIF
                       enddo
                    enddo
                 enddo
              enddo
              iflb (sna, snb) = 1
              write(*, *) "inverse=", invs(irot)
           enddo
        endif
     enddo
  enddo
    phi (:, :, :, :) = phi (:, :, :, :) / DBLE(nsymq)

  write(*, *) "RT D R with respect to the small group of q have done !"
  write(*, *) "Dynamical Matrix after symmetrizing ="
  write(*, '(6F10.5)') phi

  !
  ! DEBUG: residual after small-group symmetrization
  ! This checks whether the symmetrized phi is invariant under each
  ! operation of the small group, using exactly the same convention
  ! as the symmetrization above.
  !
  allocate(phig(3,3,nat,nat))

  norm_phi = sqrt(sum(abs(phi(:,:,:,:))**2))

  write(*, *) " "
  write(*, *) "### DEBUG symdynph_gq: residual AFTER small-group symmetrization"
  write(*, *) "DEBUG: ||phi_sym|| = ", norm_phi

  max_rel_res = 0.d0
  max_isym = 0
  max_na = 0
  max_nb = 0

  do isymq = 1, nsymq

     irot = isymq
     phig(:,:,:,:) = (0.d0, 0.d0)

     do na = 1, nat
        do nb = 1, nat

           sna = irt(irot, na)
           snb = irt(irot, nb)

           arg = 0.d0
           do ipol = 1, 3
              arg = arg + xq(ipol) * &
                   (rtau(ipol, irot, na) - rtau(ipol, irot, nb))
           enddo
           arg = arg
           fase = CMPLX(cos(arg), sin(arg), kind=DP)

           do ipol = 1, 3
              do jpol = 1, 3
                 do kpol = 1, 3
                    do lpol = 1, 3
                       if (t_rev(isymq) == 1) then
                          phig(ipol,jpol,na,nb) = phig(ipol,jpol,na,nb) + &
                               s(ipol,kpol,irot) * s(jpol,lpol,irot) * &
                               CONJG(phi(kpol,lpol,sna,snb) * fase)
                       else
                          phig(ipol,jpol,na,nb) = phig(ipol,jpol,na,nb) + &
                               s(ipol,kpol,irot) * s(jpol,lpol,irot) * &
                               phi(kpol,lpol,sna,snb) * fase
                       endif
                    enddo
                 enddo
              enddo
           enddo

        enddo
     enddo

     norm_res = sqrt(sum(abs(phig(:,:,:,:) - phi(:,:,:,:))**2))
     rel_res = 0.d0
     if (norm_phi > 0.d0) rel_res = norm_res / norm_phi

     write(*,*) "DEBUG AFTER: isymq = ", isymq, &
                " ||Tg(phi_sym)-phi_sym|| = ", norm_res, &
                " relative = ", rel_res

     if (rel_res > max_rel_res) then
        max_rel_res = rel_res
        max_isym = isymq
     endif

     do na = 1, nat
        do nb = 1, nat
           norm_pair = sqrt(sum(abs(phi(:,:,na,nb))**2))
           norm_res_pair = sqrt(sum(abs(phig(:,:,na,nb) - phi(:,:,na,nb))**2))
           rel_pair = 0.d0
           if (norm_pair > 0.d0) rel_pair = norm_res_pair / norm_pair

           write(*,*) "DEBUG AFTER block: isymq=", isymq, &
                      " na=", na, " nb=", nb, &
                      " ||block||=", norm_pair, &
                      " ||res||=", norm_res_pair, &
                      " rel=", rel_pair
        enddo
     enddo

  enddo

  write(*,*) "DEBUG AFTER: max relative residual = ", max_rel_res, &
             " at isymq = ", max_isym

  deallocate(phig)

  return
end subroutine symdynph_gq_new
