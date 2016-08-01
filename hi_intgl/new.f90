module hs_elem_mod
    type param_set
        integer :: npwr_k,npwr_g,npw
        integer,parameter :: nf = 8 
        real(rk),allocatable :: coef_g(:),coef_h(:),gpl(:),gwl(:)
        real(8) ::  beta
    end type

    type :: hs_elem
       private 
       integer,parameter ::  ndim = 3
       integer,parameter :: nedge=4
       integer tp
       real(rk),allocatable,dimension(:,:) :: nodes2D,nodes3D,nrmls3D

       type(point2D) :: src2D!,pt2D
       type(point3D) :: src3D!,pt3D

       type(param_set),dimension(:),allocatable :: params

       
   contains
       private
       procedure :: init
   end type

contains
    subroutine init(this,
        type(hs_elem) :: this

        elem_type = 8 !NCN(IELEM)!! to be changed
        this%beta = 3.
        allocate(this%nodes_m(3,8),this%nrmls_m(3,8))

    end subroutine

    ! @func : evaluate hyper singular integration)
    function evaluate_hsi(this)
        type(hs_elem) :: this
        type(edge2D) :: side,slice
        ! @var :[intg_pt] integration point
        type(point2D) :: intg_pt

        do i_edge = 1,num_edge ! ITERATE through each edge

                ks=ksb(i_edge)
                if(dabs(src_lcl(iabs(ks))-dble(ks)/dabs(dble(ks))).lt.tol) then
                    !print *,"src on current edge,Current edge skipped!"
                    goto 100
                end if
                !get the two end nodes coordinates for current edge

                ! assign slice head
                call slice%update_head(side%head)
                ! assign intg_pt part
                call side%setFixed(intg_pt,side%getFixed(slice%head)) 

                do num_converge=1,500 

                    step = get_step(side,slice)
                    if (step < -8.0d0) then
                        ! get out this side
                    end if
                   ! =====================

                    do igl = 1,iabs(ngl) !cutting each slice to ngl parts(gaussian points)

                        !pt_intg(unfixed)=seg_start(unfixed)+half_step*(1.d0+gpl(igl))
                        ! use the middle point
                        call side%setUnfixed(intg_pt,side%getUnFixed(slice%head)+0.5*step*(1.d0+gpl(igl)))
                        ! 1+gpl , is shifted gpl to start at least from 0
                        ! gpl span from -1 to 1, so intg-seg-size if half buffer point step
                        !rho_q=norm2(pt_intg-src_lcl)

                        rho_q=getDist(src2D,intg_pt)

                        !drdn_p=dabs(pt_intg(fixed)-src_lcl(fixed))/rho_q !sin(theta)

                        drdn_p=dabs(side%getFixed(intg_pt)-side%getFixed(src2D))/rho_q

                        ! calc radical integration along rho ----------------

                        !call compute_coeff_gh(num_dim,num_dim - 1,npw,elem_type,n_pwr_g,src_glb &
                                                !& ,src_lcl,pt_intg,coef_g,coef_h)
                        call this%compute_coeff_GH(intg_pt)


                        !call integrate_rho(1,ndim,nf,hi_beta,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint)
                        rint =  this%integrate_rho(1,params1,intg_pt)


                        !call integrate_rho(1,ndim,nf,hi_beta,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint)
                        !call integrate_rho(2,ndim,nf,2.d0,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint2)
                        !call integrate_rho(3,ndim,nf,1.d0,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint3)

                        ! end of radical calc along rho ------------------------

                        hiresult = hiresult + (dabs(half_step)*gwl(igl)*drdn_p/rho_q)*rint
                        str_result = str_result + (dabs(half_step)*gwl(igl)*drdn_p/rho_q)*rint2
                        wak_result = wak_result + (dabs(half_step)*gwl(igl)*drdn_p/rho_q)*rint3
                        ! equation (3-6-50)
                    end do ! igl =1,iabs(ngl)

                    !seg_start(unfixed)=seg_start(unfixed)+edge_direct*seg_step
                    ! move slice head by step size along side direction
                    call side%setUnfixed(slice%head,side%getUnfixed(slice%head)+side%dir*step)
                end do ! num_converge = 1,500

        100          end do ! i_edge = 1,num_edge

    end function

       function get_step(side,slice) result(res)
        type(edge2D) :: side,slice
        type(point2D) :: src2D

        real(rk) :: seg_step,d0,d1,d2,bool_expr
        real(rk) :: tmp
        ! distance from slice's tail to side edge's head
        d0 = this%dir*(side%getUnFixed(slice%tail)-side%getUnfixed(side%head))
        if(d0<1.0e-8) then
            !print *, "if seg_start located out of edge,exit num_converge loop"
            res = -9.0d0
            return
        end if

        d1 =           side%getFixed(src2D)  - side%getFixed(slice%head)
        d2 = side%dir*(side%getUnfixed(src2D)- side%getUnfixed(slice%head))

        ! compute step size====================

        if (d2 < 1.d-8) then !if src2D outside slice
            seg_step = fk*getDist(src2D,slice%head)
        else 
            tmp=d1**2*+d2**2 - (d1*fk)**2
            if(tmp < 0.d0) then
                seg_step=d2
            else
                seg_step=fk*(fk*d2-dsqrt(tmp))/(fk**2-1.)
            endif
        endif

        bool_expr = side%dir*(side%dir*seg_step+side%getUnfixed(slice%head))+1.d-8 

        if (bool_expr > d0) then
             if next step goes beyond the edge, use the end node
            seg_step = d0
        endif
        half_step = 0.5D0*side%dir*seg_step
        res = half_step
    end function

    function calc_wfa(beta,tolgp) result(res)
        real(rk),intent(in) :: beta,tolgp
        real(rk) :: res
        res=dsqrt(beta*2.d0/3.d0+0.4d0)*dlog(dabs(tolgp)/2.d0)  
    end function

    function calc_fk(wfa,ngl) result(res)
        integer,intent(in) :: ngl
        real(rk),intent(in) :: wfa
        res=3.d0/8.d0*(-10.d0*dble(iabs(ngl))/wfa-1.d0)**(4.d0/3.d0)
    end function



end module

