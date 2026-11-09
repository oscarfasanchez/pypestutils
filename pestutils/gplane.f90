! GPLANE projection and inverse-distance factors, following PLPROC 5 section 2.5.
! Stateless: callers own the factors and can reuse them for new property values.
integer(c_int) function calc_gplane_factors(ns,su,sv,nt,tx,ty,tz,origin,ex,ey,normal, &
    nb,bu,bv,active,width,decay,power,ratio,bearing,tol,weights,influence) &
    bind(C,name="calc_gplane_factors")
    use iso_c_binding, only: c_int,c_double
    use utilities, only: amessage
    implicit none
    integer(c_int), intent(in) :: ns,nt,nb,active(nt)
    real(c_double), intent(in) :: su(ns),sv(ns),tx(nt),ty(nt),tz(nt)
    real(c_double), intent(in) :: origin(3),ex(3),ey(3),normal(3),bu(nb),bv(nb)
    real(c_double), intent(in) :: width,decay,power,ratio,bearing,tol
    real(c_double), intent(out) :: weights(ns,nt),influence(nt)
    real(c_double) :: delta(3),u,v,d,qx,qy,dx,dy,t,edge2,best2,cx,cy
    real(c_double) :: cs,sn,du,dv,d2,mind2,total,f1,f2,along,across,dscale,rtol
    integer :: i,j,k,l,nearest
    logical :: inside

    calc_gplane_factors=1
    if(ns<1.or.nt<1.or.(nb/=0.and.nb<3).or.width<0.or.decay<=0.or. &
       power<=0.or.ratio<=0.or.tol<=0) then
        amessage='Invalid dimensions or interpolation settings in calc_gplane_factors.'
        return
    endif
    weights=0.0d0
    influence=0.0d0
    cs=cos(bearing*acos(-1.0d0)/180.0d0)
    sn=sin(bearing*acos(-1.0d0)/180.0d0)
    do i=1,nt
        if(active(i)==0) cycle
        delta=[tx(i),ty(i),tz(i)]-origin
        u=dot_product(delta,ex)
        v=dot_product(delta,ey)
        d=abs(dot_product(delta,normal))
        f1=decay_weight(d,width,decay)
        if(f1==0.0d0) cycle
        qx=u
        qy=v
        f2=1.0d0
        if(nb>0) then
            inside=.false.
            best2=huge(1.0d0)
            do j=1,nb
                k=mod(j,nb)+1
                dx=bu(k)-bu(j)
                dy=bv(k)-bv(j)
                edge2=dx*dx+dy*dy
                t=max(0.0d0,min(1.0d0,((u-bu(j))*dx+(v-bv(j))*dy)/edge2))
                cx=bu(j)+t*dx
                cy=bv(j)+t*dy
                d2=(u-cx)**2+(v-cy)**2
                if(d2<best2) then
                    best2=d2
                    qx=cx
                    qy=cy
                endif
                if((bv(j)>v).neqv.(bv(k)>v)) then
                    if(u < bu(j)+(v-bv(j))*dx/dy) inside=.not.inside
                endif
            enddo
            if(inside.or.best2<=tol*tol) then
                qx=u
                qy=v
            else
                f2=decay_weight(sqrt(best2),width,decay)
            endif
        endif
        influence(i)=f1*f2
        if(influence(i)==0.0d0) cycle
        mind2=huge(1.0d0)
        nearest=1
        ! Common distance scaling cancels in normalized IDW, and avoids d**4 overflow.
        dscale=max(maxval(abs(su-qx)),maxval(abs(sv-qy)),tol)
        rtol=tol/dscale
        if(ratio>=1.0d0) rtol=rtol/ratio
        do j=1,ns
            du=(su(j)-qx)/dscale
            dv=(sv(j)-qy)/dscale
            along=du*sn+dv*cs
            across=du*cs-dv*sn
            if(ratio>=1.0d0) then
                along=along/ratio
            else
                across=across*ratio
            endif
            ! PLPROC 5 GPLANE IDW uses sqrt(du**4 + dv**4) as distance.
            ! Keep this executable-observed convention for numeric parity.
            d2=along**4+across**4
            weights(j,i)=d2
            if(d2<mind2) then
                mind2=d2
                nearest=j
            endif
        enddo
        if(mind2<=rtol**4) then
            weights(:,i)=0.0d0
            weights(nearest,i)=1.0d0
        else
            total=0.0d0
            do l=1,ns
                weights(l,i)=(mind2/weights(l,i))**(power/2.0d0)
                total=total+weights(l,i)
            enddo
            weights(:,i)=weights(:,i)/total
        endif
    enddo
    calc_gplane_factors=0
contains
    real(c_double) function decay_weight(distance,w,a) result(f)
        real(c_double), intent(in) :: distance,w,a
        real(c_double) :: excess
        excess=max(0.0d0,distance-w/2.0d0)
        f=0.0d0
        ! Compare before dividing so tiny positive a never trips the FPE trap.
        if(excess/27.0d0<a) f=exp(-(excess/a)**2)
    end function decay_weight
end function calc_gplane_factors
