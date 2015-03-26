!#############################################################################################################################################################################################################################

subroutine PearsnR4 (x,y,n,r)

implicit none
integer n
real(4) prob,r,z,x(n),y(n),TINY
parameter (tiny=1.e-20)
integer j
real(4) ax,ay,df,sxx,sxy,syy,t,xt,betai,yt

ax=0.0
ay=0.0
DO j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
END DO
ax=ax/n
ay=ay/n
sxx=0.
syy=0.
sxy=0.
DO j=1,n
        xt=x(j)-ax
        yt=y(j)-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
END DO
r=sxy/(SQRT(sxx*syy)+TINY)
z=0.5*LOG(((1.+r)+TINY)/((1.-r)+TINY))
df=n-2
t=r*SQRT(df/(((1.-r)+TINY)*((1.+r)+TINY)))
!prob=betai(0.5*df,0.5,df/(df+t**2))
!prob=erfcc(ABS(z*SQRT(n-1.))/1.4142136)
prob=0
return

end subroutine PearsnR4

!###########################################################################################################################################################
