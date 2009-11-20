program test_nanCheck

implicit none

real :: a

a   = 5

write (*,*) 'a: ', a

a   = a / 0

write (*,*) 'a: ', a

if ( a*0 .ne. 0 ) write(*,*) 'a is inf'

a   = -5.0
a = sqrt(a)

write(*,*) 'a: ', a

if ( a .ne. a ) write (*,*) 'a is NaN'
end program test_nanCheck
