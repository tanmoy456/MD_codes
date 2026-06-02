program q1
    implicit none
    integer     :: i,j,temp,c
    integer ,dimension(10)  :: x
  !  print*, "Enter the numbers sequentially:"
    x = (/0,-20,27,27,13,15,-20,0,22,-20/)
    !x = (/9,0,5,7,2,88,3,4,6,1/)
    !do i =1,10
     !   read*,x(i)
   ! end do
   print*,"Initial x:",x
    do j=1,10
        do i=1,9
            if (x(i) > x(i+1)) then 
                temp = x(i)
                x(i) = x(i+1)
                x(i+1) = temp
            end if
      ! print*,x
        end do
      ! print*,x
    end do
    
    print*,"Final x:",x
    print*,"The largest number is ",x(10)
    print*,"The smallest number is",x(1)

end program q1
