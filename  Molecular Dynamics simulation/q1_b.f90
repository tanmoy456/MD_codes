program frequency 
    implicit none
    integer     :: i,j,temp,c,d
    integer ,dimension(10)  :: x
  !  print*, "Enter the numbers sequentially:"
    x = (/0,-20,27,27,13,15,-20,0,22,-20/)
    !x = (/9,0,5,7,2,88,3,4,6,1/)
    !do i =1,10
     !   read*,x(i)
   ! end do
    print*,"Initial x:",x
    print*,"Repeating Numbers ","Crossponding positions"
    do j =1,10
        c = 1
       ! print*,x(j),c
        do i=j+1,10
            if (x(j) == x(i)) then
               ! c = c+1
               print*,x(i),j,i
            end if
        end do 
    end do
end program frequency
