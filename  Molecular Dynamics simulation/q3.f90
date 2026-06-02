program q3
    implicit none
    integer     :: i,j,k
    integer,dimension(3,4) :: A
    integer,dimension(4,2) :: B
    integer,dimension(3,2) :: C
    integer,dimension(3,2) :: D
    !Creating Matrix A
    do i =1,3 !i =1 to 3 is row
        do j=1,4 !j=1 to 4 is column
            read*,A(i,j)
        end do
    end do
    print*,"The Matrix A :"
    do i=1,3 !3 row
        write(*,*)A(i,:)
    end do
    !Creating Matrix B
    do i =1,4 !i =1 to 4 is row
        do j=1,2 !j=1 to 2 is column
            read*,B(i,j)
        end do
    end do
    print*,"The Matrix B :"
    do i=1,4 !4 row
        write(*,*)B(i,:)
    end do
    !Matrix Multiplication
    !using matmul
    print*,"Matmul D:",matmul(A,B)
    !using do loop
    do i=1,3
        do j=1,2
            C(i,j) = 0
            do k=1,4
                C(i,j) = C(i,j) + A(i,k)*B(k,j)
            end do
        end do
    end do
    print*,"The Matrix C using do loop:"
    do i=1,3
        write(*,*)C(i,:)
    end do

end program 
