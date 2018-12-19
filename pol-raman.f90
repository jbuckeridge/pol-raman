program pol_raman
!
! This program reads in the Raman tensors from the file
! 'Raman-Tensors.yaml', which is output from Jonathan Skelton's
! python code Phonopy-Spectroscopy used to calculate Raman 
! intensities (among other stuff!). For each Raman
! active mode, the intensity vs angle theta (t) is calculated for
! incident and reflected light polarised parallel to a user-
! defined plane. That is, we model a back-scattering Raman
! experiment where the plane defined by the user is the reflection
! plane. The intensities are calculated for the polarisation
! direction varied through 360 degrees
!
! I proportional | [a,b,c] R(omega) [a] |^2
!                |                  [b] | 
!                |                  [c] |
!
! where [a,b,c] is a vector parallel to the plane (hkl)
! i.e. perpendicular to [h,k,l]. [a,b,c] is normalised.
! We can then work out values for a,b,c and set an angle
! that is varied (see my logbook 18-12-18 for the algebra).
!
! If h^2 + l^2 = 0, then we set b=0 and a=cos(t), c=sin(t)
! etc
!
! Otherwise, c = cos(t) h / sqrt(h^2 + l^2);
! 
! b = (-2hklcos(t) - sqrt( 4h^2k^2l^2cos^2(t)/(l^2+h^2) - 
!            4(k^2+h^2)h^2(cos^2(t)-1)) ) / 2(k^2+h^2);
!
! a = sqrt(1 - b^2 - c^2)
!
! The sign of a must be chosen to make sure the vector is 
! parallel to (hkl). Note that the minus sign has been chosen
! preceding the sqrt in the denominator of b.
!
! J. Buckeridge Dec. 2018
!
  implicit none
!
  integer, parameter :: nmax = 500
  integer, parameter :: maxang = 360
  real*8, parameter :: small = 1.d-12
  integer i, j, k, band(nmax)
  integer stat, nbands
  integer hi, ki, li
  real*8 pi, ram(nmax,3,3), freq(nmax), intens, theta
  real*8 hval, kval, lval, aval(nmax), bval(nmax), cval(nmax)
  real*8 test1, test2, test3, test4
  character*11 tmpstr1
  character*12 raman1, raman2
  character*10 tmpstr2
  character*1 tmpstr3
  character*2 tmpstr4
  character*3 tmpstr5
  character*20 filename
!
  pi = acos(-1.0)
!
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') 
  write(*,'(a)') "    PPP          L    RRR   " 
  write(*,'(a)') "    P   P        L    R   R  " 
  write(*,'(a)') "    PP  P        L    RR  R  " 
  write(*,'(a)') "    PPPP  OOOOO  LL   RRRR   AAAAA  MM   MM  AAAAA  NNNNN " 
  write(*,'(a)') "    PP    O   O  LL   RR  R     AA  M M M M     AA  N   N " 
  write(*,'(a)') "    PP    OO  O  LL   RR  R  AA  A  MM M  M  AA  A  NN  N " 
  write(*,'(a)') "    PP    OOOOO  LLLL RR  R  AAAAA  MM    M  AAAAA  NN  N " 
  write(*,'(a)') 
  write(*,'(a)') 
  write(*,'(a)') "j.buckeridge@ucl.ac.uk 2018"
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') 
!
  open(unit=11,file="Raman-Tensors.yaml",status="old",iostat=stat)
!
  if(stat /= 0) then
     write(*,'(a)') "ERROR: cannot find Raman-Tensors.yaml file!!"
     stop
  else
     write(*,'(a)') "Found Raman-Tensors.yaml. Reading in tensors ..."
  endif
!
! Skip first four lines of file
!
  do i=1,4
     read(11,*)
  enddo
!
! Use a loop to read in each active mode and its Raman tensor
!
  i=0
  do
     i = i+1
     read(11,*,end=101)
     read(11,*) tmpstr1, band(i)
     read(11,*) tmpstr2, freq(i)
     read(11,*)
     do j=1,3
        read(11,*) tmpstr3, tmpstr3, raman1, raman2, ram(i,j,3)
        if(raman1(1:1) == "-") then
           read(raman1(1:11),*) ram(i,j,1)
        else
           read(raman1(1:10),*) ram(i,j,1)
        endif
        if(raman2(1:1) == "-") then
           read(raman2(1:11),*) ram(i,j,2)
        else
           read(raman2(1:10),*) ram(i,j,2)
        endif
     enddo
  enddo
101 nbands = i-1
  close(unit=11)
  write(*,'(a5,x,i3,x,a21)') "Found", nbands, "Raman-active modes..."
  write(*,*)
!
! Now ask for miller indices of plane of interest
!
  write(*,'(a)') "Please enter the miller indices of the plane of interest &
       &(as 3 integers, e.g. 1 0 0 or -1 -1 2 etc):"
  read(*,*,iostat=stat) hi, ki, li
  if(stat /= 0) then
     write(*,'(a)') "ERROR reading in miller indices!"
     stop
  endif
  write(*,*)
  hval=real(hi)
  kval=real(ki)
  lval=real(li)
!
! Calculate the a,b,c parameters for t=0 first and report the polarisation
! diretion that corresponds to t=0 and t=90
!
  test1 = hval**2 + lval**2
  test2 = hval**2 + kval**2
  test3 = kval**2 + lval**2
  if( abs(test1) < small ) then
     bval = 0.d0
     do i=1,maxang+1
        theta = real(i-1) * pi / 180.d0
        aval(i) = cos(theta)
        cval(i) = -sin(theta)
     enddo
     write(*,'(a)') "Polarisation direction at theta = 0 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(1), bval(1), cval(1), "]"
     write(*,'(a)') "Polarisation direction at theta = pi/2 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(91), bval(91), cval(91), "]"
  elseif( abs(test2) < small ) then
     cval = 0.d0
     do i=1,maxang+1
        theta = real(i-1) * pi / 180.d0
        aval(i) = cos(theta)
        bval(i) = sin(theta)
     enddo
     write(*,'(a)') "Polarisation direction at theta = 0 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(1), bval(1), cval(1), "]"
     write(*,'(a)') "Polarisation direction at theta = pi/2 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(91), bval(91), cval(91), "]"
  elseif( abs(test3) < small ) then
     aval = 0.d0
     do i=1,maxang+1
        theta = real(i-1) * pi / 180.d0
        bval(i) = cos(theta)
        cval(i) = sin(theta)
     enddo
     write(*,'(a)') "Polarisation direction at theta = 0 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(1), bval(1), cval(1), "]"
     write(*,'(a)') "Polarisation direction at theta = pi/2 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(91), bval(91), cval(91), "]"
  else
     do i=1,maxang+1
        theta = real(i-1) * pi / 180.d0
        cval(i) = hval * cos(theta) / sqrt(hval**2 + lval**2)
        bval(i) = 4.d0 * hval**2 * kval**2 * lval**2 * cos(theta)**2 / (hval**2 + lval**2)
        bval(i) = bval(i) - 4.d0 * (hval**2 + kval**2) * hval**2 * (cos(theta)**2 - 1.d0)
        if(abs(bval(i)) < small) bval(i) = 0.d0
        bval(i) = sqrt(bval(i))
        bval(i) = -2.d0 * hval * kval * lval * cos(theta) / sqrt(hval**2 + lval**2) - bval(i)
        bval(i) = bval(i) / (2.d0 * (kval**2 + hval**2))
        aval(i) = 1.d0 - bval(i)**2 - cval(i)**2
        if(abs(aval(i)) < small) aval(i) = 0.d0
        aval(i) = sqrt(aval(i))
!
! Need to check that the aval is correct - it may need a minus sign
!
        test4 = hval * aval(i) + kval * bval(i) + lval * cval(i)
        if( abs(test4) > small ) then
           aval(i) = -aval(i)
        endif
        test4 = hval * aval(i) + kval * bval(i) + lval * cval(i)
        if( abs(test4) > small ) then
           write(*,'(a)') "Hmmm not sure what's going wrong here :("
           stop
        endif
     enddo
     write(*,'(a)') "Polarisation direction at theta = 0 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(1), bval(1), cval(1), "]"
     write(*,'(a)') "Polarisation direction at theta = pi/2 :"
     write(*,'(a1,2x,3(f13.6,2x),a1)') "[", aval(91), bval(91), cval(91), "]"
  endif
!
! Now calculate intensities for reflections from user-defined plane as a function
! of angle
!
  write(*,*)
  write(*,'(a)') "Writing intensity vs angle output to files titled 'modeXX.dat', where XX is the mode index..."
  do i=1,nbands
     if( band(i) < 10 ) then
        write(tmpstr3,'(i1)') band(i)
        tmpstr3 = trim(tmpstr3)
        filename = "mode"//tmpstr3//".dat"
     elseif (band(i) < 100 ) then
        write(tmpstr4,'(i2)') band(i)
        tmpstr4 = trim(tmpstr4)
        filename = "mode"//tmpstr4//".dat"
     else
        write(tmpstr5,'(i3)') band(i)
        tmpstr5 = trim(tmpstr5)
        filename = "mode"//tmpstr5//".dat"
     endif
     open(unit=11,file=filename,status="replace",iostat=stat)
     write(11,*) "# band", band(i), "freq", freq(i)
     write(11,*) "# theta (deg), intensity"
     do j=1,maxang+1
        intens = aval(j)**2 * ram(i,1,1) + aval(j) * bval(j) * ram(i,2,1) + aval(j) * cval(j) * ram(i,3,1) &
             + aval(j) * bval(j) * ram(i,1,2) + bval(j)**2 * ram(i,2,2) + bval(j) * cval(j) * ram(i,3,2) &
             + aval(j) * cval(j) * ram(i,1,3) + bval(j) * cval(j) * ram(i,2,3) + cval(j)**2 * ram(i,3,3)
        intens = intens * intens
        write(11,'(f12.5,3x,e20.13)') real(j-1), intens
     enddo
     close(unit=11)
  enddo
  write(*,'(a)') "...done"
!
end program pol_raman
