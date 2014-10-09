      program freqanal
!
!-----Created by Bin Liu.
!-----Created on August 6 2010.
!-----Get Vib Freq info from VASP OUTCAR
!-----Output file can be open with Molden. 
!
      implicit none
!
      integer i, j, k, na, namax, ntype, ntypemax, atnumtemp, nfreq 
      integer nfreqmax
      parameter (namax = 100)
      parameter (ntypemax = 10) 
      integer iatm(namax), numatom(ntypemax), atnum(100)
      real  lattice, zero
      real  c(namax,3),cell(3,3),ctemp(namax,3),disp(namax,3)
      character*1 fix(namax,3), qX, qY, qZ, qdX, qdY, qdZ 
      character*2 atname(namax), typename(ntypemax)
      character*6 qformat
      character*11 qfreq(3*namax)
      character*20 mode
      character*25 title
      character*30 fnameout
      character*80 qline
!   
!-----Initialization.
!      
      zero = 0.0
      do i = 1, 3
         do j = 1, 3
            cell(i,j) = 0.0
         enddo
      enddo      
      do j = 1, namax
         do k = 1, 3
            disp(j,k) = 0.0
         enddo
      enddo
      write(*,*) 'Specify the output file:'
      read (*,*) fnameout
!
!-----Read CONTCAR.      
!      
      open(1,file='CONTCAR',status='old')
      write(*,*) 'How many types of atoms in the system?'
      read(*,*) ntype
      write(*,*) 'What are they?'
      do i = 1, ntype
         read(*,*)typename(i)
      enddo  
      read(1,'(a8)') title
      read(1,*) lattice
      do i = 1, 3
         read(1,*) (cell(i,j), j = 1, 3)
      enddo
      read(1,*) (numatom(i), i = 1, ntype)
      read(1,*) mode
      read(1,*) qformat
      na = 0
      do i = 1, ntype
         na = na + numatom(i)
      enddo
      
      write(*,*) 'Total number of atoms:', na
      do i = 1, na
         read(1,*) c(i,1),c(i,2),c(i,3),fix(i,1),fix(i,2),fix(i,3)    
      enddo
   
      close(1)
!
!----Assign atomic numbers to atoms read from CONTCAR.
!
      i = 1 
      do j = 1, ntype
         do k = 1, numatom(j)
           if (typename(j).eq.'H') then
              atnumtemp = 1
           elseif (typename(j).eq.'C') then
              atnumtemp = 6
           elseif (typename(j).eq.'O') then
              atnumtemp = 8
           elseif (typename(j).eq.'PT'.or.typename(j).eq.'Pt')then
              atnumtemp = 78
           endif 
           atnum(i) = atnumtemp
           i = i + 1 
         enddo 
      enddo                    
      
!
!-----Convert fraction coordinates into xyz.
!
      do i = 1, 3
        do j = 1, 3
           cell(i,j) = lattice * cell(i,j)
         enddo
      enddo
      do i = 1, 3
         write(*,*) (cell(i,j), j = 1, 3)
      enddo
      c = matmul(c,cell)
!
!----Read OUTCAR.
!      
      open(3,file=fnameout,status='replace') 
!
!-----Write output file header.
!
      write(3,'(a46)')' Entering Gaussian System, Link 0=/usr/g98/g98'
      write(3,'(a17)')' Initial command:'
      write(3,'(a48,a42)')
     $                ' /usr/g98/l1.exe /scratch/g98scratch/hhc/Gau-',
     $                '28628.inp -scrdir=/scratch/g98scratch/hhc/'
      write(3,'(a36,a16)')' Entering Link 1 = /usr/g98/l1.exe',
     $                    ' PID=     28630.'
      write(3,'(a30)')'                                       '
      write(3,'(a40,a19)')' Copyright (c) 1988,1990,1992,1993,1995,',
     $                    '1998 Gaussian, Inc.'
      write(3,'(a38)')'                  All Rights Reserved.'
      write(3,'(a41)')' ***************************************'
      write(3,'(a41)')' Gaussian 98:  SGI-G98RevA.6 19-Oct-1998'
      write(3,'(a41)')'                 7-Jan-2000             '
      write(3,'(a41)')' ***************************************'
      write(3,'(a41)')'                                        '
      write(3,'(a41)')'                                        '
      write(3,'(a46)')'                        Standard orientation:'
      write(3,'(a30,a40)')' -----------------------------',
     $                    '----------------------------------------'
      write(3,'(a29,a37)')' Center     Atomic     Atomic',
     $            '              Coordinates (Angstroms)'
      write(3,'(a29,a38)')' Number     Number      Type',
     $                    '             X           Y           Z'
      write(3,'(a30,a40)')' -----------------------------',
     $                    '----------------------------------------'
      do i = 1, na
         write(3,'(i5,i11,i14,4x,3f12.6)')i,atnum(i),0,(c(i,j),j = 1, 3)
      enddo
            write(3,'(a30,a40)')' -----------------------------',
     $                    '----------------------------------------'
      write(3,'(a41)')'                                        '
      write(3,'(a22,a29)')'   102 basis functions',
     $                    '      192 primitive gaussians'
      write(3,'(a46)')'    21 alpha electrons       21 beta electrons'
      write(3,'(a41)')'                                        '
      write(3,'(a41)')'                                        '
      write(3,'(a31,a26)')' Harmonic frequencies (cm**-1),',
     $                    ' IR intensities (KM/Mole),'
      write(3,'(a40,a29)')' Raman scattering activities (A**4/AMU),',
     $                    ' Raman depolarization ratios,'
      write(3,'(a38,a34)')' reduced masses (AMU), force constants',
     $                    ' (mDyne/A) and normal coordinates:'
!
!----Read frequencies.
! 
      open(2,file='OUTCAR',status='old')    
      nfreq = 0 
      nfreqmax = 0
      do i = 1, na
         do j = 1, 3
            if (fix(i,j).eq.'T') nfreqmax = nfreqmax + 1
         enddo
      enddo
      write(*,'(a10,i3,a20)')'There are ',nfreqmax,
     $                       ' degrees of freedom.'
   51 read(2,'(a80)') qline
      if (qline(59:60).eq.'cm') then
         nfreq = nfreq + 1
         qfreq(nfreq) = qline(47:57) 
         read(2,*) qX, qY, qZ, qdX, qdY, qdZ
!
!-----Write frequency.
! 
         write(3,'(a18,i3)')'                  ',i
         write(3,'(a22)')'                    A'''
         write(3,'(a17,a9)')' Frequencies --  ', qfreq(nfreq)
         write(3,'(a15,a11)')' Red. masses --','     0.0000'
         write(3,'(a15,a11)')' Frc consts  --','     0.0000'
         write(3,'(a15,a11)')' IR Inten    --','     0.0000'
         write(3,'(a15,a11)')' Raman Activ --','     0.0000'
         write(3,'(a15,a11)')' Depolar     --','     0.0000'
         write(3,'(a29)')' Atom AN      X      Y      Z'
         do i = 1, na        
            read (2,*) ctemp(i,1),ctemp(i,2),ctemp(i,3),
     $                disp(i,1),disp(i,2),disp(i,3) 
         enddo
      else
         go to 51 
      endif
      do j = 1, na
         write(3,'(2i4,2x,3f7.2)') j, atnum(j), (disp(j,k),k=1,3)
      enddo
      if (nfreq.lt.nfreqmax) then
         go to 51
      else
         go to 52
      endif
   52 continue
      write(3,'(a35)')' Normal termination of Gaussian 98.'
      close(2)
      close(3)
      stop 'Program successfully completed!'
      end program freqanal
