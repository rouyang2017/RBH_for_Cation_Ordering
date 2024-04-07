PROGRAM RBH
! Recommendation-based Basin Hopping for thermodynaimc global optimization of cation ordering

implicit none
integer nstep,i,j,k,l,natom_tot,natom_s1,natom_s2,natom_s3,nswap,ntype,istep,nstr,nMLAB_max,&
        nhop_ml,besttry,nvA,nvB,nvO,nlatt
real*8 prob,rn,dE,temperature,beta,markovE,latt(3,3),energy,emin,ethis
character line*100,command*100,header(9)*100
logical converged
real*8,allocatable:: coord(:,:),oldcoord(:,:)
character,allocatable:: subdir(:)*10,dynamics(:)*9
integer,allocatable:: siteID(:),atype(:),natype(:)

!------------------------------------------------------------------
! Swap between atoms at the same site_i, e.g.
! coord_1    coord_1    coord_1    T   T   T   site_1
! coord_2    coord_2    coord_2    T   T   T   site_1
! There can be maximally 3 types of site: site_1, site_2, and
! site_3 in the POSCAR.
!-----------------------------------------------------------------
command='srun vasp_gam-6.3.2'
natom_tot=300  ! total number of atom in POSCAR
ntype=5        ! number of atom type in POSCAR
nstep=10000  ! total number of MC steps
temperature=2500  ! for Metropolis MC
beta=11604.8/temperature    ! 1/kT
nswap=1    ! number atom-swap per mc-step
nMLAB_max=60  ! max number of structures in ML_AB
nhop_ml=100   ! hopping optimization via ML
nvA=0         ! number of A-site vacancies
nvB=0         ! number of B-site vacancies
nvO=0         ! number of O-site vacancies
!------------------------------------------------------------

allocate(subdir(nstep))
allocate(coord(natom_tot+nvA+nvB+nvO,3))
allocate(oldcoord(natom_tot+nvA+nvB+nvO,3))
allocate(siteID(natom_tot+nvA+nvB+nvO))
allocate(dynamics(natom_tot+nvA+nvB+nvO))
allocate(atype(natom_tot+nvA+nvB+nvO))
allocate(natype(ntype+1))

call random_seed()
open(1,file='RBH.out',status='replace')
write(1,'(a)') '! Configuration Optimization of Solid Solution with Monte Carlo &
                &powered by DFT or Machine Learning!'
open(11,file='swapinfo',status='replace')
write(11,'(3a)') 'Atoms(coordinates) swap information at each step !'

! Read the initial POSCAR to get the site information
open(2,file='POSCAR',status='old')
do i=1,9   ! skip the first 9 lines
  read(2,'(a)') header(i)
  if(i==7) read(header(i),*) natype(1:ntype)
end do

markovE=0.d0
siteID=0   ! which site
natom_s1=0  ! number of atoms for configuprobnal search at site i
natom_s2=0  ! i.e. random occupation of atoms at all the site i.
natom_s3=0
atype=-1
coord=0.d0

! initial coord, slective dynamics, swap
do i=1,natom_tot+nvA+nvB+nvO
   read(2,'(a)') line
   dynamics(i)='         '
   read(line,*) coord(i,1:3),dynamics(i)(3:3),dynamics(i)(6:6),dynamics(i)(9:9)
   if (index(line,'site_1')/=0) then
      siteID(i)=1
      natom_s1=natom_s1+1
   elseif(index(line,'site_2')/=0) then
      siteID(i)=2
      natom_s2=natom_s2+1
   elseif(index(line,'site_3')/=0) then
      siteID(i)=3
      natom_s3=natom_s3+1
   end if
   do j=1,ntype
     if( (j==1 .and. i<=natype(1)) .or. ( i>sum(natype(1:j-1)) .and. i<=sum(natype(1:j)))) then
        atype(i)=j   ! atoms
     else if (i>natom_tot) then
        atype(i)=0   ! vacancies
     end if
   end do
end do
oldcoord=coord
close(2)


!------------------------------------------
do istep=1,nstep

   write(1,'(/a,i6.6)'), 'Step ',istep
   
   !*******************************
   ! create subdir and run vasp
   !*******************************
   subdir(istep)(1:4)='step'
   write(subdir(istep)(5:10),'(i6.6)'),istep

   print *,"mystep: ", istep
   if(istep==1) then
     call system('mkdir '//subdir(istep)//'; cp INCAR_0 POSCAR KPOINTS POTCAR '//subdir(istep)//';& 
                  cd '//subdir(istep)//'; mv INCAR_0 INCAR; '//trim(adjustl(command))//'; &
                  cp ML_ABN ../ML_AB; cp ML_FFN ../ML_FF; cd ..')           
   else if(istep>1) then
     call system('mkdir '//subdir(istep)//'; cp ML_AB INCAR_1 KPOINTS POTCAR '//subdir(istep)//';& 
                  cd '//subdir(istep)//'; mv ../POSCAR_next POSCAR; mv INCAR_1 INCAR; &
                  '//trim(adjustl(command))//';  cp ML_ABN ../ML_AB; cp ML_FFN ../ML_FF; cd ..')
   end if   
   open(3,file='ML_AB',status='old')
     do j=1,4
       read(3,*)
     end do
     read(3,*) nstr
   close(3)

   if(nstr>=nMLAB_max) then
       call reduce_MLAB
       call system('mkdir newMLAB; cp INCAR_3 KPOINTS POSCAR POTCAR ML_AB newMLAB; cd newMLAB;mv INCAR_3 INCAR; &
                    '//trim(adjustl(command))//'; cp ML_ABN ../ML_AB; cd ..; rm -rf newMLAB ')
   end if
   if(istep>10) call system('rm -rf '//subdir(istep-10)//' ')
   !----------------------------------------------
   ! read the converged coordinates and energy
   !----------------------------------------------
   converged=.false.
   energy=0.d0
   latt=0.d0

   nlatt=0
    open(2,file=subdir(istep)//'/OUTCAR',status='old')
    do while(.not. eof(2))
      read(2,'(a)') line
      if(index(line,'direct lattice vectors')/=0) nlatt=nlatt+1
      if(index(line,'reached required accuracy - stopping structural energy minimisation')/=0) converged=.true.
    end do
    close(2)
   

   open(2,file=subdir(istep)//'/OUTCAR',status='old')
   k=0
   do while(.not. eof(2))
      read(2,'(a)') line
      if (index(line,'direct lattice vectors')/=0) then
          k=k+1
          if (k==nlatt) then
              do l=1,3
                 read(2,*) latt(l,:)
              end do
          end if
      end if
      if(k==nlatt .and. index(line,'POSITION')/=0) then
         read(2,*)
         do l=1,natom_tot
            read(2,*) coord(l,:)
         end do

      end if
!      if(k==nlatt .and. index(line,'energy(sigma->0) =')/=0 ) then
!        l=index(line,'energy(sigma->0) =')
!        read(line(l+18:),*) energy
!      end if
   end do
   close(2)
   

   open(2,file=subdir(istep)//'/OSZICAR',status='old')
   do while (.not. eof(2))
       read(2,'(a)') line
       k=index(line,'E0=')
       if(k/=0) read(line(k+3:),*) energy
   end do
   close(2)
   
   open(2,file=subdir(istep)//'/POSCAR',status='old')
   do k=1,9
      read(2,*)
   end do
   do k=1,natom_tot+nvA+nvB+nvO
      if(k<=natom_tot) then
          read(2,*)
      else
          read(2,*) coord(k,:)
      end if
   end do
   close(2)
   
   !-----------------------------------------------------------------
   ! output lattice vectors, coordinates, energy, and convergence
   !-----------------------------------------------------------------
   write(1,'(a)') 'Lattice vectors: '
   do j=1,3
      write(1,'(3f20.10)') latt(j,:)
   end do
   write(1,'(a)') 'Atomic coordinates: '
   do j=1,natom_tot+nvA+nvB+nvO
      if(j<=natom_tot) then
          write(1,'(3f20.10)') coord(j,:)
      else
          write(1,'(3f20.10,a)') coord(j,:),'   vacancy'
      end if
   end do
   write(1,'(a,l)') 'Structure relaxation converged? ',converged
   !*********************************
   ! Metropolis MC selection
   !*********************************
   call random_number(rn) 
   dE=energy-markovE
   prob=exp(-dE*beta)
   if(prob > rn .and. energy<0.0 ) then
         write (1,'(a,f8.4,a,f8.4)'),'Probability ',prob, ' > random r=',rn
         write(1,'(a,i6.6,a,f20.10,a)') 'Step: ',istep,'  Total energy: ',energy,'   Accepted!'
         markovE=energy
         oldcoord=coord
   else
         write (1,'(a,f8.4,a,f8.4)'),'Probability ',prob, ' < random r=',rn
         write(1,'(a,i6.6,a,f20.10,a)') 'Step: ',istep,'  Total energy: ',energy,'   Rejected!'
         coord=oldcoord
   endif

   if(istep/=nstep) then
      write(11,'(/a,i5)') 'Step ',istep+1

      write(1,'(a)') 'Hopping optimization with ML (the one with lowest energy will be adopted):'
      emin=0.d0
      besttry=0
      do j=1,nhop_ml
          write(11,'(a,i5)') 'Try ',j
          call  POSCAR_update
          call system('mkdir coord_opt; cp ML_FF POSCAR_next INCAR_2 KPOINTS POTCAR coord_opt; & 
                  cd coord_opt; mv POSCAR_next POSCAR; mv INCAR_2 INCAR; &
                  '//trim(adjustl(command))//'; cd ..')
          open(50,file='coord_opt/OSZICAR',status='old')
          do while (.not. eof(50))
              read(50,'(a)') line
              k=index(line,'E0=')
              read(line(k+3:),*) ethis
          end do
          close(50)
          if(ethis<emin) then
             emin=ethis
             besttry=j
             call system('cd coord_opt; cp POSCAR POSCAR_min; cd ..')
          end if
          write(1,'(a,i5,f20.10)') 'Hopping_try, ML-E0: ', j, ethis
          if(mod(j,20)==0 .and. emin < markovE) exit
      end do
      call system('cp coord_opt/POSCAR_min POSCAR_next; rm -rf coord_opt')

      write(1,'(a,i5)') 'Adopted try: ',besttry
      write(1,'(a)') 'New configuration generated !'
      write(1,'(a)') '-----------------------------'
   end if

end do

close(1)
close(11)

deallocate(subdir)
deallocate(coord)
deallocate(oldcoord)
deallocate(siteID)
deallocate(dynamics)
deallocate(atype)
deallocate(natype)

contains

subroutine POSCAR_update
integer i,j,k,l,nnn,pos1(natom_s1),pos2(natom_s2),pos3(natom_s3),thissite,icount,nsite
real*8 newcoord(natom_tot+nvA+nvB+nvO,3), rand, coord_tmp(3)
logical doswap
logical,allocatable:: usable_s1(:),usable_s2(:),usable_s3(:)

nsite=0
! initialization
if(natom_s1>0) then
  allocate(usable_s1(natom_s1))
  usable_s1=.true.
  nsite=nsite+1
end if
if(natom_s2>0) then
   allocate(usable_s2(natom_s2))
   usable_s2=.true. 
   nsite=nsite+1
end if
if(natom_s3>0) then
   allocate(usable_s3(natom_s3))
   usable_s3=.true.
   nsite=nsite+1
end if

! positions of the atoms at each site
newcoord=coord
j=0; k=0; l=0
do i=1,natom_tot+nvA+nvB+nvO
    if(siteID(i)==1) then
        j=j+1
        pos1(j)=i   ! position of this atom in the POSCAR
    else if(siteID(i)==2) then
        k=k+1
        pos2(k)=i
    else if(siteID(i)==3) then
        l=l+1
        pos3(l)=i
    endif
end do

! cation swap
nnn=0  
icount=0

do while (nnn<nswap)
    icount=icount+1
    if(icount>100*(natom_tot+nvA+nvB+nvO)) exit ! avoid endless loop
    call random_number(rand)
    thissite=max(ceiling(3*rand),1)

    if(thissite==1 ) then  
        do while(any(usable_s1))   ! random selection of A @ site 1
           call random_number(rand)
           k=max(1,nint(rand*natom_s1))
           if(usable_s1(k)) then  ! selected
              usable_s1(k)=.false.
              exit
           end if
        end do

        doswap=.false.
        do l=1,natom_s1
           if(atype(pos1(l))/=atype(pos1(k)) .and. usable_s1(l)) then
           doswap=.true.
           exit
           end if
        end do

        do while(any(usable_s1) .and. doswap)      ! random selection of A' at site 1
           call random_number(rand)
           j=max(1,nint(rand*natom_s1))  
           if(usable_s1(j) .and. atype(pos1(j))/=atype(pos1(k))) then  ! selected
               coord_tmp=newcoord(pos1(k),:)
               newcoord(pos1(k),:)=newcoord(pos1(j),:)
               newcoord(pos1(j),:)=coord_tmp
               usable_s1(j)=.false.
               write(11,'(i5,a,i5)') pos1(k),'  <-->',pos1(j)
               nnn=nnn+1
               exit
           end if
        end do

    elseif(thissite==2 ) then  !site2
        do while(any(usable_s2))  
           call random_number(rand)
           k=max(1,nint(rand*natom_s2))
           if(usable_s2(k)) then
               usable_s2(k)=.false.
               exit
           end if
        end do
  
        doswap=.false.
         do l=1,natom_s2
            if(atype(pos2(l))/=atype(pos2(k)) .and. usable_s2(l)) then
                doswap=.true.
                exit
            end if
         end do

        do while(any(usable_s2) .and. doswap )   
           call random_number(rand)
           j=max(1,nint(rand*natom_s2))
           if(usable_s2(j) .and. j/=k .and. atype(pos2(j))/=atype(pos2(k)) ) then
               coord_tmp=newcoord(pos2(k),:)
               newcoord(pos2(k),:)=newcoord(pos2(j),:)
               newcoord(pos2(j),:)=coord_tmp
               usable_s2(j)=.false.
               write(11,'(i5,a,i5)') pos2(k),'  <-->',pos2(j)
               nnn=nnn+1
               exit
           end if
        end do

    elseif(thissite==3 ) then  !site3
        do while(any(usable_s3))
           call random_number(rand)
           k=max(1,nint(rand*natom_s3))
           if(usable_s3(k)) then
               usable_s3(k)=.false.
               exit
           end if
        end do

        doswap=.false.
         do l=1,natom_s3
            if(atype(pos3(l))/=atype(pos3(k)) .and. usable_s3(l)) then
                doswap=.true.
                exit
            end if
         end do

        do while(any(usable_s3) .and. doswap )
           call random_number(rand)
           j=max(1,nint(rand*natom_s3))
           if(usable_s3(j) .and. j/=k .and. atype(pos3(j))/=atype(pos3(k)) ) then
               coord_tmp=newcoord(pos3(k),:)
               newcoord(pos3(k),:)=newcoord(pos3(j),:)
               newcoord(pos3(j),:)=coord_tmp
               usable_s3(j)=.false.
               write(11,'(i5,a,i5)') pos3(k),'  <-->',pos3(j)
               nnn=nnn+1
               exit
           end if
        end do

    end if
end do


! create the new POSCAR with the newcoord
open(10,file='POSCAR_next',status='replace')
do i=1,9
   write(10,'(a)') trim(header(i))
end do

do i=1,natom_tot+nvA+nvB+nvO
  write(10,'(3f20.10,2a)') newcoord(i,:),'  ',trim(dynamics(i))
end do
close(10)


if(natom_s1>0) then
  deallocate(usable_s1)
end if
if(natom_s2>0) then
   deallocate(usable_s2)
end if
if(natom_s3>0) then
   deallocate(usable_s3)
end if

end subroutine POSCAR_update


subroutine reduce_MLAB
! discard the first half of structures in ML_AB 
integer i,j,istr
character line*200


open(4,file='ML_AB',status='old')
open(5,file='ML_AB_new',status='replace')

istr=0
do while( .not. eof(4))
   read(4,'(a)') line
   
   if(index(line,'Configuration num.')/=0) then
      istr=istr+1
      write(5,'(a,i7)') 'Configuration num.',istr
   else
      write(5,'(a)') trim(line)
   end if

   if(index(line,'The number of configurations')/=0) then
     read(4,'(a)') line
     write(5,'(a)') trim(line)
     read(4,'(a)') line
     write(5,'(i10)') int(nstr/2.0)
   end if

   if(index(line,'The numbers of basis sets per atom type')/=0) then
       write(5,'(a)') '--------------------------------------------------'
       i=1
       do while(i<=ntype)
          if(i+2<=ntype) then
             write(5,'(3i8)') 1, 1, 1
          else if(i+1<=ntype) then
             write(5,'(2i8)') 1, 1
          else if(i<=ntype) then
             write(5,'(i8)') 1
          end if
          i=i+3
       end do
       write(5,'(a)') '**************************************************'

       do j=1,ntype
         do while(index(line,'Basis set for')==0 .and. (.not. eof(4)))
            read(4,'(a)') line
         end do   
         write(5,'(a)') trim(line)      
         write(5,'(a)') '--------------------------------------------------'
         write(5,'(2i8)') 1, 1
         write(5,'(a)') '**************************************************'
         line=''
       end do

       i=0
       do while((.not. eof(4)) .and. i<=(nstr-int(nstr/2.0)) )
          read(4,'(a)') line
          if(index(line,'Configuration num.')/=0) i=i+1
       end do
       istr=1
       write(5,'(a,i7)') '     Configuration num.',istr
    end if
end do
close(4)
close(5)

call system('mv ML_AB_new ML_AB')

end subroutine reduce_MLAB


end


