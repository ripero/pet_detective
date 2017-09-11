PROGRAM pd

   use, intrinsic :: iso_fortran_env, only: stderr=>error_unit, &
                                            stdout=>output_unit
   implicit none

#ifdef DEBUG
   logical :: debug = .true.
#else
   logical :: debug = .false.
#endif

   type POINT
      integer :: x
      integer :: y
      integer, allocatable :: neighbours(:) ! ids of neighbouring nodes
      integer, allocatable :: dists(:)  ! distances to all nodes
   end type

   integer :: num_points
   type(POINT), allocatable :: graph(:)

   integer, allocatable :: full_lattice(:,:)

   integer :: len_x ! dimension along x
   integer :: len_y ! dimension along y

   ! Points missing in the full lattice.
   integer :: num_missing_points
   integer, allocatable :: missing_points(:,:)

   ! Links missing in the full lattice.
   integer :: num_missing_links
   integer, allocatable :: missing_links(:,:)

   ! Auxiliary variables.
   integer :: i, j, k, dx, dy, num_links
   integer :: aux_links(4)
   logical :: is_link
   integer, allocatable :: aux_dists(:)
   integer, allocatable :: aux_dists_2D(:,:)

   ! Car
   integer, allocatable :: node_sequence(:)
   integer :: car_max_dist
   integer :: car_dist
   integer :: num_pets_in_car
   integer :: car(4)
   integer :: num_steps
   integer :: step
   integer :: num_missing_pets
   integer, allocatable :: missing_pets(:)

   integer, save :: total_calls_DFS = 0
   logical, save :: solution_found

   ! List of origin/destination nodes.
   ! This includes the car origin and all nodes with pets or homes.
   integer :: car_origin
   integer :: xp, yp, xh, yh, xc, yc
   integer :: num_pets
   integer, allocatable :: pets(:)
   integer, allocatable :: homes(:)

   ! I/O auxiliary variables.
   integer :: fid

   ! Determine size of lattice.
   open(newunit=fid, file='lengths.dat', &
        status='old', action='read', form='formatted')
   read(fid,*) len_x, len_y
   close(fid)

   ! Determine missing points
   open(newunit=fid, file='missing_points.dat', &
        status='old', action='read', form='formatted')
   read(fid,*) num_missing_points
   allocate(missing_points(2,num_missing_points))
   do i = 1, num_missing_points
      read(fid,*) missing_points(:,i)
   end do
   close(fid)

   ! Create full lattice.
   allocate(full_lattice(len_x,len_y))
   full_lattice = 1
   do i = 1, num_missing_points
      full_lattice(missing_points(1,i), missing_points(2,i)) = 0
   end do

   ! List of missing points is no longer needed
   deallocate(missing_points)

   ! Create graph
   num_points = len_x * len_y - num_missing_points
   allocate(graph(num_points))
   k = 1
   do j = 1, len_y
      do i = 1, len_x
         if (full_lattice(i,j).ne.0) then
            graph(k)%x = i
            graph(k)%y = j
            full_lattice(i,j) = k
            k = k + 1
         end if
      end do
   end do

   ! Determine missing links
   open(newunit=fid, file='missing_links.dat', &
        status='old', action='read', form='formatted')
   read(fid,*) num_missing_links
   allocate(missing_links(4,num_missing_links))
   do i = 1, num_missing_links
      ! Missing link is between (xa,ya) and (xb,yb)
      ! missing_links(1,num_missing_links) = xa
      ! missing_links(2,num_missing_links) = ya
      ! missing_links(3,num_missing_links) = xb
      ! missing_links(4,num_missing_links) = yb
      read(fid,*) missing_links(:,i)
   end do
   close(fid)

   ! Create lists of neighbours
   do i = 1, num_points
      num_links = 0
      aux_links(:) = 0
      do j = 1, num_points
         if (i==j) cycle
         dx = abs(graph(i)%x-graph(j)%x)
         if (dx.le.1) then
            dy = abs(graph(i)%y-graph(j)%y)
            if ((dy.le.1).and.(dx.ne.dy)) then
               ! Points i and j are contiguous in the lattice.
               is_link = .true.
               do k = 1, num_missing_links
                  if ((missing_links(1,k) == graph(i)%x .and. &
                       missing_links(2,k) == graph(i)%y .and. &
                       missing_links(3,k) == graph(j)%x .and. &
                       missing_links(4,k) == graph(j)%y) .or. &
                      (missing_links(1,k) == graph(j)%x .and. &
                       missing_links(2,k) == graph(j)%y .and. &
                       missing_links(3,k) == graph(i)%x .and. &
                       missing_links(4,k) == graph(i)%y)) then
                     is_link = .false.
                     exit
                  end if
               end do
               if (is_link) then
                  num_links = num_links + 1
                  aux_links(num_links) = j
               end if
            end if
         end if
      end do
      allocate(graph(i)%neighbours(num_links))
      graph(i)%neighbours(:) = aux_links(1:num_links)
   end do

   if (debug) then
      write(*,*) full_lattice
      do i = 1, num_points
         write(*,*) graph(i)%x, graph(i)%y, graph(i)%neighbours
      end do
   end if

   ! Check consistency of neighbours
   do i = 1, num_points
      do j = 1, num_points
         if ( any(j==graph(i)%neighbours(:)) .and. .not.&
              any(i==graph(j)%neighbours(:)) ) then
            STOP 'ERROR: INCONSISTENT NEIGHBOURS'
         end if
      end do
   end do

   deallocate(missing_links)

   ! Print lattice
   call print_lattice(graph, full_lattice, full_lattice)

   ! Compute distances between nodes
   write(stdout,'(a)', advance='no') 'Computing distances between nodes...'
   if (debug) write(*,'(a)') ''
   allocate(aux_dists(num_points))
   do i = 1, num_points
      call Dijkstra(graph, i, aux_dists)
      graph(i)%dists = aux_dists ! this will allocate on-the-fly the LHS array
   end do
   deallocate(aux_dists)

   write(stdout,'(a)') '... done'

   ! Check consistency of distances
   do i = 1, num_points
      do j = 1, num_points
         if (graph(i)%dists(j) .ne. graph(j)%dists(i)) then
            STOP 'ERROR: INCONSISTENT DISTANCES'
         end if
      end do
   end do

   ! Write graph with distances
   if (debug) then
      do i = 1, num_points
         write(stdout,'(a,i0,a,i0,1X,i0)') 'Distances to point ', i, &
              ' with coordinates ', graph(i)%x, graph(i)%y
         aux_dists_2D = full_lattice
         do j = 1, num_points
            aux_dists_2D(graph(j)%x,graph(j)%y) = graph(i)%dists(j)
         end do
         call print_lattice(graph, full_lattice, aux_dists_2D)
      end do
      deallocate(aux_dists_2D)
   end if

   ! Read pets
   open(newunit=fid, file='pets.dat', &
        status='old', action='read', form='formatted')
   read(fid,*) num_pets
   allocate(pets(num_pets))
   allocate(homes(num_pets))
   do i = 1, num_pets
      read(fid,*) xp, yp, xh, yh
      pets(i) = full_lattice(xp,yp)
      homes(i) = full_lattice(xh,yh)
   end do
   if (any(pets==0)) STOP 'ERROR: PET CANNOT BE IN INEXISTING NODE.'
   if (any(homes==0)) STOP 'ERROR: HOME CANNOT BE IN INEXISTING NODE.'
   close(fid)

   ! Read car
   open(newunit=fid, file='car.dat', &
        status='old', action='read', form='formatted')
   read(fid,*) car_max_dist
   read(fid,*) xc, yc
   car_origin = full_lattice(xc, yc)
   if (any(homes==0)) STOP 'ERROR: CAR ORIGIN CANNOT BE IN INEXISTING NODE.'
   close(fid)

   deallocate(full_lattice)

   ! Initialisations
   num_steps = 2*num_pets + 1
   allocate(node_sequence(num_steps))

   allocate(missing_pets(num_pets))
   missing_pets = pets
   num_missing_pets = num_pets

   node_sequence(1) = car_origin
   step = 1
   car(:) = 0
   car_dist = 0
   num_pets_in_car = 0
   solution_found = .false.

   ! Call DFS (recursive subroutine).
   call DFS(node_sequence, step, car, car_dist, num_pets_in_car, &
        missing_pets, num_missing_pets)

   ! Output final result.
   write(stdout,'(a,i0)') 'Number of calls to DFS: ', total_calls_DFS
   if (solution_found) then
      write(stdout,'(a)') 'SOLUTION FOUND:'
      write(stdout,'(100000(i0, 2X))') node_sequence
   else
      STOP 'ERROR: COULD NOT FIND A SOLUTION.'
   end if

CONTAINS

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE DFS(seq, step, old_car, car_dist, old_num_pets_in_car, &
        old_missing_pets, old_num_missing_pets)

      integer, intent(inout) :: seq(:)
      integer, intent(inout) :: step
      integer, intent(in) :: car_dist
      integer, intent(in) :: old_num_pets_in_car
      integer, intent(in) :: old_num_missing_pets
      integer, intent(in) :: old_car(:)
      integer, intent(in) :: old_missing_pets(:)
      integer :: new_car_dist
      integer :: num_pets_in_car
      integer :: num_missing_pets
      integer :: car(4)
      integer :: missing_pets(size(old_missing_pets))

      ! Local variables.
      integer :: num_missing_dests
      integer :: num_possible_dests
      integer :: possible_dests(num_pets)
      integer :: i

      total_calls_DFS = total_calls_DFS + 1

      ! Update lists after last step.
      if (step > 1) then
         if (any(seq(step)==pets)) then
            ! in the last step we picked up a pet
            num_pets_in_car = old_num_pets_in_car + 1
            ! add the home of the pet to the list of addresses of the car
            i = myfindloc(pets, seq(step))
            car(1:old_num_pets_in_car) = old_car(1:old_num_pets_in_car)
            car(num_pets_in_car) = homes(i)
            ! update the list of pets that still need to be picked
            i = myfindloc(old_missing_pets, seq(step))
            missing_pets(1:i-1) = old_missing_pets(1:i-1)
            num_missing_pets = old_num_missing_pets - 1
            missing_pets(i:num_missing_pets) = old_missing_pets(i+1:old_num_missing_pets)
         else
            ! in the last step we dropped a pet at a home
            ! update the list of pets in car
            i = myfindloc(old_car, seq(step))
            car(1:i-1) = old_car(1:i-1)
            num_pets_in_car = old_num_pets_in_car - 1
            car(i:num_pets_in_car) = old_car(i+1:old_num_pets_in_car)
            ! copy list of pets - they are the same
            num_missing_pets = old_num_missing_pets
            missing_pets(1:num_missing_pets) = old_missing_pets(1:num_missing_pets)
         end if
      else
         num_pets_in_car = old_num_pets_in_car
         num_missing_pets = old_num_missing_pets
         missing_pets(1:num_missing_pets) = old_missing_pets(1:num_missing_pets)
      end if

      step = step + 1

      ! Create list of possible destinations.
      if (num_pets_in_car==4) then
         num_possible_dests = 4
         possible_dests(1:4) = car(1:4)
      else
         possible_dests(1:num_pets_in_car) = car(1:num_pets_in_car)
         num_possible_dests = num_pets_in_car + num_missing_pets
         possible_dests(num_pets_in_car+1:num_possible_dests) = missing_pets(1:num_missing_pets)
      end if

      num_missing_dests = num_steps - step

      ! Loop over possible destinations
      do i = 1, num_possible_dests
         seq(step) = possible_dests(i)
         new_car_dist = car_dist + graph(seq(step-1))%dists(seq(step))
         if (new_car_dist + num_missing_dests > car_max_dist) then
            ! Discard: more missing destinations than possible steps.
            cycle
         else if (num_missing_dests == 0) then
            solution_found = .true.
            return
         else
            call DFS(seq, step, car, new_car_dist, num_pets_in_car, &
                 missing_pets, num_missing_pets)
            if (solution_found) return
         end if
      end do
      seq(step) = 0
      step = step - 1

   END SUBROUTINE DFS

!-------------------------------------------------------------------------------

   integer FUNCTION myfindloc(array, val)
      implicit none

      integer, intent(in) :: array(:)
      integer, intent(in) :: val

      integer :: i

      do i = lbound(array,1), ubound(array,1)
         if (val == array(i)) then
            myfindloc = i
            return
         end if
      end do
      write(stderr,*) 'VAL', val
      write(stderr,*) 'ARRAY', array
      STOP 'ERROR: MYFINDLOC COULD NOT FIND VAL IN ARRAY.'
   END FUNCTION myfindloc

!-------------------------------------------------------------------------------

!   FUNCTION car_distance(seq, step, g) result(d)
!      implicit none
!
!      integer, intent(in) :: seq(:)
!      integer, intent(in) :: step
!      type(POINT), intent(in) :: g(:)
!
!      integer :: d
!      integer :: i
!
!      d = 0
!      do i = 0, step-1
!         d = d + g(seq(i))%dists(j)
!      end do
!   END FUNCTION car_distance

!-------------------------------------------------------------------------------

   SUBROUTINE print_lattice(g, lattice, values_to_print)

      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      type(POINT), intent(in) :: g(:)
      integer,     intent(in) :: lattice(:,:)
      integer,     intent(in) :: values_to_print(:,:)

      integer :: lx, ly
      integer :: i, j

      integer :: node_fmt_len
      character(len=4) :: node_fmt  ! node format
      character(len=4) :: space_fmt ! space format
      character(len=10) :: v_link_fmt ! vertical link format
      character(len=3) :: aux_h_link
      character(len=3), parameter :: h_link    = '---'
      character(len=3), parameter :: h_no_link = '   '
      character(len=1) :: aux_v_link
      character(len=1), parameter :: v_link    = '|'
      character(len=1), parameter :: v_no_link = ' '

      node_fmt_len = int(log10(real(maxval(values_to_print),kind=real64))) + 1

      write(node_fmt, '(a,i0,a)') '(i', node_fmt_len, ')'
      write(space_fmt, '(a,i0,a)') '(a', node_fmt_len, ')'
      if ((node_fmt_len-1)/2 > 0) then
         write(v_link_fmt, '(a,i0,a,i0,a)') '(', &
              node_fmt_len-(1+(node_fmt_len-1)/2), 'X,a1,', &
              (node_fmt_len-1)/2, 'X)'
      else if (node_fmt_len-(1+(node_fmt_len-1)/2) > 0) then
         write(v_link_fmt, '(a,i0,a)') '(', &
              node_fmt_len-(1+(node_fmt_len-1)/2), 'X,a1)'
      else
         v_link_fmt = '(a1)'
      end if

      lx = size(lattice,1)
      ly = size(lattice,2)

      ! Print lattice
      do j = 1, ly
         do i = 1, lx
            if (lattice(i,j).ne.0) then
               write(*,node_fmt, advance='no') values_to_print(i,j)
            else
               write(*,space_fmt, advance='no') ''
            end if
            aux_h_link = h_no_link
            if (i < lx) then
               if (lattice(i,j).ne.0) then
                  if (lattice(i+1,j).ne.0) then
                     if (any(g(lattice(i,j))%neighbours == lattice(i+1,j))) then
                        aux_h_link = h_link
                     end if
                  end if
               end if
               write(*,'(a3)', advance='no') aux_h_link
            else ! i == lx
               write(*,'(a)') ''
            end if
         end do
         if (j < ly) then
            do i = 1, lx
               aux_v_link = v_no_link
               if (lattice(i,j).ne.0) then
                  if (lattice(i,j+1).ne.0) then
                     if (any(g(lattice(i,j))%neighbours == lattice(i,j+1))) then
                        aux_v_link = v_link
                     end if
                  end if
               end if
               write(*,v_link_fmt, advance='no') aux_v_link
               if (i < lx) then
                  write(*,'(a3)', advance='no') h_no_link
               else ! i == lx
                  write(*,'(a)') ''
               end if
            end do
         end if
      end do

   END SUBROUTINE print_lattice

!-------------------------------------------------------------------------------

   SUBROUTINE Dijkstra(g, origin, dist)

      implicit none

      type(POINT), intent(in)  :: g(:)   ! Graph
      integer,     intent(in)  :: origin
      integer,     intent(out) :: dist(size(g))

      integer :: previous(size(g))
      logical :: still_active(size(g))

      integer :: u, v, iv, alt

      dist = HUGE(dist)
      previous(i) = 0

      dist(origin) = 0
      still_active = .true.

      do while(any(still_active))
         u = minloc(dist,dim=1,mask=still_active)
         still_active(u) = .false.
         do iv = 1, ubound(graph(u)%neighbours,1)
            v = graph(u)%neighbours(iv)
            if (still_active(v)) then
               alt = dist(u) + 1 ! in the general case, the 1 would be replaced
                                 ! by the weight of the link (u,v)
               if (alt < dist(v)) then
                  dist(v) = alt
                  previous(v) = u
               end if
            end if
         end do
      end do

      if (debug) then
         write(*,*) 'origin', origin
         write(*,*) 'previous', previous
         write(*,*) 'dist', dist
      end if

   END SUBROUTINE Dijkstra

END PROGRAM pd
