# reads a SU2 native mesh file (.su2) and extract the node coordinates

# the SU2 mesh must be structured, and nodes must be listed in increasing order,
# x, then y,
# for example, for a domain x: [0, 1], y: [0, 2], the nodes with y = 0, x from 
# 0 to 1 would be listed first, then the nodes with y = delta_y and x from
# 0 to 1, etc.

# this script wil verify the sortedness and print the number of nodes in 
# each direction

# 2D meshes only (for now)

# ARGS[1]  = name of .su2 file (including extension)

if length(ARGS) != 1
  println(STDERR, "Usage: extract_su2.jl fname.su2")
  exit(1)
end

function getCoords(fname)

  f = open(fname, "r")
  reading_points = false
  coords = Array(Float64, 0, 0)
  curr_row = 1  # current row of coords
  npoints = 0  # number of points


  for line in eachline(f)
    # find the line with NPOIN

    if startswith(line, "NPOIN")  # start of coordinate section
      # parse the line 
      idx = findfirst(line, '=')

      # TODO: check that there is not a comment on this line
      npoints = parse(Int, line[(idx+1):end])
      reading_points=true
      coords = zeros(Float64, npoints, 2)

      println("npoints = ", npoints)

    elseif reading_points  # this line contains the coords of a point
      data = split(line)
      @assert length(data) == 3  # x, y, index
      coords[curr_row, 1] = parse(Float64, data[1])
      coords[curr_row, 2] = parse(Float64, data[2])
      curr_row += 1

      if curr_row > npoints
        reading_points = false
      end
    end

  end  # end for loop

  close(f)
  check_sorted(coords)
  @assert curr_row == npoints + 1

  return coords
end


function check_sorted(coords)

  # verify that there are the same number of points in each section of the array

  # find all crossover points
  idx = 0

  cross_pts = Array(Int, 0)
  for i=2:size(coords, 1)
    if coords[i, 1] < coords[i-1, 1]
      push!(cross_pts, i)
    end
  end

  # verify the crossover points are equally spaced
  dist1 = cross_pts[2] - cross_pts[1]
  for i=2:length(cross_pts)
    disti = cross_pts[i] - cross_pts[i-1]
    @assert disti == dist1
  end

  # check that the y coordinates are sorted
  for i=2:length(cross_pts)
    idx1 = cross_pts[i-1]
    idx = cross_pts[i]
    @assert coords[idx, 2] > coords[idx1, 2]
  end

end



coords = getCoords(ARGS[1])
check_sorted(coords)
writedlm("coords.dat", coords)


