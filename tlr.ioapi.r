# convenience methods for IOAPI files

library(ncdf4)

# get plot-relevant global attributes from IOAPI file, return list
get.plot.addrs.from.IOAPI <- function(
  ioapi.file       # class=ncdf4 from nc_open(...)
) {
  # Note: presumes lengths in meters!
  x.orig.m <- ncatt_get(ioapi.file, varid=0, attname="XORIG")$value
  x.orig.km <- x.orig.m/1000
  y.orig.m <- ncatt_get(ioapi.file, varid=0, attname="YORIG")$value
  y.orig.km <- y.orig.m/1000
  x.cell.size.m <- ncatt_get(ioapi.file, varid=0, attname="XCELL")$value
  x.cell.size.km <- x.cell.size.m/1000
  y.cell.size.m <- ncatt_get(ioapi.file, varid=0, attname="YCELL")$value
  y.cell.size.km <- y.cell.size.m/1000
  x.n.grid <- ncatt_get(ioapi.file, varid=0, attname="NCOLS")$value
  y.n.grid <- ncatt_get(ioapi.file, varid=0, attname="NROWS")$value

  attrs.list <- list()
  attrs.list$x.cell.centers.km <-
    seq(from=x.orig.km + x.cell.size.km/2, length=x.n.grid, by=x.cell.size.km)
  attrs.list$y.cell.centers.km <-
    seq(from=y.orig.km + y.cell.size.km/2, length=y.n.grid, by=y.cell.size.km)
  # start debugging for Doug Nychka Mon, 13 Feb 2012 21:33:36 -0700
  # print(paste('class(x.cell.centers.km)==', class(attrs.list$x.cell.centers.km), sep=""))
  # print(paste('class(y.cell.centers.km)==', class(attrs.list$y.cell.centers.km), sep=""))
  #   end debugging for Doug Nychka Mon, 13 Feb 2012 21:33:36 -0700
  attrs.list
} # end function get.plot.addrs.from.IOAPI

# if value of IOAPI global attribute with name given does not match
# desired value, make it so
fix.IOAPI.global.attribute <- function(
  ioapi.file,    # class=ncdf4 from nc_open(...)
  attr.name,     # string name of global attribute
  value.correct
) {
  # TODO: only call once, cache list(foo$hasatt, foo$value)
  if (ncatt_get(ioapi.file, varid=0, attname=attr.name)$hasatt) {
    value.current <- ncatt_get(ioapi.file, varid=0, attname=attr.name)$value
# start debugging
# this dies with error:  cat(sprintf('%s : value.current=%i, value.correct=%i\n', attr.name, value.current, value.correct))
#    cat(sprintf('%s : value.current=%f, value.correct=%f\n', attr.name, value.current, value.correct))
#   end debugging

# this is not sufficient for, e.g., vector testing:
#    if (value.current != value.correct) {
# gets warning message that
# > the condition has length > 1 and only the first element will be used
# so the test succeeds when it should fail.
    if (
        (class(value.current) != class(value.correct)) ||
        (length(value.current) != length(value.correct)) ||
        (value.current != value.correct)
       ) {
      ncatt_put(ioapi.file, varid=0, attname=attr.name, attval=value.correct,
        verbose=FALSE)
    }
  } # TODO: else, throw error (on not having the attribute)
  # return
  ioapi.file
} # end function fix.IOAPI.global.attribute

# if IOAPI global attribute name=NLAYS does not match actual layer count,
# make it so
fix.NLAYS <- function(
  ioapi.file,    # class=ncdf4 from nc_open(...)
  nlays.correct  # value
) {
  ioapi.file <- fix.IOAPI.global.attribute(
    ioapi.file, attr.name="NLAYS", value.correct=nlays.correct)
  # return
  ioapi.file
} # end function fix.NLAYS

# if IOAPI global attribute name=NVARS does not match actual variable count,
# make it so
fix.NVARS <- function(
  ioapi.file,    # class=ncdf4 from nc_open(...)
  nvars.correct  # value
) {
  ioapi.file <- fix.IOAPI.global.attribute(
    ioapi.file, attr.name="NVARS", value.correct=nvars.correct)
  # return
  ioapi.file
} # end function fix.NVARS

# if IOAPI global attribute name=VAR-LIST does not contain the "real" variable name(s)
# (minus the IOAPI accounting variable=TFLAG), make it so
fix.VARdashLIST <- function(
  source.file,  # class=ncdf4 from nc_open(...)
  names.vec,    # character vector of names of datavars we want to keep,
                # in the order we want to keep them!
  target.file
) {
  # TODO: more testing of inputs
  stopifnot(!is.element("TFLAG", names.vec)) # ASSERT: TFLAG not a member

  correct.names.string <- ""
  for (name in names.vec) {
    # IOAPI wants names with spaces appended to length=16
#    correct.name <- sprintf("%-16s", name)
#    correct.names.string <- correct.names.string + correct.name
    correct.names.string <- sprintf("%s%-16s", correct.names.string, name)
  }
  target.file <- fix.IOAPI.global.attribute(
    target.file, attr.name="VAR-LIST", value.correct=correct.names.string)
  # return
  target.file
} # end function fix.VARdashLIST

# Note VGLVLS seems to want to have
# * length(VGLVLS) == n(layers) +1
# * first member == 1.0 (or "1.f" in `ncdump`)
create.VGLVLS <- function(
  layers.n  # number of layers the created vector should reflect
) {
  # constants we'll use in calculation, suggested by values in 5yravg.test.nc:VGLVLS
  vglvls.first <- 1.0
  vglvls.from <- 0.1
  vglvls.by   <- 0.01
  correct.VGLVLS <- double(0)
  # seq(from, to, by) -> ((to - from) / by) +1 == length(seq) == layers.n
  # (to - from) / by == layers.n - 1
  # (to - from) == by * (layers.n - 1)
  # to == (by * (layers.n - 1)) + from
  if        (layers.n < 1) {
    stop(sprintf('ioapi.r: (layers.n==%i) < 1\n',
      layers.n))
  } else if (layers.n == 1) {
    correct.VGLVLS <- c(vglvls.first, vglvls.from)
  } else { # layers.n > 1
    vglvls.to <- (vglvls.by * (layers.n - 1)) + vglvls.from
    correct.VGLVLS <- c(vglvls.first,
                        seq(from=vglvls.from, to=vglvls.to, by=vglvls.by))
# start debugging
# cat(sprintf('ioapi.r:create.VGLVLS: from==%f, to==%f, by==%f\n',
#   vglvls.from, vglvls.to, vglvls.by))
#       print('                       correct.VGLVLS==')
# print(correct.VGLVLS)
#   end debugging

  } # end testing layers.n

# debugging
# print('ioapi.r:create.VGLVLS: correct.VGLVLS==') ; print(correct.VGLVLS)

  if (length(correct.VGLVLS) != (layers.n +1)) {
    stop(sprintf('ioapi.r: (length(correct.VGLVLS)==%i) != ((layers.n +1)==%i)\ncorrect.VGLVLS==%s\n',
      length(correct.VGLVLS), (layers.n +1), correct.VGLVLS))
  }
  # else, return
  correct.VGLVLS
} # end function create.VGLVLS

# testcases for create.VGLVLS:
# > source('./ioapi.r')
# > create.VGLVLS(-1)
# Error in create.VGLVLS(-1) : ioapi.r: (layers.n==-1) < 1
# > create.VGLVLS(0)
# Error in create.VGLVLS(0) : ioapi.r: (layers.n==0) < 1
# > create.VGLVLS(1)
# [1] 1.0 0.1
# > create.VGLVLS(2)
# [1] 1.00 0.10 0.11
# > create.VGLVLS(3)
# [1] 1.00 0.10 0.11 0.12
# > length(create.VGLVLS(42))
# [1] 43

# if IOAPI global attribute name=VGLVLS does not contain as many members
# as we think it should, make it so
fix.VGLVLS <- function(
  source.file,  # class=ncdf4 from nc_open(...)
  layers.n,     # number of layers this vector should reflect
  target.file
) {
  # TODO: more testing of inputs

  target.file <- fix.IOAPI.global.attribute(
    target.file, attr.name="VGLVLS", value.correct=create.VGLVLS(layers.n))
  # return
  target.file
} # end function fix.VGLVLS

# if NCO does not copy the FillValue from the source datavar, make it so
fix.FillValue <- function(
  source.file,  # class=ncdf4 from nc_open(...)
  names.vec,    # character vector of names of datavars we want to keep,
                # in the order we want to keep them!
  target.file
) {
  # TODO: more testing of inputs
  stopifnot(!is.element("TFLAG", names.vec)) # ASSERT: TFLAG not a member

  for (name in names.vec) {
    # check return types: doc unclear
    FillValue <- ncatt_get(source.file, varid=datavar.name, attname="_FillValue", verbose=FALSE)
    target.file <- fix.IOAPI.global.attribute(
      target.file, attr.name="_FillValue", value.correct=FillValue)
  }
  # return
  target.file
} # end function fix.FillValue

# The following function assumes that target.file has been created by
# 1 cp source.file target.file
# 2 removing data variables from target.file
# That works for netCDF, but not for I/O API. Among other things,
# IOAPI requires adjustments to its special "datavar""=TFLAG,
# which this function provides.

# IOAPI "datavar"=TFLAG has structure like:

# in `ncdump -h`:
# > dimensions:
# >     TSTEP = UNLIMITED ; // (1 currently)
# >     DATE-TIME = 2 ;
# >     LAY = 42 ;
# >     VAR = 29 ;
# >     ROW = 299 ;
# >     COL = 459 ;
# > variables:
# >     int TFLAG(TSTEP, VAR, DATE-TIME) ;
# >             TFLAG:units = "<YYYYDDD,HHMMSS>" ;
# >             TFLAG:long_name = "TFLAG           " ;
# >             TFLAG:var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;

# in source.file$var[["TFLAG"]]:
# > attr(,"class")
# > [1] "ncid4"
#
# > $name
# > [1] "TFLAG"
#
# > $ndims
# > [1] 3
#
# > $natts
# > [1] 3
#
# > $size
# > [1]  2 29  1
# ==DATE-TIME, VAR, TSTEP

# TODO: handle missing values explicitly!
fix.TFLAG <- function(
  source.file,  # class=ncdf4 from nc_open(...)
  names.vec,    # character vector of names of datavars we want to keep,
                # in the order we want to keep them! NOT including "TFLAG"
  timesteps.n,  # number of timesteps, should be same for all datavars
  target.file,
  target.fp     # kludge: debugging dies due to suicidal target.file?
) {
  # TODO: more testing of inputs
  stopifnot(!is.element("TFLAG", names.vec)) # ASSERT: TFLAG not a member

  source.datavar <- ncvar_get(source.file, varid="TFLAG")
  target.datavar <- ncvar_get(target.file, varid="TFLAG")

  # get indices of the datavars we want to keep, in source order
  source.datavar.names <- get.datavar.names.in.order(source.file)
#  # source.datavar.names includes datavar=TFLAG, which is *not* in itself!
#  # i.e., not in either *.datavar. So ...
#  stopifnot(source.datavar.names[1] == "TFLAG") # ASSERT
#  datavar.names <- source.datavar.names[2:(length(source.datavar.names))]
  datavar.names <- source.datavar.names
  # ... and from that get the indices of the datavars we want to keep
  source.indices.vec <- numeric(0)
# start debugging
#  cat(sprintf('fix.TFLAG:\n'))
#  cat(sprintf('\tnames.vec=')) ; print(names.vec)
#  cat(sprintf('\tdatavar.names=')) ; print(datavar.names)
#   end debugging
  for (name in names.vec) {
    source.indices.vec[length(source.indices.vec) +1] <-
      which(datavar.names == name, arr.ind = TRUE)
  }
  stopifnot(length(names.vec) == length(source.indices.vec)) # ASSERT

  # Now transfer the TFLAG data for the datavars we want in target.
  datavar.dims.n <- length(dim(source.datavar))
  datavar.datetimes.n <- dim(source.datavar)[1]
  source.datavar.vars.n <- dim(source.datavar)[2]
  target.datavar.vars.n <- length(names.vec)

  # Following are used in reading one timestep at a time
  # (Pierce-style read: see help(ncvar_get)#Examples)
  start <- rep(1,datavar.dims.n) # start=(1,1,1,...)
  # Remember timelike dim is always the LAST dimension!
  # for TFLAG: DATE-TIME, VAR, TSTEP
  count.source.raw <- dim(source.datavar)
  count.target.raw <- dim(target.datavar)
  count.target.raw[2] <- target.datavar.vars.n
  # but if val=1, timestep index is omitted from dim(source.datavar),
  # breaking Pierce-style read (below)
  if      (datavar.dims.n < 2) {
    # TODO: throw
    cat(sprintf('ERROR: ioapi.r: datavar.dims.n==%i\n', datavar.dims.n))
  } else if (datavar.dims.n == 2) {
    count.source <- c(count.source.raw, timesteps.n)
    count.target <- c(count.target.raw, timesteps.n)
    start <- c(start, timesteps.n)
    datavar.dims.max.vec <- count.source
    datavar.dims.n <- 3
  } else if (datavar.dims.n == 3) {
    count.source <- count.source.raw
    datavar.dims.max.vec <- count.source
    count.target <- count.target.raw
  } else {
    # TODO: throw
    cat(sprintf('ERROR: ioapi.r: datavar.dims.n==%i\n', datavar.dims.n))
  }

# start debugging
#print('fix.TFLAG: count.source==')
#print(count.source)
#print('fix.TFLAG: count.target==')
#print(count.target)
#   end debugging

  # for safety (and pedagogy), read in data one timestep at a time, dim-agnostically
  for (i.timestep in 1:timesteps.n) {
  #i.timestep <- 1

    # Initialize start and count.source to read one timestep of the source datavar.
    # start=(1,1,i), source.count=(DATE-TIME,VAR,i)
    start[datavar.dims.n] <- i.timestep
    count.source[datavar.dims.n] <- i.timestep
    count.target[datavar.dims.n] <- i.timestep

# start debugging
#cat(sprintf('fix.TFLAG: for timestep==%i, count.source==\n', i.timestep))
#print(count.source)
#cat(sprintf('fix.TFLAG: for timestep==%i, count.target==\n', i.timestep))
#print(count.target)
#   end debugging

#    source.timestep <- ncvar_get(source.file, varid=source.datavar, start=start, count=count.source)
# should work, but class(source.datavar) == "matrix" ???
    source.timestep <- ncvar_get(source.file, varid="TFLAG", start=start, count=count.source)
    # read same size (count.source) in target, later truncating
#    target.timestep.raw <- ncvar_get(target.file, varid=target.datavar, start=start, count=count.source)
    target.timestep.raw <- ncvar_get(target.file, varid="TFLAG", start=start, count=count.source)

# start debugging
#cat(sprintf('fix.TFLAG: for timestep==%i, target.timestep.raw==\n', i.timestep))
#print(target.timestep.raw)
#   end debugging

    # Now "transpose" the appropriate values for the datavars we want:
    # e.g., datavar that was 17th in source is 1st in target.
    # Iterating over length(names.vec):
    for (i.name.target in 1:target.datavar.vars.n) {
      i.name.source <- source.indices.vec[i.name.target]
      source.datetime <- source.timestep[,i.name.source]
      target.timestep.raw[,i.name.target] <- source.datetime
    }

    # now write back new'n'improved (and resized) target.timestep
    target.timestep <- target.timestep.raw
#    dim(target.timestep) <- count.target
# count.target is 3d, but a single timestep is 2d
#    dim(target.timestep) <- count.target[1:2]
# can't do that either: gotta truncate the matrix manually ?!?
    # TODO: find the Rish way to do this!
    as.vector(target.timestep) -> target.timestep.vec
    length(target.timestep.vec) <- 2
    dim(target.timestep.vec) <- c(2,1)
    stopifnot(is.matrix(target.timestep.vec)) # ASSERT
    target.timestep <- target.timestep.vec

# start debugging
cat(sprintf('fix.TFLAG: for timestep==%i, target.timestep==\n', i.timestep))
print(target.timestep)
#   end debugging

#    ncvar_put(target.file, varid=target.datavar, vals=target.timestep, start=start, count=count.target)
#    ncvar_put(target.file, varid="TFLAG", vals=target.timestep, start=start, count=count.target)
    target.file <- ncvar_put(target.file, varid="TFLAG", vals=target.timestep, start=start, count=count.target)

  } # end for timestep

# start debugging
# something's happening to target.file before next call, so kludge :-(
cat(sprintf('fix.TFLAG: kludging target.file from target.fp==%s\n', target.fp))
target.file <- nc_open(target.fp, write=FALSE, readunlim=FALSE)

target.datavar <- ncvar_get(target.file, varid="TFLAG")
print('fix.TFLAG: target.datavar==')
print(target.datavar)
print('fix.TFLAG: exiting')
#   end debugging

  # return
  target.file
} # end function fix.TFLAG

overwrite.layer.in.timestep <- function(
  timestep,               # to which the layer belongs
  value.to.write=NA,      # e.g.
  dims.vec=dim(timestep), # assumed to be COL, ROW, LAY
  i.lay=dim(timestep)[3]  # index of layer to overwrite
) {
# start debugging
#  cat(sprintf('overwrite.layer.in.timestep: layer to fix==%i\n', i.lay))
#  cat(sprintf('overwrite.layer.in.timestep: start: sum(!is.na(layer))==%i\n',
#    sum(!is.na(timestep[,,i.lay]))))
#  print('overwrite.layer.in.timestep: timestep dimensions==')
#  print(dims.vec)
#   end debugging
  layer.dims <- dims.vec[1:2] # LAY assumed to be 3rd dimension
  layer <- rep(value.to.write, prod(layer.dims))
  dim(layer) <- layer.dims
  timestep[,,i.lay] <- layer
# debugging
#  cat(sprintf('overwrite.layer.in.timestep:   end: sum(!is.na(layer))==%i\n',
#    sum(!is.na(timestep[,,i.lay]))))
  # note: explicit `return` halts execution!
  timestep
} # end function overwrite.layer.in.timestep

# get the names of the datavars defined in the file (including TFLAG),
# in the order in which they are defined, return as character vector
get.datavar.names.in.order <- function(
  ioapi.file  # class=ncdf4 from nc_open(...)
) {
  vars.n <- ioapi.file$nvars
  ret.vec <- character(0)
  for (i.var in 1:vars.n) {
    ret.vec[length(ret.vec) +1] <- ioapi.file$var[[ i.var ]]$name
  }
  # return
  ret.vec
}
